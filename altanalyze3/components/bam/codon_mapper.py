#!/usr/bin/env python3.11
"""codon_mapper.py — authoritative protein-codon -> hg38 genomic position mapping.

Source of truth: GENCODE v45 CDS blocks (MANE Select transcript preferred) + the hg38
reference FASTA (for reference-codon translation checks). NO coordinates are guessed:
every hotspot position is derived from the CDS structure and verified by translating the
reference codon to the expected wild-type amino acid.

FASTA contig naming: the reference genome.fa here is NOT chr-prefixed (1,2,...,X,Y,MT),
while GENCODE / the BAMs are chr-prefixed. We strip 'chr' for FASTA fetches only.
"""
import pysam

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
_COMP = {'A':'T','T':'A','C':'G','G':'C','N':'N'}


def _rc(s):
    return ''.join(_COMP.get(b, b) for b in reversed(s))


def parse_gff_cds(gff_path):
    """Return {gene_name: {transcript_id: {'strand','chrom','mane':bool,'cds':[(start,end,phase)]}}}."""
    genes = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            chrom, _src, feat, start, end, _score, strand, phase, attr = f[:9]
            if feat not in ('CDS', 'transcript'):
                continue
            A = {}
            for kv in attr.split(';'):
                if '=' in kv:
                    k, v = kv.split('=', 1)
                    A[k] = v
            gn = A.get('gene_name')
            tid = A.get('transcript_id')
            if not gn or not tid:
                continue
            g = genes.setdefault(gn, {})
            t = g.setdefault(tid, {'strand': strand, 'chrom': chrom, 'mane': False, 'cds': []})
            tags = A.get('tag', '')
            if 'MANE_Select' in tags:
                t['mane'] = True
            if feat == 'CDS':
                ph = int(phase) if phase in ('0', '1', '2') else 0
                t['cds'].append((int(start), int(end), ph))
    return genes


def _pick_transcript(gene_tx):
    """MANE Select if present, else the transcript with the longest CDS."""
    mane = [(tid, t) for tid, t in gene_tx.items() if t['mane'] and t['cds']]
    if mane:
        return mane[0]
    cand = [(tid, t) for tid, t in gene_tx.items() if t['cds']]
    if not cand:
        return None, None
    cand.sort(key=lambda kt: -sum(e - s + 1 for s, e, _ in kt[1]['cds']))
    return cand[0]


def _coding_positions(t):
    """Ordered list of genomic positions (1-based) in coding (5'->3') order."""
    cds = sorted(t['cds'], key=lambda x: x[0])
    pos = []
    if t['strand'] == '+':
        for s, e, _ in cds:
            pos.extend(range(s, e + 1))
    else:
        for s, e, _ in reversed(cds):
            pos.extend(range(e, s - 1, -1))
    return pos


def codon_to_genomic(gene_tx, codon_num):
    """Return (chrom, strand, [g1,g2,g3]) genomic 1-based positions of the codon (coding order)."""
    tid, t = _pick_transcript(gene_tx)
    if t is None:
        return None
    pos = _coding_positions(t)
    i0 = (codon_num - 1) * 3
    if i0 + 2 >= len(pos):
        return None
    return t['chrom'], t['strand'], pos[i0:i0 + 3], tid


def ref_codon(fasta, chrom, strand, gpos3):
    """Fetch the reference codon (coding orientation) from the FASTA. FASTA is not chr-prefixed."""
    fc = chrom[3:] if chrom.startswith('chr') else chrom
    bases = [fasta.fetch(fc, p - 1, p).upper() for p in gpos3]
    if strand == '-':
        bases = [_COMP.get(b, b) for b in bases]  # each base complemented; positions already 3'->5' order
    codon = ''.join(bases)
    return codon


if __name__ == '__main__':
    import sys
    GFF = "/Volumes/salomonis-archive/LabFiles/Nathan/Revio/MDS-AML-KINNEX-5/analysis/genomic_variants/variant_discovery/KINNEX_supervised_variant_detection/panels/gencode_v45_driver_cds.gff3"
    FASTA = "/Users/saljh8/Dropbox/Revio/Other/Variants/SNV/genome.fa"
    genes = parse_gff_cds(GFF)
    fa = pysam.FastaFile(FASTA)
    # (gene, codon, expected_wt_aa, known_hg38_pos_or_None)
    tests = [
        ("SRSF2", 95, "P", 76736877),   # ground truth from positive-control CSV
        ("U2AF1", 34, "S", 43104346),   # ground truth (MDS-U2AF1 c.101C>T = S34F)
        ("IDH1", 132, "R", None),
        ("IDH2", 140, "R", None),
        ("IDH2", 172, "R", None),
        ("DNMT3A", 882, "R", None),
        ("JAK2", 617, "V", None),
        ("KIT", 816, "D", None),
        ("FLT3", 835, "D", None),
        ("NRAS", 12, "G", None),
        ("NRAS", 61, "Q", None),
        ("KRAS", 12, "G", None),
        ("TP53", 175, "R", None),
        ("TP53", 248, "R", None),
        ("TP53", 273, "R", None),
        ("PTPN11", 76, "E", None),
        ("WT1", 1, "M", None),
    ]
    print(f"{'gene':8} {'codon':5} {'tid':20} {'strand':6} {'chrom':6} {'g-positions':26} {'refcodon':8} {'refAA':5} {'expAA':5} {'known':10} {'CHECK'}")
    for gene, codon, exp, known in tests:
        if gene not in genes:
            print(f"{gene:8} MISSING"); continue
        r = codon_to_genomic(genes[gene], codon)
        if r is None:
            print(f"{gene:8} {codon:<5} no-mapping"); continue
        chrom, strand, g3, tid = r
        rc = ref_codon(fa, chrom, strand, g3)
        aa = CODON_TABLE.get(rc, '?')
        ok = "OK" if aa == exp else "**AA-MISMATCH**"
        kk = ""
        if known is not None:
            kk = "OK" if known in g3 else f"**POS {known} not in {g3}**"
        print(f"{gene:8} {codon:<5} {tid:20} {strand:6} {chrom:6} {str(g3):26} {rc:8} {aa:5} {exp:5} {str(known):10} {ok} {kk}")
    fa.close()
