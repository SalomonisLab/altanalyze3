"""
Build a genomic HLA allele-signature database for a given genome build.

For each class-I gene (HLA-A/-B/-C) we project every IMGT/HLA allele's
exon-2/3/4 sequence onto genomic coordinates, so the runtime typer only has to
pileup the BAM at a fixed set of chr6 positions and match observed bases.

Anchoring (self-validating):
  1. Reconstruct the reference transcript's CDS from a GTF + chr6 FASTA
     (strand-aware) -> genome_cds, with a cDNA-index -> genomic-coord map.
  2. Find the single IMGT allele G whose exon-2/3/4 sequence is IDENTICAL to the
     genome's exon-2/3/4 sequence.  G is the allele carried by the reference
     genome (e.g. A*03:01 for GRCh38).  An exact match proves the coordinate
     mapping is correct; if none is found the build fails loudly.
  3. Walk the IMGT alignment columns of exon 2/3/4.  Each column at which G has
     a base corresponds to one genomic position; record its coordinate, the
     reference (forward-strand) base, and every allele's expected forward base.

Output: build/signatures_<build>.json.gz
"""

import os, json, gzip, re
from collections import OrderedDict, defaultdict
import pysam
from imgt_parser import parse_nuc_alignment, two_field

_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def rc(s):  return s.translate(_COMP)[::-1]
def comp(b): return b.translate(_COMP)

# exons to include in the signature (class-I antigen recognition domain + a3)
SIG_EXONS = (2, 3, 4)          # 1-based exon numbers
# expression-suffix null alleles to drop from typing candidates
_NULLISH = re.compile(r'[NLSCAQ]$')


def parse_gtf_cds(gtf_path, genes):
    """Return {gene: [ (transcript_id, strand, chrom, [(start,end),...cds]) ]}."""
    def attr(s, k):
        m = re.search(k + r' "([^"]+)"', s); return m.group(1) if m else None
    tx = defaultdict(lambda: {'cds': [], 'gene': None, 'strand': None, 'chrom': None, 'tags': set()})
    with open(gtf_path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'CDS':
                continue
            att = f[8]; tid = attr(att, 'transcript_id')
            tx[tid]['cds'].append((int(f[3]), int(f[4])))
            tx[tid]['gene'] = attr(att, 'gene_name'); tx[tid]['strand'] = f[6]; tx[tid]['chrom'] = f[0]
            tx[tid]['tags'] = set(re.findall(r'tag "([^"]+)"', att))
    out = defaultdict(list)
    for tid, d in tx.items():
        if d['gene'] in genes and d['cds']:
            out[d['gene']].append((tid, d['strand'], d['chrom'], sorted(d['cds']), d['tags']))
    return out


def build_genome_cds(fasta, chrom, strand, cds_intervals):
    """Concatenate CDS in transcription order -> (cds_seq, gpos) where
    gpos[i] is the 0-based genomic coordinate of cDNA base i (forward strand)."""
    seq = []; gpos = []
    ivs = sorted(cds_intervals)                      # ascending genomic
    if strand == '-':
        ivs = ivs[::-1]
    for (s, e) in ivs:                               # s,e are 1-based inclusive
        block = fasta.fetch(chrom, s - 1, e).upper() # forward strand
        if strand == '+':
            for j, b in enumerate(block):
                seq.append(b); gpos.append(s - 1 + j)
        else:
            for j in range(len(block) - 1, -1, -1):
                seq.append(comp(block[j])); gpos.append(s - 1 + j)
    return ''.join(seq), gpos


def ungapped_exon(aln_seq, exons, exon_no):
    """Return exon `exon_no` (1-based) of an aligned allele with gaps removed."""
    s, e = exons[exon_no - 1]
    return aln_seq[s:e].replace('.', '').replace('*', '')


def reference_exon_lengths(parsed):
    """Ungapped exon lengths of the alignment reference allele."""
    ref = parsed['alleles'][parsed['reference']]
    return [len(ref[s:e].replace('.', '')) for (s, e) in parsed['exons']]


def build_gene(gene, parsed, fasta, chrom_name, gtf_cds, verbose=True):
    """Build the signature dict for one gene on one genome build."""
    exons = parsed['exons']
    alleles = parsed['alleles']
    reflens = reference_exon_lengths(parsed)          # per-exon ungapped len (canonical)

    # --- reference genome exon-2/3/4 sequence (transcription orientation) ---
    # try each candidate transcript, prefer CCDS/basic, pick the one that
    # yields an exact IMGT match over the signature exons.
    cds_pos_of_exon = {}
    # cumulative ungapped CDS offsets per exon (canonical lengths)
    cum = [0]
    for L in reflens:
        cum.append(cum[-1] + L)
    best = None
    cand = sorted(gtf_cds[gene], key=lambda t: (0 if 'CCDS' in t[4] else 1, -sum(e - s + 1 for s, e in t[3])))
    genome_exonseq = {}
    for (tid, strand, chrom, cds_iv, tags) in cand:
        cds_seq, gpos = build_genome_cds(fasta, chrom_name, strand, cds_iv)
        # slice exon k from genome CDS using canonical cumulative offsets
        ok = True; exseq = {}; exgpos = {}
        for k in SIG_EXONS:
            a, b = cum[k - 1], cum[k]
            if b > len(cds_seq):
                ok = False; break
            exseq[k] = cds_seq[a:b]
            exgpos[k] = gpos[a:b]
        if not ok:
            continue
        # verify exon k lengths match canonical (no length polymorphism)
        if any(len(exseq[k]) != reflens[k - 1] for k in SIG_EXONS):
            continue
        genome_exonseq = exseq
        best = (tid, strand, chrom, exgpos)
        # will validate by IMGT exact match below; accept first that matches
        gseq = ''.join(exseq[k] for k in SIG_EXONS)
        # find matching IMGT allele over signature exons
        match = None
        for name, aseq in alleles.items():
            cat = ''.join(ungapped_exon(aseq, exons, k) for k in SIG_EXONS)
            if cat == gseq:
                match = name; break
        if match:
            best = (tid, strand, chrom, exgpos, match)
            break
    if best is None or len(best) < 5:
        raise RuntimeError(f'{gene}: no transcript yielded an exact IMGT exon-2/3/4 match '
                           f'(coordinate anchoring failed)')
    tid, strand, chrom, exgpos, G = best
    if verbose:
        print(f'  {gene}: transcript={tid} strand={strand} genome_allele={G}')

    # --- walk alignment columns of signature exons; map G's bases -> gpos ----
    positions = []            # list of (genomic_pos_0based, ref_forward_base)
    col_index = []            # alignment column for each position (for allele lookup)
    Gseq = alleles[G]
    for k in SIG_EXONS:
        s, e = exons[k - 1]
        gp = exgpos[k]; j = 0
        for c in range(s, e):
            gch = Gseq[c]
            if gch in 'ACGT':                         # G has a base here
                gpos0 = gp[j]
                mrna_base = gch
                fwd = mrna_base if strand == '+' else comp(mrna_base)
                # sanity: forward base must equal the FASTA reference base
                fa = fasta.fetch(chrom_name, gpos0, gpos0 + 1).upper()
                if fa != fwd:
                    raise RuntimeError(f'{gene} col{c}: FASTA {fa} != expected {fwd} '
                                       f'@chr6:{gpos0+1} (anchoring bug)')
                positions.append((gpos0, fwd)); col_index.append(c)
                j += 1
            elif gch == '.':
                continue                              # insertion in others; G lacks it
            else:
                # G unknown at a canonical position: still advances genome coord
                j += 1

    # --- 2-field representatives: fullest exon-2/3/4 coverage per 2-field ----
    groups = defaultdict(list)
    for name in alleles:
        groups[two_field(name)].append(name)

    def coverage(name):
        aseq = alleles[name]
        return sum(1 for c in col_index if aseq[c] in 'ACGT')

    allele_bases = OrderedDict()
    for tf, members in groups.items():
        rep = max(members, key=coverage)
        aseq = alleles[rep]
        s = []
        for c in col_index:
            ch = aseq[c]
            if ch in 'ACGT':
                s.append(ch if strand == '+' else comp(ch))
            else:
                s.append('.')                         # no-call
        allele_bases[tf] = ''.join(s)

    # --- keep only columns polymorphic among 2-field representatives ---------
    n = len(positions)
    reps = list(allele_bases.values())
    keep = []
    for i in range(n):
        seen = set(r[i] for r in reps if r[i] in 'ACGT')
        if len(seen) >= 2:
            keep.append(i)
    positions_k = [positions[i] for i in keep]
    allele_bases_k = OrderedDict(
        (tf, ''.join(s[i] for i in keep)) for tf, s in allele_bases.items())

    if verbose:
        print(f'       positions(all exon2/3/4)={n}  polymorphic={len(keep)}  '
              f'2field_alleles={len(allele_bases_k)}')

    return {
        'gene': gene, 'chrom': chrom, 'strand': strand, 'transcript': tid,
        'genome_allele': G,
        'positions': [p for p, _ in positions_k],
        'ref_base': ''.join(b for _, b in positions_k),
        'allele_bases': allele_bases_k,
    }


def build_signatures(build, chr6_fasta, gtf_path, imgt_dir, out_path, verbose=True):
    genes = ('HLA-A', 'HLA-B', 'HLA-C')
    gene_letter = {'HLA-A': 'A', 'HLA-B': 'B', 'HLA-C': 'C'}
    fasta = pysam.FastaFile(chr6_fasta)
    chrom_name = fasta.references[0]                   # 'chr6'
    gtf_cds = parse_gtf_cds(gtf_path, set(genes))
    sig = {'build': build, 'chrom': chrom_name, 'genes': {}}
    for g in genes:
        parsed = parse_nuc_alignment(os.path.join(imgt_dir, f'{gene_letter[g]}_nuc.txt'))
        sig['genes'][g] = build_gene(g, parsed, fasta, chrom_name, gtf_cds, verbose)
    with gzip.open(out_path, 'wt') as fh:
        json.dump(sig, fh)
    if verbose:
        print(f'  wrote {out_path}')
    return sig


if __name__ == '__main__':
    HERE = os.path.dirname(os.path.abspath(__file__))
    imgt = os.path.join(HERE, 'data', 'imgt', 'alignments')
    builds = [
        ('hg38', 'hg38.chr6.fa', 'hla_hg38.ensembl104.gtf'),
        ('hg19', 'hg19.chr6.fa', 'hla_hg19.ensembl87.gtf'),
    ]
    for build, fa, gtf in builds:
        fap = os.path.join(HERE, 'data', 'genome', fa)
        gtfp = os.path.join(HERE, 'build', gtf)
        if not (os.path.exists(fap) and os.path.exists(gtfp)):
            print(f'=== skipping {build} (missing {fa} or {gtf}) ===')
            continue
        print(f'=== building {build} signatures ===')
        build_signatures(build, fap, gtfp, imgt,
                         os.path.join(HERE, 'build', f'signatures_{build}.json.gz'))
