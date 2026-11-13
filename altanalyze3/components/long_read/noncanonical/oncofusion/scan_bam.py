"""AML/MDS fusion and FLT3-ITD read scan of one KINNEX BAM (2026-11-10). Ad hoc; validate within altanalyze3.

Request (Nathan Salomonis, 2026-11-10): any trans-splicing read (inter- or intra-chromosomal, no splice-site
requirement) that corresponds to a documented AML/MDS oncofusion (KMT2A/MLL, RUNX1, CBFB::MYH11, BCR::ABL1, PML::RARA ...),
and FLT3-ITD reads (tandem duplications; no exact database match required).

Fusion reads: every non-secondary alignment overlapping an AML/MDS fusion gene (samtools -L) that has an SA tag
(chimeric: other chromosome, other strand or out of order), or whose own blocks overlap 2 or more genes. The read's
segments (the record plus its SA entries) are annotated with every overlapping gene (GENCODE v45 genes.tsv).
FLT3: primary alignments with an aligned block in FLT3 exons 14-15 (ENST00000241453, MANE); every CIGAR insertion
>= 3 nt (the shortest reported ITD) inside the window, and soft clips >= 20 nt whose clip point lies in the window.
Usage: scan_bam.py <uid> <bam> <outdir> <genes.tsv.gz> <aml_genes.txt>
"""
import collections, csv, gzip, json, os, re, subprocess, sys, time

uid, bam, out, genes_path, aml_path = sys.argv[1:6]
os.makedirs(out, exist_ok=True)
FLT3_WIN = ('chr13', 28033887, 28034214)  # exon 15 start .. exon 14 end, GENCODE v45 ENST00000241453
EXONS_1415 = [(28033887, 28033991), (28034082, 28034214)]

genes = []
with gzip.open(genes_path, 'rt') as h:
    for r in csv.DictReader(h, delimiter='\t'):
        genes.append((r['chrom'], int(r['start']), int(r['end']), r['gene_name'], r['gene_id'], r['strand'], r['gene_type']))
aml = [l.strip() for l in open(aml_path) if l.strip() and not l.startswith('#')]
by_name = collections.defaultdict(list)
for g in genes:
    by_name[g[3]].append(g)
missing = [a for a in aml if a not in by_name]
assert not missing, f'AML genes not in genes.tsv: {missing}'
BIN = 100000
idx = collections.defaultdict(list)
for g in genes:
    for b in range(g[1] // BIN, g[2] // BIN + 1):
        idx[g[0], b].append(g)


def genes_at(ch, s, e):
    hit = {}
    for b in range(s // BIN, e // BIN + 1):
        for g in idx.get((ch, b), []):
            if g[1] <= e and s <= g[2]:
                hit[g[4]] = g
    return list(hit.values())


bed = os.path.join(out, 'regions.bed')
with open(bed, 'w') as h:
    for a in aml:
        for g in by_name[a]:
            if not g[0].startswith('chr') or '_' in g[0]:
                continue
            h.write(f'{g[0]}\t{g[1] - 1}\t{g[2]}\t{g[3]}\n')
    h.write(f'{FLT3_WIN[0]}\t{FLT3_WIN[1] - 1}\t{FLT3_WIN[2]}\tFLT3_ITD_window\n')

CIG = re.compile(r'(\d+)([MIDNSHP=X])')


def block_genes(ch, bl):
    return [{g[3] for g in genes_at(ch, s, e)} for s, e in bl]


def spans_genes(ch, bl):
    """True when the blocks fall in 2 or more genes and no single gene overlaps every block (excludes nested and
    overlapping genes that hold the whole alignment)."""
    bg = block_genes(ch, bl)
    allg = set().union(*bg) if bg else set()
    return len(allg) >= 2 and not any(all(g in b for b in bg) for g in allg)


def blocks(pos, cigar):
    out_, r = [], pos
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in 'M=X':
            out_.append((r, r + n - 1)); r += n
        elif op in 'DN':
            r += n
    return out_


def qspan(cigar, strand):
    ops = CIG.findall(cigar)
    lead = int(ops[0][0]) if ops[0][1] in 'SH' else 0
    trail = int(ops[-1][0]) if ops[-1][1] in 'SH' else 0
    qlen = sum(int(n) for n, op in ops if op in 'MIS=XH')
    aln = sum(int(n) for n, op in ops if op in 'MI=X')
    s = lead if strand == '+' else trail
    return s, s + aln, qlen


segs = collections.defaultdict(dict)
flt3_cov, flt3_events = set(), []
n_rec = n_sa = 0
t0 = time.time()
p = subprocess.Popen(['samtools', 'view', '-F', '0x100', '-L', bed, bam], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
for line in p.stdout:
    f = line.split('\t', 11)
    flag, ch, pos, mapq, cigar = int(f[1]), f[2], int(f[3]), int(f[4]), f[5]
    if cigar == '*':
        continue
    n_rec += 1
    strand = '-' if flag & 16 else '+'
    tags = f[11]
    sa = re.search(r'\tSA:Z:([^\t\n]+)', '\t' + tags)
    bl = blocks(pos, cigar)
    q = f[0]
    # FLT3-ITD window (primary records only)
    if not flag & 0x800 and ch == FLT3_WIN[0] and any(s <= FLT3_WIN[2] and FLT3_WIN[1] <= e for s, e in bl):
        flt3_cov.add(q)
        r, qp, seq = pos, 0, f[9]
        ops = CIG.findall(cigar)
        for k, (n, op) in enumerate(ops):
            n = int(n)
            if op in 'M=X':
                r += n; qp += n
            elif op in 'DN':
                r += n
            elif op == 'I':
                if n >= 3 and any(a <= r <= b + 1 for a, b in EXONS_1415):
                    flt3_events.append(dict(uid=uid, read=q, kind='insertion', length=n, ref_pos=r, in_frame=(n % 3 == 0),
                                            seq=seq[qp:qp + n] if seq != '*' else '', strand=strand, mapq=mapq))
                qp += n
            elif op == 'S':
                if n >= 20 and FLT3_WIN[1] <= r <= FLT3_WIN[2] + 1:
                    flt3_events.append(dict(uid=uid, read=q, kind='softclip_' + ('left' if k == 0 else 'right'), length=n, ref_pos=r,
                                            in_frame='', seq=seq[qp:qp + n][:300] if seq != '*' else '', strand=strand, mapq=mapq))
                qp += n
    # fusion candidates
    multi = spans_genes(ch, bl)
    if not sa and not multi:
        continue
    n_sa += 1
    key = (ch, pos, strand, cigar)
    segs[q][key] = dict(primary=not flag & 0x800, mapq=mapq)
    if sa:
        for ent in sa.group(1).strip().rstrip(';').split(';'):
            c2, p2, s2, cg2, mq2, nm2 = ent.split(',')
            segs[q].setdefault((c2, int(p2), s2, cg2), dict(primary=None, mapq=int(mq2)))
rc = p.wait()
assert rc == 0, f'samtools exited {rc}'

aml_set = set(aml)
rows = []
for q, d in segs.items():
    sl = []
    for (ch, pos, strand, cigar), m in d.items():
        bl = blocks(pos, cigar)
        gs = sorted({(g[3], g[5]) for s, e in bl for g in genes_at(ch, s, e)})
        qs, qe, qlen = qspan(cigar, strand)
        sl.append(dict(chrom=ch, start=bl[0][0], end=bl[-1][1], strand=strand, qstart=qs, qend=qe, mapq=m['mapq'],
                       genes=[g[0] for g in gs], gene_strands=[g[1] for g in gs], blocks=bl, bgenes=block_genes(ch, bl),
                       multi=spans_genes(ch, bl)))
    sl.sort(key=lambda x: x['qstart'])
    names = [g for s_ in sl for g in s_['genes']]
    if not any(g in aml_set for g in names):
        continue
    # gene pairs across segments, in read order; plus within-segment multi-gene blocks
    pairs = set()
    for i in range(len(sl)):
        for j in range(i + 1, len(sl)):
            for ga in sl[i]['genes']:
                for gb in sl[j]['genes']:
                    if ga != gb and (ga in aml_set or gb in aml_set):
                        pairs.add((ga, gb, 'split'))
        if sl[i]['multi']:
            bg = sl[i]['bgenes'] if sl[i]['strand'] == '+' else sl[i]['bgenes']  # genomic order; read order set below
            for k in range(len(bg)):
                for l in range(k + 1, len(bg)):
                    for ga in bg[k] - bg[l]:
                        for gb in bg[l] - bg[k]:
                            if ga in aml_set or gb in aml_set:
                                pairs.add((ga, gb, 'same_alignment_genomic_order'))
    # same-gene non-colinear (e.g. KMT2A partial tandem duplication): a later read segment lies upstream in the gene
    ptd = []
    for i in range(len(sl) - 1):
        a_, b_ = sl[i], sl[i + 1]
        common = set(a_['genes']) & set(b_['genes']) & aml_set
        if common and a_['chrom'] == b_['chrom'] and a_['strand'] == b_['strand']:
            backwards = (b_['start'] < a_['end']) if a_['strand'] == '+' else (b_['end'] > a_['start'])
            if backwards:
                ptd.extend(sorted(common))
    if not pairs and not ptd:
        continue
    rows.append(dict(uid=uid, read=q, n_segments=len(sl),
                     segments=' | '.join(f"{s_['chrom']}:{s_['start']}-{s_['end']}({s_['strand']}) q{s_['qstart']}-{s_['qend']} mq{s_['mapq']} [{','.join(s_['genes']) or 'intergenic'}]" for s_ in sl),
                     gene_pairs=';'.join(f'{a}>{b}:{k}' for a, b, k in sorted(pairs)), same_gene_noncolinear=';'.join(sorted(set(ptd)))))

with open(os.path.join(out, f'{uid}.fusion_reads.tsv'), 'w', newline='') as h:
    w = csv.DictWriter(h, fieldnames=['uid', 'read', 'n_segments', 'segments', 'gene_pairs', 'same_gene_noncolinear'], delimiter='\t')
    w.writeheader(); w.writerows(rows)
with open(os.path.join(out, f'{uid}.flt3_events.tsv'), 'w', newline='') as h:
    w = csv.DictWriter(h, fieldnames=['uid', 'read', 'kind', 'length', 'ref_pos', 'in_frame', 'seq', 'strand', 'mapq'], delimiter='\t')
    w.writeheader(); w.writerows(flt3_events)
json.dump(dict(uid=uid, bam=bam, records_in_regions=n_rec, chimeric_or_multigene_records=n_sa, reads_with_aml_gene_pair_or_noncolinear=len(rows),
               flt3_reads_covering_exons14_15=len(flt3_cov), flt3_reads_with_event=len({e['read'] for e in flt3_events}),
               seconds=round(time.time() - t0, 1)), open(os.path.join(out, f'{uid}.summary.json'), 'w'), indent=1)
print('DONE', uid, n_rec, n_sa, len(rows), len(flt3_cov), len(flt3_events), round(time.time() - t0, 1))
