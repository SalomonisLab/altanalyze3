"""Can the KINNEX reads separate a linear KMT2A-PTD transcript from a KMT2A circRNA? (2026-11-10, v2). Ad hoc;
validate within altanalyze3.

Both molecules give the same backward junction: donor exon d joined to an upstream acceptor exon a.
  circRNA  every base lies inside the circle, genomic [acceptor, donor]. A read crosses the junction more than once
           only when it is longer than the circle; a read 2x the circle length must cross it at least twice.
  PTD      the molecule is linear with 2 tandem copies of the block. Exon 1 lies upstream of the first copy and the
           exons after d lie downstream of the second copy. A read must be longer than the block (spliced) to reach
           them, and it crosses the junction once.
Test 1, the PTD-type reads of DOCUMENTED_LESION_READ_IDS.tsv:
  linear               aligned KMT2A sequence outside [acceptor, donor] by more than 15 nt
  circle (wrap)        2 or more backward joins
  uninformative        neither of the above
  circle_spliced_length_nt counts the MANE exon bases inside [acceptor, donor], partial exons included.
Test 2, every backward join in the KMT2A region reads. Each end is typed against all GENCODE v45 KMT2A exons:
  boundary (within 15 nt of an exon start for the acceptor, an exon end for the donor; MANE label when MANE),
  exon interior, or intron. An intron-to-intron join with the donor in introns 8-11 and the acceptor in introns 1-2
  has the shape of the genomic PTD breakpoint, which a pre-mRNA can carry and a back-splice cannot.
v1 (delete/kmt2a_ptd_vs_circrna_v1_superseded/) typed sites against MANE exons only, which labelled an annotated
non-MANE exon in intron 1 as intron; it left partial exons out of the circle length; and it stated that a read
longer than the circle must cross the junction twice, which is false.
v3 reads the PTD-type reads from the export named on the command line (DOCUMENTED_LESION_READ_IDS_v2.tsv, which
adds the alt-e2 acceptor); the logic is unchanged from v2 (delete/kmt2a_ptd_vs_circrna_v2_superseded/).
Usage: kmt2a_ptd_vs_circrna_v3.py <scan_dir> <export tsv> <MANE exon gff3 subset> <all KMT2A exons gff3> <out_reads.tsv> <out_joins.tsv>
"""
import collections, csv, os, re, sys

O, EXPORT, MANE, ALLEX, OUT1, OUT2 = sys.argv[1:7]
for p in (OUT1, OUT2):
    assert not os.path.exists(p), f'{p} exists; refusing to overwrite'
csv.field_size_limit(sys.maxsize)
EX = {}
for line in open(MANE):
    if 'gene_name=KMT2A;' in line and '\texon\t' in line and 'MANE_Select' in line:
        c = line.split('\t'); EX[int(re.search(r'exon_number=(\d+)', c[8]).group(1))] = (int(c[3]), int(c[4]))
assert sorted(EX) == list(range(1, 37)), sorted(EX)
starts, ends, interiors = collections.defaultdict(set), collections.defaultdict(set), []
for line in open(ALLEX):
    c = line.rstrip('\n').split('\t')
    if c[0] != 'chr11' or c[2] != 'exon' or 'gene_name=KMT2A;' not in c[8]:
        continue
    a, b = int(c[3]), int(c[4]); tx = re.search(r'transcript_name=([^;]+)', c[8]).group(1)
    starts[a].add(tx); ends[b].add(tx)
    if 'transcript_type=retained_intron' not in c[8]:
        interiors.append((a, b))
print('MANE exons 36 | GENCODE KMT2A exon starts', len(starts), 'ends', len(ends))
TOL = 15
MSTART = {a: n for n, (a, b) in EX.items()}
MEND = {b: n for n, (a, b) in EX.items()}


def exon_of(pos):
    for n, (a, b) in EX.items():
        if a - TOL <= pos <= b + TOL:
            return n
    return None


def site(pos, kind):
    ref = ends if kind == 'donor' else starts
    mane = MEND if kind == 'donor' else MSTART
    best = min(ref, key=lambda p: abs(p - pos))
    if abs(best - pos) <= TOL:
        if best in mane:
            return f'e{mane[best]} boundary', 'boundary'
        return f"non-MANE exon {'end' if kind == 'donor' else 'start'} {best} ({','.join(sorted(ref[best])[:3])})", 'boundary'
    if exon_of(pos) or any(a <= pos <= b for a, b in interiors):
        return 'exon interior', 'exon_interior'
    return 'intron', 'intron'


def intron_no(p):
    for k in range(1, 36):
        if EX[k][1] < p < EX[k + 1][0]:
            return k
    return None


SEG = re.compile(r'(\S+?):(\d+)-(\d+)\(([+-])\) q(\d+)-(\d+) mq(\d+) \[([^\]]*)\]')


def kmt2a_segments(s):
    return [dict(chrom=m[0], start=int(m[1]), end=int(m[2]), strand=m[3], qs=int(m[4]), qe=int(m[5]))
            for m in SEG.findall(s) if 'KMT2A' in m[7].split(',')]


def backward_joins(ks):
    out = []
    for a, b in zip(ks, ks[1:]):
        if a['chrom'] != b['chrom'] or a['strand'] != b['strand']:
            continue
        if a['strand'] == '+' and b['start'] < a['end']:
            out.append((a['end'], b['start']))
        elif a['strand'] == '-' and b['end'] > a['start']:
            out.append((b['end'], a['start']))
    return out


# ---- test 1
rows = [r for r in csv.DictReader(open(EXPORT), delimiter='\t') if r['lesion'] == 'KMT2A-PTD']
print('PTD-type reads in', EXPORT, len(rows))
res = []
for r in rows:
    ks = kmt2a_segments(r['segments'])
    bj = backward_joins(ks)
    assert bj, r['read_id']
    blo, bhi = min(a for _, a in bj), max(d for d, _ in bj)
    lo, hi = min(s['start'] for s in ks), max(s['end'] for s in ks)
    up, down = max(0, blo - lo), max(0, hi - bhi)
    circ = sum(max(0, min(b, bhi) - max(a, blo) + 1) for a, b in EX.values())
    qlen = max(s['qe'] for s in ks) - min(s['qs'] for s in ks)
    if up > TOL or down > TOL:
        cls = 'linear: sequence outside the joined block'
    elif len(bj) >= 2:
        cls = 'circle-compatible: 2 or more backward joins'
    else:
        cls = 'uninformative: 1 join, all sequence inside the block'
    res.append(dict(row_role=r['row_role'], sample=r['sample'], library=r['library'], group=r['group'], bam_index=r['bam_index'],
                    read_id=r['read_id'], event=r['event'], n_backward_joins=len(bj), block=f'chr11:{blo}-{bhi}',
                    circle_spliced_length_nt=circ, aligned_kmt2a_query_nt=qlen,
                    read_longer_than_block=qlen > circ, read_ge_2x_circle=qlen >= 2 * circ,
                    upstream_of_block_nt=up, upstream_end=(f'e{exon_of(lo)}' if exon_of(lo) else 'intron') if up > TOL else '',
                    downstream_of_block_nt=down, downstream_end=(f'e{exon_of(hi)}' if exon_of(hi) else 'intron') if down > TOL else '',
                    read_class=cls, segments=r['segments']))
with open(OUT1, 'w', newline='') as h:
    w = csv.DictWriter(h, fieldnames=list(res[0]), delimiter='\t'); w.writeheader(); w.writerows(res)
print('wrote', OUT1, len(res), 'rows')
grp = lambda x: 'patient (U2AF1-RAS)' if x['row_role'] == 'lesion_sample' else ('healthy' if x['group'] in ('young', 'aged') else 'other disease/iPSC')
for k, n in sorted(collections.Counter((grp(x), x['read_class'], 'longer than block' if x['read_longer_than_block'] else 'shorter than block') for x in res).items()):
    print('  ', k, n)
print('  reads >= 2x circle length:', sum(x['read_ge_2x_circle'] for x in res))
for x in res:
    if x['read_class'].startswith('linear') or x['read_longer_than_block']:
        print('   ', grp(x), x['sample'], x['read_id'], x['event'], x['read_class'][:9], 'circle', x['circle_spliced_length_nt'], 'read', x['aligned_kmt2a_query_nt'],
              'up', x['upstream_of_block_nt'], x['upstream_end'], 'down', x['downstream_of_block_nt'], x['downstream_end'])
q = sorted(x['aligned_kmt2a_query_nt'] for x in res if grp(x).startswith('patient'))
print('patient aligned KMT2A query nt: min', q[0], 'median', q[len(q) // 2], 'max', q[-1], '| e2-e10 spliced:', sum(EX[n][1] - EX[n][0] + 1 for n in range(2, 11)))

# ---- test 2
reg = [r for r in csv.DictReader(open(os.path.join(O, 'KMT2A_PTD_and_FLT3_ITD_region_noncolinear_reads.tsv')), delimiter='\t') if r['gene'] == 'KMT2A']
PAT = {'MDS-pre-U2AF1-RAS__Y2982', 'MDS-post-U2AF1-RAS__Y3554', 'MDS-pre-U2AF1-RAS-MISTRG__Y2982-E3'}
jr = []
for r in reg:
    for d, a in backward_joins(kmt2a_segments(r['segments'])):
        ds, dk = site(d, 'donor'); as_, ak = site(a, 'acceptor')
        di, ai = intron_no(d), intron_no(a)
        jr.append(dict(uid=r['uid'], group=r['group'], patient=r['uid'] in PAT, read=r['read'], donor_pos=d, donor_site=ds,
                       acceptor_pos=a, acceptor_site=as_, join_type=f'{dk}>{ak}',
                       ptd_breakpoint_shape=dk == 'intron' and ak == 'intron' and di in (8, 9, 10, 11) and ai in (1, 2),
                       segments=r['segments']))
with open(OUT2, 'w', newline='') as h:
    w = csv.DictWriter(h, fieldnames=list(jr[0]), delimiter='\t'); w.writeheader(); w.writerows(jr)
print('\nwrote', OUT2, len(jr), 'backward joins from', len({(x['uid'], x['read']) for x in jr}), 'of', len(reg), 'region reads')
g2 = lambda x: 'patient' if x['patient'] else ('healthy' if x['group'] in ('young', 'aged') else 'other')
for k, n in sorted(collections.Counter((g2(x), x['join_type']) for x in jr).items()):
    print('  ', k, n)
print('boundary>boundary joins by donor > acceptor, per group (top 8):')
for gname in ('patient', 'healthy', 'other'):
    c = collections.Counter(f"{x['donor_site']} > {x['acceptor_site']}" for x in jr if g2(x) == gname and x['join_type'] == 'boundary>boundary')
    print('  ', gname, sum(c.values()), c.most_common(8))
bp = [x for x in jr if x['ptd_breakpoint_shape']]
print('PTD breakpoint-shaped joins (intron 8-11 > intron 1-2):', len(bp), collections.Counter(g2(x) for x in bp))
