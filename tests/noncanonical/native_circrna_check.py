"""Unit tests for forked_altanalyze3/.../bam/noncolinear_rna.py (EXPERIMENTAL fork, 2026-11-05).

Synthetic reads are built from the real genome (genome.fa) on HIPK3 (+ strand) and RUNX1 (- strand) with
the exon-model subset noncolinear_circRNA_20261105/tests/ref/Hs_Ensembl_exon.txt. Each read is a pysam
record with a CIGAR and an SA tag laid out as pbmm2 writes them (transcript-oriented read, soft clips).
Expected values come from the read design, never from the module. The genome helpers are also compared
with the functions in scripts/report_stream.py (extracted by ast) on the test joints and on 300 random
positions. Usage: python test_noncolinear_rna.py [genome.fa] [out.json]
"""
import ast
import json
import random
import sys
from pathlib import Path

import pysam

X = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(X / 'forked_altanalyze3/altanalyze3/components'))
from altanalyze3.components.bam import noncolinear_rna as N  # noqa: E402

GENOME = sys.argv[1] if len(sys.argv) > 1 else '/data/salomonis-archive/LabFiles/Nathan/Revio/Hs/genome.fa'
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else X / 'noncolinear_circRNA_20261105/tests/unit_tests.json'
REF = Path(__file__).parent / 'fixtures/Hs_Ensembl_exon.txt'
fasta = pysam.FastaFile(GENOME)
COMP = str.maketrans('ACGTN', 'TGCAN')
rc = lambda s: s.translate(COMP)[::-1]
g = lambda c, s1, e1: fasta.fetch('MT' if c == 'chrM' else (c[3:] if c[3:] in fasta.references else c), s1 - 1, e1).upper()
header = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'}, 'SQ': [
    {'SN': 'chr11', 'LN': fasta.get_reference_length('11')}, {'SN': 'chr21', 'LN': fasta.get_reference_length('21')},
    {'SN': 'chrM', 'LN': fasta.get_reference_length('MT')}, {'SN': 'chr4', 'LN': fasta.get_reference_length('4')}]})
annotation = N.Annotation(REF)
genome = N.Genome(GENOME)
results, failures = [], []


def piece(q0, q1, chrom, strand, blocks):
    """blocks: genomic 1-based (start, end) in genomic order; total length must equal q1 - q0."""
    assert sum(e - s + 1 for s, e in blocks) == q1 - q0, (q0, q1, blocks)
    return dict(q0=q0, q1=q1, chrom=chrom, strand=strand, blocks=blocks)


def tsense(chrom, strand, blocks):
    seq = ''.join(g(chrom, s, e) for s, e in blocks)
    return seq if strand == '+' else rc(seq)


def cigar(p, L):
    lead, tail = (p['q0'], L - p['q1']) if p['strand'] == '+' else (L - p['q1'], p['q0'])
    ops = [f'{lead}S'] if lead else []
    for i, (s, e) in enumerate(p['blocks']):
        if i:
            ops.append(f"{s - p['blocks'][i - 1][1] - 1}N")
        ops.append(f'{e - s + 1}=')
    if tail:
        ops.append(f'{tail}S')
    return ''.join(ops)


def make_read(name, pieces, primary=0, read=None):
    seq = read if read is not None else ''.join(tsense(p['chrom'], p['strand'], p['blocks']) for p in pieces)
    L = len(seq)
    p = pieces[primary]
    r = pysam.AlignedSegment(header)
    r.query_name = f'molecule/{name}'
    r.flag = 16 if p['strand'] == '-' else 0
    r.reference_name = p['chrom']
    r.reference_start = p['blocks'][0][0] - 1
    r.mapping_quality = 60
    r.cigarstring = cigar(p, L)
    r.query_sequence = seq if p['strand'] == '+' else rc(seq)
    sa = ''.join(f"{o['chrom']},{o['blocks'][0][0]},{o['strand']},{cigar(o, L)},60,0;" for i, o in enumerate(pieces) if i != primary)
    r.set_tags([('SA', sa), ('CB', 'AAACCCAAGAAACACT'), ('zm', name), ('XM', 'UMI'), ('NM', 0)])
    return r


def check(label, got, expected):
    ok = got == expected
    results.append(dict(test=label, got=got, expected=expected, passed=ok))
    if not ok:
        failures.append(label)


def rows_of(read):
    return N.analyze_read(read, annotation, genome, 'TEST')


# ---------------------------------------------------------------- geometry on HIPK3 E7.1 (+), chr11:33307959-33309057
# E7.1 (AK293009-2) is an mRNA-only region with non-canonical flanks; it tests geometry, not biology.
c, s0, e0 = 'chr11', 33307959, 33309057
E7 = (s0, e0)

# T1 single back-splice crossing, contained
r = make_read(1, [piece(0, 499, c, '+', [(s0 + 600, e0)]), piece(499, 899, c, '+', [(s0, s0 + 399)])])
rows = rows_of(r)
check('T1 rows', len(rows), 1)
x = rows[0]
check('T1 joint', (x['joint_type'], x['donor'], x['acceptor']), ('back', e0, s0))
check('T1 annotations', (x['donor_annotation'], x['acceptor_annotation']), ('ENSG00000110422:E7.1', 'ENSG00000110422:E7.1'))
dn, an = g(c, e0 + 1, e0 + 2), g(c, s0 - 2, s0 - 1)
check('T1 dinucleotides = genome', (x['donor_dinucleotide'], x['acceptor_dinucleotide']), (dn, an))
check('T1 class', x['joint_class'], 'backsplice_registered_canonical' if (dn, an) in {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')} else 'backsplice_registered_noncanonical')
check('T1 circle', (x['circle_start'], x['circle_end'], x['circle_span_nt']), (s0, e0, 1099))
check('T1 contained/crossings/overlap', (x['outside_circle_nt'], x['crossings'], x['pair_overlap_nt']), (0, 1, 0))
check('T1 region labels', (x['donor_region'], x['acceptor_region']), ('HIPK3:E7.1', 'HIPK3:E7.1'))
check('T1 same gene, status, unaligned', (x['same_gene'], x['status'], x['unaligned_read_nt']), (True, 'kept', 0))
check('T1 read context joint', x['junction_read_context'].split('|')[0][-5:] + x['junction_read_context'].split('|')[2][:5],
      g(c, e0 - 4, e0) + g(c, s0, s0 + 4))

# T2 rolling circle: two crossings, overlapping passes
r = make_read(2, [piece(0, 499, c, '+', [(s0 + 600, e0)]), piece(499, 1598, c, '+', [E7]), piece(1598, 1998, c, '+', [(s0, s0 + 399)])], primary=1)
rows = rows_of(r)
check('T2 rows', len(rows), 2)
check('T2 joints', [(x['joint_type'], x['donor'], x['acceptor'], x['crossings']) for x in rows], [('back', e0, s0, 2)] * 2)
check('T2 pair overlaps', [x['pair_overlap_nt'] for x in rows], [499, 400])
check('T2 outside', [x['outside_circle_nt'] for x in rows], [0, 0])

# T3 overlap resolution: piece 1 runs 3 nt past the donor; the shift with 2 registered ends wins
r = make_read(3, [piece(0, 502, c, '+', [(s0 + 600, e0 + 3)]), piece(499, 899, c, '+', [(s0, s0 + 399)])],
              read=g(c, s0 + 600, e0) + g(c, s0, s0 + 399))
x = rows_of(r)[0]
check('T3 overlap/shift/donor/acceptor', (x['overlap_nt'], x['shift_nt'], x['donor'], x['acceptor']), (3, 3, e0, s0))
check('T3 trimmed a_exons', x['a_exons'], f'{s0 + 600}-{e0}')

# T4 linear molecule with a back joint (E7 tail > E7 > E8.1 > E9.1): not contained
E8, E9 = (33328510, 33328600), (33335658, 33335766)
r = make_read(4, [piece(0, 399, c, '+', [(e0 - 398, e0)]), piece(399, 399 + 1099 + 91 + 109, c, '+', [E7, E8, E9])], primary=1)
x = rows_of(r)[0]
check('T4 joint', (x['joint_type'], x['donor'], x['acceptor']), ('back', e0, s0))
check('T4 outside circle', x['outside_circle_nt'], 91 + 109)

# T5 forward joint in one gene: E2.1 > E14.1
E2, E14 = (33257587, 33257889), (33341563, 33341686)
x = rows_of(make_read(5, [piece(0, 303, c, '+', [E2]), piece(303, 427, c, '+', [E14])]))[0]
check('T5 forward', (x['joint_type'], x['joint_class'], x['genomic_distance']), ('forward', 'forward_same_gene', E14[0] - E2[1] - 1))
check('T5 circle fields empty', (x['circle_span_nt'], x['crossings']), ('', ''))

# T6 inverted overlapping (hairpin-like): E7 first 500 nt then the same bases reversed
r = make_read(6, [piece(0, 500, c, '+', [(s0, s0 + 499)]), piece(500, 1000, c, '-', [(s0, s0 + 499)])])
x = rows_of(r)[0]
check('T6 inverted', (x['joint_type'], x['joint_class']), ('inverted', 'inverted_overlapping'))

# T7 lariat-like (+): branch-point region 30 nt before the I7.1 end joined to the I7.1 start
i_s, i_e = 33309058, 33328509
r = make_read(7, [piece(0, 200, c, '+', [(i_e - 30 - 199, i_e - 30)]), piece(200, 500, c, '+', [(i_s, i_s + 299)])])
x = rows_of(r)[0]
check('T7 lariat', (x['joint_type'], x['joint_class'], x['lariat_like'], x['donor_to_intron_end_nt'], x['acceptor']),
      ('back', 'lariat_like', True, 30, i_s))

# ---------------------------------------------------------------- RUNX1 E10.1 (-), chr21:35932065-35932477
c2, a0, b0 = 'chr21', 35932065, 35932477
T = tsense(c2, '-', [(a0, b0)])
k, m = 150, 200
r = make_read(8, [piece(0, len(T) - k, c2, '-', [(a0, b0 - k)]), piece(len(T) - k, len(T) - k + m, c2, '-', [(b0 - m + 1, b0)])])
x = rows_of(r)[0]
check('T8 minus back-splice', (x['joint_type'], x['donor'], x['acceptor']), ('back', a0, b0))
check('T8 annotations', (x['donor_annotation'], x['acceptor_annotation']), ('ENSG00000159216:E10.1', 'ENSG00000159216:E10.1'))
dn2, an2 = rc(g(c2, a0 - 2, a0 - 1)), rc(g(c2, b0 + 1, b0 + 2))
check('T8 dinucleotides = genome', (x['donor_dinucleotide'], x['acceptor_dinucleotide']), (dn2, an2))
check('T8 contained', (x['circle_start'], x['circle_end'], x['outside_circle_nt'], x['crossings']), (a0, b0, 0, 1))
check('T8 read context', x['junction_read_context'].replace('|', ''), T[-20:] + T[:20])

# T9 minus lariat: I10.1 transcript start 35932064, end 35899191; donor 30 nt before the end
ts, te = 35932064, 35899191
r = make_read(9, [piece(0, 200, c2, '-', [(te + 30, te + 30 + 199)]), piece(200, 500, c2, '-', [(ts - 299, ts)])])
x = rows_of(r)[0]
check('T9 minus lariat', (x['joint_class'], x['donor'], x['acceptor'], x['donor_to_intron_end_nt']), ('lariat_like', te + 30, ts, 30))

# T10 mixed molecule: RUNX1 piece on chr21 then a HIPK3 circle crossing on chr11
r = make_read(10, [piece(0, 200, c2, '-', [(a0, a0 + 199)]), piece(200, 699, c, '+', [(s0 + 600, e0)]),
                   piece(699, 1099, c, '+', [(s0, s0 + 399)])], primary=1)
rows = rows_of(r)
check('T10 rows (cross pair skipped)', [(x['pair_index'], x['joint_type']) for x in rows], [(2, 'back')])
check('T10 n_chromosomes, outside', (rows[0]['n_chromosomes'], rows[0]['outside_circle_nt']), (2, 200))

# T11 cross-chromosome only molecule: no rows
check('T11 cross only', rows_of(make_read(11, [piece(0, 200, c2, '-', [(a0, a0 + 199)]), piece(200, 699, c, '+', [(s0 + 600, e0)])])), [])

# ---------------------------------------------------------------- circHIPK3: HIPK3 exon 2 of MANE ENST00000303296
# chr11:33286413-33287511 (1,099 nt) = Ensembl 91 regions E6.1 + E6.2 (reference/transcripts.tsv.gz)
h0, h1 = 33286413, 33287511
check('T13 circHIPK3 exon length', h1 - h0 + 1, 1099)
x = rows_of(make_read(13, [piece(0, 599, c, '+', [(h0 + 500, h1)]), piece(599, 899, c, '+', [(h0, h0 + 299)])]))[0]
check('T13 circHIPK3 joint', (x['donor'], x['acceptor'], x['donor_annotation'], x['acceptor_annotation']),
      (h1, h0, 'ENSG00000110422:E6.2', 'ENSG00000110422:E6.1'))
check('T13 circHIPK3 dinucleotides/class', (x['donor_dinucleotide'], x['acceptor_dinucleotide'], x['joint_class'], x['slide_nt']),
      ('GT', 'AG', 'backsplice_registered_canonical', 0))

# T14 breakpoint slide: real read molecule/16311915 (WF40 source001) geometry. Exon 2 tail + exon 3 piece
# ends at 33328630 and the next piece starts at 33286410; the 3' end of exon 3 (33328631-33328633)
# equals the 3 nt before exon 2 (33286410-33286412), so the joint slides to 33328633 > 33286413.
same3 = g(c, 33328631, 33328633) == g(c, 33286410, 33286412)
check('T14 precondition: 3 nt shared', same3, True)
right = 0
while g(c, 33328631 + right, 33328631 + right) == g(c, 33286410 + right, 33286410 + right) and right < 20:
    right += 1
pa = piece(0, 760, c, '+', [(33286873, 33287511), (33328510, 33328630)])
pb = piece(760, 760 + 647, c, '+', [(33286410, 33286410 + 646)])
x = rows_of(make_read(14, [pa, pb]))[0]
check('T14 slide', (x['aligned_donor'], x['aligned_acceptor'], x['slide_right_nt'], x['slide_nt'], x['donor'], x['acceptor']),
      (33328630, 33286410, right, 3, 33328633, 33286413))
check('T14 registered after slide', (x['donor_annotation'], x['acceptor_annotation'], x['joint_class']),
      ('ENSG00000110422:E8.2', 'ENSG00000110422:E6.1', 'backsplice_registered_canonical'))
check('T14 trimmed exons follow the slide', (x['a_exons'].split(',')[-1], x['b_exons']), ('33328510-33328633', f'33286413-{33286410 + 646}'))
seq14 = tsense(c, '+', pa['blocks']) + tsense(c, '+', pb['blocks'])
check('T14 read context at moved joint', x['junction_read_context'].replace('|', ''), seq14[763 - 20:763 + 20])

# T15 no slide when the read has inserted bases between the pieces
pa15 = piece(0, 760, c, '+', [(33286873, 33287511), (33328510, 33328630)])
pb15 = piece(765, 765 + 647, c, '+', [(33286410, 33286410 + 646)])
read15 = tsense(c, '+', pa15['blocks']) + 'ACGTA' + tsense(c, '+', pb15['blocks'])
x = rows_of(make_read(15, [pa15, pb15], read=read15))[0]
check('T15 gap blocks slide', (x['query_gap'], x['slide_nt'], x['donor'], x['acceptor']), (5, 0, 33328630, 33286410))

# T16 non-colinear joint between two genes (CAPRIN1 E3.3 end 34052636 > HIPK3 exon 2 start 33286413): not a circle
x = rows_of(make_read(16, [piece(0, 200, c, '+', [(34052636 - 199, 34052636)]), piece(200, 500, c, '+', [(h0, h0 + 299)])]))[0]
canon16 = (g(c, 34052637, 34052638), g(c, h0 - 2, h0 - 1)) in {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}
check('T16 different genes', (x['joint_type'], x['same_gene'], x['joint_class']),
      ('back', False, 'backjoint_different_genes' + ('_registered_canonical' if canon16 else '')))

# T17 joint at a contig start (cluster run failure 2026-11-05, ValueError start out of range): undefined dinucleotide
x = rows_of(make_read(17, [piece(0, 200, c2, '+', [(1000, 1199)]), piece(200, 500, c2, '+', [(1, 300)])]))[0]
check('T17 contig start', (x['joint_type'], x['acceptor'], x['acceptor_dinucleotide'], x['canonical_dinucleotides']),
      ('back', 1, '', False))

# T18 read across the chrM origin (16569 > 1): mitochondrial class, both dinucleotides undefined
mt = fasta.get_reference_length('MT')
x = rows_of(make_read(18, [piece(0, 170, 'chrM', '+', [(mt - 169, mt)]), piece(170, 370, 'chrM', '+', [(1, 200)])]))[0]
check('T18 chrM origin', (x['joint_class'], x['donor'], x['acceptor'], x['donor_dinucleotide'], x['acceptor_dinucleotide']),
      ('mitochondrial_back', mt, 1, '', ''))

# T19 circTET2 (hsa_circ_0070562, GRCh38 chr4:105233897-105237351 +): Ensembl 91 intron region TET2 I6.2 starts at the
# exon E6.1 acceptor; a back-splice into a registered exon start is a back-splice, not a lariat (PoC 2026-11-05 mislabel)
t0, t1 = 105233897, 105237351
x = rows_of(make_read(19, [piece(0, 500, 'chr4', '+', [(t1 - 499, t1)]), piece(500, 900, 'chr4', '+', [(t0, t0 + 399)])]))[0]
check('T19 circTET2 not lariat', (x['joint_type'], x['lariat_like'], x['donor_annotation'], x['acceptor_annotation'], x['joint_class']),
      ('back', False, 'ENSG00000168769:E6.2', 'ENSG00000168769:E6.1', 'backsplice_registered_canonical'))

# T12 genome helpers = scripts/report_stream.py functions
src = (Path(__file__).parent / 'fixtures/report_stream_baseline.py').read_text()
tree = ast.parse(src)
ns = dict(fasta=fasta, names=set(fasta.references), comp=str.maketrans('ACGTNacgtn', 'TGCANtgcan'))
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in ('dinucleotides', 'oriented', 'junction_homology')]
                        + [n for n in tree.body if isinstance(n, ast.Assign) and any(getattr(t, 'id', '') == 'fa_name' for t in n.targets)],
                        type_ignores=[]), 'report_stream_subset', 'exec'), ns)
rng = random.Random(0)
cases = [('chr11', '+', e0, '+', s0), ('chr21', '-', a0, '-', b0), ('chr11', '+', 33328479, '+', 33309058)]
cases += [(ch, st1, rng.randint(1_000_000, 40_000_000), st2, rng.randint(1_000_000, 40_000_000))
          for ch in ('chr11', 'chr21', 'chrX') for st1 in '+-' for st2 in '+-' for _ in range(25)]
mismatch = 0
for ch, st1, d, st2, a in cases:
    mine = (genome.dinucleotides(ch, st1, d, 'donor'), genome.dinucleotides(ch, st2, a, 'acceptor'),
            genome.junction_homology(ch, st1, d, ch, st2, a))
    ref = (ns['dinucleotides'](ch, st1, d, 'donor'), ns['dinucleotides'](ch, st2, a, 'acceptor'),
           ns['junction_homology'](ch, st1, d, ch, st2, a))
    mismatch += mine != ref
check(f'T12 genome helpers equal report_stream on {len(cases)} joints (mismatches)', mismatch, 0)

summary = dict(module=str(Path(N.__file__).resolve()), genome=GENOME, exon_model=str(REF), tests=len(results),
               passed=sum(r['passed'] for r in results), failed=failures, results=results)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(summary, indent=1, default=str) + '\n')
print(f"{summary['passed']} of {summary['tests']} pass; failed: {failures}")
for r in results:
    if not r['passed']:
        print('FAIL', r)
sys.exit(1 if failures else 0)
