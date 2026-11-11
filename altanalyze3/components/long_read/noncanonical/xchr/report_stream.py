"""Streaming report (one job; same logic as report.py plus junction_homology_nt): cross-chromosome joints from the stage-1 read-level tables, with the validity
criteria, cell states, cohorts, the GENCODE reading frame across the joint and, when stages 3-4
exist, the reused collapse catalog and its translation. Writes <root>/results/ only.

Criteria per joint (columns, never silent filters; every joint is reported):
  1 adjacent      >= 1 molecule whose two alignments meet on the read (query_gap <= 0)
  2 registered    donor = registered exon-region end of gene A (not its last region) and acceptor =
                  registered exon-region start of gene B (not its first region), Ensembl 91 model
  3 coding        both genes GENCODE v45 protein_coding and a GENCODE CDS->CDS joint in frame
  valid_coding = 1 and 2 and 3.
"""
import collections
import csv
import gzip
import json
import re
import statistics
import sys
import time
from pathlib import Path
from altanalyze3.components.long_read.noncanonical.xchr.common import config_path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from altanalyze3.components.long_read.noncanonical.xchr.common import X, REVIO, KINNEX, REF_GFF, GENOME, SAMPLES, BARCODE_CLUSTERS, inside_x, refuse_existing

csv.field_size_limit(sys.maxsize)
root = inside_x(Path(sys.argv[1])) if len(sys.argv) > 1 else X
out = root / (sys.argv[2] if len(sys.argv) > 2 else 'results')
refuse_existing(out / 'REPORT.json')
out.mkdir(exist_ok=True)
SYMBOLS = config_path('gene_symbols')
HGNC = config_path('hgnc')
BCR, ABL1 = 'ENSG00000186716', 'ENSG00000097007'
samples = {r['library']: r for r in csv.DictReader(open(SAMPLES), delimiter='\t')}
COHORTS = ['young', 'aged', 'MDS', 'AML']
report = dict(run_root=str(root), started=time.strftime('%Y-%m-%dT%H:%M:%S'))

# ---- per-library extraction stats; the read-level tables are streamed later, one library at a time ----
stats = {}
for lib in samples:
    proof = root / 'samples' / lib / 'STAGE1_COMPLETE.json'
    if proof.exists():
        stats[lib] = json.loads(proof.read_text())['extract_stats']
report['libraries_with_stage1'] = len(stats)
table = lambda lib: root / 'samples' / lib / f'{lib}.cross_chromosome_junctions.tsv.gz'


def cross_rows(lib):
    with gzip.open(table(lib), 'rt') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            if r['cross_chromosome'] == 'True':
                yield r


def full_record(r):
    # Joints get a full record (annotation, frame, homology) when either side is a registered exon
    # boundary or BCR/ABL1 is a partner. All other joints are counted by cohort in REPORT.json; their
    # rows stay in samples/<library>/<library>.cross_chromosome_junctions.tsv.gz.
    return bool(r['donor_annotation'] or r['acceptor_annotation']) or bool({r['gene_a'], r['gene_b']} & {BCR, ABL1})


genes = {BCR, ABL1}
for lib in stats:
    for r in cross_rows(lib):
        if full_record(r):
            genes.update((r['gene_a'], r['gene_b']))

clusters = collections.defaultdict(dict)
with open(BARCODE_CLUSTERS) as h:
    for line in h:
        key, cluster = line.rstrip('\n').split('\t')[:2]
        barcode, lib = key.split('.', 1)
        clusters[lib].setdefault(barcode, cluster)
rc = lambda b: b.split('-')[0].translate(str.maketrans('ACGTN', 'TGCAN'))[::-1] + '-' + b.split('-')[1]

# ---- gene symbols, HGNC groups, GENCODE v45 gene types and coding transcripts ----
symbols = {}
with open(SYMBOLS) as h:
    for line in h:
        p = line.rstrip('\n').split('\t')
        if len(p) > 1:
            symbols[p[0]] = p[1]
hgnc = {}
with open(HGNC) as h:
    for r in csv.DictReader(h, delimiter='\t'):
        if r['status'] == 'Approved':
            for g in filter(None, r['ensembl_gene_id'].split('|')):
                hgnc.setdefault(g, set()).update(filter(None, r['gene_group_id'].split('|')))
gene_type, tx = {}, collections.defaultdict(lambda: dict(exons=[], cds=[], tags='', gene='', chrom='', strand=''))
attr = lambda s, k: (re.search(rf'(?:^|;){k}=([^;]+)', s) or [None, ''])[1]
with open(REF_GFF) as h:
    for line in h:
        if line[0] == '#' or 'gene_id=ENSG' not in line:
            continue
        c = line.rstrip('\n').split('\t')
        if c[2] not in ('gene', 'exon', 'CDS'):
            continue
        g = attr(c[8], 'gene_id').split('.')[0]
        if g not in genes:
            continue
        if c[2] == 'gene':
            gene_type[g] = attr(c[8], 'gene_type')
            continue
        t = tx[attr(c[8], 'transcript_id').split('.')[0]]
        t.update(gene=g, chrom=c[0], strand=c[6], tags=attr(c[8], 'tag'))
        if c[2] == 'exon':
            t['exons'].append((int(c[3]), int(c[4])))
        else:
            t['cds'].append((int(c[3]), int(c[4]), int(c[7])))
by_gene = collections.defaultdict(list)
for tid, t in tx.items():
    by_gene[t['gene']].append((tid, t))


def side(gene, chrom, strand, pos, kind):
    """Per GENCODE transcript of the gene with an exon boundary at pos: ('CDS', codon position of pos,
    CDS nt before/after, tid) or (UTR/noncoding class, tid)."""
    hits = []
    for tid, t in by_gene.get(gene, []):
        if t['chrom'] != chrom or t['strand'] != strand:
            continue
        ends = {e if strand == '+' else s for s, e in t['exons']}
        starts = {s if strand == '+' else e for s, e in t['exons']}
        if pos not in (ends if kind == 'donor' else starts):
            continue
        if not t['cds']:
            hits.append(('noncoding_transcript', None, None, tid))
            continue
        cds = sorted(t['cds'], key=lambda x: x[0], reverse=strand == '-')
        total = sum(e - s + 1 for s, e, _ in cds)
        before = 0
        found = None
        for s, e, phase in cds:
            length = e - s + 1
            p5, p3 = (s, e) if strand == '+' else (e, s)
            if kind == 'donor' and pos == p3:
                found = ('CDS', (length - 1 - phase) % 3, before + length, tid)
            if kind == 'acceptor' and pos == p5:
                found = ('CDS', (-phase) % 3, total - before, tid)
            before += length
        if found:
            hits.append(found)
            continue
        lo, hi = min(s for s, _, _ in cds), max(e for _, e, _ in cds)
        upstream = pos < lo if strand == '+' else pos > hi
        hits.append(('5UTR' if upstream else '3UTR', None, None, tid))
    return hits


def frame(key):
    ga, ca, sa, donor, gb, cb, sb, acceptor = key
    ca, cb = (c if c.startswith('chr') else 'chr' + c for c in (ca, cb))
    d, a = side(ga, ca, sa, donor, 'donor'), side(gb, cb, sb, acceptor, 'acceptor')
    dc = [h for h in d if h[0] == 'CDS']; ac = [h for h in a if h[0] == 'CDS']
    pref = lambda h: ('MANE_Select' in tx[h[3]]['tags'], 'Ensembl_canonical' in tx[h[3]]['tags'])
    inframe = [(x, y) for x in dc for y in ac if y[1] == (x[1] + 1) % 3]
    if inframe:
        x, y = max(inframe, key=lambda p: (pref(p[0]), pref(p[1])))
        return dict(frame_class='CDS_to_CDS_in_frame', tx_a=x[3], tx_b=y[3], aa_from_a=x[2] // 3,
                    aa_from_b=y[2] // 3, expected_fusion_aa=(x[2] + y[2]) // 3)
    if dc and ac:
        return dict(frame_class='CDS_to_CDS_out_of_frame', tx_a=dc[0][3], tx_b=ac[0][3])
    if any(h[0] == '5UTR' for h in a):
        return dict(frame_class='A_to_B_5UTR_B_start_codon_retained', tx_b=next(h[3] for h in a if h[0] == '5UTR'))
    if ac and any(h[0] == '5UTR' for h in d):
        return dict(frame_class='A_5UTR_to_B_CDS', tx_b=ac[0][3])
    if dc:
        return dict(frame_class='A_CDS_to_B_UTR_or_noncoding', tx_a=dc[0][3])
    if not d or not a:
        return dict(frame_class='no_GENCODE_transcript_boundary_at_joint')
    return dict(frame_class='noncoding_or_UTR_both_sides')


fasta = pysam.FastaFile(str(GENOME))
names = set(fasta.references)
fa_name = lambda c: (c[3:] if c.startswith('chr') else c) if (c[3:] if c.startswith('chr') else c) in names else ('MT' if c in ('chrM', 'M') else c)
comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def dinucleotides(chrom, strand, pos, kind):
    c = fa_name(chrom)
    if c not in names:
        return ''
    if kind == 'donor':
        s = fasta.fetch(c, pos, pos + 2) if strand == '+' else fasta.fetch(c, pos - 3, pos - 1).translate(comp)[::-1]
    else:
        s = fasta.fetch(c, pos - 3, pos - 1) if strand == '+' else fasta.fetch(c, pos, pos + 2).translate(comp)[::-1]
    return s.upper()


def oriented(chrom, start0, end0, strand):
    c = fa_name(chrom)
    if c not in names:
        return ''
    piece = fasta.fetch(c, max(0, start0), max(0, end0)).upper()
    return piece if strand == '+' else piece.translate(comp)[::-1]


def junction_homology(ca, sa, d, cb, sb, a, span=20):
    """Bases the breakpoint can slide over sequence shared by both partners (left + right of the
    joint, at most span each). RT template switching leaves shared sequence at the joint."""
    a_exon_end = oriented(ca, d - span, d, '+') if sa == '+' else oriented(ca, d - 1, d - 1 + span, '-')
    a_intron = oriented(ca, d, d + span, '+') if sa == '+' else oriented(ca, d - 1 - span, d - 1, '-')
    b_intron = oriented(cb, a - 1 - span, a - 1, '+') if sb == '+' else oriented(cb, a, a + span, '-')
    b_exon = oriented(cb, a - 1, a - 1 + span, '+') if sb == '+' else oriented(cb, a - span, a, '-')
    left = 0
    while left < min(len(a_exon_end), len(b_intron)) and a_exon_end[-1 - left] == b_intron[-1 - left]:
        left += 1
    right = 0
    while right < min(len(a_intron), len(b_exon)) and a_intron[right] == b_exon[right]:
        right += 1
    return left + right


# ---- optional: reused collapse catalog + translation for kept molecules ----
protein_of = {}
gff_out = root / 'gff-output'
if (root / 'STAGE3_COMPLETE.json').exists():
    def any_open(path):
        gz = Path(str(path) + '.gz')
        return gzip.open(gz, 'rt') if gz.exists() else (open(path) if Path(path).exists() else None)
    report['stage3'] = json.loads((root / 'STAGE3_COMPLETE.json').read_text())
    s2e = {}
    with any_open(gff_out / 'FINAL_structure_to_exemplar.tsv') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            s2e[(r['gene'], r['structure'])] = r['final_isoform_id']
    proteins = {}
    h = any_open(gff_out / 'protein_summary.txt')
    if h is not None:
        with h:
            for r in csv.DictReader(h, delimiter='\t'):
                proteins[r['Transcript ID']] = (r['Protein Length'], r['NMD Status'])
    for lib in stats:
        h = any_open(root / 'samples' / lib / 'gff-output/transcript_associations.txt')
        if h is None:
            continue
        with h:
            for line in h:
                p = line.rstrip('\n').split('\t')
                ex = s2e.get((p[0], p[2]))
                if ex:
                    bare = ex.split('.', 1)[0]
                    protein_of[f'{p[0]}:{p[3]}'] = (ex,) + proteins.get(bare, ('', ''))
    report['catalog_exemplars'] = len(set(s2e.values()))
    report['kept_molecules_with_catalog_exemplar'] = len(protein_of)

# ---- aggregate per joint, streaming one library at a time ----
J = {}
unk = lambda g: 'UNK' if g.startswith('UNK') else g  # UNK<n> numbers are per-process counters
orientation, match_rates, bcr_abl1_rows, n_rows = {}, {}, [], 0
counted_only = collections.Counter()
for lib in stats:
    bcs = {r['barcode'] for r in cross_rows(lib) if r['barcode']}
    fwd = sum(b in clusters[lib] for b in bcs)
    rev = sum(rc(b) in clusters[lib] for b in bcs)
    orientation[lib] = 'reverse_complement' if rev > fwd else 'forward'
    match_rates[lib] = dict(barcodes=len(bcs), forward=fwd, reverse_complement=rev)
    s_ = samples[lib]
    for r in cross_rows(lib):
        n_rows += 1
        b = r['barcode']
        r['cell_state'] = (clusters[lib].get(rc(b) if orientation[lib] == 'reverse_complement' else b, '')) if b else ''
        r['gene_a'], r['gene_b'] = unk(r['gene_a']), unk(r['gene_b'])
        if {r['gene_a'], r['gene_b']} == {BCR, ABL1}:
            bcr_abl1_rows.append(r)
        if not full_record(r):
            counted_only['rows'] += 1
            counted_only[s_['cohort'] + '_rows'] += 1
            counted_only['status:' + r['status']] += 1
            continue
        key = (r['gene_a'], r['a_chrom'], r['a_strand'], int(r['donor']), r['gene_b'], r['b_chrom'], r['b_strand'], int(r['acceptor']))
        j = J.get(key)
        if j is None:
            j = J[key] = dict(mol=0, kept=0, adjacent=0, libs=set(), donors=set(), cohort_mol=collections.Counter(),
                              cohort_donors=collections.defaultdict(set), states=collections.Counter(), gaps=[], mapq=[],
                              overlap=0, first={k: r[k] for k in ('donor_annotation', 'acceptor_annotation',
                                                                  'donor_is_gene_last_region', 'acceptor_is_gene_first_region')},
                              example='', exemplars=collections.Counter())
        j['mol'] += 1; j['libs'].add(lib); j['donors'].add(s_['original_uid'])
        j['cohort_mol'][s_['cohort']] += 1; j['cohort_donors'][s_['cohort']].add(s_['original_uid'])
        j['gaps'].append(int(r['query_gap'])); j['mapq'].append(min(int(r['a_mapq']), int(r['b_mapq'])))
        j['overlap'] = max(j['overlap'], int(r['overlap_nt']))
        if int(r['query_gap']) <= 0:
            j['adjacent'] += 1
        if r['status'] == 'kept':
            j['kept'] += 1
            if not j['example']:
                j['example'] = r['junction_read_context']
            if r['isoform_id'] in protein_of:
                j['exemplars'][protein_of[r['isoform_id']]] += 1
        if r['cell_state']:
            j['states'][r['cell_state']] += 1
report['cross_chromosome_joint_rows'] = n_rows
report['joint_rows_without_registered_boundary_counted_only'] = dict(counted_only)
report['barcode_match'] = match_rates
out_rows = []
for key, j in J.items():
    ga, ca, sa, donor, gb, cb, sb, acceptor = key
    f = j['first']
    registered = bool(f['donor_annotation']) and bool(f['acceptor_annotation']) and \
        f['donor_is_gene_last_region'] != 'True' and f['acceptor_is_gene_first_region'] != 'True'
    fr = frame(key)
    coding = gene_type.get(ga) == 'protein_coding' and gene_type.get(gb) == 'protein_coding' and fr['frame_class'] == 'CDS_to_CDS_in_frame'
    row = dict(gene_a=ga, symbol_a=symbols.get(ga, ''), gene_b=gb, symbol_b=symbols.get(gb, ''),
               joint=f'{ca}:{donor}:{sa} > {cb}:{acceptor}:{sb}', donors=len(j['donors']), molecules=j['mol'],
               kept_molecules=j['kept'], adjacent_molecules=j['adjacent'], libraries=len(j['libs']),
               **{f'{c}_donors': len(j['cohort_donors'][c]) for c in COHORTS},
               **{f'{c}_molecules': j['cohort_mol'][c] for c in COHORTS},
               donor_ids=','.join(sorted(j['donors'])),
               cell_states='; '.join(f'{k}:{v}' for k, v in j['states'].most_common(6)),
               donor_annotation=f['donor_annotation'], acceptor_annotation=f['acceptor_annotation'],
               donor_is_gene_last_region=f['donor_is_gene_last_region'], acceptor_is_gene_first_region=f['acceptor_is_gene_first_region'],
               donor_dinucleotide=dinucleotides(ca, sa, donor, 'donor'), acceptor_dinucleotide=dinucleotides(cb, sb, acceptor, 'acceptor'),
               gene_type_a=gene_type.get(ga, ''), gene_type_b=gene_type.get(gb, ''),
               shared_HGNC_group=bool(hgnc.get(ga, set()) & hgnc.get(gb, set())),
               junction_homology_nt=junction_homology(ca, sa, donor, cb, sb, acceptor),
               max_read_overlap_nt=j['overlap'], query_gap_min=min(j['gaps']), query_gap_median=statistics.median(j['gaps']),
               min_mapq_median=statistics.median(j['mapq']), **fr,
               catalog_exemplars='; '.join(f'{e[0]} ({e[1]} aa, {e[2]}) x{n}' for e, n in j['exemplars'].most_common(3)),
               example_read_context=j['example'],
               criterion1_adjacent=j['adjacent'] > 0, criterion2_registered=registered, criterion3_coding=coding)
    row['valid_coding'] = row['criterion1_adjacent'] and registered and coding
    out_rows.append(row)
out_rows.sort(key=lambda r: (-r['valid_coding'], -r['donors'], -r['molecules'], r['joint']))
fields = list(out_rows[0]) if out_rows else []
for name, subset in (('cross_chromosome_joints_all.tsv', out_rows),
                     ('cross_chromosome_joints_valid_coding.tsv', [r for r in out_rows if r['valid_coding']]),
                     ('BCR_ABL1_joints.tsv', [r for r in out_rows if {r['gene_a'], r['gene_b']} == {BCR, ABL1}]),
                     ('BCR_or_ABL1_any_partner_joints.tsv', [r for r in out_rows if {r['gene_a'], r['gene_b']} & {BCR, ABL1}])):
    with (gzip.open(out / (name + '.gz'), 'wt') if name.endswith('_all.tsv') else open(out / name, 'w')) as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(subset)
    report[name] = len(subset)
with open(out / 'BCR_ABL1_molecules.tsv', 'w') as h:
    sub = bcr_abl1_rows
    if sub:
        w = csv.DictWriter(h, fieldnames=list(sub[0]), delimiter='\t')
        w.writeheader()
        w.writerows(sub)
    report['BCR_ABL1_molecule_rows'] = len(sub)
coh = {}
for c in COHORTS:
    libs = [l for l in stats if samples[l]['cohort'] == c]
    primary = sum(stats[l].get('primary_records', 0) for l in libs)
    xreads = sum(stats[l].get('cross_chromosome_reads', 0) for l in libs)
    kept = sum(stats[l].get('kept_molecules', 0) for l in libs)
    valid = [r for r in out_rows if r['valid_coding'] and r[f'{c}_molecules']]
    coh[c] = dict(libraries=len(libs), donors=len({samples[l]['original_uid'] for l in libs}), primary_records=primary,
                  cross_chromosome_primary_reads=xreads, cross_chromosome_per_million_primary=round(1e6 * xreads / primary, 1) if primary else None,
                  kept_molecules=kept, valid_coding_joints=len(valid),
                  valid_coding_molecules=sum(r[f'{c}_molecules'] for r in valid))
report['cohorts'] = coh
status = collections.Counter()
for s in stats.values():
    status.update({k: v for k, v in s.items() if k.startswith('status:')})
report['extraction_status_counts'] = dict(status)
report['joints'] = len(out_rows)
report['criteria_counts'] = {k: sum(bool(r[k]) for r in out_rows) for k in
                             ('criterion1_adjacent', 'criterion2_registered', 'criterion3_coding', 'valid_coding')}
report['frame_class_counts'] = dict(collections.Counter(r['frame_class'] for r in out_rows))
report['finished'] = time.strftime('%Y-%m-%dT%H:%M:%S')
(out / 'REPORT.json').write_text(json.dumps(report, indent=1, default=str) + '\n')
print(json.dumps({k: v for k, v in report.items() if k != 'barcode_match'}, indent=1, default=str))
