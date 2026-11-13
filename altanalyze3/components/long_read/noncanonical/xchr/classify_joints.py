"""Classify cross-chromosome joints (EXPERIMENTAL extra-chr-trans-splicing).

Classes (first match wins; rules in forked_altanalyze3/.../bam/cross_chromosome_rules.py, tested in
scripts/test_xchr_rules.py): artifact_blacklist, artifact_paralog_sequence, artifact_template_switch,
gene_fusion_candidate, artifact_not_splice_precise, trans_splicing_candidate, unresolved.

Evidence per joint (only this study's reads, the genome, GENCODE v45 / Ensembl 91 annotation, and the
curated cancer fusion collections of the CTAT human fusion library for database membership):
  template switching   genome and mRNA-context homology vs the reference null (99th percentile),
                       polyA/T run at the joint, splice precision, position specificity (Poisson test of
                       the molecules at this exact joint against uniform breakpoints along both mRNAs,
                       using every molecule between the same two genes)
  paralog sequence     the representative segment of one partner aligns (mappy, identity >= 0.8 over
                       >= 50 nt) to an mRNA of the other partner gene
  gene fusion          curated database pair; donor restriction (>= 90% of molecules in one donor) with
                       a reciprocal B->A molecule, an intron-to-intron molecule (both ends inside the
                       gene bodies, outside exons) in that donor, partner-exon imbalance (Fisher, joint-
                       positive cells vs negative cells of the same cell states, same library: fewer
                       reads carrying gene A exons 3' of the donor, or gene B exons 5' of the acceptor),
                       or lower expression of protein-coding genes within 1 Mb of either breakpoint.
                       Expression tests run only for precise, donor-restricted joints in >= 3 cells:
                       cell-level permutation (10,000) within library x cell state x depth quintile,
                       Benjamini-Hochberg q < 0.05 across tested joints. Donor restriction needs >= 2 cells in that donor;
                       a reciprocal molecule must be splice-precise (registered donor and acceptor,
                       alignments meet); intron-to-intron molecules must recur at one position (5-nt bins).
Arguments: run root, worker count, optional results folder name (default results), optional output folder name,
optional library -> patient map (config/patients.tsv). With the map, a donor is a patient: timepoints of one
patient count once for the >= 2 donor rule and pool for donor restriction; 'donors' holds the patient count and
'sample_donors' keeps the original_uid count.
"""
import collections
import csv
import gzip
import json
import math
import random
import sys
import time
from multiprocessing import Pool
from pathlib import Path
from altanalyze3.components.long_read.noncanonical.xchr.common import config_path

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from altanalyze3.components.long_read.noncanonical.xchr.common import X, KINNEX, GENOME, SAMPLES, inside_x, refuse_existing
from altanalyze3.components.bam import cross_chromosome_rules as R

csv.field_size_limit(sys.maxsize)
root = inside_x(Path(sys.argv[1])) if len(sys.argv) > 1 else X
nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 8
RESULTS = sys.argv[3] if len(sys.argv) > 3 else 'results'
OUT = root / (sys.argv[4] if len(sys.argv) > 4 else ('classification' if RESULTS == 'results' else f'classification_{RESULTS}'))
refuse_existing(OUT / 'CLASSIFICATION_SUMMARY.json')
OUT.mkdir(exist_ok=True)
REF = config_path('reference')
ISV = config_path('expression_dir')
P_POSITION, Q_EXPRESSION, DONOR_RESTRICTED, PARALOG_ID, PARALOG_LEN, MIN_POS_CELLS = 1e-3, 0.05, 0.9, 0.8, 50, 3
FLANK_WINDOW, N_PERM, MIN_CLONE_CELLS = 1_000_000, 10000, 2
samples = {r['library']: r for r in csv.DictReader(open(SAMPLES), delimiter='\t')}
started = time.time()
null = json.loads((REF / 'homology_null.json').read_text())
THR_G, THR_M = null['threshold_genomic_p99'], null['threshold_mrna_p99']

# ---------------------------------------------------------------- reference -------------
genes, gene_rows = {}, []
with gzip.open(REF / 'genes.tsv.gz', 'rt') as h:
    for r in csv.DictReader(h, delimiter='\t'):
        r['start'], r['end'] = int(r['start']), int(r['end'])
        genes[r['gene_id']] = r
        gene_rows.append((r['chrom'], r['start'], r['end'], (r['gene_id'], r['strand'])))
gene_index = R.IntervalIndex(gene_rows)
tx_by_gene = collections.defaultdict(list)
with gzip.open(REF / 'transcripts.tsv.gz', 'rt') as h:
    for r in csv.DictReader(h, delimiter='\t'):
        r['exons'] = [tuple(map(int, e.split('-'))) for e in r['exons'].split(',')]
        tx_by_gene[r['gene_id']].append(r)
merged_exons = {}
for g, ts in tx_by_gene.items():
    iv = sorted((min(a, b), max(a, b)) for t in ts for a, b in t['exons'])
    m = []
    for a, b in iv:
        if m and a <= m[-1][1] + 1:
            m[-1][1] = max(m[-1][1], b)
        else:
            m.append([a, b])
    merged_exons[g] = m


def in_exon(gene, pos):
    for a, b in merged_exons.get(gene, ()):
        if a <= pos <= b:
            return True
    return False


fusion_db = {}
with open(REF / 'fusion_db.tsv') as h:
    for r in csv.DictReader(h, delimiter='\t'):
        fusion_db[(r['gene_a'], r['gene_b'])] = r

# ---------------------------------------------------------------- joints ------------------
joints = {}
with gzip.open(root / RESULTS / 'cross_chromosome_joints_all.tsv.gz', 'rt') as h:
    for r in csv.DictReader(h, delimiter='\t'):
        left, right = r['joint'].split(' > ')
        ca, d, sa = left.rsplit(':', 2)
        cb, a, sb = right.rsplit(':', 2)
        key = (r['gene_a'], ca, sa, int(d), r['gene_b'], cb, sb, int(a))
        joints[key] = r
pairset = set()
for (ga, _, _, _, gb, _, _, _) in joints:
    if ga != 'UNK' and gb != 'UNK':
        pairset.add((ga, gb)); pairset.add((gb, ga))


# ---------------------------------------------------------------- scan of read tables -----
def scan(lib):
    unk = lambda g: 'UNK' if g.startswith('UNK') else g
    per_joint = {}
    pair_all = collections.Counter()
    pair_intronic = collections.Counter()
    pair_precise = collections.Counter()
    pair_intronic_sites = collections.defaultdict(collections.Counter)
    path = root / 'samples' / lib / f'{lib}.cross_chromosome_junctions.tsv.gz'
    with gzip.open(path, 'rt') as h:
        head = h.readline().rstrip('\n').split('\t')
        i = {c: k for k, c in enumerate(head)}
        for line in h:
            p = line.rstrip('\n').split('\t')
            if p[i['cross_chromosome']] != 'True':
                continue
            donor, acc = int(p[i['donor']]), int(p[i['acceptor']])
            ca, sa, cb, sb = p[i['a_chrom']], p[i['a_strand']], p[i['b_chrom']], p[i['b_strand']]
            key = (unk(p[i['gene_a']]), ca, sa, donor, unk(p[i['gene_b']]), cb, sb, acc)
            if key in joints:
                j = per_joint.setdefault(key, dict(barcodes=set(), molecules=0, concatemer=0, status=collections.Counter(),
                                                   a_exons=p[i['a_exons']], b_exons=p[i['b_exons']]))
                j['molecules'] += 1
                j['status'][p[i['status']]] += 1
                if p[i['barcode']]:
                    j['barcodes'].add(p[i['barcode']])
                if R.concatemer_signature(p[i['junction_read_context']]):
                    j['concatemer'] += 1
                if p[i['status']] == 'kept':
                    j['a_exons'], j['b_exons'] = p[i['a_exons']], p[i['b_exons']]
            ga_set = [g for g, s in gene_index.overlaps(ca, donor, donor) if s == sa]
            gb_set = [g for g, s in gene_index.overlaps(cb, acc, acc) if s == sb]
            for ga in ga_set:
                for gb in gb_set:
                    if (ga, gb) in pairset:
                        pair_all[(ga, gb)] += 1
                        if not in_exon(ga, donor) and not in_exon(gb, acc):
                            pair_intronic[(ga, gb)] += 1
                            pair_intronic_sites[(ga, gb)][(donor // 5, acc // 5)] += 1
                        if (p[i['donor_annotation']] and p[i['acceptor_annotation']] and p[i['donor_is_gene_last_region']] != 'True'
                                and p[i['acceptor_is_gene_first_region']] != 'True' and int(p[i['query_gap']]) <= 0):
                            pair_precise[(ga, gb)] += 1
    recurrent = {k: max(c.values()) for k, c in pair_intronic_sites.items()}
    return lib, per_joint, pair_all, pair_intronic, pair_precise, recurrent


libs = [l for l in samples if (root / 'samples' / l / 'STAGE1_COMPLETE.json').exists()]
with Pool(nproc) as pool:
    scanned = pool.map(scan, libs)
per_joint = collections.defaultdict(lambda: dict(barcodes={}, molecules=0, concatemer=0, status=collections.Counter(),
                                                 lib_molecules=collections.Counter(), a_exons='', b_exons=''))
pair_all_lib, pair_intronic_lib = collections.defaultdict(collections.Counter), collections.defaultdict(collections.Counter)
pair_precise_lib, pair_recurrent_lib = collections.defaultdict(collections.Counter), collections.defaultdict(collections.Counter)
for lib, pj, pa, pi, pp, pr in scanned:
    for k, n in pp.items():
        pair_precise_lib[k][lib] += n
    for k, n in pr.items():
        pair_recurrent_lib[k][lib] = n
    for key, j in pj.items():
        t = per_joint[key]
        t['molecules'] += j['molecules']; t['concatemer'] += j['concatemer']; t['status'].update(j['status'])
        t['lib_molecules'][lib] += j['molecules']
        if j['barcodes']:
            t['barcodes'][lib] = j['barcodes']
        if not t['a_exons'] or j['status'].get('kept'):
            t['a_exons'], t['b_exons'] = j['a_exons'], j['b_exons']
    for k, n in pa.items():
        pair_all_lib[k][lib] += n
    for k, n in pi.items():
        pair_intronic_lib[k][lib] += n
PATIENTS = Path(sys.argv[5]) if len(sys.argv) > 5 else None
if PATIENTS:
    patient_of = {r['library']: r['patient'] for r in csv.DictReader(open(PATIENTS), delimiter='\t')}
    missing = [l for l in libs if l not in patient_of]
    if missing:
        raise SystemExit(f'{PATIENTS} lacks libraries: {missing}')
    donor_of = lambda lib: patient_of[lib]
else:
    donor_of = lambda lib: samples[lib]['original_uid']

# ---------------------------------------------------------------- sequence helpers --------
fasta = pysam.FastaFile(str(GENOME))
fnames = set(fasta.references)
fa = lambda c: (c[3:] if c[3:] in fnames else ('MT' if c == 'chrM' else c)) if c.startswith('chr') else c
mrna_cache = {}


def mrna(t):
    k = t['transcript_id']
    if k not in mrna_cache:
        parts = [fasta.fetch(fa(t['chrom']), a - 1, b).upper() for a, b in t['exons']]
        mrna_cache[k] = ''.join(parts if t['strand'] == '+' else [R.revcomp(x) for x in parts])
    return mrna_cache[k]


def boundary_transcripts(gene, pos, kind):
    """(transcript, mRNA index of the boundary base) for transcripts with an exon boundary at pos."""
    out = []
    for t in tx_by_gene.get(gene, ()):
        idx = 0
        for a, b in t['exons']:
            five, three = (a, b) if t['strand'] == '+' else (b, a)
            if kind == 'donor' and pos == three:
                out.append((t, idx + (b - a)))
            if kind == 'acceptor' and pos == five:
                out.append((t, idx))
            idx += b - a + 1
    return out


import mappy
aligner_cache = {}


def gene_aligner(gene):
    if gene not in aligner_cache:
        seqs = [mrna(t) for t in tx_by_gene.get(gene, ())]
        aligner_cache[gene] = mappy.Aligner(seq=('N' * 50).join(seqs), preset='map-ont') if seqs else None
    return aligner_cache[gene]


def segment_seq(chrom, exons_text):
    if not exons_text or fa(chrom) not in fnames:
        return ''
    ex = [tuple(map(int, e.split('-'))) for e in exons_text.split(',')]
    return ''.join(fasta.fetch(fa(chrom), a - 1, b).upper() for a, b in sorted(ex))


def similar(seq, gene):
    al = gene_aligner(gene) if seq and gene != 'UNK' else None
    if al is None:
        return 0.0, 0
    best = (0.0, 0)
    for hit in al.map(seq):
        ident = hit.mlen / hit.blen if hit.blen else 0.0
        if hit.blen >= PARALOG_LEN and ident > best[0]:
            best = (ident, hit.blen)
    return best


def symbols(gene, report_symbol):
    out = {report_symbol} if report_symbol else set()
    if gene in genes:
        out.add(genes[gene]['gene_name'])
    return out


# ---------------------------------------------------------------- evidence per joint ------
rows = []
for key, r in joints.items():
    ga, ca, sa, d, gb, cb, sb, acc = key
    pj = per_joint.get(key)
    molecules = int(r['molecules'])
    statuses = pj['status'] if pj else collections.Counter()
    blacklisted_all = bool(statuses) and all(s.startswith('blacklist_') for s in statuses)
    # template switching
    a_tx, b_tx = boundary_transcripts(ga, d, 'donor'), boundary_transcripts(gb, acc, 'acceptor')
    mh, complete = 0, False
    for ta, ia in a_tx[:5]:
        for tb, ib in b_tx[:5]:
            h_, c_ = R.mrna_homology(mrna(ta), ia, mrna(tb), ib)
            if h_ > mh or (h_ == mh and c_):
                mh, complete = h_, c_
    gh = int(r['junction_homology_nt'])
    homology_flag = gh > THR_G or mh > THR_M
    concatemer_flag = bool(pj) and pj['concatemer'] / max(1, pj['molecules']) >= 0.5
    precise = r['criterion2_registered'] == 'True' and R.canonical_splice(r['donor_dinucleotide'], r['acceptor_dinucleotide'])
    pair_total = sum(pair_all_lib.get((ga, gb), {}).values())
    la = min((sum(b - a + 1 for a, b in t['exons']) for t, _ in a_tx), default=0)
    lb = min((sum(b - a + 1 for a, b in t['exons']) for t, _ in b_tx), default=0)
    slack = (1 + int(r['max_read_overlap_nt'])) * (1 + gh)
    if la and lb and pair_total:
        expected, p_pos = R.position_specificity(molecules, max(pair_total, molecules), la, lb, slack)
    else:
        expected, p_pos = None, None
    position_specific = molecules >= 2 and p_pos is not None and p_pos < P_POSITION
    ts_status = R.template_switch_status(precise, homology_flag, concatemer_flag, position_specific, int(r['adjacent_molecules']) > 0)
    # paralog sequence
    a_seq, b_seq = (segment_seq(ca, pj['a_exons']), segment_seq(cb, pj['b_exons'])) if pj else ('', '')
    sim_ab, sim_ba = similar(a_seq, gb), similar(b_seq, ga)
    paralog_partner = max(sim_ab[0], sim_ba[0]) >= PARALOG_ID
    # fusion database
    db_hits = [(x, y, fusion_db[(x, y)]) for x in symbols(ga, r['symbol_a']) for y in symbols(gb, r['symbol_b']) if (x, y) in fusion_db]
    db_rev = [(y, x, fusion_db[(y, x)]) for x in symbols(ga, r['symbol_a']) for y in symbols(gb, r['symbol_b']) if (y, x) in fusion_db]
    curated = any(e[2]['curated_cancer'] == 'True' for e in db_hits + db_rev)
    tags = sorted({t for e in db_hits + db_rev for t in e[2]['tags'].split(',')})
    # donor restriction, reciprocal and intronic molecules
    donor_mol = collections.Counter()
    for lib, n in (pj['lib_molecules'] if pj else {}).items():
        donor_mol[donor_of(lib)] += n
    top_donor, top_n = donor_mol.most_common(1)[0] if donor_mol else ('', 0)
    pos_cells = {lib: bcs for lib, bcs in (pj['barcodes'] if pj else {}).items()}
    top_donor_cells = len(set().union(*[b for lib, b in pos_cells.items() if donor_of(lib) == top_donor] or [set()]))
    donor_restricted = bool(donor_mol) and top_n / sum(donor_mol.values()) >= DONOR_RESTRICTED and top_donor_cells >= MIN_CLONE_CELLS
    rec = sum(n for lib, n in pair_all_lib.get((gb, ga), {}).items() if donor_of(lib) == top_donor)
    rec_precise = sum(n for lib, n in pair_precise_lib.get((gb, ga), {}).items() if donor_of(lib) == top_donor)
    intr = sum(n for lib, n in pair_intronic_lib.get((ga, gb), {}).items() if donor_of(lib) == top_donor)
    intr_recurrent = max([n for lib, n in pair_recurrent_lib.get((ga, gb), {}).items() if donor_of(lib) == top_donor] or [0])
    rows.append(dict(r, molecules_scanned=pj['molecules'] if pj else 0,
                     status_counts=';'.join(f'{k}:{v}' for k, v in statuses.most_common()),
                     blacklisted_all=blacklisted_all, mrna_homology_nt=mh, mrna_homology_complete=complete,
                     homology_threshold_genomic=THR_G, homology_threshold_mrna=THR_M, homology_flag=homology_flag,
                     concatemer_molecules=pj['concatemer'] if pj else 0, concatemer_flag=concatemer_flag, splice_precise=precise,
                     pair_molecules_all_positions=pair_total, positions_a=la, positions_b=lb, position_slack=slack,
                     position_expected=expected, position_p=p_pos, position_specific=position_specific, ts_status=ts_status,
                     paralog_identity_a_in_b=round(sim_ab[0], 3), paralog_identity_b_in_a=round(sim_ba[0], 3),
                     paralog_partner=paralog_partner, fusion_db_pairs=';'.join(f'{x}--{y}' for x, y, _ in db_hits + db_rev),
                     fusion_db_tags=','.join(tags), fusion_db_curated=curated,
                     fusion_db_rnaseq_or_normal_tags_info_only=','.join(t for t in tags if t in R.RNASEQ_DERIVED_TAGS),
                     top_donor=top_donor, top_donor_share=round(top_n / sum(donor_mol.values()), 3) if donor_mol else None,
                     top_donor_cells=top_donor_cells, donor_restricted=donor_restricted,
                     reciprocal_molecules_top_donor=rec, reciprocal_precise_molecules_top_donor=rec_precise,
                     reciprocal_precise_same_donor=rec_precise > 0,
                     intronic_molecules_top_donor=intr, intronic_max_same_site_top_donor=intr_recurrent,
                     intronic_breakpoint_recurrent_same_donor=intr_recurrent >= 2,
                     positive_cells=sum(len(v) for v in pos_cells.values()), _pos_cells=pos_cells, _key=key))
    if PATIENTS:
        ids = sorted({donor_of(lib) for lib in (pj['lib_molecules'] if pj else {})})
        rows[-1].update(sample_donors=r['donors'], donors=len(ids), patient_ids=','.join(ids))


# ---------------------------------------------------------------- expression tests --------
# Only joints where a fusion call is at stake: splice-precise, not an artifact, donor-restricted, and
# seen in >= MIN_POS_CELLS cells. Cell-level values; labels permuted within strata of library, cell
# state and depth quintile; Benjamini-Hochberg across the tested joints for each test.
candidates = [x for x in rows if not x['blacklisted_all'] and not x['paralog_partner'] and x['ts_status'] != 'template_switch_likely'
              and x['splice_precise'] and x['donor_restricted'] and x['positive_cells'] >= MIN_POS_CELLS
              and x['gene_a'] in genes and x['gene_b'] in genes]
by_lib = collections.defaultdict(list)
for x in candidates:
    for lib, bcs in x['_pos_cells'].items():
        by_lib[lib].append((x['_key'], bcs))


def expression_cells(lib):
    """Per joint: per-cell records for joint-positive cells and negative cells of the same cell states."""
    import anndata as ad
    from scipy import sparse
    a = ad.read_h5ad(ISV / f'{lib}-isoform.h5ad')
    obs = a.obs['cluster'].astype(str)
    cell_index = {b: i for i, b in enumerate(a.obs_names)}
    var_gene = [v.split(':', 1)[0] for v in a.var_names]
    gene_codes = {g: k for k, g in enumerate(sorted(set(var_gene)))}
    indicator = sparse.csr_matrix((np.ones(len(var_gene)), (np.arange(len(var_gene)), [gene_codes[g] for g in var_gene])),
                                  shape=(len(var_gene), len(gene_codes)))
    X_ = a.X.tocsr() if sparse.issparse(a.X) else sparse.csr_matrix(a.X)
    gene_X = (X_ @ indicator).tocsr()
    totals = np.asarray(X_.sum(axis=1)).ravel()
    bam = pysam.AlignmentFile(samples[lib]['bam'])
    results = {}
    for key, pos in by_lib[lib]:
        ga, ca, sa, d, gb, cb, sb, acc = key
        pos_in = [b for b in pos if b in cell_index]
        states = {obs.iloc[cell_index[b]] for b in pos_in}
        neg = [b for b, i in cell_index.items() if obs.iloc[i] in states and b not in pos]
        cells = pos_in + neg
        if not pos_in or not neg:
            results[key] = dict(library=lib, cells=[])
            continue
        where = {b: k for k, b in enumerate(cells)}
        counts = np.zeros((len(cells), 4))

        def read_side(gene, chrom, strand, bp, side, col):
            g = genes[gene]
            for read in bam.fetch(chrom, g['start'] - 1, g['end']):
                if read.is_secondary or read.is_supplementary or read.is_unmapped or read.has_tag('SA') or not read.has_tag('CB'):
                    continue
                if ('-' if read.is_reverse else '+') != strand:
                    continue
                bc = read.get_tag('CB'); bc = bc if '-' in bc else bc + '-1'
                k = where.get(bc)
                if k is None:
                    continue
                blocks = read.get_blocks()
                if side == 'after_donor':
                    hit = any((e > bp) if strand == '+' else (s0 + 1 < bp) for s0, e in blocks)
                else:
                    hit = any((s0 + 1 < bp) if strand == '+' else (e > bp) for s0, e in blocks)
                counts[k, col] += hit
                counts[k, col + 1] += 1
        read_side(ga, ca, sa, d, 'after_donor', 0)
        read_side(gb, cb, sb, acc, 'before_acceptor', 2)
        flank = [g for g, s in gene_index.overlaps(ca, d - FLANK_WINDOW, d + FLANK_WINDOW) if genes[g]['gene_type'] == 'protein_coding' and g != ga]
        flank += [g for g, s in gene_index.overlaps(cb, acc - FLANK_WINDOW, acc + FLANK_WINDOW) if genes[g]['gene_type'] == 'protein_coding' and g != gb]
        cols = [gene_codes[g] for g in set(flank) if g in gene_codes]
        idx = [cell_index[b] for b in cells]
        flank_score = (np.log1p(1e4 * np.asarray(gene_X[idx][:, cols].sum(axis=1)).ravel() / np.maximum(totals[idx], 1))
                       if cols else np.full(len(cells), np.nan))
        results[key] = dict(library=lib, flank_genes=len(set(flank)), cells=[
            (b in pos, obs.iloc[cell_index[b]], float(totals[cell_index[b]]), *counts[k], float(flank_score[k]))
            for k, b in enumerate(cells)])
    bam.close()
    return lib, results


def run_tests(item):
    key, parts = item
    cells, strata_raw = [], []
    for part in parts:
        lib = part['library']
        by_state = collections.defaultdict(list)
        for c in part['cells']:
            by_state[c[1]].append(c)
        for state, cs in by_state.items():
            depth = np.array([c[2] for c in cs])
            edges = np.quantile(depth, [0.2, 0.4, 0.6, 0.8]) if len(cs) >= 5 else []
            for c in cs:
                cells.append(c)
                strata_raw.append(f'{lib}|{state}|{int(np.searchsorted(edges, c[2], side="right"))}')
    if not cells:
        return key, {}
    pos = [c[0] for c in cells]
    col = lambda k: [c[k] for c in cells]
    out = dict(positive_cells_tested=sum(pos), negative_cells_tested=len(pos) - sum(pos),
               a_after_donor_reads_positive=sum(c[3] for c in cells if c[0]), a_reads_positive=sum(c[4] for c in cells if c[0]),
               a_after_donor_reads_negative=sum(c[3] for c in cells if not c[0]), a_reads_negative=sum(c[4] for c in cells if not c[0]),
               b_before_acceptor_reads_positive=sum(c[5] for c in cells if c[0]), b_reads_positive=sum(c[6] for c in cells if c[0]),
               b_before_acceptor_reads_negative=sum(c[5] for c in cells if not c[0]), b_reads_negative=sum(c[6] for c in cells if not c[0]))
    out['a_stat'], out['a_p'] = R.stratified_permutation_lower(col(3), col(4), pos, strata_raw, n_perm=N_PERM)
    out['b_stat'], out['b_p'] = R.stratified_permutation_lower(col(5), col(6), pos, strata_raw, n_perm=N_PERM)
    fl = col(7)
    if not any(np.isnan(fl)):
        out['flank_stat'], out['flank_p'] = R.stratified_permutation_lower(fl, None, pos, strata_raw, n_perm=N_PERM)
    return key, out


with Pool(nproc) as pool:
    expr = pool.map(expression_cells, sorted(by_lib))
parts_by_key = collections.defaultdict(list)
flank_genes = {}
for lib, res in expr:
    for key, part in res.items():
        parts_by_key[key].append(part)
        flank_genes[key] = part.get('flank_genes')
with Pool(nproc) as pool:
    tested = dict(pool.map(run_tests, sorted(parts_by_key.items())))
keys = sorted(tested)
for test in ('a', 'b', 'flank'):
    qs = R.bh_qvalues([tested[k].get(f'{test}_p') for k in keys])
    for k, q in zip(keys, qs):
        tested[k][f'{test}_q'] = q
expr_rows = []
for x in rows:
    t = tested.get(x['_key'])
    imbalance = flank_decrease = False
    if t:
        imbalance = any(t.get(f'{s}_q') is not None and t[f'{s}_q'] < Q_EXPRESSION for s in ('a', 'b'))
        flank_decrease = t.get('flank_q') is not None and t['flank_q'] < Q_EXPRESSION
        x.update({f'expression_{k}': v for k, v in t.items()}, expression_flank_genes=flank_genes.get(x['_key']))
        expr_rows.append(dict(joint=x['joint'], gene_a=x['gene_a'], gene_b=x['gene_b'], symbol_a=x['symbol_a'], symbol_b=x['symbol_b'],
                              top_donor=x['top_donor'], flank_genes=flank_genes.get(x['_key']), **t))
    x.update(expression_tested=bool(t), partner_exon_imbalance=imbalance, flanking_expression_decrease=flank_decrease)
    x['donors'] = int(x['donors'])
    x['class'] = R.classify(x)
    if x['class'] == 'gene_fusion_candidate':
        x['fusion_evidence'] = ','.join(n for n, v in (('curated_database', x['fusion_db_curated']),
                                                        ('reciprocal_precise', x['donor_restricted'] and x['reciprocal_precise_same_donor']),
                                                        ('recurrent_intronic_breakpoint', x['donor_restricted'] and x['intronic_breakpoint_recurrent_same_donor']),
                                                        ('partner_exon_loss', x['donor_restricted'] and imbalance),
                                                        ('flanking_expression_decrease', x['donor_restricted'] and flank_decrease)) if v)
    else:
        x['fusion_evidence'] = ''

# ---------------------------------------------------------------- outputs -----------------
order = ['trans_splicing_candidate', 'gene_fusion_candidate', 'unresolved', 'artifact_template_switch',
         'artifact_not_splice_precise', 'artifact_paralog_sequence', 'artifact_blacklist']
rows.sort(key=lambda x: (order.index(x['class']), -x['donors'], -int(x['molecules'])))
fields = []
for x in rows:  # union of columns in first-seen order: tested joints carry extra expression_* columns
    for k in x:
        if not k.startswith('_') and k not in fields:
            fields.append(k)
with gzip.open(OUT / 'joint_classification.tsv.gz', 'wt') as h:
    w = csv.DictWriter(h, fieldnames=fields, delimiter='\t', extrasaction='ignore')
    w.writeheader(); w.writerows(rows)
with open(OUT / 'joint_classification_candidates.tsv', 'w') as h:
    w = csv.DictWriter(h, fieldnames=fields, delimiter='\t', extrasaction='ignore')
    w.writeheader(); w.writerows(x for x in rows if x['class'] in ('trans_splicing_candidate', 'gene_fusion_candidate'))
if expr_rows:
    with open(OUT / 'expression_tests.tsv', 'w') as h:
        keys = sorted({k for e in expr_rows for k in e}, key=lambda k: (k not in ('joint', 'gene_a', 'gene_b', 'library'), k))
        w = csv.DictWriter(h, fieldnames=keys, delimiter='\t')
        w.writeheader(); w.writerows(expr_rows)
summary = dict(run_root=str(root), joints=len(rows), thresholds=dict(genomic_homology_p99=THR_G, mrna_homology_p99=THR_M,
               position_p=P_POSITION, expression_bh_q=Q_EXPRESSION, donor_restricted_share=DONOR_RESTRICTED,
               paralog_identity=PARALOG_ID, paralog_min_aligned=PARALOG_LEN, min_positive_cells_for_expression=MIN_POS_CELLS,
               min_cells_in_top_donor_for_donor_restriction=MIN_CLONE_CELLS, intronic_site_bin_nt=5,
               flank_window=FLANK_WINDOW, permutations=N_PERM),
               classes=dict(collections.Counter(x['class'] for x in rows)),
               ts_status=dict(collections.Counter(x['ts_status'] for x in rows)),
               expression_tested_joints=sum(x['expression_tested'] for x in rows),
               classes_with_2plus_donors=dict(collections.Counter(x['class'] for x in rows if x['donors'] >= 2)),
               fusion_db_curated_joints=sum(bool(x['fusion_db_curated']) for x in rows),
               seconds=round(time.time() - started, 1))
(OUT / 'CLASSIFICATION_SUMMARY.json').write_text(json.dumps(summary, indent=1, default=str) + '\n')
print(json.dumps(summary, indent=1, default=str))
