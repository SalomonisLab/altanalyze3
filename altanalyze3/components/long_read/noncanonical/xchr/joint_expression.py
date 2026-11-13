"""Expression, cell-state predominance and cohort predominance of cross-chromosome joints with at least
MIN_READS molecules in each of at least MIN_SAMPLES samples (EXPERIMENTAL extra-chr-trans-splicing).

Thresholds (user, 2026-10-30): "minimum of 50 reads in at least one sample" (50, 1) and "> 5 reads in > 4
samples" (6, 5). A read is one deduplicated molecule (isoseq groupdedup); a sample is one original_uid (split
BAMs of one library pool). Every class passes through and is labelled from the classification folder; the
threshold is the only filter.
Cells: cellHarmony-annotated cells only (isv/combined_cellHarmony_unique_barcodes.txt). Junction-table and
h5ad barcodes map to the annotation orientation by forward or reverse complement, whichever matches more.
Depth: isoform counts of the annotated cells in the standard per-library isoform h5ad (isv/<library>-isoform.h5ad).
Units for cohort tests: patient x cohort (config/patients.tsv); timepoints of one patient in one cohort pool.
1 expression: molecules per sample, positive cells, molecules per million depth, molecules per 10,000
  molecules of the donor gene and of the acceptor gene in the same cells.
2 cell state: positive cells per state against the within-library expectation (a positive cell of original
  library L falls in state s with probability depth(L,s)/depth(L)); exact upper tail of the sum of binomials;
  Benjamini-Hochberg over all joint x state tests.
3 cohort: O/E per cohort, E = sum over states of unit depth x the joint's pooled state rate (removes cell-state
  composition); contrasts aged vs young, MDS vs aged, AML vs aged; statistic = difference of pooled O/E; exact
  two-sided permutation of unit labels; Benjamini-Hochberg per contrast over testable joints.
Arguments: run root, classification folder, output folder, patients TSV, worker count, minimum reads per sample
(inclusive, default 50), minimum samples (inclusive, default 1).
"""
import collections
import csv
import gzip
import itertools
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path
from altanalyze3.components.long_read.noncanonical.xchr.common import config_path

import numpy as np
from scipy import stats as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
from altanalyze3.components.long_read.noncanonical.xchr.common import KINNEX, SAMPLES, BARCODE_CLUSTERS, inside_x, refuse_existing
from altanalyze3.components.bam import cross_chromosome_rules as R

csv.field_size_limit(sys.maxsize)
root = inside_x(Path(sys.argv[1]))
CLASS_DIR = root / sys.argv[2]
OUT = root / sys.argv[3]
PATIENTS = Path(sys.argv[4])
nproc = int(sys.argv[5])
MIN_READS = int(sys.argv[6]) if len(sys.argv) > 6 else 50
MIN_SAMPLES = int(sys.argv[7]) if len(sys.argv) > 7 else 1
refuse_existing(OUT / 'EXPRESSION_SUMMARY.json')
OUT.mkdir(exist_ok=True)
ISV = config_path('expression_dir')
COHORTS = ['young', 'aged', 'MDS', 'AML']
CONTRASTS = [('aged', 'young'), ('MDS', 'aged'), ('AML', 'aged')]
ORDER = ['HSC-1', 'HSC-2', 'MPP-1', 'MPP-2', 'MPP-MEP', 'LMPP-1', 'LMPP-1-cycling', 'Multilin-1', 'Multilin-2', 'MEP-1', 'MEP-2',
         'MEP-Eryth-1', 'MEP-Eryth-2', 'ERP-1', 'ERP-2', 'ERP-3', 'ERP-4', 'ERP-5', 'ERP-6', 'ERP-7', 'ERP-8', 'Erythroblast-1',
         'Erythroblast-2', 'Erythroblast-3', 'MKP-early', 'MKP-late', 'MK-Platelet', 'BMCP-1', 'BMCP-2', 'LMPP-2', 'CLP',
         'Pro-B-Early', 'Pro-B-1', 'Pro-B-cycling-2', 'Pro-B-Early-cycling', 'Pro-B-cycling-1', 'Transitional-B-1',
         'Transitional-B-2', 'Pro-B-2', 'Pro-B-3', 'pre-B', 'B Memory-1', 'B Memory-2', 'Plasma Cell', 'CD4 TCM',
         'T CD4 Activated', 'T CD4 Naive-1', 'T CD4 Naive-2', 'T CD8 Naive', 'CD8 TEM', 'MAIT', 'NK-Mature', 'MultiLin-GMP-1',
         'MultiLin-GMP-2', 'MultiLin-GMP-3', 'preNeu', 'immNeu-1', 'immNeu-2', 'cMOP', 'Mono-1', 'Mono-2', 'Intermediate Mono-1',
         'Intermediate Mono-2', 'Intermediate Mono-3', 'Classical-Mono', 'Non-Classical Mono-1', 'Non-Classical Mono-2',
         'Myeloid intermediate 1', 'Myeloid intermediate 2', 'Myeloid intermediate 3', 'MDP-1', 'MDP-2', 'MDP-3', 'MDP-4',
         'MDP-5', 'MDP-6', 'pre-DC-1', 'pre-DC-2', 'pre-DC-3', 'cDC1', 'cDC2-1', 'cDC2-2', 'pDC', 'ASDC', 'Mac',
         'MSC Fibroblast-1', 'MSC Fibroblast-2', 'Osteoblast', 'Stromal Vascular']
started = time.time()
samples = {r['library']: r for r in csv.DictReader(open(SAMPLES), delimiter='\t')}
patient_of = {r['library']: r['patient'] for r in csv.DictReader(open(PATIENTS), delimiter='\t')}
libs = [l for l in samples if (root / 'samples' / l / 'STAGE1_COMPLETE.json').exists()]
missing = [l for l in libs if l not in patient_of]
if missing:
    raise SystemExit(f'{PATIENTS} lacks libraries: {missing}')
uid_of = {l: samples[l]['original_uid'] for l in libs}
olib_of = {l: samples[l]['original_library'] for l in libs}
cohort_of = {l: samples[l]['cohort'] for l in libs}
unit_of = {l: f'{patient_of[l]}|{cohort_of[l]}' for l in libs}
unk = lambda g: 'UNK' if g.startswith('UNK') else g


def rc(b):
    core, _, suffix = b.partition('-')
    return core.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1] + ('-' + suffix if suffix else '')


# ---------------------------------------------------------------- joints and classes -------
joints = {}
with gzip.open(CLASS_DIR / 'joint_classification.tsv.gz', 'rt') as h:
    for r in csv.DictReader(h, delimiter='\t'):
        left, right = r['joint'].split(' > ')
        ca, d, sa = left.rsplit(':', 2)
        cb, a, sb = right.rsplit(':', 2)
        joints[(r['gene_a'], ca, sa, int(d), r['gene_b'], cb, sb, int(a))] = {
            k: r.get(k, '') for k in ('gene_a', 'symbol_a', 'gene_b', 'symbol_b', 'joint', 'class', 'donors', 'sample_donors',
                                      'molecules', 'molecules_scanned', 'frame_class', 'valid_coding', 'ts_status', 'splice_precise',
                                      'junction_homology_nt', 'mrna_homology_nt', 'paralog_partner', 'fusion_db_curated',
                                      'fusion_db_tags', 'donor_restricted', 'top_donor', 'top_donor_share')}
eligible = {k for k, j in joints.items() if int(j['molecules']) >= MIN_READS * MIN_SAMPLES}
clusters = collections.defaultdict(dict)
with open(BARCODE_CLUSTERS) as h:
    for line in h:
        key, cluster = line.rstrip('\n').split('\t')[:2]
        barcode, lib = key.split('.', 1)
        clusters[lib].setdefault(barcode, cluster)


# ---------------------------------------------------------------- one pass per library -----
def scan(lib):
    import anndata as ad
    from scipy import sparse
    ann = clusters.get(lib, {})
    counts = collections.Counter()
    rows = []
    bcs = set()
    path = root / 'samples' / lib / f'{lib}.cross_chromosome_junctions.tsv.gz'
    with gzip.open(path, 'rt') as h:
        head = h.readline().rstrip('\n').split('\t')
        i = {c: k for k, c in enumerate(head)}
        for line in h:
            p = line.rstrip('\n').split('\t')
            if p[i['cross_chromosome']] != 'True':
                continue
            key = (unk(p[i['gene_a']]), p[i['a_chrom']], p[i['a_strand']], int(p[i['donor']]),
                   unk(p[i['gene_b']]), p[i['b_chrom']], p[i['b_strand']], int(p[i['acceptor']]))
            if key not in joints:
                continue
            counts[key] += 1
            b = p[i['barcode']]
            if b:
                bcs.add(b)
            if key in eligible:
                rows.append((key, b))
    fwd, rev = sum(b in ann for b in bcs), sum(rc(b) in ann for b in bcs)
    j_orient = rc if rev > fwd else (lambda b: b)
    a = ad.read_h5ad(ISV / f'{lib}-isoform.h5ad')
    names = list(a.obs_names)
    hf, hr = sum(b in ann for b in names), sum(rc(b) in ann for b in names)
    h_orient = rc if hr > hf else (lambda b: b)
    mapped = [h_orient(b) for b in names]
    state = [ann.get(b, '') for b in mapped]
    obs_cluster = a.obs['cluster'].astype(str).tolist() if 'cluster' in a.obs else [''] * len(names)
    annotated = [k for k, s in enumerate(state) if s]
    agree = sum(state[k] == obs_cluster[k] for k in annotated)
    X_ = a.X.tocsr() if sparse.issparse(a.X) else sparse.csr_matrix(a.X)
    totals = np.asarray(X_.sum(axis=1)).ravel()
    depth, cells = collections.Counter(), collections.Counter()
    for k in annotated:
        depth[state[k]] += float(totals[k]); cells[state[k]] += 1
    cell_depth = {mapped[k]: float(totals[k]) for k in annotated}
    genes_needed = {joints[k][g] for k in eligible for g in ('gene_a', 'gene_b')} - {'UNK'}
    var_gene = [v.split(':', 1)[0] for v in a.var_names]
    gene_counts = collections.defaultdict(collections.Counter)
    cols = [c for c, g in enumerate(var_gene) if g in genes_needed]
    if cols and annotated:
        sub = X_[annotated][:, cols].tocsc()
        st_ann = [state[k] for k in annotated]
        for n, c in enumerate(cols):
            col = sub.getcol(n)
            for r_, v in zip(col.indices, col.data):
                gene_counts[var_gene[c]][st_ann[r_]] += float(v)
    per_key = collections.defaultdict(lambda: dict(total=0, annotated=0, with_depth=0, mol=collections.Counter(), cells=collections.defaultdict(set)))
    for key, b in rows:
        x = per_key[key]
        x['total'] += 1
        if not b:
            continue
        cb = j_orient(b)
        s = ann.get(cb, '')
        if not s:
            continue
        x['annotated'] += 1
        if cb not in cell_depth:
            continue
        x['with_depth'] += 1
        x['mol'][s] += 1
        x['cells'][s].add(cb)
    info = dict(library=lib, junction_barcodes=len(bcs), junction_forward=fwd, junction_reverse_complement=rev,
                h5ad_cells=len(names), h5ad_forward=hf, h5ad_reverse_complement=hr, h5ad_annotated=len(annotated),
                h5ad_cluster_agrees=agree, annotation_cells=len(ann), depth_annotated=float(sum(depth.values())))
    return lib, dict(counts), dict(per_key), dict(depth), dict(cells), {g: dict(c) for g, c in gene_counts.items()}, info


with Pool(nproc) as pool:
    scanned = pool.map(scan, libs)

# ---------------------------------------------------------------- conservation ------------
lib_counts = {lib: c for lib, c, *_ in scanned}
total_scan = collections.Counter()
for c in lib_counts.values():
    total_scan.update(c)
mismatch = [k for k, j in joints.items() if j['molecules_scanned'] != '' and total_scan.get(k, 0) != int(j['molecules_scanned'])]
if mismatch:
    raise SystemExit(f'{len(mismatch)} joints: scanned molecules differ from the classification molecules_scanned, e.g. {mismatch[:3]}')
per_sample = collections.defaultdict(collections.Counter)
per_library_max = collections.Counter()
libs_ge = collections.Counter()
for lib, c in lib_counts.items():
    for k, n in c.items():
        per_sample[k][uid_of[lib]] += n
        per_library_max[k] = max(per_library_max[k], n)
        libs_ge[k] += n >= MIN_READS
max_sample = {k: max(v.values()) for k, v in per_sample.items()}
samples_ge = {k: sum(n >= MIN_READS for n in v.values()) for k, v in per_sample.items()}
passing = sorted(k for k, m in samples_ge.items() if m >= MIN_SAMPLES)
if not set(passing) <= eligible:
    raise SystemExit('a passing joint was not eligible for the cell scan; eligibility rule is wrong')
class_dist = collections.defaultdict(collections.Counter)
for k, j in joints.items():
    class_dist[j['class']]['joints'] += 1
    class_dist[j['class']]['pass_sample'] += samples_ge.get(k, 0) >= MIN_SAMPLES
    class_dist[j['class']]['pass_library'] += libs_ge.get(k, 0) >= MIN_SAMPLES
max_by_class = {c: max((max_sample.get(k, 0) for k, j in joints.items() if j['class'] == c), default=0) for c in class_dist}
max_samples_ge_by_class = {c: max((samples_ge.get(k, 0) for k, j in joints.items() if j['class'] == c), default=0) for c in class_dist}

# ---------------------------------------------------------------- depth tables ------------
states = set()
D_lib = {}
cells_lib = {}
for lib, _, _, depth, cells, _, _ in scanned:
    D_lib[lib], cells_lib[lib] = depth, cells
    states.update(depth)
states = [s for s in ORDER if s in states] + sorted(s for s in states if s not in ORDER)
units = sorted({unit_of[l] for l in libs}, key=lambda u: (COHORTS.index(u.split('|')[1]), u))
olibs = sorted({olib_of[l] for l in libs})
sidx = {s: k for k, s in enumerate(states)}
D_unit = np.zeros((len(units), len(states)))
D_olib = np.zeros((len(olibs), len(states)))
for lib in libs:
    for s, v in D_lib[lib].items():
        D_unit[units.index(unit_of[lib]), sidx[s]] += v
        D_olib[olibs.index(olib_of[lib]), sidx[s]] += v
gene_total = collections.Counter()
for x in scanned:
    for g, c in x[5].items():
        gene_total[g] += sum(c.values())
total_depth = float(D_unit.sum())


# ---------------------------------------------------------------- tests --------------------
def sum_binom_upper(ns, ps, obs):
    pmf = np.array([1.0])
    for n, p in zip(ns, ps):
        pmf = np.convolve(pmf, st.binom.pmf(np.arange(n + 1), n, p))
    return float(min(1.0, pmf[obs:].sum())), float(np.dot(ns, ps))


def state_tests(cells_by_olib_state):
    """cells_by_olib_state: array (olibs x states) of positive cells. Returns list of (state, O, E, p)."""
    n_l = cells_by_olib_state.sum(axis=1)
    use = np.where(n_l > 0)[0]
    out = []
    if not len(use):
        return out
    P = D_olib[use] / D_olib[use].sum(axis=1, keepdims=True)
    for s in range(len(states)):
        obs = int(cells_by_olib_state[use, s].sum())
        if not (P[:, s] > 0).any():
            if obs:
                raise SystemExit(f'positive cells in state {states[s]} with zero depth in their libraries')
            continue
        p, e = sum_binom_upper(n_l[use].astype(int), P[:, s], obs)
        out.append((states[s], obs, e, p))
    return out


combos = {}
for a_, b_ in CONTRASTS:
    ua = [k for k, u in enumerate(units) if u.endswith('|' + a_)]
    ub = [k for k, u in enumerate(units) if u.endswith('|' + b_)]
    idx = ua + ub
    rows_ = [[1.0 if k in set(c) else 0.0 for k in range(len(idx))] for c in itertools.combinations(range(len(idx)), len(ua))]
    combos[(a_, b_)] = (idx, len(ua), np.array(rows_))


def cohort_test(O_us):
    """O_us: molecules (units x states). Returns per-cohort O/E and per-contrast (A O/E, B O/E, diff, p)."""
    O_u = O_us.sum(axis=1)
    r_s = np.divide(O_us.sum(axis=0), D_unit.sum(axis=0), out=np.zeros(len(states)), where=D_unit.sum(axis=0) > 0)
    E_u = D_unit @ r_s
    res = {}
    for c in COHORTS:
        m = np.array([u.endswith('|' + c) for u in units])
        res[f'{c}_units'] = int(m.sum())
        res[f'{c}_units_detected'] = int((O_u[m] > 0).sum())
        res[f'{c}_molecules'] = int(O_u[m].sum())
        res[f'{c}_per_million'] = float(O_u[m].sum() / D_unit[m].sum() * 1e6)
        res[f'{c}_obs_over_exp'] = float(O_u[m].sum() / E_u[m].sum()) if E_u[m].sum() > 0 else None
    for (a_, b_), (idx, na, Cm) in combos.items():
        name = f'{a_}_vs_{b_}'
        O, E = O_u[idx], E_u[idx]
        if O.sum() == 0:
            res.update({f'{name}_diff': None, f'{name}_p': None, f'{name}_permutations': len(Cm)})
            continue
        OA, EA = Cm @ O, Cm @ E
        OB, EB = O.sum() - OA, E.sum() - EA
        ok = (EA > 0) & (EB > 0)
        stat = np.full(len(Cm), np.nan)
        stat[ok] = OA[ok] / EA[ok] - OB[ok] / EB[ok]
        obs_row = np.zeros(len(idx)); obs_row[:na] = 1
        obs = float(obs_row @ O / (obs_row @ E) - (O.sum() - obs_row @ O) / (E.sum() - obs_row @ E))
        finite = stat[np.isfinite(stat)]
        res.update({f'{name}_diff': obs, f'{name}_p': float(np.mean(np.abs(finite) >= abs(obs) - 1e-12)),
                    f'{name}_permutations': int(len(finite)), f'{name}_undefined_permutations': int(len(Cm) - len(finite))})
    return res, O_u, E_u


# ---------------------------------------------------------------- per joint ----------------
summary_rows, state_rows, unit_rows = [], [], []
by_key = collections.defaultdict(list)
for lib, _, per_key, *_ in scanned:
    for k, x in per_key.items():
        by_key[k].append((lib, x))
dispersion = []
anomalies = []
for k in passing:
    j = joints[k]
    O_us = np.zeros((len(units), len(states)))
    C_ls = np.zeros((len(olibs), len(states)))
    cellsets = collections.defaultdict(set)
    tot = ann_n = dep_n = 0
    for lib, x in by_key.get(k, []):
        tot += x['total']; ann_n += x['annotated']; dep_n += x['with_depth']
        for s, n in x['mol'].items():
            O_us[units.index(unit_of[lib]), sidx[s]] += n
        for s, bs in x['cells'].items():
            cellsets[(olib_of[lib], s)].update(bs)
    for (ol, s), bs in cellsets.items():
        C_ls[olibs.index(ol), sidx[s]] = len(bs)
        if D_olib[olibs.index(ol), sidx[s]] <= 0:
            anomalies.append((j['joint'], ol, s))
    tests = state_tests(C_ls)
    coh, O_u, E_u = cohort_test(O_us)
    m = E_u > 0
    if O_u.sum() >= 20 and m.sum() > 1:
        dispersion.append(float(((O_u[m] - E_u[m]) ** 2 / E_u[m]).sum() / (m.sum() - 1)))
    smp = per_sample[k]
    row = dict(j, threshold_max_sample_molecules=max_sample[k],
               threshold_top_sample=max(smp, key=smp.get), samples_ge_threshold=sum(v >= MIN_READS for v in smp.values()),
               samples_detected=len(smp), molecules_all_cells=tot, molecules_annotated_cells=ann_n,
               molecules_annotated_with_depth=dep_n, positive_cells=int(C_ls.sum()),
               per_million_depth=dep_n / total_depth * 1e6,
               per_10k_gene_a=(dep_n / gene_total[j['gene_a']] * 1e4) if gene_total.get(j['gene_a']) else None,
               per_10k_gene_b=(dep_n / gene_total[j['gene_b']] * 1e4) if gene_total.get(j['gene_b']) else None,
               molecules_per_sample=';'.join(f'{u}:{n}' for u, n in smp.most_common()), **coh)
    for s, o, e, p in tests:
        pat = sorted({units[u].split('|')[0] for u in range(len(units)) if O_us[u, sidx[s]] > 0})
        state_rows.append(dict(joint=j['joint'], symbol_a=j['symbol_a'], symbol_b=j['symbol_b'], **{'class': j['class']},
                               state=s, positive_cells=o, expected_cells=e, obs_over_exp=(o / e) if e > 0 else None, p=p,
                               molecules=int(O_us[:, sidx[s]].sum()), patients_with_molecules=len(pat), patients=','.join(pat)))
    for u in range(len(units)):
        unit_rows.append(dict(joint=j['joint'], symbol_a=j['symbol_a'], symbol_b=j['symbol_b'], unit=units[u],
                              patient=units[u].split('|')[0], cohort=units[u].split('|')[1], molecules=int(O_u[u]),
                              depth=float(D_unit[u].sum()), expected=float(E_u[u]),
                              per_million=float(O_u[u] / D_unit[u].sum() * 1e6) if D_unit[u].sum() > 0 else None))
    summary_rows.append(row)

qs = R.bh_qvalues([r['p'] for r in state_rows])
for r, q in zip(state_rows, qs):
    r['q'] = q
for a_, b_ in CONTRASTS:
    name = f'{a_}_vs_{b_}'
    for r, q in zip(summary_rows, R.bh_qvalues([r[f'{name}_p'] for r in summary_rows])):
        r[f'{name}_q'] = q
for r in summary_rows:
    mine = [x for x in state_rows if x['joint'] == r['joint']]
    best = min(mine, key=lambda x: (x['p'], -(x['obs_over_exp'] or 0))) if mine else None
    r.update(top_state=best['state'] if best else '', top_state_cells=best['positive_cells'] if best else None,
             top_state_expected=best['expected_cells'] if best else None, top_state_obs_over_exp=best['obs_over_exp'] if best else None,
             top_state_p=best['p'] if best else None, top_state_q=best['q'] if best else None,
             top_state_patients=best['patients_with_molecules'] if best else None,
             states_q_lt_0_05=';'.join(f"{x['state']}({x['positive_cells']}/{x['expected_cells']:.1f})" for x in mine if x['q'] < 0.05),
             state_cells=';'.join(f"{x['state']}:{x['positive_cells']}" for x in sorted(mine, key=lambda x: -x['positive_cells'])[:8] if x['positive_cells']))
summary_rows.sort(key=lambda r: (-r['samples_ge_threshold'], -r['molecules_annotated_with_depth']))

# ---------------------------------------------------------------- calibration and controls -
rng = np.random.default_rng(0)
null_p = []
for k in passing:
    C_ls = np.zeros((len(olibs), len(states)))
    for lib, x in by_key.get(k, []):
        for s, bs in x['cells'].items():
            C_ls[olibs.index(olib_of[lib]), sidx[s]] += len(bs)
    n_l = C_ls.sum(axis=1)
    for rep in range(200):
        sim = np.zeros_like(C_ls)
        for L in np.where(n_l > 0)[0]:
            sim[L] = rng.multinomial(int(n_l[L]), D_olib[L] / D_olib[L].sum())
        null_p.extend(p for _, _, _, p in state_tests(sim))
null_p = np.array(null_p)
young_units = [u for u in range(len(units)) if units[u].endswith('|young')]
aged_units = [u for u in range(len(units)) if units[u].endswith('|aged')]
plant = np.zeros((len(units), len(states)))
for u in aged_units:
    plant[u] = rng.multinomial(10, D_unit[u] / D_unit[u].sum())
plant_res, _, _ = cohort_test(plant)
null_cohort = []
for rep in range(200):
    sim = np.zeros((len(units), len(states)))
    tot = rng.multinomial(60, D_unit.sum(axis=1) / D_unit.sum())
    for u in range(len(units)):
        sim[u] = rng.multinomial(int(tot[u]), D_unit[u] / D_unit[u].sum())
    res, _, _ = cohort_test(sim)
    null_cohort.append([res[f'{a}_vs_{b}_p'] for a, b in CONTRASTS])
null_cohort = np.array(null_cohort, dtype=float)
# planted state enrichment: 8 positive cells of one library, all in its smallest non-zero-depth state
state_plant = np.zeros((len(olibs), len(states)))
L0 = 0
share = np.where(D_olib[L0] > 0, D_olib[L0] / D_olib[L0].sum(), np.inf)
planted_state = states[int(np.argmin(share))]
state_plant[L0, sidx[planted_state]] = 8
planted_state_p = [p for s, o, e, p in state_tests(state_plant) if s == planted_state][0]

# ---------------------------------------------------------------- outputs -----------------
def write(name, rows):
    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with open(OUT / name, 'w') as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter='\t')
        w.writeheader(); w.writerows(rows)


write('joints_min_reads_summary.tsv', summary_rows)
write('joints_min_reads_state_tests.tsv', state_rows)
write('joints_min_reads_unit_molecules.tsv', unit_rows)
write('library_state_depth.tsv', [dict(library=l, original_library=olib_of[l], unit=unit_of[l], state=s, cells=cells_lib[l].get(s, 0), depth=v)
                                  for l in libs for s, v in sorted(D_lib[l].items())])
write('library_barcode_mapping.tsv', [x[6] for x in scanned])
with gzip.open(OUT / 'max_molecules_per_sample_all_joints.tsv.gz', 'wt') as h:
    w = csv.writer(h, delimiter='\t')
    w.writerow(['joint', 'symbol_a', 'symbol_b', 'class', 'molecules', 'max_sample_molecules', 'max_sample', 'max_library_molecules',
                f'samples_ge_{MIN_READS}', f'libraries_ge_{MIN_READS}', 'passes'])
    for k, j in joints.items():
        smp = per_sample.get(k, {})
        w.writerow([j['joint'], j['symbol_a'], j['symbol_b'], j['class'], j['molecules'], max_sample.get(k, 0),
                    max(smp, key=smp.get) if smp else '', per_library_max.get(k, 0), samples_ge.get(k, 0), libs_ge.get(k, 0),
                    samples_ge.get(k, 0) >= MIN_SAMPLES])
info = [x[6] for x in scanned]
summary = dict(run_root=str(root), classification=str(CLASS_DIR), patients=str(PATIENTS), annotation=str(BARCODE_CLUSTERS),
               threshold=dict(min_reads_per_sample=MIN_READS, min_samples=MIN_SAMPLES,
                              unit='molecules in one original_uid (split BAMs pooled); both bounds inclusive'),
               max_samples_meeting_min_reads_by_class=max_samples_ge_by_class,
               joints=len(joints), eligible_by_total=len(eligible), passing=len(passing),
               classes=dict((c, dict(v)) for c, v in class_dist.items()), max_sample_molecules_by_class=max_by_class,
               conservation_scanned_equals_classification=len(mismatch) == 0,
               units=units, unit_counts={c: sum(u.endswith('|' + c) for u in units) for c in COHORTS},
               permutations={f'{a}_vs_{b}': len(combos[(a, b)][2]) for a, b in CONTRASTS},
               min_attainable_p={f'{a}_vs_{b}': 1 / len(combos[(a, b)][2]) for a, b in CONTRASTS},
               state_tests=len(state_rows), state_tests_q_lt_0_05=sum(r['q'] < 0.05 for r in state_rows),
               cohort_p_lt_0_05={f'{a}_vs_{b}': sum(r[f'{a}_vs_{b}_p'] is not None and r[f'{a}_vs_{b}_p'] < 0.05 for r in summary_rows) for a, b in CONTRASTS},
               cohort_q_lt_0_05={f'{a}_vs_{b}': sum(r[f'{a}_vs_{b}_q'] is not None and r[f'{a}_vs_{b}_q'] < 0.05 for r in summary_rows) for a, b in CONTRASTS},
               cohort_untestable={f'{a}_vs_{b}': sum(r[f'{a}_vs_{b}_p'] is None for r in summary_rows) for a, b in CONTRASTS},
               dispersion_pearson_chi2_per_df=dict(joints=len(dispersion), median=float(np.median(dispersion)) if dispersion else None,
                                                   above_1_5=int(sum(d > 1.5 for d in dispersion))),
               depth_anomalies_positive_cell_in_zero_depth_state=len(anomalies),
               barcode_mapping=dict(junction_reverse_complement_libraries=sum(x['junction_reverse_complement'] > x['junction_forward'] for x in info),
                                    h5ad_reverse_complement_libraries=sum(x['h5ad_reverse_complement'] > x['h5ad_forward'] for x in info),
                                    h5ad_cells=sum(x['h5ad_cells'] for x in info), h5ad_annotated=sum(x['h5ad_annotated'] for x in info),
                                    h5ad_cluster_agrees=sum(x['h5ad_cluster_agrees'] for x in info)),
               retention=dict(molecules_all_cells=sum(r['molecules_all_cells'] for r in summary_rows),
                              annotated=sum(r['molecules_annotated_cells'] for r in summary_rows),
                              annotated_with_depth=sum(r['molecules_annotated_with_depth'] for r in summary_rows)),
               calibration=dict(state_null_tests=int(len(null_p)), state_null_fraction_p_lt_0_05=float(np.mean(null_p < 0.05)) if len(null_p) else None,
                                cohort_null_replicates=200,
                                cohort_null_fraction_p_lt_0_05={f'{a}_vs_{b}': float(np.nanmean(null_cohort[:, n] < 0.05)) for n, (a, b) in enumerate(CONTRASTS)},
                                planted_state=dict(library=olibs[L0], state=planted_state, depth_share=float(share.min()), cells=8,
                                                   p=planted_state_p),
                                planted_aged_only_10_molecules_per_aged_unit={f'{a}_vs_{b}_p': plant_res[f'{a}_vs_{b}_p'] for a, b in CONTRASTS}),
               seconds=round(time.time() - started, 1))
(OUT / 'EXPRESSION_SUMMARY.json').write_text(json.dumps(summary, indent=1, default=str) + '\n')
print(json.dumps(summary, indent=1, default=str))
