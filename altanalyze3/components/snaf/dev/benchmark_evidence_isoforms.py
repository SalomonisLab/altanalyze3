"""Reproducible gene-held-out evaluation of sample-supported isoform inference.

Long-read chains/proteins are read only by this evaluator, never by the generator.
The primary structural label is the previously associated long-read isoform. Also
report structural agreement with ANY target-containing long-read chain in the gene.
This is a pooled-cohort benchmark, not matched-sample validation or molecular phasing.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys
import time
from collections import defaultdict

os.environ.setdefault('SNAF_OFFLINE', '1')
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np
import pysam
from Bio.Align import PairwiseAligner
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from altanalyze3.components.snaf.surface import predict_isoform as P
from altanalyze3.components.snaf.surface import evidence_isoform as E

EXAMPLES = {
    'ENSG00000103257:E3.1-I3.1_87843819': ('SLC7A5', 'chr16', (87843819, 87851724)),
    'ENSG00000197965:I6.1_167744987-E10.3': ('MPZL1', 'chr1', (167744987, 167765583)),
    'ENSG00000140090:E3.2-I3.1_92324060': ('SLC24A4', 'chr14', (92323960, 92324060)),
}


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def chain(text):
    return sorted(tuple(map(int, x.split('-'))) for x in text.split(',') if x)


def write_tsv(path, records):
    if not records:
        return
    with open(path, 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(records[0]), delimiter='\t')
        w.writeheader()
        w.writerows(records)


def fold(gene, n=5):
    return int(hashlib.sha256(gene.encode()).hexdigest()[:12], 16) % n


def structure_scores(predicted, truth):
    p, t = set(predicted), set(truth)
    shared = len(p & t)
    f1 = 2 * shared / (len(p) + len(t)) if p or t else 1.0
    return f1, float(tuple(predicted) == tuple(truth))


ALIGNER = PairwiseAligner(mode='global', match_score=1, mismatch_score=-1,
                          open_gap_score=-2, extend_gap_score=-0.5)


def identity(a, b):
    if not a or not b:
        return 0.0
    if a == b:
        return 1.0
    al = ALIGNER.align(a, b)[0]
    matched = 0
    for (a0, a1), (b0, b1) in zip(*al.aligned):
        matched += sum(x == y for x, y in zip(a[a0:a1], b[b0:b1]))
    return matched / max(len(a), len(b))


def measure(pred, target, truth, truth_protein, epitope, lr_chains, memo):
    if pred is None:
        return dict(protein_identity=0., intron_f1=0., exact_chain=0.,
                    any_lr_exact=0., epitope=0., prediction=0., nmd=0.)
    pi = E.introns(pred.chain)
    key = (pred.protein, truth_protein)
    if key not in memo:
        memo[key] = identity(*key)
    f1, exact = structure_scores(pi, truth)
    return dict(protein_identity=memo[key], intron_f1=f1, exact_chain=exact,
                any_lr_exact=float(pi in lr_chains),
                epitope=float(bool(epitope) and epitope in pred.protein),
                prediction=1., nmd=float(pred.nmd))


def model_spec(name):
    if name == 'pairwise':
        return E.LinearQueryRanker(alpha=10)
    if name == 'ridge':
        return make_pipeline(StandardScaler(), Ridge(alpha=100))
    if name == 'trees':
        return ExtraTreesRegressor(n_estimators=160, min_samples_leaf=12,
                                   max_features=0.8, random_state=193, n_jobs=4)
    return HistGradientBoostingRegressor(max_iter=120, max_leaf_nodes=7,
                                        l2_regularization=10, random_state=193)


def train_model(name, data, indices, objective='joint'):
    x, y, weights = [], [], []
    for i in indices:
        d = data[i]
        if not d.get('truth_eligible', True):
            continue
        start = len(x)
        for c, s in zip(d['candidates'], d['metrics']):
            x.append(c.features)
            y.append(s['intron_f1'] if objective == 'structure' else
                     (s['intron_f1'] + s['protein_identity']) / 2)
            weights.append(1 / max(1, len(d['candidates'])))
        if name == 'pairwise' and len(x) > start:
            x[start:] = (np.asarray(x[start:]) - np.mean(x[start:], axis=0)).tolist()
            y[start:] = (np.asarray(y[start:]) - np.mean(y[start:])).tolist()
    model = model_spec(name)
    if name == 'ridge':
        model.fit(x, y, ridge__sample_weight=weights)
    else:
        model.fit(x, y, sample_weight=weights)
    return model


def choose(model, d):
    if not d['candidates']:
        return -1
    return int(np.argmax(model.predict([c.features for c in d['candidates']])))


def summary(records):
    keys = ['protein_identity', 'intron_f1', 'exact_chain', 'any_lr_exact',
            'epitope', 'prediction', 'nmd']
    return {k: float(np.mean([r[k] for r in records])) for k in keys}


def generate(args):
    root = Path(args.project)
    out = Path(args.out)
    truth_rows = rows(root / 'results/03_truth.tsv')
    old_rows = {r['NeoJunction']: r for r in rows(root / 'results/predicted_final_multi/predicted_isoforms.tsv')}
    samples, counts = pickle.load(open(args.evidence, 'rb'))
    models = pickle.load(open(root / 'results/index/ens91_models.pkl', 'rb'))
    annotation = E.load_exon_annotation(args.exon_annotation, set(models)) if args.exon_annotation else {}
    lr = pickle.load(open(root / 'results/index/longread_chains.pkl', 'rb'))
    fa = pysam.FastaFile(args.genome)
    all_data = []
    start = time.time()
    for i, tr in enumerate(truth_rows):
        if args.limit and i >= args.limit:
            break
        uid, gene = tr['uid'], tr['gene']
        target = tuple(sorted((int(tr['junction_start']), int(tr['junction_end']))))
        # Trans-spliced calls cannot be represented by one genomic interval.
        chrom = tr['chrom'] if tr['chrom'] in fa.references else tr['chrom'].removeprefix('chr')
        valid = uid.count('ENSG') == 1 and target[0] > 0 and chrom in fa.references
        gm = models.get(gene, {})
        ev = E.Evidence({j: v['counts'] for j, v in counts.get(gene, {}).items()},
                        min_reads=args.min_reads)
        fetch = lambda s, e: fa.fetch(chrom, s - 1, e)
        truth_exons = chain(tr['lr_chain'])
        truth_source = 'longread'
        truth_chrom = tr['lr_chrom'].removeprefix('chr')
        # Some historical "long-read" associations actually name an ENST reference
        # transcript. Recover its annotated chain and report its provenance separately.
        if not truth_exons and tr['iso_id'] in gm:
            truth_exons = gm[tr['iso_id']][1]
            truth_source = 'reference_annotation'
            truth_chrom = chrom
        ti = E.introns(truth_exons)
        eligible = valid and E.contains(truth_exons, target) and chrom == truth_chrom
        reason = ('unsupported_fusion' if uid.count('ENSG') != 1 else
                  'missing_truth_chain' if not truth_exons else
                  'chromosome_mismatch' if chrom != truth_chrom else
                  'target_absent_from_truth_chain' if not E.contains(truth_exons, target)
                  else 'verified_' + truth_source)
        refchains = {E.introns(ch) for _, ch in gm.values()}
        # 03_truth.truth_prot is a segment, not the complete long-read protein. Use
        # neoisoform_protein_seq from the actual association table for full-length scoring.
        truth_protein = args.full_truth[uid]['neoisoform_protein_seq']
        epitope = args.full_truth[uid]['novel_insert_peptide']
        d = dict(uid=uid, gene=gene, symbol=tr['symbol'], target=target,
                 truth_eligible=eligible, truth_audit=reason,
                 truth_source=truth_source,
                 truth_chain=ti, truth_protein=truth_protein, epitope=epitope,
                 truth_kind='known_chain' if ti in refchains else 'novel_chain',
                 junction_kind='annotated' if target in {j for c in refchains for j in c} else 'novel',
                 target_samples=ev.recurrence((target,)), candidates=[], metrics=[], baseline={})
        memo = {}
        def metric(p):
            return measure(p, target, ti, truth_protein, epitope, lr.get(gene, {}), memo)
        baseline_j = [tuple(sorted((int(t['junction_start']), int(t['junction_end']))))
                      for t in truth_rows if t['gene'] == gene and t['uid'] != uid]
        for label, others in [('legacy_single', []), ('legacy_pair', baseline_j)]:
            pred = P.predict_isoform(gm, *target, fetch, co_junctions=others) if valid else None
            d['baseline'][label] = metric(pred)
        candidate_union = {}
        for edits in (1, 2, 3, 4):
            candidates = E.generate_candidates(gm, target, fetch, ev, max_edits=edits,
                                               junction_label=uid,
                                               exon_annotation=annotation.get(gene)) if valid else []
            clean = [c.prediction for c in candidates if not c.prediction.nmd]
            best = max(clean, key=lambda p: len(p.protein), default=None)
            d['baseline']['local%d_longest' % edits] = metric(best)
            ranked = E.rank_candidates(candidates)
            d['baseline']['local%d_evidence' % edits] = metric(ranked[0][1].prediction if ranked else None)
            ranked_clean = E.rank_candidates([c for c in candidates if not c.prediction.nmd])
            d['baseline']['local%d_evidence_non_nmd' % edits] = metric(
                ranked_clean[0][1].prediction if ranked_clean else None)
            for c in candidates:
                key = (tuple(c.prediction.chain), c.prediction.cds_start, c.prediction.protein)
                candidate_union.setdefault(key, c)
        d['candidates'] = list(candidate_union.values())
        # A larger beam must not evict the better two-event hypotheses. Recompute
        # query-relative length after merging candidates from all search depths.
        longest = max((len(c.prediction.protein) for c in d['candidates']), default=1)
        for c in d['candidates']:
            f = list(c.features)
            f[1] = len(c.prediction.protein) / longest
            c.features = tuple(f)
        d['metrics'] = [metric(c.prediction) for c in d['candidates']]
        ranked = E.rank_candidates(d['candidates'])
        d['baseline']['union_evidence'] = metric(ranked[0][1].prediction if ranked else None)
        d['empty'] = metric(None)
        all_data.append(d)
        if i % 5 == 0:
            print(i + 1, '/', len(truth_rows), uid, 'candidates', len(d['candidates']),
                  'seconds', round(time.time() - start), flush=True)
    with open(out / 'candidates.pkl', 'wb') as fh:
        pickle.dump(all_data, fh)
    fa.close()
    return all_data


def evaluate(args, data):
    out = Path(args.out)
    records, choices = [], []
    names = ['ridge', 'trees', 'boost', 'pairwise']
    for d in data:
        for method, scores in d['baseline'].items():
            records.append(dict(uid=d['uid'], gene=d['gene'], truth_kind=d['truth_kind'],
                                junction_kind=d['junction_kind'], method=method, **scores))
        if d['metrics']:
            oracle = max(d['metrics'], key=lambda s: (s['protein_identity'] + s['intron_f1']) / 2)
        else:
            oracle = d['empty']
        records.append(dict(uid=d['uid'], gene=d['gene'], truth_kind=d['truth_kind'],
                            junction_kind=d['junction_kind'], method='oracle_joint', **oracle))
    for f in range(5):
        train = [i for i, d in enumerate(data) if fold(d['gene']) != f and d['candidates']
                 and d.get('truth_eligible', True)]
        test = [i for i, d in enumerate(data) if fold(d['gene']) == f]
        assert {data[i]['gene'] for i in train}.isdisjoint({data[i]['gene'] for i in test})
        # One inner gene split selects model/objective without inspecting outer test labels.
        inner_val = [i for i in train if fold(data[i]['gene'], 3) == 0]
        inner_train = [i for i in train if i not in inner_val]
        selection = []
        for name in names:
            for objective in ('joint', 'structure'):
                m = train_model(name, data, inner_train, objective)
                scores = [data[i]['metrics'][choose(m, data[i])] for i in inner_val]
                score = float(np.mean([(s['protein_identity'] + s['intron_f1']) / 2 for s in scores]))
                selection.append((score, name, objective))
        _, selected, selected_objective = max(selection)
        choices.append({'fold': f, 'model': selected, 'objective': selected_objective,
                        'inner_scores': selection, 'test_genes': sorted({data[i]['gene'] for i in test})})
        print('fold', f, 'chosen', selected, selected_objective, flush=True)
        for name, objective in [(n, 'joint') for n in names] + [(selected, selected_objective)]:
            method = 'learned_' + name
            m = train_model(name, data, train, objective)
            is_selected = name == selected and objective == selected_objective
            if is_selected:
                with open(out / ('fold%d_ranker.pkl' % f), 'wb') as fh:
                    pickle.dump({'model': m, 'feature_names': E.FEATURE_NAMES,
                                 'training_genes': sorted({data[i]['gene'] for i in train})}, fh)
            for i in test:
                d = data[i]
                k = choose(m, d)
                scores = d['metrics'][k] if k >= 0 else d['empty']
                if objective == 'joint':
                    row = dict(uid=d['uid'], gene=d['gene'], truth_kind=d['truth_kind'],
                               junction_kind=d['junction_kind'], method=method, **scores)
                    if not any(r['uid'] == d['uid'] and r['method'] == method for r in records):
                        records.append(row)
                if is_selected:
                    row = dict(uid=d['uid'], gene=d['gene'], truth_kind=d['truth_kind'],
                               junction_kind=d['junction_kind'], method='nested_learned', **scores)
                    if not any(r['uid'] == d['uid'] and r['method'] == 'nested_learned' for r in records):
                        records.append(row)
                        d['selected'] = k
    eligible_uids = {d['uid'] for d in data if d.get('truth_eligible', True)}
    truth_sources = {d['uid']: d.get('truth_source', 'unknown') for d in data}
    for r in records:
        r['truth_verified'] = int(r['uid'] in eligible_uids)
        r['truth_source'] = truth_sources[r['uid']]
    write_tsv(out / 'metrics_per_junction.tsv', records)
    write_tsv(out / 'truth_audit.tsv', [{'junction_id': d['uid'], 'gene': d['gene'],
               'status': d.get('truth_audit', 'not_audited'),
               'used_as_supervised_truth': d.get('truth_eligible', True)} for d in data])
    summaries = []
    for method in sorted({r['method'] for r in records}):
        for subset in ('all', 'verified_truth', 'verified_longread', 'reference_annotation',
                       'known_chain', 'novel_chain',
                       'annotated_junction', 'novel_junction'):
            rr = [r for r in records if r['method'] == method and
                  (subset == 'all' or (subset == 'verified_truth' and r['truth_verified']) or
                   (subset == 'verified_longread' and r['truth_verified'] and r['truth_source'] == 'longread') or
                   (subset == 'reference_annotation' and r['truth_source'] == 'reference_annotation') or
                   r['truth_kind'] == subset or
                   r['junction_kind'] + '_junction' == subset)]
            if rr:
                summaries.append(dict(method=method, subset=subset, n=len(rr), **summary(rr)))
    write_tsv(out / 'summary.tsv', summaries)
    (out / 'folds.json').write_text(json.dumps(choices, indent=2))
    # Gene bootstrap uncertainty: resample genes, retaining their junction clusters.
    base = {r['uid']: r for r in records if r['method'] == 'legacy_pair'}
    selected = {r['uid']: r for r in records if r['method'] == 'nested_learned'}
    genes = sorted({d['gene'] for d in data})
    rng = np.random.default_rng(193)
    intervals = {}
    for metric in ('protein_identity', 'intron_f1', 'exact_chain', 'epitope'):
        by_gene = {g: [selected[d['uid']][metric] - base[d['uid']][metric]
                       for d in data if d['gene'] == g] for g in genes}
        means = [np.mean([v for g in rng.choice(genes, len(genes), replace=True) for v in by_gene[g]])
                 for _ in range(2000)]
        intervals[metric] = dict(delta=float(np.mean([v for vals in by_gene.values() for v in vals])),
                                ci95=np.quantile(means, [.025, .975]).tolist())
    (out / 'uncertainty.json').write_text(json.dumps(intervals, indent=2))
    predictions = []
    for d in data:
        k = d['selected']
        p = d['candidates'][k].prediction if k >= 0 else None
        row = {'junction_id': d['uid'], 'gene': d['gene'], 'symbol': d['symbol'],
               'truth_kind': d['truth_kind'], 'target_samples': d['target_samples'],
               'candidate_count': len(d['candidates'])}
        row.update({c: '' for c in P.PREDICTION_COLUMNS} if p is None else P.prediction_row(p, d['uid']))
        predictions.append(row)
    write_tsv(out / 'heldout_predictions.tsv', predictions)
    # A deployable model trained after evaluation, using the most frequently inner-selected spec.
    from collections import Counter
    spec = Counter((c['model'], c['objective']) for c in choices).most_common(1)[0][0]
    fitting = [i for i, d in enumerate(data) if d['candidates'] and d.get('truth_eligible', True)]
    fitted = train_model(*spec[:1], data, fitting, objective=spec[1])
    with open(out / 'ranker.pkl', 'wb') as fh:
        pickle.dump({'model': fitted, 'feature_names': E.FEATURE_NAMES, 'spec': spec,
                     'training_genes': sorted({data[i]['gene'] for i in fitting}),
                     'min_reads': args.min_reads}, fh)
    for r in summaries:
        if r['subset'] == 'all':
            print(r['method'], 'identity', round(r['protein_identity'], 4), 'F1', round(r['intron_f1'], 4),
                  'exact', round(r['exact_chain'], 4), 'epitope', round(r['epitope'], 4), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project', required=True)
    ap.add_argument('--evidence', required=True)
    ap.add_argument('--genome', required=True)
    ap.add_argument('--associations', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--min-reads', type=float, default=3)
    ap.add_argument('--exon-annotation')
    ap.add_argument('--reuse', action='store_true')
    ap.add_argument('--limit', type=int, default=0)
    args = ap.parse_args()
    Path(args.out).mkdir(parents=True, exist_ok=True)
    manifest = vars(args).copy()
    manifest['feature_names'] = E.FEATURE_NAMES
    manifest['code_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in (Path(__file__), Path(E.__file__), Path(P.__file__))}
    args.full_truth = {r['NeoJunction']: r for r in rows(args.associations)}
    (Path(args.out) / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    data = pickle.load(open(Path(args.out) / 'candidates.pkl', 'rb')) if args.reuse else generate(args)
    if args.limit:
        print('Pilot generation complete; CV requires the full benchmark.', flush=True)
        return
    evaluate(args, data)


if __name__ == '__main__':
    main()
