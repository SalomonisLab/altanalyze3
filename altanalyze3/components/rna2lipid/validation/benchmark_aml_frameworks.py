"""Compare AML ridge and the lung sparse ElasticNet on identical CPTAC folds.

External studies never enter feature selection, scaling, or model selection.
Run with Python 3.11; source_root is the original Human-MS-impute workspace.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import pickle
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import spearmanr, rankdata
from sklearn.linear_model import ElasticNetCV, Ridge
from sklearn.model_selection import KFold
from threadpoolctl import threadpool_limits


def fit_target(ti, xs, ys, xt, order, prior, pool, full=False):
    """Lung candidate rule, including its original training sparsity score.

    Inner ElasticNetCV sees features selected on the outer training set; only
    the outer folds are used to report held-out performance.
    """
    y = ys[:, ti]
    best = None
    candidates = []
    with threadpool_limits(limits=1), warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for n in (25, 50, 100, 200):
            ix = order[ti, :n]
            m = ElasticNetCV(l1_ratio=[.7, .9, .95, 1.],
                             alphas=np.logspace(-3, 1, 30), cv=3,
                             max_iter=50000, n_jobs=1, random_state=1,
                             selection='cyclic').fit(xs[:, ix], y)
            fitted = m.predict(xs[:, ix])
            r2 = 1 - np.sum((y-fitted)**2) / max(np.sum((y-y.mean())**2), 1e-30)
            nnz = int(np.count_nonzero(m.coef_))
            score = r2 - .002 * nnz
            candidates.append((n, float(m.alpha_), float(m.l1_ratio_), nnz, float(score)))
            if best is None or score > best[0]:
                best = (score, ix, m)
        _, ix, m = best
        # Original ridge selection: prior genes ranked by training Spearman,
        # followed by expressed genes. All choices are in the outer train fold.
        pri, corr = prior
        chosen = sorted(pri, key=lambda j: -corr[j])[:120]
        taken = set(chosen)
        for j in pool[np.argsort(-corr[pool])]:
            if len(chosen) == 120:
                break
            if j not in taken:
                chosen.append(int(j)); taken.add(j)
        rm = Ridge(alpha=100).fit(xs[:, chosen], y)
        return ti, m.predict(xt[:, ix]), rm.predict(xt[:, chosen]), {
            'ix': ix, 'coef': m.coef_, 'intercept': float(m.intercept_),
            'top_n': len(ix), 'alpha': float(m.alpha_),
            'l1_ratio': float(m.l1_ratio_), 'nonzero': int(np.count_nonzero(m.coef_)),
            'candidates': candidates,
        }


def metrics(y, p):
    rows = []
    for i, t in enumerate(y.columns):
        a, b = y.iloc[:, i].values, p[:, i]
        sp = float(spearmanr(a, b).statistic) if np.std(b) > 1e-12 else np.nan
        rows.append({'target': t, 'spearman': sp,
                     'r2': 1-np.sum((a-b)**2)/np.sum((a-a.mean())**2),
                     'rmse': np.sqrt(np.mean((a-b)**2))})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source-root', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--jobs', type=int, default=6)
    args = ap.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location('original_aml_evaluation', args.source_root/'code/evaluate_imputation.py')
    original = importlib.util.module_from_spec(spec); spec.loader.exec_module(original)
    original.ROOT = str(args.source_root)
    x, y, pri = original.load('lipid')
    pd.DataFrame({'case': x.index, 'fold': -1}).to_csv(args.output/'cases.csv', index=False)
    x.to_csv(args.output/'training_rna_cp10k_log1p.csv.gz', compression='gzip')
    y.to_csv(args.output/'training_lipids_native.csv.gz', compression='gzip')
    print(f'{x.shape=} {y.shape=}', flush=True)
    gidx = {g: i for i, g in enumerate(x.columns)}
    prior = [[gidx[g] for g in pri[t] if g in gidx] for t in y.columns]
    xv, yv = x.values, y.values
    assert np.isfinite(xv).all() and np.isfinite(yv).all()
    # Retain original expressed-pool definition for an exact original baseline.
    pool = np.where(xv.mean(0) > np.quantile(xv.mean(0), .5))[0]
    splits = list(KFold(5, shuffle=True, random_state=0).split(x))
    predictions = {name: np.full(y.shape, np.nan) for name in ['sparse_elasticnet', 'original_ridge', 'mean']}
    fold_cases = []
    for fold, (tr, te) in enumerate(splits + [(np.arange(len(x)), np.arange(len(x)))]):
        print(f'Start fold {fold+1}/6 (6 = full training)', flush=True)
        mu, sd = xv[tr].mean(0), xv[tr].std(0); sd[sd < 1e-8] = 1
        ym, yd = yv[tr].mean(0), yv[tr].std(0); yd[yd < 1e-8] = 1
        xs, xt = (xv[tr]-mu)/sd, (xv[te]-mu)/sd
        ys = (yv[tr]-ym)/yd
        with threadpool_limits(limits=1):
            corr = np.nan_to_num((ys.T @ xs)/len(tr))
            order = np.argsort(-np.abs(corr), axis=1, kind='stable')[:, :200]
            sr = np.abs(original.col_spearman(yv[tr], xv[tr]))
        cache = args.output/f'fold_{fold}.pkl.gz'
        if cache.exists():
            with gzip.open(cache, 'rb') as f: results = pickle.load(f)
        else:
            results = Parallel(n_jobs=args.jobs, backend='loky', verbose=5)(
                delayed(fit_target)(ti, xs, ys, xt, order, (prior[ti], sr[ti]), pool)
                for ti in range(y.shape[1]))
            with gzip.open(cache, 'wb') as f: pickle.dump(results, f)
        if fold < 5:
            for ti, ep, rp, _ in results:
                predictions['sparse_elasticnet'][te, ti] = ep*yd[ti]+ym[ti]
                predictions['original_ridge'][te, ti] = rp*yd[ti]+ym[ti]
            predictions['mean'][te] = ym
            fold_cases.extend({'case': x.index[i], 'fold': fold} for i in te)
        else:
            results.sort(key=lambda r: r[0]); models = [r[3] for r in results]
            b = {'X_columns': x.columns.tolist(), 'Y_columns': y.columns.tolist(),
                 'mu': mu, 'sd': sd, 'sel_idx': [m['ix'] for m in models],
                 'coef': [m['coef']*yd[i] for i,m in enumerate(models)],
                 'intercept': np.array([m['intercept']*yd[i]+ym[i] for i,m in enumerate(models)]),
                 'metadata': {'modality': 'lipid', 'estimator': 'sparse per-lipid ElasticNetCV',
                              'normalization': 'CP10k+log1p; training z-score',
                              'n_train_cases': len(x), 'source': str(args.source_root),
                              'target_scale': 'log2(intensity+1), per-sample median centered separately within ion mode'}}
            with gzip.open(args.output/'sparse_elasticnet_aml_bundle.pkl.gz','wb') as f:pickle.dump(b,f)
            pd.DataFrame([{k:v for k,v in m.items() if k not in ['ix','coef','candidates']} | {'target':y.columns[i]} for i,m in enumerate(models)]).to_csv(args.output/'sparse_model_parameters.csv',index=False)
        print(f'Finished fold {fold+1}/6', flush=True)
    pd.DataFrame(fold_cases).to_csv(args.output/'cases.csv', index=False)
    summary = {'cases': len(x), 'genes': x.shape[1], 'lipids': y.shape[1],
               'split': '5-fold case-disjoint KFold(shuffle=True, random_state=0)',
               'primary_selection_metric': 'median per-lipid out-of-fold Spearman; constant predictions score zero',
               'limitations': ['Gene universe and ridge expressed-fill pool retained from original full cohort.',
                               'Hyperparameter menus were fixed before external evaluation; this is a development benchmark, not a fresh prospective CPTAC test.',
                               'Lung inner CV ranks genes on its outer training fold; outer test cases are untouched.']}
    source_files=[args.source_root/'rnaseq_counts/star_counts_lipid_matched.txt.gz',
                  args.source_root/'rnaseq_counts/star_tpm_lipid_matched.txt.gz',
                  args.source_root/'rnaseq_counts/gene_annotation.tsv',
                  args.source_root/'data/lipidomics_unique_matched.txt.gz',
                  args.source_root/'data/lipidomics_unique_annotation.tsv',
                  args.source_root/'metabolic_model/lipid_class_gene_prior.json',
                  args.source_root/'metabolic_model/lung_rna2lipid_genes.txt',
                  args.source_root/'code/evaluate_imputation.py']
    summary['training_source_files']=[{'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size} for f in source_files]
    tables = {}
    for name, p in predictions.items():
        pd.DataFrame(p,index=x.index,columns=y.columns).to_csv(args.output/f'{name}_oof.csv.gz',compression='gzip')
        tables[name] = metrics(y,p);tables[name].to_csv(args.output/f'{name}_metrics.csv',index=False)
        sp = tables[name].spearman.fillna(0)
        summary[name] = {'median_spearman_constant_zero':float(sp.median()),
                         'median_spearman_nonconstant':float(tables[name].spearman.median()),
                         'sp_gt_0p3':int((sp>.3).sum()),'sp_gt_0p5':int((sp>.5).sum()),
                         'constant_predictions':int(tables[name].spearman.isna().sum()),
                         'median_r2':float(tables[name].r2.median())}
    old = np.load(args.source_root/'reports/impute_lipid_oof.npy')
    summary['original_oof_max_absolute_difference'] = float(np.nanmax(np.abs(old-predictions['original_ridge'])))
    assert summary['original_oof_max_absolute_difference'] < 1e-6, 'Original OOF reproduction failed'
    candidates = ['sparse_elasticnet','original_ridge']
    summary['selected_model'] = max(candidates, key=lambda n: summary[n]['median_spearman_constant_zero'])
    delta=tables['sparse_elasticnet'].spearman.fillna(0)-tables['original_ridge'].spearman.fillna(0)
    pd.DataFrame({'target':y.columns,'sparse_minus_ridge_spearman':delta}).to_csv(args.output/'paired_model_comparison.csv',index=False)
    summary['sparse_wins_per_lipid'] = int((delta>0).sum())
    # Resample cases jointly across targets and models to retain lipid dependence.
    def median_sp(a,b):
        ar,br=rankdata(a,axis=0),rankdata(b,axis=0)
        ar-=ar.mean(0);br-=br.mean(0)
        den=np.sqrt(np.sum(ar**2,0)*np.sum(br**2,0))
        r=np.divide(np.sum(ar*br,0),den,out=np.zeros(a.shape[1]),where=den>1e-20)
        return float(np.median(r))
    rng=np.random.default_rng(0);boot=[]
    for _ in range(500):
        ix=rng.integers(0,len(x),len(x))
        boot.append(median_sp(yv[ix],predictions['sparse_elasticnet'][ix])-median_sp(yv[ix],predictions['original_ridge'][ix]))
    summary['paired_case_bootstrap_sparse_minus_ridge_median_spearman_95ci']=np.quantile(boot,[.025,.975]).tolist()
    summary['bootstrap_note']='500 paired case resamples of pooled out-of-fold predictions; preserves target dependence but does not refit training models'
    with gzip.open(args.output/'sparse_elasticnet_aml_bundle.pkl.gz','rb') as f: b=pickle.load(f)
    b['metadata']['per_target']={'heldout_spearman':dict(zip(y.columns,tables['sparse_elasticnet'].spearman.fillna(0)))}
    b['metadata']['heldout_median_spearman']=summary['sparse_elasticnet']['median_spearman_constant_zero']
    b['metadata']['selected_for_external_validation']=summary['selected_model']=='sparse_elasticnet'
    b['metadata']['target_scale']='log2(intensity+1), per-sample median centered separately within ion mode'
    with gzip.open(args.output/'sparse_elasticnet_aml_bundle.pkl.gz','wb') as f:pickle.dump(b,f)
    (args.output/'benchmark_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == '__main__':
    main()
