"""Reproducible, checkpointed ND20_167 HSC/MPP benchmark of the scALABLE interface.

Run from the repository with .venv/bin/python -m
altanalyze3.components.ambient_rna.evaluation.benchmark_nd20_167.
Original matrices are read only. Outputs live beside this script in benchmarking/.
"""
from pathlib import Path
import gc
import hashlib
import json
import time
import os
import sys
import shutil

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, io
from sklearn.neighbors import NearestNeighbors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from altanalyze3.components.cellHarmony import cellHarmony_lite as lite
from altanalyze3.components.cellHarmony import cellHarmony_differential as diff
from altanalyze3.components.ambient_rna import ambient_subtract as ambient

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'benchmarking' / 'ND20_167_20261002'
BASE = Path('/Volumes/salomonis2/Grimes/RNA/scRNA-Seq/10x-Genomics/X202SC21094303-Z02-F001')
REF = ROOT.parent / 'cellHarmony/flask/references/Human/BoneMarrow/Zhang-2024/Hs-MarrowAtlas-L3M.txt'
KEY = REF.stem
METHODS = ['uncorrected', 'python_auto', 'python_0.10', 'python_0.20', 'python_0.30',
           'SoupX_0.10', 'SoupX_0.15', 'SoupX_0.20', 'SoupX_0.25']

def total(x):
    return float(x.sum(dtype=np.float64))

def load_inputs(method):
    cache = OUT / 'inputs' / (method + '.h5ad')
    if cache.exists():
        return ad.read_h5ad(cache)
    captures = []
    for cap in ['HSC', 'MPP']:
        path = BASE / f'ND20_167_{cap}' / 'outs'
        if method == 'uncorrected':
            a = sc.read_10x_h5(path / 'filtered_feature_bc_matrix.h5')
        else:
            path = path / ('soupX-contamination-fraction-' + method.split('_')[1])
            genes = pd.read_csv(path / 'genes.tsv', sep='\t', header=None)
            bars = pd.read_csv(path / 'barcodes.tsv', sep='\t', header=None)[0].astype(str)
            matrix_path = path / 'matrix.mtx'
            if method == 'SoupX_0.25':
                # Stage a byte copy for the fast multithreaded Matrix Market reader.
                # A direct read from the mounted volume failed despite valid integer tokens.
                staged = OUT/'inputs'/f'{method}_{cap}_matrix.mtx'
                staged.parent.mkdir(parents=True,exist_ok=True)
                if not staged.exists():
                    partial = staged.with_suffix('.partial')
                    shutil.copyfile(matrix_path,partial)
                    partial.replace(staged)
                matrix_path = staged
            a = ad.AnnData(io.mmread(matrix_path).T.tocsr().astype(np.float32),
                           obs=pd.DataFrame(index=pd.Index(bars.to_numpy(),name=None)),
                           var=pd.DataFrame(index=pd.Index(genes[1].astype(str).to_numpy(),name=None)))
        a.var_names_make_unique()
        a.obs['barcode'] = a.obs_names.astype(str)
        a.obs['Library'] = cap
        a.obs['patient'] = 'ND20_167'
        a.obs_names = a.obs_names.astype(str) + '.' + cap
        captures.append(a)
        print('LOADED', method, cap, a.shape, total(a.X), flush=True)
    assert captures[0].var_names.equals(captures[1].var_names)
    a = ad.concat(captures, join='inner', merge='same')
    cache.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(cache)
    return a

def embed(a, seed):
    genes = pd.read_csv(REF, sep='\t', index_col=0).index.astype(str)
    genes = [g for g in genes if g in a.var_names]
    # Exact cellHarmony-lite umap_features='reference' pathway; no integration.
    b = ad.AnnData(a[:, genes].X.copy(), obs=a.obs[['Library', 'fixed_population']].copy())
    sc.pp.pca(b, n_comps=min(50, b.n_vars - 1, b.n_obs - 1), random_state=seed)
    sc.pp.neighbors(b, random_state=seed)
    sc.tl.umap(b, random_state=seed, init_pos='spectral')
    return b.obsm['X_umap'], b.obsm['X_pca'], len(genes)

def distances(a, xy, pcs, method, seed, states):
    rows = []
    global_scale = np.sqrt(np.mean(np.sum((xy - xy.mean(0))**2, axis=1)))
    for state in states:
        idx = np.flatnonzero(a.obs.fixed_population.to_numpy() == state)
        lib = a.obs.Library.to_numpy()[idx]
        x = xy[idx][lib == 'HSC']; y = xy[idx][lib == 'MPP']
        delta = x.mean(0) - y.mean(0)
        within = np.sqrt((np.mean(np.sum((x-x.mean(0))**2, axis=1)) +
                          np.mean(np.sum((y-y.mean(0))**2, axis=1))) / 2)
        p = pcs[idx]; px = p[lib == 'HSC']; py = p[lib == 'MPP']
        pw = np.sqrt((np.mean(np.sum((px-px.mean(0))**2, axis=1)) +
                      np.mean(np.sum((py-py.mean(0))**2, axis=1))) / 2)
        k = min(30, len(idx)-1)
        nn = NearestNeighbors(n_neighbors=k+1).fit(p).kneighbors(p, return_distance=False)[:, 1:]
        cross = float(np.mean(lib[nn] != lib[:, None]))
        expected = 2 * len(x)*len(y)/(len(idx)*(len(idx)-1))
        rows.append(dict(method=method, seed=seed, population=state, n_HSC=len(x), n_MPP=len(y),
                         delta_UMAP1=float(delta[0]), delta_UMAP2=float(delta[1]),
                         UMAP_distance=float(np.linalg.norm(delta)),
                         UMAP_distance_global_scaled=float(np.linalg.norm(delta)/global_scale),
                         UMAP_distance_within_scaled=float(np.linalg.norm(delta)/max(within, 1e-9)),
                         PCA_distance_within_scaled=float(np.linalg.norm(px.mean(0)-py.mean(0))/max(pw, 1e-9)),
                         PCA_neighbor_mixing=cross/expected))
    return rows

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    np.random.seed(0)
    raw = load_inputs('uncorrected')
    # Freeze QC from the original filtered counts, before seeing method results.
    totals = np.asarray(raw.X.sum(1)).ravel()
    detected = np.asarray((raw.X > 0).sum(1)).ravel()
    mito = np.asarray(raw.X[:, raw.var_names.str.startswith('MT-')].sum(1)).ravel()
    keep = (totals >= 500) & (detected >= 200) & (mito / np.maximum(totals,1) < .10)
    cohort = raw.obs_names[keep]
    pd.DataFrame({'cell':raw.obs_names, 'Library':raw.obs.Library.values,
                  'UMIs':totals, 'genes':detected, 'mito_fraction':mito/np.maximum(totals,1),
                  'included':keep}).to_csv(OUT/'qc.tsv',sep='\t',index=False)
    metrics = []; losses = []; embeds = []; assignments = []
    fixed = None; predominant = None; shared = None
    for method in METHODS:
        dest = OUT / method
        dest.mkdir(exist_ok=True)
        processed = dest/'aligned_fixed_cohort.h5ad'
        if processed.exists():
            a = ad.read_h5ad(processed)
            counts = pd.read_csv(dest/'count_loss.tsv',sep='\t')
        else:
            inp = raw.copy() if method.startswith('python') or method == 'uncorrected' else load_inputs(method)
            assert inp.obs_names.is_unique and inp.var_names.is_unique
            assert set(inp.obs_names) == set(raw.obs_names), 'Barcode mismatch'
            if method.startswith('SoupX') and not set(inp.var_names) == set(raw.var_names):
                # R make.unique uses '.1', whereas AnnData uses '-1'. Confirm every
                # positional mismatch is exactly a duplicate-symbol suffix difference;
                # never strip dots from legitimate gene symbols such as RP11-34P13.7.
                assert inp.n_vars == raw.n_vars
                mapping = []
                seen = set()
                for canonical, supplied in zip(raw.var_names, inp.var_names):
                    if canonical != supplied:
                        stem, suffix = canonical.rsplit('-',1)
                        assert suffix.isdigit() and stem in seen and supplied == stem+'.'+suffix, (canonical,supplied)
                        mapping.append(dict(SoupX_symbol=supplied, canonical_symbol=canonical))
                    seen.add(canonical)
                pd.DataFrame(mapping).to_csv(dest/'duplicate_symbol_mapping.tsv',sep='\t',index=False)
                inp.var_names = raw.var_names.copy()
            assert set(inp.var_names) == set(raw.var_names), 'Feature mismatch'
            inp = inp[raw.obs_names, raw.var_names].copy()
            cutoff = ('auto' if method == 'python_auto' else method.split('_')[1]) if method.startswith('python') else None
            timings = []
            original = ambient.process_anndata
            def timed(*args, **kwargs):
                start = time.perf_counter()
                result = original(*args, **kwargs)
                timings.append(time.perf_counter()-start)
                return result
            ambient.process_anndata = timed
            start = time.perf_counter()
            try:
                _, aligned = lite.combine_and_align_h5([], str(REF), adata=inp,
                    output_dir=str(dest), min_genes=0, min_cells=None, min_counts=0,
                    mit_percent=101, generate_umap=False, save_adata=False,
                    unsupervised_cluster=False, alignment_mode='cosine', min_alignment_score=None,
                    ambient_correct_cutoff=cutoff, ambient_memory_efficient=True, return_adata=True)
            finally:
                ambient.process_anndata = original
            elapsed = time.perf_counter()-start
            assert set(aligned.obs_names) == set(raw.obs_names)
            aligned = aligned[raw.obs_names].copy()
            # Verify scale handling; X must be log(CP10k+1), counts must be corrected counts.
            assert 'counts' in aligned.layers, 'Counts were misclassified by expression-scale detection'
            assert np.isfinite(aligned.X.data).all()
            test = aligned.layers['counts'][:10].copy()
            check = ad.AnnData(test); lite.normalize_adata(check)
            assert np.max(np.abs((check.X-aligned.X[:10]).data), initial=0) < 1e-4
            records = []
            for cap in ['HSC','MPP']:
                rows = (raw.obs.Library == cap).to_numpy()
                qrows = rows & keep
                original_total = total(raw.X[rows]); corrected_total = total(aligned.layers['counts'][rows])
                original_q = total(raw.X[qrows]); corrected_q = total(aligned.layers['counts'][qrows])
                per_cell_raw = np.asarray(raw.X[qrows].sum(1)).ravel()
                per_cell_cor = np.asarray(aligned.layers['counts'][qrows].sum(1)).ravel()
                records.append(dict(method=method, capture=cap, filtered_cells=int(rows.sum()),
                    QC_cells=int(qrows.sum()), original_UMIs=original_total, retained_UMIs=corrected_total,
                    loss_percent=100*(1-corrected_total/original_total), QC_loss_percent=100*(1-corrected_q/original_q),
                    median_cell_loss_percent=float(np.median(100*(1-per_cell_cor/per_cell_raw))),
                    correction_seconds=sum(timings) if timings else np.nan, alignment_seconds=elapsed))
            counts = pd.DataFrame(records); counts.to_csv(dest/'count_loss.tsv',sep='\t',index=False)
            a = aligned[cohort].copy()
            if method == 'uncorrected':
                score = pd.read_csv(dest/'cellHarmony_lite_assignments.txt',sep='\t').set_index('CellBarcode').AlignmentScore
                a = a[score.reindex(a.obs_names).to_numpy() >= .1].copy()
            a.raw = None
            for layer in list(a.layers):
                if layer != 'counts': del a.layers[layer]
            a.write_h5ad(processed)
            del inp, aligned; gc.collect()
        losses.extend(counts.to_dict('records'))
        if method == 'uncorrected':
            cohort = a.obs_names.copy()
            fixed = a.obs[KEY].astype(str).copy()
            tab = pd.crosstab(fixed, a.obs.Library).reindex(columns=['HSC','MPP'],fill_value=0)
            tab.to_csv(OUT/'population_counts.tsv',sep='\t')
            shared = tab.index[(tab.min(axis=1) >= 20) & (tab.index != 'Unassigned')].tolist()
            predominant = tab.index[(tab.min(axis=1) >= 100) & (tab.index != 'Unassigned')].tolist()
            assert predominant
            fixed.to_csv(OUT/'fixed_population_assignments.tsv',sep='\t')
            print('PREDOMINANT', predominant, 'SHARED', shared, flush=True)
        a.obs['fixed_population'] = fixed.reindex(a.obs_names).values
        assignments.append(dict(method=method, population_reassignment_percent=100*float(np.mean(a.obs[KEY].astype(str).values != a.obs.fixed_population.values))))
        de_file = dest/'DE_summary.tsv'
        if not de_file.exists():
            # Use the actual cellHarmony differential interface, including its count-based FC.
            b = a[a.obs.fixed_population.isin(shared)].copy()
            store = diff.run_de_for_comparisons(b, 'fixed_population','Library','MPP','HSC',
                        method='wilcoxon',alpha=.05,fc_thresh=1.2,min_cells_per_group=20,use_rawp=False)
            summary_table = store['summary_per_population']
            # cellHarmony currently omits a summary row when a tested state has zero DEGs.
            # Recover that zero from its complete statistics rather than calling it untested.
            missing = []
            for state in shared:
                if state not in set(summary_table.population):
                    s = store['per_population_deg'][state]
                    n_deg = int(((s.fdr < .05) & (s.log2fc.abs() > np.log2(1.2))).sum())
                    assert n_deg == 0
                    labels = b.obs.loc[b.obs.fixed_population == state, 'Library']
                    missing.append(dict(population=state,n_case=int((labels=='MPP').sum()),
                                        n_control=int((labels=='HSC').sum()),num_DEG=0,tested_genes=len(s)))
            summary_table = pd.concat([summary_table,pd.DataFrame(missing)],ignore_index=True)
            summary_table.to_csv(de_file,sep='\t',index=False)
            store['detailed_deg'].to_csv(dest/'DEGs.tsv.gz',sep='\t',index=False)
            # Save continuous statistics for every tested feature and fold threshold sensitivity.
            all_stats = pd.concat([v.reset_index().assign(population=k) for k,v in store['per_population_deg'].items()], ignore_index=True)
            all_stats.to_csv(dest/'all_gene_statistics.tsv.gz',sep='\t',index=False)
            del b, store, all_stats; gc.collect()
        summary = pd.read_csv(de_file,sep='\t')
        pred = summary[summary.population.isin(predominant)]
        stats = pd.read_csv(dest/'all_gene_statistics.tsv.gz',sep='\t')
        primary = stats[stats.population.isin(predominant)]
        record = dict(method=method, predominant_DEG_gene_state_pairs=int(pred.num_DEG.sum()),
            shared_DEG_gene_state_pairs=int(summary.num_DEG.sum()),
            predominant_unique_DEGs=int(primary.loc[(primary.fdr<.05)&(primary.log2fc.abs()>np.log2(1.2)),'gene'].nunique()),
            median_abs_log2FC=float(primary.log2fc.abs().median()))
        for fold in [1.1,1.5,2.0]:
            record['DEG_FC_'+str(fold)] = int(((primary.fdr<.05)&(primary.log2fc.abs()>np.log2(fold))).sum())
        metrics.append(record)
        del stats, primary
        for seed in [0,1,2]:
            emb_file = dest/f'embedding_seed{seed}.npz'
            if emb_file.exists():
                e = np.load(emb_file); xy=e['umap']; pcs=e['pca']
            else:
                xy, pcs, n_genes = embed(a,seed)
                np.savez_compressed(emb_file,umap=xy,pca=pcs,cell=a.obs_names.to_numpy(dtype=str))
            embeds.extend(distances(a,xy,pcs,method,seed,shared))
            if seed == 0:
                pd.DataFrame(dict(cell=a.obs_names, capture=a.obs.Library, population=a.obs.fixed_population,
                                 UMAP1=xy[:,0],UMAP2=xy[:,1])).to_csv(dest/'umap_coordinates.tsv',sep='\t',index=False)
        print('DONE', record, flush=True)
        del a; gc.collect()
        pd.DataFrame(metrics).to_csv(OUT/'method_metrics.tsv',sep='\t',index=False)
        pd.DataFrame(losses).to_csv(OUT/'transcript_loss.tsv',sep='\t',index=False)
        pd.DataFrame(embeds).to_csv(OUT/'embedding_distances.tsv',sep='\t',index=False)
        pd.DataFrame(assignments).to_csv(OUT/'assignment_stability.tsv',sep='\t',index=False)
    print('ALL METHODS COMPLETE',flush=True)

if __name__ == '__main__':
    main()
