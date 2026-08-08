"""
Prove the three memory changes reproduce the labels, and measure what they cost.

The changes:
  A  metrics._shared_index_csr    shares src's index arrays instead of copying
                                  the whole matrix, at the t-test squaring and
                                  the tf-idf nonzero indicator
  B  metrics.check_filter_single_cluster  returns the object itself when no
                                  cluster holds one cell, instead of a view that
                                  rebuilds X on every attribute access
  C  cli.py                       frees layers that --layer does not name

C is measured separately; it does nothing on an input with no layers, which is
the case for both files here.

A and B are compared by restoring the OLD code with a monkeypatch and running
the same pipeline twice in the same interpreter. Peak memory comes from
tracemalloc, which counts python allocations only, so treat it as the delta
between the two runs rather than as the process high-water mark.

Writes /Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/
sctriangulate/benchmarks/three_changes_validate.json
"""

import os
import sys
import gc
import json
import time
import tracemalloc

import numpy as np
import matplotlib
matplotlib.use('Agg')
import scanpy as sc
from scipy.sparse import csr_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))

from altanalyze3.components.sctriangulate import ScTriangulate        # noqa: E402
from altanalyze3.components.sctriangulate import metrics as M         # noqa: E402
from altanalyze3.components.sctriangulate import main_class as MC     # noqa: E402

NEW_SHARED = M._shared_index_csr
NEW_FILTER = M.check_filter_single_cluster


def old_shared_index_csr(src, data):
    """What the two sites did before: copy everything, then replace the data."""
    out = src.copy()
    out.data = data
    return out


def old_check_filter_single_cluster(adata, key):
    """Upstream: always a view, even when nothing is excluded."""
    vc = adata.obs[key].value_counts()
    exclude_clusters = vc.loc[vc == 1].index
    truth = np.logical_not(adata.obs[key].isin(exclude_clusters).values)
    return adata[truth, :]


def install(which):
    if which == 'old':
        M._shared_index_csr = old_shared_index_csr
        M.check_filter_single_cluster = old_check_filter_single_cluster
        MC.check_filter_single_cluster = old_check_filter_single_cluster
    else:
        M._shared_index_csr = NEW_SHARED
        M.check_filter_single_cluster = NEW_FILTER
        MC.check_filter_single_cluster = NEW_FILTER


def run(h5ad, keys, outdir, which):
    install(which)
    M.RUN_ENRICHMENT = False
    adata = sc.read(h5ad)
    gc.collect()
    tracemalloc.start()
    t0 = time.perf_counter()
    s = ScTriangulate(dir=outdir, adata=adata, query=keys, predict_doublet=False)
    s.run_enrichment = False
    # one process: tracemalloc cannot see into workers
    s.lazy_run(compute_metrics_parallel=False, compute_shapley_parallel=False)
    secs = time.perf_counter() - t0
    _cur, pk = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    obs = s.adata.obs.copy()
    del s, adata
    gc.collect()
    return obs, pk / 1e6, secs


def compare(h5ad, keys, label, out):
    print('\n=== {} ==='.format(label))
    old_obs, old_pk, old_t = run(h5ad, keys, os.path.join(HERE, 'out_3chg_old'), 'old')
    print('  old code: peak {:8.1f} MB  {:6.1f} s'.format(old_pk, old_t))
    new_obs, new_pk, new_t = run(h5ad, keys, os.path.join(HERE, 'out_3chg_new'), 'new')
    print('  new code: peak {:8.1f} MB  {:6.1f} s'.format(new_pk, new_t))

    row = {'dataset': label, 'cells': int(old_obs.shape[0]),
           'old_peak_MB': old_pk, 'new_peak_MB': new_pk,
           'peak_ratio': old_pk / new_pk if new_pk else float('nan'),
           'old_seconds': old_t, 'new_seconds': new_t,
           'seconds_ratio': old_t / new_t if new_t else float('nan')}
    for col in ['raw', 'pruned', 'final_annotation', 'confidence']:
        if col not in old_obs.columns:
            continue
        a = old_obs[col].astype(str).values
        b = new_obs.loc[old_obs.index, col].astype(str).values
        n_same = int((a == b).sum())
        row['identical_' + col] = n_same == len(a)
        row['agree_' + col] = '{} of {}'.format(n_same, len(a))
        print('  {:<18s} {} of {} identical'.format(col, n_same, len(a)))
    out.append(row)


def main():
    out = []
    compare(os.path.join(HERE, 'data', 'demo_pbmc3k.h5ad'),
            ['sctri_rna_leiden_1', 'sctri_rna_leiden_2', 'sctri_rna_leiden_3'],
            'demo_2700x32738', out)

    sim = os.path.join(HERE, 'data', 'sim_8k_20k.h5ad')
    if os.path.exists(sim):
        a = sc.read(sim, backed='r')
        simkeys = [c for c in a.obs.columns if c.startswith('sctri_')]
        del a
        compare(sim, simkeys, 'sim_8000x20000', out)

    path = os.path.join(HERE, 'three_changes_validate.json')
    json.dump(out, open(path, 'w'), indent=1)
    print('\nwrote {}'.format(path))


if __name__ == '__main__':
    main()
