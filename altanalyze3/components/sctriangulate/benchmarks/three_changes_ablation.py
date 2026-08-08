"""
Which of the two changes moved the time, and is the move real?

three_changes_validate.py ran old-then-new and reported the new code slower.
Two explanations fit: one of the changes really costs time, or the second run in
an interpreter is slower for reasons that have nothing to do with the change.
This runs all four combinations of A and B, twice each, in both orders, with
tracemalloc off so its overhead cannot confound the timing.

  A = _shared_index_csr shares the index arrays
  B = check_filter_single_cluster returns the object when nothing is dropped

Writes /Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/
sctriangulate/benchmarks/three_changes_ablation.json
"""

import os
import sys
import gc
import json
import time
import itertools

import numpy as np
import matplotlib
matplotlib.use('Agg')
import scanpy as sc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))

from altanalyze3.components.sctriangulate import ScTriangulate        # noqa: E402
from altanalyze3.components.sctriangulate import metrics as M         # noqa: E402
from altanalyze3.components.sctriangulate import main_class as MC     # noqa: E402

NEW_SHARED = M._shared_index_csr
NEW_FILTER = M.check_filter_single_cluster
DEMO = os.path.join(HERE, 'data', 'demo_pbmc3k.h5ad')
KEYS = ['sctri_rna_leiden_1', 'sctri_rna_leiden_2', 'sctri_rna_leiden_3']


def old_shared_index_csr(src, data):
    out = src.copy()
    out.data = data
    return out


def old_check_filter_single_cluster(adata, key):
    vc = adata.obs[key].value_counts()
    exclude_clusters = vc.loc[vc == 1].index
    truth = np.logical_not(adata.obs[key].isin(exclude_clusters).values)
    return adata[truth, :]


def install(a_new, b_new):
    M._shared_index_csr = NEW_SHARED if a_new else old_shared_index_csr
    f = NEW_FILTER if b_new else old_check_filter_single_cluster
    M.check_filter_single_cluster = f
    MC.check_filter_single_cluster = f


def run(a_new, b_new, outdir):
    install(a_new, b_new)
    M.RUN_ENRICHMENT = False
    adata = sc.read(DEMO)
    gc.collect()
    t0 = time.perf_counter()
    s = ScTriangulate(dir=outdir, adata=adata, query=KEYS, predict_doublet=False)
    s.run_enrichment = False
    s.lazy_run(compute_metrics_parallel=False, compute_shapley_parallel=False)
    secs = time.perf_counter() - t0
    labels = s.adata.obs['final_annotation'].astype(str).values.copy()
    del s, adata
    gc.collect()
    return secs, labels


def main():
    combos = list(itertools.product([False, True], [False, True]))
    order = combos + combos[::-1]          # each combination twice, both orders
    timings = {c: [] for c in combos}
    ref = None
    outdir = os.path.join(HERE, 'out_3chg_ablation')

    for i, (a, b) in enumerate(order):
        secs, labels = run(a, b, outdir)
        timings[(a, b)].append(secs)
        if ref is None:
            ref = labels
            agree = len(labels)
        else:
            agree = int((labels == ref).sum())
        print('run {}/{}  A={:<5} B={:<5}  {:6.1f} s   labels {} of {}'.format(
            i + 1, len(order), str(a), str(b), secs, agree, len(ref)))
        assert agree == len(ref), 'A={} B={} changed final_annotation'.format(a, b)

    out = []
    print('\n{:<22s} {:>10s} {:>10s} {:>10s}'.format('combination', 'run 1', 'run 2', 'median'))
    for (a, b), ts in timings.items():
        name = 'A={} B={}'.format('new' if a else 'old', 'new' if b else 'old')
        print('{:<22s} {:>9.1f}s {:>9.1f}s {:>9.1f}s'.format(name, ts[0], ts[1],
                                                             float(np.median(ts))))
        out.append({'A_new': a, 'B_new': b, 'seconds': ts,
                    'median_seconds': float(np.median(ts))})

    first = [r['seconds'][0] for r in out]
    second = [r['seconds'][1] for r in out]
    print('\nposition in the interpreter, not the change:')
    print('  mean of runs 1-4 (first pass):  {:.1f} s'.format(float(np.mean(first))))
    print('  mean of runs 5-8 (second pass): {:.1f} s'.format(float(np.mean(second))))

    path = os.path.join(HERE, 'three_changes_ablation.json')
    json.dump({'combinations': out,
               'mean_first_pass': float(np.mean(first)),
               'mean_second_pass': float(np.mean(second))}, open(path, 'w'), indent=1)
    print('\nwrote {}'.format(path))


if __name__ == '__main__':
    main()
