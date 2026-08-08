"""
Two questions about --downsample.

1. How many cells does the shared whitelist save against the naive scheme, where
   every annotation draws its own quota independently?
2. Does downsampling change the answer? Run the full pipeline with and without
   it and compare the labels on the cells both runs analysed.

Writes /Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/
sctriangulate/benchmarks/downsample_validate.json
"""

import os
import sys
import json
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))

from altanalyze3.components.sctriangulate import ScTriangulate            # noqa: E402
from altanalyze3.components.sctriangulate import metrics as M             # noqa: E402
from altanalyze3.components.sctriangulate.prefilter import (              # noqa: E402
    downsample_indices, _group_positions)

DEMO = os.path.join(HERE, 'data', 'demo_pbmc3k.h5ad')
SIM = os.path.join(HERE, 'data', 'sim_100k_20k.h5ad')


def independent_union(obs, keys, cap, seed=0):
    """Naive scheme: each annotation draws its own cap per cluster, blind to the
    others. Report the union, which is what you would actually keep."""
    rng = np.random.default_rng(seed)
    keep = np.zeros(obs.shape[0], dtype=bool)
    quota = 0
    for key in keys:
        _names, groups = _group_positions(obs[key].astype(str).to_numpy())
        for idx in groups:
            take = idx if len(idx) <= cap else rng.choice(idx, size=cap, replace=False)
            quota += len(take)
            keep[take] = True
    return int(keep.sum()), int(quota)


def coverage(obs_full, obs_kept, keys, cap):
    """Smallest ratio of kept to min(size, cap) over every cluster. 1.0 means no
    cluster of any annotation fell short of its cap."""
    worst = np.inf
    for key in keys:
        full = obs_full[key].value_counts()
        got = obs_kept[key].value_counts()
        for cluster, size in full.items():
            want = min(int(size), cap)
            have = int(got.get(cluster, 0))
            worst = min(worst, have / want)
    return float(worst)


def sizes_table(path, keys, caps, label):
    adata = sc.read(path)
    obs = adata.obs
    rows = []
    for cap in caps:
        t0 = time.perf_counter()
        idx, info = downsample_indices(obs, keys, cap, seed=0, report=False)
        secs = time.perf_counter() - t0
        naive_union, naive_quota = independent_union(obs, keys, cap, seed=0)
        rows.append({
            'dataset': label, 'n_cells': int(obs.shape[0]), 'cap': cap,
            'shared_whitelist': int(len(idx)),
            'independent_union': naive_union,
            'naive_quota_sum': naive_quota,
            'saved_vs_independent': naive_union - int(len(idx)),
            'fraction_of_full': len(idx) / obs.shape[0],
            'min_cluster_coverage': coverage(obs, obs.iloc[idx], keys, cap),
            'seconds': secs,
        })
        print('{:>10s} cap {:>5d}  shared {:>7d}  independent {:>7d}  '
              'naive quota {:>8d}  coverage {:.3f}  {:.2f}s'.format(
                  label, cap, len(idx), naive_union, naive_quota,
                  rows[-1]['min_cluster_coverage'], secs))
    del adata
    return rows


def run(cap, outdir):
    M.RUN_ENRICHMENT = False
    adata = sc.read(DEMO)
    keys = ['sctri_rna_leiden_1', 'sctri_rna_leiden_2', 'sctri_rna_leiden_3']
    t0 = time.perf_counter()
    s = ScTriangulate(dir=outdir, adata=adata, query=keys, predict_doublet=False)
    s.run_enrichment = False
    info = s.downsample_cells(max_per_cluster=cap) if cap else None
    s.lazy_run(compute_metrics_parallel=True, compute_shapley_parallel=True)
    return s.adata.obs.copy(), time.perf_counter() - t0, info


def main():
    keys = ['sctri_rna_leiden_1', 'sctri_rna_leiden_2', 'sctri_rna_leiden_3']
    out = {'size_scaling': [], 'label_effect': []}

    print('--- whitelist size: shared vs independent ---')
    out['size_scaling'] += sizes_table(DEMO, keys, [50, 100, 250, 500, 1000], 'demo_2700')
    if os.path.exists(SIM):
        sim = sc.read(SIM, backed='r')
        simkeys = [c for c in sim.obs.columns if c.startswith('sctri_')]
        del sim
        out['size_scaling'] += sizes_table(SIM, simkeys, [100, 250, 500, 1000], 'sim_100k')

    print('\n--- does it change the answer? demo, full pipeline ---')
    ref_obs, ref_t, _ = run(0, os.path.join(HERE, 'out_downsample_full'))
    print('full: {} cells, {:.1f}s'.format(ref_obs.shape[0], ref_t))

    for cap in [1000, 500, 250, 100]:
        obs, secs, info = run(cap, os.path.join(HERE, 'out_downsample_%d' % cap))
        shared = obs.index.intersection(ref_obs.index)
        row = {'cap': cap, 'cells': int(obs.shape[0]), 'shared_cells': int(len(shared)),
               'seconds': secs, 'full_seconds': ref_t,
               'n_pruned_full': int(ref_obs.loc[shared, 'pruned'].nunique()),
               'n_pruned_down': int(obs.loc[shared, 'pruned'].nunique())}
        for col in ['raw', 'pruned', 'final_annotation']:
            a = ref_obs.loc[shared, col].astype(str).values
            b = obs.loc[shared, col].astype(str).values
            row['ari_' + col] = float(adjusted_rand_score(a, b))
            row['nmi_' + col] = float(normalized_mutual_info_score(a, b))
            row['exact_' + col] = float((a == b).mean())
        out['label_effect'].append(row)
        print('cap {:>5d}  cells {:>5d}  ARI raw {:.4f}  pruned {:.4f}  final {:.4f}  '
              '{:.1f}s'.format(cap, obs.shape[0], row['ari_raw'], row['ari_pruned'],
                               row['ari_final_annotation'], secs))

    path = os.path.join(HERE, 'downsample_validate.json')
    json.dump(out, open(path, 'w'), indent=1)
    print('\nwrote {}'.format(path))


if __name__ == '__main__':
    main()
