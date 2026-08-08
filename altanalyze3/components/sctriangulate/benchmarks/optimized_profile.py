"""Stage-level timing of the AltAnalyze3-integrated scTriangulate.

Mirrors benchmarks/baseline_profile.py stage for stage so the two logs and the
two JSON dumps line up. Sequential and single-process, same as the baseline.

Usage:
  python3.11 optimized_profile.py <input.h5ad> <outdir> [key1,key2,...]
"""

import os
import sys
import json
import time

import numpy as np
import scanpy as sc

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    '..', '..', '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from altanalyze3.components.sctriangulate import ScTriangulate                      # noqa: E402
from altanalyze3.components.sctriangulate.metrics import (                          # noqa: E402
    check_filter_single_cluster, marker_gene, reassign_score,
    tf_idf10_for_cluster, tf_idf5_for_cluster, SCCAF_score, doublet_compute,
)
from altanalyze3.components.sctriangulate.main_class import run_shapley, run_assign  # noqa: E402


class Timer:
    def __init__(self):
        self.t = {}

    def __call__(self, name):
        self._name = name
        return self

    def __enter__(self):
        self._t0 = time.perf_counter()
        return self

    def __exit__(self, *a):
        dt = time.perf_counter() - self._t0
        self.t[self._name] = self.t.get(self._name, 0.0) + dt
        print('    {:<28s} {:8.2f} s'.format(self._name, dt), flush=True)
        return False


def main():
    h5ad = sys.argv[1]
    outdir = sys.argv[2]
    keys = sys.argv[3].split(',') if len(sys.argv) > 3 else \
        ['sctri_rna_leiden_1', 'sctri_rna_leiden_2', 'sctri_rna_leiden_3']
    os.makedirs(outdir, exist_ok=True)

    adata = sc.read(h5ad)
    print('adata: {} cells x {} genes'.format(adata.shape[0], adata.shape[1]), flush=True)

    timer = Timer()
    with timer('00_construct_object'):
        sctri = ScTriangulate(dir=outdir, adata=adata, query=keys, predict_doublet=False)
    sctri._to_dense()

    reference = {}
    for key in keys:
        print('  --- key {} ---'.format(key), flush=True)
        adata_c = check_filter_single_cluster(sctri.adata, key)
        with timer('01_marker_gene[{}]'.format(key)):
            mg = marker_gene(adata_c, key, sctri.species, sctri.criterion, sctri.dir)
        with timer('02_reassign[{}]'.format(key)):
            c_reassign, conf_reassign = reassign_score(adata_c, key, mg)
        with timer('03_tfidf10[{}]'.format(key)):
            c_tfidf10, exclusive = tf_idf10_for_cluster(adata_c, key, sctri.species, sctri.criterion)
        with timer('04_tfidf5[{}]'.format(key)):
            c_tfidf5 = tf_idf5_for_cluster(adata_c, key, sctri.species, sctri.criterion)
        with timer('05_sccaf[{}]'.format(key)):
            c_sccaf, conf_sccaf = SCCAF_score(adata_c, key, sctri.species, sctri.criterion, False)
        with timer('06_doublet[{}]'.format(key)):
            c_doublet = doublet_compute(adata_c, key)

        reference[key] = {
            'marker_genes': {c: list(mg.loc[c, 'whole_marker_genes']) for c in mg.index},
            'purify': {c: list(mg.loc[c, 'purify']) for c in mg.index},
            'enrichr': {c: mg.loc[c, 'enrichr'] for c in mg.index},
            'gsea': {c: {k2: list(v2) for k2, v2 in mg.loc[c, 'gsea'].items()} for c in mg.index},
            'reassign': {str(k2): float(v2) for k2, v2 in c_reassign.items()},
            'tfidf10': {str(k2): float(v2) for k2, v2 in c_tfidf10.items()},
            'tfidf5': {str(k2): float(v2) for k2, v2 in c_tfidf5.items()},
            'SCCAF': {str(k2): float(v2) for k2, v2 in c_sccaf.items()},
            'doublet': {str(k2): float(v2) for k2, v2 in c_doublet.items()},
            'exclusive_genes_head': {c: list(exclusive[c].keys())[:50] for c in exclusive.index},
            'confusion_reassign': conf_reassign.values.tolist(),
            'confusion_sccaf': conf_sccaf.values.tolist(),
        }
        for metric, d in [('reassign', c_reassign), ('tfidf10', c_tfidf10),
                          ('SCCAF', c_sccaf), ('doublet', c_doublet), ('tfidf5', c_tfidf5)]:
            sctri.adata.obs['{}@{}'.format(metric, key)] = \
                sctri.adata.obs[key].astype('str').map(d).fillna(0).values

    sctri.total_metrics = ['reassign', 'tfidf10', 'SCCAF', 'doublet', 'tfidf5']
    score_colname = [m for m in sctri.total_metrics if m != 'doublet']
    obs = sctri.adata.obs
    data = np.empty([len(keys), obs.shape[0], len(score_colname)])
    for i, key in enumerate(keys):
        data[i, :, :] = obs[[n + '@' + key for n in score_colname]].values

    with timer('07_shapley_all_or_none'):
        final, intermediate = run_shapley(obs, keys, sctri.reference, sctri.size_dict,
                                          data, 'shapley_all_or_none', 0.01)
    reference['shapley_all_or_none'] = {'final': list(final),
                                        'values': np.asarray(intermediate).tolist()}

    with timer('08_shapley'):
        f2, i2 = run_shapley(obs, keys, sctri.reference, sctri.size_dict, data, 'shapley', 0.01)
    reference['shapley'] = {'final': list(f2), 'values': np.asarray(i2).tolist()}

    with timer('09_rank'):
        f3, i3 = run_shapley(obs, keys, sctri.reference, sctri.size_dict, data, 'rank', 0.01)
    reference['rank'] = {'final': list(f3), 'values': np.asarray(i3).tolist()}

    with timer('10_rank_all_or_none'):
        f4, i4 = run_shapley(obs, keys, sctri.reference, sctri.size_dict, data,
                             'rank_all_or_none', 0.01)
    reference['rank_all_or_none'] = {'final': list(f4), 'values': np.asarray(i4).tolist()}

    obs['final_annotation'] = final
    with timer('11_run_assign'):
        obs2 = run_assign(obs.copy())
    reference['raw'] = list(obs2['raw'].values)

    reference['_timings'] = timer.t
    reference['_shape'] = list(adata.shape)
    reference['_keys'] = keys
    reference['_data_matrix'] = data.tolist()

    out = os.path.join(outdir, 'optimized_reference.json')
    with open(out, 'w') as f:
        json.dump(reference, f)
    print('\nwrote {}'.format(out), flush=True)
    print('\nTOTAL {:.2f} s'.format(sum(timer.t.values())), flush=True)
    for k in sorted(timer.t):
        print('  {:<32s} {:8.2f} s'.format(k, timer.t[k]))


if __name__ == '__main__':
    main()
