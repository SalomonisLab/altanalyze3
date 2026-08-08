"""Compare the optimized SCCAF against the legacy one: score, speed, memory, ARI.

Two levels, because a change to one metric only matters if it survives to the
labels the user actually reads.

LEVEL 1, the metric itself. Runs both modes on each annotation and reports the
per-cluster accuracy difference, the wall clock and the peak memory of the
matrix handed to the solver. It also runs the legacy mode several times to
re-measure upstream's own spread, so the legacy-versus-optimized difference can
be read against the right yardstick.

LEVEL 2, the labels. Runs the whole pipeline through `lazy_run` once per mode
and scores the resulting labels against each other with exact agreement, ARI
and NMI. ARI ignores cluster names and compares the partitions.

Usage:
  python3.11 sccaf_compare.py [--h5ad data/demo_pbmc3k.h5ad] [--keys a,b,c]
  python3.11 sccaf_compare.py --child <mode> <h5ad> <outdir> <keys>
"""

import os
import sys
import json
import time
import shutil
import argparse
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)
PY = sys.executable
N_LEGACY_REPEATS = 3


def matrix_mb(X):
    from scipy.sparse import issparse as _sp
    if _sp(X):
        return (X.data.nbytes + X.indices.nbytes + X.indptr.nbytes) / 1e6
    return X.nbytes / 1e6


def level1(h5ad, keys):
    import scanpy as sc
    from altanalyze3.components.sctriangulate import ScTriangulate
    from altanalyze3.components.sctriangulate import metrics as M

    adata = sc.read(h5ad)
    sctri = ScTriangulate(dir=os.path.join(HERE, 'out_sccaf_cmp'), adata=adata,
                          query=keys, predict_doublet=False, verbose=0)
    # measure the shipped path: compute_metrics parks the sparse copy before
    # densifying, so SCCAF never pays a dense-to-sparse conversion
    sctri._install_sparse_sidecar()
    sctri._to_dense()

    print('LEVEL 1 -- the SCCAF metric itself')
    print()
    print('{:<26s} {:>5s} {:>9s} {:>9s} {:>8s} {:>9s} {:>9s} {:>11s} {:>11s}'.format(
        'annotation', 'cls', 'legacy s', 'optim s', 'speedup',
        'legacy MB', 'optim MB', 'mean|dAcc|', 'max|dAcc|'))
    print('-' * 108)

    rows = []
    for key in keys:
        adata_c = M.check_filter_single_cluster(sctri.adata, key)
        artifact = M._artifact_gene_set('human', 2)
        keep = np.asarray(~adata_c.var_names.isin(artifact))
        mb_dense = matrix_mb(M._as_dense(adata_c[:, keep].X))
        mb_sparse = matrix_mb(M._csr_from_columns(adata_c.X, keep))

        t0 = time.perf_counter()
        acc_leg, _ = M.SCCAF_score(adata_c, key, 'human', 2, False, mode='legacy')
        t_leg = time.perf_counter() - t0

        t0 = time.perf_counter()
        acc_opt, _ = M.SCCAF_score(adata_c, key, 'human', 2, False, mode='optimized')
        t_opt = time.perf_counter() - t0

        clusters = sorted(acc_leg)
        d = np.array([abs(float(acc_leg[c]) - float(acc_opt[c])) for c in clusters])

        # upstream's own spread on the same input, for scale
        repeats = [acc_leg]
        for _ in range(N_LEGACY_REPEATS - 1):
            a, _c = M.SCCAF_score(adata_c, key, 'human', 2, False, mode='legacy')
            repeats.append(a)
        legacy_spread = max(max(float(r[c]) for r in repeats) - min(float(r[c]) for r in repeats)
                            for c in clusters)
        opt_repeat, _ = M.SCCAF_score(adata_c, key, 'human', 2, False, mode='optimized')
        opt_stable = all(float(opt_repeat[c]) == float(acc_opt[c]) for c in clusters)

        print('{:<26s} {:>5d} {:>9.2f} {:>9.2f} {:>7.2f}x {:>9.0f} {:>9.0f} {:>11.4f} {:>11.4f}'.format(
            key, len(clusters), t_leg, t_opt, t_leg / t_opt, mb_dense, mb_sparse, d.mean(), d.max()))
        print('      legacy repeats identical: {} ; optimized repeats identical: {} ; '
              'legacy own spread over {} runs: {:.4f}'.format(
                  all(r == acc_leg for r in repeats), opt_stable, N_LEGACY_REPEATS, legacy_spread))

        rows.append({'key': key, 'n_clusters': len(clusters),
                     'legacy_seconds': t_leg, 'optimized_seconds': t_opt,
                     'legacy_matrix_mb': mb_dense, 'optimized_matrix_mb': mb_sparse,
                     'mean_abs_delta': float(d.mean()), 'max_abs_delta': float(d.max()),
                     'legacy_self_spread': float(legacy_spread),
                     'optimized_reproducible': bool(opt_stable),
                     'legacy_reproducible': bool(all(r == acc_leg for r in repeats))})
    print()
    return rows


CHILD = r'''
import os, sys, json, time


def main():
    mode, h5ad, outdir, keys = sys.argv[1:5]
    sys.path.insert(0, {repo!r})
    import matplotlib; matplotlib.use('Agg')
    import scanpy as sc
    from altanalyze3.components.sctriangulate import ScTriangulate
    adata = sc.read(h5ad)
    t0 = time.perf_counter()
    sctri = ScTriangulate(dir=outdir, adata=adata, query=keys.split(','), predict_doublet=False)
    sctri.sccaf_mode = mode
    sctri.lazy_run()
    wall = time.perf_counter() - t0
    obs = sctri.adata.obs
    json.dump({{'mode': mode, 'wall_seconds': wall,
               'pruned': list(obs['pruned'].astype(str).values),
               'raw': list(obs['raw'].astype(str).values),
               'final_annotation': list(obs['final_annotation'].astype(str).values)}},
              open(os.path.join(outdir, 'sccaf_mode_result.json'), 'w'))
    print('@@OK@@')


if __name__ == '__main__':
    main()
'''


def level2(h5ad, keys):
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    print('LEVEL 2 -- the labels the pipeline ends up with')
    print()
    results = {}
    for mode in ['legacy', 'optimized']:
        outdir = os.path.join(HERE, 'out_sccaf_mode_{}'.format(mode))
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.makedirs(outdir, exist_ok=True)
        script = os.path.join(outdir, '_child.py')
        with open(script, 'w') as f:
            f.write(CHILD.format(repo=REPO))
        env = dict(os.environ)
        env.setdefault('OMP_NUM_THREADS', '8')
        env['PYTHONHASHSEED'] = '1000'
        out = subprocess.run([PY, script, mode, os.path.abspath(h5ad), outdir, ','.join(keys)],
                             capture_output=True, text=True, env=env, cwd=HERE)
        if '@@OK@@' not in out.stdout:
            print(out.stdout[-3000:]); print(out.stderr[-3000:])
            raise RuntimeError('lazy_run failed for mode {}'.format(mode))
        results[mode] = json.load(open(os.path.join(outdir, 'sccaf_mode_result.json')))
        print('  lazy_run [{:<9s}] {:.2f} s'.format(mode, results[mode]['wall_seconds']))

    print()
    print('  {:<20s} {:>14s} {:>9s} {:>9s}'.format('COLUMN', 'EXACT', 'ARI', 'NMI'))
    print('  ' + '-' * 56)
    rows = []
    a, b = results['legacy'], results['optimized']
    for col in ['final_annotation', 'raw', 'pruned']:
        x, y = a[col], b[col]
        n = sum(1 for p, q in zip(x, y) if p == q)
        ari = adjusted_rand_score(x, y)
        nmi = normalized_mutual_info_score(x, y)
        print('  {:<20s} {:>7d}/{:<6d} {:>9.4f} {:>9.4f}'.format(col, n, len(x), ari, nmi))
        rows.append({'column': col, 'exact': n, 'n': len(x), 'ari': ari, 'nmi': nmi})
    print()
    print('  cluster counts: legacy {} | optimized {}'.format(
        {c: len(set(a[c])) for c in ['raw', 'pruned']},
        {c: len(set(b[c])) for c in ['raw', 'pruned']}))
    return {'wall': {m: results[m]['wall_seconds'] for m in results}, 'agreement': rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--h5ad', default=os.path.join(HERE, 'data', 'demo_pbmc3k.h5ad'))
    ap.add_argument('--keys', default='sctri_rna_leiden_1,sctri_rna_leiden_2,sctri_rna_leiden_3')
    ap.add_argument('--skip-level2', action='store_true')
    args = ap.parse_args()
    keys = args.keys.split(',')

    out = {'h5ad': os.path.abspath(args.h5ad), 'keys': keys}
    out['level1'] = level1(os.path.abspath(args.h5ad), keys)
    if not args.skip_level2:
        out['level2'] = level2(os.path.abspath(args.h5ad), keys)

    path = os.path.join(HERE, 'sccaf_compare_{}.json'.format(
        os.path.basename(args.h5ad).replace('.h5ad', '')))
    with open(path, 'w') as f:
        json.dump(out, f, indent=1)
    print('\nwrote {}'.format(path))


if __name__ == '__main__':
    main()
