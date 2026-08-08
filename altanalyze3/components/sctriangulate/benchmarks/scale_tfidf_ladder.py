"""Peak memory of the tf-idf step against cell count, upstream vs this version.

tf-idf is the step that decides whether a large run survives. Upstream builds a
dense cells-by-genes DataFrame and then takes a full-matrix `.values` copy for
every cluster in turn. This version counts nonzeros in fixed-size chunks.

Each measurement runs in its own process so `ru_maxrss` is that step alone.
A ladder of cell counts is used instead of one run at 100,000, because the
upstream allocation at that size is large enough to disturb the machine.

Usage: python3.11 scale_tfidf_ladder.py [--genes 20000]
"""
import os, sys, json, time, resource, argparse, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
PY = sys.executable


def child(impl, n_cells, n_genes, n_clusters):
    sys.path.insert(0, HERE); sys.path.insert(0, REPO)
    import numpy as np, pandas as pd, anndata as ad
    from scipy.sparse import csr_matrix
    rng = np.random.default_rng(0)
    # float32 directly and in row blocks: rng.random(shape) would build a float64
    # array first, which at 40,000 cells is 6.4 GB and would dominate ru_maxrss
    X = np.zeros((n_cells, n_genes), dtype=np.float32)
    for st in range(0, n_cells, 2000):
        sp = min(st + 2000, n_cells)
        blk = rng.random((sp - st, n_genes), dtype=np.float32)
        blk[blk < 0.94] = 0.0                             # ~6% nonzero
        X[st:sp] = blk
        del blk
    labels = rng.integers(0, n_clusters, n_cells).astype(str)
    adata = ad.AnnData(X=X,
                       obs=pd.DataFrame({'k': pd.Categorical(labels)},
                                        index=['c%d' % i for i in range(n_cells)]),
                       var=pd.DataFrame(index=['G%05d' % i for i in range(n_genes)]))
    import gc; gc.collect()
    base = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    t0 = time.perf_counter()
    if impl == 'upstream':
        import _upstream_shim  # noqa: F401
        from sctriangulate.metrics import tf_idf_bare_compute
        df = pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names)
        df['cluster'] = adata.obs['k'].astype('str').values
        for c in adata.obs['k'].cat.categories:
            tf_idf_bare_compute(df, c)
    else:
        from altanalyze3.components.sctriangulate import metrics as M
        M._TFIDF_CACHE.update({'adata': None, 'matrix': None, 'ranked': None})
        M._tfidf_matrix(adata, 'k')
    dt = time.perf_counter() - t0
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9   # macOS: bytes
    print('@@R@@' + json.dumps({'impl': impl, 'cells': n_cells, 'seconds': dt,
                                'peak_gb': peak, 'base_gb': base,
                                'step_gb': max(peak - base, 0.0)}))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--genes', type=int, default=20000)
    ap.add_argument('--clusters', type=int, default=40)
    ap.add_argument('--ladder', default='5000,10000,20000,40000')
    a = ap.parse_args()
    print('tf-idf step only, {} genes, {} clusters, ~6% nonzero\n'.format(a.genes, a.clusters))
    print('{:>8s} {:>12s} {:>12s} {:>10s} {:>12s} {:>12s} {:>9s}'.format(
        'cells', 'upstream +GB', 'this +GB', 'GB ratio', 'upstream s', 'this s', 'speedup'))
    print('-' * 82)
    rows = []
    for n in [int(x) for x in a.ladder.split(',')]:
        got = {}
        for impl in ['upstream', 'optimized']:
            env = dict(os.environ, PYTHONPATH=os.path.join(HERE, '_shimpath'))
            out = subprocess.run([PY, os.path.abspath(__file__), '--child', impl,
                                  str(n), str(a.genes), str(a.clusters)],
                                 capture_output=True, text=True, env=env, cwd=HERE)
            line = [l for l in out.stdout.splitlines() if l.startswith('@@R@@')]
            if not line:
                print(out.stderr[-1500:]); got[impl] = None; continue
            got[impl] = json.loads(line[0][5:])
        u, o = got.get('upstream'), got.get('optimized')
        if u and o:
            print('{:>8d} {:>12.2f} {:>12.2f} {:>9.1f}x {:>12.1f} {:>12.1f} {:>8.1f}x'.format(
                n, u['step_gb'], o['step_gb'], u['step_gb'] / max(o['step_gb'],1e-9),
                u['seconds'], o['seconds'], u['seconds'] / o['seconds']))
            rows.append({'cells': n, 'upstream': u, 'optimized': o})
    json.dump(rows, open(os.path.join(HERE, 'scale_tfidf_ladder.json'), 'w'), indent=1)
    if len(rows) >= 2:
        f, l = rows[0], rows[-1]
        du = (l['upstream']['step_gb'] - f['upstream']['step_gb']) / (l['cells'] - f['cells'])
        do = (l['optimized']['step_gb'] - f['optimized']['step_gb']) / (l['cells'] - f['cells'])
        print('\nmeasured slope, GB per 1,000 cells: upstream {:.3f}, this version {:.3f}'.format(
            du * 1000, do * 1000))
        print('additional GB projected to 100,000 cells: upstream {:.1f}, this version {:.1f}'.format(
            f['upstream']['step_gb'] + du * (100000 - f['cells']),
            f['optimized']['step_gb'] + do * (100000 - f['cells'])))


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--child':
        child(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))
    else:
        main()
