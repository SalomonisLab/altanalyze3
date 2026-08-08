"""End-to-end wall clock through the real entry point, upstream vs integrated.

Runs `ScTriangulate(...).lazy_run()` -- the documented top-level function, not a
hand-rolled driver -- once for each implementation, in a fresh subprocess, and
reports wall clock plus the agreement of the final labels.

Settings timed for each implementation:
  sequential : compute_metrics_parallel=False   (pure algorithmic comparison)
  parallel   : compute_metrics_parallel=True    (the default a user gets)

Each sequential setting runs TWICE, under different PYTHONHASHSEED values. The
repeat measures how well an implementation reproduces ITSELF, which is the only
fair yardstick for how well the port reproduces upstream: two upstream runs
disagree with each other because upstream leaves two sklearn estimators
unseeded (see benchmarks/determinism_audit.py).

Usage:
  python3.11 end_to_end.py <input.h5ad> [key1,key2,...]
  python3.11 end_to_end.py --child <impl> <parallel> <h5ad> <outdir> <keys>
"""

import os
import sys
import json
import time
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
PY = sys.executable
TIMEOUT_SECONDS = 1500   # a config that exceeds this is reported as not finishing


def child(impl, parallel, h5ad, outdir, keys):
    parallel = (parallel == 'parallel')
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir, exist_ok=True)

    if impl == 'upstream':
        sys.path.insert(0, HERE)
        import _upstream_shim  # noqa: F401
        from sctriangulate import ScTriangulate
    else:
        sys.path.insert(0, REPO)
        from altanalyze3.components.sctriangulate import ScTriangulate

    import scanpy as sc
    import matplotlib
    matplotlib.use('Agg')

    adata = sc.read(h5ad)
    t0 = time.perf_counter()
    sctri = ScTriangulate(dir=outdir, adata=adata, query=keys, predict_doublet=False)
    sctri.lazy_run(compute_metrics_parallel=parallel,
                   compute_shapley_parallel=parallel)
    wall = time.perf_counter() - t0

    obs = sctri.adata.obs
    result = {
        'impl': impl,
        'parallel': parallel,
        'wall_seconds': wall,
        'n_cells': int(obs.shape[0]),
        'pruned': list(obs['pruned'].astype(str).values),
        'raw': list(obs['raw'].astype(str).values),
        'final_annotation': list(obs['final_annotation'].astype(str).values),
    }
    with open(os.path.join(outdir, 'end_to_end_result.json'), 'w') as f:
        json.dump(result, f)
    print('@@WALL@@{:.3f}'.format(wall))


def agreement(a, b):
    return sum(1 for x, y in zip(a, b) if x == y), len(a)


def driver(h5ad, keys):
    runs = {}
    configs = [('upstream', 'sequential', 1000), ('upstream', 'sequential', 2000),
               ('upstream', 'parallel', 1000),
               ('optimized', 'sequential', 1000), ('optimized', 'sequential', 2000),
               ('optimized', 'parallel', 1000)]
    for impl, mode, seed in configs:
        if True:
            tag = '{}_{}_{}'.format(impl, mode, seed)
            outdir = os.path.join(HERE, 'out_e2e_{}'.format(tag))
            cmd = [PY, os.path.abspath(__file__), '--child', impl, mode,
                   os.path.abspath(h5ad), outdir, ','.join(keys)]
            env = dict(os.environ)
            env.setdefault('OMP_NUM_THREADS', '8')
            env['PYTHONHASHSEED'] = str(seed)
            if impl == 'upstream':
                # Reach multiprocessing 'spawn' children. Without it, upstream's
                # workers re-import sctriangulate -> umap -> TensorFlow, die
                # during bootstrap, and the parent's pool.join() never returns.
                # Shimming them keeps the comparison about the parallel design.
                shim = os.path.join(HERE, '_shimpath')
                env['PYTHONPATH'] = shim + os.pathsep + env.get('PYTHONPATH', '')
            t0 = time.perf_counter()
            try:
                out = subprocess.run(cmd, capture_output=True, text=True, env=env,
                                     cwd=HERE, timeout=TIMEOUT_SECONDS)
            except subprocess.TimeoutExpired:
                print('{:<28s} DID NOT FINISH within {} s'.format(tag, TIMEOUT_SECONDS))
                runs[(impl, mode, seed)] = None
                continue
            elapsed = time.perf_counter() - t0
            if out.returncode != 0:
                print('--- {} FAILED ---'.format(tag))
                print(out.stdout[-4000:])
                print(out.stderr[-4000:])
                runs[(impl, mode, seed)] = None
                continue
            res = json.load(open(os.path.join(outdir, 'end_to_end_result.json')))
            res['process_seconds'] = elapsed
            res['tag'] = tag
            runs[(impl, mode, seed)] = res
            print('{:<28s} lazy_run {:8.2f} s   (process {:.2f} s)'.format(
                tag, res['wall_seconds'], elapsed))

    def get(impl, mode, seed=1000):
        return runs.get((impl, mode, seed))

    print()
    print('{:<52s} {:>10s}'.format('WALL CLOCK (lazy_run)', 'SPEEDUP'))
    print('-' * 66)
    for mode in ['sequential', 'parallel']:
        u, o = get('upstream', mode), get('optimized', mode)
        if u and o:
            print('{:<52s} {:>9.2f}x'.format(
                'upstream {0} -> optimized {0}'.format(mode),
                u['wall_seconds'] / o['wall_seconds']))
    u_seq, o_par = get('upstream', 'sequential'), get('optimized', 'parallel')
    if u_seq and o_par:
        print('{:<52s} {:>9.2f}x'.format('upstream sequential -> optimized parallel',
                                         u_seq['wall_seconds'] / o_par['wall_seconds']))

    print()
    print('REPRODUCIBILITY: does an implementation reproduce ITSELF?')
    print('(same input, two runs, different PYTHONHASHSEED)')
    print('-' * 66)
    for impl in ['upstream', 'optimized']:
        a, b = get(impl, 'sequential', 1000), get(impl, 'sequential', 2000)
        if not (a and b):
            continue
        for col in ['final_annotation', 'raw', 'pruned']:
            m, t = agreement(a[col], b[col])
            print('  {:<10s} run A vs run B  {:<18s} {:>6d}/{:<6d}  {:6.2f}%'.format(
                impl, col, m, t, 100.0 * m / t))
    a, b = get('optimized', 'sequential'), get('optimized', 'parallel')
    if a and b:
        for col in ['final_annotation', 'raw', 'pruned']:
            m, t = agreement(a[col], b[col])
            print('  {:<10s} seq vs parallel {:<18s} {:>6d}/{:<6d}  {:6.2f}%'.format(
                'optimized', col, m, t, 100.0 * m / t))

    print()
    print('PORT FIDELITY: optimized vs upstream, against that yardstick')
    print('-' * 66)
    u, o = get('upstream', 'sequential'), get('optimized', 'sequential')
    if u and o:
        for col in ['final_annotation', 'raw', 'pruned']:
            m, t = agreement(u[col], o[col])
            print('  optimized vs upstream    {:<18s} {:>6d}/{:<6d}  {:6.2f}%'.format(
                col, m, t, 100.0 * m / t))

    with open(os.path.join(HERE, 'end_to_end_summary.json'), 'w') as f:
        json.dump({'|'.join(map(str, k)): (None if v is None else
                                           {kk: vv for kk, vv in v.items()
                                            if kk not in ('pruned', 'raw', 'final_annotation')})
                   for k, v in runs.items()}, f, indent=1)
    print('\nwrote {}'.format(os.path.join(HERE, 'end_to_end_summary.json')))


if __name__ == '__main__':
    if sys.argv[1] == '--child':
        child(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6].split(','))
    else:
        ks = sys.argv[2].split(',') if len(sys.argv) > 2 else \
            ['sctri_rna_leiden_1', 'sctri_rna_leiden_2', 'sctri_rna_leiden_3']
        driver(sys.argv[1], ks)
