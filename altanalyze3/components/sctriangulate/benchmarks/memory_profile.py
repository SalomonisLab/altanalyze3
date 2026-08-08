"""Measure peak memory of upstream scTriangulate against this version.

WHY THIS MATTERS
----------------
Runtime is not the only reason a job fails on a laptop. Upstream's tf-idf step
builds a dense cells-by-genes pandas DataFrame and then takes full-matrix copies
of it once per cluster, so peak memory grows with cells times genes times a
constant that has nothing to do with the answer.

WHAT IS MEASURED
----------------
Peak resident set size of the whole process tree, sampled every 0.1 s with `ps`,
plus the kernel's own high-water mark from `/usr/bin/time -l` for the single
process case. The sampler catches concurrent workers, which the kernel figure
does not sum.

Sampling can miss a spike shorter than the interval, so the sampled figure is a
lower bound. The `time -l` figure is exact for the parent but reports only the
largest single child, not the sum, so neither number alone is complete. They are
reported side by side.

Usage:
  python3.11 memory_profile.py                       # demo, both modes
  python3.11 memory_profile.py --h5ad data/sim_8k_20k.h5ad \
      --keys sctri_sim_leiden_10,sctri_sim_leiden_20,sctri_sim_leiden_30
"""

import os
import re
import sys
import time
import json
import shutil
import argparse
import threading
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
PY = sys.executable
SAMPLE_SECONDS = 0.25   # each sample spawns a `ps`; 0.1 s added measurable load


def tree_rss_kb(root_pid):
    """Summed RSS in KB over root_pid and every descendant, via one ps call."""
    try:
        out = subprocess.run(['ps', '-eo', 'pid,ppid,rss'], capture_output=True,
                             text=True, timeout=5).stdout
    except Exception:
        return 0
    children, rss = {}, {}
    for line in out.splitlines()[1:]:
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            pid, ppid, r = int(parts[0]), int(parts[1]), int(parts[2])
        except ValueError:
            continue
        children.setdefault(ppid, []).append(pid)
        rss[pid] = r
    total, stack = 0, [root_pid]
    seen = set()
    while stack:
        p = stack.pop()
        if p in seen:
            continue
        seen.add(p)
        total += rss.get(p, 0)
        stack.extend(children.get(p, []))
    return total


class Sampler(threading.Thread):
    def __init__(self, pid):
        super().__init__(daemon=True)
        self.pid = pid
        self.peak_kb = 0
        self.stop_flag = False

    def run(self):
        while not self.stop_flag:
            self.peak_kb = max(self.peak_kb, tree_rss_kb(self.pid))
            time.sleep(SAMPLE_SECONDS)


# The __main__ guard is not optional. Under the 'spawn' start method every
# worker re-imports the main module, so without the guard each worker would
# re-run the whole analysis and the pool would never return. Leaving it out
# once already produced a fake "upstream parallel hangs" reading.
CHILD_SCRIPT = r'''
import os, sys, time


def main():
    impl, parallel, h5ad, outdir, keys = sys.argv[1:6]
    parallel = (parallel == 'parallel')
    keys = keys.split(',')
    if impl == 'upstream':
        sys.path.insert(0, {here!r})
        import _upstream_shim
        from sctriangulate import ScTriangulate
    else:
        sys.path.insert(0, {repo!r})
        from altanalyze3.components.sctriangulate import ScTriangulate
    import matplotlib; matplotlib.use('Agg')
    import scanpy as sc
    adata = sc.read(h5ad)
    t0 = time.perf_counter()
    sctri = ScTriangulate(dir=outdir, adata=adata, query=keys, predict_doublet=False)
    sctri.lazy_run(compute_metrics_parallel=parallel, compute_shapley_parallel=parallel)
    print('@@SECONDS@@{{:.2f}}'.format(time.perf_counter() - t0))


if __name__ == '__main__':
    main()
'''


def run_one(impl, mode, h5ad, keys, outdir, timeout=None):
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir, exist_ok=True)
    script = os.path.join(outdir, '_child.py')
    with open(script, 'w') as f:
        f.write(CHILD_SCRIPT.format(here=HERE, repo=REPO))

    env = dict(os.environ)
    env.setdefault('OMP_NUM_THREADS', '8')
    env['PYTHONHASHSEED'] = '1000'
    if impl == 'upstream':
        env['PYTHONPATH'] = os.path.join(HERE, '_shimpath') + os.pathsep + env.get('PYTHONPATH', '')

    cmd = ['/usr/bin/time', '-l', PY, script, impl, mode, h5ad, outdir, ','.join(keys)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, env=env, cwd=HERE)
    sampler = Sampler(proc.pid)
    sampler.start()
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        out, err = proc.communicate()
        sampler.stop_flag = True
        sampler.join(timeout=2)
        print('  {} {} did not finish within {} s; peak tree RSS at kill: {:.0f} MB'.format(
            impl, mode, timeout, sampler.peak_kb / 1024.0))
        return None
    sampler.stop_flag = True
    sampler.join(timeout=2)

    if proc.returncode != 0:
        print(out[-2500:]); print(err[-2500:])
        return None

    m = re.search(r'^\s*(\d+)\s+maximum resident set size', err, re.M)
    kernel_bytes = int(m.group(1)) if m else 0
    s = re.search(r'@@SECONDS@@([\d.]+)', out)
    return {'impl': impl, 'mode': mode,
            'kernel_peak_mb': kernel_bytes / 1e6,
            'sampled_tree_peak_mb': sampler.peak_kb / 1024.0,
            'seconds': float(s.group(1)) if s else None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--h5ad', default=os.path.join(HERE, 'data', 'demo_pbmc3k.h5ad'))
    ap.add_argument('--keys', default='sctri_rna_leiden_1,sctri_rna_leiden_2,sctri_rna_leiden_3')
    ap.add_argument('--modes', default='sequential,parallel')
    ap.add_argument('--impls', default='upstream,optimized')
    ap.add_argument('--timeout', type=float, default=None, help='seconds per run')
    args = ap.parse_args()

    keys = args.keys.split(',')
    import scanpy as sc
    shape = sc.read(args.h5ad).shape
    print('input: {}  ({} cells x {} genes, {} annotations)'.format(
        os.path.basename(args.h5ad), shape[0], shape[1], len(keys)))
    print()
    print('{:<12s} {:<12s} {:>14s} {:>18s} {:>10s}'.format(
        'impl', 'mode', 'kernel peak', 'sampled tree peak', 'seconds'))
    print('-' * 72)

    rows = []
    for mode in args.modes.split(','):
        for impl in args.impls.split(','):
            outdir = os.path.join(HERE, 'out_mem_{}_{}'.format(impl, mode))
            r = run_one(impl, mode, os.path.abspath(args.h5ad), keys, outdir, args.timeout)
            if r is None:
                print('{:<12s} {:<12s} {:>14s}'.format(impl, mode, 'FAILED'))
                continue
            rows.append(r)
            print('{:<12s} {:<12s} {:>11.0f} MB {:>15.0f} MB {:>9.1f}s'.format(
                impl, mode, r['kernel_peak_mb'], r['sampled_tree_peak_mb'], r['seconds']))
        u = [x for x in rows if x['mode'] == mode and x['impl'] == 'upstream']
        o = [x for x in rows if x['mode'] == mode and x['impl'] == 'optimized']
        if u and o:
            print('{:<12s} {:<12s} {:>11.2f}x {:>15.2f}x   <- reduction'.format(
                '', mode,
                u[0]['kernel_peak_mb'] / max(o[0]['kernel_peak_mb'], 1e-9),
                u[0]['sampled_tree_peak_mb'] / max(o[0]['sampled_tree_peak_mb'], 1e-9)))
        print()

    out = os.path.join(HERE, 'memory_profile_{}.json'.format(
        os.path.basename(args.h5ad).replace('.h5ad', '')))
    with open(out, 'w') as f:
        json.dump(rows, f, indent=1)
    print('wrote {}'.format(out))


if __name__ == '__main__':
    main()
