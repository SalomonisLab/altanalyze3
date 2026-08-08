"""Quantify upstream scTriangulate's run-to-run variation in its stability metrics.

WHY THIS EXISTS
---------------
`reassign_score` and `SCCAF_score` are the two stability metrics whose optimized
values do not reproduce the upstream values bit for bit. This script shows why:
neither upstream function returns the same number twice on identical input.

`reassign_score` has two unseeded sources, measured separately:

  (a) `PCA(n_components=30)` with `random_state=None`. sklearn's `svd_solver='auto'`
      picks the randomized solver at this shape, so repeated calls differ.
      Measured by repeating the call inside ONE process.

  (b) `pool = list(set(pool))`. Python randomises string hashing per interpreter,
      so the gene column order fed to PCA differs between processes -- including
      between the worker processes upstream spawns in `compute_metrics`.
      Measured by re-running under different PYTHONHASHSEED values.

`SCCAF_score` has one: `LogisticRegression(solver='liblinear')` with
`random_state=None`. sklearn uses that argument to shuffle the data for the
liblinear solver.

The optimized implementations sort the pool and seed both estimators, so each
returns one value. The question this script answers is whether those values sit
inside the range upstream itself produces.

Marker genes are computed once with the optimized `marker_gene`, whose output is
proven bit-identical to upstream's by `compare_outputs.py` (all
marker/purify/exclusive gene lists match). Upstream's own `reassign_score` is
then called on exactly that input.

Usage:
  # child mode (one process, N repeats, prints JSON)
  python3.11 determinism_audit.py --child <h5ad> <key> <n_repeats>
  # driver
  python3.11 determinism_audit.py <h5ad> [key1,key2,...]
"""

import os
import sys
import json
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
PY = sys.executable
N_REPEATS_IN_PROCESS = 10
N_PROCESSES = 5


def child(h5ad, key, n_repeats):
    sys.path.insert(0, HERE)
    import _upstream_shim  # noqa: F401
    sys.path.insert(0, REPO)

    import scanpy as sc
    from altanalyze3.components.sctriangulate import ScTriangulate
    from altanalyze3.components.sctriangulate.metrics import (
        check_filter_single_cluster, marker_gene,
        reassign_score as reassign_optimized,
        SCCAF_score as sccaf_optimized,
    )
    from sctriangulate.metrics import (
        reassign_score as reassign_upstream,
        SCCAF_score as sccaf_upstream,
    )

    adata = sc.read(h5ad)
    sctri = ScTriangulate(dir=os.path.join(HERE, 'out_determinism'), adata=adata,
                          query=[key], predict_doublet=False, verbose=0)
    sctri._to_dense()
    adata_c = check_filter_single_cluster(sctri.adata, key)
    mg = marker_gene(adata_c, key, 'human', 2, sctri.dir, run_enrichment=False)

    out = {'hashseed': os.environ.get('PYTHONHASHSEED')}

    runs = []
    for _ in range(n_repeats):
        d, _conf = reassign_upstream(adata_c, key, mg)
        runs.append({str(k): float(v) for k, v in d.items()})
    out['reassign_upstream'] = runs
    runs = []
    for _ in range(3):
        d, _conf = reassign_optimized(adata_c, key, mg)
        runs.append({str(k): float(v) for k, v in d.items()})
    out['reassign_optimized'] = runs

    runs = []
    for _ in range(n_repeats):
        d, _conf = sccaf_upstream(adata_c, key, 'human', 2, False)
        runs.append({str(k): float(v) for k, v in d.items()})
    out['SCCAF_upstream'] = runs
    runs = []
    for _ in range(3):
        d, _conf = sccaf_optimized(adata_c, key, 'human', 2, False)
        runs.append({str(k): float(v) for k, v in d.items()})
    out['SCCAF_optimized'] = runs

    print('@@JSON@@' + json.dumps(out))


def _summarise(metric, key, per_process, n_processes, n_repeats):
    up_key, opt_key = metric + '_upstream', metric + '_optimized'
    clusters = sorted(per_process[0][up_key][0])

    within = {c: [] for c in clusters}   # one process: solver seed only
    across = {c: [] for c in clusters}   # all processes: solver seed + hash order
    for p in per_process:
        for r in p[up_key]:
            for c in clusters:
                across[c].append(r[c])
    for r in per_process[0][up_key]:
        for c in clusters:
            within[c].append(r[c])

    first = per_process[0][opt_key][0]
    opt_stable = all(r == first for p in per_process for r in p[opt_key])
    n_ident = sum(1 for p in per_process for r in p[up_key] if r == per_process[0][up_key][0])
    n_total = n_processes * n_repeats

    within_spread = {c: max(within[c]) - min(within[c]) for c in clusters}
    across_spread = {c: max(across[c]) - min(across[c]) for c in clusters}
    inside = {c: (min(across[c]) - 1e-12) <= first[c] <= (max(across[c]) + 1e-12)
              for c in clusters}
    n_inside = sum(inside.values())

    print('  [{}]'.format(metric))
    print('    upstream runs identical to the first : {}/{}'.format(n_ident, n_total))
    print('    optimized identical across {}x3 runs  : {}'.format(n_processes, opt_stable))
    print('    upstream spread, one process         : max {:.6f}  mean {:.6f}'.format(
        max(within_spread.values()), float(np.mean(list(within_spread.values())))))
    print('    upstream spread, across processes    : max {:.6f}  mean {:.6f}'.format(
        max(across_spread.values()), float(np.mean(list(across_spread.values())))))
    print('    optimized inside upstream range      : {}/{} clusters'.format(
        n_inside, len(clusters)))
    for c in clusters:
        if not inside[c]:
            print('      OUTSIDE: cluster {} optimized={:.6f} upstream range=[{:.6f},{:.6f}]'.format(
                c, first[c], min(across[c]), max(across[c])))

    return {'metric': metric, 'key': key, 'clusters': clusters,
            'optimized': first,
            'optimized_stable': bool(opt_stable),
            'upstream_identical_runs': n_ident, 'upstream_total_runs': n_total,
            'upstream_min': {c: min(across[c]) for c in clusters},
            'upstream_max': {c: max(across[c]) for c in clusters},
            'within_process_spread': within_spread,
            'across_process_spread': across_spread,
            'n_inside_upstream_range': n_inside,
            'n_clusters': len(clusters),
            'all_inside': bool(n_inside == len(clusters))}


def driver(h5ad, keys):
    all_rows = []
    for key in keys:
        per_process = []
        for i in range(N_PROCESSES):
            env = dict(os.environ, PYTHONHASHSEED=str(1000 + i))
            cmd = [PY, os.path.abspath(__file__), '--child', h5ad, key,
                   str(N_REPEATS_IN_PROCESS)]
            out = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=HERE)
            line = [l for l in out.stdout.splitlines() if l.startswith('@@JSON@@')]
            if not line:
                print(out.stdout[-3000:]); print(out.stderr[-3000:])
                raise RuntimeError('child failed for {}'.format(key))
            per_process.append(json.loads(line[0][len('@@JSON@@'):]))

        print('\n=== {} ({} clusters) ==='.format(
            key, len(per_process[0]['reassign_upstream'][0])))
        for metric in ['reassign', 'SCCAF']:
            all_rows.append(_summarise(metric, key, per_process,
                                       N_PROCESSES, N_REPEATS_IN_PROCESS))

    out_path = os.path.join(HERE, 'out_determinism', 'determinism_audit.json')
    with open(out_path, 'w') as f:
        json.dump(all_rows, f, indent=1)

    ok = all(r['all_inside'] for r in all_rows)
    stable = all(r['optimized_stable'] for r in all_rows)
    print('\nOPTIMIZED DETERMINISTIC EVERYWHERE      : {}'.format(stable))
    print('OPTIMIZED INSIDE UPSTREAM RANGE ALWAYS  : {}'.format(ok))
    print('wrote {}'.format(out_path))
    return 0 if (ok and stable) else 1


if __name__ == '__main__':
    os.makedirs(os.path.join(HERE, 'out_determinism'), exist_ok=True)
    if sys.argv[1] == '--child':
        child(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    else:
        keys = sys.argv[2].split(',') if len(sys.argv) > 2 else \
            ['sctri_rna_leiden_1', 'sctri_rna_leiden_2', 'sctri_rna_leiden_3']
        sys.exit(driver(os.path.abspath(sys.argv[1]), keys))
