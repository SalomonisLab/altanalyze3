"""Shapley engine: upstream vs batched, as the number of annotations grows.

The bundled demo has 3 annotations, where the coalition enumeration is only
2**2 = 4 terms per player and the upstream loop looks cheap. scTriangulate is
designed for more than that (upstream switches to 'rank' mode above 15
annotations), and the enumeration is exponential, so this script measures the
scaling on simulated score matrices.

The simulation is deliberately realistic in the way that matters for the engine:
scores are per-cluster constants broadcast to cells, so cells sharing a cluster
tuple share a score matrix -- exactly the structure the deduplication exploits.
Values are drawn on a coarse grid so exact ties occur, which is the branch that
`cheat_add_bonus` decides.

Correctness is checked against the verbatim upstream functions on every
configuration, not just timed.

Usage:
  python3.11 shapley_scaling.py [--cells 20000]
"""

import os
import sys
import json
import time
import argparse
import importlib.util

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
COMPONENT = os.path.dirname(HERE)
REPO = os.path.abspath(os.path.join(COMPONENT, '..', '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from altanalyze3.components.sctriangulate import shapley as fast   # noqa: E402


def load_reference():
    path = os.path.join(COMPONENT, '_reference', 'upstream_shapley.py')
    spec = importlib.util.spec_from_file_location('upstream_shapley', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


REF = load_reference()
BONUS = 0.01


def simulate(n_annotations, n_cells, n_metrics, n_clusters, seed=0):
    """Per-cluster metric values broadcast to cells, as the real pipeline builds them."""
    rng = np.random.default_rng(seed)
    data = np.empty((n_annotations, n_cells, n_metrics))
    obs = {}
    for a in range(n_annotations):
        # coarse grid -> genuine ties between annotations
        table = np.round(rng.random((n_clusters, n_metrics)), 2)
        assign = rng.integers(0, n_clusters, n_cells)
        data[a] = table[assign]
        obs['ann{}'.format(a)] = assign.astype(str)
    return data, pd.DataFrame(obs)


def upstream_scores(data, mode, n_cells_limit):
    n_players, _, _ = data.shape
    out = np.empty((n_cells_limit, n_players))
    for i in range(n_cells_limit):
        layer = data[:, i, :]
        for j in range(n_players):
            out[i, j] = REF.wrapper_shapley(j, layer, mode, BONUS)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cells', type=int, default=20000)
    ap.add_argument('--metrics', type=int, default=4)
    ap.add_argument('--clusters', type=int, default=12)
    ap.add_argument('--verify-cells', type=int, default=300,
                    help='how many cells to also score with the upstream loop for verification')
    args = ap.parse_args()

    mode = 'shapley_all_or_none'
    rows = []
    print('mode={}  cells={}  metrics={}  clusters/annotation={}'.format(
        mode, args.cells, args.metrics, args.clusters))
    print()
    print('{:>5s} {:>12s} {:>14s} {:>14s} {:>10s} {:>9s} {:>10s}'.format(
        'ann', 'coalitions', 'upstream (s)', 'batched (s)', 'speedup', 'unique', 'verified'))
    print('-' * 82)

    for n_ann in [3, 4, 5, 6, 8, 10, 12]:
        data, obs = simulate(n_ann, args.cells, args.metrics, args.clusters)

        # upstream: time a subset, extrapolate linearly (it is strictly per-cell)
        n_probe = min(args.verify_cells, args.cells)
        t0 = time.perf_counter()
        ref_scores = upstream_scores(data, mode, n_probe)
        t_probe = time.perf_counter() - t0
        t_upstream_full = t_probe * args.cells / n_probe

        t0 = time.perf_counter()
        got = fast.shapley_batch(data, mode=mode, bonus=BONUS)
        t_fast = time.perf_counter() - t0

        ok = np.allclose(got[:n_probe], ref_scores, rtol=0, atol=1e-9)
        n_unique = len(np.unique(
            np.transpose(data, (1, 0, 2)).reshape(args.cells, -1), axis=0))

        n_coalitions = n_ann * (2 ** (n_ann - 1) - 1)
        print('{:>5d} {:>12d} {:>14.2f} {:>14.3f} {:>9.1f}x {:>9d} {:>10s}'.format(
            n_ann, n_coalitions, t_upstream_full, t_fast, t_upstream_full / t_fast,
            n_unique, 'yes' if ok else 'NO'))
        rows.append({'n_annotations': n_ann, 'coalitions_per_cell': n_coalitions,
                     'upstream_seconds_extrapolated': t_upstream_full,
                     'batched_seconds': t_fast,
                     'speedup': t_upstream_full / t_fast,
                     'unique_score_matrices': n_unique,
                     'n_cells': args.cells, 'verified_on_cells': n_probe,
                     'verified': bool(ok)})
        if not ok:
            raise SystemExit('MISMATCH against upstream at {} annotations'.format(n_ann))

    print()
    print('upstream times for {} cells are extrapolated from the first {} cells '
          '(the upstream loop is strictly per-cell, so this is linear).'.format(
              args.cells, args.verify_cells))

    out = os.path.join(HERE, 'shapley_scaling.json')
    with open(out, 'w') as f:
        json.dump(rows, f, indent=1)
    print('wrote {}'.format(out))


if __name__ == '__main__':
    main()
