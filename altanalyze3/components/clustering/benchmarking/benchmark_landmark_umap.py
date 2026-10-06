"""Isolated UMAP fit/landmark-transform benchmark with explicit identity gates.

Requires a verified, prepared float32 UMAP input and its ordered cell/feature lists.
Does not reconstruct an absent feature panel or modify any production output.
Run each fit in a fresh process, then compare saved full and landmark outputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import resource
import sys
import time
from pathlib import Path

import numpy as np


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b''):
            value.update(block)
    return value.hexdigest()


def load_verified_input(root):
    manifest = json.loads((root / 'manifest.json').read_text())
    if manifest.get('baseline_verified') is not True:
        raise ValueError('Baseline source/algorithm verification is required before fitting.')
    required = {'X.npy', 'states.npy', 'cells.txt', 'features.txt'}
    if not required.issubset(manifest.get('sha256', {})):
        raise ValueError('Hashes for the matrix and complete ordered rosters are required.')
    for name, expected in manifest['sha256'].items():
        if digest(root / name) != expected:
            raise ValueError(f'Input provenance hash differs: {name}')
    x = np.load(root / 'X.npy', mmap_mode='r')
    cells = (root / 'cells.txt').read_text().splitlines()
    genes = (root / 'features.txt').read_text().splitlines()
    labels = np.load(root / 'states.npy', allow_pickle=False)
    if x.shape != (len(cells), len(genes)) or len(labels) != len(cells):
        raise ValueError('Complete ordered matrix, cell, feature and state rosters must agree.')
    if len(set(cells)) != len(cells) or len(set(genes)) != len(genes):
        raise ValueError('Duplicate cell or feature identities are unresolved.')
    if x.dtype != np.float32 or not x.flags.c_contiguous:
        raise ValueError('Use the exact C-contiguous float32 matrix consumed by baseline UMAP.')
    for start in range(0, len(cells), 4096):
        if not np.isfinite(x[start:start + 4096]).all():
            raise ValueError('Nonfinite UMAP inputs must be resolved before comparison.')
    return x, cells, genes, labels, manifest


def select_landmarks(labels, count, seed, strategy='state-stratified', min_per_state=200):
    n = len(labels)
    if not 2 < count < n:
        raise ValueError('Landmark count must be between 3 and n_cells - 1.')
    if strategy == 'state-stratified':
        from altanalyze3.components.clustering.umap_fit import select_landmarks as production_selection
        return production_selection(labels, count, seed, min_per_state)
    rng = np.random.default_rng(seed)
    chosen = []
    remaining = np.ones(n, dtype=bool)
    remaining[chosen] = False
    extra = rng.choice(np.flatnonzero(remaining), count - len(chosen), replace=False)
    return np.sort(np.concatenate([np.asarray(chosen, dtype=np.int64), extra]))


def peak_gib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024**3 if sys.platform == 'darwin' else 1024**2)


def run_fit(args):
    from importlib.metadata import version
    import umap
    root = args.input.resolve()
    start = time.perf_counter()
    x, cells, genes, states, manifest = load_verified_input(root)
    verification_seconds = time.perf_counter() - start
    if args.mode != 'full' and manifest.get('landmark_benchmark_authorized') is not True:
        raise ValueError('An explicit authorization record is required for the landmark comparison.')
    args.output.mkdir(parents=True, exist_ok=False)
    parameters = manifest['umap_parameters']
    model = umap.UMAP(**parameters, verbose=True)
    selected = np.arange(len(cells)) if args.mode == 'full' else select_landmarks(
        states, args.landmarks, int(parameters['random_state']), args.strategy, args.min_per_state)
    started = time.perf_counter()
    coordinates = np.empty((len(cells), 2), dtype=np.float32)
    if args.mode == 'landmark' and args.strategy == 'state-stratified':
        # Benchmark the same bounded fit/transform implementation the app uses.
        from altanalyze3.components.clustering.umap_fit import fit_umap
        coordinates, info, production_selected = fit_umap(
            model, lambda rows: np.asarray(x[rows]), len(cells), mode='landmark', labels=states,
            max_fit_cells=args.landmarks, batch_cells=args.batch_cells,
            min_per_state=args.min_per_state, seed=int(parameters['random_state']), log=print)
        if not np.array_equal(selected, production_selected):
            raise ValueError('Benchmark and production landmark identities differ.')
        fit_seconds, transform_seconds = info['fit_seconds'], info['transform_seconds']
    else:
        coordinates[selected] = model.fit_transform(x if args.mode == 'full' else np.asarray(x[selected]))
        fit_seconds = time.perf_counter() - started
        transform_seconds = 0.
    if args.mode != 'full' and args.strategy != 'state-stratified':
        missing = np.ones(len(cells), dtype=bool)
        missing[selected] = False
        rows = np.flatnonzero(missing)
        started = time.perf_counter()
        # Balanced chunks avoid a tiny tail changing UMAP's automatic transform
        # epochs. Batch size and model defaults are recorded, not assumed equal.
        for block in np.array_split(rows, max(1, math.ceil(len(rows) / args.batch_cells))):
            coordinates[block] = model.transform(np.asarray(x[block]))
            print(f'Transformed {len(block):,} cells; peak RSS {peak_gib():.2f} GiB', flush=True)
        transform_seconds = time.perf_counter() - started
    if not np.isfinite(coordinates).all():
        raise ValueError('UMAP returned missing coordinates; the comparison is incomplete.')
    np.save(args.output / 'coordinates.npy', coordinates)
    np.save(args.output / 'landmark_rows.npy', selected)
    if args.mode == 'full' and getattr(model, '_knn_indices', None) is not None:
        np.save(args.output / 'input_neighbors.npy', model._knn_indices)
    report = dict(mode=args.mode, strategy=args.strategy if args.mode != 'full' else 'all cells',
                  cells=len(cells), features=len(genes), landmarks=len(selected),
                  manifest_sha256=digest(root / 'manifest.json'), parameters=model.get_params(),
                  verification_seconds=verification_seconds, fit_seconds=fit_seconds,
                  transform_seconds=transform_seconds, embedding_seconds=fit_seconds + transform_seconds,
                  peak_rss_gib=peak_gib(), batch_cells=args.batch_cells,
                  min_per_state=args.min_per_state, finite_coordinates=True,
                  versions={name:version(name) for name in ('umap-learn','pynndescent','numpy','scipy','numba','scikit-learn')})
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report), flush=True)


def coordinate_neighbors(coords, anchors, k):
    from sklearn.neighbors import NearestNeighbors
    raw = NearestNeighbors(n_neighbors=k + 1, algorithm='kd_tree', n_jobs=1).fit(coords).kneighbors(
        coords[anchors], return_distance=False)
    return np.array([[j for j in row if j != anchor][:k] for anchor, row in zip(anchors, raw)])


def compare(args):
    x, cells, genes, states, manifest = load_verified_input(args.input)
    a = np.load(args.baseline / 'coordinates.npy', mmap_mode='r')
    b = np.load(args.candidate / 'coordinates.npy', mmap_mode='r')
    reports = [json.loads((path / 'report.json').read_text()) for path in (args.baseline, args.candidate)]
    for report, coords in zip(reports, (a,b)):
        if report['manifest_sha256'] != digest(args.input / 'manifest.json') or coords.shape != (len(cells),2):
            raise ValueError('Runs do not contain the same verified input and full cell roster.')
    rng = np.random.default_rng(23)
    anchors = np.sort(np.concatenate([rng.choice(np.flatnonzero(states == state),
                              min(50, int(np.sum(states == state))), replace=False)
                              for state in np.unique(states)]))
    k = min(15, len(cells)-1)
    an = coordinate_neighbors(a, anchors, k)
    bn = coordinate_neighbors(b, anchors, k)
    high_path = args.baseline / 'input_neighbors.npy'
    if high_path.exists():
        high_raw = np.load(high_path, mmap_mode='r')[anchors]
    else:
        from sklearn.neighbors import NearestNeighbors
        high_raw = NearestNeighbors(n_neighbors=k+1, metric=manifest['umap_parameters']['metric'],
                                    algorithm='brute', n_jobs=1).fit(x).kneighbors(x[anchors],return_distance=False)
    high = np.array([[j for j in row if j >= 0 and j != anchor][:k]
                     for anchor,row in zip(anchors, high_raw)])
    overlap = lambda p,q: np.array([len(set(left)&set(right))/k for left,right in zip(p,q)])
    ap = np.mean(states[an] == states[anchors,None], axis=1)
    bp = np.mean(states[bn] == states[anchors,None], axis=1)
    state_metrics = []
    ar,br,retention = overlap(an,high),overlap(bn,high),overlap(an,bn)
    for state in np.unique(states):
        m = states[anchors] == state
        state_metrics.append(dict(state=str(state), cells=int(np.sum(states == state)),
            evaluated_cells=int(np.sum(m)), baseline_input_neighbor_recall=float(ar[m].mean()),
            landmark_input_neighbor_recall=float(br[m].mean()),
            coordinate_neighbor_retention=float(retention[m].mean()),
            baseline_state_purity=float(ap[m].mean()), landmark_state_purity=float(bp[m].mean())))
    output = dict(cells=len(cells), features=len(genes), evaluated_cells=len(anchors), k=k,
                  note='Per-state stratified evaluation. State purity is an existing-cluster separation diagnostic, not an independent biological truth or proof of UMAP fidelity.',
                  baseline_input_neighbor_recall=float(ar.mean()), landmark_input_neighbor_recall=float(br.mean()),
                  coordinate_neighbor_retention=float(retention.mean()), states=state_metrics,
                  baseline_embedding_seconds=reports[0]['embedding_seconds'],
                  landmark_embedding_seconds=reports[1]['embedding_seconds'],
                  speedup=reports[0]['embedding_seconds']/reports[1]['embedding_seconds'],
                  baseline_peak_rss_gib=reports[0]['peak_rss_gib'],
                  landmark_peak_rss_gib=reports[1]['peak_rss_gib'])
    args.output.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k != 'states'}),flush=True)


def main():
    parser = argparse.ArgumentParser(__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    fit = sub.add_parser('fit')
    fit.add_argument('input', type=Path)
    fit.add_argument('output', type=Path)
    fit.add_argument('--mode', choices=['full','landmark'], required=True)
    fit.add_argument('--landmarks', type=int, default=30000)
    fit.add_argument('--strategy', choices=['random','state-stratified'], default='state-stratified')
    fit.add_argument('--min-per-state', type=int, default=200)
    fit.add_argument('--batch-cells', type=int, default=50000)
    comparison = sub.add_parser('compare')
    comparison.add_argument('input',type=Path)
    comparison.add_argument('baseline',type=Path)
    comparison.add_argument('candidate',type=Path)
    comparison.add_argument('output',type=Path)
    args = parser.parse_args()
    (run_fit if args.command=='fit' else compare)(args)


if __name__=='__main__': main()
