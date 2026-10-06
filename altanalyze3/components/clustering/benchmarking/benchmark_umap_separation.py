"""Paired, UMAP-only separation comparison with complete identity gates.

Cluster separation is descriptive, not biological validation. No cluster labels
are supplied to UMAP's objective, and no post-fit cluster movement is performed.
"""
import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

from .benchmark_landmark_umap import digest, load_verified_input


def input_neighbor_reference(matrix, states, k=15, per_state=5):
    """Exact correlation-space neighbors for a fixed stratified diagnostic set."""
    matrix = np.array(matrix, dtype=np.float32, copy=True)
    matrix -= matrix.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(matrix, axis=1)
    np.divide(matrix, norm[:, None], out=matrix, where=norm[:, None] > 0)
    states = np.asarray(states, dtype=str)
    rng = np.random.default_rng(29)
    anchors = np.sort(np.concatenate([rng.choice(np.flatnonzero(states == name),
        min(per_state, np.sum(states == name)), replace=False) for name in np.unique(states)]))
    k = min(k, len(matrix) - 1)
    neighbors = []
    for start in range(0, len(anchors), 16):
        selected = anchors[start:start + 16]
        distances = 1 - matrix[selected] @ matrix.T
        for i, row in enumerate(selected):
            if norm[row] == 0:
                distances[i, norm == 0] = 0
            distances[i, row] = np.inf
        neighbors.extend(np.argpartition(distances, k - 1, axis=1)[:, :k])
    return anchors, np.asarray(neighbors)


def input_neighbor_recall(coordinates, anchors, reference):
    from sklearn.neighbors import NearestNeighbors
    k = reference.shape[1]
    raw = NearestNeighbors(n_neighbors=k + 1, algorithm="kd_tree").fit(coordinates).kneighbors(
        coordinates[anchors], return_distance=False)
    neighbors = [[j for j in row if j != a][:k] for a, row in zip(anchors, raw)]
    return float(np.mean([len(set(a) & set(b)) / k for a, b in zip(neighbors, reference)]))


def separation(coordinates, states):
    """Scale-invariant nearest-centroid distances relative to RMS state radius."""
    from scipy.spatial.distance import cdist
    from sklearn.neighbors import NearestNeighbors
    coordinates = np.asarray(coordinates, dtype=np.float64)
    states = np.asarray(states, dtype=str)
    if coordinates.shape != (len(states), 2) or not np.isfinite(coordinates).all():
        raise ValueError("Finite coordinates for every state-labeled cell are required.")
    names = np.unique(states)
    if len(names) < 2:
        raise ValueError("Separation comparison requires at least two states.")
    centers, radii, anchors = [], [], []
    rng = np.random.default_rng(23)
    for name in names:
        rows = np.flatnonzero(states == name)
        center = coordinates[rows].mean(axis=0)
        centers.append(center)
        radii.append(np.sqrt(np.mean(np.sum((coordinates[rows] - center) ** 2, axis=1))))
        anchors.extend(rng.choice(rows, min(50, len(rows)), replace=False))
    centers, radii = np.asarray(centers), np.asarray(radii)
    distances = cdist(centers, centers)
    relative = np.divide(distances, radii[:, None] + radii[None, :],
                         out=np.full_like(distances, np.inf),
                         where=(radii[:, None] + radii[None, :]) > 0)
    np.fill_diagonal(distances, np.inf)
    np.fill_diagonal(relative, np.inf)
    nearest = relative.min(axis=1)
    anchors = np.asarray(anchors)
    k = min(15, len(states) - 1)
    raw = NearestNeighbors(n_neighbors=k + 1, algorithm="kd_tree").fit(coordinates).kneighbors(
        coordinates[anchors], return_distance=False)
    neighbors = np.array([[j for j in row if j != a][:k] for a, row in zip(anchors, raw)])
    purity = np.mean(states[neighbors] == states[anchors, None], axis=1)
    per_state = [{"state": name, "cells": int(np.sum(states == name)),
                  "rms_radius": float(radii[i]), "nearest_centroid_distance": float(distances[i].min()),
                  "nearest_relative_centroid_distance": float(nearest[i]),
                  "neighbor_state_purity": float(purity[states[anchors] == name].mean())}
                 for i, name in enumerate(names)]
    return {"median_nearest_relative_centroid_distance": float(np.median(nearest)),
            "p10_nearest_relative_centroid_distance": float(np.percentile(nearest, 10)),
            "macro_neighbor_state_purity": float(np.mean([s["neighbor_state_purity"] for s in per_state])),
            "per_state": per_state,
            "scope": "Existing-cluster geometry diagnostic; not biological validation. All states included."}


def fit(args):
    import scanpy as sc
    import anndata as ad
    import scipy.sparse as sp
    from altanalyze3.components.clustering.ICGS import _import_umap_with_local_retry
    from altanalyze3.components.clustering.umap_fit import fit_umap
    x, cells, genes, labels, manifest = load_verified_input(args.input)
    if not manifest.get("separation_benchmark_authorization"):
        raise ValueError("An actual user request to compare UMAP separation is required.")
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    if args.pcs:
        graph = ad.AnnData(sp.csr_matrix(x))
        sc.pp.pca(graph, n_comps=min(args.pcs, len(genes) - 1), zero_center=True,
                  svd_solver="arpack", random_state=0)
        representation = graph.obsm["X_pca"]
    else:
        representation = x
    model = _import_umap_with_local_retry().UMAP(
        n_neighbors=args.neighbors, min_dist=args.min_dist,
        metric=args.metric, random_state=0, repulsion_strength=args.repulsion_strength)
    coordinates, info, selected = fit_umap(
        model, lambda rows: np.asarray(representation[rows]), len(cells), mode=args.mode,
        labels=labels, max_fit_cells=30000, batch_cells=50000,
        min_per_state=args.min_per_state, seed=0, log=print)
    elapsed = time.perf_counter() - started
    np.save(args.output / "coordinates.npy", coordinates)
    np.save(args.output / "selected_rows.npy", selected)
    report = {"embedding_seconds": elapsed,
              "peak_rss_gib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss /
                  (1024**3 if sys.platform == "darwin" else 1024**2),
              "input_manifest_sha256": digest(args.input / "manifest.json"),
              "cells": len(cells), "features": len(genes), "parameters": dict(
                  info, n_pcs=args.pcs, metric=args.metric, n_neighbors=args.neighbors,
                  min_dist=args.min_dist, repulsion_strength=args.repulsion_strength),
              "separation": separation(coordinates, labels)}
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k:v for k,v in report.items() if k != "separation"}), flush=True)


def compile_comparison(root, output):
    """Evaluate every saved candidate against one exact input-neighbor reference."""
    from importlib.metadata import version
    from altanalyze3.components.clustering.umap_fit import select_landmarks
    root = Path(root)
    x, cells, genes, states, manifest = load_verified_input(root)
    anchors, reference = input_neighbor_reference(x, states)
    reports = {}
    for path in sorted(root.glob("*/report.json")):
        report = json.loads(path.read_text())
        if report["input_manifest_sha256"] != digest(root / "manifest.json"):
            raise ValueError("Candidate input provenance differs from the baseline.")
        coords = np.load(path.parent / "coordinates.npy")
        selected = np.load(path.parent / "selected_rows.npy")
        params = report["parameters"]
        expected = (np.arange(len(cells)) if params['fit_mode'] == 'full' else
                    select_landmarks(states, 30000, 0, params['min_landmarks_per_state']))
        if not np.array_equal(selected, expected):
            raise ValueError("Candidate fit identities differ from the recorded sampling rule.")
        report['separation'] = separation(coords, states)
        report['input_neighbor_recall'] = input_neighbor_recall(coords, anchors, reference)
        report['fit_identity_check'] = True
        reports[path.parent.name] = report
    saved = {}
    for name in ['current', 'previous']:
        coords = np.load(root / (name + '_coordinates.npy'))
        saved[name] = {'separation': separation(coords, states),
                       'input_neighbor_recall': input_neighbor_recall(coords, anchors, reference)}
    result = {'cells': len(cells), 'features': len(genes), 'states': len(np.unique(states)),
              'source_manifest': manifest, 'source_manifest_sha256': digest(root / 'manifest.json'),
              'versions': {name: version(name) for name in ['numpy', 'scipy', 'scanpy', 'anndata',
                                                           'umap-learn', 'scikit-learn']},
              'neighbor_evaluation': {'anchors': len(anchors), 'k': 15, 'max_per_state': 5,
                                      'seed': 29, 'metric': 'exact original-feature correlation'},
              'candidates': reports, 'saved_baselines': saved,
              'scope': 'UMAP-only comparison on identical final cells and marker panel; no biological accuracy claim.'}
    Path(output).write_text(json.dumps(result, indent=2) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--pcs", type=int, default=0)
    parser.add_argument("--neighbors", type=int, default=50)
    parser.add_argument("--min-dist", type=float, default=0.75)
    parser.add_argument("--metric", choices=["correlation", "euclidean"], default="correlation")
    parser.add_argument("--mode", choices=["full", "landmark"], default="landmark")
    parser.add_argument("--min-per-state", type=int, default=200)
    parser.add_argument("--repulsion-strength", type=float, default=1.0)
    fit(parser.parse_args())


if __name__ == "__main__":
    main()
