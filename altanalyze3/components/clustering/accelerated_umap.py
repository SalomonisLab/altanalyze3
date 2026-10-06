"""UMAP-only PCA projection; independent of clustering, markers and expression."""
from __future__ import annotations

import time

import anndata as ad
import numpy as np
import scanpy as sc


def pca_umap(matrix, *, method="scanpy", labels=None, n_pcs=50, n_neighbors=15,
             min_dist=0.75, seed=0, fit_cells=30000, batch_cells=50000, log=None):
    if method not in {"scanpy", "landmark"}:
        raise ValueError(f"Unknown accelerated UMAP method: {method}")
    cells, features = matrix.shape
    if cells < 3 or features < 2:
        raise ValueError("PCA UMAP requires at least three cells and two input features.")
    components = min(n_pcs, cells - 1, features - 1)
    neighbors = min(n_neighbors, cells - 1)
    if log:
        log(f"UMAP PCA: {cells:,} cells x {features:,} features to {components} components")
    graph = ad.AnnData(X=matrix)
    started = time.perf_counter()
    # Standard centered Scanpy PCA. No scaling, HVG filtering or expression mutation.
    sc.pp.pca(graph, n_comps=components, zero_center=True, svd_solver="arpack",
              random_state=seed)
    pca_seconds = time.perf_counter() - started
    if log:
        log(f"UMAP PCA completed in {pca_seconds:.1f} s")
    started = time.perf_counter()
    if method == "scanpy":
        if log:
            log(f"UMAP fitting PCA neighbors for all {cells:,} cells")
        sc.pp.neighbors(graph, n_neighbors=neighbors, use_rep="X_pca",
                        metric="euclidean", random_state=seed)
        neighbors_seconds = time.perf_counter() - started
        started = time.perf_counter()
        sc.tl.umap(graph, min_dist=min_dist, random_state=seed)
        coordinates = graph.obsm["X_umap"]
        info = {"fit_mode": "full", "fit_cells": cells, "total_cells": cells,
                "neighbors_seconds": neighbors_seconds,
                "fit_seconds": time.perf_counter() - started, "transform_seconds": 0.0}
        selected = np.arange(cells)
    else:
        from .ICGS import _import_umap_with_local_retry
        from .umap_fit import fit_umap
        umap = _import_umap_with_local_retry()
        model = umap.UMAP(n_neighbors=neighbors, min_dist=min_dist,
                          metric="euclidean", random_state=seed)
        pcs = graph.obsm["X_pca"]
        coordinates, info, selected = fit_umap(
            model, lambda rows: pcs[rows], cells, mode="landmark", labels=labels,
            max_fit_cells=fit_cells, batch_cells=batch_cells, seed=seed, log=log)
    if coordinates.shape != (cells, 2) or not np.isfinite(coordinates).all():
        raise ValueError("Accelerated UMAP must provide finite coordinates for every input cell.")
    info.update(source=f"scanpy_pca_{method}", pca_seconds=pca_seconds, n_pcs=components,
                pca_solver="arpack", pca_zero_center=True, n_features=features,
                n_neighbors=neighbors, min_dist=min_dist, metric="euclidean", random_state=seed)
    if log:
        log(f"UMAP fit completed in {info['fit_seconds']:.1f} s; transform in {info['transform_seconds']:.1f} s")
    return np.asarray(coordinates, dtype=np.float32), info, selected
