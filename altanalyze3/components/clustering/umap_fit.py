"""Opt-in landmark UMAP with bounded input blocks and coordinates for every cell."""
from __future__ import annotations

import math
import time

import numpy as np


def select_landmarks(labels, budget, seed, min_per_state=200):
    """Stratify within states, preserving minimum coverage and rare states."""
    labels = np.asarray(labels, dtype=str)
    rng = np.random.default_rng(seed)
    _, inverse, counts = np.unique(labels, return_inverse=True, return_counts=True)
    # Sorting once avoids scanning the entire roster for each state.
    order = np.argsort(inverse, kind="stable")
    groups = np.split(order, np.cumsum(counts)[:-1])
    quota = np.minimum(counts, min_per_state)
    budget = min(len(labels), max(int(budget), int(quota.sum())))
    available = counts - quota
    extra = budget - int(quota.sum())
    if extra:
        shares = available * (extra / int(available.sum()))
        allocation = np.floor(shares).astype(np.int64)
        left = extra - int(allocation.sum())
        # Largest-remainder allocation fills the exact budget without exceeding
        # a state's population. Sampling remains within each state, not global.
        order = np.argsort(-(shares - allocation), kind="stable")
        allocation[order[:left]] += 1
        quota += allocation
    selected = np.concatenate([rng.choice(rows, int(n), replace=False)
                               for rows, n in zip(groups, quota)])
    return np.sort(selected)


def fit_umap(model, load_rows, n_cells, *, mode="full", labels=None,
             max_fit_cells=30000, batch_cells=50000, min_per_state=200, seed=0, log=None):
    """Keep the full-fit call unchanged; landmark mode only changes the embedding.

    load_rows receives ordered integer row positions and returns the baseline's
    feature matrix for those cells. No feature selection or clustering occurs here.
    """
    if mode not in {"full", "landmark"}:
        raise ValueError(f"Unknown UMAP fit mode: {mode}")
    if max_fit_cells < 3 or batch_cells < 1 or min_per_state < 1:
        raise ValueError("UMAP landmark and batch limits must be positive (at least 3 landmarks).")
    rows = np.arange(n_cells)
    selected = rows
    if mode == "landmark" and n_cells > max_fit_cells:
        if labels is None or len(labels) != n_cells:
            raise ValueError("Landmark UMAP requires a state label for every input cell.")
        budget = max(max_fit_cells, int(model.n_neighbors) + 1)
        selected = select_landmarks(labels, budget, seed, min_per_state)
    effective_mode = "landmark" if len(selected) < n_cells else "full"
    if log:
        log(f"UMAP fit mode={effective_mode} (requested={mode}); fitting {len(selected):,} "
            f"of {n_cells:,} cells; all cells receive coordinates")
        if effective_mode == "landmark":
            log(f"UMAP landmarks: minimum {min_per_state} cells per state; all cells in smaller "
                "states; remaining places sampled proportionally within larger states")
    started = time.perf_counter()
    fitted = model.fit_transform(load_rows(selected))
    fit_seconds = time.perf_counter() - started
    if log:
        log(f"UMAP fitting finished in {fit_seconds:.1f} s")
    transform_seconds = 0.0
    if effective_mode == "full":
        coordinates = fitted
    else:
        coordinates = np.empty((n_cells, 2), dtype=np.float32)
        coordinates[selected] = fitted
        remaining = np.ones(n_cells, dtype=bool)
        remaining[selected] = False
        to_transform = rows[remaining]
        started = time.perf_counter()
        # Balance blocks so a tiny last batch does not get extra UMAP epochs.
        for block in np.array_split(to_transform, math.ceil(len(to_transform) / batch_cells)):
            if log:
                log(f"UMAP mapping {len(block):,} remaining cells")
            coordinates[block] = model.transform(load_rows(block))
            if log:
                log(f"UMAP transformed {len(block):,} remaining cells")
        transform_seconds = time.perf_counter() - started
    if coordinates.shape != (n_cells, 2) or not np.isfinite(coordinates).all():
        raise ValueError("UMAP must return finite coordinates for the complete input cell roster.")
    info = {"requested_fit_mode": mode, "fit_mode": effective_mode,
            "fit_cells": len(selected), "total_cells": n_cells,
            "max_fit_cells": max_fit_cells, "transform_batch_cells": batch_cells,
            "min_landmarks_per_state": min_per_state, "random_state": seed,
            "landmark_strategy": "cluster_stratified_proportional" if effective_mode == "landmark" else "all_cells",
            "fit_seconds": fit_seconds, "transform_seconds": transform_seconds}
    if log:
        log(f"UMAP fit completed in {fit_seconds:.1f} s; transform in {transform_seconds:.1f} s")
    return coordinates, info, selected
