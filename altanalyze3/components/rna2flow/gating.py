"""MarkerFinder-derived virtual flow gates, then iterative optimization.

A FlowJo gate is drawn by eye on two channels. Here the two channels are CHOSEN by
MarkerFinder: for a target population, the channels with the highest point-biserial rho are
the ones that separate it, which is the same statistic altanalyze3 uses to call markers.

The gate is then optimized rather than guessed. Coordinate ascent moves one rectangle bound at
a time to maximize F1 against population membership, which is the quantity a gate is actually
for: capture the population (recall) without admitting others (precision).

This extends the published KDE pipeline. That pipeline maps and transfers labels; it stops
there. Turning a transferred label into a reproducible gate, and scoring that gate, is the
step added here.
"""
from __future__ import annotations

import numpy as np

from ..cellHarmony.markerFinder import marker_finder

__all__ = ["rank_channels_for_population", "propose_gate", "optimize_gate", "gate_population",
           "gate_all_populations"]


def _f1(sel, truth):
    tp = float(np.sum(sel & truth))
    if tp == 0:
        return 0.0, 0.0, 0.0
    prec = tp / float(np.sum(sel))
    rec = tp / float(np.sum(truth))
    return 2 * prec * rec / (prec + rec), prec, rec


def rank_channels_for_population(X, labels, channels, population, min_cells=20):
    """MarkerFinder rho for every channel against this population's 0/1 indicator."""
    labels = np.asarray(labels).astype(str)
    keep = {p for p, n in zip(*np.unique(labels, return_counts=True)) if n >= min_cells}
    m = np.isin(labels, list(keep))
    if population not in keep or m.sum() == 0:
        return None
    r_df, _ = marker_finder(np.asarray(X, dtype=np.float64)[m], list(labels[m]),
                            gene_names=list(channels), validate_scaling=False)
    if population not in r_df.columns:
        return None
    return r_df[population].sort_values(ascending=False)


def propose_gate(x, y, truth, q=(10.0, 90.0)):
    """A starting rectangle: the central quantile box of the target population."""
    return [float(np.percentile(x[truth], q[0])), float(np.percentile(x[truth], q[1])),
            float(np.percentile(y[truth], q[0])), float(np.percentile(y[truth], q[1]))]


def optimize_gate(x, y, truth, box, rounds=6, steps=12):
    """Coordinate ascent on the four bounds, maximizing F1 against population membership."""
    box = list(box)
    best = _f1((x >= box[0]) & (x <= box[1]) & (y >= box[2]) & (y <= box[3]), truth)[0]
    for _ in range(rounds):
        improved = False
        for i, arr in ((0, x), (1, x), (2, y), (3, y)):
            lo, hi = np.percentile(arr, 0.5), np.percentile(arr, 99.5)
            for cand in np.linspace(lo, hi, steps):
                trial = list(box)
                trial[i] = float(cand)
                if trial[0] >= trial[1] or trial[2] >= trial[3]:
                    continue
                sel = (x >= trial[0]) & (x <= trial[1]) & (y >= trial[2]) & (y <= trial[3])
                f = _f1(sel, truth)[0]
                if f > best + 1e-6:
                    best, box, improved = f, trial, True
        if not improved:
            break
    return box, best


def gate_population(X, labels, channels, population, min_cells=20, rounds=6):
    """Pick two channels by MarkerFinder, propose a gate, optimize it, score it."""
    rho = rank_channels_for_population(X, labels, channels, population, min_cells)
    if rho is None or len(rho) < 2:
        return None
    cx, cy = rho.index[0], rho.index[1]
    ix, iy = channels.index(cx), channels.index(cy)
    x, y = np.asarray(X[:, ix], float), np.asarray(X[:, iy], float)
    truth = np.asarray(labels).astype(str) == population
    if truth.sum() < min_cells:
        return None
    box0 = propose_gate(x, y, truth)
    f0 = _f1((x >= box0[0]) & (x <= box0[1]) & (y >= box0[2]) & (y <= box0[3]), truth)[0]
    box, f1 = optimize_gate(x, y, truth, box0, rounds=rounds)
    sel = (x >= box[0]) & (x <= box[1]) & (y >= box[2]) & (y <= box[3])
    f, prec, rec = _f1(sel, truth)
    return {"population": population, "x_channel": cx, "y_channel": cy,
            "rho_x": round(float(rho.iloc[0]), 4), "rho_y": round(float(rho.iloc[1]), 4),
            "gate": [round(v, 5) for v in box], "f1_initial": round(float(f0), 4),
            "f1": round(float(f), 4), "precision": round(float(prec), 4),
            "recall": round(float(rec), 4), "n_population": int(truth.sum()),
            "n_gated": int(sel.sum()),
            "polygon": [[box[0], box[2]], [box[1], box[2]], [box[1], box[3]], [box[0], box[3]]]}


def gate_all_populations(X, labels, channels, min_cells=20, rounds=6):
    labels = np.asarray(labels).astype(str)
    out = []
    for p, n in zip(*np.unique(labels, return_counts=True)):
        if n < min_cells or p == "unassigned":
            continue
        g = gate_population(X, labels, channels, p, min_cells, rounds)
        if g:
            out.append(g)
    return sorted(out, key=lambda d: -d["f1"])
