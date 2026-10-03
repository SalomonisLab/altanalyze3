"""Cross-platform feature normalization for CITE-seq to flow-cytometry mapping.

The published method (MultiLin_Project_Code_Repository, 3-CITE-seq_InfinityFlow_Integration,
kde_mapping_normalization.ipynb) maps each feature by rank onto a KDE-smoothed reference CDF:

  1. scale_feature        percentile-clip, then MinMax to [0, 1]
  2. reference spline     gaussian_kde -> integrate to a CDF -> monotone k=1 spline
  3. map by rank          rank the query, evaluate the spline at (r + 1) / (n + 2)

Compatibility with the notebook requires integration='published' and ties='published'.
The defaults use trapezoid integration and stable ranks: these are approximations, especially
for tied detector floors and clipped values. The older reference helper uses stable ranks
and is not an independent reproduction of the notebook's tie behavior.

Why the reference implementation is slow, measured rather than assumed:
  * the spline is evaluated one element at a time in a Python list comprehension
  * ranks come from an argsort plus a pandas Series sort, two O(n log n) passes and a copy
  * `integrate_box` is called once per grid point, each a fresh quadrature
"""
from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.interpolate import InterpolatedUnivariateSpline

__all__ = [
    "scale_feature",
    "reference_spline",
    "kde_quantile_map_reference",
    "kde_quantile_map",
    "rank_plotting_positions",
]


def scale_feature(x, min_pct: float = 0.0, max_pct: float = 100.0) -> np.ndarray:
    """Percentile-clip then MinMax to [0, 1]. Matches the published scale_feature."""
    a = np.asarray(x, dtype=np.float64).ravel().copy()
    if a.size == 0:
        return a
    lo = np.percentile(a, min_pct)
    hi = np.percentile(a, max_pct)
    np.clip(a, lo, hi, out=a)
    span = a.max() - a.min()
    if span <= 0:
        # A constant feature carries no rank information. Return zeros rather than NaN, and
        # let the caller decide whether to drop it; silently emitting NaN would poison the
        # downstream distance metric.
        return np.zeros_like(a)
    return (a - a.min()) / span


def rank_plotting_positions(x, ties: str = "first") -> np.ndarray:
    """(rank + 1) / (n + 2) for each element.

    ties='published' uses numpy's default argsort, as in Kyle's notebook. Its ordering
                   among ties depends on the NumPy version and is not stable.
    ties='first'   equal values receive DIFFERENT ranks, ordered by position,
                   so two cells with an identical marker value map to
                   different outputs. Percentile clipping in `scale_feature` manufactures such
                   ties at both bounds, and a flow detector floor manufactures many more.
    ties='average' gives equal values the same rank, so identical inputs map identically.

    Inverting the argsort permutation avoids the notebook's extra pandas sort.
    Use ties='published' to preserve its default (unstable) sorting convention.
    """
    a = np.asarray(x).ravel()
    n = a.size
    if ties == "average":
        from scipy.stats import rankdata
        ranks = rankdata(a, method="average") - 1.0
    elif ties in {"first", "published"}:
        order = np.argsort(a, kind="stable") if ties == "first" else np.argsort(a)
        ranks = np.empty(n, dtype=np.float64)
        ranks[order] = np.arange(n, dtype=np.float64)
    else:
        raise ValueError("ties must be 'first', 'published' or 'average', got %r" % (ties,))
    return (ranks + 1.0) / (n + 2.0)


def reference_spline(reference_signal, n_grid: int = 100, integration: str = "trapezoid"):
    """Monotone percentile -> value spline over a gaussian_kde of the reference.

    `reference_signal` must already sit in [0, 1] (i.e. be scale_feature output), because the
    integration bounds are 0 and 1, exactly as the published code assumes.
    """
    ref = np.asarray(reference_signal, dtype=np.float64).ravel()
    kde = stats.gaussian_kde(ref)
    grid = np.linspace(0.0, 1.0, n_grid)
    total = kde.integrate_box_1d(0.0, 1.0)
    # Trapezoid integration is an approximation; published mode preserves exact integrals.
    if integration == "published":
        # Exact one-dimensional Gaussian integrals, as used by the notebook.
        cum = np.array([kde.integrate_box_1d(0.0, g) for g in grid])
    elif integration == "trapezoid":
        dens = kde(grid)
        cum = np.concatenate(([0.0], np.cumsum((dens[1:] + dens[:-1]) * 0.5 * np.diff(grid))))
    else:
        raise ValueError("integration must be 'published' or 'trapezoid'")
    pct = cum / total if total > 0 else cum
    keep = (pct > 0) & (pct < 1)
    xs, ys = pct[keep], grid[keep]
    _, uniq = np.unique(xs, return_index=True)
    xs, ys = xs[uniq], ys[uniq]
    xs = np.concatenate(([0.0], xs, [1.0]))
    ys = np.concatenate(([0.0], ys, [1.0]))
    return InterpolatedUnivariateSpline(xs, ys, k=1, bbox=[0.0, 1.0], ext=3)


def kde_quantile_map_reference(query_signal, reference_signal, n_grid: int = 100) -> np.ndarray:
    """Exact KDE integrals with stable ranks; differs from the notebook at ties."""
    kde = stats.gaussian_kde(np.asarray(reference_signal, dtype=np.float64).ravel())
    total = kde.integrate_box_1d(0.0, 1.0)
    grid = np.linspace(0.0, 1.0, n_grid)
    pct = np.array([kde.integrate_box_1d(0.0, g) / total for g in grid])
    keep = (pct > 0) & (pct < 1)
    xs, ys = pct[keep], grid[keep]
    _, uniq = np.unique(xs, return_index=True)
    xs, ys = xs[uniq], ys[uniq]
    xs = np.concatenate(([0.0], xs, [1.0]))
    ys = np.concatenate(([0.0], ys, [1.0]))
    spline = InterpolatedUnivariateSpline(xs, ys, k=1, bbox=[0.0, 1.0], ext=3)
    pos = rank_plotting_positions(query_signal)
    return np.array([float(spline(p)) for p in pos])


def kde_quantile_map(query_signal, reference_signal, n_grid: int = 100,
                     spline=None, ties: str = "first", integration: str = "trapezoid") -> np.ndarray:
    """Rank-based KDE mapping, with explicit integration and tie conventions.

    Pass `spline` to reuse one reference fit across many query batches.
    """
    if spline is None:
        spline = reference_spline(reference_signal, n_grid=n_grid, integration=integration)
    return np.asarray(spline(rank_plotting_positions(query_signal, ties=ties)), dtype=np.float64)
