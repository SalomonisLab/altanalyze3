"""CITE-seq to flow label transfer. Several methods, one interface, so they can be compared.

Direction is always CITE-seq -> flow: CITE-seq carries the labels, flow events receive them.

Each method takes (cite X, flow X) over the SHARED markers plus the CITE labels, and returns a
label per flow event. Nothing here decides which method is best; `evaluate.py` does that from
measurements.
"""
from __future__ import annotations

import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import LabelEncoder

from .normalize import scale_feature, reference_spline, kde_quantile_map

__all__ = ["METHODS", "kde_knn", "zscore_knn", "minmax_knn", "harmony_knn", "xgboost_transfer"]


def _scale_both(cite, flow, lo=1.0, hi=99.0):
    """Percentile-clip + MinMax each feature of each platform independently."""
    C = np.column_stack([scale_feature(cite[:, j], lo, hi) for j in range(cite.shape[1])])
    F = np.column_stack([scale_feature(flow[:, j], lo, hi) for j in range(flow.shape[1])])
    return C, F


def kde_knn(cite, flow, labels, k=15, ties="first", lo=1.0, hi=99.0, seed=0, n_ref=20000):
    """The published method: KDE quantile-map CITE onto the flow distribution, then kNN.

    The reference distribution is the FLOW channel, because the labels must land in flow space.
    """
    C, F = _scale_both(cite, flow, lo, hi)
    rng = np.random.default_rng(seed)
    Cm = np.empty_like(C)
    for j in range(C.shape[1]):
        ref = F[:, j]
        if ref.size > n_ref:                      # gaussian_kde is O(n^2) in evaluation
            ref = ref[rng.choice(ref.size, n_ref, replace=False)]
        spl = reference_spline(ref)
        Cm[:, j] = kde_quantile_map(C[:, j], None, spline=spl, ties=ties)
    clf = KNeighborsClassifier(n_neighbors=k, n_jobs=-1).fit(Cm, labels)
    return clf.predict(F)


def zscore_knn(cite, flow, labels, k=15, **kw):
    """Per-platform z-score, then kNN. The simplest scale alignment that can work."""
    z = lambda M: (M - M.mean(0)) / np.where(M.std(0) == 0, 1, M.std(0))
    clf = KNeighborsClassifier(n_neighbors=k, n_jobs=-1).fit(z(cite), labels)
    return clf.predict(z(flow))


def minmax_knn(cite, flow, labels, k=15, lo=1.0, hi=99.0, **kw):
    """Percentile-clipped MinMax only, no distribution matching. The control for KDE."""
    C, F = _scale_both(cite, flow, lo, hi)
    clf = KNeighborsClassifier(n_neighbors=k, n_jobs=-1).fit(C, labels)
    return clf.predict(F)


def harmony_knn(cite, flow, labels, k=15, lo=1.0, hi=99.0, seed=0, **kw):
    """Harmony batch-corrects the two platforms jointly, then kNN in the corrected space."""
    import harmonypy as hm
    import pandas as pd
    C, F = _scale_both(cite, flow, lo, hi)
    X = np.vstack([C, F])
    meta = pd.DataFrame({"platform": ["cite"] * len(C) + ["flow"] * len(F)})
    ho = hm.run_harmony(X, meta, ["platform"], max_iter_harmony=10, random_state=seed)
    Z = np.asarray(ho.Z_corr)
    # harmonypy returns (d, N) in some versions and (N, d) in others. Pick the orientation
    # whose row count equals the number of cells rather than assuming either one.
    if Z.shape[0] != X.shape[0]:
        Z = Z.T
    if Z.shape[0] != X.shape[0]:
        raise ValueError("harmony returned %s for %d cells" % (Z.shape, X.shape[0]))
    clf = KNeighborsClassifier(n_neighbors=k, n_jobs=-1).fit(Z[: len(C)], labels)
    return clf.predict(Z[len(C):])


def xgboost_transfer(cite, flow, labels, lo=1.0, hi=99.0, seed=0, **kw):
    """Train a classifier on CITE-seq, predict flow. No distribution matching at all."""
    from xgboost import XGBClassifier
    C, F = _scale_both(cite, flow, lo, hi)
    le = LabelEncoder().fit(labels)
    # n_jobs=1: xgboost's OpenMP runtime clashes with the one already loaded by numpy/sklearn
    # on macOS and the process dies without a traceback. Single-threaded is slower and survives.
    clf = XGBClassifier(n_estimators=200, max_depth=6, learning_rate=0.2, tree_method="hist",
                        n_jobs=1, random_state=seed, verbosity=0)
    clf.fit(C, le.transform(labels))
    return le.inverse_transform(clf.predict(F))


def cellharmony_transfer(*args, **kwargs):
    from .cellharmony import kde_cellharmony
    return kde_cellharmony(*args, **kwargs)


METHODS = {
    'kde_cellharmony': cellharmony_transfer,
    "kde_knn": kde_knn,
    "kde_knn_avgties": lambda *a, **k: kde_knn(*a, ties="average", **k),
    "zscore_knn": zscore_knn,
    "minmax_knn": minmax_knn,
    "harmony_knn": harmony_knn,
    "xgboost": xgboost_transfer,
}
