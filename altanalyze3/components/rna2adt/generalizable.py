"""Cluster-free per-ADT predictors designed to generalize aberrant expression.

The centroid-based models cannot detect a cell that aberrantly expresses a
surface protein because they predict each cell's ADT as the mean of normal
cells in its assigned state. The models in this module deliberately avoid
any cluster-conditional shortcut: they use per-cell RNA only, after a
within-dataset kNN-smoothing step that mitigates RNA dropout while
preserving per-cell deviations from the local manifold.

Three head variants are provided, all sharing the same preprocessing trunk:

* ``GeneralizableRidge``    -- per-ADT Ridge regression on smoothed canonical
                                partner genes (linear, fast).
* ``GeneralizableSplineGAM`` -- per-ADT Ridge on quantile-rank features +
                                 second-order Hermite expansion of the most
                                 informative partner gene; captures saturation.
* ``GeneralizableXGBoost``  -- per-ADT XGBoost regression on smoothed partner
                                genes; captures gene interactions.

All three predict per-cell ADT from per-cell RNA without any reference cell
state, so a Pro-B cell aberrantly expressing CD7 RNA will produce a high
predicted CD7 protein value; a normal Pro-B with no CD7 RNA will not. This
is the property required for surfacing disease biology.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import QuantileTransformer

from .adt_rna_map import audit_against_genes, load_curated_adt_rna_map


# --------------------------------------------------------------------------- #
# Preprocessing                                                                #
# --------------------------------------------------------------------------- #


@dataclass
class _SmoothingState:
    pca_components: Optional[np.ndarray]    # (n_pca, n_genes) used for query projection
    pca_mean: Optional[np.ndarray]          # (n_genes,)
    train_pca: Optional[np.ndarray]         # (n_train_cells, n_pca) — reference index
    train_smoothed: Optional[np.ndarray]    # (n_train_cells, n_genes) — for fit-time
    n_neighbors: int
    panel_indices: List[int]                # which columns of full RNA are smoothed
    panel_gene_names: List[str]


def _build_smoothing_state(
    X_train_full: np.ndarray,
    rna_genes: Sequence[str],
    panel_gene_names: Sequence[str],
    *,
    n_pca: int = 30,
    n_neighbors: int = 25,
) -> _SmoothingState:
    """Panel-only PCA + kNN smoothing of the panel-partner subset.

    Both the PCA basis and the smoothed gene subset are restricted to the
    canonical RNA partners of the ADT panel (~136 genes). This means:

      * Aberrant expression is preserved: a Pro-B cell with high CD7 RNA
        moves toward CD7+ neighbors in the panel-PCA space because the
        basis is dominated by panel signal, so its smoothed CD7 stays high.
        A whole-RNA basis would sort the cell back to Pro-B-like neighbors
        via 3000 non-panel genes and dilute its CD7 signal.
      * No leakage from non-panel RNA: the model truly only depends on
        the curated partner-gene relationships, in keeping with the
        "use ADT-to-gene lookups only" requirement.
    """
    rna_idx = {g: i for i, g in enumerate(rna_genes)}
    panel_indices = [rna_idx[g] for g in panel_gene_names if g in rna_idx]
    if not panel_indices:
        raise ValueError("No panel genes found in rna_genes")

    panel_train = X_train_full[:, panel_indices].astype(np.float32, copy=False)
    print(f"[smooth] n_train={X_train_full.shape[0]} n_panel={len(panel_indices)} "
          f"pca_on=panel pca={n_pca} k={n_neighbors}", flush=True)

    n_components = min(n_pca, panel_train.shape[1] - 1, panel_train.shape[0] - 1)
    n_components = max(2, n_components)
    svd = TruncatedSVD(n_components=n_components, random_state=0)
    Z_train = svd.fit_transform(panel_train).astype(np.float32)
    norms = np.linalg.norm(Z_train, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    Z_train_normed = Z_train / norms

    nn = NearestNeighbors(n_neighbors=n_neighbors, algorithm="auto", metric="euclidean")
    nn.fit(Z_train_normed)
    _, train_idx = nn.kneighbors(Z_train_normed, n_neighbors=n_neighbors)

    smoothed = np.zeros_like(panel_train)
    chunk = 4096
    for start in range(0, panel_train.shape[0], chunk):
        stop = min(start + chunk, panel_train.shape[0])
        nbrs = train_idx[start:stop]
        smoothed[start:stop] = panel_train[nbrs].mean(axis=1).astype(np.float32)

    return _SmoothingState(
        pca_components=svd.components_.astype(np.float32),  # (n_pca, n_panel)
        pca_mean=getattr(svd, "mean_", None) if hasattr(svd, "mean_") else None,
        train_pca=Z_train_normed,
        train_smoothed=smoothed,
        n_neighbors=n_neighbors,
        panel_indices=panel_indices,
        panel_gene_names=list(panel_gene_names),
    )


def _smooth_query(
    state: _SmoothingState,
    X_query_full: np.ndarray,
    *,
    self_smoothing_k: int = 5,
) -> np.ndarray:
    """Return the kNN-smoothed panel-gene expression of a query batch.

    The query is smoothed *against itself*, NOT against the training set,
    so that disease-driven aberrant expression is preserved: if a Pro-B
    cell aberrantly expresses CD7, smoothing among its query-set neighbors
    (likely also aberrant Pro-B cells in a diseased sample) keeps the CD7
    signal; smoothing against a normal-bone-marrow reference would average
    it out.

    Smoothing uses the same PCA basis the model was trained on, so the
    neighbor space is comparable to training. ``self_smoothing_k`` is set
    smaller than training-time k to keep variance in 'rare aberrant cell'
    estimates instead of regressing them to the local mean.
    """
    if state.pca_components is None:
        raise RuntimeError("smoothing state is unfit")
    Xq = X_query_full.astype(np.float32, copy=False)
    panel_q = Xq[:, state.panel_indices].astype(np.float32, copy=False)
    if panel_q.shape[0] <= 1 or self_smoothing_k <= 1:
        return panel_q

    # Project the panel subset through the panel-trained PCA basis, so the
    # query lives in the same neighbor space the model was fit on.
    Z = panel_q @ state.pca_components.T
    norms = np.linalg.norm(Z, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    Zn = Z / norms

    k = min(self_smoothing_k, max(2, Xq.shape[0]))
    nn = NearestNeighbors(n_neighbors=k, algorithm="auto", metric="euclidean")
    nn.fit(Zn)
    _, idx = nn.kneighbors(Zn, n_neighbors=k)
    out = np.zeros_like(panel_q)
    chunk = 4096
    for start in range(0, Xq.shape[0], chunk):
        stop = min(start + chunk, Xq.shape[0])
        nbrs = idx[start:stop]
        out[start:stop] = panel_q[nbrs].mean(axis=1).astype(np.float32)
    return out


# --------------------------------------------------------------------------- #
# Per-ADT heads                                                                #
# --------------------------------------------------------------------------- #


@dataclass
class _GenericHead:
    adt_name: str
    feature_indices: List[int]      # positions in the panel-partner matrix
    feature_names: List[str]
    model: object                   # sklearn-like with .predict(X)
    fallback_value: float


class _GeneralizableBase:
    """Common predictor scaffolding shared by Ridge / spline / XGBoost heads."""

    head_kind: str = "base"

    def __init__(
        self,
        *,
        rna_genes: Sequence[str],
        adt_names: Sequence[str],
        n_pca: int = 50,
        n_neighbors: int = 25,
    ) -> None:
        self.rna_genes = list(rna_genes)
        self.adt_names = list(adt_names)
        self.n_pca = int(n_pca)
        self.n_neighbors = int(n_neighbors)
        self.heads: List[_GenericHead] = []
        self._smoothing: Optional[_SmoothingState] = None
        self._panel_partner_names: List[str] = []
        self._adt_to_partners: Dict[str, List[str]] = {}
        self.n_features_in_: Optional[int] = None

    # --- Subclasses override --- #

    def _fit_head(self, Xj: np.ndarray, yj: np.ndarray):  # pragma: no cover
        raise NotImplementedError

    # --- Common pipeline --- #

    def _build_panel_union(self) -> List[str]:
        gene_idx = set(self.rna_genes)
        adt_to_partners = load_curated_adt_rna_map()
        audit = audit_against_genes(self.adt_names, self.rna_genes, map_obj=adt_to_partners)
        seen: set[str] = set()
        union: List[str] = []
        per_adt: Dict[str, List[str]] = {}
        for adt in self.adt_names:
            partners = [g for g in audit[adt]["matched"] if g in gene_idx]  # type: ignore[index]
            per_adt[adt] = partners
            for g in partners:
                if g not in seen:
                    seen.add(g)
                    union.append(g)
        self._adt_to_partners = per_adt
        return union

    def fit(self, X: np.ndarray, Y: np.ndarray) -> "_GeneralizableBase":
        self.n_features_in_ = int(X.shape[1])
        panel_union = self._build_panel_union()
        if not panel_union:
            raise ValueError("Empty panel-partner union; check curated map vs rna_genes")
        self._panel_partner_names = panel_union
        self._smoothing = _build_smoothing_state(
            X, self.rna_genes, panel_union,
            n_pca=self.n_pca, n_neighbors=self.n_neighbors,
        )
        Xs = self._smoothing.train_smoothed
        if Xs is None:
            raise RuntimeError("smoothing failed")

        gene_to_panel_idx = {g: i for i, g in enumerate(panel_union)}
        self.heads = []
        for j, adt in enumerate(self.adt_names):
            partner_names = self._adt_to_partners.get(adt, [])
            partner_idx = [gene_to_panel_idx[g] for g in partner_names if g in gene_to_panel_idx]
            yj = np.asarray(Y[:, j], dtype=np.float64)
            fallback = float(np.mean(yj))
            if not partner_idx:
                self.heads.append(_GenericHead(
                    adt_name=adt,
                    feature_indices=[],
                    feature_names=[],
                    model=_MeanFallback(fallback),
                    fallback_value=fallback,
                ))
                continue
            Xj = Xs[:, partner_idx]
            model = self._fit_head(Xj.astype(np.float64, copy=False), yj)
            self.heads.append(_GenericHead(
                adt_name=adt,
                feature_indices=list(partner_idx),
                feature_names=[panel_union[i] for i in partner_idx],
                model=model,
                fallback_value=fallback,
            ))
        return self

    def predict(self, X: np.ndarray, cluster_labels=None) -> np.ndarray:
        if self._smoothing is None:
            raise RuntimeError("model is unfit")
        Xs = _smooth_query(self._smoothing, np.asarray(X, dtype=np.float32))
        out = np.zeros((Xs.shape[0], len(self.adt_names)), dtype=np.float32)
        for j, head in enumerate(self.heads):
            if not head.feature_indices:
                out[:, j] = head.fallback_value
                continue
            Xj = Xs[:, head.feature_indices]
            preds = head.model.predict(Xj.astype(np.float64, copy=False))
            out[:, j] = np.asarray(preds, dtype=np.float32).ravel()
        return out


class _MeanFallback:
    def __init__(self, value: float) -> None:
        self.value = float(value)

    def predict(self, X) -> np.ndarray:
        return np.full(X.shape[0], self.value, dtype=np.float32)


class GeneralizableRidge(_GeneralizableBase):
    head_kind = "ridge"

    def __init__(self, *, alpha: float = 1.0, **kwargs) -> None:
        super().__init__(**kwargs)
        self.alpha = float(alpha)

    def _fit_head(self, Xj: np.ndarray, yj: np.ndarray):
        m = Ridge(alpha=self.alpha, random_state=0)
        m.fit(Xj, yj)
        return m


class GeneralizableElasticNet(_GeneralizableBase):
    """Per-ADT ElasticNet on smoothed canonical partner genes.

    Sparsity (L1) lets each head zero out partner genes that don't carry
    signal for that particular ADT — useful for multi-subunit ADTs where
    only one of the listed partners is actually informative (e.g. CD3D for
    Hu.CD3) and lets us hope to land between Ridge (overshoot) and Spline
    (sometimes undershoots saturation).
    """

    head_kind = "elastic"

    def __init__(self, *, alpha: float = 0.001, l1_ratio: float = 0.5, **kwargs) -> None:
        super().__init__(**kwargs)
        self.alpha = float(alpha)
        self.l1_ratio = float(l1_ratio)

    def _fit_head(self, Xj: np.ndarray, yj: np.ndarray):
        m = ElasticNet(alpha=self.alpha, l1_ratio=self.l1_ratio,
                       max_iter=2000, random_state=0)
        m.fit(Xj, yj)
        return m


class _RankElastic:
    """ElasticNet on per-feature quantile-rank transforms.

    Each input feature is mapped to its empirical CDF (rank in [0,1]) using
    the training distribution. ElasticNet then fits a linear model on the
    rank features. Because rank is bounded, even when the input RNA is
    artificially boosted to extreme values the rank input is clipped to 1,
    which kills the linear-extrapolation overshoot we saw with raw-feature
    Ridge/ElasticNet while keeping a faithful monotone response.
    """

    def __init__(self, alpha: float = 0.001, l1_ratio: float = 0.5,
                 n_quantiles: int = 1000) -> None:
        self.alpha = float(alpha)
        self.l1_ratio = float(l1_ratio)
        self.n_quantiles = int(n_quantiles)
        self._qt: Optional[QuantileTransformer] = None
        self._enet: Optional[ElasticNet] = None

    def fit(self, X: np.ndarray, y: np.ndarray) -> "_RankElastic":
        X = np.asarray(X, dtype=np.float64)
        n_q = min(self.n_quantiles, max(2, X.shape[0]))
        self._qt = QuantileTransformer(
            n_quantiles=n_q, output_distribution="uniform",
            subsample=min(100_000, X.shape[0]), random_state=0,
        )
        Xr = self._qt.fit_transform(X)
        self._enet = ElasticNet(alpha=self.alpha, l1_ratio=self.l1_ratio,
                                max_iter=2000, random_state=0)
        self._enet.fit(Xr, y)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        X = np.asarray(X, dtype=np.float64)
        Xr = self._qt.transform(X)
        return self._enet.predict(Xr)


class GeneralizableElasticRank(_GeneralizableBase):
    head_kind = "elastic_rank"

    def __init__(self, *, alpha: float = 0.001, l1_ratio: float = 0.5, **kwargs) -> None:
        super().__init__(**kwargs)
        self.alpha = float(alpha)
        self.l1_ratio = float(l1_ratio)

    def _fit_head(self, Xj: np.ndarray, yj: np.ndarray):
        m = _RankElastic(alpha=self.alpha, l1_ratio=self.l1_ratio)
        m.fit(Xj, yj)
        return m


class _SplineElastic:
    """ElasticNet on the spline-expanded feature basis.

    Reuses the same Hermite-style expansion as ``_SplineRidge`` (linear,
    soft-saturation tanh, quadratic) but fits with ElasticNet so L1 can
    deselect basis functions that don't carry signal. The tanh term keeps
    extrapolation bounded.
    """

    def __init__(self, alpha: float = 0.001, l1_ratio: float = 0.5) -> None:
        self.alpha = float(alpha)
        self.l1_ratio = float(l1_ratio)
        self._scales: Optional[np.ndarray] = None
        self._means: Optional[np.ndarray] = None
        self._enet: Optional[ElasticNet] = None

    def _expand(self, X: np.ndarray) -> np.ndarray:
        if self._scales is None or self._means is None:
            raise RuntimeError("spline_elastic model unfit")
        x_centered = X - self._means
        sat = np.tanh(X / self._scales)
        sq = x_centered ** 2
        return np.concatenate([X, sat, sq], axis=1)

    def fit(self, X: np.ndarray, y: np.ndarray) -> "_SplineElastic":
        X = np.asarray(X, dtype=np.float64)
        self._scales = np.maximum(np.std(X, axis=0), 1e-6).astype(np.float64)
        self._means = np.mean(X, axis=0).astype(np.float64)
        Xp = self._expand(X)
        self._enet = ElasticNet(alpha=self.alpha, l1_ratio=self.l1_ratio,
                                max_iter=2000, random_state=0)
        self._enet.fit(Xp, y)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        X = np.asarray(X, dtype=np.float64)
        return self._enet.predict(self._expand(X))


class GeneralizableElasticSpline(_GeneralizableBase):
    head_kind = "elastic_spline"

    def __init__(self, *, alpha: float = 0.001, l1_ratio: float = 0.5, **kwargs) -> None:
        super().__init__(**kwargs)
        self.alpha = float(alpha)
        self.l1_ratio = float(l1_ratio)

    def _fit_head(self, Xj: np.ndarray, yj: np.ndarray):
        m = _SplineElastic(alpha=self.alpha, l1_ratio=self.l1_ratio)
        m.fit(Xj, yj)
        return m


class _SplineRidge:
    """Ridge on a Hermite-style spline expansion of each input feature.

    Adds for each feature x:
        x, x^2 (centred), and a soft saturation term tanh(x / s)
    where s is the per-feature scale. Captures saturation that's typical
    of RNA->protein dose-response (high RNA + plateau in protein expression)
    while staying linear in the spline coefficients.
    """

    def __init__(self, alpha: float = 1.0) -> None:
        self.alpha = float(alpha)
        self._scales: Optional[np.ndarray] = None
        self._means: Optional[np.ndarray] = None
        self._ridge: Optional[Ridge] = None

    def _expand(self, X: np.ndarray) -> np.ndarray:
        if self._scales is None or self._means is None:
            raise RuntimeError("spline model unfit")
        x_centered = X - self._means
        sat = np.tanh(X / self._scales)
        sq = x_centered ** 2
        return np.concatenate([X, sat, sq], axis=1)

    def fit(self, X: np.ndarray, y: np.ndarray) -> "_SplineRidge":
        X = np.asarray(X, dtype=np.float64)
        self._scales = np.maximum(np.std(X, axis=0), 1e-6).astype(np.float64)
        self._means = np.mean(X, axis=0).astype(np.float64)
        Xp = self._expand(X)
        self._ridge = Ridge(alpha=self.alpha, random_state=0)
        self._ridge.fit(Xp, y)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        X = np.asarray(X, dtype=np.float64)
        return self._ridge.predict(self._expand(X))


class GeneralizableSplineGAM(_GeneralizableBase):
    head_kind = "spline_gam"

    def __init__(self, *, alpha: float = 1.0, **kwargs) -> None:
        super().__init__(**kwargs)
        self.alpha = float(alpha)

    def _fit_head(self, Xj: np.ndarray, yj: np.ndarray):
        m = _SplineRidge(alpha=self.alpha)
        m.fit(Xj, yj)
        return m


class GeneralizableXGBoost(_GeneralizableBase):
    head_kind = "xgboost"

    def __init__(
        self,
        *,
        n_estimators: int = 200,
        max_depth: int = 4,
        learning_rate: float = 0.08,
        min_child_weight: float = 5.0,
        subsample: float = 0.8,
        colsample_bytree: float = 0.8,
        reg_lambda: float = 1.0,
        n_jobs: int = 4,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs)
        self.n_estimators = int(n_estimators)
        self.max_depth = int(max_depth)
        self.learning_rate = float(learning_rate)
        self.min_child_weight = float(min_child_weight)
        self.subsample = float(subsample)
        self.colsample_bytree = float(colsample_bytree)
        self.reg_lambda = float(reg_lambda)
        self.n_jobs = int(n_jobs)

    def _fit_head(self, Xj: np.ndarray, yj: np.ndarray):
        from xgboost import XGBRegressor
        # On macOS arm64 the hist method's thread pool can segfault when many
        # boosters are fit back to back. Force single-threaded fit and the
        # exact tree method, which is acceptably fast for ~25k rows by ~4
        # features per ADT head.
        m = XGBRegressor(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            learning_rate=self.learning_rate,
            min_child_weight=self.min_child_weight,
            subsample=self.subsample,
            colsample_bytree=self.colsample_bytree,
            reg_lambda=self.reg_lambda,
            tree_method="exact",
            n_jobs=1,
            random_state=0,
            objective="reg:squarederror",
            verbosity=0,
        )
        m.fit(Xj, yj)
        return m
