"""Per-ADT targeted regression using only canonical RNA partners.

Each ADT gets its own small regressor trained on a tight set of mapped genes
(typically 1–4). With so few features the model cannot memorize donor batch
effects through irrelevant signal, which is the dominant failure mode of the
global multi-output models on leave-one-donor-out splits.

Two approaches are supported:

* ``ridge``  — closed-form per-protein Ridge regression. Fastest and the
  recommended default.
* ``elastic`` — per-protein ElasticNet with a fixed ``alpha``. Slightly more
  conservative for genes with co-expressed partners.

Either way the model conforms to the ``rna2adt.api.Rna2AdtBundle`` protocol:
it exposes a ``.predict(X)`` method that takes a (n_cells, n_input_genes)
matrix aligned to ``X_columns``. ``X_columns`` here is the union of all
matched RNA partners across the panel (so the bundle's input alignment logic
in ``api`` works unchanged), but the per-protein head only consumes the
columns it was fit on.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNet, Ridge


@dataclass
class _ProteinHead:
    adt_name: str
    feature_indices: List[int]      # positions in the unioned X_columns
    feature_names: List[str]
    coef: np.ndarray                # shape (len(feature_indices),)
    intercept: float
    fallback_value: float           # used if all features missing or isotype


class TargetedAdtRegressor:
    """Per-ADT linear regressor on a sparse (ADT -> RNA partners) graph."""

    def __init__(
        self,
        *,
        feature_names: Sequence[str],
        adt_names: Sequence[str],
        adt_to_features: Dict[str, List[str]],
        method: str = "ridge",
        alpha: float = 1.0,
        l1_ratio: float = 0.5,
        share_panel_context: bool = False,
    ) -> None:
        self.feature_names = list(feature_names)
        self._feature_idx = {name: i for i, name in enumerate(self.feature_names)}
        self.adt_names = list(adt_names)
        self.adt_to_features = {k: list(v) for k, v in adt_to_features.items()}
        self.method = method
        self.alpha = float(alpha)
        self.l1_ratio = float(l1_ratio)
        self.share_panel_context = bool(share_panel_context)
        self.heads: List[_ProteinHead] = []
        self.n_features_in_: Optional[int] = None

    def fit(self, X: np.ndarray, Y: np.ndarray) -> "TargetedAdtRegressor":
        if X.shape[1] != len(self.feature_names):
            raise ValueError(
                f"X has {X.shape[1]} columns, expected {len(self.feature_names)} (matches feature_names)"
            )
        if Y.shape[1] != len(self.adt_names):
            raise ValueError(
                f"Y has {Y.shape[1]} columns, expected {len(self.adt_names)} (matches adt_names)"
            )
        self.n_features_in_ = int(X.shape[1])
        self.heads = []
        for j, adt in enumerate(self.adt_names):
            partner_names = self.adt_to_features.get(adt, [])
            partner_idx = [self._feature_idx[n] for n in partner_names if n in self._feature_idx]
            if self.share_panel_context:
                extra_idx = [i for i in range(len(self.feature_names)) if i not in set(partner_idx)]
                fit_idx = partner_idx + extra_idx
            else:
                fit_idx = partner_idx
            yj = np.asarray(Y[:, j], dtype=np.float64)
            fallback = float(np.mean(yj))
            if not partner_idx:
                # No partners present — head returns the training mean.
                self.heads.append(_ProteinHead(
                    adt_name=adt,
                    feature_indices=[],
                    feature_names=[],
                    coef=np.zeros(0, dtype=np.float64),
                    intercept=fallback,
                    fallback_value=fallback,
                ))
                continue
            Xj = X[:, fit_idx].astype(np.float64, copy=False)
            if self.method == "elastic":
                model = ElasticNet(alpha=self.alpha, l1_ratio=self.l1_ratio,
                                   max_iter=2000, random_state=0)
            else:
                model = Ridge(alpha=self.alpha, random_state=0)
            model.fit(Xj, yj)
            coef = np.asarray(model.coef_, dtype=np.float64).ravel()
            intercept = float(model.intercept_)
            kept_names = [self.feature_names[i] for i in fit_idx]
            self.heads.append(_ProteinHead(
                adt_name=adt,
                feature_indices=list(fit_idx),
                feature_names=kept_names,
                coef=coef,
                intercept=intercept,
                fallback_value=fallback,
            ))
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        X = np.asarray(X)
        n = X.shape[0]
        out = np.zeros((n, len(self.heads)), dtype=np.float32)
        for j, head in enumerate(self.heads):
            if not head.feature_indices:
                out[:, j] = head.fallback_value
                continue
            Xj = X[:, head.feature_indices].astype(np.float64, copy=False)
            preds = Xj @ head.coef + head.intercept
            out[:, j] = preds.astype(np.float32)
        return out

    def head_summary(self) -> pd.DataFrame:
        rows = []
        for head in self.heads:
            rows.append({
                "adt_name": head.adt_name,
                "n_features": len(head.feature_indices),
                "features": ",".join(head.feature_names) if head.feature_names else "(fallback=mean)",
                "intercept": head.intercept,
                "fallback_value": head.fallback_value,
                "coef_l1": float(np.abs(head.coef).sum()),
            })
        return pd.DataFrame(rows)
