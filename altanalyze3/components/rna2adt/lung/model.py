"""The human lung rna2adt model class.

Kept in its own module on purpose. ``train_lung`` runs as ``python -m ...``,
which makes it ``__main__``; a class defined there pickles as
``__main__.HumanLungRna2AdtModel`` and then fails to unpickle anywhere else,
including inside cellHarmony-web. Defining it here pickles the class as
``altanalyze3.components.rna2adt.lung.model.HumanLungRna2AdtModel``, which
loads from any process.
"""

from __future__ import annotations

from typing import List

import numpy as np
from sklearn.linear_model import ElasticNet
from sklearn.preprocessing import StandardScaler


class HumanLungRna2AdtModel:
    """Per-ADT ElasticNet on the panel-union of human lung whitelist genes.

    Architecturally identical to ``components.rna2adt.mouse.train_mouse.
    MouseRna2AdtModel`` and to the human bone marrow
    ``Rna2LipidArchPanelPerProtein(feature_source="whitelist_union")``.
    """

    head_kind = "rna2adt_hs_lung_panel_per_protein_whitelist_union"

    def __init__(self, *, alpha: float = 0.01, l1_ratio: float = 0.5,
                 max_iter: int = 2000) -> None:
        self.alpha = float(alpha)
        self.l1_ratio = float(l1_ratio)
        self.max_iter = int(max_iter)
        self.rna_genes: List[str] = []
        self.adt_names: List[str] = []
        self.scaler_x = None
        self.scaler_y = None
        self.models: List[ElasticNet] = []
        self._coef_matrix = None
        self._intercept_vector = None
        self.n_features_in_ = None

    def fit(self, X: np.ndarray, Y: np.ndarray, *, rna_genes, adt_names) -> "HumanLungRna2AdtModel":
        self.rna_genes = list(rna_genes)
        self.adt_names = list(adt_names)
        self.n_features_in_ = int(X.shape[1])
        Xr = np.asarray(X, dtype=np.float32)
        Yr = np.asarray(Y, dtype=np.float32)
        self.scaler_x = StandardScaler()
        Xs = self.scaler_x.fit_transform(Xr)
        self.scaler_y = StandardScaler()
        Ys = self.scaler_y.fit_transform(Yr)
        print(f"[train] fitting {Ys.shape[1]} per-ADT ElasticNets alpha={self.alpha} "
              f"l1_ratio={self.l1_ratio} n_features={Xs.shape[1]} n_cells={Xs.shape[0]}",
              flush=True)
        self.models = []
        self._coef_matrix = None
        self._intercept_vector = None
        for j in range(Ys.shape[1]):
            head = ElasticNet(alpha=self.alpha, l1_ratio=self.l1_ratio,
                              max_iter=self.max_iter, tol=1e-3, random_state=0)
            head.fit(Xs, Ys[:, j])
            self.models.append(head)
        print(f"[train] done {len(self.models)} heads", flush=True)
        return self

    def predict(self, X: np.ndarray, cluster_labels=None) -> np.ndarray:
        if not self.models:
            raise RuntimeError("model unfit")
        Xs = self.scaler_x.transform(np.asarray(X, dtype=np.float32))
        coef_matrix = getattr(self, "_coef_matrix", None)
        intercept_vector = getattr(self, "_intercept_vector", None)
        if coef_matrix is None or intercept_vector is None:
            coef_matrix = np.column_stack(
                [np.asarray(m.coef_, dtype=np.float32).ravel() for m in self.models]
            ).astype(np.float32, copy=False)
            intercept_vector = np.asarray([float(m.intercept_) for m in self.models],
                                          dtype=np.float32)
            self._coef_matrix = coef_matrix
            self._intercept_vector = intercept_vector
        ys = Xs @ coef_matrix + intercept_vector
        return np.asarray(self.scaler_y.inverse_transform(ys), dtype=np.float32)
