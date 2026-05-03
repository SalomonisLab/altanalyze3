"""Cell-state centroid + per-protein RNA-residual ADT predictor.

A query cell's predicted ADT is::

    pred_adt = state_centroid_adt[cluster_label]
             + residual_model(rna_residual)

where ``state_centroid_adt`` is the mean ADT profile of each reference cell
state (computed once at training time from atlas cells), and
``rna_residual`` is the cell's RNA expression after centering on its assigned
state's mean RNA.

Why this generalises better than global linear models on small donor cohorts:

* The dominant signal driving ADT (cell-state identity) is captured by the
  centroid lookup, which is *shared across donors*.
* The remaining within-state RNA→ADT relationship is usually small in
  magnitude, dominated by biology rather than donor batch, and only kicks
  in once a state assignment is in place.

Inference requires the caller to supply each cell's assigned cluster (the
cellHarmony pipeline already does this upstream). When the cluster label is
unknown, ``predict_with_rna`` falls back to a kNN-in-PCA centroid lookup so
the model still works on raw uploads.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import Ridge
from sklearn.neighbors import NearestNeighbors

from .adt_rna_map import audit_against_genes, load_curated_adt_rna_map


@dataclass
class _StateCentroid:
    cluster: str
    rna_mean: np.ndarray   # (n_rna,)
    adt_mean: np.ndarray   # (n_adt,)
    n_cells: int


class CentroidAdtRegressor:
    """Centroid-based ADT predictor with optional RNA residual head.

    Two modes:
      * ``predict(X, cluster_labels)`` — caller supplies cluster labels (from
        cellHarmony) and we look up centroids directly.
      * ``predict(X)`` — no labels supplied; we use the bundled kNN-in-PCA
        index to assign each query cell to its nearest training cluster.

    Each per-protein RNA residual head is a small Ridge regression on the
    RNA-centered-within-state matrix. Setting ``use_residual=False`` disables
    the residual entirely (centroid-only baseline).
    """

    def __init__(
        self,
        *,
        rna_genes: Sequence[str],
        adt_names: Sequence[str],
        residual_alpha: float = 10.0,
        use_residual: bool = True,
        residual_kind: str = "global_ridge",
        residual_mlp_hidden: Tuple[int, ...] = (128, 64),
        residual_mlp_dropout: float = 0.3,
        residual_mlp_epochs: int = 20,
        residual_mlp_lr: float = 1e-3,
        residual_mlp_batch: int = 1024,
        knn_n_components: int = 50,
        knn_n_neighbors: int = 25,
    ) -> None:
        self.rna_genes = list(rna_genes)
        self.adt_names = list(adt_names)
        self.residual_alpha = float(residual_alpha)
        self.use_residual = bool(use_residual)
        if residual_kind not in {"global_ridge", "targeted", "mlp"}:
            raise ValueError(f"unknown residual_kind: {residual_kind!r}")
        self.residual_kind = residual_kind
        self.residual_mlp_hidden = tuple(int(h) for h in residual_mlp_hidden)
        self.residual_mlp_dropout = float(residual_mlp_dropout)
        self.residual_mlp_epochs = int(residual_mlp_epochs)
        self.residual_mlp_lr = float(residual_mlp_lr)
        self.residual_mlp_batch = int(residual_mlp_batch)
        self.knn_n_components = int(knn_n_components)
        self.knn_n_neighbors = int(knn_n_neighbors)
        self.centroids: Dict[str, _StateCentroid] = {}
        self._cluster_order: List[str] = []
        self._fallback_adt: Optional[np.ndarray] = None
        self._fallback_rna: Optional[np.ndarray] = None
        self._residual_model: Optional[object] = None
        # Targeted-residual heads, indexed in order of self.adt_names
        self._targeted_heads: Optional[List[Tuple[List[int], np.ndarray, float]]] = None
        # MLP residual: torch.nn.Sequential + chosen feature subset
        self._mlp_module = None
        self._mlp_device = None
        self._mlp_feature_indices: Optional[List[int]] = None
        self._mlp_x_mean: Optional[np.ndarray] = None
        self._mlp_x_std: Optional[np.ndarray] = None
        self._mlp_y_mean: Optional[np.ndarray] = None
        self._mlp_y_std: Optional[np.ndarray] = None
        self._svd: Optional[TruncatedSVD] = None
        self._nn: Optional[NearestNeighbors] = None
        self._nn_centroid_clusters: Optional[List[str]] = None
        self.n_features_in_: Optional[int] = None

    def fit(
        self,
        X: np.ndarray,
        Y: np.ndarray,
        cluster_labels: Sequence[str],
    ) -> "CentroidAdtRegressor":
        cluster_labels = np.asarray([str(c) for c in cluster_labels])
        if len(cluster_labels) != X.shape[0]:
            raise ValueError("cluster_labels length must match X rows")
        self.n_features_in_ = int(X.shape[1])
        self.centroids = {}
        order: List[str] = []
        for state in pd.unique(cluster_labels):
            mask = cluster_labels == state
            if mask.sum() < 1:
                continue
            xc = np.asarray(X[mask], dtype=np.float64)
            yc = np.asarray(Y[mask], dtype=np.float64)
            self.centroids[str(state)] = _StateCentroid(
                cluster=str(state),
                rna_mean=xc.mean(axis=0).astype(np.float32),
                adt_mean=yc.mean(axis=0).astype(np.float32),
                n_cells=int(mask.sum()),
            )
            order.append(str(state))
        self._cluster_order = order
        self._fallback_adt = np.asarray(Y, dtype=np.float64).mean(axis=0).astype(np.float32)
        self._fallback_rna = np.asarray(X, dtype=np.float64).mean(axis=0).astype(np.float32)

        if self.use_residual:
            X_resid = np.zeros_like(X, dtype=np.float32)
            Y_resid = np.zeros_like(Y, dtype=np.float32)
            for state in self._cluster_order:
                mask = cluster_labels == state
                X_resid[mask] = X[mask].astype(np.float32) - self.centroids[state].rna_mean
                Y_resid[mask] = Y[mask].astype(np.float32) - self.centroids[state].adt_mean

            if self.residual_kind == "global_ridge":
                self._residual_model = Ridge(alpha=self.residual_alpha, random_state=0)
                self._residual_model.fit(X_resid, Y_resid)
            elif self.residual_kind == "targeted":
                self._fit_targeted_residual(X_resid, Y_resid)
            elif self.residual_kind == "mlp":
                self._fit_mlp_residual(X_resid, Y_resid)

        # Build a kNN index on per-state mean RNA so callers without labels
        # can still get a centroid lookup.
        n_states = len(self._cluster_order)
        if n_states >= 2:
            centroid_matrix = np.stack(
                [self.centroids[s].rna_mean for s in self._cluster_order], axis=0
            )
            n_components = min(self.knn_n_components, centroid_matrix.shape[1] - 1, max(centroid_matrix.shape[0] - 1, 1))
            n_components = max(2, n_components)
            self._svd = TruncatedSVD(n_components=n_components, random_state=0)
            Z = self._svd.fit_transform(centroid_matrix)
            norms = np.linalg.norm(Z, axis=1, keepdims=True)
            norms[norms == 0] = 1.0
            self._nn = NearestNeighbors(n_neighbors=1, algorithm="auto", metric="euclidean")
            self._nn.fit(Z / norms)
            self._nn_centroid_clusters = list(self._cluster_order)

        return self

    def predict(self, X: np.ndarray, cluster_labels: Optional[Sequence[str]] = None) -> np.ndarray:
        X = np.asarray(X, dtype=np.float32)
        if cluster_labels is None:
            cluster_labels = self._assign_clusters_via_knn(X)
        else:
            cluster_labels = np.asarray([str(c) for c in cluster_labels])
        if len(cluster_labels) != X.shape[0]:
            raise ValueError("cluster_labels length must match X rows")

        n = X.shape[0]
        out = np.zeros((n, len(self.adt_names)), dtype=np.float32)
        rna_means = np.zeros_like(X)
        for i, state in enumerate(cluster_labels):
            cent = self.centroids.get(str(state))
            if cent is None:
                out[i] = self._fallback_adt
                rna_means[i] = self._fallback_rna
            else:
                out[i] = cent.adt_mean
                rna_means[i] = cent.rna_mean

        if self.use_residual:
            X_resid = X - rna_means
            residual = self._predict_residual(X_resid)
            if residual is not None:
                out = out + residual.astype(np.float32)
        return out

    # ------------------------------------------------------------------ #
    # Residual variants                                                   #
    # ------------------------------------------------------------------ #

    def _fit_targeted_residual(self, X_resid: np.ndarray, Y_resid: np.ndarray) -> None:
        """Per-ADT Ridge on the canonical RNA partner subset only.

        We fit the residual ``ADT - centroid_ADT`` against
        ``RNA[partners] - centroid_RNA[partners]`` for each ADT independently.
        With ~1-4 features per protein the model cannot memorise donor batch
        effects through irrelevant genes; its coefficients capture
        within-state biology that should travel across donors.
        """
        gene_idx = {g: i for i, g in enumerate(self.rna_genes)}
        adt_to_partners = load_curated_adt_rna_map()
        audit = audit_against_genes(self.adt_names, self.rna_genes, map_obj=adt_to_partners)
        heads: List[Tuple[List[int], np.ndarray, float]] = []
        for j, adt in enumerate(self.adt_names):
            partner_names = audit[adt]["matched"]  # type: ignore[index]
            partner_idx = [gene_idx[g] for g in partner_names if g in gene_idx]
            if not partner_idx:
                heads.append(([], np.zeros(0, dtype=np.float64), 0.0))
                continue
            Xj = X_resid[:, partner_idx].astype(np.float64, copy=False)
            yj = Y_resid[:, j].astype(np.float64, copy=False)
            model = Ridge(alpha=self.residual_alpha, random_state=0)
            model.fit(Xj, yj)
            heads.append((partner_idx, np.asarray(model.coef_, dtype=np.float64).ravel(),
                          float(model.intercept_)))
        self._targeted_heads = heads

    def _fit_mlp_residual(self, X_resid: np.ndarray, Y_resid: np.ndarray) -> None:
        """Small MLP residual on the panel-partner gene subset.

        Subset because feeding all 3000 RNA genes to an MLP just lets it
        re-learn donor signature; restricting to the canonical-partner union
        (~136 genes) caps the donor-overfit attack surface while still giving
        the network enough features to capture non-linear within-state
        structure (saturation, OR-of-genes for multi-subunit proteins).
        """
        try:
            import torch
            from torch import nn
        except ImportError as exc:
            raise RuntimeError("torch required for residual_kind='mlp'") from exc

        # Build the panel-partner union from the curated map.
        adt_to_partners = load_curated_adt_rna_map()
        audit = audit_against_genes(self.adt_names, self.rna_genes, map_obj=adt_to_partners)
        gene_idx = {g: i for i, g in enumerate(self.rna_genes)}
        seen: set[str] = set()
        feat_idx: List[int] = []
        for adt in self.adt_names:
            for g in audit[adt]["matched"]:  # type: ignore[index]
                if g not in seen and g in gene_idx:
                    seen.add(g)
                    feat_idx.append(gene_idx[g])
        if not feat_idx:
            raise ValueError("No canonical partner genes available for MLP residual.")
        self._mlp_feature_indices = feat_idx

        Xs = X_resid[:, feat_idx].astype(np.float32, copy=False)
        Ys = Y_resid.astype(np.float32, copy=False)
        x_mean = Xs.mean(0); x_std = Xs.std(0); x_std[x_std < 1e-6] = 1.0
        y_mean = Ys.mean(0); y_std = Ys.std(0); y_std[y_std < 1e-6] = 1.0
        self._mlp_x_mean = x_mean; self._mlp_x_std = x_std
        self._mlp_y_mean = y_mean; self._mlp_y_std = y_std
        Xn = (Xs - x_mean) / x_std
        Yn = (Ys - y_mean) / y_std

        torch.manual_seed(0)
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self._mlp_device = device
        layers: List[nn.Module] = []
        prev = Xn.shape[1]
        for h in self.residual_mlp_hidden:
            layers += [nn.Linear(prev, h), nn.ReLU(), nn.Dropout(self.residual_mlp_dropout)]
            prev = h
        layers += [nn.Linear(prev, Yn.shape[1])]
        net = nn.Sequential(*layers).to(device)
        opt = torch.optim.Adam(net.parameters(), lr=self.residual_mlp_lr, weight_decay=1e-4)
        loss_fn = nn.MSELoss()

        Xt = torch.from_numpy(Xn)
        Yt = torch.from_numpy(Yn)
        ds = torch.utils.data.TensorDataset(Xt, Yt)
        loader = torch.utils.data.DataLoader(ds, batch_size=self.residual_mlp_batch, shuffle=True)
        for epoch in range(self.residual_mlp_epochs):
            net.train()
            running = 0.0; nb = 0
            for xb, yb in loader:
                xb = xb.to(device); yb = yb.to(device)
                opt.zero_grad()
                pred = net(xb)
                loss = loss_fn(pred, yb)
                loss.backward(); opt.step()
                running += float(loss.item()) * xb.shape[0]; nb += xb.shape[0]
            print(f"[centroid-mlp] epoch {epoch+1}/{self.residual_mlp_epochs} loss={running/max(nb,1):.4f}",
                  flush=True)
        self._mlp_module = net

    def _predict_residual(self, X_resid: np.ndarray) -> Optional[np.ndarray]:
        if self.residual_kind == "global_ridge":
            if self._residual_model is None:
                return None
            return self._residual_model.predict(X_resid)
        if self.residual_kind == "targeted":
            if self._targeted_heads is None:
                return None
            n = X_resid.shape[0]
            out = np.zeros((n, len(self.adt_names)), dtype=np.float32)
            for j, (idx, coef, intercept) in enumerate(self._targeted_heads):
                if not idx:
                    continue
                Xj = X_resid[:, idx].astype(np.float64, copy=False)
                out[:, j] = (Xj @ coef + intercept).astype(np.float32)
            return out
        if self.residual_kind == "mlp":
            import torch
            if self._mlp_module is None or self._mlp_feature_indices is None:
                return None
            Xs = X_resid[:, self._mlp_feature_indices].astype(np.float32, copy=False)
            Xn = (Xs - self._mlp_x_mean) / self._mlp_x_std
            self._mlp_module.eval()
            with torch.no_grad():
                yt = self._mlp_module(torch.from_numpy(Xn).to(self._mlp_device))
            preds = yt.cpu().numpy()
            preds = preds * self._mlp_y_std + self._mlp_y_mean
            return preds
        return None

    def _assign_clusters_via_knn(self, X: np.ndarray) -> np.ndarray:
        if self._svd is None or self._nn is None or self._nn_centroid_clusters is None:
            if not self._cluster_order:
                return np.array(["__fallback__"] * X.shape[0])
            # Single cluster — return it.
            return np.array([self._cluster_order[0]] * X.shape[0])
        Z = self._svd.transform(np.asarray(X, dtype=np.float32))
        norms = np.linalg.norm(Z, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        Z = Z / norms
        _, idx = self._nn.kneighbors(Z, n_neighbors=1, return_distance=True)
        return np.asarray([self._nn_centroid_clusters[int(i)] for i in idx.ravel()])
