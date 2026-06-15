"""rna2grn.model — the GRN activity predictor.

``RegulonEdgeRegressor`` is the shipped model. It predicts each TF->target edge's
activity score from the **expression of that edge's target gene, the TF, and the
TF's regulon mean** via a per-edge ridge regression. This is a supervised,
learned form of regulon-activity inference (cf. VIPER / AUCell / decoupleR): it
**extrapolates** — when a regulon's target genes shift (e.g. a mutation changes a
TF's activity while its mRNA stays flat), the predicted edge activities move
beyond any reference profile, and the perturbed TF's edges rise specifically.

This is the key behavioural difference from a retrieval model: a nearest-neighbor
predictor can only return a blend of existing reference profiles and therefore
cannot reveal a perturbation that is absent from the reference. In an in-silico
regulon-perturbation test the per-edge regressor ranks the perturbed TF at median
rank 2/217 (kNN: 74); it also has higher absolute LOSO accuracy (R²=0.76 vs 0.63).

The model is numpy/scipy-only and tiny (per-edge coefficients + standardization
stats + a sparse regulon-averaging matrix), so the bundle is well under 1 MB.

``KnnGrnRegressor`` is retained as an alternative for pure profile reconstruction,
but it is NOT suitable for perturbation / differential-activity detection.

Input contract for ``predict``: X is (n_pseudobulks, n_feature_genes) in
**CP10k + log1p** space, column-aligned to ``self.feature_genes``.
"""
from __future__ import annotations

from typing import Dict, List, Sequence

import numpy as np
import scipy.sparse as sp


class RegulonEdgeRegressor:
    def __init__(
        self,
        *,
        feature_genes: Sequence[str],
        edge_tf: Sequence[str],
        edge_gene: Sequence[str],
        tf_col: np.ndarray,          # (E,) col index of each edge's TF
        target_col: np.ndarray,      # (E,) col index of each edge's target gene
        tf_of_edge: np.ndarray,      # (E,) index of the edge's TF in tf_list
        regulon_matrix: sp.csr_matrix,  # (n_tf, n_feat) row-averaging over regulon targets
        coefs: np.ndarray,           # (E, 4): intercept, target, TF, regulon
        feat_mu: np.ndarray,         # (E, 3): standardization mean
        feat_sd: np.ndarray,         # (E, 3): standardization std
        tf_list: Sequence[str],
    ) -> None:
        self.feature_genes: List[str] = [str(g) for g in feature_genes]
        self.edge_tf = np.asarray(edge_tf, dtype=object)
        self.edge_gene = np.asarray(edge_gene, dtype=object)
        self.tf_col = np.asarray(tf_col, dtype=np.int64)
        self.target_col = np.asarray(target_col, dtype=np.int64)
        self.tf_of_edge = np.asarray(tf_of_edge, dtype=np.int64)
        self.regulon_matrix = regulon_matrix.tocsr()
        self.coefs = np.asarray(coefs, dtype=np.float64)
        self.feat_mu = np.asarray(feat_mu, dtype=np.float64)
        self.feat_sd = np.asarray(feat_sd, dtype=np.float64)
        self.tf_list = [str(t) for t in tf_list]
        self.n_features_in_ = len(self.feature_genes)
        self.n_edges_ = self.coefs.shape[0]

    def _edge_features(self, X: np.ndarray):
        X = np.asarray(X, dtype=np.float64)
        t = X[:, self.target_col]                       # (n, E) target expression
        f = X[:, self.tf_col]                           # (n, E) TF expression
        reg_all = X @ self.regulon_matrix.T             # (n, n_tf) regulon means
        r = np.asarray(reg_all)[:, self.tf_of_edge]     # (n, E)
        return t, f, r

    def predict(self, X: np.ndarray, *, return_neighbors: bool = False) -> np.ndarray:
        t, f, r = self._edge_features(X)
        sd = self.feat_sd.copy(); sd[sd < 1e-8] = 1.0
        ts = (t - self.feat_mu[:, 0]) / sd[:, 0]
        fs = (f - self.feat_mu[:, 1]) / sd[:, 1]
        rs = (r - self.feat_mu[:, 2]) / sd[:, 2]
        pred = (self.coefs[:, 0]
                + self.coefs[:, 1] * ts
                + self.coefs[:, 2] * fs
                + self.coefs[:, 3] * rs)
        # true connection scores are non-negative; clip small negative extrapolations
        np.clip(pred, 0.0, None, out=pred)
        if return_neighbors:        # not applicable to this model; keep API parity
            return pred, None, None
        return pred

    # convenience: aggregate predicted edges to a per-TF activity vector
    def tf_activity(self, edge_pred: np.ndarray) -> Dict[str, np.ndarray]:
        out = {}
        for ti, tf in enumerate(self.tf_list):
            mask = self.tf_of_edge == ti
            out[tf] = edge_pred[:, mask].mean(axis=1)
        return out


class KnnGrnRegressor:
    """Retrieval baseline (kept for profile reconstruction; NOT for perturbation)."""

    def __init__(self, *, feature_genes, ref_X_std, ref_Y, mean, std, k=15, metric="cosine"):
        self.feature_genes = [str(g) for g in feature_genes]
        self.ref_Y = np.asarray(ref_Y, dtype=np.float32)
        self.mean = np.asarray(mean, dtype=np.float32)
        self.std = np.asarray(std, dtype=np.float32)
        self.k = int(k); self.metric = str(metric)
        ref = np.asarray(ref_X_std, dtype=np.float32)
        if self.metric == "cosine":
            norm = np.linalg.norm(ref, axis=1, keepdims=True); norm[norm < 1e-12] = 1.0
            self._ref = (ref / norm).astype(np.float32)
        else:
            self._ref = ref
        self.n_features_in_ = len(self.feature_genes)

    def _standardize(self, X):
        sd = self.std.copy(); sd[sd < 1e-8] = 1.0
        return (np.asarray(X, dtype=np.float32) - self.mean) / sd

    def predict(self, X, *, return_neighbors=False):
        Xs = self._standardize(X); k = min(self.k, self._ref.shape[0])
        if self.metric == "cosine":
            norm = np.linalg.norm(Xs, axis=1, keepdims=True); norm[norm < 1e-12] = 1.0
            sim = (Xs / norm) @ self._ref.T
            idx = np.argpartition(-sim, kth=k - 1, axis=1)[:, :k]
            row = np.arange(sim.shape[0])[:, None]; s = sim[row, idx]
            w = 1.0 / ((1.0 - s) + 1e-6)
        else:
            d2 = (Xs ** 2).sum(1, keepdims=True) - 2 * (Xs @ self._ref.T) + (self._ref ** 2).sum(1)[None, :]
            idx = np.argpartition(d2, kth=k - 1, axis=1)[:, :k]
            row = np.arange(d2.shape[0])[:, None]
            w = 1.0 / (np.sqrt(np.clip(d2[row, idx], 0, None)) + 1e-6)
        w = w / w.sum(1, keepdims=True)
        pred = np.einsum("ij,ijk->ik", w.astype(np.float32), self.ref_Y[idx])
        if return_neighbors:
            return pred, idx, w
        return pred
