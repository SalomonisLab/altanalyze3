"""Train and benchmark RNA -> ADT imputation models on a CITE-seq atlas.

Usage::

    python -m altanalyze3.components.rna2adt.training \
        --train-h5ad /path/to/atlas.h5ad \
        --heldout-h5ad /path/to/atlas.heldout.h5ad \
        --output-dir components/rna2adt/artifacts \
        --bundle-path components/rna2adt/rna2adt_bm_bundle.pkl

Trains three candidate models on the bone marrow atlas, evaluates them on
both a cell-level holdout and a leave-one-donor-out split, prints a per-model
report (mean per-protein Pearson r, Spearman rho, RMSE), and writes the
winning model out as a bundle compatible with ``rna2adt.api.Rna2AdtBundle``.
"""

from __future__ import annotations

import argparse
import json
import pickle
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import MultiTaskElasticNet, Ridge
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

from .adt_rna_map import audit_against_genes, load_curated_adt_rna_map
from .centroid import CentroidAdtRegressor
from .generalizable import (
    GeneralizableElasticNet,
    GeneralizableElasticRank,
    GeneralizableElasticSpline,
    GeneralizableRidge,
    GeneralizableSplineGAM,
    GeneralizableXGBoost,
)
from .rna2lipid_arch import Rna2LipidArchFull, Rna2LipidArchPanel, Rna2LipidArchPanelPerProtein
from .targeted import TargetedAdtRegressor


_ADT_PREFIX = "Hu."


def _is_adt(name: str) -> bool:
    return str(name).startswith(_ADT_PREFIX)


def _to_dense(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        return np.asarray(matrix.todense(), dtype=np.float32)
    return np.asarray(matrix, dtype=np.float32)


@dataclass
class SplitData:
    X_train: np.ndarray
    Y_train: np.ndarray
    X_test: np.ndarray
    Y_test: np.ndarray
    train_meta: pd.DataFrame
    test_meta: pd.DataFrame
    rna_genes: List[str]
    adt_names: List[str]
    label: str = "split"
    cluster_col: Optional[str] = None
    train_clusters: Optional[np.ndarray] = None
    test_clusters: Optional[np.ndarray] = None


def load_split_arrays(
    train_path: Path,
    *,
    test_path: Optional[Path] = None,
    leave_out_donor: Optional[str] = None,
    donor_col: str = "Donor",
    cluster_col: str = "Level 3 Multimodal",
    max_train_cells: Optional[int] = None,
    max_test_cells: Optional[int] = None,
    rna_top_n: int = 4000,
    extra_rna_genes: Optional[Sequence[str]] = None,
    seed: int = 0,
    require_cluster_col: bool = True,
) -> SplitData:
    """Load RNA / ADT matrices and produce a train/test split.

    Either ``test_path`` (a separate h5ad) or ``leave_out_donor`` must be
    provided. RNA features are reduced to the union of the top-variance genes
    (across the training cells) and any ``extra_rna_genes`` the caller passes
    in (e.g. surface marker partners of the ADT panel).
    """
    rng = np.random.default_rng(seed)

    print(f"[load] reading train atlas: {train_path}")
    full = ad.read_h5ad(train_path)
    var_names = np.array([str(v) for v in full.var_names])
    adt_mask = np.array([_is_adt(v) for v in var_names])
    rna_mask = ~adt_mask
    rna_names = var_names[rna_mask].tolist()
    adt_names = var_names[adt_mask].tolist()
    print(f"[load] {full.n_obs} cells, {rna_mask.sum()} RNA vars, {adt_mask.sum()} ADT vars")

    if leave_out_donor is not None:
        if donor_col not in full.obs.columns:
            raise KeyError(f"Donor column {donor_col!r} not in obs")
        train_idx = np.where(full.obs[donor_col].astype(str).to_numpy() != str(leave_out_donor))[0]
        test_idx = np.where(full.obs[donor_col].astype(str).to_numpy() == str(leave_out_donor))[0]
        if len(test_idx) == 0:
            raise ValueError(f"Donor {leave_out_donor} not present in {donor_col}")
        train_obs = full.obs.iloc[train_idx]
        test_obs = full.obs.iloc[test_idx]
        train_X_full = full.X[train_idx]
        test_X_full = full.X[test_idx]
        label = f"LODO[{leave_out_donor}]"
    elif test_path is not None:
        print(f"[load] reading heldout atlas: {test_path}")
        held = ad.read_h5ad(test_path)
        held_var = np.array([str(v) for v in held.var_names])
        if not np.array_equal(held_var, var_names):
            common, train_pos, test_pos = np.intersect1d(var_names, held_var, return_indices=True)
            order = np.argsort(np.concatenate([np.where(adt_mask & np.isin(var_names, common))[0],
                                              np.where(rna_mask & np.isin(var_names, common))[0]]))
            mask_keep_train = np.isin(var_names, common)
            train_X_full = full.X[:, mask_keep_train]
            mask_keep_test = np.isin(held_var, common)
            test_X_full = held.X[:, mask_keep_test]
            var_names = var_names[mask_keep_train]
            adt_mask = np.array([_is_adt(v) for v in var_names])
            rna_mask = ~adt_mask
            rna_names = var_names[rna_mask].tolist()
            adt_names = var_names[adt_mask].tolist()
        else:
            train_X_full = full.X
            test_X_full = held.X
        train_obs = full.obs
        test_obs = held.obs
        label = f"cell_holdout[{Path(test_path).name}]"
    else:
        raise ValueError("Must provide either test_path or leave_out_donor")

    if max_train_cells is not None and train_obs.shape[0] > max_train_cells:
        sel = rng.choice(train_obs.shape[0], size=max_train_cells, replace=False)
        sel.sort()
        train_obs = train_obs.iloc[sel]
        train_X_full = train_X_full[sel]
    if max_test_cells is not None and test_obs.shape[0] > max_test_cells:
        sel = rng.choice(test_obs.shape[0], size=max_test_cells, replace=False)
        sel.sort()
        test_obs = test_obs.iloc[sel]
        test_X_full = test_X_full[sel]

    print(f"[load] split={label} train={train_obs.shape[0]} test={test_obs.shape[0]}")

    print("[load] densifying ADT block")
    Y_train = _to_dense(train_X_full[:, adt_mask])
    Y_test = _to_dense(test_X_full[:, adt_mask])

    print(f"[load] selecting RNA features (top {rna_top_n} variance + companions)")
    rna_idx_full = np.where(rna_mask)[0]
    rna_train = train_X_full[:, rna_mask]
    var = _approx_variance(rna_train)
    var = np.asarray(var).ravel()
    top_idx = np.argsort(-var)[:rna_top_n]
    keep_genes = {rna_names[i] for i in top_idx.tolist()}
    if extra_rna_genes:
        keep_genes.update(g for g in extra_rna_genes if g)
    keep_set = keep_genes
    keep_positions = np.array([i for i, g in enumerate(rna_names) if g in keep_set], dtype=int)
    keep_genes_list = [rna_names[i] for i in keep_positions.tolist()]
    print(f"[load] kept {len(keep_genes_list)} RNA genes")

    X_train = _to_dense(rna_train[:, keep_positions])
    X_test = _to_dense(test_X_full[:, rna_mask][:, keep_positions])

    train_clusters = None
    test_clusters = None
    if require_cluster_col:
        if cluster_col in train_obs.columns:
            train_clusters = train_obs[cluster_col].astype(str).to_numpy()
        if cluster_col in test_obs.columns:
            test_clusters = test_obs[cluster_col].astype(str).to_numpy()
        if train_clusters is None or test_clusters is None:
            print(f"[load] WARNING: cluster col {cluster_col!r} missing from one split; "
                  "centroid model will fall back to kNN assignment.")

    return SplitData(
        X_train=X_train,
        Y_train=Y_train,
        X_test=X_test,
        Y_test=Y_test,
        train_meta=train_obs.reset_index().rename(columns={"index": "cell_id"}),
        test_meta=test_obs.reset_index().rename(columns={"index": "cell_id"}),
        rna_genes=keep_genes_list,
        adt_names=adt_names,
        label=label,
        cluster_col=cluster_col,
        train_clusters=train_clusters,
        test_clusters=test_clusters,
    )


def _approx_variance(matrix, *, sample_rows: int = 20000, seed: int = 0) -> np.ndarray:
    n = matrix.shape[0]
    if n > sample_rows:
        rng = np.random.default_rng(seed)
        idx = rng.choice(n, size=sample_rows, replace=False)
        idx.sort()
        sub = matrix[idx]
    else:
        sub = matrix
    sub = _to_dense(sub)
    return sub.var(axis=0)


# ---------------------------------------------------------------------------
# kNN-in-PCA regressor
# ---------------------------------------------------------------------------


class KnnAdtRegressor:
    """Predicts ADTs by averaging the ADT vectors of the k nearest reference cells.

    Stored state: PCA components (TruncatedSVD), reference PCA embedding, the
    reference ADT matrix, and the fit ``NearestNeighbors`` index. Inference
    transforms an input matrix into PCA space, queries the index, and gathers
    a (weighted) average of the reference ADTs.
    """

    def __init__(self, *, n_components: int = 50, n_neighbors: int = 25, weights: str = "distance") -> None:
        self.n_components = int(n_components)
        self.n_neighbors = int(n_neighbors)
        self.weights = weights
        self._svd: Optional[TruncatedSVD] = None
        self._nn: Optional[NearestNeighbors] = None
        self._Y_ref: Optional[np.ndarray] = None
        self.n_features_in_: Optional[int] = None

    def fit(self, X: np.ndarray, Y: np.ndarray) -> "KnnAdtRegressor":
        self.n_features_in_ = int(X.shape[1])
        self._svd = TruncatedSVD(n_components=min(self.n_components, X.shape[1] - 1, X.shape[0] - 1), random_state=0)
        Z = self._svd.fit_transform(X)
        Z = Z.astype(np.float32, copy=False)
        # L2-normalize for cosine-like distances
        norms = np.linalg.norm(Z, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        Z = Z / norms
        self._nn = NearestNeighbors(n_neighbors=self.n_neighbors, algorithm="auto", metric="euclidean")
        self._nn.fit(Z)
        self._Y_ref = np.asarray(Y, dtype=np.float32)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        if self._svd is None or self._nn is None or self._Y_ref is None:
            raise RuntimeError("KnnAdtRegressor must be fit before predict()")
        Z = self._svd.transform(np.asarray(X, dtype=np.float32))
        Z = Z.astype(np.float32, copy=False)
        norms = np.linalg.norm(Z, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        Z = Z / norms
        dist, idx = self._nn.kneighbors(Z, return_distance=True)
        if self.weights == "distance":
            inv = 1.0 / np.maximum(dist, 1e-6)
            w = inv / inv.sum(axis=1, keepdims=True)
        else:
            w = np.full_like(dist, 1.0 / self.n_neighbors, dtype=np.float32)
        # Gather: for each query row, weighted average of k reference Y rows
        Y_neighbors = self._Y_ref[idx]  # (n_query, k, n_adts)
        out = (w[..., None] * Y_neighbors).sum(axis=1)
        return out.astype(np.float32)


# ---------------------------------------------------------------------------
# MLP regressor (lightweight, optional torch dependency)
# ---------------------------------------------------------------------------


class MlpAdtRegressor:
    def __init__(self, *, hidden: Tuple[int, ...] = (512, 256), dropout: float = 0.1,
                 lr: float = 1e-3, batch_size: int = 1024, epochs: int = 25, seed: int = 0) -> None:
        self.hidden = tuple(int(h) for h in hidden)
        self.dropout = float(dropout)
        self.lr = float(lr)
        self.batch_size = int(batch_size)
        self.epochs = int(epochs)
        self.seed = int(seed)
        self._net = None
        self._device = None
        self.n_features_in_: Optional[int] = None

    def fit(self, X: np.ndarray, Y: np.ndarray) -> "MlpAdtRegressor":
        try:
            import torch
            from torch import nn
        except ImportError as exc:
            raise RuntimeError("torch is required to fit MlpAdtRegressor") from exc
        torch.manual_seed(self.seed)
        self._device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.n_features_in_ = int(X.shape[1])

        layers: List[nn.Module] = []
        prev = X.shape[1]
        for h in self.hidden:
            layers += [nn.Linear(prev, h), nn.ReLU(), nn.Dropout(self.dropout)]
            prev = h
        layers += [nn.Linear(prev, Y.shape[1])]
        self._net = nn.Sequential(*layers).to(self._device)
        opt = torch.optim.Adam(self._net.parameters(), lr=self.lr)
        loss_fn = nn.MSELoss()

        X_t = torch.from_numpy(np.asarray(X, dtype=np.float32))
        Y_t = torch.from_numpy(np.asarray(Y, dtype=np.float32))
        ds = torch.utils.data.TensorDataset(X_t, Y_t)
        loader = torch.utils.data.DataLoader(ds, batch_size=self.batch_size, shuffle=True, drop_last=False)

        for epoch in range(self.epochs):
            t0 = time.time()
            self._net.train()
            running = 0.0
            n = 0
            for xb, yb in loader:
                xb = xb.to(self._device)
                yb = yb.to(self._device)
                opt.zero_grad()
                pred = self._net(xb)
                loss = loss_fn(pred, yb)
                loss.backward()
                opt.step()
                running += float(loss.item()) * xb.shape[0]
                n += xb.shape[0]
            print(f"[mlp] epoch {epoch+1}/{self.epochs} loss={running/max(n,1):.4f} ({time.time()-t0:.1f}s)")
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        import torch
        if self._net is None:
            raise RuntimeError("MlpAdtRegressor must be fit before predict()")
        self._net.eval()
        X_t = torch.from_numpy(np.asarray(X, dtype=np.float32)).to(self._device)
        with torch.no_grad():
            preds = self._net(X_t).cpu().numpy()
        return preds.astype(np.float32)


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------


@dataclass
class EvalReport:
    model_name: str
    split_label: str
    mean_pearson: float
    median_pearson: float
    mean_spearman: float
    rmse: float
    per_protein: pd.DataFrame
    fit_seconds: float = 0.0
    predict_seconds: float = 0.0


def _evaluate(model_name: str, split_label: str, Y_true: np.ndarray, Y_pred: np.ndarray,
              adt_names: Sequence[str], *, fit_seconds: float = 0.0, predict_seconds: float = 0.0) -> EvalReport:
    rows = []
    for j, name in enumerate(adt_names):
        yt = Y_true[:, j]
        yp = Y_pred[:, j]
        if np.std(yt) < 1e-9 or np.std(yp) < 1e-9:
            pearson = float("nan")
            spear = float("nan")
        else:
            pearson = float(np.corrcoef(yt, yp)[0, 1])
            try:
                spear = float(spearmanr(yt, yp).correlation)
            except Exception:
                spear = float("nan")
        rmse_j = float(np.sqrt(np.mean((yt - yp) ** 2)))
        rows.append({"adt": name, "pearson": pearson, "spearman": spear, "rmse": rmse_j,
                     "true_mean": float(yt.mean()), "pred_mean": float(yp.mean())})
    per = pd.DataFrame(rows)
    return EvalReport(
        model_name=model_name,
        split_label=split_label,
        mean_pearson=float(np.nanmean(per["pearson"])),
        median_pearson=float(np.nanmedian(per["pearson"])),
        mean_spearman=float(np.nanmean(per["spearman"])),
        rmse=float(np.sqrt(np.mean((Y_true - Y_pred) ** 2))),
        per_protein=per.sort_values("pearson", ascending=False).reset_index(drop=True),
        fit_seconds=fit_seconds,
        predict_seconds=predict_seconds,
    )


def _summary_row(r: EvalReport) -> Dict[str, object]:
    return {
        "model": r.model_name,
        "split": r.split_label,
        "mean_pearson": round(r.mean_pearson, 4),
        "median_pearson": round(r.median_pearson, 4),
        "mean_spearman": round(r.mean_spearman, 4),
        "rmse": round(r.rmse, 4),
        "fit_s": round(r.fit_seconds, 1),
        "predict_s": round(r.predict_seconds, 1),
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


@dataclass
class TrainedCandidate:
    name: str
    model: object
    scaler_x: Optional[object]
    scaler_y: Optional[object]
    target_scaling_mode: str = "none"
    extra_metadata: Dict[str, object] = field(default_factory=dict)


def _fit_knn(split: SplitData, *, n_components: int, n_neighbors: int) -> TrainedCandidate:
    print(f"[fit] kNN n_components={n_components} k={n_neighbors}")
    scaler_x = StandardScaler(with_mean=True, with_std=True)
    Xtr = scaler_x.fit_transform(split.X_train)
    model = KnnAdtRegressor(n_components=n_components, n_neighbors=n_neighbors, weights="distance")
    model.fit(Xtr, split.Y_train)
    return TrainedCandidate(
        name=f"knn_pca{n_components}_k{n_neighbors}",
        model=model,
        scaler_x=scaler_x,
        scaler_y=None,
        target_scaling_mode="none",
        extra_metadata={"approach": "knn_in_pca", "n_components": n_components, "n_neighbors": n_neighbors},
    )


def _fit_elastic(split: SplitData, *, alpha: float, l1_ratio: float) -> TrainedCandidate:
    print(f"[fit] ElasticNet alpha={alpha} l1_ratio={l1_ratio}", flush=True)
    scaler_x = StandardScaler(with_mean=True, with_std=True)
    scaler_y = StandardScaler(with_mean=True, with_std=True)
    Xtr = scaler_x.fit_transform(split.X_train)
    Ytr = scaler_y.fit_transform(split.Y_train)
    model = MultiTaskElasticNet(alpha=alpha, l1_ratio=l1_ratio, max_iter=500, tol=1e-3, random_state=0)
    model.fit(Xtr, Ytr)
    return TrainedCandidate(
        name=f"elastic_a{alpha}_l1{l1_ratio}",
        model=model,
        scaler_x=scaler_x,
        scaler_y=scaler_y,
        target_scaling_mode="standard",
        extra_metadata={"approach": "multi_task_elastic_net", "alpha": alpha, "l1_ratio": l1_ratio},
    )


def _fit_rna2lipid_arch_per_protein(
    split: SplitData,
    *,
    alpha: float = 0.01,
    l1_ratio: float = 0.5,
    max_iter: int = 2000,
    feature_source: str = "curated",
    whitelist_path: Optional[Path] = None,
    label_suffix: str = "",
) -> TrainedCandidate:
    print(f"[fit] rna2lipid-arch-per-protein alpha={alpha} l1_ratio={l1_ratio} "
          f"feature_source={feature_source}", flush=True)
    model = Rna2LipidArchPanelPerProtein(
        rna_genes=split.rna_genes,
        adt_names=split.adt_names,
        alpha=alpha,
        l1_ratio=l1_ratio,
        max_iter=max_iter,
        feature_source=feature_source,
        whitelist_path=whitelist_path,
    )
    name = f"rna2lipid_arch_panel_per_protein_a{alpha}_l1{l1_ratio}_{feature_source}"
    if label_suffix:
        name = name + "_" + label_suffix
    model.fit(split.X_train, split.Y_train)
    return TrainedCandidate(
        name=name,
        model=model,
        scaler_x=None,
        scaler_y=None,
        target_scaling_mode="none",
        extra_metadata={
            "approach": "rna2lipid_arch_panel_per_protein",
            "feature_subset": "panel",
            "feature_source": feature_source,
            "alpha": alpha,
            "l1_ratio": l1_ratio,
            "_lipidarch_model": True,
        },
    )


def _fit_rna2lipid_arch(
    split: SplitData,
    *,
    feature_subset: str,
    cv: int = 5,
    max_iter: int = 1500,
    n_jobs: int = 4,
    l1_ratio: Sequence[float] | float = (0.1, 0.3, 0.5, 0.7, 0.9),
    alphas: Optional[Sequence[float]] = None,
    use_cv: bool = True,
    fixed_alpha: float = 0.01,
    label_suffix: str = "",
) -> TrainedCandidate:
    print(f"[fit] rna2lipid-arch feature_subset={feature_subset} cv={cv} use_cv={use_cv}", flush=True)
    common = dict(rna_genes=split.rna_genes, adt_names=split.adt_names,
                  l1_ratio=l1_ratio, alphas=alphas, cv=cv, max_iter=max_iter, n_jobs=n_jobs,
                  use_cv=use_cv, fixed_alpha=fixed_alpha)
    if feature_subset == "panel":
        model = Rna2LipidArchPanel(**common)
        base = "rna2lipid_arch_panel"
    else:
        model = Rna2LipidArchFull(**common)
        base = "rna2lipid_arch_full"
    if not use_cv:
        suffix = f"_a{fixed_alpha}_l1{float(l1_ratio[0]) if isinstance(l1_ratio, (list, tuple, np.ndarray)) else float(l1_ratio)}"
    else:
        suffix = f"_cv{cv}"
    if label_suffix:
        suffix = suffix + "_" + label_suffix
    name = base + suffix
    model.fit(split.X_train, split.Y_train)
    return TrainedCandidate(
        name=name,
        model=model,
        scaler_x=None,
        scaler_y=None,
        target_scaling_mode="none",
        extra_metadata={
            "approach": f"rna2lipid_arch_{feature_subset}",
            "feature_subset": feature_subset,
            "cv": cv,
            "use_cv": use_cv,
            "fixed_alpha": fixed_alpha,
            "_lipidarch_model": True,
        },
    )


def _fit_generalizable(
    split: SplitData,
    *,
    head_kind: str,
    alpha: float = 1.0,
    l1_ratio: float = 0.5,
    n_pca: int = 50,
    n_neighbors: int = 25,
    xgb_n_estimators: int = 200,
    xgb_max_depth: int = 4,
    xgb_lr: float = 0.08,
) -> TrainedCandidate:
    print(f"[fit] Generalizable head_kind={head_kind} alpha={alpha} n_pca={n_pca} k={n_neighbors}",
          flush=True)
    common = dict(rna_genes=split.rna_genes, adt_names=split.adt_names,
                  n_pca=n_pca, n_neighbors=n_neighbors)
    if head_kind == "ridge":
        model = GeneralizableRidge(alpha=alpha, **common)
        name = f"genz_ridge_a{alpha}_k{n_neighbors}"
    elif head_kind == "elastic":
        model = GeneralizableElasticNet(alpha=alpha, l1_ratio=l1_ratio, **common)
        name = f"genz_elastic_a{alpha}_l1{l1_ratio}_k{n_neighbors}"
    elif head_kind == "elastic_rank":
        model = GeneralizableElasticRank(alpha=alpha, l1_ratio=l1_ratio, **common)
        name = f"genz_elastic_rank_a{alpha}_l1{l1_ratio}_k{n_neighbors}"
    elif head_kind == "elastic_spline":
        model = GeneralizableElasticSpline(alpha=alpha, l1_ratio=l1_ratio, **common)
        name = f"genz_elastic_spline_a{alpha}_l1{l1_ratio}_k{n_neighbors}"
    elif head_kind == "spline":
        model = GeneralizableSplineGAM(alpha=alpha, **common)
        name = f"genz_spline_a{alpha}_k{n_neighbors}"
    elif head_kind == "xgboost":
        model = GeneralizableXGBoost(
            n_estimators=xgb_n_estimators,
            max_depth=xgb_max_depth,
            learning_rate=xgb_lr,
            **common,
        )
        name = f"genz_xgb_n{xgb_n_estimators}_d{xgb_max_depth}_k{n_neighbors}"
    else:
        raise ValueError(f"unknown head_kind: {head_kind}")
    model.fit(split.X_train, split.Y_train)
    return TrainedCandidate(
        name=name,
        model=model,
        scaler_x=None,
        scaler_y=None,
        target_scaling_mode="none",
        extra_metadata={
            "approach": f"generalizable_{head_kind}",
            "head_kind": head_kind,
            "alpha": alpha,
            "n_pca": n_pca,
            "n_neighbors": n_neighbors,
            "_generalizable_model": True,
        },
    )


def _fit_centroid(
    split: SplitData,
    *,
    use_residual: bool,
    residual_alpha: float,
    residual_kind: str = "global_ridge",
    mlp_hidden: Tuple[int, ...] = (128, 64),
    mlp_dropout: float = 0.3,
    mlp_epochs: int = 20,
) -> TrainedCandidate:
    if split.train_clusters is None:
        raise ValueError("Centroid model needs cluster labels in training split.")
    print(f"[fit] Centroid use_residual={use_residual} residual_kind={residual_kind} "
          f"residual_alpha={residual_alpha}", flush=True)
    model = CentroidAdtRegressor(
        rna_genes=split.rna_genes,
        adt_names=split.adt_names,
        residual_alpha=residual_alpha,
        use_residual=use_residual,
        residual_kind=residual_kind,
        residual_mlp_hidden=mlp_hidden,
        residual_mlp_dropout=mlp_dropout,
        residual_mlp_epochs=mlp_epochs,
    )
    model.fit(split.X_train, split.Y_train, split.train_clusters)
    if not use_residual:
        suffix = "only"
    elif residual_kind == "global_ridge":
        suffix = f"resid_a{residual_alpha}"
    elif residual_kind == "targeted":
        suffix = f"resid_targeted_a{residual_alpha}"
    elif residual_kind == "mlp":
        suffix = f"resid_mlp_{'x'.join(str(h) for h in mlp_hidden)}_d{mlp_dropout}"
    else:
        suffix = residual_kind
    return TrainedCandidate(
        name=f"centroid_{suffix}",
        model=model,
        scaler_x=None,
        scaler_y=None,
        target_scaling_mode="none",
        extra_metadata={
            "approach": (
                "centroid_only" if not use_residual
                else f"centroid_plus_{residual_kind}_residual"
            ),
            "residual_alpha": residual_alpha,
            "use_residual": use_residual,
            "residual_kind": residual_kind,
            "_centroid_model": True,
        },
    )


def _fit_targeted(
    split: SplitData,
    *,
    method: str,
    alpha: float,
    l1_ratio: float = 0.5,
    share_panel_context: bool = False,
) -> TrainedCandidate:
    """Per-ADT regression on canonical RNA partners only."""
    print(f"[fit] Targeted method={method} alpha={alpha} share_panel_context={share_panel_context}", flush=True)
    map_obj = load_curated_adt_rna_map()
    audit = audit_against_genes(split.adt_names, split.rna_genes, map_obj=map_obj)
    # Union of all partner genes that exist in this RNA matrix.
    partner_union: List[str] = []
    seen: set[str] = set()
    for adt in split.adt_names:
        for g in audit[adt]["matched"]:
            if g not in seen:
                seen.add(g)
                partner_union.append(g)
    rna_idx = [split.rna_genes.index(g) for g in partner_union]
    X_train_targeted = split.X_train[:, rna_idx]
    print(f"[fit] targeted feature union: {len(partner_union)} unique RNA partners", flush=True)
    n_unmapped = sum(1 for adt, info in audit.items() if info.get("unmapped"))
    n_isotype = sum(1 for adt, info in audit.items() if info.get("isotype"))
    n_no_match = sum(1 for adt, info in audit.items()
                     if not info["matched"] and not info["unmapped"] and not info["isotype"])
    print(f"[fit] targeted audit: unmapped={n_unmapped} no_match={n_no_match} isotype={n_isotype}", flush=True)

    model = TargetedAdtRegressor(
        feature_names=partner_union,
        adt_names=split.adt_names,
        adt_to_features={adt: list(audit[adt]["matched"]) for adt in split.adt_names},
        method=method,
        alpha=alpha,
        l1_ratio=l1_ratio,
        share_panel_context=share_panel_context,
    )
    model.fit(X_train_targeted, split.Y_train)
    # Stash the partner-union slice so _predict can re-slice the test matrix.
    extra: Dict[str, object] = {
        "approach": f"targeted_{method}",
        "alpha": alpha,
        "n_partner_union": len(partner_union),
        "n_unmapped": int(n_unmapped),
        "n_isotype": int(n_isotype),
        "n_no_match": int(n_no_match),
        "_targeted_partner_union": list(partner_union),
        "_targeted_partner_indices": list(rna_idx),
    }
    if method == "elastic":
        extra["l1_ratio"] = l1_ratio
    return TrainedCandidate(
        name=f"targeted_{method}_a{alpha}",
        model=model,
        scaler_x=None,
        scaler_y=None,
        target_scaling_mode="none",
        extra_metadata=extra,
    )


def _fit_ridge(split: SplitData, *, alpha: float) -> TrainedCandidate:
    """Closed-form multi-output Ridge regression on standardized inputs.

    Much faster than MultiTaskElasticNet for our shape (30k cells, 3k genes,
    129 proteins) — a single SVD vs. iterative coordinate descent — and tends
    to match it on dense ADT prediction since we already use a strong RNA
    feature filter.
    """
    print(f"[fit] Ridge alpha={alpha}", flush=True)
    scaler_x = StandardScaler(with_mean=True, with_std=True)
    scaler_y = StandardScaler(with_mean=True, with_std=True)
    Xtr = scaler_x.fit_transform(split.X_train)
    Ytr = scaler_y.fit_transform(split.Y_train)
    model = Ridge(alpha=alpha, solver="auto", random_state=0)
    model.fit(Xtr, Ytr)
    return TrainedCandidate(
        name=f"ridge_a{alpha}",
        model=model,
        scaler_x=scaler_x,
        scaler_y=scaler_y,
        target_scaling_mode="standard",
        extra_metadata={"approach": "ridge", "alpha": alpha},
    )


def _fit_mlp(split: SplitData, *, hidden: Tuple[int, ...], epochs: int) -> TrainedCandidate:
    print(f"[fit] MLP hidden={hidden} epochs={epochs}")
    scaler_x = StandardScaler(with_mean=True, with_std=True)
    Xtr = scaler_x.fit_transform(split.X_train)
    model = MlpAdtRegressor(hidden=hidden, epochs=epochs)
    model.fit(Xtr, split.Y_train)
    return TrainedCandidate(
        name=f"mlp_{'x'.join(str(h) for h in hidden)}_e{epochs}",
        model=model,
        scaler_x=scaler_x,
        scaler_y=None,
        target_scaling_mode="none",
        extra_metadata={"approach": "mlp", "hidden": list(hidden), "epochs": epochs},
    )


def _predict(candidate: TrainedCandidate, X_test: np.ndarray, Y_test: np.ndarray,
             *, test_clusters: Optional[np.ndarray] = None) -> Tuple[np.ndarray, float]:
    t0 = time.time()
    targeted_idx = candidate.extra_metadata.get("_targeted_partner_indices")
    is_centroid = bool(candidate.extra_metadata.get("_centroid_model"))
    if targeted_idx is not None:
        X = X_test[:, targeted_idx]
    elif candidate.scaler_x is not None:
        X = candidate.scaler_x.transform(X_test)
    else:
        X = X_test
    if is_centroid:
        pred = candidate.model.predict(X, cluster_labels=test_clusters)
    else:
        pred = candidate.model.predict(X)
    if candidate.target_scaling_mode == "standard" and candidate.scaler_y is not None:
        pred = candidate.scaler_y.inverse_transform(pred)
    return np.asarray(pred, dtype=np.float32), time.time() - t0


def save_bundle(
    candidate: TrainedCandidate,
    *,
    bundle_path: Path,
    rna_genes: Sequence[str],
    adt_names: Sequence[str],
    metadata: Dict[str, object],
) -> Path:
    # For targeted models, swap the bundle's X_columns to the partner-only
    # union so Rna2AdtBundle.predict feeds the model exactly its expected
    # columns (and we can drop the StandardScaler entirely).
    targeted_partner_union = candidate.extra_metadata.get("_targeted_partner_union")
    if targeted_partner_union:
        x_columns = list(targeted_partner_union)
    else:
        x_columns = list(rna_genes)
    extra_meta = {k: v for k, v in candidate.extra_metadata.items() if not k.startswith("_")}
    bundle = {
        "model": candidate.model,
        "scaler_x": candidate.scaler_x,
        "scaler_y": candidate.scaler_y,
        "X_columns": x_columns,
        "Y_columns": list(adt_names),
        "metadata": {
            **metadata,
            "approach": candidate.name,
            "target_scaling": {"mode": candidate.target_scaling_mode},
            **extra_meta,
        },
    }
    bundle_path.parent.mkdir(parents=True, exist_ok=True)
    with bundle_path.open("wb") as fh:
        pickle.dump(bundle, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"[save] bundle -> {bundle_path}")
    return bundle_path


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--train-h5ad", required=True, type=Path)
    p.add_argument("--heldout-h5ad", type=Path, default=None)
    p.add_argument("--output-dir", type=Path, default=Path("artifacts"))
    p.add_argument("--bundle-path", type=Path, default=Path("rna2adt_bm_bundle.pkl"))
    p.add_argument("--targeted-bundle-path", type=Path, default=None,
                   help="If set, always write the targeted candidate to this path (in addition to "
                        "the winning bundle at --bundle-path).")
    p.add_argument("--centroid-bundle-path", type=Path, default=None,
                   help="If set, always write the FIRST centroid candidate to this path "
                        "(usually the global_ridge variant, kept for evaluation).")
    p.add_argument("--centroid-targeted-bundle-path", type=Path, default=None,
                   help="If set, write the centroid+targeted-residual candidate.")
    p.add_argument("--centroid-mlp-bundle-path", type=Path, default=None,
                   help="If set, write the centroid+MLP-residual candidate.")
    p.add_argument("--centroid-only-bundle-path", type=Path, default=None,
                   help="If set, write the centroid-only (no residual) candidate.")
    p.add_argument("--genz-kinds", default="",
                   help="Comma-separated generalizable heads to fit. "
                        "Choices: ridge, elastic, spline, xgboost. Empty disables.")
    p.add_argument("--genz-alpha", type=float, default=1.0)
    p.add_argument("--genz-elastic-l1", type=float, default=0.5)
    p.add_argument("--genz-n-pca", type=int, default=50)
    p.add_argument("--genz-n-neighbors", type=int, default=25)
    p.add_argument("--genz-xgb-n-estimators", type=int, default=200)
    p.add_argument("--genz-xgb-max-depth", type=int, default=4)
    p.add_argument("--genz-xgb-lr", type=float, default=0.08)
    p.add_argument("--genz-ridge-bundle-path", type=Path, default=None)
    p.add_argument("--genz-elastic-bundle-path", type=Path, default=None)
    p.add_argument("--genz-spline-bundle-path", type=Path, default=None)
    p.add_argument("--genz-xgb-bundle-path", type=Path, default=None)
    p.add_argument("--max-train-cells", type=int, default=40000)
    p.add_argument("--max-test-cells", type=int, default=15000)
    p.add_argument("--rna-top-n", type=int, default=4000)
    p.add_argument("--lodo-donor", type=str, default=None,
                   help="If set, also report leave-one-donor-out using this donor.")
    p.add_argument("--donor-col", default="Donor")
    p.add_argument("--cluster-col", default="Level 3 Multimodal")
    p.add_argument("--skip-mlp", action="store_true")
    p.add_argument("--skip-elastic", action="store_true")
    p.add_argument("--skip-knn", action="store_true")
    p.add_argument("--skip-ridge", action="store_true")
    p.add_argument("--skip-targeted", action="store_true")
    p.add_argument("--skip-centroid", action="store_true")
    p.add_argument("--centroid-residual-alpha", type=float, default=10.0)
    p.add_argument("--centroid-no-residual", action="store_true",
                   help="Use centroid lookup only, no per-protein RNA residual head.")
    p.add_argument("--centroid-kinds", default="global_ridge",
                   help="Comma-separated list of residual kinds to fit. "
                        "Choices: global_ridge, targeted, mlp, none. Use 'none' to fit a "
                        "centroid-only candidate. Example: 'global_ridge,targeted,mlp'.")
    p.add_argument("--centroid-mlp-hidden", default="128,64")
    p.add_argument("--centroid-mlp-dropout", type=float, default=0.3)
    p.add_argument("--centroid-mlp-epochs", type=int, default=20)
    p.add_argument("--targeted-method", choices=("ridge", "elastic"), default="ridge")
    p.add_argument("--targeted-alpha", type=float, default=1.0)
    p.add_argument("--targeted-l1", type=float, default=0.5)
    p.add_argument("--targeted-share-panel-context", action="store_true",
                   help="Each per-ADT head sees the full panel-partner union (with own genes prioritised) "
                        "rather than only its own canonical partners.")
    p.add_argument("--knn-k", type=int, default=25)
    p.add_argument("--knn-pca", type=int, default=50)
    p.add_argument("--elastic-alpha", type=float, default=0.001)
    p.add_argument("--elastic-l1", type=float, default=0.5)
    p.add_argument("--ridge-alpha", type=float, default=1.0)
    p.add_argument("--mlp-hidden", type=str, default="512,256")
    p.add_argument("--mlp-epochs", type=int, default=20)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Make sure every gene referenced by the curated ADT->RNA map survives
    # the variance filter, even if it happens to fall below the top-N cutoff.
    curated_map = load_curated_adt_rna_map()
    curated_partners = sorted({g for genes in curated_map.values() for g in genes})
    extra_genes = sorted(set(_surface_marker_companions()) | set(curated_partners))

    splits: List[SplitData] = []
    if args.heldout_h5ad is not None:
        splits.append(load_split_arrays(
            args.train_h5ad,
            test_path=args.heldout_h5ad,
            max_train_cells=args.max_train_cells,
            max_test_cells=args.max_test_cells,
            rna_top_n=args.rna_top_n,
            extra_rna_genes=extra_genes,
            seed=args.seed,
        ))
    if args.lodo_donor is not None:
        splits.append(load_split_arrays(
            args.train_h5ad,
            leave_out_donor=args.lodo_donor,
            donor_col=args.donor_col,
            max_train_cells=args.max_train_cells,
            max_test_cells=args.max_test_cells,
            rna_top_n=args.rna_top_n,
            extra_rna_genes=extra_genes,
            seed=args.seed,
        ))
    if not splits:
        raise SystemExit("Provide --heldout-h5ad and/or --lodo-donor")

    # Use the first split (typically the cell-level holdout) as the deciding
    # split for picking which model to ship in the bundle.
    primary = splits[0]
    rna_genes = primary.rna_genes
    adt_names = primary.adt_names

    candidates: List[TrainedCandidate] = []
    candidate_fit_secs: Dict[str, float] = {}
    if not args.skip_targeted:
        t = time.time()
        candidates.append(_fit_targeted(
            primary,
            method=args.targeted_method,
            alpha=args.targeted_alpha,
            l1_ratio=args.targeted_l1,
            share_panel_context=args.targeted_share_panel_context,
        ))
        candidate_fit_secs[candidates[-1].name] = time.time() - t
    genz_kinds = [k.strip() for k in args.genz_kinds.split(",") if k.strip()]
    for kind in genz_kinds:
        t = time.time()
        # Each kind takes the global --genz-alpha by default, but elastic
        # benefits from a much smaller alpha because it's a per-protein
        # model fit on 1-4 features.
        kind_alpha = 0.001 if kind in ("elastic", "elastic_rank", "elastic_spline") else args.genz_alpha
        candidates.append(_fit_generalizable(
            primary,
            head_kind=kind,
            alpha=kind_alpha,
            l1_ratio=args.genz_elastic_l1,
            n_pca=args.genz_n_pca,
            n_neighbors=args.genz_n_neighbors,
            xgb_n_estimators=args.genz_xgb_n_estimators,
            xgb_max_depth=args.genz_xgb_max_depth,
            xgb_lr=args.genz_xgb_lr,
        ))
        candidate_fit_secs[candidates[-1].name] = time.time() - t

    if not args.skip_centroid:
        kinds = [k.strip() for k in args.centroid_kinds.split(",") if k.strip()]
        mlp_hidden = tuple(int(x) for x in args.centroid_mlp_hidden.split(",") if x.strip())
        for kind in kinds:
            t = time.time()
            if kind == "none":
                candidates.append(_fit_centroid(
                    primary,
                    use_residual=False,
                    residual_alpha=args.centroid_residual_alpha,
                ))
            else:
                candidates.append(_fit_centroid(
                    primary,
                    use_residual=not args.centroid_no_residual,
                    residual_alpha=args.centroid_residual_alpha,
                    residual_kind=kind,
                    mlp_hidden=mlp_hidden,
                    mlp_dropout=args.centroid_mlp_dropout,
                    mlp_epochs=args.centroid_mlp_epochs,
                ))
            candidate_fit_secs[candidates[-1].name] = time.time() - t
    if not args.skip_knn:
        t = time.time()
        candidates.append(_fit_knn(primary, n_components=args.knn_pca, n_neighbors=args.knn_k))
        candidate_fit_secs[candidates[-1].name] = time.time() - t
    if not args.skip_ridge:
        t = time.time()
        candidates.append(_fit_ridge(primary, alpha=args.ridge_alpha))
        candidate_fit_secs[candidates[-1].name] = time.time() - t
    if not args.skip_elastic:
        t = time.time()
        candidates.append(_fit_elastic(primary, alpha=args.elastic_alpha, l1_ratio=args.elastic_l1))
        candidate_fit_secs[candidates[-1].name] = time.time() - t
    if not args.skip_mlp:
        try:
            t = time.time()
            hidden = tuple(int(x) for x in args.mlp_hidden.split(",") if x.strip())
            candidates.append(_fit_mlp(primary, hidden=hidden, epochs=args.mlp_epochs))
            candidate_fit_secs[candidates[-1].name] = time.time() - t
        except Exception as exc:
            print(f"[mlp] skipped due to: {exc}")

    # Evaluate every candidate on every split.
    all_reports: List[EvalReport] = []
    for split in splits:
        for cand in candidates:
            pred, predict_secs = _predict(cand, split.X_test, split.Y_test,
                                          test_clusters=split.test_clusters)
            report = _evaluate(cand.name, split.label, split.Y_test, pred, split.adt_names,
                               fit_seconds=candidate_fit_secs.get(cand.name, 0.0),
                               predict_seconds=predict_secs)
            all_reports.append(report)
            per_path = args.output_dir / f"per_protein__{cand.name}__{_safe(split.label)}.tsv"
            report.per_protein.to_csv(per_path, sep="\t", index=False)
            print(f"[eval] {cand.name} on {split.label}: "
                  f"mean_pearson={report.mean_pearson:.3f} median={report.median_pearson:.3f} "
                  f"rmse={report.rmse:.3f}")

    summary_df = pd.DataFrame([_summary_row(r) for r in all_reports])
    summary_path = args.output_dir / "summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print("\n=== Summary ===")
    print(summary_df.to_string(index=False))
    print(f"[summary] {summary_path}")

    # Pick best candidate by mean_pearson on the primary split.
    primary_reports = [r for r in all_reports if r.split_label == primary.label]
    if not primary_reports:
        raise SystemExit("No primary-split reports were produced.")
    best = max(primary_reports, key=lambda r: r.mean_pearson)
    print(f"[winner] {best.model_name} mean_pearson={best.mean_pearson:.3f} on {best.split_label}")
    winning = next(c for c in candidates if c.name == best.model_name)

    metadata = {
        "atlas": str(args.train_h5ad),
        "heldout": str(args.heldout_h5ad) if args.heldout_h5ad else None,
        "donor_col": args.donor_col,
        "cluster_col": args.cluster_col,
        "n_train_cells": int(primary.X_train.shape[0]),
        "n_test_cells": int(primary.X_test.shape[0]),
        "n_rna_genes": int(len(rna_genes)),
        "n_adts": int(len(adt_names)),
        "primary_split": primary.label,
        "evaluation_summary": summary_df.to_dict(orient="records"),
        "winning_split_mean_pearson": best.mean_pearson,
    }
    save_bundle(
        winning,
        bundle_path=args.bundle_path,
        rna_genes=rna_genes,
        adt_names=adt_names,
        metadata=metadata,
    )

    def _find_centroid(*, residual_kind: Optional[str], use_residual: bool):
        for c in candidates:
            if not c.extra_metadata.get("_centroid_model"):
                continue
            if c.extra_metadata.get("use_residual") != use_residual:
                continue
            if use_residual and c.extra_metadata.get("residual_kind") != residual_kind:
                continue
            return c
        return None

    def _find_genz(kind: str):
        for c in candidates:
            if c.extra_metadata.get("head_kind") == kind:
                return c
        return None

    bundle_writes = [
        (args.centroid_bundle_path, _find_centroid(residual_kind="global_ridge", use_residual=True),
         "centroid_plus_global_ridge_residual"),
        (args.centroid_targeted_bundle_path, _find_centroid(residual_kind="targeted", use_residual=True),
         "centroid_plus_targeted_residual"),
        (args.centroid_mlp_bundle_path, _find_centroid(residual_kind="mlp", use_residual=True),
         "centroid_plus_mlp_residual"),
        (args.centroid_only_bundle_path, _find_centroid(residual_kind=None, use_residual=False),
         "centroid_only"),
        (args.genz_ridge_bundle_path, _find_genz("ridge"), "generalizable_ridge"),
        (args.genz_elastic_bundle_path, _find_genz("elastic"), "generalizable_elastic"),
        (args.genz_spline_bundle_path, _find_genz("spline"), "generalizable_spline"),
        (args.genz_xgb_bundle_path, _find_genz("xgboost"), "generalizable_xgboost"),
    ]
    for path, candidate, label in bundle_writes:
        if path is None:
            continue
        if candidate is None:
            print(f"[bundle:{label}] no candidate trained; skipping write to {path}.")
            continue
        save_bundle(
            candidate,
            bundle_path=path,
            rna_genes=rna_genes,
            adt_names=adt_names,
            metadata={
                **metadata,
                "selected_strategy": label,
            },
        )

    if args.targeted_bundle_path is not None:
        targeted_candidate = next(
            (c for c in candidates if c.extra_metadata.get("approach", "").startswith("targeted_")),
            None,
        )
        if targeted_candidate is None:
            print("[targeted-bundle] no targeted candidate trained; skipping write.")
        else:
            save_bundle(
                targeted_candidate,
                bundle_path=args.targeted_bundle_path,
                rna_genes=rna_genes,
                adt_names=adt_names,
                metadata={
                    **metadata,
                    "selected_for_disease_generalization": True,
                },
            )

    meta_json_path = args.output_dir / "training_metadata.json"
    with meta_json_path.open("w") as fh:
        json.dump(metadata, fh, indent=2, default=str)
    print(f"[meta] {meta_json_path}")


def _safe(text: str) -> str:
    return "".join(ch if (ch.isalnum() or ch in {"-", "_"}) else "_" for ch in str(text))


def _surface_marker_companions() -> List[str]:
    """Hard-coded list of canonical RNA partners of common ADT panel proteins.

    Including these guarantees the model has the most informative gene for
    every classical surface marker even if a particular gene happens to fall
    below the variance cutoff in the current atlas.
    """
    return sorted({
        # T/NK/B/myeloid lineage
        "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "CD19", "MS4A1", "CD22",
        "CD27", "CD28", "CD38", "CD79A", "CD79B", "IGHD", "IGHM",
        # Activation / costim
        "ICOS", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "TNFRSF4", "TNFRSF9",
        # NK / cytotoxic
        "NCAM1", "FCGR3A", "KLRD1", "KLRB1", "KLRG1", "GNLY",
        # Monocyte / DC
        "CD14", "FCGR1A", "LYZ", "ITGAX", "ITGAM", "CST3", "HLA-DRA", "HLA-DRB1",
        "CLEC4C", "IL3RA", "CLEC9A", "CD1C",
        # Stem / progenitor
        "CD34", "KIT", "PROM1", "FLT3", "MPO", "ELANE", "MME", "CD33",
        # Erythroid / megakaryocyte
        "GYPA", "GYPB", "ITGA2B", "GP1BA", "GP9", "PF4",
        # Adhesion / homing
        "SELL", "CCR7", "CXCR4", "CXCR3", "ITGB1", "ITGB7", "ITGAE", "S1PR1",
        # MHC / Treg / others
        "FOXP3", "IL7R", "IL2RA", "CD69", "ENTPD1", "NT5E", "B3GAT1", "THY1",
    })


if __name__ == "__main__":
    main()
