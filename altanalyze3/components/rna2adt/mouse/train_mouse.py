"""Train the mouse rna2adt bundle: per-protein ElasticNet on whitelist union.

Mirrors the production human pipeline (whitelist_union variant of
``Rna2LipidArchPanelPerProtein``) for the murine atlas. Saves the trained
bundle to ``rna2adt/mouse/rna2adt_mm_bundle.pkl`` for cellHarmony-web.
"""

from __future__ import annotations

import argparse
import pickle
import time
from pathlib import Path
from typing import List

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr
from sklearn.linear_model import ElasticNet
from sklearn.preprocessing import StandardScaler

from .adt_mgi_map import strip_prefix


_ADT_PREFIX = "ADT-"
_ISOTYPE_TOKENS = ("Isotype_Ctrl", "isotype_Ctrl", "isotype_ctrl")


def _is_isotype(name: str) -> bool:
    return any(tok in name for tok in _ISOTYPE_TOKENS)


def _load_whitelist(path: Path) -> dict:
    df = pd.read_csv(path, sep="\t")
    out: dict = {}
    for _, r in df.iterrows():
        adt = str(r["adt_raw"])
        genes = [g.strip() for g in str(r["feature_genes"]).split(",") if g.strip()]
        out[adt] = genes
    return out


class MouseRna2AdtModel:
    """Per-ADT ElasticNet on the panel-union of mouse whitelist genes.

    Architecturally identical to the human ``Rna2LipidArchPanelPerProtein``
    with ``feature_source="whitelist_union"``: every per-ADT head sees the
    full panel-union as input; L1 sparsity inside each head selects the
    informative subset for that ADT.

    Attributes set after ``.fit``:
        rna_genes        - input gene order (the full panel-union)
        adt_names        - output ADT order (with "ADT-" prefix preserved)
        scaler_x         - fitted StandardScaler over rna_genes
        scaler_y         - fitted StandardScaler over adt_names
        models           - list of ElasticNet, len == len(adt_names)
    """

    head_kind = "rna2adt_mm_panel_per_protein_whitelist_union"

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

    def fit(self, X: np.ndarray, Y: np.ndarray,
            *, rna_genes, adt_names) -> "MouseRna2AdtModel":
        self.rna_genes = list(rna_genes)
        self.adt_names = list(adt_names)
        self.n_features_in_ = int(X.shape[1])
        Xr = np.asarray(X, dtype=np.float32)
        Yr = np.asarray(Y, dtype=np.float32)
        self.scaler_x = StandardScaler()
        Xs = self.scaler_x.fit_transform(Xr)
        self.scaler_y = StandardScaler()
        Ys = self.scaler_y.fit_transform(Yr)
        print(f"[train] fitting {Ys.shape[1]} per-ADT ElasticNets "
              f"alpha={self.alpha} l1_ratio={self.l1_ratio} "
              f"n_features={Xs.shape[1]}", flush=True)
        self.models = []
        self._coef_matrix = None
        self._intercept_vector = None
        for j in range(Ys.shape[1]):
            m = ElasticNet(alpha=self.alpha, l1_ratio=self.l1_ratio,
                           max_iter=self.max_iter, tol=1e-3, random_state=0)
            m.fit(Xs, Ys[:, j])
            self.models.append(m)
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
            intercept_vector = np.asarray(
                [float(m.intercept_) for m in self.models],
                dtype=np.float32,
            )
            self._coef_matrix = coef_matrix
            self._intercept_vector = intercept_vector
        ys = Xs @ coef_matrix + intercept_vector
        y = self.scaler_y.inverse_transform(ys)
        return np.asarray(y, dtype=np.float32)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--adata", type=Path,
                   default=Path('/Users/saljh8/Dropbox/Manuscripts/InProgress/ALRG_Paper_2022/Azimuth-Tests/adata_combined_60k_rna_adt.h5ad'))
    p.add_argument("--whitelist", type=Path,
                   default=Path(__file__).parent / "configs" / "empirical_whitelist.tsv")
    p.add_argument("--bundle-out", type=Path,
                   default=Path(__file__).parent / "rna2adt_mm_bundle.pkl")
    p.add_argument("--max-train-cells", type=int, default=40000)
    p.add_argument("--max-test-cells", type=int, default=15000)
    p.add_argument("--alpha", type=float, default=0.01)
    p.add_argument("--l1-ratio", type=float, default=0.5)
    p.add_argument("--max-iter", type=int, default=2000)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    T0 = time.time()
    def step(msg: str) -> None:
        print(f"[{time.time()-T0:6.1f}s] {msg}", flush=True)

    step(f"loading {args.adata}")
    a = ad.read_h5ad(args.adata)
    print(f"  {a.n_obs} cells, {a.n_vars} vars")

    var_names = np.array([str(v) for v in a.var_names])
    is_adt = np.array([s.startswith(_ADT_PREFIX) for s in var_names])
    is_isotype = np.array([_is_isotype(s) for s in var_names])
    is_adt_kept = is_adt & ~is_isotype
    rna_idx_full = np.where(~is_adt)[0]
    adt_idx = np.where(is_adt_kept)[0]
    adt_names_raw = var_names[is_adt_kept].tolist()
    rna_full = var_names[~is_adt].tolist()
    print(f"  {len(adt_names_raw)} ADTs (kept), {len(rna_full)} RNA genes")

    step("loading whitelist")
    wl = _load_whitelist(args.whitelist)
    panel_union: List[str] = []
    seen: set = set()
    n_no_entry = 0
    for adt in adt_names_raw:
        feats = wl.get(adt) or []
        if not feats:
            n_no_entry += 1
        for g in feats:
            if g not in seen:
                seen.add(g)
                panel_union.append(g)
    rna_set = set(rna_full)
    panel_union = [g for g in panel_union if g in rna_set]
    print(f"  {len(panel_union)} panel-union genes (no_entry ADTs={n_no_entry})")

    rna_pos = {g: i for i, g in enumerate(rna_full)}
    panel_idx = [rna_pos[g] for g in panel_union]

    step("subsampling and densifying matrices")
    rng = np.random.default_rng(args.seed)
    n_keep = min(a.n_obs, args.max_train_cells + args.max_test_cells)
    sel = np.sort(rng.choice(a.n_obs, n_keep, replace=False))
    n_train = min(args.max_train_cells, n_keep - 1)
    train_sel = sel[:n_train]
    test_sel = sel[n_train:]

    rna_full_idx = rna_idx_full
    panel_full_idx = rna_full_idx[panel_idx]

    def _slice(cells, cols):
        sub = a.X[cells][:, cols]
        if sp.issparse(sub):
            sub = np.asarray(sub.todense())
        return np.asarray(sub, dtype=np.float32)

    X_train = _slice(train_sel, panel_full_idx)
    Y_train = _slice(train_sel, adt_idx)
    X_test  = _slice(test_sel,  panel_full_idx)
    Y_test  = _slice(test_sel,  adt_idx)
    print(f"  train: X {X_train.shape} Y {Y_train.shape}")
    print(f"  test:  X {X_test.shape}  Y {Y_test.shape}")

    step("fitting model")
    model = MouseRna2AdtModel(alpha=args.alpha, l1_ratio=args.l1_ratio,
                              max_iter=args.max_iter)
    model.fit(X_train, Y_train, rna_genes=panel_union, adt_names=adt_names_raw)

    step("evaluating on cell-holdout test split")
    pred = model.predict(X_test)
    pearsons = []
    spearmans = []
    for j in range(Y_test.shape[1]):
        yt, yp = Y_test[:, j], pred[:, j]
        if yt.std() < 1e-9 or yp.std() < 1e-9:
            pearsons.append(np.nan); spearmans.append(np.nan); continue
        pearsons.append(float(np.corrcoef(yt, yp)[0, 1]))
        try:
            spearmans.append(float(spearmanr(yt, yp).correlation))
        except Exception:
            spearmans.append(np.nan)
    pearsons = np.array(pearsons); spearmans = np.array(spearmans)
    print(f"\n=== cell-holdout (n_test={Y_test.shape[0]}, n_adts={Y_test.shape[1]}) ===")
    print(f"mean Pearson:  {np.nanmean(pearsons):.3f}   median: {np.nanmedian(pearsons):.3f}")
    print(f"mean Spearman: {np.nanmean(spearmans):.3f}  median: {np.nanmedian(spearmans):.3f}")
    print(f"valid ADTs: {np.isfinite(pearsons).sum()} / {len(pearsons)}")

    step("saving bundle")
    # Strip the ADT- prefix and anti_<species>_<species>_ wrapper from output
    # names so downstream consumers see clean marker names ("CD19", "Ly_6C")
    # instead of "ADT-CD19", "ADT-anti_mouse_human_CD11b". The model's internal
    # adt_names is also patched so model.predict returns the clean names.
    adt_names_clean_out = [strip_prefix(n) for n in adt_names_raw]
    if len(set(adt_names_clean_out)) != len(adt_names_clean_out):
        raise RuntimeError("ADT name collision after stripping prefixes; investigate")
    model.adt_names = list(adt_names_clean_out)

    bundle = {
        "model": model,
        "scaler_x": None,  # handled inside model
        "scaler_y": None,
        "X_columns": panel_union,
        "Y_columns": adt_names_clean_out,
        "metadata": {
            "approach": "rna2adt_mm_panel_per_protein_whitelist_union",
            "feature_source": "whitelist_union",
            "label": "rna2adt_mm_bundle",
            "adt_name_format": "clean (ADT- prefix and anti_<species> wrapper stripped)",
            "alpha": args.alpha,
            "l1_ratio": args.l1_ratio,
            "n_train_cells": int(X_train.shape[0]),
            "n_test_cells": int(X_test.shape[0]),
            "n_panel_union": int(len(panel_union)),
            "n_adts": int(len(adt_names_raw)),
            "evaluation_summary": {
                "cell_holdout": {
                    "mean_pearson": float(np.nanmean(pearsons)),
                    "median_pearson": float(np.nanmedian(pearsons)),
                    "mean_spearman": float(np.nanmean(spearmans)),
                    "median_spearman": float(np.nanmedian(spearmans)),
                }
            },
        },
    }
    args.bundle_out.parent.mkdir(parents=True, exist_ok=True)
    with args.bundle_out.open("wb") as fh:
        pickle.dump(bundle, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"  saved {args.bundle_out}")

    # Also write per-ADT metrics table
    metrics_path = args.bundle_out.parent / "rna2adt_mm_per_adt_metrics.tsv"
    pd.DataFrame({
        "adt": adt_names_raw,
        "adt_clean": [strip_prefix(n) for n in adt_names_raw],
        "pearson": pearsons,
        "spearman": spearmans,
    }).to_csv(metrics_path, sep="\t", index=False)
    print(f"  per-ADT metrics: {metrics_path}")

    step("done")


if __name__ == "__main__":
    main()
