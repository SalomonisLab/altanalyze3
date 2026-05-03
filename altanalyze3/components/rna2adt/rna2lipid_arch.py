"""rna2lipid-style sibling for ADT prediction.

Faithful port of the architecture used by ``components.rna2lipid``:
    StandardScaler(X) -> MultiTaskElasticNetCV -> StandardScaler.inverse_transform(Y)

Two variants:

* ``Rna2LipidArchFull``  -- input is the full RNA matrix the pipeline kept
  (top-variance genes + curated companions); the multi-task fit selects
  which features matter via L1 sparsity. This is the literal sibling.
* ``Rna2LipidArchPanel`` -- input is restricted to the curated ADT panel
  partner genes (typically ~136). Closest to the "ADT-to-gene lookups only"
  requirement. Worth testing because rna2lipid's CV-tuned multi-task
  shrinkage may behave differently from the per-ADT ElasticNet we tried.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNet, MultiTaskElasticNet, MultiTaskElasticNetCV
from sklearn.preprocessing import StandardScaler

from .adt_rna_map import audit_against_genes, load_curated_adt_rna_map


class _BaseLipidArch:
    head_kind: str = "rna2lipid_arch_base"

    def __init__(
        self,
        *,
        rna_genes: Sequence[str],
        adt_names: Sequence[str],
        l1_ratio: Sequence[float] | float = (0.1, 0.3, 0.5, 0.7, 0.9),
        alphas: Optional[Sequence[float]] = None,
        cv: int = 5,
        max_iter: int = 1500,
        n_jobs: int = 4,
        feature_subset: str = "full",       # "full" or "panel"
        use_cv: bool = True,
        fixed_alpha: float = 0.01,
    ) -> None:
        self.rna_genes = list(rna_genes)
        self.adt_names = list(adt_names)
        self.l1_ratio = l1_ratio
        self.alphas = (None if alphas is None else np.asarray(alphas, dtype=float))
        self.cv = int(cv)
        self.max_iter = int(max_iter)
        self.n_jobs = int(n_jobs)
        self.feature_subset = str(feature_subset)
        self.use_cv = bool(use_cv)
        self.fixed_alpha = float(fixed_alpha)
        self._scaler_x: Optional[StandardScaler] = None
        self._scaler_y: Optional[StandardScaler] = None
        self._model: Optional[MultiTaskElasticNetCV] = None
        self._panel_indices: Optional[List[int]] = None
        self._panel_gene_names: List[str] = []
        self.n_features_in_: Optional[int] = None

    def _resolve_panel_indices(self) -> Tuple[List[int], List[str]]:
        adt_to_partners = load_curated_adt_rna_map()
        audit = audit_against_genes(self.adt_names, self.rna_genes, map_obj=adt_to_partners)
        gene_idx = {g: i for i, g in enumerate(self.rna_genes)}
        seen: set[str] = set()
        names: List[str] = []
        idx: List[int] = []
        for adt in self.adt_names:
            for g in audit[adt]["matched"]:  # type: ignore[index]
                if g in gene_idx and g not in seen:
                    seen.add(g)
                    names.append(g)
                    idx.append(gene_idx[g])
        return idx, names

    def _restrict(self, X: np.ndarray) -> np.ndarray:
        if self.feature_subset == "panel":
            if self._panel_indices is None:
                raise RuntimeError("panel indices not resolved")
            return X[:, self._panel_indices]
        return X

    def fit(self, X: np.ndarray, Y: np.ndarray) -> "_BaseLipidArch":
        self.n_features_in_ = int(X.shape[1])
        if self.feature_subset == "panel":
            idx, names = self._resolve_panel_indices()
            self._panel_indices = idx
            self._panel_gene_names = names
            print(f"[lipidarch] panel feature subset: {len(idx)} genes", flush=True)
        Xr = self._restrict(np.asarray(X, dtype=np.float32))
        Yr = np.asarray(Y, dtype=np.float32)
        self._scaler_x = StandardScaler()
        Xs = self._scaler_x.fit_transform(Xr)
        self._scaler_y = StandardScaler()
        Ys = self._scaler_y.fit_transform(Yr)

        if self.use_cv:
            l1_ratios = list(self.l1_ratio) if isinstance(self.l1_ratio, (list, tuple, np.ndarray)) else [float(self.l1_ratio)]
            print(f"[lipidarch] fitting MultiTaskElasticNetCV cv={self.cv} l1_ratios={l1_ratios} "
                  f"alphas={'auto' if self.alphas is None else 'grid'} n_features={Xs.shape[1]}",
                  flush=True)
            kwargs: Dict[str, object] = dict(
                l1_ratio=l1_ratios,
                cv=self.cv,
                max_iter=self.max_iter,
                n_jobs=self.n_jobs,
                tol=1e-3,
            )
            if self.alphas is not None:
                kwargs["alphas"] = self.alphas
            self._model = MultiTaskElasticNetCV(**kwargs)
            self._model.fit(Xs, Ys)
            print(f"[lipidarch] selected alpha={self._model.alpha_:.4g} l1_ratio={self._model.l1_ratio_:.2f}",
                  flush=True)
        else:
            # Fixed-hyperparameter MultiTaskElasticNet (no CV) — much faster.
            l1_ratio = float(self.l1_ratio[0]) if isinstance(self.l1_ratio, (list, tuple, np.ndarray)) else float(self.l1_ratio)
            print(f"[lipidarch] fitting MultiTaskElasticNet (no CV) alpha={self.fixed_alpha} "
                  f"l1_ratio={l1_ratio} n_features={Xs.shape[1]}", flush=True)
            self._model = MultiTaskElasticNet(
                alpha=self.fixed_alpha,
                l1_ratio=l1_ratio,
                max_iter=self.max_iter,
                tol=1e-3,
                random_state=0,
            )
            self._model.fit(Xs, Ys)
        return self

    def predict(self, X: np.ndarray, cluster_labels=None) -> np.ndarray:
        if self._model is None or self._scaler_x is None or self._scaler_y is None:
            raise RuntimeError("model unfit")
        Xr = self._restrict(np.asarray(X, dtype=np.float32))
        Xs = self._scaler_x.transform(Xr)
        ys = self._model.predict(Xs)
        y = self._scaler_y.inverse_transform(ys)
        return np.asarray(y, dtype=np.float32)


class Rna2LipidArchFull(_BaseLipidArch):
    head_kind = "rna2lipid_arch_full"

    def __init__(self, **kwargs) -> None:
        kwargs.setdefault("feature_subset", "full")
        super().__init__(**kwargs)


class Rna2LipidArchPanel(_BaseLipidArch):
    head_kind = "rna2lipid_arch_panel"

    def __init__(self, **kwargs) -> None:
        kwargs.setdefault("feature_subset", "panel")
        super().__init__(**kwargs)


class Rna2LipidArchPanelPerProtein:
    """rna2lipid scaling architecture but with per-ADT independent ElasticNet.

    Same as ``Rna2LipidArchPanel`` for everything except the model: instead
    of one MultiTaskElasticNet with a joint penalty across all 129 outputs,
    we fit 129 independent ElasticNet regressors in sequence. This removes
    the conservative averaging effect of multi-task shrinkage that was
    suppressing aberrant amplitude in the strict and loose variants.

    Two feature-resolution modes:

    * ``feature_source="curated"`` (default) -- each head sees the full
      panel-union (the union of curated partner genes from
      ``adt_rna_map.py``, ~136 genes), and L1 sparsity inside each head
      selects which partners to use for that protein.
    * ``feature_source="whitelist"`` -- each head sees ONLY its own
      whitelist genes from ``configs/empirical_whitelist.tsv`` (typically
      1-3 genes per head). The panel-union (input vector to ``predict``)
      is the union of whitelist genes across all ADTs (~222 genes). This
      is more biologically targeted and faster: each head fits on a
      tiny feature set that already passed the empirical correlation
      filter.
    """

    head_kind = "rna2lipid_arch_panel_per_protein"

    def __init__(
        self,
        *,
        rna_genes: Sequence[str],
        adt_names: Sequence[str],
        alpha: float = 0.01,
        l1_ratio: float = 0.5,
        max_iter: int = 2000,
        feature_source: str = "curated",
        whitelist_path: Optional[Path] = None,
    ) -> None:
        self.rna_genes = list(rna_genes)
        self.adt_names = list(adt_names)
        self.alpha = float(alpha)
        self.l1_ratio = float(l1_ratio)
        self.max_iter = int(max_iter)
        if feature_source not in {"curated", "whitelist", "whitelist_union"}:
            raise ValueError(
                f"feature_source must be 'curated', 'whitelist', or 'whitelist_union', got {feature_source!r}"
            )
        self.feature_source = feature_source
        self.whitelist_path = (whitelist_path if whitelist_path is not None
                               else Path(__file__).parent / "configs" / "empirical_whitelist.tsv")
        self._scaler_x: Optional[StandardScaler] = None
        self._scaler_y: Optional[StandardScaler] = None
        self._models: List[ElasticNet] = []
        self._panel_indices: Optional[List[int]] = None
        self._panel_gene_names: List[str] = []
        # For whitelist mode: per-ADT indices into the panel-union slice
        self._head_local_indices: Optional[List[List[int]]] = None
        self._coef_matrix: Optional[np.ndarray] = None
        self._intercept_vector: Optional[np.ndarray] = None
        self.n_features_in_: Optional[int] = None

    def _resolve_panel_indices_curated(self) -> Tuple[List[int], List[str]]:
        adt_to_partners = load_curated_adt_rna_map()
        audit = audit_against_genes(self.adt_names, self.rna_genes, map_obj=adt_to_partners)
        gene_idx = {g: i for i, g in enumerate(self.rna_genes)}
        seen: set[str] = set()
        names: List[str] = []
        idx: List[int] = []
        for adt in self.adt_names:
            for g in audit[adt]["matched"]:  # type: ignore[index]
                if g in gene_idx and g not in seen:
                    seen.add(g)
                    names.append(g)
                    idx.append(gene_idx[g])
        return idx, names

    def _resolve_features_whitelist(self) -> Tuple[List[int], List[str], List[List[str]]]:
        """Read empirical_whitelist.tsv and return (panel_idx, panel_gene_names, per_adt_genes)."""
        import pandas as pd
        if not self.whitelist_path.exists():
            raise FileNotFoundError(f"whitelist not found at {self.whitelist_path}")
        df = pd.read_csv(self.whitelist_path, sep="\t")
        # Map adt_raw -> feature_genes
        adt_to_features: Dict[str, List[str]] = {}
        for _, row in df.iterrows():
            adt_raw = str(row["adt_raw"])
            genes = [g.strip() for g in str(row["feature_genes"]).split(",") if g.strip()]
            adt_to_features[adt_raw] = genes
        # Also support clean-name lookup as a fallback
        adt_clean_to_features = {str(r["adt_clean"]): [g.strip() for g in str(r["feature_genes"]).split(",") if g.strip()]
                                 for _, r in df.iterrows()}

        gene_idx = {g: i for i, g in enumerate(self.rna_genes)}
        seen: set[str] = set()
        panel_names: List[str] = []
        panel_idx: List[int] = []
        per_adt_genes: List[List[str]] = []
        missing_genes: set[str] = set()
        n_no_entry = 0
        for adt in self.adt_names:
            features = adt_to_features.get(adt) or adt_clean_to_features.get(adt)
            if features is None:
                # No whitelist row found — fall back to empty (head will fit on a single-gene
                # surrogate in the panel union; the head won't have meaningful features).
                n_no_entry += 1
                per_adt_genes.append([])
                continue
            kept: List[str] = []
            for g in features:
                if g in gene_idx:
                    kept.append(g)
                    if g not in seen:
                        seen.add(g)
                        panel_names.append(g)
                        panel_idx.append(gene_idx[g])
                else:
                    missing_genes.add(g)
            per_adt_genes.append(kept)
        print(f"[whitelist] {len(panel_idx)} unique features across {len(self.adt_names)} ADTs "
              f"({n_no_entry} ADTs had no whitelist entry, "
              f"{len(missing_genes)} whitelist genes missing from RNA matrix)",
              flush=True)
        return panel_idx, panel_names, per_adt_genes

    def fit(self, X: np.ndarray, Y: np.ndarray) -> "Rna2LipidArchPanelPerProtein":
        self.n_features_in_ = int(X.shape[1])
        if self.feature_source == "whitelist":
            panel_idx, panel_names, per_adt_genes = self._resolve_features_whitelist()
            # For each head we also need its local-index slice into the panel-union vector.
            panel_pos = {g: i for i, g in enumerate(panel_names)}
            self._head_local_indices = [
                [panel_pos[g] for g in feats] for feats in per_adt_genes
            ]
        elif self.feature_source == "whitelist_union":
            # Same panel union as 'whitelist' (208 genes), but every head sees the
            # full union as input — L1 inside the per-head ElasticNet decides
            # which features to use. Recovers the dropout-resilience of the
            # curated variant while broadening the input set to empirically
            # informative genes, not just curated partners.
            panel_idx, panel_names, _ = self._resolve_features_whitelist()
            self._head_local_indices = None
        else:
            panel_idx, panel_names = self._resolve_panel_indices_curated()
            self._head_local_indices = None
        self._panel_indices = panel_idx
        self._panel_gene_names = panel_names
        print(f"[lipidarch-per-protein] feature_source={self.feature_source} "
              f"panel size: {len(panel_idx)} genes", flush=True)
        Xr = np.asarray(X, dtype=np.float32)[:, panel_idx]
        Yr = np.asarray(Y, dtype=np.float32)
        self._scaler_x = StandardScaler()
        Xs = self._scaler_x.fit_transform(Xr)
        self._scaler_y = StandardScaler()
        Ys = self._scaler_y.fit_transform(Yr)

        print(f"[lipidarch-per-protein] fitting {Ys.shape[1]} per-ADT ElasticNets "
              f"alpha={self.alpha} l1_ratio={self.l1_ratio} max_iter={self.max_iter}",
              flush=True)
        self._models = []
        self._coef_matrix = None
        self._intercept_vector = None
        for j in range(Ys.shape[1]):
            m = ElasticNet(alpha=self.alpha, l1_ratio=self.l1_ratio,
                           max_iter=self.max_iter, tol=1e-3, random_state=0)
            if self._head_local_indices is not None:
                local_idx = self._head_local_indices[j]
                if not local_idx:
                    # No features assigned — fit on the global mean (intercept only).
                    fallback = np.zeros((Xs.shape[0], 1), dtype=np.float32)
                    m.fit(fallback, Ys[:, j])
                else:
                    m.fit(Xs[:, local_idx], Ys[:, j])
            else:
                m.fit(Xs, Ys[:, j])
            self._models.append(m)
        print(f"[lipidarch-per-protein] done {len(self._models)} heads", flush=True)
        return self

    def predict(self, X: np.ndarray, cluster_labels=None) -> np.ndarray:
        if not self._models or self._scaler_x is None or self._scaler_y is None:
            raise RuntimeError("model unfit")
        Xr = np.asarray(X, dtype=np.float32)[:, self._panel_indices]
        Xs = self._scaler_x.transform(Xr)
        if self._head_local_indices is None:
            coef_matrix = getattr(self, "_coef_matrix", None)
            intercept_vector = getattr(self, "_intercept_vector", None)
            if coef_matrix is None or intercept_vector is None:
                coef_matrix = np.column_stack(
                    [np.asarray(m.coef_, dtype=np.float32).ravel() for m in self._models]
                ).astype(np.float32, copy=False)
                intercept_vector = np.asarray(
                    [float(m.intercept_) for m in self._models],
                    dtype=np.float32,
                )
                self._coef_matrix = coef_matrix
                self._intercept_vector = intercept_vector
            ys = Xs @ coef_matrix + intercept_vector
        else:
            ys = np.zeros((Xs.shape[0], len(self._models)), dtype=np.float32)
            for j, m in enumerate(self._models):
                local_idx = self._head_local_indices[j]
                if not local_idx:
                    ys[:, j] = m.predict(np.zeros((Xs.shape[0], 1), dtype=np.float32)).astype(np.float32)
                else:
                    ys[:, j] = m.predict(Xs[:, local_idx]).astype(np.float32)
        y = self._scaler_y.inverse_transform(ys)
        return np.asarray(y, dtype=np.float32)
