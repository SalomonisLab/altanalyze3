"""AML RNA -> mass-spec abundance imputation engine (per-molecule regularized linear model).

Shared, modality-agnostic core for the AML ``rna2metabolite`` and ``rna2lipid.aml``
components. Each component ships its own copy of this file so it is fully self-contained
and independent of the (non-AML) lung ``rna2lipid`` code.

A bundle is an (optionally gzip-compressed) pickle **dict** (plain numpy arrays only — no
custom classes to unpickle)::

    X_columns : list[str]            model gene symbols (z-score / weight order)
    Y_columns : list[str]            molecule names (targets)
    mu, sd    : float[n_genes]       per-gene training mean / std
    sel_idx   : list[int[]]          per molecule, indices into X_columns
    coef      : list[float[]]        per molecule, ridge weights (z-score space)
    intercept : float[n_targets]
    metadata  : dict                 provenance (estimator, normalization, held-out perf)

Prediction is per-molecule and linear::

    z      = (CP10k+log1p(counts) - mu) / sd          # genes absent from the input -> z=0
    yhat_t = intercept_t + coef_t . z[sel_idx_t]

i.e. one independent model per molecule, evaluated on the sample's own expression. There
is no per-sample neighbour search; on held-out CV this beats k-NN retrieval on the same
genes (metabolite 0.268 vs 0.223, lipid 0.355 vs 0.308 median Spearman).
"""
from __future__ import annotations

import gzip
import pickle
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PredictionResult:
    predictions: pd.DataFrame          # (n_samples, n_molecules)
    summary: Dict[str, object]


def _clean(values) -> List[str]:
    return [str(v).strip() for v in values]


def cp10k_log1p(counts: np.ndarray, totals: Optional[np.ndarray] = None) -> np.ndarray:
    """CP10k + log1p of a (rows x genes) count matrix. ``totals`` overrides the per-row
    library size (use the all-input-gene library size, matching training)."""
    counts = np.asarray(counts, dtype=np.float64)
    if totals is None:
        totals = counts.sum(axis=1)
    tot = np.asarray(totals, dtype=np.float64).reshape(-1, 1).copy()
    tot[tot == 0] = 1.0
    return np.log1p(counts / tot * 1e4)


def _open(path: Path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rb") if magic == b"\x1f\x8b" else open(path, "rb")


class PerTargetImputeBundle:
    """Loads a bundle dict and imputes molecule abundance from RNA counts/expression."""

    def __init__(self, *, bundle_path, X_columns, Y_columns, mu, sd,
                 sel_idx, coef, intercept, metadata=None):
        self.bundle_path = Path(bundle_path)
        self.input_genes = tuple(_clean(X_columns))
        self.targets = tuple(_clean(Y_columns))
        self.mu = np.asarray(mu, np.float64)
        self.sd = np.asarray(sd, np.float64).copy(); self.sd[self.sd < 1e-8] = 1.0
        self.sel_idx = [np.asarray(s, np.int64).ravel() for s in sel_idx]
        self.coef = [np.asarray(c, np.float64).ravel() for c in coef]
        self.intercept = np.asarray(intercept, np.float64)
        self.metadata = dict(metadata or {})
        self._gene_to_idx = {g: i for i, g in enumerate(self.input_genes)}
        self._W = None
        self._valid = None

    # ------------------------------------------------------------------ load
    @classmethod
    def load(cls, bundle_path) -> "PerTargetImputeBundle":
        bundle_path = Path(bundle_path)
        with _open(bundle_path) as fh:
            b = pickle.load(fh)
        req = {"X_columns", "Y_columns", "mu", "sd", "sel_idx", "coef", "intercept"}
        miss = req.difference(b)
        if miss:
            raise ValueError(f"imputation bundle is missing keys: {sorted(miss)}")
        return cls(bundle_path=bundle_path, X_columns=b["X_columns"], Y_columns=b["Y_columns"],
                   mu=b["mu"], sd=b["sd"], sel_idx=b["sel_idx"], coef=b["coef"],
                   intercept=b["intercept"], metadata=b.get("metadata"))

    def model_info(self) -> Dict[str, object]:
        info = {"bundle_path": str(self.bundle_path), "n_input_genes": len(self.input_genes),
                "n_molecules": len(self.targets), "model_class": "PerTargetRidge"}
        for k in ("modality", "estimator", "normalization", "gene_filter", "nfeat", "alpha",
                  "n_train_cases", "heldout_median_spearman", "n_imputable_sp_gt_0p3"):
            if k in self.metadata:
                info[k] = self.metadata[k]
        return info

    # --------------------------------------------------------------- weights
    def _weight_matrix(self):
        if self._W is None:
            from scipy.sparse import csr_matrix
            rows, cols, data = [], [], []
            valid = np.zeros(len(self.targets), bool)
            for ti, (idx, c) in enumerate(zip(self.sel_idx, self.coef)):
                if idx.size and np.isfinite(self.intercept[ti]) and idx.size == c.size:
                    rows.append(np.full(idx.size, ti)); cols.append(idx); data.append(c)
                    valid[ti] = True
            self._W = csr_matrix(
                (np.concatenate(data) if data else np.zeros(0),
                 (np.concatenate(rows) if rows else np.zeros(0, int),
                  np.concatenate(cols) if cols else np.zeros(0, int))),
                shape=(len(self.targets), len(self.input_genes)))
            self._valid = valid
        return self._W, self._valid

    # --------------------------------------------------------------- align
    def _align(self, matrix, gene_labels):
        """Align a (rows x genes_in) matrix to the model gene order. Duplicate input genes
        are summed; the CP10k library size is taken over ALL input genes (training
        convention). Returns aligned counts, per-row totals, present-gene mask, n matched."""
        import scipy.sparse as _sp
        pos = defaultdict(list)
        for i, g in enumerate(_clean(gene_labels)):
            if g:
                pos[g].append(i)
        issparse = _sp.issparse(matrix)
        M = matrix.tocsc() if issparse else np.asarray(matrix, np.float64)   # avoid full densify
        totals = np.asarray(M.sum(axis=1)).reshape(-1).astype(np.float64)
        aligned = np.zeros((M.shape[0], len(self.input_genes)), np.float64)
        present = np.zeros(len(self.input_genes), bool); matched = 0
        for g, mi in self._gene_to_idx.items():
            src = pos.get(g)
            if not src:
                continue
            matched += 1; present[mi] = True
            col = M[:, src]
            aligned[:, mi] = (np.asarray(col.sum(axis=1)).reshape(-1) if issparse
                              else (col.sum(axis=1) if len(src) > 1 else col[:, 0]))
        return aligned, totals, present, matched

    def _predict_from_aligned(self, aligned, totals, present, *, normalized):
        Xn = aligned if normalized else cp10k_log1p(aligned, totals)
        Z = (Xn - self.mu[None, :]) / self.sd[None, :]
        Z[:, ~present] = 0.0                                   # genes absent in input -> training mean
        W, valid = self._weight_matrix()
        Yhat = np.asarray((W @ Z.T).T) + self.intercept[None, :]
        Yhat[:, ~valid] = np.nan
        return Yhat, valid

    # ------------------------------------------------------------- dataframe
    def predict_from_dataframe(self, expression: pd.DataFrame, *, normalized: bool = False) -> PredictionResult:
        """expression: rows = samples/pseudobulks, columns = gene symbols."""
        if expression.empty:
            raise ValueError("expression dataframe is empty")
        df = expression.copy()
        df.index = _clean(df.index); df.columns = _clean(df.columns)
        df = df.T.groupby(level=0).sum().T
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        aligned, totals, present, matched = self._align(df.to_numpy(), df.columns)
        if matched == 0:
            raise ValueError("no model feature genes were found in the dataframe")
        Yhat, valid = self._predict_from_aligned(aligned, totals, present, normalized=normalized)
        return self._wrap(Yhat, pd.Index(df.index), matched, "dataframe", valid)

    # ----------------------------------------------------------------- adata
    def _group_counts(self, adata, groupby, layer, gene_symbol_col, min_cells):
        gene_labels = (adata.var[gene_symbol_col].astype(str).tolist()
                       if gene_symbol_col else _clean(adata.var_names))
        X = adata.layers[layer] if layer else adata.X
        kept = None
        if groupby is None:
            counts = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
            index = pd.Index(_clean(adata.obs_names))
            obs = adata.obs.copy()
        else:
            if groupby not in adata.obs.columns:
                raise KeyError(f"obs column {groupby!r} not found")
            labels = adata.obs[groupby].astype(str).values
            groups, counts, kept = [], [], {}
            for g in pd.unique(labels):
                m = labels == g
                if m.sum() < min_cells:
                    continue
                sub = X[m]; sub = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
                counts.append(sub.sum(0)); groups.append(g); kept[str(g)] = int(m.sum())
            counts = np.vstack(counts); index = pd.Index(groups, name=groupby)
            obs = pd.DataFrame(index=index)
        return counts, gene_labels, index, obs, kept

    def predict_from_adata(self, adata, *, groupby=None, layer=None, gene_symbol_col=None,
                           min_cells=1, normalized=False) -> PredictionResult:
        counts, gene_labels, index, _obs, kept = self._group_counts(
            adata, groupby, layer, gene_symbol_col, min_cells)
        aligned, totals, present, matched = self._align(counts, gene_labels)
        if matched == 0:
            raise ValueError("no model feature genes were found in the AnnData")
        Yhat, valid = self._predict_from_aligned(aligned, totals, present, normalized=normalized)
        res = self._wrap(Yhat, index, matched, "adata", valid)
        if kept is not None:
            res.summary["n_cells_per_group"] = kept
        return res

    def predict_from_h5ad(self, h5ad_path, **kw) -> PredictionResult:
        import anndata as ad
        return self.predict_from_adata(ad.read_h5ad(str(h5ad_path)), **kw)

    def impute_anndata(self, adata, *, groupby=None, layer=None, gene_symbol_col=None,
                       min_cells=1, normalized=False, carry_obs=True):
        """Return an AnnData of imputed abundance (samples x molecules); obs is copied from
        the source when ``groupby is None`` (rows are already pseudobulks)."""
        import anndata as ad
        counts, gene_labels, index, obs, kept = self._group_counts(
            adata, groupby, layer, gene_symbol_col, min_cells)
        aligned, totals, present, matched = self._align(counts, gene_labels)
        if matched == 0:
            raise ValueError("no model feature genes were found in the AnnData")
        Yhat, valid = self._predict_from_aligned(aligned, totals, present, normalized=normalized)
        keep = np.where(valid)[0]
        var = pd.DataFrame(index=pd.Index([self.targets[i] for i in keep], name="molecule"))
        for col in ("heldout_spearman", "heldout_r2", "feature_genes"):
            vals = self.metadata.get("per_target", {}).get(col)
            if vals is not None:
                var[col] = [vals.get(self.targets[i], np.nan) for i in keep]
        out = ad.AnnData(X=Yhat[:, keep].astype(np.float32),
                         obs=(obs if carry_obs else pd.DataFrame(index=index)), var=var)
        out.obs_names = [str(x) for x in index]
        out.uns["imputation"] = {**self.model_info(),
                                  "matched_genes": int(matched), "n_molecules": int(valid.sum())}
        if kept is not None:
            out.obs["n_cells_group"] = [kept.get(str(x), np.nan) for x in index]
        return out

    # ----------------------------------------------------------------- wrap
    def _wrap(self, Yhat, index, matched, kind, valid) -> PredictionResult:
        df = pd.DataFrame(np.asarray(Yhat), index=index, columns=list(self.targets))
        summary = {
            "bundle_path": str(self.bundle_path), "input_kind": kind,
            "n_samples": int(df.shape[0]), "matched_genes": int(matched),
            "missing_genes": len(self.input_genes) - int(matched),
            "model_gene_count": len(self.input_genes),
            "n_molecules": len(self.targets), "n_molecules_imputed": int(np.asarray(valid).sum()),
            "modality": self.metadata.get("modality"),
        }
        return PredictionResult(predictions=df, summary=summary)


def load_bundle(bundle_path) -> PerTargetImputeBundle:
    return PerTargetImputeBundle.load(bundle_path)
