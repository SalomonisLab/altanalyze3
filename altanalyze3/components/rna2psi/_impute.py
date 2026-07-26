"""AML RNA -> splicing PSI imputation engine (per-event regularized linear model).

Self-contained inference core for the ``rna2psi`` component, mirroring
``components.rna2metabolite._impute`` so the cellHarmony-web plumbing can register a
PSI-imputation modality on the same rails. One independent model per splicing event
(AltAnalyze MultiPath-PSI ``UID``); each model uses at most a handful of predictor genes.

Bundle = an (optionally gzip-compressed) pickle **dict** (plain numpy arrays only)::

    X_columns : list[str]            model feature genes (Ensembl gene IDs, z-score order)
    Y_columns : list[str]            splicing-event UIDs (targets)
    mu, sd    : float[n_genes]       per-gene training mean / std (on cp10k+log1p features)
    sel_idx   : list[int[]]          per event, <=5 indices into X_columns
    coef      : list[float[]]        per event, EN/Ridge weights (z-score space)
    intercept : float[n_targets]
    metadata  : dict                 provenance (per-event estimator + held-out CV perf)

Prediction is per-event and linear, then bounded to the PSI range::

    z      = (cp10k+log1p(counts) - mu) / sd          # genes absent from input -> z=0
    psi_e  = clip(intercept_e + coef_e . z[sel_idx_e], 0, 1)

Feature genes are matched to the input by **Ensembl gene ID** (the native key of the
Leucegene Kallisto counts). Inputs start from counts and are normalized here exactly as at
training (cp10k+log1p on the all-gene library size, then z-score with training mu/sd), so
raw counts (or cellHarmony-lite scaled values via ``normalized=True``) can be passed
straight through.
"""
from __future__ import annotations

import gzip
import pickle
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PredictionResult:
    predictions: pd.DataFrame          # (n_samples, n_events)  PSI in [0,1], NaN where unimputable
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


class PerEventPSIBundle:
    """Loads a bundle dict and imputes splicing-event PSI from RNA counts/expression."""

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
    def load(cls, bundle_path) -> "PerEventPSIBundle":
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
                "n_events": len(self.targets), "model_class": "PerEventElasticNetOrRidge"}
        for k in ("modality", "gene_id", "normalization", "gene_filter", "max_features",
                  "n_candidates", "n_train_samples", "heldout_median_spearman",
                  "n_imputable_sp_gt_0p3"):
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
        are summed; the cp10k library size is taken over ALL input genes (training
        convention). Returns aligned counts, per-row totals, present-gene mask, n matched."""
        import scipy.sparse as _sp
        pos = defaultdict(list)
        for i, g in enumerate(_clean(gene_labels)):
            if g:
                pos[g].append(i)
        issparse = _sp.issparse(matrix)
        M = matrix.tocsc() if issparse else np.asarray(matrix, np.float64)
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
        np.clip(Yhat, 0.0, 1.0, out=Yhat)                     # PSI is a bounded ratio
        Yhat[:, ~valid] = np.nan
        return Yhat, valid

    # ------------------------------------------------------------- dataframe
    def predict_from_dataframe(self, expression: pd.DataFrame, *, gene_axis: str = "rows",
                               normalized: bool = False) -> PredictionResult:
        """expression: an RNA matrix. ``gene_axis='rows'`` (default) means genes are the
        index and samples the columns (native AltAnalyze / Kallisto counts layout);
        ``gene_axis='columns'`` means the transpose (samples x genes)."""
        if expression.empty:
            raise ValueError("expression dataframe is empty")
        df = expression.copy()
        if gene_axis == "rows":
            df = df.T
        elif gene_axis != "columns":
            raise ValueError("gene_axis must be 'rows' or 'columns'")
        df.index = _clean(df.index); df.columns = _clean(df.columns)   # rows=samples, cols=genes
        df = df.T.groupby(level=0).sum().T
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        aligned, totals, present, matched = self._align(df.to_numpy(), df.columns)
        if matched == 0:
            raise ValueError("no model feature genes were found in the input")
        Yhat, valid = self._predict_from_aligned(aligned, totals, present, normalized=normalized)
        return self._wrap(Yhat, pd.Index(df.index), matched, "dataframe", valid)

    def predict_from_file(self, path, *, sep="\t", gene_axis="rows", normalized=False) -> PredictionResult:
        df = pd.read_csv(path, sep=sep, index_col=0)
        return self.predict_from_dataframe(df, gene_axis=gene_axis, normalized=normalized)

    # ----------------------------------------------------------------- adata
    def predict_from_adata(self, adata, *, groupby=None, layer=None, gene_symbol_col=None,
                           min_cells=1, normalized=False) -> PredictionResult:
        gene_labels = (adata.var[gene_symbol_col].astype(str).tolist()
                       if gene_symbol_col else _clean(adata.var_names))
        X = adata.layers[layer] if layer else adata.X
        if groupby is None:
            counts = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
            index = pd.Index(_clean(adata.obs_names))
        else:
            if groupby not in adata.obs.columns:
                raise KeyError(f"obs column {groupby!r} not found")
            labels = adata.obs[groupby].astype(str).values
            groups, rows = [], []
            for g in pd.unique(labels):
                m = labels == g
                if m.sum() < min_cells:
                    continue
                sub = X[m]; sub = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
                rows.append(sub.sum(0)); groups.append(g)
            counts = np.vstack(rows); index = pd.Index(groups, name=groupby)
        aligned, totals, present, matched = self._align(counts, gene_labels)
        if matched == 0:
            raise ValueError("no model feature genes were found in the AnnData")
        Yhat, valid = self._predict_from_aligned(aligned, totals, present, normalized=normalized)
        return self._wrap(Yhat, index, matched, "adata", valid)

    # ----------------------------------------------------------------- wrap
    def _wrap(self, Yhat, index, matched, kind, valid) -> PredictionResult:
        df = pd.DataFrame(np.asarray(Yhat), index=index, columns=list(self.targets))
        summary = {
            "bundle_path": str(self.bundle_path), "input_kind": kind,
            "n_samples": int(df.shape[0]), "matched_genes": int(matched),
            "missing_genes": len(self.input_genes) - int(matched),
            "model_gene_count": len(self.input_genes),
            "n_events": len(self.targets), "n_events_imputed": int(np.asarray(valid).sum()),
            "modality": self.metadata.get("modality", "psi"),
        }
        return PredictionResult(predictions=df, summary=summary)


def load_bundle(bundle_path) -> PerEventPSIBundle:
    return PerEventPSIBundle.load(bundle_path)
