"""Inference API for predicting gene-regulatory-network (GRN) connection scores
from RNA, mirroring ``components.rna2adt.api`` / ``components.rna2lipid.api`` so
the cellHarmony-web pipeline can register a GRN-imputation modality on the same
plumbing.

A bundle is a (optionally gzip-compressed) pickle with::

    model      : KnnGrnRegressor with ``.predict(X)`` -> (n_pseudobulks, n_edges)
    scaler_x   : None (the kNN model standardizes internally)
    scaler_y   : None
    X_columns  : RNA feature gene symbols (TF u target genes)
    Y_columns  : GRN edge ids "TF|Gene"
    metadata   : training provenance + per-reference-row labels

A GRN is a *pseudobulk-level* quantity (it summarizes regulation across many
cells of one cell state), so — unlike the per-cell ADT/lipid models — inference
aggregates **expression** for each group into a pseudobulk, normalizes it
(CP10k + log1p), and predicts one GRN vector per group. Because the model is a
retrieval average, that is the statistically correct unit of prediction.
"""
from __future__ import annotations

import gzip
import pickle
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

BUNDLE_DIR = Path(__file__).parent

#: Named references shipped with the module. Each maps a tissue/atlas name to the
#: bundle trained on that atlas's paired (RNA pseudobulk, GRN) tables. ``leukemia``
#: stays the default so existing callers are unaffected.
REFERENCE_BUNDLES: Dict[str, Path] = {
    "leukemia": BUNDLE_DIR / "rna2grn_bundle.pkl.gz",
    "lung": BUNDLE_DIR / "rna2grn_lung_bundle.pkl.gz",
}
DEFAULT_REFERENCE = "leukemia"
DEFAULT_BUNDLE_PATH = REFERENCE_BUNDLES[DEFAULT_REFERENCE]

#: How a bundle's training pseudobulks were formed from single cells. A bundle
#: declares its own in ``metadata["pseudobulk_statistic"]``; inference must build
#: the query pseudobulk the same way or the model sees an input space it never saw
#: in training.
#:   sum_counts                    log1p(CP10k of the SUMMED counts)  [leukemia]
#:   mean_over_cells_of_log1p_cp10k  MEAN over cells of log1p(CP10k)   [lung]
PSEUDOBULK_STATISTICS = ("sum_counts", "mean_over_cells_of_log1p_cp10k")
DEFAULT_PSEUDOBULK_STATISTIC = "sum_counts"


def available_references() -> Dict[str, str]:
    """Named references and their bundle paths (present ones only are usable)."""
    return {name: str(path) for name, path in REFERENCE_BUNDLES.items()}


def resolve_bundle_path(
    bundle_path: Path | str | None = None, reference: Optional[str] = None
) -> Path:
    """Resolve an explicit bundle path or a named reference to a file path.

    Passing both is an error rather than a silent precedence rule.
    """
    if bundle_path is not None and reference is not None:
        raise ValueError("pass either bundle_path or reference, not both")
    if bundle_path is not None:
        return Path(bundle_path)
    name = DEFAULT_REFERENCE if reference is None else str(reference)
    if name not in REFERENCE_BUNDLES:
        raise KeyError(
            f"unknown GRN reference {name!r}; available: {sorted(REFERENCE_BUNDLES)}"
        )
    path = REFERENCE_BUNDLES[name]
    if not path.exists():
        raise FileNotFoundError(
            f"reference {name!r} is registered but its bundle is missing: {path}"
        )
    return path


def _clean_labels(values) -> List[str]:
    return [str(v).strip() for v in values]


@dataclass(frozen=True)
class PredictionResult:
    predictions: pd.DataFrame          # (n_groups, n_edges)
    summary: Dict[str, object]


def cp10k_log1p(counts: np.ndarray) -> np.ndarray:
    """CP10k + log1p normalization of a (rows x genes) count matrix."""
    counts = np.asarray(counts, dtype=np.float64)
    tot = counts.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1.0
    return np.log1p(counts / tot * 1e4).astype(np.float32)


def _open_bundle(path: Path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":   # gzip
        return gzip.open(path, "rb")
    return open(path, "rb")


class Rna2GrnBundle:
    def __init__(
        self,
        *,
        bundle_path: Path,
        model,
        input_genes: Sequence[str],
        output_edges: Sequence[str],
        metadata: Optional[Dict[str, object]] = None,
    ) -> None:
        self.bundle_path = Path(bundle_path)
        self.model = model
        self.input_genes = tuple(_clean_labels(input_genes))
        self.output_edges = tuple(_clean_labels(output_edges))
        self.metadata = metadata or {}
        self._gene_to_idx = {g: i for i, g in enumerate(self.input_genes)}

    # ------------------------------------------------------------------ load
    @classmethod
    def load(
        cls,
        bundle_path: Path | str | None = None,
        *,
        reference: Optional[str] = None,
    ) -> "Rna2GrnBundle":
        """Load a bundle by path, or by named reference (``leukemia``/``lung``).

        With neither argument the default reference (``leukemia``) is loaded, so
        existing callers keep their behaviour.
        """
        bundle_path = resolve_bundle_path(bundle_path, reference)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with _open_bundle(bundle_path) as handle:
                bundle = pickle.load(handle)
        required = {"model", "X_columns", "Y_columns"}
        missing = required.difference(bundle)
        if missing:
            raise ValueError(f"GRN bundle is missing required keys: {sorted(missing)}")
        return cls(
            bundle_path=bundle_path,
            model=bundle["model"],
            input_genes=bundle["X_columns"],
            output_edges=bundle["Y_columns"],
            metadata=bundle.get("metadata"),
        )

    @property
    def pseudobulk_statistic(self) -> str:
        """The statistic this bundle's training pseudobulks were built with."""
        value = str(self.metadata.get("pseudobulk_statistic") or "").strip()
        return value or DEFAULT_PSEUDOBULK_STATISTIC

    def model_info(self) -> Dict[str, object]:
        info = {
            "bundle_path": str(self.bundle_path),
            "n_input_genes": len(self.input_genes),
            "n_output_edges": len(self.output_edges),
            "model_class": type(self.model).__name__,
        }
        for key in ("reference", "tissue", "approach", "features", "ridge_lambda",
                    "normalization", "pseudobulk_statistic", "input_already_normalized",
                    "internal_validation", "n_train_pseudobulks", "n_edges", "n_tfs",
                    "k", "metric", "include_controls", "grn_source", "rna_source"):
            if key in self.metadata:
                info[key] = self.metadata[key]
        return info

    # --------------------------------------------------------------- helpers
    def _align_matrix(self, matrix: np.ndarray, gene_labels: Sequence[str]) -> Tuple[np.ndarray, int]:
        """Align a (rows x genes_in) matrix to the model feature order.

        Duplicate input genes are summed (counts); missing model genes are 0.
        """
        from collections import defaultdict
        pos_by_gene: Dict[str, List[int]] = defaultdict(list)
        for i, g in enumerate(_clean_labels(gene_labels)):
            if g:
                pos_by_gene[g].append(i)
        aligned = np.zeros((matrix.shape[0], len(self.input_genes)), dtype=np.float64)
        matched = 0
        M = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
        # Library size over ALL input genes — the CP10k denominator used at training
        # (and by cellHarmony). When the input is already cellHarmony-normalized
        # expression, use predict_*(normalized=True) and this path is skipped.
        totals = np.asarray(M.sum(axis=1)).reshape(-1).astype(np.float64)
        for g, model_idx in self._gene_to_idx.items():
            src = pos_by_gene.get(g)
            if not src:
                continue
            matched += 1
            aligned[:, model_idx] = M[:, src].sum(axis=1) if len(src) > 1 else M[:, src[0]]
        return aligned, matched, totals

    def _warn_if_counts_path(self) -> None:
        """Warn when a bundle trained on already-normalized pseudobulks is fed counts.

        The lung reference was trained on the MEAN over cells of log1p(CP10k); the
        counts path here computes log1p(CP10k of the SUMMED counts). Both are
        CP10k+log1p, but they are different statistics, so the caller is told
        rather than silently served a shifted input space.
        """
        if self.metadata.get("input_already_normalized"):
            warnings.warn(
                f"bundle {self.metadata.get('reference', '?')!r} was trained on "
                f"{self.metadata.get('pseudobulk_statistic', 'normalized')} pseudobulks; "
                "you passed counts, which are normalized here as log1p(CP10k of summed "
                "counts). Pass an already-normalized matrix with normalized=True to "
                "match the training statistic.",
                RuntimeWarning, stacklevel=3,
            )

    def _predict_aligned_counts(self, aligned_counts, totals, *, return_neighbors=False):
        self._warn_if_counts_path()
        tot = np.asarray(totals, dtype=np.float64).reshape(-1, 1)
        tot[tot == 0] = 1.0
        Xn = np.log1p(aligned_counts / tot * 1e4).astype(np.float32)
        return self.model.predict(Xn, return_neighbors=return_neighbors)

    def _neighbor_frame(self, idx: np.ndarray, w: np.ndarray, index) -> pd.DataFrame:
        rs = self.metadata.get("reference_sample", [])
        rc = self.metadata.get("reference_cell_state", [])
        rows = []
        for r, (nbrs, wts) in enumerate(zip(idx, w)):
            order = np.argsort(-wts)
            tops = [f"{rs[j] if rs else j}|{rc[j] if rc else ''}({wts[o]:.2f})"
                    for o, j in zip(order, np.asarray(nbrs)[order])]
            rows.append("; ".join(tops[:5]))
        return pd.DataFrame({"top_neighbors": rows}, index=index)

    # ------------------------------------------------------------- dataframe
    def predict_from_dataframe(
        self, expression: pd.DataFrame, *, normalized: bool = False,
        return_neighbors: bool = False,
    ) -> PredictionResult:
        """Predict GRNs from a pseudobulk-by-gene matrix.

        expression : rows = pseudobulks (sample x cell-state), columns = genes.
        normalized : set True if the matrix is already CP10k+log1p (then the
                     normalization step is skipped; alignment still applies).
        """
        if expression.empty:
            raise ValueError("Expression dataframe is empty")
        df = expression.copy()
        df.index = _clean_labels(df.index)
        df.columns = _clean_labels(df.columns)
        df = df.T.groupby(level=0).sum().T
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        aligned, matched, totals = self._align_matrix(df.to_numpy(), df.columns)
        if matched == 0:
            raise ValueError("No model feature genes were found in the dataframe")
        if normalized:
            pred = self.model.predict(aligned, return_neighbors=return_neighbors)
        else:
            pred = self._predict_aligned_counts(aligned, totals, return_neighbors=return_neighbors)
        return self._wrap(pred, df.index, matched, "dataframe", return_neighbors)

    # ----------------------------------------------------------------- adata
    def predict_from_adata(
        self, adata, *, groupby: Optional[str] = None, layer: Optional[str] = None,
        gene_symbol_col: Optional[str] = None, min_cells: int = 1,
        return_neighbors: bool = False, pseudobulk_statistic: Optional[str] = None,
    ) -> PredictionResult:
        """Predict GRNs from an AnnData of single cells (or pseudobulks).

        groupby : obs column whose groups become pseudobulks.
                  If None, each obs row is treated as an already-formed pseudobulk.
        pseudobulk_statistic :
            ``"sum_counts"`` (default for bundles that declare nothing) sums each
            group's counts, then applies CP10k + log1p — the leukemia convention.
            ``"mean_over_cells_of_log1p_cp10k"`` averages each group's ALREADY
            CP10k+log1p-normalized rows — the lung convention. When omitted, the
            bundle's own declared statistic is used, so a caller that passes nothing
            still gets the input space the model was trained on.
        """
        stat = str(pseudobulk_statistic or self.pseudobulk_statistic)
        if stat not in PSEUDOBULK_STATISTICS:
            raise ValueError(
                f"unknown pseudobulk_statistic {stat!r}; expected one of "
                f"{list(PSEUDOBULK_STATISTICS)}")
        normalized_input = stat == "mean_over_cells_of_log1p_cp10k"

        gene_labels = (adata.var[gene_symbol_col].astype(str).tolist()
                       if gene_symbol_col else _clean_labels(adata.var_names))
        X = adata.layers[layer] if layer else adata.X

        if normalized_input:
            # Fail loudly rather than average a count matrix as if it were log1p(CP10k).
            peak = float(X.max()) if hasattr(X, "max") else float(np.asarray(X).max())
            if peak > 20.0:
                raise ValueError(
                    f"pseudobulk_statistic={stat!r} needs a CP10k+log1p matrix, but the "
                    f"selected matrix peaks at {peak:.1f} (log1p(CP10k) cannot exceed "
                    f"{float(np.log1p(1e4)):.2f}). Pass layer=None so adata.X is used, or "
                    "pass pseudobulk_statistic='sum_counts' for a count matrix.")

        if groupby is None:
            grouped = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
            index = pd.Index(_clean_labels(adata.obs_names))
            kept = None
        else:
            if groupby not in adata.obs.columns:
                raise KeyError(f"obs column {groupby!r} not found")
            labels = adata.obs[groupby].astype(str).values
            groups, rows, kept = [], [], []
            for g in pd.unique(labels):
                m = labels == g
                if m.sum() < min_cells:
                    continue
                sub = X[m]
                sub = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
                rows.append(sub.mean(0) if normalized_input else sub.sum(0))
                groups.append(g); kept.append(int(m.sum()))
            if not rows:
                raise ValueError(f"no group in {groupby!r} reached min_cells={min_cells}")
            grouped = np.vstack(rows)
            index = pd.Index(groups, name=groupby)

        aligned, matched, totals = self._align_matrix(grouped, gene_labels)
        if matched == 0:
            raise ValueError("No model feature genes were found in the AnnData")
        if normalized_input:
            pred = self.model.predict(aligned, return_neighbors=return_neighbors)
        else:
            pred = self._predict_aligned_counts(aligned, totals, return_neighbors=return_neighbors)
        res = self._wrap(pred, index, matched, "adata", return_neighbors)
        res.summary["pseudobulk_statistic"] = stat
        if groupby is not None:
            res.summary["n_cells_per_group"] = dict(zip([str(x) for x in index], kept))
        return res

    def predict_from_h5ad(self, h5ad_path, **kw) -> PredictionResult:
        import anndata as ad
        return self.predict_from_adata(ad.read_h5ad(str(h5ad_path)), **kw)

    def predict_from_10x_h5(self, h5_path, **kw) -> PredictionResult:
        """Read a 10x ``filtered_feature_bc_matrix.h5`` and predict.

        With no ``groupby`` this returns one GRN for the whole-sample pseudobulk;
        pass ``groupby`` after annotating cells (e.g. cellHarmony cell states) for
        per-cell-state GRNs.
        """
        import scanpy as sc
        adata = sc.read_10x_h5(str(h5_path))
        adata.var_names_make_unique()
        # A raw 10x .h5 is cell-level. Without a grouping column, predict one GRN
        # for the whole-sample pseudobulk rather than one per cell.
        if kw.get("groupby") is None:
            adata.obs["_whole_sample"] = "whole_sample"
            kw["groupby"] = "_whole_sample"
        return self.predict_from_adata(adata, **kw)

    # ------------------------------------------------- TF activity / differential
    def tf_activity(self, predictions: pd.DataFrame) -> pd.DataFrame:
        """Aggregate per-edge predictions into a per-TF activity matrix
        (groups x TFs), the mean activity of each TF's outgoing edges."""
        edge_tf = np.asarray(self.metadata.get("edge_tf"))
        if edge_tf is None or len(edge_tf) != predictions.shape[1]:
            raise ValueError("bundle metadata lacks edge_tf aligned to outputs")
        vals = predictions.to_numpy()
        tfs = pd.unique(edge_tf)
        out = {tf: vals[:, edge_tf == tf].mean(1) for tf in tfs}
        return pd.DataFrame(out, index=predictions.index)

    def differential_tf_ranking(
        self, query_pred: pd.DataFrame, control_pred: pd.DataFrame,
    ) -> pd.DataFrame:
        """Rank TFs by increased regulatory activity in query vs control.

        Both inputs are per-edge prediction frames (e.g. one row per cell state)
        aligned on their index; returns TFs sorted by mean (query-control) edge
        activity — the TFs whose regulons are predicted MORE active in the query.
        This is the core 'identify the perturbed TF' output.
        """
        common = query_pred.index.intersection(control_pred.index)
        q = self.tf_activity(query_pred.loc[common])
        c = self.tf_activity(control_pred.loc[common])
        diff = (q - c).mean(axis=0).sort_values(ascending=False)
        return diff.rename("delta_activity").to_frame().reset_index().rename(
            columns={"index": "TF"})

    # ----------------------------------------------------------------- wrap
    def _wrap(self, pred, index, matched, kind, return_neighbors) -> PredictionResult:
        if return_neighbors and isinstance(pred, tuple):
            P, idx, w = pred
            nbr = self._neighbor_frame(idx, w, index) if idx is not None else None
        else:
            P, nbr = (pred[0] if isinstance(pred, tuple) else pred), None
        df = pd.DataFrame(np.asarray(P), index=index, columns=list(self.output_edges))
        summary = {
            "bundle_path": str(self.bundle_path),
            "input_kind": kind,
            "n_groups": int(df.shape[0]),
            "matched_genes": int(matched),
            "missing_genes": len(self.input_genes) - int(matched),
            "model_gene_count": len(self.input_genes),
            "n_output_edges": len(self.output_edges),
            "model_class": type(self.model).__name__,
        }
        if nbr is not None:
            summary["neighbors"] = nbr
        return PredictionResult(predictions=df, summary=summary)


def load_bundle(
    bundle_path: Path | str | None = None, *, reference: Optional[str] = None
) -> Rna2GrnBundle:
    """Load a GRN bundle by path or by named reference (default: ``leukemia``)."""
    return Rna2GrnBundle.load(bundle_path, reference=reference)
