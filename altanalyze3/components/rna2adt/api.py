"""Inference API for predicting CITE-seq ADT (antibody-derived tag) values from RNA.

This module mirrors the structure of ``components.rna2lipid.api`` so the
cellHarmony-web pipeline can register an ADT imputation modality alongside
lipids using the same plumbing. A bundle is a pickle with the keys::

    model      : object with ``.predict(X)`` returning (n_cells, n_adts)
    scaler_x   : sklearn-style transformer with ``.transform`` (or None)
    scaler_y   : sklearn-style transformer used to inverse-transform predictions
                 if the model was trained on standardized targets (or None)
    X_columns  : sequence of RNA gene symbols expected as input
    Y_columns  : sequence of ADT names (e.g. ``"Hu.CD4"``) produced as output
    metadata   : optional dict of training metadata

The model object can be any class exposing ``.predict``; for the kNN-in-PCA
baseline see ``training.KnnAdtRegressor``.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import pickle
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

try:
    from sklearn.exceptions import InconsistentVersionWarning
    from sklearn.utils.validation import check_is_fitted
except ImportError:  # pragma: no cover
    InconsistentVersionWarning = Warning

    def check_is_fitted(estimator, attributes=None) -> None:
        if not hasattr(estimator, "n_features_in_"):
            raise AttributeError("Estimator does not appear to be fitted")


DEFAULT_BUNDLE_PATH = Path(__file__).with_name("rna2adt_bm_bundle.pkl")


def _clean_labels(values: Iterable[object]) -> List[str]:
    return [str(value).strip() for value in values]


def _make_unique(values: Sequence[str]) -> List[str]:
    counts: Dict[str, int] = {}
    unique: List[str] = []
    for value in values:
        if value not in counts:
            counts[value] = 1
            unique.append(value)
        else:
            suffix = counts[value]
            unique.append(f"{value}.{suffix}")
            counts[value] += 1
    return unique


@dataclass(frozen=True)
class PredictionResult:
    predictions: pd.DataFrame
    summary: Dict[str, object]


def _transformer_is_fitted(transformer) -> bool:
    if transformer is None:
        return False
    fitted_state_attrs = ("scale_", "mean_", "var_", "n_samples_seen_")
    if any(hasattr(transformer, attr) for attr in fitted_state_attrs):
        return True
    try:
        check_is_fitted(transformer, attributes=fitted_state_attrs)
    except Exception:
        return False
    return True


def _model_output_scaling_mode(metadata: Optional[Dict[str, object]], scaler_y) -> str:
    scaling_info = (metadata or {}).get("target_scaling")
    if isinstance(scaling_info, dict):
        mode = str(scaling_info.get("mode", "")).strip().lower()
        if mode in {"none", "standard"}:
            return mode
    return "standard" if _transformer_is_fitted(scaler_y) else "none"


class Rna2AdtBundle:
    def __init__(
        self,
        *,
        bundle_path: Path,
        model,
        scaler_x,
        scaler_y,
        input_genes: Sequence[str],
        output_adts: Sequence[str],
        metadata: Optional[Dict[str, object]] = None,
    ) -> None:
        self.bundle_path = Path(bundle_path)
        self.model = model
        self.scaler_x = scaler_x
        self.scaler_y = scaler_y
        self.input_genes = tuple(_clean_labels(input_genes))
        self.output_adts = tuple(_clean_labels(output_adts))
        self.metadata = metadata or {}
        self.target_scaling_mode = _model_output_scaling_mode(self.metadata, scaler_y)
        self._input_gene_to_index = {gene: idx for idx, gene in enumerate(self.input_genes)}

    @classmethod
    def load(cls, bundle_path: Path | str = DEFAULT_BUNDLE_PATH) -> "Rna2AdtBundle":
        bundle_path = Path(bundle_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", InconsistentVersionWarning)
            with bundle_path.open("rb") as handle:
                bundle = pickle.load(handle)

        required = {"model", "X_columns", "Y_columns"}
        missing = required.difference(bundle)
        if missing:
            raise ValueError(f"Model bundle is missing required keys: {sorted(missing)}")

        return cls(
            bundle_path=bundle_path,
            model=bundle["model"],
            scaler_x=bundle.get("scaler_x"),
            scaler_y=bundle.get("scaler_y"),
            input_genes=bundle["X_columns"],
            output_adts=bundle["Y_columns"],
            metadata=bundle.get("metadata"),
        )

    def model_info(self) -> Dict[str, object]:
        info = {
            "bundle_path": str(self.bundle_path),
            "n_input_genes": len(self.input_genes),
            "n_output_adts": len(self.output_adts),
            "model_class": type(self.model).__name__,
            "scaler_x_class": type(self.scaler_x).__name__ if self.scaler_x is not None else None,
            "scaler_y_class": type(self.scaler_y).__name__ if self.scaler_y is not None else None,
            "target_scaling_mode": self.target_scaling_mode,
        }
        if self.metadata:
            info["metadata"] = self.metadata
        return info

    def predict_from_dataframe(self, expression: pd.DataFrame) -> PredictionResult:
        aligned, matched_genes = self._align_dataframe(expression)
        predicted = self._predict_aligned_matrix(aligned)
        summary = self._build_summary(
            input_rows=int(aligned.shape[0]),
            matched_genes=matched_genes,
            input_kind="dataframe",
        )
        predicted.index = aligned.index
        summary["output_rows"] = int(predicted.shape[0])
        return PredictionResult(predictions=predicted, summary=summary)

    def predict_from_h5ad(
        self,
        h5ad_path: Path | str,
        *,
        layer: Optional[str] = None,
        gene_symbol_col: Optional[str] = None,
        groupby: Optional[str] = None,
        chunk_size: int = 4096,
    ) -> PredictionResult:
        h5ad_path = Path(h5ad_path)
        adata = ad.read_h5ad(h5ad_path)
        return self.predict_from_adata(
            adata,
            layer=layer,
            gene_symbol_col=gene_symbol_col,
            groupby=groupby,
            chunk_size=chunk_size,
        )

    def predict_from_adata(
        self,
        adata: ad.AnnData,
        *,
        layer: Optional[str] = None,
        gene_symbol_col: Optional[str] = None,
        groupby: Optional[str] = None,
        chunk_size: int = 4096,
        cluster_obs_col: Optional[str] = None,
    ) -> PredictionResult:
        obs_names = _clean_labels(adata.obs_names)
        if len(set(obs_names)) != len(obs_names):
            obs_names = _make_unique(obs_names)
        # avoid copying the full RNA matrix; only patch obs index
        adata_local = adata
        if list(adata_local.obs_names.astype(str)) != obs_names:
            adata_local = adata.copy()
            adata_local.obs_names = pd.Index(obs_names)

        if groupby is not None and groupby not in adata_local.obs.columns:
            raise KeyError(f"Column {groupby!r} was not found in adata.obs")

        gene_labels = self._extract_gene_labels(adata_local, gene_symbol_col=gene_symbol_col)
        matched_positions = self._matched_model_positions(gene_labels)
        if not matched_positions:
            raise ValueError("No model genes were found in the provided AnnData")

        cluster_labels: Optional[np.ndarray] = None
        if cluster_obs_col is not None:
            if cluster_obs_col not in adata_local.obs.columns:
                raise KeyError(f"obs column {cluster_obs_col!r} not found")
            cluster_labels = adata_local.obs[cluster_obs_col].astype(str).to_numpy()

        predictions: List[pd.DataFrame] = []
        for start in range(0, adata_local.n_obs, chunk_size):
            stop = min(start + chunk_size, adata_local.n_obs)
            matrix = self._read_matrix_chunk(adata_local, start=start, stop=stop, layer=layer)
            aligned = self._align_matrix_chunk(matrix, matched_positions)
            chunk_clusters = cluster_labels[start:stop] if cluster_labels is not None else None
            predicted = self._predict_aligned_matrix(aligned, cluster_labels=chunk_clusters)
            predicted.index = pd.Index(obs_names[start:stop])
            predictions.append(predicted)

        prediction_df = pd.concat(predictions, axis=0)
        summary = self._build_summary(
            input_rows=int(adata_local.n_obs),
            matched_genes=len(matched_positions),
            input_kind="adata",
        )
        summary["layer"] = layer or "X"
        summary["gene_symbol_source"] = gene_symbol_col or "var_names"

        if groupby is not None:
            groups = adata_local.obs[groupby].astype(str).str.strip()
            prediction_df = prediction_df.groupby(groups, dropna=False).mean()
            prediction_df.index.name = groupby

        summary["output_rows"] = int(prediction_df.shape[0])
        return PredictionResult(predictions=prediction_df, summary=summary)

    def _align_dataframe(self, expression: pd.DataFrame) -> Tuple[pd.DataFrame, int]:
        if expression.empty:
            raise ValueError("Expression dataframe is empty")
        df = expression.copy()
        df.index = _clean_labels(df.index)
        df.columns = _clean_labels(df.columns)
        df = df.T.groupby(level=0).mean().T
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        matched_genes = int(df.columns.isin(self.input_genes).sum())
        if matched_genes == 0:
            raise ValueError("No model genes were found in the provided dataframe")
        aligned = df.reindex(columns=self.input_genes, fill_value=0.0)
        return aligned, matched_genes

    def _extract_gene_labels(self, adata: ad.AnnData, *, gene_symbol_col: Optional[str]) -> List[str]:
        if gene_symbol_col is None:
            return _clean_labels(adata.var_names)
        if gene_symbol_col not in adata.var.columns:
            raise KeyError(f"Column {gene_symbol_col!r} was not found in adata.var")
        return _clean_labels(adata.var[gene_symbol_col].tolist())

    def _matched_model_positions(self, gene_labels: Sequence[str]) -> List[Tuple[int, List[int]]]:
        positions_by_gene: Dict[str, List[int]] = defaultdict(list)
        for idx, gene in enumerate(gene_labels):
            if gene:
                positions_by_gene[gene].append(idx)
        matched: List[Tuple[int, List[int]]] = []
        for gene, source_positions in positions_by_gene.items():
            model_idx = self._input_gene_to_index.get(gene)
            if model_idx is not None:
                matched.append((model_idx, source_positions))
        return matched

    def _read_matrix_chunk(
        self,
        adata: ad.AnnData,
        *,
        start: int,
        stop: int,
        layer: Optional[str],
    ):
        if layer is None:
            return adata.X[start:stop]
        if layer not in adata.layers:
            raise KeyError(f"Layer {layer!r} was not found in adata.layers")
        return adata.layers[layer][start:stop]

    def _align_matrix_chunk(
        self,
        matrix,
        matched_positions: Sequence[Tuple[int, Sequence[int]]],
    ) -> np.ndarray:
        n_rows = matrix.shape[0]
        aligned = np.zeros((n_rows, len(self.input_genes)), dtype=np.float32)
        single_model_idx: List[int] = []
        single_source_idx: List[int] = []
        multi_source_positions: List[Tuple[int, Sequence[int]]] = []
        for model_idx, source_positions in matched_positions:
            if len(source_positions) == 1:
                single_model_idx.append(model_idx)
                single_source_idx.append(int(source_positions[0]))
            else:
                multi_source_positions.append((model_idx, source_positions))

        if single_source_idx:
            chunk = matrix[:, single_source_idx]
            if sparse.issparse(chunk):
                chunk = chunk.toarray()
            else:
                chunk = np.asarray(chunk)
            aligned[:, single_model_idx] = np.asarray(chunk, dtype=np.float32)

        for model_idx, source_positions in multi_source_positions:
            chunk = matrix[:, source_positions]
            if sparse.issparse(chunk):
                chunk = chunk.toarray()
            else:
                chunk = np.asarray(chunk)
            if chunk.ndim == 1:
                aligned[:, model_idx] = chunk.astype(np.float32, copy=False)
            else:
                aligned[:, model_idx] = np.asarray(chunk, dtype=np.float32).mean(axis=1)
        return aligned

    def _predict_aligned_matrix(self, aligned_matrix, *, cluster_labels=None) -> pd.DataFrame:
        if isinstance(aligned_matrix, pd.DataFrame):
            aligned_values = aligned_matrix.to_numpy(dtype=np.float32, copy=False)
        else:
            aligned_values = np.asarray(aligned_matrix, dtype=np.float32)
        if self.scaler_x is not None and _transformer_is_fitted(self.scaler_x):
            transformed = self.scaler_x.transform(aligned_values)
        else:
            transformed = aligned_values
        # Centroid models accept a kwarg; everything else uses the standard
        # one-arg .predict.
        try:
            predicted = self.model.predict(transformed, cluster_labels=cluster_labels)
        except TypeError:
            predicted = self.model.predict(transformed)
        if self.target_scaling_mode == "standard":
            if self.scaler_y is None or not _transformer_is_fitted(self.scaler_y):
                raise ValueError(
                    f"Bundle {self.bundle_path} declares standardized targets but scaler_y is unavailable."
                )
            predicted = self.scaler_y.inverse_transform(predicted)
        else:
            predicted = np.asarray(predicted, dtype=float)
        return pd.DataFrame(predicted, columns=list(self.output_adts))

    def _build_summary(self, *, input_rows: int, matched_genes: int, input_kind: str) -> Dict[str, object]:
        return {
            "bundle_path": str(self.bundle_path),
            "input_kind": input_kind,
            "input_rows": input_rows,
            "matched_genes": matched_genes,
            "missing_genes": len(self.input_genes) - matched_genes,
            "model_gene_count": len(self.input_genes),
            "output_adt_count": len(self.output_adts),
            "target_scaling_mode": self.target_scaling_mode,
            "model_class": type(self.model).__name__,
        }


def load_bundle(bundle_path: Path | str = DEFAULT_BUNDLE_PATH) -> Rna2AdtBundle:
    return Rna2AdtBundle.load(bundle_path=bundle_path)
