from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

try:
    from sklearn.exceptions import InconsistentVersionWarning
except ImportError:  # pragma: no cover
    InconsistentVersionWarning = Warning

import pickle


DEFAULT_BUNDLE_PATH = Path(__file__).with_name("otherelastic_multitask_try.pkl")


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


class Rna2LipidBundle:
    def __init__(
        self,
        *,
        bundle_path: Path,
        model,
        scaler_x,
        scaler_y,
        input_genes: Sequence[str],
        output_lipids: Sequence[str],
    ) -> None:
        self.bundle_path = Path(bundle_path)
        self.model = model
        self.scaler_x = scaler_x
        self.scaler_y = scaler_y
        self.input_genes = tuple(_clean_labels(input_genes))
        self.output_lipids = tuple(_clean_labels(output_lipids))
        self._input_gene_to_index = {gene: idx for idx, gene in enumerate(self.input_genes)}

    @classmethod
    def load(cls, bundle_path: Path | str = DEFAULT_BUNDLE_PATH) -> "Rna2LipidBundle":
        bundle_path = Path(bundle_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", InconsistentVersionWarning)
            with bundle_path.open("rb") as handle:
                bundle = pickle.load(handle)

        required = {"model", "scaler_x", "scaler_y", "X_columns", "Y_columns"}
        missing = required.difference(bundle)
        if missing:
            raise ValueError(f"Model bundle is missing required keys: {sorted(missing)}")

        return cls(
            bundle_path=bundle_path,
            model=bundle["model"],
            scaler_x=bundle["scaler_x"],
            scaler_y=bundle["scaler_y"],
            input_genes=bundle["X_columns"],
            output_lipids=bundle["Y_columns"],
        )

    def model_info(self) -> Dict[str, object]:
        return {
            "bundle_path": str(self.bundle_path),
            "n_input_genes": len(self.input_genes),
            "n_output_lipids": len(self.output_lipids),
            "model_class": type(self.model).__name__,
            "scaler_x_class": type(self.scaler_x).__name__,
            "scaler_y_class": type(self.scaler_y).__name__,
        }

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
        chunk_size: int = 2048,
    ) -> PredictionResult:
        h5ad_path = Path(h5ad_path)
        adata = ad.read_h5ad(h5ad_path)

        obs_names = _clean_labels(adata.obs_names)
        if len(set(obs_names)) != len(obs_names):
            obs_names = _make_unique(obs_names)
        adata.obs_names = pd.Index(obs_names)

        if groupby is not None and groupby not in adata.obs.columns:
            raise KeyError(f"Column {groupby!r} was not found in adata.obs")

        gene_labels = self._extract_gene_labels(adata, gene_symbol_col=gene_symbol_col)
        matched_positions = self._matched_model_positions(gene_labels)
        if not matched_positions:
            raise ValueError("No model genes were found in the provided h5ad")
        predictions: List[pd.DataFrame] = []

        for start in range(0, adata.n_obs, chunk_size):
            stop = min(start + chunk_size, adata.n_obs)
            matrix = self._read_matrix_chunk(adata, start=start, stop=stop, layer=layer)
            aligned = self._align_matrix_chunk(matrix, matched_positions)
            predicted = self._predict_aligned_matrix(aligned)
            predicted.index = adata.obs_names[start:stop]
            predictions.append(predicted)

        prediction_df = pd.concat(predictions, axis=0)
        summary = self._build_summary(
            input_rows=int(adata.n_obs),
            matched_genes=len(matched_positions),
            input_kind="h5ad",
        )
        summary["h5ad_path"] = str(h5ad_path)
        summary["layer"] = layer or "X"
        summary["gene_symbol_source"] = gene_symbol_col or "var_names"

        if groupby is not None:
            groups = adata.obs[groupby].astype(str).str.strip()
            prediction_df = prediction_df.groupby(groups, dropna=False).mean()
            prediction_df.index.name = groupby
            summary["groupby"] = groupby
            summary["output_rows"] = int(prediction_df.shape[0])
        else:
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
        aligned = np.zeros((n_rows, len(self.input_genes)), dtype=float)

        for model_idx, source_positions in matched_positions:
            chunk = matrix[:, source_positions]
            if sparse.issparse(chunk):
                chunk = chunk.toarray()
            else:
                chunk = np.asarray(chunk)

            if chunk.ndim == 1:
                aligned[:, model_idx] = chunk.astype(float, copy=False)
            else:
                aligned[:, model_idx] = np.asarray(chunk, dtype=float).mean(axis=1)

        return aligned

    def _predict_aligned_matrix(self, aligned_matrix) -> pd.DataFrame:
        if not isinstance(aligned_matrix, pd.DataFrame):
            aligned_matrix = pd.DataFrame(aligned_matrix, columns=self.input_genes)
        transformed = self.scaler_x.transform(aligned_matrix)
        predicted = self.model.predict(transformed)
        predicted = self.scaler_y.inverse_transform(predicted)
        return pd.DataFrame(predicted, columns=self.output_lipids)

    def _build_summary(self, *, input_rows: int, matched_genes: int, input_kind: str) -> Dict[str, object]:
        return {
            "bundle_path": str(self.bundle_path),
            "input_kind": input_kind,
            "input_rows": input_rows,
            "matched_genes": matched_genes,
            "missing_genes": len(self.input_genes) - matched_genes,
            "model_gene_count": len(self.input_genes),
            "output_lipid_count": len(self.output_lipids),
        }


def load_bundle(bundle_path: Path | str = DEFAULT_BUNDLE_PATH) -> Rna2LipidBundle:
    return Rna2LipidBundle.load(bundle_path=bundle_path)
