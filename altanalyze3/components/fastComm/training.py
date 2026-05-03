from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ResponseTrainingParams:
    input: Path
    output: Path
    manifest: Optional[Path] = None
    id_col: str = "signature"
    gene_col: str = "gene"
    score_col: str = "score"
    input_format: str = "auto"
    top_genes: int = 200
    min_abs_score: float = 0.0
    l2_normalize: bool = True


def _read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python")


def _matrix_from_long(df: pd.DataFrame, *, id_col: str, gene_col: str, score_col: str) -> pd.DataFrame:
    missing = {id_col, gene_col, score_col}.difference(df.columns)
    if missing:
        raise ValueError(f"Long response table is missing columns: {sorted(missing)}")
    work = df[[id_col, gene_col, score_col]].copy()
    work[id_col] = work[id_col].astype(str).str.strip()
    work[gene_col] = work[gene_col].astype(str).str.strip()
    work[score_col] = pd.to_numeric(work[score_col], errors="coerce")
    work = work.dropna(subset=[score_col])
    work = work.loc[work[id_col].ne("") & work[gene_col].ne("")]
    if work.empty:
        raise ValueError("No valid response rows remain after cleaning")
    return work.pivot_table(index=id_col, columns=gene_col, values=score_col, aggfunc="mean", fill_value=0.0)


def _matrix_from_wide(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep=None, engine="python", index_col=0)
    matrix.index = matrix.index.astype(str).str.strip()
    matrix.columns = matrix.columns.astype(str).str.strip()
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    matrix = matrix.loc[matrix.index.to_series().ne("").to_numpy()]
    matrix = matrix.loc[:, matrix.columns.to_series().ne("").to_numpy()]
    return matrix


def load_response_training_input(params: ResponseTrainingParams) -> pd.DataFrame:
    input_format = params.input_format.strip().lower()
    if input_format not in {"auto", "long", "wide"}:
        raise ValueError("input_format must be one of: auto, long, wide")
    if input_format == "wide":
        return _matrix_from_wide(params.input)

    raw = _read_table(params.input)
    long_columns = {params.id_col, params.gene_col, params.score_col}
    if input_format == "long" or long_columns.issubset(raw.columns):
        return _matrix_from_long(
            raw,
            id_col=params.id_col,
            gene_col=params.gene_col,
            score_col=params.score_col,
        )
    return _matrix_from_wide(params.input)


def prune_response_matrix(
    matrix: pd.DataFrame,
    *,
    top_genes: int = 200,
    min_abs_score: float = 0.0,
    l2_normalize: bool = True,
) -> pd.DataFrame:
    if matrix.empty:
        raise ValueError("Response matrix is empty")
    pruned = matrix.copy().astype(float)
    if min_abs_score > 0:
        pruned = pruned.where(pruned.abs() >= min_abs_score, 0.0)
    if top_genes > 0:
        rows = []
        for _, row in pruned.iterrows():
            keep = row.abs().sort_values(ascending=False).head(top_genes).index
            filtered = pd.Series(0.0, index=pruned.columns)
            filtered.loc[keep] = row.loc[keep]
            rows.append(filtered)
        pruned = pd.DataFrame(rows, index=pruned.index, columns=pruned.columns)
    if l2_normalize:
        norms = np.sqrt((pruned * pruned).sum(axis=1))
        pruned = pruned.div(norms.replace(0.0, 1.0), axis=0)
    nonzero_cols = pruned.columns[(pruned != 0).any(axis=0)]
    return pruned.loc[:, nonzero_cols]


def build_response_matrix(params: ResponseTrainingParams) -> Dict[str, object]:
    matrix = load_response_training_input(params)
    pruned = prune_response_matrix(
        matrix,
        top_genes=params.top_genes,
        min_abs_score=params.min_abs_score,
        l2_normalize=params.l2_normalize,
    )
    params.output.parent.mkdir(parents=True, exist_ok=True)
    pruned.to_csv(params.output, sep="\t")
    manifest = {
        "input": str(params.input),
        "output": str(params.output),
        "input_shape": [int(matrix.shape[0]), int(matrix.shape[1])],
        "output_shape": [int(pruned.shape[0]), int(pruned.shape[1])],
        "top_genes": int(params.top_genes),
        "min_abs_score": float(params.min_abs_score),
        "l2_normalize": bool(params.l2_normalize),
    }
    manifest_path = params.manifest
    if manifest_path is not None:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest
