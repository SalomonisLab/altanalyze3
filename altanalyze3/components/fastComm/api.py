from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd

from .scoring import (
    add_receiver_response_scores,
    clean_labels,
    finalize_scores,
    load_response_matrix,
    make_receiver_delta,
    make_state_pseudobulk,
    prepare_ligand_receptor_table,
    score_ligand_receptor_expression,
    split_complex,
    to_dense_frame,
)


PACKAGE_DIR = Path(__file__).resolve().parent
DEFAULT_LR_TABLE = PACKAGE_DIR / "resources" / "seed_ligand_receptor.tsv"
DEFAULT_RESPONSE_MATRIX = PACKAGE_DIR / "resources" / "seed_response_signatures.tsv"


@dataclass(frozen=True)
class FastCommParams:
    expression: Optional[Path] = None
    metadata: Optional[Path] = None
    h5ad: Optional[Path] = None
    output: Path = Path("fastcomm_scores.tsv")
    lr_table: Path = DEFAULT_LR_TABLE
    response_matrix: Optional[Path] = DEFAULT_RESPONSE_MATRIX
    state_pair_output: Optional[Path] = None
    state_expression_output: Optional[Path] = None
    state_key: str = "cell_state"
    layer: Optional[str] = None
    gene_symbol_col: Optional[str] = None
    min_cells: int = 1
    min_ligand_expr: float = 0.01
    min_receptor_expr: float = 0.01
    min_lr_expression_score: float = 0.001
    include_self_edges: bool = True
    baseline_state: Optional[str] = None
    top_n: Optional[int] = None


@dataclass(frozen=True)
class FastCommResult:
    scores: pd.DataFrame
    state_expression: pd.DataFrame
    state_sizes: pd.Series
    summary: Dict[str, object]


def _read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python", index_col=0)


def _materialize_matrix(matrix) -> np.ndarray:
    if hasattr(matrix, "toarray"):
        return matrix.toarray()
    return np.asarray(matrix)


def _deduplicate_columns(df: pd.DataFrame) -> pd.DataFrame:
    if not pd.Index(df.columns).duplicated().any():
        return df
    return df.T.groupby(level=0, sort=False).mean().T


def _gene_match_key(value: object) -> str:
    return str(value or "").strip().upper()


def _required_gene_map(required_genes: Optional[Sequence[str]]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for gene in clean_labels(required_genes or []):
        if gene:
            mapping.setdefault(_gene_match_key(gene), gene)
    return mapping


def _select_gene_positions(
    genes: Sequence[str],
    required_genes: Optional[Sequence[str]],
) -> tuple[list[int], list[str], Dict[str, object]]:
    if not required_genes:
        return list(range(len(genes))), clean_labels(genes), {
            "n_required_genes": 0,
            "n_matched_required_genes": len(set(clean_labels(genes))),
            "n_missing_required_genes": 0,
            "missing_required_genes": [],
        }

    required_by_key = _required_gene_map(required_genes)
    matched_keys: set[str] = set()
    positions: list[int] = []
    canonical_columns: list[str] = []
    for idx, gene in enumerate(genes):
        key = _gene_match_key(gene)
        canonical = required_by_key.get(key)
        if canonical is None:
            continue
        positions.append(idx)
        canonical_columns.append(canonical)
        matched_keys.add(key)

    missing = [
        gene
        for key, gene in required_by_key.items()
        if key not in matched_keys
    ]
    diagnostics = {
        "n_required_genes": int(len(required_by_key)),
        "n_matched_required_genes": int(len(matched_keys)),
        "n_missing_required_genes": int(len(missing)),
        "missing_required_genes": missing[:50],
    }
    return positions, canonical_columns, diagnostics


def _expression_from_h5ad(
    params: FastCommParams,
    *,
    required_genes: Optional[Sequence[str]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    if params.h5ad is None:
        raise ValueError("h5ad path is required")
    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover - exercised only in minimal envs
        raise ImportError("AnnData input requires the 'anndata' package") from exc

    adata = ad.read_h5ad(params.h5ad, backed="r")
    if params.gene_symbol_col:
        if params.gene_symbol_col not in adata.var.columns:
            raise KeyError(f"Column {params.gene_symbol_col!r} was not found in adata.var")
        genes = adata.var[params.gene_symbol_col].astype(str).tolist()
    else:
        genes = adata.var_names.astype(str).tolist()

    gene_positions, selected_genes, diagnostics = _select_gene_positions(genes, required_genes)
    if required_genes and not gene_positions:
        raise ValueError("None of the requested ligand/receptor/response genes were found in the h5ad")

    matrix_source = adata.layers[params.layer] if params.layer else adata.X
    matrix = _materialize_matrix(matrix_source[:, gene_positions])
    expression = to_dense_frame(matrix, index=adata.obs_names.astype(str).tolist(), columns=selected_genes)
    expression = _deduplicate_columns(expression)
    metadata = adata.obs.copy()
    adata.file.close()
    diagnostics["n_loaded_columns"] = int(expression.shape[1])
    return expression, metadata, diagnostics


def _load_inputs(
    params: FastCommParams,
    *,
    required_genes: Optional[Sequence[str]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    if params.h5ad is not None:
        return _expression_from_h5ad(params, required_genes=required_genes)
    if params.expression is None or params.metadata is None:
        raise ValueError("Provide either --h5ad or both --expression and --metadata")
    expression = _read_table(params.expression)
    original_columns = clean_labels(expression.columns)
    expression.columns = original_columns
    diagnostics: Dict[str, object] = {
        "n_required_genes": 0,
        "n_matched_required_genes": int(len(set(original_columns))),
        "n_missing_required_genes": 0,
        "missing_required_genes": [],
    }
    if required_genes:
        required_by_key = _required_gene_map(required_genes)
        rename: Dict[str, str] = {}
        keep: list[str] = []
        matched_keys: set[str] = set()
        for column in original_columns:
            canonical = required_by_key.get(_gene_match_key(column))
            if canonical is None:
                continue
            keep.append(column)
            rename[column] = canonical
            matched_keys.add(_gene_match_key(column))
        if not keep:
            raise ValueError("None of the requested ligand/receptor/response genes were found in the expression matrix")
        missing = [gene for key, gene in required_by_key.items() if key not in matched_keys]
        expression = expression.loc[:, pd.Index(keep).drop_duplicates()].rename(columns=rename).copy()
        expression = _deduplicate_columns(expression)
        diagnostics = {
            "n_required_genes": int(len(required_by_key)),
            "n_matched_required_genes": int(len(matched_keys)),
            "n_missing_required_genes": int(len(missing)),
            "missing_required_genes": missing[:50],
        }
    metadata = _read_table(params.metadata)
    diagnostics["n_loaded_columns"] = int(expression.shape[1])
    return expression, metadata, diagnostics


def _required_genes(lr_table: pd.DataFrame, response_matrix: Optional[pd.DataFrame]) -> list[str]:
    genes: list[str] = []
    lr = prepare_ligand_receptor_table(lr_table)
    for edge in lr.itertuples(index=False):
        genes.extend(split_complex(edge.ligand))
        genes.extend(split_complex(edge.receptor))
    if response_matrix is not None:
        genes.extend(response_matrix.columns.astype(str).tolist())
    seen = set()
    ordered: list[str] = []
    for gene in genes:
        gene = str(gene).strip()
        if gene and gene not in seen:
            seen.add(gene)
            ordered.append(gene)
    return ordered


def summarize_state_pairs(scores: pd.DataFrame) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame(
            columns=[
                "sender_state",
                "receiver_state",
                "n_edges",
                "max_fastcomm_score",
                "sum_fastcomm_score",
                "mean_receiver_response_score",
                "top_interaction",
            ]
        )

    rows = []
    for (sender, receiver), group in scores.groupby(["sender_state", "receiver_state"], sort=False):
        top = group.sort_values("fastcomm_score", ascending=False).iloc[0]
        rows.append(
            {
                "sender_state": sender,
                "receiver_state": receiver,
                "n_edges": int(group.shape[0]),
                "max_fastcomm_score": float(group["fastcomm_score"].max()),
                "sum_fastcomm_score": float(group["fastcomm_score"].sum()),
                "mean_receiver_response_score": float(group["receiver_response_score"].mean()),
                "top_interaction": f"{top['ligand']}->{top['receptor']}",
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["max_fastcomm_score", "sum_fastcomm_score"],
        ascending=[False, False],
    )


def run_fastcomm(params: FastCommParams) -> FastCommResult:
    lr_table = pd.read_csv(params.lr_table, sep=None, engine="python")
    response_matrix = load_response_matrix(str(params.response_matrix)) if params.response_matrix else None
    required_genes = _required_genes(lr_table, response_matrix)
    expression, metadata, gene_diagnostics = _load_inputs(params, required_genes=required_genes)
    state = make_state_pseudobulk(
        expression,
        metadata,
        state_key=params.state_key,
        min_cells=params.min_cells,
    )
    edges = score_ligand_receptor_expression(
        state,
        lr_table,
        min_ligand_expr=params.min_ligand_expr,
        min_receptor_expr=params.min_receptor_expr,
        include_self_edges=params.include_self_edges,
    )
    if params.min_lr_expression_score > 0 and not edges.empty:
        edges = edges.loc[edges["lr_expression_score"] >= params.min_lr_expression_score].copy()

    receiver_delta = make_receiver_delta(state.expression, baseline_state=params.baseline_state)
    edges = add_receiver_response_scores(edges, receiver_delta, response_matrix)
    scores = finalize_scores(edges)
    if params.top_n is not None:
        scores = scores.head(params.top_n).copy()

    params.output.parent.mkdir(parents=True, exist_ok=True)
    scores.to_csv(params.output, sep="\t", index=False)
    if params.state_pair_output:
        params.state_pair_output.parent.mkdir(parents=True, exist_ok=True)
        summarize_state_pairs(scores).to_csv(params.state_pair_output, sep="\t", index=False)
    if params.state_expression_output:
        params.state_expression_output.parent.mkdir(parents=True, exist_ok=True)
        state.expression.to_csv(params.state_expression_output, sep="\t")

    summary = {
        "output": str(params.output),
        "n_cells": int(expression.shape[0]),
        "n_genes": int(expression.shape[1]),
        "n_states": int(state.expression.shape[0]),
        "n_lr_edges": int(lr_table.shape[0]),
        "n_scored_edges": int(scores.shape[0]),
        "n_loaded_genes": int(expression.shape[1]),
        **gene_diagnostics,
        "min_ligand_expr": float(params.min_ligand_expr),
        "min_receptor_expr": float(params.min_receptor_expr),
        "min_lr_expression_score": float(params.min_lr_expression_score),
        "state_key": params.state_key,
        "response_matrix": str(params.response_matrix) if params.response_matrix else None,
        "state_pair_output": str(params.state_pair_output) if params.state_pair_output else None,
        "state_expression_output": str(params.state_expression_output) if params.state_expression_output else None,
    }
    return FastCommResult(
        scores=scores,
        state_expression=state.expression,
        state_sizes=state.sizes,
        summary=summary,
    )
