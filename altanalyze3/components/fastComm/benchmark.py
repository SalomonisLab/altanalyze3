from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import time
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np

from .api import (
    DEFAULT_LR_TABLE,
    DEFAULT_RESPONSE_MATRIX,
    FastCommParams,
    _matrix_inputs_from_h5ad,
    _deduplicate_columns,
    _filter_lr_sources,
    _required_genes,
    _resolve_default_resource_paths,
)
from .scoring import (
    add_receiver_response_scores,
    finalize_scores,
    load_response_matrix,
    limit_lr_candidates_per_state_pair,
    make_receiver_delta,
    make_state_pseudobulk,
    PseudobulkState,
    clean_labels,
    to_dense_frame,
    score_ligand_receptor_expression,
)


@dataclass(frozen=True)
class FastCommBenchmarkParams:
    h5ad: Path
    output_dir: Path
    state_key: str = "cell_state"
    split_key: str = "Donor"
    lr_table: Path = DEFAULT_LR_TABLE
    response_matrix: Optional[Path] = DEFAULT_RESPONSE_MATRIX
    lr_sources: Optional[List[str]] = None
    layer: Optional[str] = None
    gene_symbol_col: Optional[str] = None
    species: Optional[str] = None
    min_cells: int = 20
    min_ligand_expr: float = 0.01
    min_receptor_expr: float = 0.01
    min_lr_expression_score: float = 0.001
    max_lr_candidates_per_state_pair: int = 25
    include_self_edges: bool = False
    top_n_stability: int = 100


def _edge_key_frame(scores: pd.DataFrame) -> pd.Series:
    if scores.empty:
        return pd.Series(dtype=str)
    return (
        scores["sender_state"].astype(str)
        + "|"
        + scores["receiver_state"].astype(str)
        + "|"
        + scores["ligand"].astype(str)
        + "|"
        + scores["receptor"].astype(str)
    )


def _score_subset(
    expression,
    metadata: pd.DataFrame,
    *,
    lr_table: pd.DataFrame,
    response_matrix: Optional[pd.DataFrame],
    params: FastCommBenchmarkParams,
    matrix_index=None,
    matrix_columns=None,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    started = time.perf_counter()
    if matrix_index is None:
        state = make_state_pseudobulk(expression, metadata,
                                      state_key=params.state_key, min_cells=params.min_cells)
        n_cells, n_genes = expression.shape
    else:
        state = _batched_state(expression, matrix_index, matrix_columns, metadata,
                               state_key=params.state_key, min_cells=params.min_cells)
        n_cells = len(pd.Index(matrix_index).intersection(metadata.index.astype(str)))
        n_genes = len(set(clean_labels(matrix_columns)))
    edges = score_ligand_receptor_expression(
        state,
        lr_table,
        min_ligand_expr=params.min_ligand_expr,
        min_receptor_expr=params.min_receptor_expr,
        min_lr_expression_score=params.min_lr_expression_score,
        include_self_edges=params.include_self_edges,
    )
    if params.min_lr_expression_score > 0 and not edges.empty:
        edges = edges.loc[edges["lr_expression_score"] >= params.min_lr_expression_score].copy()
    edges = limit_lr_candidates_per_state_pair(edges, params.max_lr_candidates_per_state_pair)
    receiver_delta = make_receiver_delta(state.expression)
    edges = add_receiver_response_scores(edges, receiver_delta, response_matrix)
    scores = finalize_scores(edges)
    elapsed = time.perf_counter() - started
    summary = {
        "n_cells": int(n_cells),
        "n_loaded_genes": int(n_genes),
        "n_states": int(state.expression.shape[0]),
        "n_edges": int(scores.shape[0]),
        "elapsed_seconds": round(elapsed, 4),
    }
    if not scores.empty:
        top = scores.iloc[0]
        summary.update(
            {
                "top_interaction": f"{top['sender_state']}->{top['receiver_state']}:{top['ligand']}->{top['receptor']}",
                "top_fastcomm_score": float(top["fastcomm_score"]),
            }
        )
    return scores, summary


def _batched_state(matrix, index, columns, metadata, *, state_key, min_cells, block_genes=128):
    """Original pandas means/detection, processing bounded feature panels.

    Every retained cell enters each calculation in its original order. Unlike
    adding partial row sums, this preserves pandas' floating-point accumulation
    and near-tied interaction ranks. Duplicate symbols are averaged per cell
    before detection, exactly as in the original dense benchmark.
    """
    if state_key not in metadata:
        raise KeyError(f"State column {state_key!r} was not found in metadata")
    obs_index = pd.Index(clean_labels(index))
    metadata = metadata.copy()
    metadata.index = clean_labels(metadata.index)
    common = obs_index.intersection(metadata.index)
    if common.empty:
        raise ValueError("No shared cells between expression index and metadata index")
    states = metadata.loc[common, state_key].astype(str).str.strip()
    valid = states.ne("") & states.notna()
    common, states = common[valid], states.loc[valid]
    sizes = states.value_counts().sort_index()
    keep = sizes.index[sizes >= min_cells]
    if keep.empty:
        raise ValueError(f"No states passed min_cells={min_cells}")
    valid = states.isin(keep)
    common, states = common[valid], states.loc[valid]
    positions = obs_index.get_indexer(common)
    column_labels = pd.Index(clean_labels(columns))
    names = column_labels.drop_duplicates()
    # Sparse rows stay sparse. No all-feature dense single-cell table exists.
    subset = matrix[positions]
    means, detections = [], []
    if not len(names):
        empty = pd.DataFrame(index=keep, columns=names, dtype=float)
        return PseudobulkState(empty, empty.copy(), sizes.loc[keep])
    for start in range(0, len(names), block_genes):
        selected = names[start:start + block_genes]
        gene_positions = np.flatnonzero(column_labels.isin(selected))
        block = _deduplicate_columns(to_dense_frame(subset[:, gene_positions],
                                                     index=common, columns=column_labels[gene_positions]))
        means.append(block.groupby(states, sort=True).mean().astype(float))
        detections.append(block.gt(0).groupby(states, sort=True).mean().astype(float))
    return PseudobulkState(pd.concat(means, axis=1).reindex(columns=names),
                           pd.concat(detections, axis=1).reindex(columns=names), sizes.loc[keep])


def _compare_to_full(full_scores: pd.DataFrame, split_scores: pd.DataFrame, *, top_n: int) -> Dict[str, object]:
    if full_scores.empty or split_scores.empty:
        return {
            "top_jaccard_vs_full": 0.0,
            "shared_edges_vs_full": 0,
            "score_corr_vs_full": 0.0,
        }

    full = full_scores.copy()
    split = split_scores.copy()
    full["edge_key"] = _edge_key_frame(full)
    split["edge_key"] = _edge_key_frame(split)

    full_top = set(full.head(top_n)["edge_key"])
    split_top = set(split.head(top_n)["edge_key"])
    union = full_top | split_top
    jaccard = len(full_top & split_top) / len(union) if union else 0.0

    merged = full[["edge_key", "fastcomm_score"]].merge(
        split[["edge_key", "fastcomm_score"]],
        on="edge_key",
        suffixes=("_full", "_split"),
    )
    corr = 0.0
    if merged.shape[0] >= 2:
        corr_value = merged["fastcomm_score_full"].corr(merged["fastcomm_score_split"], method="spearman")
        corr = 0.0 if pd.isna(corr_value) else float(corr_value)
    return {
        "top_jaccard_vs_full": float(jaccard),
        "shared_edges_vs_full": int(merged.shape[0]),
        "score_corr_vs_full": corr,
    }


def run_benchmark(params: FastCommBenchmarkParams) -> Dict[str, object]:
    params.output_dir.mkdir(parents=True, exist_ok=True)
    resolved_lr_table, resolved_response_path, inferred_species = _resolve_default_resource_paths(
        FastCommParams(
            h5ad=params.h5ad,
            lr_table=params.lr_table,
            response_matrix=params.response_matrix,
            layer=params.layer,
            gene_symbol_col=params.gene_symbol_col,
            species=params.species,
        )
    )
    lr_table = _filter_lr_sources(pd.read_csv(resolved_lr_table, sep="\t"), params.lr_sources)
    response_matrix = load_response_matrix(str(resolved_response_path)) if resolved_response_path else None
    required_genes = _required_genes(lr_table, response_matrix)
    expression, matrix_index, matrix_columns, metadata, gene_diagnostics = _matrix_inputs_from_h5ad(
        FastCommParams(
            h5ad=params.h5ad,
            lr_table=resolved_lr_table,
            response_matrix=resolved_response_path,
            output=params.output_dir / "_unused.tsv",
            state_key=params.state_key,
            layer=params.layer,
            gene_symbol_col=params.gene_symbol_col,
            species=inferred_species,
        ),
        required_genes=required_genes,
    )

    full_scores, full_summary = _score_subset(
        expression,
        metadata,
        lr_table=lr_table,
        response_matrix=response_matrix,
        params=params,
        matrix_index=matrix_index,
        matrix_columns=matrix_columns,
    )
    full_scores.to_csv(params.output_dir / "full_scores.tsv", sep="\t", index=False)

    split_rows: List[Dict[str, object]] = []
    split_score_frames: List[pd.DataFrame] = []
    if params.split_key not in metadata.columns:
        raise KeyError(f"Split column {params.split_key!r} was not found in h5ad obs")

    for split_name, split_metadata in metadata.groupby(metadata[params.split_key].astype(str), sort=True):
        if len(split_metadata) < params.min_cells:
            continue
        try:
            split_scores, split_summary = _score_subset(
                expression,
                split_metadata,
                lr_table=lr_table,
                response_matrix=response_matrix,
                params=params,
                matrix_index=matrix_index,
                matrix_columns=matrix_columns,
            )
        except ValueError as exc:
            split_rows.append(
                {
                    "split": split_name,
                    "status": "skipped",
                    "reason": str(exc),
                    "n_cells": int(len(split_metadata)),
                }
            )
            continue
        split_scores.to_csv(params.output_dir / f"split_{split_name}_scores.tsv", sep="\t", index=False)
        if not split_scores.empty:
            split_frame = split_scores.copy()
            split_frame.insert(0, "split", str(split_name))
            split_score_frames.append(split_frame)
        row = {"split": split_name, "status": "ok"}
        row.update(split_summary)
        row.update(_compare_to_full(full_scores, split_scores, top_n=params.top_n_stability))
        split_rows.append(row)

    split_summary_df = pd.DataFrame(split_rows)
    split_summary_df.to_csv(params.output_dir / "split_stability.tsv", sep="\t", index=False)
    split_scores_long_path = params.output_dir / "split_scores_long.tsv"
    if split_score_frames:
        pd.concat(split_score_frames, ignore_index=True).to_csv(split_scores_long_path, sep="\t", index=False)
    else:
        pd.DataFrame(columns=["split"]).to_csv(split_scores_long_path, sep="\t", index=False)

    summary = {
        "h5ad": str(params.h5ad),
        "output_dir": str(params.output_dir),
        "split_scores_long_tsv": str(split_scores_long_path),
        "state_key": params.state_key,
        "split_key": params.split_key,
        "species": inferred_species,
        "lr_table": str(resolved_lr_table),
        "response_matrix": str(resolved_response_path) if resolved_response_path else None,
        "top_n_stability": int(params.top_n_stability),
        "full": full_summary,
        "n_splits": int(split_summary_df.shape[0]),
        "loaded_genes": int(len(set(clean_labels(matrix_columns)))),
        **gene_diagnostics,
    }
    (params.output_dir / "benchmark_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary
