from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import time
from typing import Dict, List, Optional, Tuple

import pandas as pd

from .api import (
    DEFAULT_LR_TABLE,
    DEFAULT_RESPONSE_MATRIX,
    FastCommParams,
    _expression_from_h5ad,
    _required_genes,
    _resolve_default_resource_paths,
)
from .scoring import (
    add_receiver_response_scores,
    finalize_scores,
    load_response_matrix,
    make_receiver_delta,
    make_state_pseudobulk,
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
    layer: Optional[str] = None
    gene_symbol_col: Optional[str] = None
    species: Optional[str] = None
    min_cells: int = 20
    min_ligand_expr: float = 0.01
    min_receptor_expr: float = 0.01
    min_lr_expression_score: float = 0.001
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
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    lr_table: pd.DataFrame,
    response_matrix: Optional[pd.DataFrame],
    params: FastCommBenchmarkParams,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    started = time.perf_counter()
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
    receiver_delta = make_receiver_delta(state.expression)
    edges = add_receiver_response_scores(edges, receiver_delta, response_matrix)
    scores = finalize_scores(edges)
    elapsed = time.perf_counter() - started
    summary = {
        "n_cells": int(expression.shape[0]),
        "n_loaded_genes": int(expression.shape[1]),
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
    lr_table = pd.read_csv(resolved_lr_table, sep=None, engine="python")
    response_matrix = load_response_matrix(str(resolved_response_path)) if resolved_response_path else None
    required_genes = _required_genes(lr_table, response_matrix)
    expression, metadata, gene_diagnostics = _expression_from_h5ad(
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
    )
    full_scores.to_csv(params.output_dir / "full_scores.tsv", sep="\t", index=False)

    split_rows: List[Dict[str, object]] = []
    split_score_frames: List[pd.DataFrame] = []
    if params.split_key not in metadata.columns:
        raise KeyError(f"Split column {params.split_key!r} was not found in h5ad obs")

    for split_name, split_metadata in metadata.groupby(metadata[params.split_key].astype(str), sort=True):
        split_cells = split_metadata.index.astype(str)
        split_expression = expression.loc[split_cells]
        if split_expression.shape[0] < params.min_cells:
            continue
        try:
            split_scores, split_summary = _score_subset(
                split_expression,
                split_metadata,
                lr_table=lr_table,
                response_matrix=response_matrix,
                params=params,
            )
        except ValueError as exc:
            split_rows.append(
                {
                    "split": split_name,
                    "status": "skipped",
                    "reason": str(exc),
                    "n_cells": int(split_expression.shape[0]),
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
        "loaded_genes": int(expression.shape[1]),
        **gene_diagnostics,
    }
    (params.output_dir / "benchmark_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary
