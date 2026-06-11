from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd

from .scoring import (
    add_receiver_response_scores,
    clean_labels,
    finalize_scores,
    limit_lr_candidates_per_state_pair,
    load_response_matrix,
    make_receiver_delta,
    make_state_pseudobulk,
    make_state_pseudobulk_from_matrix,
    prepare_ligand_receptor_table,
    score_ligand_receptor_expression,
    split_complex,
    to_dense_frame,
)
from .upstream_resources import bundle_paths_for_species


PACKAGE_DIR = Path(__file__).resolve().parent
DEFAULT_LR_TABLE = PACKAGE_DIR / "resources" / "__auto__ligand_receptor.tsv"
DEFAULT_RESPONSE_MATRIX = PACKAGE_DIR / "resources" / "__auto__response_signatures.tsv"


@dataclass(frozen=True)
class FastCommParams:
    expression: Optional[Path] = None
    metadata: Optional[Path] = None
    h5ad: Optional[Path] = None
    adata: Optional[object] = None
    output: Path = Path("fastcomm_scores.tsv")
    lr_table: Path = DEFAULT_LR_TABLE
    response_matrix: Optional[Path] = DEFAULT_RESPONSE_MATRIX
    lr_sources: Optional[Sequence[str]] = None
    state_pair_output: Optional[Path] = None
    state_expression_output: Optional[Path] = None
    state_key: str = "cell_state"
    layer: Optional[str] = None
    gene_symbol_col: Optional[str] = None
    species: Optional[str] = None
    min_cells: int = 1
    min_ligand_expr: float = 0.01
    min_receptor_expr: float = 0.01
    min_lr_expression_score: float = 0.001
    max_lr_candidates_per_state_pair: int = 25
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


def _looks_human_symbol(gene: str) -> bool:
    letters = re.sub(r"[^A-Za-z]", "", gene)
    return bool(letters) and letters.isupper()


def _looks_mouse_symbol(gene: str) -> bool:
    letters = re.sub(r"[^A-Za-z]", "", gene)
    return len(letters) >= 2 and letters[0].isupper() and letters[1:].islower()


def infer_species_from_gene_symbols(genes: Sequence[str]) -> Optional[str]:
    human = 0
    mouse = 0
    examined = 0
    for gene in genes[:500]:
        gene = str(gene).strip()
        if _looks_human_symbol(gene):
            human += 1
            examined += 1
        elif _looks_mouse_symbol(gene):
            mouse += 1
            examined += 1
    if human and human == examined:
        return "human"
    if mouse and mouse == examined:
        return "mouse"
    if human >= 10 and human >= mouse * 2:
        return "human"
    if mouse >= 10 and mouse >= human * 2:
        return "mouse"
    return None


def _peek_input_genes(params: FastCommParams) -> list[str]:
    if params.adata is not None:
        adata = params.adata
        if params.gene_symbol_col and params.gene_symbol_col in adata.var.columns:
            return adata.var[params.gene_symbol_col].astype(str).tolist()
        return adata.var_names.astype(str).tolist()
    if params.h5ad is not None:
        try:
            import anndata as ad
        except ImportError:
            return []
        adata = ad.read_h5ad(params.h5ad, backed="r")
        try:
            if params.gene_symbol_col and params.gene_symbol_col in adata.var.columns:
                genes = adata.var[params.gene_symbol_col].astype(str).tolist()
            else:
                genes = adata.var_names.astype(str).tolist()
        finally:
            adata.file.close()
        return genes
    if params.expression is not None:
        header = pd.read_csv(params.expression, sep=None, engine="python", nrows=0, index_col=0)
        return [str(column) for column in header.columns]
    return []


def _resolve_default_resource_paths(params: FastCommParams) -> tuple[Path, Optional[Path], Optional[str]]:
    using_default_lr = params.lr_table == DEFAULT_LR_TABLE
    using_default_response = params.response_matrix == DEFAULT_RESPONSE_MATRIX
    if not using_default_lr and not using_default_response:
        return params.lr_table, params.response_matrix, params.species

    species = (params.species or "").strip().lower() or infer_species_from_gene_symbols(_peek_input_genes(params))
    if species not in {"human", "mouse"}:
        if using_default_lr:
            raise ValueError(
                "No default fastComm ligand-receptor resource is bundled. "
                "Run `fastComm build-upstream-resources --species human|mouse` or supply --lr-table explicitly."
            )
        return params.lr_table, None if using_default_response else params.response_matrix, species or None

    bundle = bundle_paths_for_species(species)
    if using_default_lr:
        if not bundle["lr_table"].exists():
            raise ValueError(
                f"No built fastComm ligand-receptor bundle was found for species={species!r}. "
                "Run `fastComm build-upstream-resources` first or supply --lr-table explicitly."
            )
        lr_table = bundle["lr_table"]
    else:
        lr_table = params.lr_table
    response_matrix = bundle["response_matrix"] if using_default_response and bundle["response_matrix"].exists() else params.response_matrix
    if using_default_response and response_matrix == DEFAULT_RESPONSE_MATRIX:
        response_matrix = None
    return lr_table, response_matrix, species


def _materialize_matrix(matrix) -> np.ndarray:
    if hasattr(matrix, "toarray"):
        return matrix.toarray()
    return np.asarray(matrix)


def _slice_matrix_from_adata(adata, params: FastCommParams, gene_positions: Sequence[int]):
    matrix_source = adata.layers[params.layer] if params.layer else adata.X
    return matrix_source[:, gene_positions]


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


def _filter_lr_sources(lr_table: pd.DataFrame, lr_sources: Optional[Sequence[str]]) -> pd.DataFrame:
    if not lr_sources:
        return lr_table
    if "source" not in lr_table.columns:
        raise ValueError("LR source filtering requires a 'source' column in the ligand-receptor table")
    allowed = {str(source).strip().lower() for source in lr_sources if str(source).strip()}
    if not allowed:
        return lr_table
    keep = lr_table["source"].astype(str).str.strip().str.lower().isin(allowed)
    filtered = lr_table.loc[keep].copy()
    if filtered.empty:
        raise ValueError(f"No ligand-receptor rows matched lr_sources={sorted(allowed)}")
    return filtered


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


def _matrix_inputs_from_h5ad(
    params: FastCommParams,
    *,
    required_genes: Optional[Sequence[str]] = None,
) -> tuple[object, list[str], list[str], pd.DataFrame, Dict[str, object]]:
    if params.h5ad is None:
        raise ValueError("h5ad path is required")
    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover - exercised only in minimal envs
        raise ImportError("AnnData input requires the 'anndata' package") from exc

    adata = ad.read_h5ad(params.h5ad, backed="r")
    try:
        if params.gene_symbol_col:
            if params.gene_symbol_col not in adata.var.columns:
                raise KeyError(f"Column {params.gene_symbol_col!r} was not found in adata.var")
            genes = adata.var[params.gene_symbol_col].astype(str).tolist()
        else:
            genes = adata.var_names.astype(str).tolist()

        gene_positions, selected_genes, diagnostics = _select_gene_positions(genes, required_genes)
        if required_genes and not gene_positions:
            raise ValueError("None of the requested ligand/receptor/response genes were found in the h5ad")

        index = adata.obs_names.astype(str).tolist()
        metadata = adata.obs.copy()
        try:
            matrix = _slice_matrix_from_adata(adata, params, gene_positions)
        except AttributeError as exc:
            if "_validate_indices" not in str(exc):
                raise
            matrix_source = adata.layers[params.layer] if params.layer else adata.X
            if hasattr(matrix_source, "to_memory"):
                matrix = matrix_source.to_memory()[:, gene_positions]
            else:
                adata.file.close()
                adata = ad.read_h5ad(params.h5ad)
                matrix = _slice_matrix_from_adata(adata, params, gene_positions)
    finally:
        if getattr(adata, "isbacked", False):
            adata.file.close()
    diagnostics["n_loaded_columns"] = int(len(selected_genes))
    return matrix, index, selected_genes, metadata, diagnostics


def _matrix_inputs_from_adata(
    params: FastCommParams,
    *,
    required_genes: Optional[Sequence[str]] = None,
) -> tuple[object, list[str], list[str], pd.DataFrame, Dict[str, object]]:
    if params.adata is None:
        raise ValueError("adata is required")
    adata = params.adata
    if params.gene_symbol_col:
        if params.gene_symbol_col not in adata.var.columns:
            raise KeyError(f"Column {params.gene_symbol_col!r} was not found in adata.var")
        genes = adata.var[params.gene_symbol_col].astype(str).tolist()
    else:
        genes = adata.var_names.astype(str).tolist()

    gene_positions, selected_genes, diagnostics = _select_gene_positions(genes, required_genes)
    if required_genes and not gene_positions:
        raise ValueError("None of the requested ligand/receptor/response genes were found in the adata")

    matrix = _slice_matrix_from_adata(adata, params, gene_positions)
    diagnostics["n_loaded_columns"] = int(len(selected_genes))
    return (
        matrix,
        adata.obs_names.astype(str).tolist(),
        selected_genes,
        adata.obs.copy(),
        diagnostics,
    )


def _expression_from_h5ad(
    params: FastCommParams,
    *,
    required_genes: Optional[Sequence[str]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    matrix, index, selected_genes, metadata, diagnostics = _matrix_inputs_from_h5ad(
        params,
        required_genes=required_genes,
    )
    expression = to_dense_frame(_materialize_matrix(matrix), index=index, columns=selected_genes)
    expression = _deduplicate_columns(expression)
    return expression, metadata, diagnostics


def _load_inputs(
    params: FastCommParams,
    *,
    required_genes: Optional[Sequence[str]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    if params.adata is not None:
        matrix, index, selected_genes, metadata, diagnostics = _matrix_inputs_from_adata(
            params,
            required_genes=required_genes,
        )
        expression = to_dense_frame(_materialize_matrix(matrix), index=index, columns=selected_genes)
        expression = _deduplicate_columns(expression)
        return expression, metadata, diagnostics
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
    resolved_lr_table, resolved_response_path, inferred_species = _resolve_default_resource_paths(params)
    lr_table = _filter_lr_sources(pd.read_csv(resolved_lr_table, sep="\t"), params.lr_sources)
    response_matrix = load_response_matrix(str(resolved_response_path)) if resolved_response_path else None
    required_genes = _required_genes(lr_table, response_matrix)
    if params.adata is not None or params.h5ad is not None:
        if params.adata is not None:
            matrix, index, columns, metadata, gene_diagnostics = _matrix_inputs_from_adata(
                params,
                required_genes=required_genes,
            )
        else:
            matrix, index, columns, metadata, gene_diagnostics = _matrix_inputs_from_h5ad(
                params,
                required_genes=required_genes,
            )
        state = make_state_pseudobulk_from_matrix(
            matrix,
            index=index,
            columns=columns,
            metadata=metadata,
            state_key=params.state_key,
            min_cells=params.min_cells,
        )
        state = type(state)(
            expression=_deduplicate_columns(state.expression),
            detection=_deduplicate_columns(state.detection),
            sizes=state.sizes,
        )
        n_input_cells = len(index)
        n_input_genes = len(columns)
    else:
        expression, metadata, gene_diagnostics = _load_inputs(params, required_genes=required_genes)
        state = make_state_pseudobulk(
            expression,
            metadata,
            state_key=params.state_key,
            min_cells=params.min_cells,
        )
        n_input_cells = int(expression.shape[0])
        n_input_genes = int(expression.shape[1])
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
        "n_cells": n_input_cells,
        "n_genes": n_input_genes,
        "n_states": int(state.expression.shape[0]),
        "n_lr_edges": int(lr_table.shape[0]),
        "n_scored_edges": int(scores.shape[0]),
        "n_loaded_genes": n_input_genes,
        **gene_diagnostics,
        "min_ligand_expr": float(params.min_ligand_expr),
        "min_receptor_expr": float(params.min_receptor_expr),
        "min_lr_expression_score": float(params.min_lr_expression_score),
        "max_lr_candidates_per_state_pair": int(params.max_lr_candidates_per_state_pair),
        "state_key": params.state_key,
        "species": inferred_species,
        "lr_table": str(resolved_lr_table),
        "lr_sources": list(params.lr_sources) if params.lr_sources else None,
        "response_matrix": str(resolved_response_path) if resolved_response_path else None,
        "state_pair_output": str(params.state_pair_output) if params.state_pair_output else None,
        "state_expression_output": str(params.state_expression_output) if params.state_expression_output else None,
    }
    return FastCommResult(
        scores=scores,
        state_expression=state.expression,
        state_sizes=state.sizes,
        summary=summary,
    )
