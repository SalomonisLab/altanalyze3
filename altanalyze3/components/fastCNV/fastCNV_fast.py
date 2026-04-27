"""Fast variant of fastCNV: global preprocess + per-chromosome vectorized scoring.

Same scoring model as ``main.py`` (median baseline, 1.4826 * MAD scale, rolling
window mean of residuals, run-length interval calling), but reorganized so the
expensive chunk extractions and normalizations happen once per chromosome
instead of once per (state x chromosome). Designed for runs that take hours
in the per-state implementation.

Differences vs ``main.py``:
  * Requires a control h5ad. Internal-anchor mode is not implemented here -
    use ``main.py`` if you need it.
  * Single-pass scoring (no burden anchor refinement; the control already
    provides the WT baseline).
  * Skips clone discovery and the PDF by default. The npz score matrix lets
    you re-run those downstream.
"""

from __future__ import annotations

import argparse
import logging
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from altanalyze3.components.fastCNV.main import (
    RESOURCE_FILES,
    Window,
    _call_intervals_for_cell,
    _chr_slices,
    _cnv_burden,
    _format_duration,
    _window_frame,
    bundled_gene_coordinates,
    build_windows,
    load_gene_coordinates,
)


LOGGER = logging.getLogger("fastCNV.fast")


@dataclass
class FastParams:
    h5ad: Path
    control_h5ad: Path
    gene_coordinates: Optional[Path]
    output_prefix: Path
    state_key: str
    sample_key: Optional[str] = None
    control_state_key: Optional[str] = None
    layer: str = "auto"
    input_normalized: bool = False
    window_genes: int = 41
    stride_genes: int = 7
    min_chr_genes: int = 25
    min_state_cells: int = 30
    high_threshold: float = 2.6
    low_threshold: float = 1.6
    min_run_windows: int = 3
    min_interval_genes: int = 60
    min_mean_score: float = 1.8
    burden_quantile: float = 0.95
    cnv_burden_threshold: float = 1.8
    write_h5ad: bool = False


def _matrix(adata: ad.AnnData, layer: str) -> sp.spmatrix | np.ndarray:
    if layer == "auto":
        return adata.layers["counts"] if "counts" in adata.layers else adata.X
    if layer == "X":
        return adata.X
    if layer not in adata.layers:
        raise KeyError(f"Layer '{layer}' not found in AnnData.")
    return adata.layers[layer]


def _row_sums(matrix: sp.spmatrix | np.ndarray) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32)


def _normalize_chunk(
    matrix: sp.spmatrix | np.ndarray,
    col_indices: np.ndarray,
    library_sizes: np.ndarray,
    input_normalized: bool,
) -> np.ndarray:
    """Densify all rows for the given column indices and log-normalize in place."""
    sub = matrix[:, col_indices]
    if sp.issparse(sub):
        sub = sub.toarray()
    sub = np.asarray(sub, dtype=np.float32)
    if input_normalized:
        return sub
    factors = (10000.0 / np.maximum(library_sizes.astype(np.float32), 1.0)).astype(np.float32)
    sub *= factors[:, None]
    np.log1p(sub, out=sub)
    return sub


def _rolling_window_mean(values: np.ndarray, windows: Sequence[Window]) -> np.ndarray:
    if values.shape[1] == 0 or not windows:
        return np.zeros((values.shape[0], 0), dtype=np.float32)
    prefix = np.empty((values.shape[0], values.shape[1] + 1), dtype=np.float32)
    prefix[:, 0] = 0.0
    np.cumsum(values, axis=1, dtype=np.float32, out=prefix[:, 1:])
    starts = np.asarray([w.start_offset for w in windows], dtype=np.int64)
    ends = np.asarray([w.end_offset for w in windows], dtype=np.int64)
    lengths = np.maximum(ends - starts, 1).astype(np.float32)
    return ((prefix[:, ends] - prefix[:, starts]) / lengths).astype(np.float32, copy=False)


def _mad(values: np.ndarray, axis: int = 0) -> np.ndarray:
    med = np.nanmedian(values, axis=axis)
    return np.nanmedian(np.abs(values - np.expand_dims(med, axis)), axis=axis)


def _control_var_map(query_var: pd.Index, control_var: pd.Index) -> np.ndarray:
    lookup = pd.Series(np.arange(len(control_var), dtype=np.int64), index=control_var.astype(str))
    return lookup.reindex(query_var.astype(str)).to_numpy()


def run_fast(params: FastParams) -> Dict[str, Path]:
    run_start = time.perf_counter()
    params.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    LOGGER.info("Reading query AnnData: %s", params.h5ad)
    t0 = time.perf_counter()
    query = ad.read_h5ad(params.h5ad)
    LOGGER.info("Query loaded: %d cells x %d genes (%s)", query.n_obs, query.n_vars, _format_duration(time.perf_counter() - t0))

    available = ", ".join(query.obs.columns.astype(str)) or "<none>"
    if params.state_key not in query.obs.columns:
        raise KeyError(f"State column '{params.state_key}' not found. Available: {available}")
    if params.sample_key and params.sample_key not in query.obs.columns:
        raise KeyError(f"Sample column '{params.sample_key}' not found. Available: {available}")

    LOGGER.info("Reading control AnnData: %s", params.control_h5ad)
    t0 = time.perf_counter()
    control = ad.read_h5ad(params.control_h5ad)
    LOGGER.info("Control loaded: %d cells x %d genes (%s)", control.n_obs, control.n_vars, _format_duration(time.perf_counter() - t0))

    if params.control_state_key and params.control_state_key not in control.obs.columns:
        raise KeyError(
            f"Control state column '{params.control_state_key}' not found. "
            f"Available: {', '.join(control.obs.columns.astype(str))}"
        )

    if params.gene_coordinates is None:
        raise ValueError("gene_coordinates must be provided.")
    coords = load_gene_coordinates(params.gene_coordinates, query.var_names)
    windows = build_windows(coords, params.window_genes, params.stride_genes, params.min_chr_genes)
    slices_by_chr = _chr_slices(windows)
    LOGGER.info("Built %d windows across %d chromosomes.", len(windows), len(slices_by_chr))

    query_matrix = _matrix(query, params.layer)
    if sp.issparse(query_matrix):
        query_matrix = query_matrix.tocsc()
    control_matrix = _matrix(control, params.layer)
    if sp.issparse(control_matrix):
        control_matrix = control_matrix.tocsc()

    LOGGER.info("Computing library sizes.")
    t0 = time.perf_counter()
    query_lib = np.ones(query.n_obs, dtype=np.float32) if params.input_normalized else _row_sums(query_matrix.tocsr() if sp.issparse(query_matrix) else query_matrix)
    control_lib = np.ones(control.n_obs, dtype=np.float32) if params.input_normalized else _row_sums(control_matrix.tocsr() if sp.issparse(control_matrix) else control_matrix)
    LOGGER.info("Library sizes done (%s).", _format_duration(time.perf_counter() - t0))

    var_map = _control_var_map(query.var_names, control.var_names)

    state_values = query.obs[params.state_key]
    states = [str(s) for s in pd.Index(state_values[state_values.notna()]).unique().tolist() if str(s)]
    state_index_map: Dict[str, np.ndarray] = {
        state: np.flatnonzero(state_values.astype(str).to_numpy() == state).astype(np.int64)
        for state in states
    }
    eligible_states = [s for s in states if state_index_map[s].size >= params.min_state_cells]
    LOGGER.info("States: %d total, %d eligible (>= %d cells).", len(states), len(eligible_states), params.min_state_cells)

    control_state_pool: Optional[Dict[str, np.ndarray]] = None
    if params.control_state_key:
        control_state_values = control.obs[params.control_state_key].astype(str).to_numpy()
        control_state_pool = {}
        for state in eligible_states:
            mask = control_state_values == state
            control_state_pool[state] = np.flatnonzero(mask).astype(np.int64) if mask.any() else np.arange(control.n_obs, dtype=np.int64)
    all_control_rows = np.arange(control.n_obs, dtype=np.int64)

    all_scores = np.full((query.n_obs, len(windows)), np.nan, dtype=np.float32)
    LOGGER.info(
        "Allocated score matrix: %d x %d float32 (~%.1f GB).",
        query.n_obs, len(windows), all_scores.nbytes / 1e9,
    )

    chrom_groups = list(coords.groupby("chr", sort=False))
    LOGGER.info("Scoring %d chromosomes.", len(chrom_groups))

    for chr_index, (chrom, chr_coords) in enumerate(chrom_groups):
        chrom_str = str(chrom)
        chr_slice = slices_by_chr.get(chrom_str)
        if chr_slice is None:
            continue
        chr_windows = windows[chr_slice]
        if not chr_windows:
            continue
        chr_t0 = time.perf_counter()

        query_cols = chr_coords["var_index"].to_numpy(dtype=np.int64)
        query_expr = _normalize_chunk(query_matrix, query_cols, query_lib, params.input_normalized)

        mapped = var_map[query_cols]
        valid_mask = ~pd.isna(mapped)
        valid_control_cols = mapped[valid_mask].astype(np.int64) if valid_mask.any() else np.empty(0, dtype=np.int64)
        if valid_control_cols.size == 0:
            LOGGER.warning("Chromosome %s: no overlapping control genes; skipping.", chrom_str)
            continue
        control_expr_full = np.zeros((control.n_obs, query_cols.size), dtype=np.float32)
        control_sub = _normalize_chunk(control_matrix, valid_control_cols, control_lib, params.input_normalized)
        control_expr_full[:, valid_mask] = control_sub
        del control_sub

        if control_state_pool is None:
            baseline = np.nanmedian(control_expr_full[all_control_rows], axis=0).astype(np.float32)
            residual_query = query_expr - baseline[None, :]
            residual_control = control_expr_full - baseline[None, :]
            window_query = _rolling_window_mean(residual_query, chr_windows)
            window_control = _rolling_window_mean(residual_control, chr_windows)
            center = np.nanmedian(window_control, axis=0).astype(np.float32)
            scale = (1.4826 * _mad(window_control, axis=0)).astype(np.float32)
            scale = np.where(scale < 0.10, 0.10, scale).astype(np.float32)
            scores_chr = ((window_query - center[None, :]) / scale[None, :]).astype(np.float32)
            for state in eligible_states:
                rows = state_index_map[state]
                all_scores[rows[:, None], np.arange(chr_slice.start, chr_slice.stop)[None, :]] = scores_chr[rows]
        else:
            chr_window_count = len(chr_windows)
            chr_abs_indices = np.arange(chr_slice.start, chr_slice.stop)
            for state in eligible_states:
                rows = state_index_map[state]
                control_rows = control_state_pool.get(state, all_control_rows)
                if control_rows.size == 0:
                    control_rows = all_control_rows
                baseline = np.nanmedian(control_expr_full[control_rows], axis=0).astype(np.float32)
                state_residual = query_expr[rows] - baseline[None, :]
                control_residual = control_expr_full[control_rows] - baseline[None, :]
                window_state = _rolling_window_mean(state_residual, chr_windows)
                window_control = _rolling_window_mean(control_residual, chr_windows)
                center = np.nanmedian(window_control, axis=0).astype(np.float32)
                scale = (1.4826 * _mad(window_control, axis=0)).astype(np.float32)
                scale = np.where(scale < 0.10, 0.10, scale).astype(np.float32)
                state_scores = ((window_state - center[None, :]) / scale[None, :]).astype(np.float32)
                all_scores[rows[:, None], chr_abs_indices[None, :]] = state_scores

        del query_expr, control_expr_full
        LOGGER.info(
            "Chromosome %s [%d/%d]: %d windows scored (%s).",
            chrom_str, chr_index + 1, len(chrom_groups), len(chr_windows),
            _format_duration(time.perf_counter() - chr_t0),
        )

    LOGGER.info("Calling per-cell intervals.")
    interval_t0 = time.perf_counter()
    cell_records: List[Dict[str, object]] = []
    interval_records: List[Dict[str, object]] = []
    sample_series = query.obs[params.sample_key].astype(str) if params.sample_key else None
    state_series = query.obs[params.state_key].astype(str)

    for state in eligible_states:
        rows = state_index_map[state]
        scores_block = all_scores[rows]
        burden = _cnv_burden(scores_block, params.burden_quantile)
        for local_index, row in enumerate(rows):
            cell_scores = scores_block[local_index]
            calls = _call_intervals_for_cell(
                cell_scores, windows, slices_by_chr,
                high_threshold=params.high_threshold,
                low_threshold=params.low_threshold,
                min_run_windows=params.min_run_windows,
                min_interval_genes=params.min_interval_genes,
                min_mean_score=params.min_mean_score,
            )
            status = "CNV" if calls and burden[local_index] >= params.cnv_burden_threshold else "WT"
            barcode = str(query.obs_names[row])
            sample = sample_series.iloc[row] if sample_series is not None else ""
            cell_records.append({
                "CellBarcode": barcode,
                "cell_state": state,
                "sample": sample,
                "cnv_status": status,
                "cnv_burden": float(burden[local_index]),
                "baseline_source": "control",
                "n_cnv_intervals": len(calls) if status == "CNV" else 0,
            })
            if status == "CNV":
                for call in calls:
                    interval_records.append({
                        "CellBarcode": barcode,
                        "cell_state": state,
                        "sample": sample,
                        **call,
                    })

    skipped_states = [s for s in states if s not in set(eligible_states)]
    for state in skipped_states:
        for row in state_index_map[state]:
            cell_records.append({
                "CellBarcode": str(query.obs_names[row]),
                "cell_state": state,
                "sample": sample_series.iloc[row] if sample_series is not None else "",
                "cnv_status": "low_power",
                "cnv_burden": float("nan"),
                "baseline_source": "none",
                "n_cnv_intervals": 0,
            })
    LOGGER.info("Interval calling done (%s).", _format_duration(time.perf_counter() - interval_t0))

    cells_path = params.output_prefix.with_suffix(".cnv_cells.tsv")
    intervals_path = params.output_prefix.with_suffix(".cnv_intervals.tsv")
    windows_path = params.output_prefix.with_suffix(".cnv_windows.tsv")
    scores_path = params.output_prefix.with_suffix(".cnv_window_scores.npz")
    h5ad_path = params.output_prefix.with_suffix(".fastcnv.h5ad")

    pd.DataFrame(cell_records).to_csv(cells_path, sep="\t", index=False)
    interval_columns = [
        "CellBarcode", "cell_state", "sample", "call", "chr", "start", "end",
        "start_gene", "end_gene", "n_windows", "n_genes", "mean_score", "max_score", "confidence",
    ]
    pd.DataFrame(interval_records, columns=interval_columns).to_csv(intervals_path, sep="\t", index=False)
    _window_frame(windows).to_csv(windows_path, sep="\t", index=False)
    np.savez(
        scores_path,
        scores=all_scores,
        cell_barcodes=np.asarray(query.obs_names.astype(str), dtype=object),
        states=state_series.to_numpy(dtype=object),
        window_ids=np.asarray([w.window_id for w in windows], dtype=object),
    )

    outputs = {
        "cells": cells_path,
        "intervals": intervals_path,
        "windows": windows_path,
        "scores": scores_path,
    }
    if params.write_h5ad:
        cell_df = pd.DataFrame(cell_records).set_index("CellBarcode")
        for col in ("cnv_status", "cnv_burden"):
            query.obs[f"fastcnv_{col}"] = cell_df.reindex(query.obs_names.astype(str))[col].values
        query.write_h5ad(h5ad_path)
        outputs["h5ad"] = h5ad_path

    LOGGER.info("fastCNV (fast) total runtime: %s", _format_duration(time.perf_counter() - run_start))
    return outputs


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Fast variant of fastCNV: requires a control h5ad; vectorized per-chromosome scoring.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--h5ad", required=True, type=Path)
    p.add_argument("--control-h5ad", required=True, type=Path)
    p.add_argument("--gene-coordinates", default=None, type=Path)
    p.add_argument("--species", choices=sorted(RESOURCE_FILES), default=None)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--state-key", required=True)
    p.add_argument("--sample-key", default=None)
    p.add_argument("--control-state-key", default=None)
    p.add_argument("--layer", default="auto")
    p.add_argument("--input-normalized", action="store_true")
    p.add_argument("--window-genes", type=int, default=41)
    p.add_argument("--stride-genes", type=int, default=7)
    p.add_argument("--min-chr-genes", type=int, default=25)
    p.add_argument("--min-state-cells", type=int, default=30)
    p.add_argument("--high-threshold", type=float, default=2.6)
    p.add_argument("--low-threshold", type=float, default=1.6)
    p.add_argument("--min-run-windows", type=int, default=3)
    p.add_argument("--min-interval-genes", type=int, default=60)
    p.add_argument("--min-mean-score", type=float, default=1.8)
    p.add_argument("--burden-quantile", type=float, default=0.95)
    p.add_argument("--cnv-burden-threshold", type=float, default=1.8)
    p.add_argument("--write-h5ad", action="store_true")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args(argv)


def params_from_args(args: argparse.Namespace) -> FastParams:
    gene_coords = Path(args.gene_coordinates) if args.gene_coordinates else None
    if gene_coords is None and args.species:
        gene_coords = bundled_gene_coordinates(args.species)
    return FastParams(
        h5ad=Path(args.h5ad),
        control_h5ad=Path(args.control_h5ad),
        gene_coordinates=gene_coords,
        output_prefix=Path(args.output),
        state_key=args.state_key,
        sample_key=args.sample_key,
        control_state_key=args.control_state_key,
        layer=args.layer,
        input_normalized=bool(args.input_normalized),
        window_genes=int(args.window_genes),
        stride_genes=int(args.stride_genes),
        min_chr_genes=int(args.min_chr_genes),
        min_state_cells=int(args.min_state_cells),
        high_threshold=float(args.high_threshold),
        low_threshold=float(args.low_threshold),
        min_run_windows=int(args.min_run_windows),
        min_interval_genes=int(args.min_interval_genes),
        min_mean_score=float(args.min_mean_score),
        burden_quantile=float(args.burden_quantile),
        cnv_burden_threshold=float(args.cnv_burden_threshold),
        write_h5ad=bool(args.write_h5ad),
    )


def main(argv: Optional[Sequence[str]] = None) -> Dict[str, Path]:
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s | %(message)s")
    outputs = run_fast(params_from_args(args))
    for name, path in outputs.items():
        LOGGER.info("%s: %s", name, path)
    return outputs


if __name__ == "__main__":
    main()
