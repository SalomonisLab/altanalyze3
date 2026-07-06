"""Per-gene expression-weighted candidate fastCNV caller.

This module is deliberately isolated from the current production fastCNV code.
The design follows the acceptance notes in `/Users/saljh8/Dropbox/Revio/fastCNV`:

* compare every query cell to a cell-state-matched healthy reference centroid;
* use per-gene log2 fold-change in CP10k space, not an equal-weight mean of
  log-normalized expression;
* ignore genes that are silent or nearly silent in that reference state;
* aggregate along genomic windows using reference-expression weights;
* emit both cell-level regional calls and sample/state-level global summaries.

The thresholds here are general fold-change/fraction thresholds, not region- or
LOY-specific parameters. They are intentionally exposed so benchmark gates can
be run before any production integration.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import time
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from altanalyze3.components.fastCNV.main import (
    Window,
    build_windows,
    bundled_gene_coordinates,
    load_gene_coordinates,
)


@dataclass(frozen=True)
class CandidateParams:
    h5ad: Path
    control_h5ad: Path
    output_prefix: Path
    state_key: str
    control_state_key: str
    gene_coordinates: Optional[Path] = None
    sample_key: Optional[str] = None
    control_sample_key: Optional[str] = None
    chromosomes: Optional[Tuple[str, ...]] = None
    control_filter_field: Optional[str] = None
    control_filter_value: Optional[str] = None
    obs_names_file: Optional[Path] = None
    layer: str = "X"
    control_layer: str = "X"
    input_normalized: bool = True
    control_input_normalized: bool = True
    window_genes: int = 101
    stride_genes: int = 25
    min_chr_genes: int = 50
    min_ref_cells: int = 30
    max_ref_cells_per_state: int = 600
    chunk_cells: int = 512
    min_ref_detect_frac: float = 0.05
    min_ref_cp10k: float = 0.05
    min_weighted_genes: float = 20.0
    min_regulated_weight_frac: float = 0.45
    min_abs_mean_log2: float = 0.45
    loss_log2: float = -1.0
    gain_log2: float = 1.0
    amp_log2: float = 2.0
    ratio_pseudocount: float = 1.0
    empirical_alpha: float = 0.001
    empirical_tolerance: float = 0.05
    min_copy_z: float = 4.0
    copy_margin_z: float = 3.0
    min_null_sd: float = 0.08
    haploid_chromosomes: Tuple[str, ...] = ("chrY",)
    mixture_calibration: bool = True
    mixture_max_windows: int = 5
    mixture_min_cells: int = 80
    mixture_min_separation_z: float = 2.5
    mixture_min_component_frac: float = 0.01
    near_zero_cp10k: float = 0.03
    center_per_cell: bool = True
    min_run_windows: int = 3
    refine_intervals: bool = True
    min_interval_support_genes: int = 5
    max_cells: Optional[int] = None
    random_state: int = 0


@dataclass
class ReferenceState:
    state: str
    n_cells: int
    rows: np.ndarray
    mean_cp10k: np.ndarray
    detect_frac: np.ndarray
    weights: np.ndarray


@dataclass
class ChromNull:
    null_center: np.ndarray
    null_sd: np.ndarray
    low_cut: np.ndarray
    high_cut: np.ndarray
    expected: Dict[int, np.ndarray]


@dataclass
class SparseChromBuffer:
    mean_log2: np.ndarray
    loss_frac: np.ndarray
    zero_frac: np.ndarray
    gain_frac: np.ndarray
    amp_frac: np.ndarray
    confidence: np.ndarray


def _format_seconds(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    return f"{seconds / 60:.1f}m"


def _matrix_for_adata(adata: ad.AnnData, layer: str):
    if layer == "X":
        if adata.X is None:
            raise ValueError("AnnData X is empty; specify a populated layer with --layer.")
        return adata.X
    if layer == "auto":
        return adata.layers["counts"] if "counts" in adata.layers else adata.X
    if layer not in adata.layers:
        raise KeyError(f"Layer {layer!r} is not present in AnnData.")
    return adata.layers[layer]


def _resolved_layer_name(adata: ad.AnnData, layer: str) -> str:
    if layer == "auto":
        return "counts" if "counts" in adata.layers else "X"
    return layer


def _row_sums(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        return np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32)
    if hasattr(matrix, "group") and getattr(matrix, "format", None) == "csr":
        indptr = np.asarray(matrix.indptr, dtype=np.int64)
        data_ds = matrix.group["data"]
        sums = np.zeros(matrix.shape[0], dtype=np.float32)
        block_rows = 4096
        for start_row in range(0, matrix.shape[0], block_rows):
            end_row = min(start_row + block_rows, matrix.shape[0])
            data_start = int(indptr[start_row])
            data_end = int(indptr[end_row])
            if data_end <= data_start:
                continue
            block_data = np.asarray(data_ds[data_start:data_end], dtype=np.float32)
            counts = np.diff(indptr[start_row:end_row + 1])
            offsets = np.r_[0, np.cumsum(counts[:-1], dtype=np.int64)]
            sums[start_row:end_row] = np.add.reduceat(block_data, offsets)
        return sums
    return np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32)


def _slice_to_dense(matrix, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    chunk = matrix[rows][:, cols]
    if sp.issparse(chunk):
        chunk = chunk.toarray()
    return np.asarray(chunk, dtype=np.float32)


def _cp10k_chunk(
    matrix,
    rows: np.ndarray,
    cols: np.ndarray,
    library_sizes: Optional[np.ndarray],
    *,
    input_normalized: bool,
) -> np.ndarray:
    chunk = _slice_to_dense(matrix, rows, cols)
    if input_normalized:
        np.expm1(chunk, out=chunk)
        np.maximum(chunk, 0.0, out=chunk)
        return chunk
    if library_sizes is None:
        raise ValueError("library_sizes are required when input_normalized=False")
    scale = 10000.0 / np.maximum(library_sizes[rows].astype(np.float32), 1.0)
    return chunk * scale[:, None]


def _ordered_gene_table(adata: ad.AnnData, coord_path: Path) -> pd.DataFrame:
    coords = load_gene_coordinates(coord_path, adata.var_names)
    return coords.reset_index(drop=True)


def _normalize_chrom(chrom: str) -> str:
    value = str(chrom).strip()
    if not value:
        return value
    return value if value.startswith("chr") else f"chr{value}"


def _parse_csv_values(value: Optional[str]) -> Optional[Tuple[str, ...]]:
    if not value:
        return None
    parsed = tuple(item.strip() for item in value.split(",") if item.strip())
    return parsed or None


def _row_subset_from_obs_names(adata: ad.AnnData, obs_names_file: Optional[Path]) -> Optional[np.ndarray]:
    if obs_names_file is None:
        return None
    names = pd.read_csv(obs_names_file, sep="\t")
    if "CellBarcode" in names.columns:
        wanted = names["CellBarcode"].astype(str)
    else:
        wanted = names.iloc[:, 0].astype(str)
    lookup = pd.Series(np.arange(adata.n_obs, dtype=np.int64), index=adata.obs_names.astype(str))
    rows = lookup.reindex(wanted).dropna().astype(np.int64).to_numpy()
    if rows.size == 0:
        raise ValueError(f"No obs names from {obs_names_file} matched {adata.filename}.")
    return np.sort(rows)


def _backed_csr_subset(matrix, rows: Optional[np.ndarray], cols: np.ndarray):
    if not hasattr(matrix, "group") or getattr(matrix, "format", None) != "csr":
        raise TypeError(f"Unsupported backed sparse matrix type: {type(matrix)!r}")
    cols = np.asarray(cols, dtype=np.int64)
    rows = np.arange(matrix.shape[0], dtype=np.int64) if rows is None else np.asarray(rows, dtype=np.int64)
    col_map = np.full(matrix.shape[1], -1, dtype=np.int64)
    col_map[cols] = np.arange(cols.size, dtype=np.int64)
    group = matrix.group
    data_ds = group["data"]
    indices_ds = group["indices"]
    indptr = np.asarray(matrix.indptr, dtype=np.int64)
    all_columns = cols.size == matrix.shape[1] and np.array_equal(cols, np.arange(matrix.shape[1], dtype=np.int64))
    if all_columns:
        out_data: List[np.ndarray] = []
        out_indices: List[np.ndarray] = []
        out_indptr = np.zeros(rows.size + 1, dtype=np.int64)
        nnz = 0
        start_pos = 0
        while start_pos < rows.size:
            end_pos = start_pos + 1
            while end_pos < rows.size and rows[end_pos] == rows[end_pos - 1] + 1:
                end_pos += 1
            row_block = rows[start_pos:end_pos]
            data_start = int(indptr[row_block[0]])
            data_end = int(indptr[row_block[-1] + 1])
            block_data = np.asarray(data_ds[data_start:data_end])
            block_indices = np.asarray(indices_ds[data_start:data_end], dtype=np.int32)
            block_counts = np.diff(indptr[row_block])
            out_data.append(block_data)
            out_indices.append(block_indices)
            out_indptr[start_pos + 1:end_pos + 1] = nnz + np.cumsum(block_counts, dtype=np.int64)
            nnz += int(block_counts.sum())
            start_pos = end_pos
        data = np.concatenate(out_data) if out_data else np.array([], dtype=matrix.dtype)
        indices = np.concatenate(out_indices) if out_indices else np.array([], dtype=np.int32)
        return sp.csr_matrix((data, indices, out_indptr.astype(np.int32, copy=False)), shape=(rows.size, cols.size))
    out_data: List[np.ndarray] = []
    out_indices: List[np.ndarray] = []
    out_indptr = np.zeros(rows.size + 1, dtype=np.int64)
    nnz = 0
    start_pos = 0
    while start_pos < rows.size:
        end_pos = start_pos + 1
        while end_pos < rows.size and rows[end_pos] == rows[end_pos - 1] + 1:
            end_pos += 1
        row_block = rows[start_pos:end_pos]
        data_start = int(indptr[row_block[0]])
        data_end = int(indptr[row_block[-1] + 1])
        block_data = np.asarray(data_ds[data_start:data_end])
        block_indices = np.asarray(indices_ds[data_start:data_end], dtype=np.int64)
        for offset, row in enumerate(row_block):
            rel_start = int(indptr[row] - data_start)
            rel_end = int(indptr[row + 1] - data_start)
            if rel_end > rel_start:
                row_indices = block_indices[rel_start:rel_end]
                mapped = col_map[row_indices]
                keep = mapped >= 0
                if np.any(keep):
                    out_data.append(block_data[rel_start:rel_end][keep])
                    out_indices.append(mapped[keep].astype(np.int32, copy=False))
                    nnz += int(keep.sum())
            out_indptr[start_pos + offset + 1] = nnz
        start_pos = end_pos
    data = np.concatenate(out_data) if out_data else np.array([], dtype=matrix.dtype)
    indices = np.concatenate(out_indices) if out_indices else np.array([], dtype=np.int32)
    return sp.csr_matrix((data, indices, out_indptr.astype(np.int32, copy=False)), shape=(rows.size, cols.size))


def _subset_matrix(matrix, rows: Optional[np.ndarray], cols: np.ndarray):
    cols = np.asarray(cols, dtype=np.int64)
    if sp.issparse(matrix):
        return matrix[:, cols].copy() if rows is None else matrix[rows][:, cols].copy()
    if hasattr(matrix, "group"):
        return _backed_csr_subset(matrix, rows, cols)
    subset = matrix[:, cols] if rows is None else matrix[rows][:, cols]
    if hasattr(subset, "to_memory"):
        subset = subset.to_memory()
    return np.asarray(subset)


def _materialize_gene_subset(
    adata: ad.AnnData,
    genes: Sequence[str],
    layer: str,
    rows: Optional[np.ndarray] = None,
) -> ad.AnnData:
    """Load only the required genes into memory.

    Backed AnnData slicing materializes populated layers without reading every
    gene. This keeps chromosome-restricted benchmarks fast and avoids loading
    the full 80k x 36k LOY query matrix when only chrY is being scored.
    """
    gene_index = pd.Index(adata.var_names.astype(str))
    cols = gene_index.get_indexer(pd.Index(list(genes), dtype=str))
    if np.any(cols < 0):
        missing = [gene for gene, idx in zip(genes, cols) if idx < 0][:5]
        raise ValueError(f"Requested genes are absent from AnnData: {missing}")
    actual_layer = _resolved_layer_name(adata, layer)
    matrix = _matrix_for_adata(adata, layer)
    subset_matrix = _subset_matrix(matrix, rows, cols)
    obs = adata.obs.copy() if rows is None else adata.obs.iloc[rows].copy()
    subset = ad.AnnData(
        X=subset_matrix if actual_layer == "X" else None,
        obs=obs,
        var=adata.var.iloc[cols].copy(),
    )
    if actual_layer != "X":
        subset.layers[actual_layer] = subset_matrix
    return subset


def _state_values(adata: ad.AnnData, key: str) -> pd.Series:
    if key not in adata.obs.columns:
        raise KeyError(f"State key {key!r} is not present in obs.")
    return adata.obs[key].astype(str).str.strip().replace({"nan": "", "None": ""})


def _sample_values(adata: ad.AnnData, key: Optional[str]) -> pd.Series:
    if key and key in adata.obs.columns:
        return adata.obs[key].astype(str).str.strip().replace({"nan": "", "None": ""})
    return pd.Series(["sample"] * adata.n_obs, index=adata.obs_names, dtype=str)


def _choose_reference_rows(
    state_series: pd.Series,
    state: str,
    eligible_mask: np.ndarray,
    *,
    min_ref_cells: int,
    max_ref_cells: int,
    rng: np.random.Generator,
) -> np.ndarray:
    state_values = state_series.to_numpy(dtype=str)
    rows = np.flatnonzero((state_values == str(state)) & eligible_mask)
    if rows.size < min_ref_cells:
        rows = np.flatnonzero(eligible_mask)
    if rows.size > max_ref_cells > 0:
        rows = np.sort(rng.choice(rows, size=max_ref_cells, replace=False))
    return rows.astype(np.int64)


def build_reference_states(
    control: ad.AnnData,
    coords: pd.DataFrame,
    params: CandidateParams,
    library_sizes: Optional[np.ndarray] = None,
) -> Dict[str, ReferenceState]:
    matrix = _matrix_for_adata(control, params.control_layer)
    if not params.control_input_normalized and library_sizes is None:
        library_sizes = _row_sums(matrix)
    states = _state_values(control, params.control_state_key)
    cols = coords["var_index"].to_numpy(dtype=np.int64)
    rng = np.random.default_rng(params.random_state)
    eligible = np.ones(control.n_obs, dtype=bool)
    if params.control_filter_field or params.control_filter_value:
        if not params.control_filter_field or params.control_filter_field not in control.obs.columns:
            raise KeyError(f"Control filter field {params.control_filter_field!r} is not present in obs.")
        wanted = str(params.control_filter_value or "").strip()
        eligible &= control.obs[params.control_filter_field].astype(str).str.strip().to_numpy() == wanted
        if not np.any(eligible):
            raise ValueError(
                f"Control filter {params.control_filter_field}={wanted!r} selected zero cells."
            )

    ref: Dict[str, ReferenceState] = {}
    for state in sorted(value for value in pd.unique(states) if value):
        rows = _choose_reference_rows(
            states,
            state,
            eligible,
            min_ref_cells=params.min_ref_cells,
            max_ref_cells=params.max_ref_cells_per_state,
            rng=rng,
        )
        cp = _cp10k_chunk(
            matrix,
            rows,
            cols,
            library_sizes,
            input_normalized=params.control_input_normalized,
        )
        mean_cp = cp.mean(axis=0, dtype=np.float64).astype(np.float32)
        detect = (cp > 0).mean(axis=0, dtype=np.float64).astype(np.float32)
        weights = np.log1p(mean_cp).astype(np.float32)
        weights[(detect < params.min_ref_detect_frac) | (mean_cp < params.min_ref_cp10k)] = 0.0
        ref[state] = ReferenceState(
            state=state,
            n_cells=int(rows.size),
            rows=rows,
            mean_cp10k=mean_cp,
            detect_frac=detect,
            weights=weights,
        )
    return ref


def _windows_by_chrom(windows: Sequence[Window]) -> Dict[str, List[Window]]:
    by_chrom: Dict[str, List[Window]] = {}
    for window in windows:
        by_chrom.setdefault(window.chrom, []).append(window)
    return by_chrom


def _prefix_sum(values: np.ndarray) -> np.ndarray:
    return np.concatenate([[0.0], np.cumsum(values, dtype=np.float64)])


def _window_sums(prefix: np.ndarray, windows: Sequence[Window]) -> np.ndarray:
    starts = np.fromiter((w.start_offset for w in windows), dtype=np.int64, count=len(windows))
    ends = np.fromiter((w.end_offset for w in windows), dtype=np.int64, count=len(windows))
    return prefix[ends] - prefix[starts]


def _weighted_window_mean(values: np.ndarray, weights: np.ndarray, windows: Sequence[Window]) -> np.ndarray:
    weighted_prefix = np.concatenate(
        [
            np.zeros((values.shape[0], 1), dtype=np.float64),
            np.cumsum(values * weights[None, :], axis=1, dtype=np.float64),
        ],
        axis=1,
    )
    starts = np.fromiter((w.start_offset for w in windows), dtype=np.int64, count=len(windows))
    ends = np.fromiter((w.end_offset for w in windows), dtype=np.int64, count=len(windows))
    weight_prefix = _prefix_sum(weights)
    wsum = _window_sums(weight_prefix, windows)
    return (weighted_prefix[:, ends] - weighted_prefix[:, starts]) / np.maximum(wsum[None, :], 1e-6)


def _normal_copy_number(chrom: str, params: CandidateParams) -> int:
    return 1 if _normalize_chrom(chrom) in set(params.haploid_chromosomes) else 2


def _copy_state_expectations(
    ref_cp: np.ndarray,
    weights: np.ndarray,
    windows: Sequence[Window],
    chrom: str,
    params: CandidateParams,
) -> Dict[int, np.ndarray]:
    normal = _normal_copy_number(chrom, params)
    copy_states = (0, 1, 2, 4) if normal == 1 else (0, 1, 2, 3, 4)
    expected: Dict[int, np.ndarray] = {}
    for copies in copy_states:
        factor = copies / max(float(normal), 1.0)
        gene_expected = np.log2(
            (ref_cp * factor + params.ratio_pseudocount)
            / (ref_cp + params.ratio_pseudocount)
        ).astype(np.float32)
        prefix = _prefix_sum(gene_expected * weights)
        weight_prefix = _prefix_sum(weights)
        wsum = _window_sums(weight_prefix, windows)
        expected[copies] = _window_sums(prefix, windows) / np.maximum(wsum, 1e-6)
    return expected


def _mad(values: np.ndarray, axis: int = 0) -> np.ndarray:
    med = np.nanmedian(values, axis=axis, keepdims=True)
    return np.nanmedian(np.abs(values - med), axis=axis)


def _build_chrom_nulls(
    control_matrix,
    control_library_sizes: Optional[np.ndarray],
    coords: pd.DataFrame,
    chrom_windows: Dict[str, List[Window]],
    chrom_positions: Dict[str, np.ndarray],
    ref_state: ReferenceState,
    params: CandidateParams,
) -> Dict[str, ChromNull]:
    ref_cp = np.maximum(ref_state.mean_cp10k, 0.0)
    weights = ref_state.weights.astype(np.float32)
    valid = weights > 0
    ref_cells = _cp10k_chunk(
        control_matrix,
        ref_state.rows,
        coords["var_index"].to_numpy(dtype=np.int64),
        control_library_sizes,
        input_normalized=params.control_input_normalized,
    )
    ratios = np.log2(
        (ref_cells + params.ratio_pseudocount) / (ref_cp[None, :] + params.ratio_pseudocount)
    ).astype(np.float32)
    ratios[:, ~valid] = 0.0

    autosome_mask = valid & coords["chr"].astype(str).str.match(r"^chr([1-9]|1[0-9]|2[0-2])$").to_numpy()
    if params.center_per_cell and np.any(autosome_mask):
        center_weights = weights[autosome_mask]
        centers = (
            (ratios[:, autosome_mask] * center_weights[None, :]).sum(axis=1)
            / max(float(center_weights.sum()), 1e-6)
        ).astype(np.float32)
        ratios[:, valid] = ratios[:, valid] - centers[:, None]

    nulls: Dict[str, ChromNull] = {}
    for chrom, windows in chrom_windows.items():
        positions = chrom_positions[chrom]
        chrom_weights = weights[positions]
        if float(chrom_weights.sum()) <= 0:
            continue
        ref_win = _weighted_window_mean(ratios[:, positions], chrom_weights, windows)
        center = np.nanmedian(ref_win, axis=0)
        sd = np.maximum(1.4826 * _mad(ref_win, axis=0), params.min_null_sd)
        low_cut = np.nanquantile(ref_win, params.empirical_alpha, axis=0)
        high_cut = np.nanquantile(ref_win, 1.0 - params.empirical_alpha, axis=0)
        expected = _copy_state_expectations(ref_cp[positions], chrom_weights, windows, chrom, params)
        nulls[chrom] = ChromNull(
            null_center=center.astype(np.float32),
            null_sd=sd.astype(np.float32),
            low_cut=low_cut.astype(np.float32),
            high_cut=high_cut.astype(np.float32),
            expected={key: value.astype(np.float32) for key, value in expected.items()},
        )
    return nulls


def _call_code(mean_log2: np.ndarray, loss_frac: np.ndarray, zero_frac: np.ndarray, gain_frac: np.ndarray, amp_frac: np.ndarray, params: CandidateParams) -> np.ndarray:
    calls = np.zeros(mean_log2.shape, dtype=np.int8)
    loss_like = np.maximum(loss_frac, zero_frac)
    calls[(loss_like >= params.min_regulated_weight_frac) & (mean_log2 <= -params.min_abs_mean_log2)] = -1
    calls[(zero_frac >= params.min_regulated_weight_frac) & (mean_log2 <= -params.min_abs_mean_log2)] = -2
    calls[(gain_frac >= params.min_regulated_weight_frac) & (mean_log2 >= params.min_abs_mean_log2)] = 1
    calls[(amp_frac >= params.min_regulated_weight_frac) & (mean_log2 >= max(params.min_abs_mean_log2, 1.0))] = 2
    return calls


def _copy_state_calls(
    mean_log2: np.ndarray,
    null: ChromNull,
    chrom: str,
    loss_frac: np.ndarray,
    zero_frac: np.ndarray,
    gain_frac: np.ndarray,
    amp_frac: np.ndarray,
    params: CandidateParams,
) -> Tuple[np.ndarray, np.ndarray]:
    normal = _normal_copy_number(chrom, params)
    z_normal = np.abs((mean_log2 - null.expected[normal][None, :]) / null.null_sd[None, :])
    best_dist = z_normal.copy()
    best_copy = np.full(mean_log2.shape, normal, dtype=np.int16)
    for copies, expected in null.expected.items():
        if copies == normal:
            continue
        dist = np.abs((mean_log2 - expected[None, :]) / null.null_sd[None, :])
        better = dist < best_dist
        best_dist[better] = dist[better]
        best_copy[better] = int(copies)

    calls = np.zeros(mean_log2.shape, dtype=np.int8)
    confident = (z_normal >= params.min_copy_z) & ((z_normal - best_dist) >= params.copy_margin_z)
    empirical_loss = mean_log2 <= (null.low_cut[None, :] + params.empirical_tolerance)
    empirical_gain = mean_log2 >= (null.high_cut[None, :] - params.empirical_tolerance)
    loss_like = np.maximum(loss_frac, zero_frac)
    gain_like = np.maximum(gain_frac, amp_frac)
    min_frac = max(params.min_regulated_weight_frac * 0.5, 0.05)

    loss = confident & empirical_loss & (best_copy < normal) & (loss_like >= min_frac)
    gain = confident & empirical_gain & (best_copy > normal) & (gain_like >= min_frac)
    calls[loss] = -1
    calls[gain] = 1
    calls[loss & (best_copy == 0)] = -2
    calls[gain & (best_copy >= normal * 2)] = 2
    confidence = z_normal - best_dist
    return calls, confidence.astype(np.float32)


def _two_means(values: np.ndarray) -> Optional[Tuple[np.ndarray, np.ndarray, float]]:
    values = np.asarray(values, dtype=np.float32)
    valid = np.isfinite(values)
    if int(valid.sum()) < 10:
        return None
    x = values[valid]
    if float(np.nanstd(x)) <= 1e-6:
        return None
    centers = np.array([np.nanquantile(x, 0.25), np.nanquantile(x, 0.75)], dtype=np.float32)
    labels = np.zeros(x.shape[0], dtype=np.int8)
    for _ in range(25):
        distances = np.abs(x[:, None] - centers[None, :])
        new_labels = np.argmin(distances, axis=1).astype(np.int8)
        if np.array_equal(new_labels, labels):
            break
        labels = new_labels
        for idx in (0, 1):
            if np.any(labels == idx):
                centers[idx] = float(np.nanmean(x[labels == idx]))
    order = np.argsort(centers)
    centers = centers[order]
    relabel = np.zeros_like(labels)
    relabel[labels == order[0]] = 0
    relabel[labels == order[1]] = 1
    full_labels = np.full(values.shape[0], -1, dtype=np.int8)
    full_labels[valid] = relabel
    pooled_sd = max(float(np.nanstd(x)), 1e-6)
    separation = float((centers[1] - centers[0]) / pooled_sd)
    return full_labels, centers.astype(np.float32), separation


def _mixture_sparse_calls(
    buffer: SparseChromBuffer,
    chrom: str,
    params: CandidateParams,
) -> Tuple[np.ndarray, str]:
    score = np.nanmean(buffer.mean_log2, axis=1)
    fit = _two_means(score)
    if fit is None:
        return np.zeros(score.shape[0], dtype=np.int8), ""
    labels, centers, separation = fit
    if separation < params.mixture_min_separation_z:
        return np.zeros(score.shape[0], dtype=np.int8), ""
    calls = np.zeros(score.shape[0], dtype=np.int8)
    n_valid = int(np.count_nonzero(labels >= 0))
    min_component = max(int(np.ceil(n_valid * params.mixture_min_component_frac)), 3)
    loss_like = np.maximum(np.nanmean(buffer.loss_frac, axis=1), np.nanmean(buffer.zero_frac, axis=1))
    gain_like = np.maximum(np.nanmean(buffer.gain_frac, axis=1), np.nanmean(buffer.amp_frac, axis=1))
    min_frac = max(params.min_regulated_weight_frac * 0.5, 0.05)

    low = labels == 0
    high = labels == 1
    if int(low.sum()) >= min_component and centers[0] <= -params.min_abs_mean_log2:
        calls[low & (loss_like >= min_frac)] = -1
    if int(high.sum()) >= min_component and centers[1] >= params.min_abs_mean_log2:
        calls[high & (gain_like >= min_frac)] = 1
    label = f"mixture_sep={separation:.3f};centers={centers[0]:.3f},{centers[1]:.3f}"
    return calls, label


def _merge_cell_runs(
    cell: str,
    state: str,
    sample: str,
    chrom_windows: Sequence[Window],
    calls: np.ndarray,
    mean_log2: np.ndarray,
    loss_frac: np.ndarray,
    zero_frac: np.ndarray,
    gain_frac: np.ndarray,
    amp_frac: np.ndarray,
    confidence: np.ndarray,
    gene_ratios: np.ndarray,
    gene_cp10k: np.ndarray,
    chrom_weights: np.ndarray,
    chrom_gene_table: pd.DataFrame,
    params: CandidateParams,
) -> List[dict]:
    rows: List[dict] = []
    start = 0
    while start < len(chrom_windows):
        code = int(calls[start])
        if code == 0:
            start += 1
            continue
        end = start + 1
        while end < len(chrom_windows) and int(calls[end]) == code:
            end += 1
        min_run_windows = min(params.min_run_windows, max(len(chrom_windows), 1))
        if end - start >= min_run_windows:
            block = chrom_windows[start:end]
            label = {-2: "homozygous_loss", -1: "loss", 1: "gain", 2: "amplification"}[code]
            refined = _refine_interval_bounds(
                block,
                code,
                gene_ratios,
                gene_cp10k,
                chrom_weights,
                chrom_gene_table,
                params,
            )
            rows.append(
                {
                    "CellBarcode": cell,
                    "cell_state": state,
                    "sample": sample,
                    "call": label,
                    "chr": block[0].chrom,
                    "start": int(refined["start"]),
                    "end": int(refined["end"]),
                    "start_gene": refined["start_gene"],
                    "end_gene": refined["end_gene"],
                    "n_windows": int(end - start),
                    "n_genes": int(sum(w.n_genes for w in block)),
                    "n_informative_genes": int(refined["n_informative_genes"]),
                    "n_evidence_genes": int(refined["n_evidence_genes"]),
                    "evidence_gene_fraction": float(refined["evidence_gene_fraction"]),
                    "boundary_source": refined["boundary_source"],
                    "mean_log2": float(np.nanmean(mean_log2[start:end])),
                    "max_abs_log2": float(np.nanmax(np.abs(mean_log2[start:end]))),
                    "loss_weight_frac": float(np.nanmean(loss_frac[start:end])),
                    "zero_weight_frac": float(np.nanmean(zero_frac[start:end])),
                    "gain_weight_frac": float(np.nanmean(gain_frac[start:end])),
                    "amp_weight_frac": float(np.nanmean(amp_frac[start:end])),
                    "mean_confidence_z": float(np.nanmean(confidence[start:end])),
                    "max_confidence_z": float(np.nanmax(confidence[start:end])),
                    "caller": "copy_state",
                }
            )
        start = end
    return rows


def _refine_interval_bounds(
    block: Sequence[Window],
    code: int,
    gene_ratios: np.ndarray,
    gene_cp10k: np.ndarray,
    chrom_weights: np.ndarray,
    chrom_gene_table: pd.DataFrame,
    params: CandidateParams,
) -> dict:
    """Refine a called window run to the altered gene-supported span.

    Detection is still done on robust sliding windows for speed/noise control,
    but reported CNV intervals should not default to chromosome/arm spans when
    only a subset of expressed genes supports the event. This refinement uses
    only per-gene evidence inside the called run and falls back to the window
    span when there are too few informative altered genes.
    """
    run_start = min(window.start_offset for window in block)
    run_end = max(window.end_offset for window in block)
    gene_frame = chrom_gene_table.iloc[run_start:run_end].reset_index(drop=True)
    ratios = np.asarray(gene_ratios[run_start:run_end], dtype=np.float32)
    cp10k = np.asarray(gene_cp10k[run_start:run_end], dtype=np.float32)
    weights = np.asarray(chrom_weights[run_start:run_end], dtype=np.float32)
    informative = np.isfinite(ratios) & (weights > 0)
    if code < 0:
        evidence = informative & ((ratios <= params.loss_log2) | (cp10k <= params.near_zero_cp10k))
    elif code > 0:
        threshold = params.amp_log2 if code == 2 else params.gain_log2
        evidence = informative & (ratios >= threshold)
    else:
        evidence = np.zeros(ratios.shape, dtype=bool)

    n_info = int(np.count_nonzero(informative))
    n_evidence = int(np.count_nonzero(evidence))
    evidence_fraction = n_evidence / max(n_info, 1)
    if params.refine_intervals and n_evidence >= params.min_interval_support_genes:
        offsets = np.flatnonzero(evidence)
        start_idx = int(offsets[0])
        end_idx = int(offsets[-1])
        start_row = gene_frame.iloc[start_idx]
        end_row = gene_frame.iloc[end_idx]
        return {
            "start": int(min(start_row["start"], start_row["end"])),
            "end": int(max(end_row["start"], end_row["end"])),
            "start_gene": str(start_row["gene"]),
            "end_gene": str(end_row["gene"]),
            "n_informative_genes": n_info,
            "n_evidence_genes": n_evidence,
            "evidence_gene_fraction": evidence_fraction,
            "boundary_source": "gene_evidence",
        }

    return {
        "start": int(min(window.start for window in block)),
        "end": int(max(window.end for window in block)),
        "start_gene": block[0].start_gene,
        "end_gene": block[-1].end_gene,
        "n_informative_genes": n_info,
        "n_evidence_genes": n_evidence,
        "evidence_gene_fraction": evidence_fraction,
        "boundary_source": "window_span",
    }


def call_candidate(params: CandidateParams) -> dict:
    started = time.perf_counter()
    query = ad.read_h5ad(params.h5ad, backed="r")
    control = ad.read_h5ad(params.control_h5ad, backed="r")
    try:
        coord_path = params.gene_coordinates or bundled_gene_coordinates("human")
        coords = _ordered_gene_table(query, Path(coord_path))
        control_var = pd.Index(control.var_names.astype(str))
        coords = coords.loc[coords["gene"].astype(str).isin(control_var)].reset_index(drop=True)
        if params.chromosomes:
            wanted_chroms = {_normalize_chrom(chrom) for chrom in params.chromosomes}
            coords = coords.loc[coords["chr"].astype(str).map(_normalize_chrom).isin(wanted_chroms)].reset_index(drop=True)
        if coords.empty:
            raise ValueError("No query genes with coordinates were present in the control AnnData.")
        query_full_library_sizes = None
        control_full_library_sizes = None
        if not params.input_normalized:
            query_full_library_sizes = _row_sums(_matrix_for_adata(query, params.layer))
        if not params.control_input_normalized:
            control_full_library_sizes = _row_sums(_matrix_for_adata(control, params.control_layer))
        query_rows = _row_subset_from_obs_names(query, params.obs_names_file)
        if params.max_cells and params.max_cells > 0 and params.max_cells < query.n_obs:
            rng = np.random.default_rng(params.random_state)
            base_rows = np.arange(query.n_obs, dtype=np.int64) if query_rows is None else query_rows
            query_rows = np.sort(rng.choice(base_rows, size=min(params.max_cells, base_rows.size), replace=False))
        if query_rows is not None and query_full_library_sizes is not None:
            query_full_library_sizes = query_full_library_sizes[query_rows]
        genes = coords["gene"].astype(str).tolist()
        query_memory = _materialize_gene_subset(query, genes, params.layer, query_rows)
        control_memory = _materialize_gene_subset(control, genes, params.control_layer)
        if getattr(query, "isbacked", False):
            query.file.close()
        if getattr(control, "isbacked", False):
            control.file.close()
        query = query_memory
        control = control_memory
        coords["var_index"] = np.arange(len(genes), dtype=np.int64)
        coords["control_var_index"] = np.arange(len(genes), dtype=np.int64)
        windows = build_windows(coords, params.window_genes, params.stride_genes, params.min_chr_genes)
        chrom_windows = _windows_by_chrom(windows)
        chrom_positions = {
            str(chrom): np.flatnonzero(coords["chr"].astype(str).to_numpy() == str(chrom)).astype(np.int64)
            for chrom in chrom_windows
        }

        control_coords = coords.copy()
        control_coords["var_index"] = control_coords["control_var_index"]
        reference = build_reference_states(control, control_coords, params, control_full_library_sizes)
        control_matrix = _matrix_for_adata(control, params.control_layer)
        control_library_sizes = None if params.control_input_normalized else control_full_library_sizes

        q_matrix = _matrix_for_adata(query, params.layer)
        q_library_sizes = None if params.input_normalized else query_full_library_sizes
        q_states = _state_values(query, params.state_key)
        q_samples = _sample_values(query, params.sample_key)
        all_rows = np.arange(query.n_obs, dtype=np.int64)

        interval_rows: List[dict] = []
        cell_rows: List[dict] = []
        for state in sorted(pd.unique(q_states.iloc[all_rows])):
            state = str(state)
            if not state:
                continue
            ref_state = reference.get(state)
            if ref_state is None:
                ref_state = reference.get(next(iter(reference)))
            state_rows = all_rows[q_states.iloc[all_rows].to_numpy(dtype=str) == state]
            if state_rows.size == 0:
                continue
            ref_cp = np.maximum(ref_state.mean_cp10k, 0.0)
            weights = ref_state.weights.astype(np.float32)
            valid = weights > 0
            autosome_mask = valid & coords["chr"].astype(str).str.match(r"^chr([1-9]|1[0-9]|2[0-2])$").to_numpy()
            chrom_nulls = _build_chrom_nulls(
                control_matrix,
                control_library_sizes,
                control_coords,
                chrom_windows,
                chrom_positions,
                ref_state,
                params,
            )
            state_interval_count = np.zeros(state_rows.size, dtype=np.int32)
            row_to_local = {int(row): idx for idx, row in enumerate(state_rows)}
            state_called_chroms: List[set] = [set() for _ in range(state_rows.size)]
            sparse_buffers: Dict[str, SparseChromBuffer] = {}
            if params.mixture_calibration and state_rows.size >= params.mixture_min_cells:
                for chrom, chr_windows in chrom_windows.items():
                    if len(chr_windows) <= params.mixture_max_windows:
                        shape = (state_rows.size, len(chr_windows))
                        sparse_buffers[chrom] = SparseChromBuffer(
                            mean_log2=np.full(shape, np.nan, dtype=np.float32),
                            loss_frac=np.full(shape, np.nan, dtype=np.float32),
                            zero_frac=np.full(shape, np.nan, dtype=np.float32),
                            gain_frac=np.full(shape, np.nan, dtype=np.float32),
                            amp_frac=np.full(shape, np.nan, dtype=np.float32),
                            confidence=np.full(shape, np.nan, dtype=np.float32),
                        )

            for chunk_start in range(0, state_rows.size, params.chunk_cells):
                rows = state_rows[chunk_start:chunk_start + params.chunk_cells]
                cp = _cp10k_chunk(
                    q_matrix,
                    rows,
                    coords["var_index"].to_numpy(dtype=np.int64),
                    q_library_sizes,
                    input_normalized=params.input_normalized,
                )
                ratios = np.log2((cp + params.ratio_pseudocount) / (ref_cp[None, :] + params.ratio_pseudocount)).astype(np.float32)
                ratios[:, ~valid] = 0.0
                if params.center_per_cell and np.any(autosome_mask):
                    center_weights = weights[autosome_mask]
                    centers = (
                        (ratios[:, autosome_mask] * center_weights[None, :]).sum(axis=1)
                        / max(float(center_weights.sum()), 1e-6)
                    ).astype(np.float32)
                    ratios[:, valid] = ratios[:, valid] - centers[:, None]

                for chrom, chr_windows in chrom_windows.items():
                    if not chr_windows:
                        continue
                    positions = chrom_positions[chrom]
                    chrom_weights = weights[positions]
                    if float(chrom_weights.sum()) <= 0:
                        continue
                    local_windows = list(chr_windows)
                    weight_prefix = _prefix_sum(chrom_weights)
                    wsum = _window_sums(weight_prefix, local_windows)
                    usable_threshold = min(params.min_weighted_genes, max(float(chrom_weights.sum()) * 0.8, 1.0))
                    usable = wsum >= usable_threshold
                    if not np.any(usable):
                        continue

                    r = ratios[:, positions]
                    c = cp[:, positions]
                    weighted_ratio_prefix = np.concatenate(
                        [np.zeros((r.shape[0], 1), dtype=np.float64), np.cumsum(r * chrom_weights[None, :], axis=1, dtype=np.float64)],
                        axis=1,
                    )
                    starts = np.fromiter((w.start_offset for w in local_windows), dtype=np.int64, count=len(local_windows))
                    ends = np.fromiter((w.end_offset for w in local_windows), dtype=np.int64, count=len(local_windows))
                    mean_log2 = (weighted_ratio_prefix[:, ends] - weighted_ratio_prefix[:, starts]) / np.maximum(wsum[None, :], 1e-6)
                    loss_prefix = np.concatenate(
                        [np.zeros((r.shape[0], 1), dtype=np.float64), np.cumsum((r <= params.loss_log2) * chrom_weights[None, :], axis=1, dtype=np.float64)],
                        axis=1,
                    )
                    zero_prefix = np.concatenate(
                        [np.zeros((r.shape[0], 1), dtype=np.float64), np.cumsum((c <= params.near_zero_cp10k) * chrom_weights[None, :], axis=1, dtype=np.float64)],
                        axis=1,
                    )
                    gain_prefix = np.concatenate(
                        [np.zeros((r.shape[0], 1), dtype=np.float64), np.cumsum((r >= params.gain_log2) * chrom_weights[None, :], axis=1, dtype=np.float64)],
                        axis=1,
                    )
                    amp_prefix = np.concatenate(
                        [np.zeros((r.shape[0], 1), dtype=np.float64), np.cumsum((r >= params.amp_log2) * chrom_weights[None, :], axis=1, dtype=np.float64)],
                        axis=1,
                    )
                    loss_frac = (loss_prefix[:, ends] - loss_prefix[:, starts]) / np.maximum(wsum[None, :], 1e-6)
                    zero_frac = (zero_prefix[:, ends] - zero_prefix[:, starts]) / np.maximum(wsum[None, :], 1e-6)
                    gain_frac = (gain_prefix[:, ends] - gain_prefix[:, starts]) / np.maximum(wsum[None, :], 1e-6)
                    amp_frac = (amp_prefix[:, ends] - amp_prefix[:, starts]) / np.maximum(wsum[None, :], 1e-6)
                    null = chrom_nulls.get(chrom)
                    if null is not None:
                        calls, confidence = _copy_state_calls(
                            mean_log2,
                            null,
                            chrom,
                            loss_frac,
                            zero_frac,
                            gain_frac,
                            amp_frac,
                            params,
                        )
                    else:
                        calls = _call_code(mean_log2, loss_frac, zero_frac, gain_frac, amp_frac, params)
                        confidence = np.abs(mean_log2).astype(np.float32)
                    calls[:, ~usable] = 0
                    confidence[:, ~usable] = 0.0
                    buffer = sparse_buffers.get(chrom)
                    if buffer is not None:
                        local_indices = np.fromiter(
                            (row_to_local[int(row_idx)] for row_idx in rows),
                            dtype=np.int64,
                            count=rows.size,
                        )
                        buffer.mean_log2[local_indices, :] = mean_log2.astype(np.float32)
                        buffer.loss_frac[local_indices, :] = loss_frac.astype(np.float32)
                        buffer.zero_frac[local_indices, :] = zero_frac.astype(np.float32)
                        buffer.gain_frac[local_indices, :] = gain_frac.astype(np.float32)
                        buffer.amp_frac[local_indices, :] = amp_frac.astype(np.float32)
                        buffer.confidence[local_indices, :] = confidence.astype(np.float32)
                    chrom_gene_table = coords.iloc[positions].reset_index(drop=True)
                    for i, row_idx in enumerate(rows):
                        local_row = row_to_local[int(row_idx)]
                        merged = _merge_cell_runs(
                            str(query.obs_names[row_idx]),
                            state,
                            str(q_samples.iloc[row_idx] or "sample"),
                            chr_windows,
                            calls[i],
                            mean_log2[i],
                            loss_frac[i],
                            zero_frac[i],
                            gain_frac[i],
                            amp_frac[i],
                            confidence[i],
                            r[i],
                            c[i],
                            chrom_weights,
                            chrom_gene_table,
                            params,
                        )
                        if merged:
                            state_interval_count[local_row] += len(merged)
                            for item in merged:
                                state_called_chroms[local_row].add(str(item["chr"]))
                            interval_rows.extend(merged)

            for chrom, buffer in sparse_buffers.items():
                mix_calls, mix_label = _mixture_sparse_calls(buffer, chrom, params)
                if not np.any(mix_calls):
                    continue
                chr_windows = chrom_windows[chrom]
                for local_idx, code in enumerate(mix_calls):
                    code = int(code)
                    if code == 0 or chrom in state_called_chroms[local_idx]:
                        continue
                    block = chr_windows
                    label = {-1: "loss", 1: "gain"}[code]
                    row_idx = int(state_rows[local_idx])
                    row = {
                        "CellBarcode": str(query.obs_names[row_idx]),
                        "cell_state": state,
                        "sample": str(q_samples.iloc[row_idx] or "sample"),
                        "call": label,
                        "chr": block[0].chrom,
                        "start": int(min(w.start for w in block)),
                        "end": int(max(w.end for w in block)),
                        "start_gene": block[0].start_gene,
                        "end_gene": block[-1].end_gene,
                        "n_windows": int(len(block)),
                        "n_genes": int(sum(w.n_genes for w in block)),
                        "n_informative_genes": int(sum(w.n_genes for w in block)),
                        "n_evidence_genes": np.nan,
                        "evidence_gene_fraction": np.nan,
                        "boundary_source": "sparse_chromosome_mixture",
                        "mean_log2": float(np.nanmean(buffer.mean_log2[local_idx])),
                        "max_abs_log2": float(np.nanmax(np.abs(buffer.mean_log2[local_idx]))),
                        "loss_weight_frac": float(np.nanmean(buffer.loss_frac[local_idx])),
                        "zero_weight_frac": float(np.nanmean(buffer.zero_frac[local_idx])),
                        "gain_weight_frac": float(np.nanmean(buffer.gain_frac[local_idx])),
                        "amp_weight_frac": float(np.nanmean(buffer.amp_frac[local_idx])),
                        "mean_confidence_z": float(np.nanmean(buffer.confidence[local_idx])),
                        "max_confidence_z": float(np.nanmax(buffer.confidence[local_idx])),
                        "caller": mix_label,
                    }
                    state_interval_count[local_idx] += 1
                    state_called_chroms[local_idx].add(chrom)
                    interval_rows.append(row)

            for row_idx, n_intervals in zip(state_rows, state_interval_count):
                cell_rows.append(
                    {
                        "CellBarcode": str(query.obs_names[row_idx]),
                        "cell_state": state,
                        "sample": str(q_samples.iloc[row_idx] or "sample"),
                        "cnv_status": "CNV" if n_intervals else "WT",
                        "n_cnv_intervals": int(n_intervals),
                    }
                )

        intervals = pd.DataFrame(interval_rows)
        cells = pd.DataFrame(cell_rows)
        if intervals.empty:
            intervals = pd.DataFrame(
                columns=[
                    "CellBarcode", "cell_state", "sample", "call", "chr", "start", "end",
                    "start_gene", "end_gene", "n_windows", "n_genes",
                    "n_informative_genes", "n_evidence_genes", "evidence_gene_fraction",
                    "boundary_source", "mean_log2",
                    "max_abs_log2", "loss_weight_frac", "zero_weight_frac",
                    "gain_weight_frac", "amp_weight_frac",
                    "mean_confidence_z", "max_confidence_z",
                    "caller",
                ]
            )
        global_summary = _summarize_global(intervals, cells)
        params.output_prefix.parent.mkdir(parents=True, exist_ok=True)
        cells.to_csv(f"{params.output_prefix}.candidate_cnv_cells.tsv", sep="\t", index=False)
        intervals.to_csv(f"{params.output_prefix}.candidate_cnv_regions.tsv", sep="\t", index=False)
        global_summary.to_csv(f"{params.output_prefix}.candidate_global_summary.tsv", sep="\t", index=False)
        reference_summary = pd.DataFrame(
            [
                {
                    "state": item.state,
                    "n_reference_cells": item.n_cells,
                    "n_weighted_genes": int(np.count_nonzero(item.weights > 0)),
                    "median_weight": float(np.median(item.weights[item.weights > 0])) if np.any(item.weights > 0) else 0.0,
                }
                for item in reference.values()
            ]
        )
        reference_summary.to_csv(f"{params.output_prefix}.candidate_reference_summary.tsv", sep="\t", index=False)
        return {
            "n_cells": int(cells.shape[0]),
            "n_cnv_cells": int((cells["cnv_status"] == "CNV").sum()) if not cells.empty else 0,
            "n_regions": int(intervals.shape[0]),
            "n_global_rows": int(global_summary.shape[0]),
            "elapsed_seconds": round(time.perf_counter() - started, 3),
            "output_prefix": str(params.output_prefix),
        }
    finally:
        if getattr(query, "isbacked", False):
            query.file.close()
        if getattr(control, "isbacked", False):
            control.file.close()


def _summarize_global(intervals: pd.DataFrame, cells: pd.DataFrame) -> pd.DataFrame:
    if intervals.empty:
        return pd.DataFrame(columns=["sample", "cell_state", "call", "chr", "region", "n_cells", "fraction_cells", "mean_log2"])
    total = cells.groupby(["sample", "cell_state"], observed=False).size().rename("n_total_cells").reset_index()
    frame = intervals.copy()
    frame["region"] = frame["chr"].astype(str) + ":" + frame["start"].astype(str) + "-" + frame["end"].astype(str)
    grouped = (
        frame.groupby(["sample", "cell_state", "call", "chr", "region"], observed=False)
        .agg(n_cells=("CellBarcode", "nunique"), mean_log2=("mean_log2", "mean"))
        .reset_index()
    )
    grouped = grouped.merge(total, on=["sample", "cell_state"], how="left")
    grouped["fraction_cells"] = grouped["n_cells"] / grouped["n_total_cells"].replace(0, np.nan)
    return grouped.sort_values(["fraction_cells", "n_cells"], ascending=[False, False]).reset_index(drop=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Experimental per-gene fastCNV candidate caller.")
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--control-h5ad", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--state-key", required=True)
    parser.add_argument("--control-state-key", required=True)
    parser.add_argument("--gene-coordinates", type=Path)
    parser.add_argument("--sample-key")
    parser.add_argument("--control-sample-key")
    parser.add_argument("--chromosomes", help="Comma-separated chromosome subset, e.g. chrY or chr5,chr7,chr8.")
    parser.add_argument("--control-filter-field", help="Optional obs column used to restrict control cells.")
    parser.add_argument("--control-filter-value", help="Value required in --control-filter-field.")
    parser.add_argument("--obs-names-file", type=Path, help="Optional TSV/list of query CellBarcode values to evaluate.")
    parser.add_argument("--layer", default="X")
    parser.add_argument("--control-layer", default="X")
    parser.add_argument("--counts-input", action="store_true", help="Input h5ad X/layer is counts, not log1p(CP10k).")
    parser.add_argument("--counts-control", action="store_true", help="Control h5ad X/layer is counts, not log1p(CP10k).")
    parser.add_argument("--window-genes", type=int, default=101)
    parser.add_argument("--stride-genes", type=int, default=25)
    parser.add_argument("--min-chr-genes", type=int, default=50)
    parser.add_argument("--chunk-cells", type=int, default=512)
    parser.add_argument("--min-ref-detect-frac", type=float, default=0.05)
    parser.add_argument("--min-ref-cp10k", type=float, default=0.05)
    parser.add_argument("--min-weighted-genes", type=float, default=20.0)
    parser.add_argument("--min-regulated-weight-frac", type=float, default=0.45)
    parser.add_argument("--min-abs-mean-log2", type=float, default=0.45)
    parser.add_argument("--empirical-alpha", type=float, default=0.001)
    parser.add_argument("--empirical-tolerance", type=float, default=0.05)
    parser.add_argument("--min-copy-z", type=float, default=4.0)
    parser.add_argument("--copy-margin-z", type=float, default=3.0)
    parser.add_argument("--min-null-sd", type=float, default=0.08)
    parser.add_argument("--haploid-chromosomes", default="chrY")
    parser.add_argument("--no-mixture-calibration", action="store_true")
    parser.add_argument("--mixture-max-windows", type=int, default=5)
    parser.add_argument("--mixture-min-cells", type=int, default=80)
    parser.add_argument("--mixture-min-separation-z", type=float, default=2.5)
    parser.add_argument("--mixture-min-component-frac", type=float, default=0.01)
    parser.add_argument("--min-run-windows", type=int, default=3)
    parser.add_argument("--no-refine-intervals", action="store_true")
    parser.add_argument("--min-interval-support-genes", type=int, default=5)
    parser.add_argument("--max-cells", type=int)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    params = CandidateParams(
        h5ad=args.h5ad,
        control_h5ad=args.control_h5ad,
        output_prefix=args.output_prefix,
        state_key=args.state_key,
        control_state_key=args.control_state_key,
        gene_coordinates=args.gene_coordinates,
        sample_key=args.sample_key,
        control_sample_key=args.control_sample_key,
        chromosomes=tuple(_normalize_chrom(chrom) for chrom in _parse_csv_values(args.chromosomes) or ()),
        control_filter_field=args.control_filter_field,
        control_filter_value=args.control_filter_value,
        obs_names_file=args.obs_names_file,
        layer=args.layer,
        control_layer=args.control_layer,
        input_normalized=not args.counts_input,
        control_input_normalized=not args.counts_control,
        window_genes=args.window_genes,
        stride_genes=args.stride_genes,
        min_chr_genes=args.min_chr_genes,
        chunk_cells=args.chunk_cells,
        min_ref_detect_frac=args.min_ref_detect_frac,
        min_ref_cp10k=args.min_ref_cp10k,
        min_weighted_genes=args.min_weighted_genes,
        min_regulated_weight_frac=args.min_regulated_weight_frac,
        min_abs_mean_log2=args.min_abs_mean_log2,
        empirical_alpha=args.empirical_alpha,
        empirical_tolerance=args.empirical_tolerance,
        min_copy_z=args.min_copy_z,
        copy_margin_z=args.copy_margin_z,
        min_null_sd=args.min_null_sd,
        haploid_chromosomes=tuple(_normalize_chrom(chrom) for chrom in _parse_csv_values(args.haploid_chromosomes) or ()),
        mixture_calibration=not args.no_mixture_calibration,
        mixture_max_windows=args.mixture_max_windows,
        mixture_min_cells=args.mixture_min_cells,
        mixture_min_separation_z=args.mixture_min_separation_z,
        mixture_min_component_frac=args.mixture_min_component_frac,
        min_run_windows=args.min_run_windows,
        refine_intervals=not args.no_refine_intervals,
        min_interval_support_genes=args.min_interval_support_genes,
        max_cells=args.max_cells,
    )
    summary = call_candidate(params)
    print(pd.Series(summary).to_string())
    print(f"elapsed={_format_seconds(float(summary['elapsed_seconds']))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
