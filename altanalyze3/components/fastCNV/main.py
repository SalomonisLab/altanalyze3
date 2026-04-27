"""Fast state-aware first-pass RNA CNV caller.

The base fastCNV caller operates on a cellHarmony-produced AnnData object. It
builds Ensembl-coordinate adaptive gene windows, learns WT-like anchors within
each assigned cell state, and emits cell-level CNV intervals plus a continuous
window score matrix for downstream sparse NMF clone discovery.
"""

from __future__ import annotations

import argparse
import logging
import math
import time
from importlib import resources
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


LOGGER = logging.getLogger(__name__)
RESOURCE_FILES = {
    "human": "Hs_Ensembl_GRCh38_genes.tsv",
    "mouse": "Mm_Ensembl_GRCm39_genes.tsv",
}

_GENE_COLUMNS = ("gene", "gene_id", "gene_name", "symbol", "Gene", "GeneID")
_CHR_COLUMNS = ("chr", "chromosome", "seqname", "chrom", "Chr")
_START_COLUMNS = ("start", "tx_start", "gene_start", "Start")
_END_COLUMNS = ("end", "tx_end", "gene_end", "End")
_MAIN_CHR_ORDER = {str(i): i for i in range(1, 23)}
_MAIN_CHR_ORDER.update({"X": 23, "Y": 24, "M": 25, "MT": 25})


@dataclass
class FastCNVParams:
    """Parameters controlling first-pass fastCNV cell-level calling."""

    h5ad: Path
    gene_coordinates: Optional[Path]
    output_prefix: Path
    state_key: str
    sample_key: Optional[str] = None
    layer: str = "auto"
    input_normalized: bool = False
    control_h5ad: Optional[Path] = None
    control_state_key: Optional[str] = None
    control_sample_key: Optional[str] = None
    min_control_cells: int = 30
    window_genes: int = 41
    stride_genes: int = 7
    min_chr_genes: int = 25
    min_state_cells: int = 30
    anchor_fraction: float = 0.25
    min_anchor_cells: int = 20
    high_threshold: float = 2.6
    low_threshold: float = 1.6
    min_run_windows: int = 3
    min_interval_genes: int = 60
    min_mean_score: float = 1.8
    burden_quantile: float = 0.95
    cnv_burden_threshold: float = 1.8
    max_clones_per_state: int = 10
    max_global_clones: int = 10
    min_clone_cells: int = 5
    clone_similarity_threshold: float = 0.88
    clone_consensus_fraction: float = 0.45
    nmf_max_iter: int = 300
    skip_clones: bool = False
    skip_pdf: bool = False
    write_h5ad: bool = False
    max_cells_per_state: Optional[int] = None
    random_state: int = 0


@dataclass
class Window:
    window_id: str
    chrom: str
    start: int
    end: int
    start_gene: str
    end_gene: str
    start_offset: int
    end_offset: int
    n_genes: int


def bundled_gene_coordinates(species: str) -> Path:
    """Return the bundled fastCNV gene-coordinate table for a species."""
    if species not in RESOURCE_FILES:
        raise ValueError(f"Unsupported species '{species}'. Expected one of: {', '.join(sorted(RESOURCE_FILES))}")
    resource = resources.files("altanalyze3.components.fastCNV").joinpath("resources", RESOURCE_FILES[species])
    path = Path(str(resource))
    if not path.exists():
        raise FileNotFoundError(f"Bundled fastCNV coordinate resource is missing: {path}")
    return path


def _first_existing(columns: Iterable[str], candidates: Sequence[str], label: str) -> str:
    available = set(columns)
    for candidate in candidates:
        if candidate in available:
            return candidate
    raise ValueError(f"Gene coordinate file is missing a {label} column. Tried: {', '.join(candidates)}")


def _chr_sort_key(chrom: object) -> Tuple[int, str]:
    text = str(chrom).strip()
    clean = text[3:] if text.lower().startswith("chr") else text
    clean_upper = clean.upper()
    return (_MAIN_CHR_ORDER.get(clean_upper, 1000), text)


def _normalize_chr(chrom: object) -> str:
    text = str(chrom).strip()
    if not text:
        return text
    return text if text.lower().startswith("chr") else f"chr{text}"


def load_gene_coordinates(path: Path, var_names: pd.Index) -> pd.DataFrame:
    """Load and align an Ensembl-style gene coordinate table to AnnData genes."""
    coords = pd.read_csv(path, sep=None, engine="python")
    gene_col = _first_existing(coords.columns, _GENE_COLUMNS, "gene")
    chr_col = _first_existing(coords.columns, _CHR_COLUMNS, "chromosome")
    start_col = _first_existing(coords.columns, _START_COLUMNS, "start")
    end_col = _first_existing(coords.columns, _END_COLUMNS, "end")

    coords = coords[[gene_col, chr_col, start_col, end_col]].copy()
    coords.columns = ["gene", "chr", "start", "end"]
    coords["gene"] = coords["gene"].astype(str)
    coords["chr"] = coords["chr"].map(_normalize_chr)
    coords["start"] = pd.to_numeric(coords["start"], errors="coerce")
    coords["end"] = pd.to_numeric(coords["end"], errors="coerce")
    coords = coords.dropna(subset=["gene", "chr", "start", "end"])
    coords["start"] = coords["start"].astype(np.int64)
    coords["end"] = coords["end"].astype(np.int64)
    coords["midpoint"] = ((coords["start"] + coords["end"]) // 2).astype(np.int64)

    var_lookup = pd.Series(np.arange(len(var_names), dtype=np.int64), index=var_names.astype(str))
    coords = coords[coords["gene"].isin(var_lookup.index)].copy()
    coords["var_index"] = coords["gene"].map(var_lookup).astype(np.int64)
    coords = coords.drop_duplicates(subset=["gene"], keep="first")
    coords["_chr_rank"] = coords["chr"].map(lambda value: _chr_sort_key(value)[0])
    coords["_chr_text"] = coords["chr"].astype(str)
    coords = coords.sort_values(["_chr_rank", "_chr_text", "midpoint", "gene"]).reset_index(drop=True)
    coords = coords.drop(columns=["_chr_rank", "_chr_text"])
    if coords.empty:
        raise ValueError("No AnnData var_names matched the provided gene coordinate file.")
    return coords


def build_windows(coords: pd.DataFrame, window_genes: int, stride_genes: int, min_chr_genes: int) -> List[Window]:
    """Build overlapping adaptive gene windows from ordered gene coordinates."""
    if window_genes < 3:
        raise ValueError("--window-genes must be at least 3.")
    if stride_genes < 1:
        raise ValueError("--stride-genes must be at least 1.")

    windows: List[Window] = []
    for chrom, chr_df in coords.groupby("chr", sort=False):
        chr_df = chr_df.reset_index(drop=True)
        n_genes = len(chr_df)
        if n_genes < min_chr_genes:
            continue
        size = min(window_genes, n_genes)
        starts = list(range(0, max(n_genes - size + 1, 1), stride_genes))
        final_start = max(n_genes - size, 0)
        if starts[-1] != final_start:
            starts.append(final_start)
        for start_offset in starts:
            end_offset = min(start_offset + size, n_genes)
            block = chr_df.iloc[start_offset:end_offset]
            window_id = f"{chrom}:{int(block['start'].min())}-{int(block['end'].max())}:{len(windows)}"
            windows.append(
                Window(
                    window_id=window_id,
                    chrom=str(chrom),
                    start=int(block["start"].min()),
                    end=int(block["end"].max()),
                    start_gene=str(block.iloc[0]["gene"]),
                    end_gene=str(block.iloc[-1]["gene"]),
                    start_offset=int(start_offset),
                    end_offset=int(end_offset),
                    n_genes=int(end_offset - start_offset),
                )
            )
    if not windows:
        raise ValueError("No genome windows were built. Check gene coordinates and minimum chromosome gene settings.")
    return windows


def _matrix_for_adata(adata: ad.AnnData, layer: str) -> sp.spmatrix | np.ndarray:
    if layer == "auto":
        return adata.layers["counts"] if "counts" in adata.layers else adata.X
    if layer == "X":
        return adata.X
    if layer not in adata.layers:
        raise KeyError(f"Layer '{layer}' not found in AnnData.")
    return adata.layers[layer]


def _row_sums(matrix: sp.spmatrix | np.ndarray) -> np.ndarray:
    if sp.issparse(matrix):
        return np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32)
    return np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32)


def _extract_normalized_chunk(
    matrix: sp.spmatrix | np.ndarray,
    row_indices: np.ndarray,
    col_indices: np.ndarray,
    library_sizes: np.ndarray,
    input_normalized: bool,
) -> np.ndarray:
    chunk = matrix[row_indices][:, col_indices]
    if sp.issparse(chunk):
        chunk = chunk.toarray()
    else:
        chunk = np.asarray(chunk)
    chunk = chunk.astype(np.float32, copy=False)
    if input_normalized:
        return chunk
    factors = 10000.0 / np.maximum(library_sizes[row_indices].astype(np.float32), 1.0)
    chunk = chunk * factors[:, None]
    np.log1p(chunk, out=chunk)
    return chunk


def _rolling_window_mean(values: np.ndarray, windows: Sequence[Window]) -> np.ndarray:
    if values.shape[1] == 0 or not windows:
        return np.zeros((values.shape[0], 0), dtype=np.float32)
    prefix = np.concatenate(
        [np.zeros((values.shape[0], 1), dtype=np.float32), np.cumsum(values, axis=1, dtype=np.float32)],
        axis=1,
    )
    starts = np.asarray([window.start_offset for window in windows], dtype=np.int64)
    ends = np.asarray([window.end_offset for window in windows], dtype=np.int64)
    lengths = np.maximum(ends - starts, 1).astype(np.float32)
    return ((prefix[:, ends] - prefix[:, starts]) / lengths).astype(np.float32, copy=False)


def _mad(values: np.ndarray, axis: int = 0) -> np.ndarray:
    med = np.nanmedian(values, axis=axis)
    return np.nanmedian(np.abs(values - np.expand_dims(med, axis)), axis=axis)


@dataclass
class ControlReference:
    """External control AnnData aligned to the query gene space."""

    matrix: sp.spmatrix | np.ndarray
    library_sizes: np.ndarray
    obs: pd.DataFrame
    var_index_map: np.ndarray
    input_normalized: bool


def load_control_reference(
    path: Path,
    query_var_names: pd.Index,
    input_normalized: bool,
    layer: str,
    state_key: Optional[str],
    sample_key: Optional[str],
) -> ControlReference:
    """Load a control AnnData and align its var space to the query genes."""
    LOGGER.info("Reading control AnnData: %s", path)
    control = ad.read_h5ad(path)
    available_control_obs = ", ".join(control.obs.columns.astype(str)) or "<none>"
    if state_key and state_key not in control.obs.columns:
        raise KeyError(
            f"Control state column '{state_key}' not found in control AnnData obs. Available obs columns: {available_control_obs}"
        )
    if sample_key and sample_key not in control.obs.columns:
        raise KeyError(
            f"Control sample column '{sample_key}' not found in control AnnData obs. Available obs columns: {available_control_obs}"
        )

    query_set = set(query_var_names.astype(str))
    control_var = control.var_names.astype(str)
    shared_mask = control_var.isin(query_set)
    n_shared = int(shared_mask.sum())
    overlap_fraction = n_shared / max(len(query_var_names), 1)
    if n_shared == 0:
        raise ValueError("Control AnnData shares no genes with the query AnnData.")
    if overlap_fraction < 0.50:
        raise ValueError(
            f"Control AnnData overlaps only {overlap_fraction:.1%} of query genes "
            f"({n_shared}/{len(query_var_names)}); refusing to use as baseline."
        )
    if overlap_fraction < 0.95:
        LOGGER.warning(
            "Control AnnData covers %.1f%% of query genes (%d/%d). Missing genes will be ignored in baseline.",
            overlap_fraction * 100.0,
            n_shared,
            len(query_var_names),
        )

    matrix = _matrix_for_adata(control, layer)
    if sp.issparse(matrix):
        matrix = matrix.tocsr()
    library_sizes = np.ones(control.n_obs, dtype=np.float32) if input_normalized else _row_sums(matrix)

    control_lookup = pd.Series(np.arange(control.n_vars, dtype=np.int64), index=control_var)
    var_index_map = control_lookup.reindex(query_var_names.astype(str)).to_numpy()
    return ControlReference(
        matrix=matrix,
        library_sizes=library_sizes,
        obs=control.obs.copy(),
        var_index_map=var_index_map,
        input_normalized=input_normalized,
    )


def _control_rows_for_state(
    control: ControlReference,
    state: str,
    sample: str,
    state_key: Optional[str],
    sample_key: Optional[str],
) -> np.ndarray:
    rows = np.arange(control.obs.shape[0], dtype=np.int64)
    if state_key and state_key in control.obs.columns:
        state_match = control.obs[state_key].astype(str).to_numpy() == state
        if state_match.any():
            rows = rows[state_match]
        else:
            LOGGER.warning("No control cells match state '%s'; falling back to all controls.", state)
            rows = np.arange(control.obs.shape[0], dtype=np.int64)
    if sample_key and sample and sample_key in control.obs.columns:
        sample_match = control.obs.iloc[rows][sample_key].astype(str).to_numpy() == sample
        if sample_match.any():
            rows = rows[sample_match]
    return rows


def _control_chunk(
    control: ControlReference,
    control_rows: np.ndarray,
    query_col_indices: np.ndarray,
) -> np.ndarray:
    """Return a normalized control expression chunk aligned to query columns."""
    mapped = control.var_index_map[query_col_indices]
    valid_mask = ~pd.isna(mapped)
    valid_indices = mapped[valid_mask].astype(np.int64)
    chunk = np.zeros((control_rows.size, query_col_indices.size), dtype=np.float32)
    if valid_indices.size == 0:
        return chunk
    sub = control.matrix[control_rows][:, valid_indices]
    if sp.issparse(sub):
        sub = sub.toarray()
    sub = np.asarray(sub, dtype=np.float32)
    if not control.input_normalized:
        factors = 10000.0 / np.maximum(control.library_sizes[control_rows].astype(np.float32), 1.0)
        sub = sub * factors[:, None]
        np.log1p(sub, out=sub)
    chunk[:, valid_mask] = sub
    return chunk


def _window_scores_for_state(
    matrix: sp.spmatrix | np.ndarray,
    coords: pd.DataFrame,
    windows_by_chr: Dict[str, List[Window]],
    state_rows: np.ndarray,
    library_sizes: np.ndarray,
    input_normalized: bool,
    anchor_mask: Optional[np.ndarray],
    control: Optional[ControlReference] = None,
    control_rows: Optional[np.ndarray] = None,
) -> np.ndarray:
    use_control = control is not None and control_rows is not None and control_rows.size > 0
    state_scores: List[np.ndarray] = []
    for chrom, chr_windows in windows_by_chr.items():
        chr_coords = coords[coords["chr"] == chrom]
        col_indices = chr_coords["var_index"].to_numpy(dtype=np.int64)
        expr = _extract_normalized_chunk(matrix, state_rows, col_indices, library_sizes, input_normalized)
        if use_control:
            control_expr = _control_chunk(control, control_rows, col_indices)
            baseline = np.nanmedian(control_expr, axis=0).astype(np.float32)
            residual = expr - baseline[None, :]
            raw_window_scores = _rolling_window_mean(residual, chr_windows)
            control_residual = control_expr - baseline[None, :]
            control_window_scores = _rolling_window_mean(control_residual, chr_windows)
            center = np.nanmedian(control_window_scores, axis=0).astype(np.float32)
            scale = (1.4826 * _mad(control_window_scores, axis=0)).astype(np.float32)
        else:
            baseline_source = expr if anchor_mask is None else expr[anchor_mask]
            baseline = np.nanmedian(baseline_source, axis=0).astype(np.float32)
            residual = expr - baseline[None, :]
            raw_window_scores = _rolling_window_mean(residual, chr_windows)
            scale_source = raw_window_scores if anchor_mask is None else raw_window_scores[anchor_mask]
            center = np.nanmedian(scale_source, axis=0).astype(np.float32)
            scale = (1.4826 * _mad(scale_source, axis=0)).astype(np.float32)
        scale = np.where(scale < 0.10, 0.10, scale).astype(np.float32)
        scores = ((raw_window_scores - center[None, :]) / scale[None, :]).astype(np.float32)
        state_scores.append(scores)
    return np.concatenate(state_scores, axis=1) if state_scores else np.zeros((len(state_rows), 0), dtype=np.float32)


def _cnv_burden(scores: np.ndarray, quantile: float) -> np.ndarray:
    if scores.shape[1] == 0:
        return np.zeros(scores.shape[0], dtype=np.float32)
    quantile = min(max(float(quantile), 0.50), 0.99)
    k = max(1, int(math.ceil((1.0 - quantile) * scores.shape[1])))
    abs_scores = np.abs(scores)
    top = np.partition(abs_scores, scores.shape[1] - k, axis=1)[:, -k:]
    return np.nanmean(top, axis=1).astype(np.float32)


def _iter_runs(mask: np.ndarray) -> Iterable[Tuple[int, int]]:
    padded = np.concatenate([[False], mask.astype(bool), [False]])
    changes = np.flatnonzero(padded[1:] != padded[:-1])
    for start, end in zip(changes[0::2], changes[1::2]):
        yield int(start), int(end)


def _call_intervals_for_cell(
    scores: np.ndarray,
    windows: Sequence[Window],
    chr_slices: Dict[str, slice],
    *,
    high_threshold: float,
    low_threshold: float,
    min_run_windows: int,
    min_interval_genes: int,
    min_mean_score: float,
) -> List[Dict[str, object]]:
    calls: List[Dict[str, object]] = []
    for chrom, chr_slice in chr_slices.items():
        chr_scores = scores[chr_slice]
        chr_windows = windows[chr_slice]
        for direction, sign in (("gain", 1.0), ("loss", -1.0)):
            signed = chr_scores * sign
            low_mask = signed >= low_threshold
            high_mask = signed >= high_threshold
            for start, end in _iter_runs(low_mask):
                if end - start < min_run_windows or not np.any(high_mask[start:end]):
                    continue
                run_scores = chr_scores[start:end]
                run_windows = chr_windows[start:end]
                n_genes = int(run_windows[-1].end_offset - run_windows[0].start_offset)
                mean_signed = float(np.nanmean(run_scores * sign))
                if n_genes < min_interval_genes or mean_signed < min_mean_score:
                    continue
                max_signed = float(np.nanmax(run_scores * sign))
                calls.append(
                    {
                        "call": direction,
                        "chr": chrom,
                        "start": int(run_windows[0].start),
                        "end": int(run_windows[-1].end),
                        "start_gene": run_windows[0].start_gene,
                        "end_gene": run_windows[-1].end_gene,
                        "n_windows": int(end - start),
                        "n_genes": n_genes,
                        "mean_score": mean_signed,
                        "max_score": max_signed,
                        "confidence": mean_signed * math.sqrt(max(end - start, 1)),
                    }
                )
    return calls


def _window_frame(windows: Sequence[Window]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "window_id": [window.window_id for window in windows],
            "chr": [window.chrom for window in windows],
            "start": [window.start for window in windows],
            "end": [window.end for window in windows],
            "start_gene": [window.start_gene for window in windows],
            "end_gene": [window.end_gene for window in windows],
            "n_genes": [window.n_genes for window in windows],
        }
    )


def _chr_slices(windows: Sequence[Window]) -> Dict[str, slice]:
    slices: Dict[str, slice] = {}
    start = 0
    while start < len(windows):
        chrom = windows[start].chrom
        end = start + 1
        while end < len(windows) and windows[end].chrom == chrom:
            end += 1
        slices[chrom] = slice(start, end)
        start = end
    return slices


def _state_indices(obs: pd.DataFrame, state_key: str, state: str) -> np.ndarray:
    return np.flatnonzero(obs[state_key].astype(str).to_numpy() == state).astype(np.int64)


def _subsample_indices(indices: np.ndarray, max_cells: Optional[int], seed: int) -> np.ndarray:
    if max_cells is None or len(indices) <= max_cells:
        return indices
    rng = np.random.default_rng(seed)
    selected = rng.choice(indices, size=max_cells, replace=False)
    return np.sort(selected.astype(np.int64))


def _cnv_features(scores: np.ndarray) -> np.ndarray:
    finite = np.nan_to_num(scores.astype(np.float32, copy=False), nan=0.0, posinf=0.0, neginf=0.0)
    clipped = np.clip(finite, -8.0, 8.0)
    return np.concatenate([np.maximum(clipped, 0.0), np.maximum(-clipped, 0.0)], axis=1).astype(np.float32)


def _estimate_nmf_rank(features: np.ndarray, max_rank: int) -> int:
    if features.shape[0] < 2 or features.shape[1] < 2:
        return 1
    centered = features - np.nanmean(features, axis=0, keepdims=True)
    try:
        singular_values = np.linalg.svd(centered, compute_uv=False)
    except np.linalg.LinAlgError:
        return 1
    if singular_values.size == 0 or float(np.sum(singular_values**2)) <= 0:
        return 1
    variance = (singular_values**2) / np.sum(singular_values**2)
    informative = int(np.sum(variance >= 0.05))
    cumulative = int(np.searchsorted(np.cumsum(variance), 0.90) + 1)
    rank = max(1, min(informative if informative > 0 else cumulative, cumulative, max_rank))
    return int(min(rank, features.shape[0], features.shape[1], max_rank))


def _discover_state_clones(
    scores: np.ndarray,
    cnv_mask: np.ndarray,
    anchor_mask: np.ndarray,
    params: FastCNVParams,
    seed: int,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Assign WT and state-local clone IDs from the continuous CNV score matrix."""
    labels = np.full(scores.shape[0], "WT", dtype=object)
    candidates = np.flatnonzero(cnv_mask)
    if params.skip_clones or candidates.size < params.min_clone_cells:
        return labels, pd.DataFrame()

    features = _cnv_features(scores[candidates])
    useful = np.nanmax(features, axis=0) > 0.25
    features = features[:, useful]
    if features.shape[1] == 0:
        return labels, pd.DataFrame()

    max_rank = max(1, min(params.max_clones_per_state, candidates.size // max(params.min_clone_cells, 1)))
    rank = _estimate_nmf_rank(features, max_rank)
    if rank <= 1:
        component_labels = np.zeros(candidates.size, dtype=int)
        component_profiles = np.nanmean(scores[candidates], axis=0, keepdims=True)
    else:
        from sklearn.decomposition import NMF

        model = NMF(
            n_components=rank,
            init="nndsvda",
            random_state=seed,
            max_iter=params.nmf_max_iter,
            solver="cd",
            beta_loss="frobenius",
        )
        weights = model.fit_transform(np.maximum(features, 0.0))
        component_labels = np.argmax(weights, axis=1).astype(int)
        component_profiles = np.vstack(
            [
                np.nanmean(scores[candidates[component_labels == component]], axis=0)
                if np.any(component_labels == component)
                else np.zeros(scores.shape[1], dtype=np.float32)
                for component in range(rank)
            ]
        )

    component_rows = []
    component_order = sorted(
        set(component_labels.tolist()),
        key=lambda component: int(np.sum(component_labels == component)),
        reverse=True,
    )
    clone_number = 1
    for component in component_order:
        local = np.flatnonzero(component_labels == component)
        support = int(local.size)
        if support < params.min_clone_cells:
            continue
        clone_id = f"clone{clone_number}"
        clone_number += 1
        labels[candidates[local]] = clone_id
        profile = component_profiles[component]
        component_rows.append(
            {
                "state_clone_id": clone_id,
                "n_cells": support,
                "mean_burden": float(np.nanmean(_cnv_burden(scores[candidates[local]], params.burden_quantile))),
                "max_abs_score": float(np.nanmax(np.abs(profile))) if profile.size else 0.0,
            }
        )
        if clone_number > params.max_clones_per_state:
            break

    return labels, pd.DataFrame(component_rows)


def _cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    a = np.nan_to_num(a.astype(np.float32, copy=False), nan=0.0)
    b = np.nan_to_num(b.astype(np.float32, copy=False), nan=0.0)
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0:
        return 0.0
    return float(np.dot(a, b) / denom)


def _merge_global_clones(
    clone_profiles: Dict[str, np.ndarray],
    clone_sizes: Dict[str, int],
    params: FastCNVParams,
) -> Dict[str, str]:
    """Merge state-local clones into at most max_global_clones global clone IDs."""
    global_profiles: List[np.ndarray] = []
    global_sizes: List[int] = []
    assignments: Dict[str, str] = {"WT": "WT"}
    ordered = sorted(clone_profiles, key=lambda key: clone_sizes.get(key, 0), reverse=True)
    for key in ordered:
        profile = clone_profiles[key]
        if not global_profiles:
            global_profiles.append(profile.copy())
            global_sizes.append(clone_sizes.get(key, 1))
            assignments[key] = "clone1"
            continue
        similarities = [_cosine_similarity(profile, existing) for existing in global_profiles]
        best_index = int(np.argmax(similarities))
        should_merge = similarities[best_index] >= params.clone_similarity_threshold
        if should_merge or len(global_profiles) >= params.max_global_clones:
            assignments[key] = f"clone{best_index + 1}"
            size = clone_sizes.get(key, 1)
            total = global_sizes[best_index] + size
            global_profiles[best_index] = ((global_profiles[best_index] * global_sizes[best_index]) + (profile * size)) / max(total, 1)
            global_sizes[best_index] = total
        else:
            global_profiles.append(profile.copy())
            global_sizes.append(clone_sizes.get(key, 1))
            assignments[key] = f"clone{len(global_profiles)}"
    return assignments


def _clone_profile_frame(
    cell_df: pd.DataFrame,
    all_scores: np.ndarray,
    barcode_to_index: Dict[str, int],
    clone_col: str,
) -> tuple[pd.DataFrame, Dict[str, np.ndarray], Dict[str, int]]:
    rows = []
    profiles: Dict[str, np.ndarray] = {}
    sizes: Dict[str, int] = {}
    for clone_id, group in cell_df.groupby(clone_col, sort=True):
        if clone_id == "WT" or not str(clone_id):
            continue
        indices = np.asarray([barcode_to_index[barcode] for barcode in group["CellBarcode"] if barcode in barcode_to_index], dtype=int)
        if indices.size == 0:
            continue
        profile = np.nanmean(all_scores[indices], axis=0)
        profiles[str(clone_id)] = profile
        sizes[str(clone_id)] = int(indices.size)
        rows.append(
            {
                clone_col: clone_id,
                "n_cells": int(indices.size),
                "mean_burden": float(np.nanmean(group["cnv_burden"].astype(float))),
                "max_abs_score": float(np.nanmax(np.abs(profile))),
            }
        )
    return pd.DataFrame(rows), profiles, sizes


def _call_clone_consensus_intervals(
    cell_df: pd.DataFrame,
    all_scores: np.ndarray,
    windows: Sequence[Window],
    slices_by_chr: Dict[str, slice],
    barcode_to_index: Dict[str, int],
    params: FastCNVParams,
) -> pd.DataFrame:
    records = []
    for clone_id, group in cell_df.groupby("global_clone_id", sort=True):
        if clone_id == "WT":
            continue
        indices = np.asarray([barcode_to_index[barcode] for barcode in group["CellBarcode"] if barcode in barcode_to_index], dtype=int)
        if indices.size == 0:
            continue
        clone_scores = all_scores[indices]
        profile = np.nanmean(clone_scores, axis=0)
        calls = _call_intervals_for_cell(
            profile,
            windows,
            slices_by_chr,
            high_threshold=params.high_threshold,
            low_threshold=params.low_threshold,
            min_run_windows=params.min_run_windows,
            min_interval_genes=params.min_interval_genes,
            min_mean_score=params.min_mean_score,
        )
        for call in calls:
            chrom_slice = slices_by_chr[call["chr"]]
            chr_windows = windows[chrom_slice]
            window_mask = np.array(
                [
                    window.start >= call["start"] and window.end <= call["end"]
                    for window in chr_windows
                ],
                dtype=bool,
            )
            absolute_mask = np.zeros(len(windows), dtype=bool)
            absolute_indices = np.arange(chrom_slice.start, chrom_slice.stop)
            absolute_mask[absolute_indices[window_mask]] = True
            sign = 1.0 if call["call"] == "gain" else -1.0
            support = np.nanmean(np.nanmean(clone_scores[:, absolute_mask] * sign, axis=1) >= params.low_threshold)
            if float(support) < params.clone_consensus_fraction:
                continue
            records.append(
                {
                    "global_clone_id": clone_id,
                    "n_cells": int(indices.size),
                    "support_fraction": float(support),
                    **call,
                }
            )
    return pd.DataFrame(records)


def _write_clone_pdf(
    path: Path,
    cell_df: pd.DataFrame,
    all_scores: np.ndarray,
    windows: Sequence[Window],
    barcode_to_index: Dict[str, int],
    low_threshold: float,
) -> Optional[Path]:
    clone_ids = [clone for clone in sorted(cell_df["global_clone_id"].dropna().unique()) if clone != "WT"]
    if not clone_ids:
        return None
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    logging.getLogger("fontTools").setLevel(logging.WARNING)
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    x = np.arange(len(windows))
    matrix = []
    labels = []
    for clone_id in clone_ids:
        group = cell_df[cell_df["global_clone_id"] == clone_id]
        indices = np.asarray([barcode_to_index[barcode] for barcode in group["CellBarcode"] if barcode in barcode_to_index], dtype=int)
        if indices.size == 0:
            continue
        profile = np.nanmean(all_scores[indices], axis=0)
        calls = np.zeros(profile.size, dtype=int)
        calls[profile >= low_threshold] = 1
        calls[profile <= -low_threshold] = -1
        matrix.append(calls)
        labels.append(f"{clone_id} ({indices.size})")
    if not matrix:
        return None

    call_matrix = np.vstack(matrix)
    fig_height = max(2.5, 0.35 * len(labels) + 1.5)
    fig, ax = plt.subplots(figsize=(14, fig_height))
    cmap = ListedColormap(["#2b66c3", "#ffffff", "#d62f27"])
    ax.imshow(call_matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=-1, vmax=1)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Genome windows")
    ax.set_title("fastCNV clone-level genome calls")
    chrom_boundaries = []
    chrom_labels = []
    last_chrom = None
    for index, window in enumerate(windows):
        if window.chrom != last_chrom:
            chrom_boundaries.append(index)
            chrom_labels.append(window.chrom.replace("chr", ""))
            last_chrom = window.chrom
    for boundary in chrom_boundaries:
        ax.axvline(boundary - 0.5, color="#cccccc", linewidth=0.4)
    ax.set_xticks(chrom_boundaries)
    ax.set_xticklabels(chrom_labels, fontsize=7, rotation=90)
    ax.tick_params(axis="x", length=0)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def _format_duration(seconds: float) -> str:
    seconds = max(float(seconds), 0.0)
    if seconds < 60.0:
        return f"{seconds:.2f}s"
    minutes, remainder = divmod(seconds, 60.0)
    if minutes < 60.0:
        return f"{int(minutes)}m{remainder:05.2f}s"
    hours, minutes = divmod(minutes, 60.0)
    return f"{int(hours)}h{int(minutes):02d}m{remainder:05.2f}s"


def run_fastcnv(params: FastCNVParams) -> Dict[str, Path]:
    """Run state-aware first-pass CNV calling and write output artifacts."""
    run_start = time.perf_counter()
    params.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Reading AnnData: %s", params.h5ad)
    adata = ad.read_h5ad(params.h5ad)
    available_obs = ", ".join(adata.obs.columns.astype(str)) or "<none>"
    if params.state_key not in adata.obs.columns:
        raise KeyError(
            f"State column '{params.state_key}' not found in AnnData obs. Available obs columns: {available_obs}"
        )
    if params.sample_key and params.sample_key not in adata.obs.columns:
        raise KeyError(
            f"Sample column '{params.sample_key}' not found in AnnData obs. Available obs columns: {available_obs}"
        )

    matrix = _matrix_for_adata(adata, params.layer)
    if sp.issparse(matrix):
        matrix = matrix.tocsr()
    library_sizes = np.ones(adata.n_obs, dtype=np.float32) if params.input_normalized else _row_sums(matrix)

    if params.gene_coordinates is None:
        raise ValueError("gene_coordinates must be resolved before running fastCNV.")
    coords = load_gene_coordinates(params.gene_coordinates, adata.var_names)
    windows = build_windows(coords, params.window_genes, params.stride_genes, params.min_chr_genes)

    control_reference: Optional[ControlReference] = None
    if params.control_h5ad is not None:
        control_reference = load_control_reference(
            params.control_h5ad,
            adata.var_names,
            input_normalized=params.input_normalized,
            layer=params.layer,
            state_key=params.control_state_key,
            sample_key=params.control_sample_key,
        )
    windows_by_chr: Dict[str, List[Window]] = {}
    for window in windows:
        windows_by_chr.setdefault(window.chrom, []).append(window)
    slices_by_chr = _chr_slices(windows)
    window_ids = np.asarray([window.window_id for window in windows], dtype=object)

    all_scores = np.full((adata.n_obs, len(windows)), np.nan, dtype=np.float32)
    cell_records: List[Dict[str, object]] = []
    interval_records: List[Dict[str, object]] = []

    state_values = adata.obs[params.state_key]
    states = [str(state) for state in pd.Index(state_values[state_values.notna()]).unique().tolist() if str(state)]
    LOGGER.info("Calling CNVs across %d states and %d genome windows.", len(states), len(windows))

    states_loop_start = time.perf_counter()
    for state_number, state in enumerate(states):
        state_start = time.perf_counter()
        all_state_rows = _state_indices(adata.obs, params.state_key, state)
        state_rows = _subsample_indices(all_state_rows, params.max_cells_per_state, params.random_state + state_number)
        skipped = len(state_rows) < params.min_state_cells
        if skipped:
            LOGGER.info("Skipping state %s: %d cells below minimum %d.", state, len(state_rows), params.min_state_cells)
            for row in state_rows:
                cell_records.append(
                    {
                        "CellBarcode": adata.obs_names[row],
                        "cell_state": state,
                        "sample": adata.obs.iloc[row][params.sample_key] if params.sample_key else "",
                        "cnv_status": "low_power",
                        "cnv_burden": np.nan,
                        "is_wt_anchor": False,
                        "baseline_source": "none",
                        "n_cnv_intervals": 0,
                    }
                )
            continue

        control_rows_for_state: Optional[np.ndarray] = None
        if control_reference is not None:
            sample_for_state = ""
            if params.sample_key:
                state_samples = adata.obs.iloc[state_rows][params.sample_key].astype(str)
                if not state_samples.empty:
                    sample_for_state = str(state_samples.mode().iat[0])
            control_rows_for_state = _control_rows_for_state(
                control_reference,
                state=state,
                sample=sample_for_state,
                state_key=params.control_state_key,
                sample_key=params.control_sample_key,
            )
            if control_rows_for_state.size < params.min_control_cells:
                LOGGER.warning(
                    "State %s: only %d matched control cells (< %d); falling back to internal anchors.",
                    state,
                    int(control_rows_for_state.size),
                    params.min_control_cells,
                )
                control_rows_for_state = None

        if control_rows_for_state is not None:
            anchor_mask = np.zeros(len(state_rows), dtype=bool)
            scores = _window_scores_for_state(
                matrix,
                coords,
                windows_by_chr,
                state_rows,
                library_sizes,
                params.input_normalized,
                anchor_mask=None,
                control=control_reference,
                control_rows=control_rows_for_state,
            )
        else:
            first_scores = _window_scores_for_state(
                matrix,
                coords,
                windows_by_chr,
                state_rows,
                library_sizes,
                params.input_normalized,
                anchor_mask=None,
            )
            first_burden = _cnv_burden(first_scores, params.burden_quantile)
            n_anchor = max(params.min_anchor_cells, int(math.ceil(len(state_rows) * params.anchor_fraction)))
            n_anchor = min(max(n_anchor, 1), len(state_rows))
            anchor_order = np.argsort(first_burden, kind="mergesort")
            anchor_mask = np.zeros(len(state_rows), dtype=bool)
            anchor_mask[anchor_order[:n_anchor]] = True

            scores = _window_scores_for_state(
                matrix,
                coords,
                windows_by_chr,
                state_rows,
                library_sizes,
                params.input_normalized,
                anchor_mask=anchor_mask,
            )
        burden = _cnv_burden(scores, params.burden_quantile)
        all_scores[state_rows, :] = scores

        baseline_source = "control" if control_rows_for_state is not None else "internal_anchor"
        for local_index, row in enumerate(state_rows):
            cell_scores = scores[local_index]
            calls = _call_intervals_for_cell(
                cell_scores,
                windows,
                slices_by_chr,
                high_threshold=params.high_threshold,
                low_threshold=params.low_threshold,
                min_run_windows=params.min_run_windows,
                min_interval_genes=params.min_interval_genes,
                min_mean_score=params.min_mean_score,
            )
            status = "CNV" if calls and burden[local_index] >= params.cnv_burden_threshold else "WT"
            barcode = str(adata.obs_names[row])
            sample = str(adata.obs.iloc[row][params.sample_key]) if params.sample_key else ""
            cell_records.append(
                {
                    "CellBarcode": barcode,
                    "cell_state": state,
                    "sample": sample,
                    "cnv_status": status,
                    "cnv_burden": float(burden[local_index]),
                    "is_wt_anchor": bool(anchor_mask[local_index]),
                    "baseline_source": baseline_source,
                    "n_cnv_intervals": len(calls) if status == "CNV" else 0,
                }
            )
            if status != "CNV":
                continue
            for call in calls:
                interval_records.append(
                    {
                        "CellBarcode": barcode,
                        "cell_state": state,
                        "sample": sample,
                        **call,
                    }
                )
        LOGGER.info(
            "State %s: cells=%d anchors=%d CNV=%d (%s)",
            state,
            len(state_rows),
            int(anchor_mask.sum()),
            sum(1 for record in cell_records if record["cell_state"] == state and record["cnv_status"] == "CNV"),
            _format_duration(time.perf_counter() - state_start),
        )
    LOGGER.info("All %d states scored in %s.", len(states), _format_duration(time.perf_counter() - states_loop_start))

    cells_path = params.output_prefix.with_suffix(".cnv_cells.tsv")
    intervals_path = params.output_prefix.with_suffix(".cnv_intervals.tsv")
    windows_path = params.output_prefix.with_suffix(".cnv_windows.tsv")
    scores_path = params.output_prefix.with_suffix(".cnv_window_scores.npz")

    cell_df = pd.DataFrame(cell_records)
    if not cell_df.empty:
        cell_df["state_clone_id"] = "WT"
        cell_df["global_clone_id"] = "WT"
    barcode_to_index = {str(barcode): index for index, barcode in enumerate(adata.obs_names.astype(str))}

    state_clone_rows: List[Dict[str, object]] = []
    if not params.skip_clones and not cell_df.empty:
        clone_group_columns = ["cell_state"]
        if params.sample_key and "sample" in cell_df.columns:
            clone_group_columns.append("sample")
        for state_number, (group_key, state_cell_df) in enumerate(cell_df.groupby(clone_group_columns, sort=False)):
            if isinstance(group_key, tuple):
                state = str(group_key[0])
                sample = str(group_key[1]) if len(group_key) > 1 else ""
            else:
                state = str(group_key)
                sample = ""
            valid_state_df = state_cell_df[state_cell_df["cnv_status"].isin(["CNV", "WT"])]
            if valid_state_df.empty:
                continue
            indices = np.asarray(
                [barcode_to_index[barcode] for barcode in valid_state_df["CellBarcode"] if barcode in barcode_to_index],
                dtype=int,
            )
            if indices.size == 0:
                continue
            state_scores = all_scores[indices]
            cnv_mask = valid_state_df["cnv_status"].to_numpy(dtype=object) == "CNV"
            anchor_mask = valid_state_df["is_wt_anchor"].astype(bool).to_numpy()
            local_labels, local_summary = _discover_state_clones(
                state_scores,
                cnv_mask=cnv_mask,
                anchor_mask=anchor_mask,
                params=params,
                seed=params.random_state + state_number,
            )
            resolved_labels = np.asarray(
                ["WT" if label == "WT" else f"{state}::{sample}::{label}" for label in local_labels],
                dtype=object,
            )
            cell_df.loc[valid_state_df.index, "state_clone_id"] = resolved_labels
            if not local_summary.empty:
                local_summary.insert(0, "cell_state", state)
                local_summary.insert(1, "sample", sample)
                local_summary["state_clone_id"] = local_summary["state_clone_id"].map(lambda label: f"{state}::{sample}::{label}")
                state_clone_rows.extend(local_summary.to_dict("records"))

        _, state_profiles, state_sizes = _clone_profile_frame(cell_df, all_scores, barcode_to_index, "state_clone_id")
        global_assignments = _merge_global_clones(state_profiles, state_sizes, params)
        cell_df["global_clone_id"] = cell_df["state_clone_id"].map(lambda clone: global_assignments.get(str(clone), "WT"))
        for row in state_clone_rows:
            row["global_clone_id"] = global_assignments.get(str(row["state_clone_id"]), "WT")

    clone_intervals_df = (
        _call_clone_consensus_intervals(cell_df, all_scores, windows, slices_by_chr, barcode_to_index, params)
        if not params.skip_clones and not cell_df.empty
        else pd.DataFrame()
    )
    state_clones_path = params.output_prefix.with_suffix(".state_clones.tsv")
    clone_intervals_path = params.output_prefix.with_suffix(".clone_intervals.tsv")
    clone_pdf_path = params.output_prefix.with_suffix(".clone_genome.pdf")
    updated_h5ad_path = params.output_prefix.with_suffix(".fastcnv.h5ad")

    cell_df.to_csv(cells_path, sep="\t", index=False)
    interval_columns = [
        "CellBarcode",
        "cell_state",
        "sample",
        "call",
        "chr",
        "start",
        "end",
        "start_gene",
        "end_gene",
        "n_windows",
        "n_genes",
        "mean_score",
        "max_score",
        "confidence",
    ]
    pd.DataFrame(interval_records, columns=interval_columns).to_csv(intervals_path, sep="\t", index=False)
    pd.DataFrame(state_clone_rows).to_csv(state_clones_path, sep="\t", index=False)
    clone_interval_columns = [
        "global_clone_id",
        "n_cells",
        "support_fraction",
        "call",
        "chr",
        "start",
        "end",
        "start_gene",
        "end_gene",
        "n_windows",
        "n_genes",
        "mean_score",
        "max_score",
        "confidence",
    ]
    pd.DataFrame(clone_intervals_df, columns=clone_interval_columns).to_csv(clone_intervals_path, sep="\t", index=False)
    _window_frame(windows).to_csv(windows_path, sep="\t", index=False)
    np.savez(
        scores_path,
        scores=all_scores,
        cell_barcodes=np.asarray(adata.obs_names.astype(str), dtype=object),
        states=adata.obs[params.state_key].astype(str).to_numpy(dtype=object),
        window_ids=window_ids,
    )
    pdf_path = None
    if not params.skip_pdf and not cell_df.empty:
        pdf_path = _write_clone_pdf(
            clone_pdf_path,
            cell_df,
            all_scores,
            windows,
            barcode_to_index,
            low_threshold=params.low_threshold,
        )

    if params.write_h5ad and not cell_df.empty:
        obs_updates = cell_df.set_index("CellBarcode")
        for column in ("cnv_status", "cnv_burden", "state_clone_id", "global_clone_id"):
            adata.obs[f"fastcnv_{column}"] = obs_updates.reindex(adata.obs_names.astype(str))[column].values
        adata.write_h5ad(updated_h5ad_path)

    outputs = {
        "cells": cells_path,
        "intervals": intervals_path,
        "windows": windows_path,
        "scores": scores_path,
        "state_clones": state_clones_path,
        "clone_intervals": clone_intervals_path,
    }
    if pdf_path is not None:
        outputs["clone_pdf"] = pdf_path
    if params.write_h5ad:
        outputs["h5ad"] = updated_h5ad_path
    LOGGER.info("fastCNV total runtime: %s", _format_duration(time.perf_counter() - run_start))
    return outputs


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run first-pass state-aware fastCNV calls from a cellHarmony h5ad.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--h5ad", required=True, type=Path, help="Input cellHarmony AnnData file.")
    parser.add_argument("--gene-coordinates", default=None, type=Path, help="Gene coordinate TSV/CSV.")
    parser.add_argument("--species", choices=sorted(RESOURCE_FILES), default=None, help="Use a bundled coordinate resource.")
    parser.add_argument("--output", required=True, type=Path, help="Output prefix.")
    parser.add_argument("--state-key", required=True, help="AnnData obs column with cellHarmony cell-state labels.")
    parser.add_argument("--sample-key", default=None, help="Optional AnnData obs sample column.")
    parser.add_argument("--layer", default="auto", help="'auto', 'X', or a layer name. Auto prefers layers['counts'].")
    parser.add_argument("--input-normalized", action="store_true", help="Treat selected matrix values as already log-normalized.")
    parser.add_argument("--control-h5ad", default=None, type=Path, help="Optional external control AnnData used as the WT baseline (overrides per-state internal anchors).")
    parser.add_argument("--control-state-key", default=None, help="Optional obs column in the control AnnData used to match control cells to query cell-state.")
    parser.add_argument("--control-sample-key", default=None, help="Optional obs column in the control AnnData used to match control cells to query sample.")
    parser.add_argument("--min-control-cells", type=int, default=30, help="If fewer than this many control cells match a query state, fall back to internal anchors for that state.")
    parser.add_argument("--window-genes", type=int, default=41)
    parser.add_argument("--stride-genes", type=int, default=7)
    parser.add_argument("--min-chr-genes", type=int, default=25)
    parser.add_argument("--min-state-cells", type=int, default=30)
    parser.add_argument("--anchor-fraction", type=float, default=0.25)
    parser.add_argument("--min-anchor-cells", type=int, default=20)
    parser.add_argument("--high-threshold", type=float, default=2.6)
    parser.add_argument("--low-threshold", type=float, default=1.6)
    parser.add_argument("--min-run-windows", type=int, default=3)
    parser.add_argument("--min-interval-genes", type=int, default=60)
    parser.add_argument("--min-mean-score", type=float, default=1.8)
    parser.add_argument("--burden-quantile", type=float, default=0.95)
    parser.add_argument("--cnv-burden-threshold", type=float, default=1.8)
    parser.add_argument("--max-clones-per-state", type=int, default=10)
    parser.add_argument("--max-global-clones", type=int, default=10)
    parser.add_argument("--min-clone-cells", type=int, default=5)
    parser.add_argument("--clone-similarity-threshold", type=float, default=0.88)
    parser.add_argument("--clone-consensus-fraction", type=float, default=0.45)
    parser.add_argument("--nmf-max-iter", type=int, default=300)
    parser.add_argument("--skip-clones", action="store_true", help="Only write first-pass cell-level CNV calls.")
    parser.add_argument("--skip-pdf", action="store_true", help="Skip clone-level genome PDF export.")
    parser.add_argument("--write-h5ad", action="store_true", help="Write an h5ad copy with fastCNV obs columns.")
    parser.add_argument("--max-cells-per-state", type=int, default=None)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def params_from_args(args: argparse.Namespace) -> FastCNVParams:
    gene_coordinates = Path(args.gene_coordinates) if args.gene_coordinates else None
    if gene_coordinates is None and getattr(args, "species", None):
        gene_coordinates = bundled_gene_coordinates(args.species)
    return FastCNVParams(
        h5ad=Path(args.h5ad),
        gene_coordinates=gene_coordinates,
        output_prefix=Path(args.output),
        state_key=args.state_key,
        sample_key=args.sample_key,
        layer=args.layer,
        input_normalized=bool(args.input_normalized),
        control_h5ad=Path(args.control_h5ad) if args.control_h5ad else None,
        control_state_key=args.control_state_key,
        control_sample_key=args.control_sample_key,
        min_control_cells=int(args.min_control_cells),
        window_genes=int(args.window_genes),
        stride_genes=int(args.stride_genes),
        min_chr_genes=int(args.min_chr_genes),
        min_state_cells=int(args.min_state_cells),
        anchor_fraction=float(args.anchor_fraction),
        min_anchor_cells=int(args.min_anchor_cells),
        high_threshold=float(args.high_threshold),
        low_threshold=float(args.low_threshold),
        min_run_windows=int(args.min_run_windows),
        min_interval_genes=int(args.min_interval_genes),
        min_mean_score=float(args.min_mean_score),
        burden_quantile=float(args.burden_quantile),
        cnv_burden_threshold=float(args.cnv_burden_threshold),
        max_clones_per_state=int(args.max_clones_per_state),
        max_global_clones=int(args.max_global_clones),
        min_clone_cells=int(args.min_clone_cells),
        clone_similarity_threshold=float(args.clone_similarity_threshold),
        clone_consensus_fraction=float(args.clone_consensus_fraction),
        nmf_max_iter=int(args.nmf_max_iter),
        skip_clones=bool(args.skip_clones),
        skip_pdf=bool(args.skip_pdf),
        write_h5ad=bool(args.write_h5ad),
        max_cells_per_state=args.max_cells_per_state,
        random_state=int(args.random_state),
    )


def main(argv: Optional[Sequence[str]] = None) -> Dict[str, Path]:
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s | %(message)s")
    outputs = run_fastcnv(params_from_args(args))
    for name, path in outputs.items():
        LOGGER.info("%s: %s", name, path)
    return outputs


def run_from_altanalyze_args(args: argparse.Namespace) -> Dict[str, Path]:
    return run_fastcnv(params_from_args(args))


if __name__ == "__main__":
    main()
