"""fastCNV: global preprocess + per-chromosome vectorized scoring.

Same scoring model as ``main.py`` (median baseline, 1.4826 * MAD scale, rolling
window mean of residuals, run-length interval calling), but reorganized so the
expensive chunk extractions and normalizations happen once per chromosome
instead of once per (state x chromosome). Includes KMeans-based clone
discovery, cached normalized expression for fast log2-ratio plots, and
zygosity annotations on per-cell and per-clone interval outputs.

Differences vs ``main.py``:
  * Requires a control h5ad. Internal-anchor mode is not implemented here -
    use ``main.py`` if you need it.
  * Single-pass scoring (no burden anchor refinement; the control already
    provides the WT baseline).
"""

from __future__ import annotations

import argparse
import logging
import math
import time
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from altanalyze3.components.fastCNV.main import (
    RESOURCE_FILES,
    FastCNVParams,
    Window,
    _call_clone_consensus_intervals,
    _call_intervals_for_cell,
    _chr_slices,
    _clone_profile_frame,
    _cnv_burden,
    _discover_state_clones,
    _format_duration,
    _merge_global_clones,
    _window_frame,
    bundled_gene_coordinates,
    build_windows,
    load_gene_coordinates,
)


LOGGER = logging.getLogger("fastCNV")


CHRY_Y_ONLY_MARKERS: Tuple[str, ...] = (
    "RPS4Y1", "DDX3Y", "EIF1AY", "KDM5D", "UTY", "USP9Y", "NLGN4Y", "TMSB4Y",
    "PRKY", "ZFY", "TSPY1", "BCORP1", "PRY", "TBL1Y", "TXLNGY", "AMELY",
)
SEX_DETECTION_CHRY_PCT_DEFAULT: float = 5.0


def _infer_sample_sex(
    adata: ad.AnnData,
    sample_labels: np.ndarray,
    threshold_pct: float,
    chry_markers: Sequence[str] = CHRY_Y_ONLY_MARKERS,
) -> pd.DataFrame:
    """Infer per-sample biological sex from chrY-Y-only marker expression.

    For each sample (a unique value in ``sample_labels``), count the fraction of
    cells with >= 1 UMI on any of the chrY-Y-only markers present in
    ``adata.var_names``. Samples with > ``threshold_pct``% chrY-positive cells
    are called male, else female. This is robust to LOY: even in a sample with
    80% LOY cells the unaffected 20% will still express chrY, far above the
    ~1% background expected from a true female sample.

    Returns a DataFrame indexed by sample with columns
    n_cells, n_chrY_pos, pct_chrY_pos, inferred_sex, n_markers_present.
    """
    present = [g for g in chry_markers if g in adata.var_names]
    if len(present) == 0:
        return pd.DataFrame({
            "sample": pd.unique(sample_labels),
            "n_cells": [int((sample_labels == s).sum()) for s in pd.unique(sample_labels)],
            "n_chrY_pos": 0,
            "pct_chrY_pos": 0.0,
            "inferred_sex": "unknown",
            "n_markers_present": 0,
        }).set_index("sample")

    sub = adata[:, present].X
    if sp.issparse(sub):
        sub = sub.tocsr()
        n_pos_per_cell = np.asarray((sub > 0).sum(axis=1)).ravel()
    else:
        n_pos_per_cell = (np.asarray(sub) > 0).sum(axis=1)
    cell_chry_pos = n_pos_per_cell > 0

    rows = []
    for sample in pd.unique(sample_labels):
        mask = sample_labels == sample
        n_cells = int(mask.sum())
        n_pos = int(cell_chry_pos[mask].sum())
        pct = 100.0 * n_pos / max(n_cells, 1)
        rows.append({
            "sample": sample,
            "n_cells": n_cells,
            "n_chrY_pos": n_pos,
            "pct_chrY_pos": pct,
            "inferred_sex": "male" if pct > threshold_pct else "female",
            "n_markers_present": len(present),
        })
    return pd.DataFrame(rows).set_index("sample")


CANCER_DRIVER_GENES: Tuple[str, ...] = (
    "MCL1", "MYCN", "GLI2", "PDGFRA", "FGFR3", "KIT", "EGFR", "CDK6", "MET",
    "FGFR1", "MYC", "CDKN2A", "CDKN2B", "PTCH1", "JAK2", "PTEN", "FAS",
    "CCND1", "ATM", "KRAS", "CDK4", "MDM2", "RB1", "BRCA2", "FOXO1", "NF1",
    "TP53", "ERBB2", "BRCA1", "STK11", "CCNE1", "AURKA", "RUNX1", "AKT1",
    "BCL2", "SMAD4", "APC", "NOTCH1", "FBXW7", "ARID1A", "TET2", "DNMT3A",
    "ASXL1", "EZH2", "SF3B1", "U2AF1", "SRSF2", "IDH1", "IDH2", "FLT3",
    "NPM1", "WT1", "GATA2", "SETBP1", "BCOR", "STAG2", "PHF6", "RAD21",
    "CALR", "MPL",
    # Sex-chromosome dosage markers (always shown on the driver-gene strip plot
    # so LOY appears as a strong negative on RPS4Y1 and X-inactivation/loss
    # appears on XIST):
    "XIST", "RPS4Y1",
)


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
    max_clones_per_state: int = 10
    max_global_clones: int = 10
    min_clone_cells: int = 5
    clone_similarity_threshold: float = 0.88
    clone_consensus_fraction: float = 0.45
    nmf_max_iter: int = 100
    clone_min_active_fraction: float = 0.05
    clone_max_features: int = 400
    clone_min_cnv_fraction: float = 0.015
    clone_min_cells_confident: int = 30
    zygosity_mode: str = "relative"
    skip_clones: bool = False
    skip_pdf: bool = False
    write_h5ad: bool = False
    random_state: int = 0
    n_jobs: int = -1
    pdf_smooth_genes: int = 50
    pdf_y_clip: float = 1.0
    pdf_label_threshold: float = 0.25
    pdf_score_y_clip: float = 4.0
    sex_chrom_mode: str = "absolute_log2"
    sex_chrom_log2_unit: float = 0.040
    sex_detection_threshold_pct: float = 5.0
    control_sample_key: Optional[str] = None
    sex_chrom_het_loss: float = -0.6
    sex_chrom_hom_loss: float = -1.5
    sex_chrom_het_gain: float = 0.4
    sex_chrom_hom_gain: float = 0.7
    sex_chroms: Tuple[str, ...] = ("chrY",)


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


def _absolute_log2_window_scores(
    query_state_log1p: np.ndarray,
    control_state_log1p: np.ndarray,
    chr_windows: Sequence[Window],
    log2_unit: float,
) -> np.ndarray:
    """Per-cell sex-chrom window scores from absolute expression difference (log1p space).

    Bypasses MAD-standardization, which collapses real LOY signal because the male
    reference's chrY MAD is large (donor-level age-related sub-clinical LOY).

    For each cell we compute a per-gene difference (query_log1p - control_state_mean_log1p),
    rolling-window mean it, then divide by ``log2_unit`` so the result lives in the same
    standardized-score units (~|score| in [-3, +3]) used by the autosomal interval caller.
    Working in log1p space avoids the divide-by-zero / pseudocount-tuning problem of
    CP10K log2 ratios on sparse sex-chromosome counts. The mean log1p difference is
    monotonic with biological log2 fold-change for moderately expressed genes.
    """
    if query_state_log1p.size == 0 or control_state_log1p.size == 0:
        return np.zeros((query_state_log1p.shape[0], len(chr_windows)), dtype=np.float32)
    control_mean_log1p = np.nanmean(control_state_log1p, axis=0).astype(np.float32)
    diff_per_gene = (query_state_log1p.astype(np.float32) - control_mean_log1p[None, :])
    windowed = _rolling_window_mean(diff_per_gene, chr_windows)
    unit = max(float(log2_unit), 1e-3)
    return (windowed / unit).astype(np.float32)


def _control_var_map(query_var: pd.Index, control_var: pd.Index) -> np.ndarray:
    lookup = pd.Series(np.arange(len(control_var), dtype=np.int64), index=control_var.astype(str))
    return lookup.reindex(query_var.astype(str)).to_numpy()


def _kmeans_state_clones(
    scores: np.ndarray,
    cnv_mask: np.ndarray,
    max_clones: int,
    min_cells: int,
    seed: int,
) -> Tuple[np.ndarray, pd.DataFrame]:
    """KMeans replacement for NMF-based state-local clone discovery.

    Operates only on CNV-positive cells. Returns per-cell labels (WT or clone1..N)
    and a summary frame with size + mean burden per clone.
    """
    labels = np.full(scores.shape[0], "WT", dtype=object)
    candidates = np.flatnonzero(cnv_mask)
    if candidates.size < min_cells:
        return labels, pd.DataFrame()

    feature_block = np.nan_to_num(scores[candidates].astype(np.float32, copy=False), nan=0.0, posinf=0.0, neginf=0.0)
    feature_block = np.clip(feature_block, -8.0, 8.0)
    if feature_block.shape[1] == 0:
        return labels, pd.DataFrame()

    n_clusters = max(1, min(int(max_clones), candidates.size // max(min_cells, 1)))
    if n_clusters <= 1:
        clone_labels = np.zeros(candidates.size, dtype=int)
    else:
        from sklearn.cluster import MiniBatchKMeans
        model = MiniBatchKMeans(
            n_clusters=n_clusters,
            random_state=seed,
            n_init=3,
            batch_size=min(1024, max(64, candidates.size)),
            max_iter=100,
        )
        clone_labels = model.fit_predict(feature_block).astype(int)

    rows: List[Dict[str, object]] = []
    order = sorted(set(int(c) for c in clone_labels), key=lambda c: int(np.sum(clone_labels == c)), reverse=True)
    next_id = 1
    for cluster in order:
        member_local = np.flatnonzero(clone_labels == cluster)
        if member_local.size < min_cells:
            continue
        clone_id = f"clone{next_id}"
        next_id += 1
        labels[candidates[member_local]] = clone_id
        cluster_scores = feature_block[member_local]
        rows.append({
            "state_clone_id": clone_id,
            "n_cells": int(member_local.size),
            "mean_burden": float(np.nanmean(np.abs(cluster_scores).max(axis=1))),
            "max_abs_score": float(np.nanmax(np.abs(cluster_scores))) if cluster_scores.size else 0.0,
        })
        if next_id > max_clones:
            break
    return labels, pd.DataFrame(rows)


def _smooth_and_center_log2_ratios(
    ratios: np.ndarray,
    coords: pd.DataFrame,
    smooth_genes: int,
) -> np.ndarray:
    """Smooth log2 ratios across a gene-neighborhood window per chromosome and center each chromosome on its median."""
    if ratios.size == 0:
        return ratios
    smoothed = np.full_like(ratios, np.nan, dtype=np.float32)
    for chrom, chr_coords in coords.groupby("chr", sort=False):
        idx = chr_coords.index.to_numpy(dtype=np.int64)
        if idx.size == 0:
            continue
        block = ratios[:, idx]
        valid = ~np.isnan(block)
        block_filled = np.where(valid, block, 0.0).astype(np.float32, copy=False)
        valid_f = valid.astype(np.float32)
        n_genes = block.shape[1]
        if smooth_genes > 1 and n_genes >= 3:
            half = max(1, smooth_genes // 2)
            sum_prefix = np.concatenate([np.zeros((block.shape[0], 1), dtype=np.float32), np.cumsum(block_filled, axis=1, dtype=np.float32)], axis=1)
            count_prefix = np.concatenate([np.zeros((block.shape[0], 1), dtype=np.float32), np.cumsum(valid_f, axis=1, dtype=np.float32)], axis=1)
            starts = np.maximum(np.arange(n_genes) - half, 0)
            ends = np.minimum(np.arange(n_genes) + half + 1, n_genes)
            window_sum = sum_prefix[:, ends] - sum_prefix[:, starts]
            window_count = count_prefix[:, ends] - count_prefix[:, starts]
            chrom_smoothed = np.where(window_count > 0, window_sum / np.maximum(window_count, 1.0), np.nan).astype(np.float32)
        else:
            chrom_smoothed = block.astype(np.float32, copy=True)
        median = np.nanmedian(chrom_smoothed, axis=1, keepdims=True)
        median = np.nan_to_num(median, nan=0.0)
        chrom_smoothed = chrom_smoothed - median
        smoothed[:, idx] = chrom_smoothed
    return smoothed


_FIXED_ZYGOSITY_THRESHOLDS = {
    "homozygous_loss": -1.30,
    "heterozygous_loss": -0.55,
    "heterozygous_gain": 0.40,
    "homozygous_gain": 0.85,
}


def _zygosity_state(log2_ratio: float, thresholds: Optional[Dict[str, float]] = None) -> str:
    """Map a smoothed log2 ratio to a zygosity-like categorical state."""
    if thresholds is None:
        thresholds = _FIXED_ZYGOSITY_THRESHOLDS
    if not np.isfinite(log2_ratio):
        return "low_signal"
    if log2_ratio <= thresholds["homozygous_loss"]:
        return "homozygous_loss"
    if log2_ratio <= thresholds["heterozygous_loss"]:
        return "heterozygous_loss"
    if log2_ratio >= thresholds["homozygous_gain"]:
        return "homozygous_gain"
    if log2_ratio >= thresholds["heterozygous_gain"]:
        return "heterozygous_gain"
    return "low_signal"


def _calibrate_zygosity_thresholds(
    interval_records: Sequence[Dict[str, object]],
    wt_indices: np.ndarray,
    cached_query_expr: np.ndarray,
    control_gene_means: np.ndarray,
    coords: pd.DataFrame,
    seed: int = 0,
    cells_per_interval: int = 50,
    max_intervals: int = 2000,
) -> Optional[Dict[str, float]]:
    """Estimate empirical zygosity thresholds from WT cells' log2 ratios in called intervals.

    Returns a dict of thresholds, or None if there's not enough data to calibrate.
    """
    if cached_query_expr is None or wt_indices.size < 100 or len(interval_records) == 0:
        return None
    rng = np.random.default_rng(seed)
    coords_chr = coords["chr"].astype(str).to_numpy()
    coords_start = coords["start"].to_numpy()
    coords_end = coords["end"].to_numpy()
    chr_to_indices: Dict[str, np.ndarray] = {
        chrom: np.flatnonzero(coords_chr == chrom).astype(np.int64)
        for chrom in pd.unique(coords_chr)
    }
    pseudocount = 1e-3

    if len(interval_records) > max_intervals:
        sampled = rng.choice(len(interval_records), size=max_intervals, replace=False)
        sampled_records = [interval_records[i] for i in sampled]
    else:
        sampled_records = list(interval_records)

    noise_means: List[float] = []
    for record in sampled_records:
        chrom = str(record["chr"])
        start = int(record["start"])
        end = int(record["end"])
        chr_idx = chr_to_indices.get(chrom)
        if chr_idx is None or chr_idx.size == 0:
            continue
        in_interval = (coords_start[chr_idx] >= start) & (coords_end[chr_idx] <= end)
        gene_indices = chr_idx[in_interval]
        if gene_indices.size == 0:
            continue
        sample_size = min(cells_per_interval, wt_indices.size)
        wt_sample = rng.choice(wt_indices, size=sample_size, replace=False)
        cell_expr = cached_query_expr[wt_sample][:, gene_indices]
        ctrl_expr = control_gene_means[gene_indices]
        per_cell_means = np.nanmean(np.log2((cell_expr + pseudocount) / (ctrl_expr + pseudocount)), axis=1)
        noise_means.extend(per_cell_means.tolist())

    if len(noise_means) < 200:
        return None
    arr = np.asarray(noise_means, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size < 200:
        return None

    q_low = float(np.quantile(arr, 0.005))
    q_lo = float(np.quantile(arr, 0.05))
    q_hi = float(np.quantile(arr, 0.95))
    q_high = float(np.quantile(arr, 0.995))
    thresholds = {
        "homozygous_loss": min(q_low, -1.0),
        "heterozygous_loss": min(q_lo, -0.4),
        "heterozygous_gain": max(q_hi, 0.4),
        "homozygous_gain": max(q_high, 0.7),
    }
    return thresholds


def _annotate_intervals_with_zygosity(
    interval_records: List[Dict[str, object]],
    cached_query_expr: Optional[np.ndarray],
    control_gene_means: np.ndarray,
    coords: pd.DataFrame,
    barcode_to_index: Dict[str, int],
    thresholds: Optional[Dict[str, float]] = None,
    sex_chroms: Sequence[str] = (),
    sex_chrom_thresholds: Optional[Dict[str, float]] = None,
) -> List[Dict[str, object]]:
    """Attach mean_log2_ratio + zygosity_state to per-cell interval records.

    For autosomal intervals, the per-gene log2 ratio uses the cached log1p(CP10K)
    expression with a tiny pseudocount (legacy behavior, valid where both sides
    are non-trivial).

    For intervals on a chromosome in ``sex_chroms``, ``mean_log2_ratio`` is the
    mean log1p(CP10K) DIFFERENCE (query - control_state_mean), not a log2 ratio.
    The autosomal log2-of-log1p formulation explodes when control expression
    approaches zero (chrY in female-mixed or sparse references). The log1p-diff
    is monotonic with biological log2 fold-change for moderately expressed genes
    and is well-defined at any expression level. ``sex_chrom_thresholds`` must be
    calibrated for log1p-diff units (defaults: het_loss<=-0.04, hom_loss<=-0.08).
    """
    if not interval_records or cached_query_expr is None:
        return interval_records
    pseudocount = 1e-3
    coords_chr = coords["chr"].astype(str).to_numpy()
    coords_start = coords["start"].to_numpy()
    coords_end = coords["end"].to_numpy()
    chr_to_indices: Dict[str, np.ndarray] = {
        chrom: np.flatnonzero(coords_chr == chrom).astype(np.int64)
        for chrom in pd.unique(coords_chr)
    }
    sex_chrom_set = {str(c) for c in sex_chroms}
    annotated: List[Dict[str, object]] = []
    for record in interval_records:
        chrom = str(record["chr"])
        start = int(record["start"])
        end = int(record["end"])
        chr_idx = chr_to_indices.get(chrom)
        if chr_idx is None or chr_idx.size == 0:
            record["mean_log2_ratio"] = float("nan")
            record["zygosity_state"] = "low_signal"
            annotated.append(record)
            continue
        in_interval_local = (coords_start[chr_idx] >= start) & (coords_end[chr_idx] <= end)
        gene_indices = chr_idx[in_interval_local]
        cell_idx = barcode_to_index.get(str(record["CellBarcode"]))
        if cell_idx is None or gene_indices.size == 0:
            record["mean_log2_ratio"] = float("nan")
            record["zygosity_state"] = "low_signal"
            annotated.append(record)
            continue
        cell_expr = cached_query_expr[cell_idx, gene_indices]
        ctrl_expr = control_gene_means[gene_indices]
        is_sex_chrom = chrom in sex_chrom_set and sex_chrom_thresholds is not None
        per_gene_log2 = np.log2((cell_expr + pseudocount) / (ctrl_expr + pseudocount))
        mean_metric = float(np.nanmean(per_gene_log2))
        record["mean_log2_ratio"] = mean_metric
        if is_sex_chrom:
            record["zygosity_state"] = _zygosity_state(mean_metric, sex_chrom_thresholds)
        else:
            record["zygosity_state"] = _zygosity_state(mean_metric, thresholds)
        annotated.append(record)
    return annotated


def _compute_clone_log2_ratios_from_cache(
    cached_query_expr: np.ndarray,
    control_gene_means: np.ndarray,
    clone_to_query_rows: Dict[str, np.ndarray],
) -> Tuple[np.ndarray, List[str]]:
    """Per-clone per-gene log2(mean_clone / mean_control) using a precomputed normalized matrix."""
    clone_order = list(clone_to_query_rows.keys())
    n_genes = cached_query_expr.shape[1]
    ratios = np.full((len(clone_order), n_genes), np.nan, dtype=np.float32)
    pseudocount = 1e-3
    for i, clone_id in enumerate(clone_order):
        rows = clone_to_query_rows[clone_id]
        if rows.size == 0:
            continue
        clone_mean = cached_query_expr[rows].mean(axis=0).astype(np.float32)
        ratios[i] = np.log2((clone_mean + pseudocount) / (control_gene_means + pseudocount))
    return ratios, clone_order


def _compute_clone_log2_ratios(
    query_matrix: sp.spmatrix | np.ndarray,
    query_lib: np.ndarray,
    coords: pd.DataFrame,
    control_gene_means: np.ndarray,
    clone_to_query_rows: Dict[str, np.ndarray],
    input_normalized: bool,
) -> Tuple[np.ndarray, List[str]]:
    """Compute per-clone log2(mean_clone / mean_control) for every coord-mapped gene.

    Returns (ratios, clone_order) with ratios of shape (n_clones, n_coord_genes).
    """
    clone_order = list(clone_to_query_rows.keys())
    n_genes = coords.shape[0]
    ratios = np.full((len(clone_order), n_genes), np.nan, dtype=np.float32)
    if not clone_order:
        return ratios, clone_order

    pseudocount = 1e-3
    for chrom, chr_coords in coords.groupby("chr", sort=False):
        query_cols = chr_coords["var_index"].to_numpy(dtype=np.int64)
        coord_rows = chr_coords.index.to_numpy(dtype=np.int64)
        for i, clone_id in enumerate(clone_order):
            rows = clone_to_query_rows[clone_id]
            if rows.size == 0:
                continue
            sub = query_matrix[rows][:, query_cols]
            if sp.issparse(sub):
                sub = sub.toarray()
            sub = np.asarray(sub, dtype=np.float32)
            if not input_normalized:
                factors = (10000.0 / np.maximum(query_lib[rows].astype(np.float32), 1.0)).astype(np.float32)
                sub *= factors[:, None]
                np.log1p(sub, out=sub)
            clone_mean = np.nanmean(sub, axis=0).astype(np.float32)
            control_mean = control_gene_means[coord_rows]
            ratios[i, coord_rows] = np.log2((clone_mean + pseudocount) / (control_mean + pseudocount))
    return ratios, clone_order


def _genome_x_positions(coords: pd.DataFrame) -> Tuple[np.ndarray, List[str], List[Tuple[float, float]]]:
    """Return per-gene cumulative genome position, ordered chromosome list, and per-chromosome (start, end) bounds."""
    chr_to_offset: Dict[str, float] = {}
    chr_order: List[str] = []
    bounds: List[Tuple[float, float]] = []
    cumulative = 0.0
    for chrom, chr_coords in coords.groupby("chr", sort=False):
        size = float(chr_coords["end"].max() - chr_coords["start"].min())
        chr_to_offset[chrom] = cumulative - float(chr_coords["start"].min())
        chr_order.append(str(chrom))
        bounds.append((cumulative, cumulative + size))
        cumulative += size + size * 0.005  # 0.5% gap between chromosomes
    midpoints = ((coords["start"] + coords["end"]) // 2).to_numpy(dtype=np.float64)
    chr_offsets = coords["chr"].map(chr_to_offset).to_numpy(dtype=np.float64)
    return midpoints + chr_offsets, chr_order, bounds


def _write_clone_genome_pdf(
    path: Path,
    coords: pd.DataFrame,
    clone_score_profiles: np.ndarray,
    windows: Sequence[Window],
    clone_order: Sequence[str],
    clone_sizes: Dict[str, int],
    clone_intervals_df: pd.DataFrame,
    score_y_clip: float,
    score_threshold: float,
) -> Optional[Path]:
    """Genome-wide CNV-score visualization (one page per clone).

    Plots per-clone mean window score along the genome rather than per-gene log2(clone/control).
    The score is the same standardized residual the interval caller uses, so what's drawn
    tracks the CNV calls directly and is not biased by transcriptional differential expression
    that's unrelated to copy number (a major confounder in AML/MDS where blast vs. healthy
    cell-state programs dominate per-gene log2).
    """
    if clone_score_profiles.shape[0] == 0 or not windows:
        return None
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LinearSegmentedColormap

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    # Window midpoints positioned along a continuous genome x-axis using the same chrom-bound
    # layout as the per-gene plot, so chromosome boundaries align with what users expect.
    _, chr_order, chr_bounds = _genome_x_positions(coords)
    chr_centers = [(s + e) / 2.0 for s, e in chr_bounds]
    chr_labels = [c.replace("chr", "") for c in chr_order]
    chr_to_bounds = dict(zip(chr_order, chr_bounds))
    chr_min_start_lookup = {
        chrom: float(coords.loc[coords["chr"] == chrom, "start"].min())
        for chrom in chr_order
    }

    window_x = np.full(len(windows), np.nan, dtype=np.float64)
    for i, w in enumerate(windows):
        chrom = str(w.chrom)
        if chrom not in chr_to_bounds:
            continue
        chr_x_start, _ = chr_to_bounds[chrom]
        chr_min_start = chr_min_start_lookup[chrom]
        midpoint_genomic = (float(w.start) + float(w.end)) / 2.0
        window_x[i] = chr_x_start + (midpoint_genomic - chr_min_start)
    valid_x_mask = np.isfinite(window_x)

    cmap = LinearSegmentedColormap.from_list(
        "fastcnv_diverging",
        [(0.0, "#c8252c"), (0.5, "#cccccc"), (1.0, "#1f863d")],
    )

    intervals_by_clone: Dict[str, pd.DataFrame] = {}
    if not clone_intervals_df.empty:
        intervals_by_clone = {
            str(cid): grp for cid, grp in clone_intervals_df.groupby("global_clone_id", sort=False)
        }

    # Color saturation point: anything beyond ±score_threshold is full red/green. Below it,
    # colors fade toward grey so windows with sub-threshold (call-failed) noise don't draw the eye.
    color_extent = max(float(score_threshold), 1e-6)

    with PdfPages(path) as pdf:
        for clone_idx, clone_id in enumerate(clone_order):
            # Skip clones that produced no called consensus intervals — NMF placed cells in
            # this clone but the clone-level signal didn't survive the interval caller, so the
            # genome plot would be uninformative (no stars, just noise).
            clone_intervals = intervals_by_clone.get(str(clone_id))
            if clone_intervals is None or clone_intervals.empty:
                continue

            clone_scores_chr = clone_score_profiles[clone_idx]
            # Recenter on the genome-wide median so each clone's baseline sits at zero. The
            # standardized residual is centered per-window-per-state during scoring, but a
            # global library/protocol mismatch between the query cohort and the reference
            # leaves a residual offset (typically -0.5 to -1.0) that's identical across the
            # whole genome — not CNV. Subtracting the per-clone median for display exposes
            # the contiguous deviations that actually correspond to copy number.
            finite = np.isfinite(clone_scores_chr)
            clone_offset = float(np.nanmedian(clone_scores_chr[finite])) if finite.any() else 0.0
            clone_scores_centered = clone_scores_chr - clone_offset
            scores_clipped = np.clip(np.nan_to_num(clone_scores_centered, nan=0.0), -score_y_clip, score_y_clip)

            fig, ax = plt.subplots(figsize=(14, 4.5))
            ax.scatter(
                window_x[valid_x_mask], scores_clipped[valid_x_mask],
                c=scores_clipped[valid_x_mask], cmap=cmap, vmin=-color_extent, vmax=color_extent,
                s=12, linewidths=0, alpha=0.95, rasterized=True,
            )
            ax.axhline(0.0, color="#888888", linewidth=0.5)
            ax.axhline(score_threshold, color="#bbbbbb", linewidth=0.4, linestyle=":")
            ax.axhline(-score_threshold, color="#bbbbbb", linewidth=0.4, linestyle=":")
            for chrom in chr_order:
                start_x, _ = chr_to_bounds[chrom]
                ax.axvline(start_x, color="#cccccc", linewidth=0.4, linestyle="--")
            ax.set_xticks(chr_centers)
            ax.set_xticklabels(chr_labels, fontsize=8)
            ax.set_xlim(chr_bounds[0][0], chr_bounds[-1][1])
            ax.set_ylim(-score_y_clip - 0.1, score_y_clip + 0.1)
            ax.set_ylabel("CNV score  (standardized residual)")
            ax.set_xlabel("Chromosome")
            n_cells = int(clone_sizes.get(str(clone_id), 0))
            ax.set_title(f"{clone_id}  (n={n_cells} cells)")

            # Overlay called intervals as black horizontal segments at the interval mean score
            # (in the same recentered units as the scatter), with a black star at the midpoint.
            window_chr = np.array([str(w.chrom) for w in windows])
            window_start = np.array([w.start for w in windows], dtype=np.int64)
            window_end = np.array([w.end for w in windows], dtype=np.int64)
            for _, interval_row in clone_intervals.iterrows():
                chrom = str(interval_row["chr"])
                if chrom not in chr_to_bounds:
                    continue
                chr_x_start, _ = chr_to_bounds[chrom]
                chr_min_start = chr_min_start_lookup[chrom]
                seg_x0 = chr_x_start + (float(interval_row["start"]) - chr_min_start)
                seg_x1 = chr_x_start + (float(interval_row["end"]) - chr_min_start)
                in_segment = (
                    (window_chr == chrom)
                    & (window_start >= int(interval_row["start"]))
                    & (window_end <= int(interval_row["end"]))
                )
                if not in_segment.any():
                    continue
                seg_y = float(np.nanmean(clone_scores_centered[in_segment]))
                seg_y_clipped = float(np.clip(seg_y, -score_y_clip, score_y_clip))
                ax.plot(
                    [seg_x0, seg_x1], [seg_y_clipped, seg_y_clipped],
                    color="black", linewidth=2.5, solid_capstyle="butt", zorder=6,
                )
                ax.scatter(
                    [(seg_x0 + seg_x1) / 2.0], [seg_y_clipped],
                    marker="*", s=140, color="black", edgecolor="white", linewidth=0.6, zorder=7,
                )
            fig.tight_layout()
            pdf.savefig(fig, dpi=200)
            plt.close(fig)
    return path


def _write_clone_pdf_hardened(
    path: Path,
    cell_df: pd.DataFrame,
    all_scores: np.ndarray,
    windows: Sequence[Window],
    barcode_to_index: Dict[str, int],
    low_threshold: float,
) -> Optional[Path]:
    """PDF writer that forces the Agg backend so the file is portable across viewers."""
    def _clone_sort_key(clone_id: str) -> Tuple[int, str]:
        text = str(clone_id)
        if text.startswith("clone") and text[5:].isdigit():
            return (0, f"{int(text[5:]):06d}")
        return (1, text)

    clone_ids = [c for c in sorted(cell_df["global_clone_id"].dropna().unique(), key=_clone_sort_key) if c != "WT"]
    if not clone_ids:
        return None
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import ListedColormap

    logging.getLogger("fontTools").setLevel(logging.WARNING)
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    matrix_rows: List[np.ndarray] = []
    labels: List[str] = []
    for clone_id in clone_ids:
        group = cell_df[cell_df["global_clone_id"] == clone_id]
        indices = np.asarray(
            [barcode_to_index[b] for b in group["CellBarcode"] if b in barcode_to_index],
            dtype=int,
        )
        if indices.size == 0:
            continue
        profile = np.nanmean(all_scores[indices], axis=0)
        calls = np.zeros(profile.size, dtype=int)
        calls[profile >= low_threshold] = 1
        calls[profile <= -low_threshold] = -1
        matrix_rows.append(calls)
        labels.append(f"{clone_id} (n={indices.size})")
    if not matrix_rows:
        return None
    call_matrix = np.vstack(matrix_rows)

    chrom_boundaries: List[int] = []
    chrom_labels: List[str] = []
    last_chrom = None
    for index, window in enumerate(windows):
        if window.chrom != last_chrom:
            chrom_boundaries.append(index)
            chrom_labels.append(window.chrom.replace("chr", ""))
            last_chrom = window.chrom

    fig_height = max(2.5, 0.35 * len(labels) + 1.5)
    cmap = ListedColormap(["#2b66c3", "#ffffff", "#d62f27"])
    n_rows, n_cols = call_matrix.shape
    with PdfPages(path) as pdf:
        fig, ax = plt.subplots(figsize=(14, fig_height))
        im = ax.imshow(
            call_matrix,
            aspect="auto",
            interpolation="nearest",
            cmap=cmap,
            vmin=-1,
            vmax=1,
            origin="upper",
            extent=(-0.5, n_cols - 0.5, n_rows - 0.5, -0.5),
        )
        im.set_rasterized(True)
        ax.set_yticks(np.arange(n_rows))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("Genome windows")
        ax.set_title("fastCNV clone-level genome calls")
        for boundary in chrom_boundaries:
            ax.axvline(boundary - 0.5, color="#cccccc", linewidth=0.4)
        ax.set_xticks(chrom_boundaries)
        ax.set_xticklabels(chrom_labels, fontsize=7, rotation=90)
        ax.tick_params(axis="x", length=0)
        fig.tight_layout()
        pdf.savefig(fig, dpi=200)
        plt.close(fig)
    return path


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

    # ---- Automatic per-sample sex detection ----------------------------------
    # ChrY scoring is meaningless when the query or control samples are female,
    # and pooling male+female controls flattens the chrY MAD/baseline so LOY
    # in male queries gets squashed. Detect sex per sample (>{threshold}% cells
    # with >=1 chrY-Y-only marker UMI -> male, else female), then gate chrY
    # scoring: only male query cells get scored against a male-only control pool.
    if params.sample_key and params.sample_key in query.obs.columns:
        query_sample_labels = query.obs[params.sample_key].astype(str).to_numpy()
    else:
        query_sample_labels = np.full(query.n_obs, "all", dtype=object)
    query_sex_df = _infer_sample_sex(query, query_sample_labels, params.sex_detection_threshold_pct)
    sample_to_query_sex = query_sex_df["inferred_sex"].to_dict()
    query_cell_sex = np.array(
        [sample_to_query_sex.get(s, "unknown") for s in query_sample_labels], dtype=object,
    )
    male_query_mask = (query_cell_sex == "male")

    if params.control_sample_key and params.control_sample_key in control.obs.columns:
        control_sample_labels = control.obs[params.control_sample_key].astype(str).to_numpy()
    else:
        control_sample_labels = np.full(control.n_obs, "all", dtype=object)
    control_sex_df = _infer_sample_sex(control, control_sample_labels, params.sex_detection_threshold_pct)
    sample_to_control_sex = control_sex_df["inferred_sex"].to_dict()
    control_cell_sex = np.array(
        [sample_to_control_sex.get(s, "unknown") for s in control_sample_labels], dtype=object,
    )
    male_control_rows = np.flatnonzero(control_cell_sex == "male").astype(np.int64)

    LOGGER.info(
        "Sex detection (>%g%% chrY-marker-positive cells = male):",
        params.sex_detection_threshold_pct,
    )
    for src_label, df in (("query", query_sex_df), ("control", control_sex_df)):
        for sample, row in df.iterrows():
            LOGGER.info(
                "  %s sample=%s n=%d chrY+ %d/%d (%.1f%%) -> %s",
                src_label, sample, int(row["n_cells"]),
                int(row["n_chrY_pos"]), int(row["n_cells"]),
                float(row["pct_chrY_pos"]), str(row["inferred_sex"]),
            )

    n_male_q = int(male_query_mask.sum())
    n_male_c = int(male_control_rows.size)
    chry_scoring_enabled = (
        params.sex_chrom_mode == "absolute_log2"
        and "chrY" in {str(c) for c in params.sex_chroms}
        and n_male_q > 0 and n_male_c > 0
    )
    if not chry_scoring_enabled and "chrY" in {str(c) for c in params.sex_chroms}:
        if n_male_q == 0:
            LOGGER.warning("No male samples in query; chrY scoring disabled.")
        elif n_male_c == 0:
            LOGGER.warning("No male cells in control; chrY scoring disabled (LOY cannot be detected).")

    all_scores = np.full((query.n_obs, len(windows)), np.nan, dtype=np.float32)
    control_gene_means = np.full(coords.shape[0], np.nan, dtype=np.float32)
    coords_var_lookup = pd.Series(
        np.arange(coords.shape[0], dtype=np.int64), index=coords["var_index"].to_numpy()
    )
    cache_query_expr = not params.skip_pdf
    if cache_query_expr:
        cached_query_expr = np.zeros((query.n_obs, coords.shape[0]), dtype=np.float32)
        LOGGER.info(
            "Caching normalized query expression: %d x %d float32 (~%.1f GB).",
            query.n_obs, coords.shape[0], cached_query_expr.nbytes / 1e9,
        )
    else:
        cached_query_expr = None
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
        if cached_query_expr is not None:
            cached_query_expr[:, chr_coords.index.to_numpy(dtype=np.int64)] = query_expr

        mapped = var_map[query_cols]
        valid_mask = ~pd.isna(mapped)
        valid_control_cols = mapped[valid_mask].astype(np.int64) if valid_mask.any() else np.empty(0, dtype=np.int64)
        if valid_control_cols.size == 0:
            LOGGER.warning("Chromosome %s: no overlapping control genes; skipping.", chrom_str)
            continue
        control_expr_full = np.zeros((control.n_obs, query_cols.size), dtype=np.float32)
        control_sub = _normalize_chunk(control_matrix, valid_control_cols, control_lib, params.input_normalized)
        control_expr_full[:, valid_mask] = control_sub
        coord_rows_for_chr = coords_var_lookup.reindex(query_cols).to_numpy()

        is_sex_chrom = (
            params.sex_chrom_mode == "absolute_log2"
            and chrom_str in params.sex_chroms
        )
        # For sex chromosomes (specifically chrY when scoring is active), the per-gene
        # control mean must be computed from male donors only — otherwise female zeros drag
        # the baseline toward zero and the per-clone log2(clone/control) used by the PDF
        # collapses to ~0 in LOY cells (where clone_mean is also ~0). The result was a
        # visually flat chrY band even though LOY is biologically dramatic. Restricting
        # to males here makes the PDF show the real ~log2(-2 to -5) chrY drop in LOY clones.
        if is_sex_chrom and chrom_str == "chrY" and chry_scoring_enabled and male_control_rows.size > 0:
            gene_chunk_means = np.nanmean(control_expr_full[male_control_rows], axis=0).astype(np.float32)
        else:
            gene_chunk_means = np.nanmean(control_expr_full, axis=0).astype(np.float32)
        control_gene_means[coord_rows_for_chr] = gene_chunk_means
        del control_sub
        # ChrY scoring requires both male query cells AND male control cells. If either is
        # missing for this chromosome, fall back to standard MAD-standardized scoring (which
        # will produce noise / no calls — but won't fabricate a LOY signal in females).
        chry_active = is_sex_chrom and chrom_str == "chrY" and chry_scoring_enabled
        # Pool of control rows used for sex-chrom baseline: males only when scoring chrY.
        sex_chrom_control_rows = male_control_rows if chry_active else all_control_rows
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
            if is_sex_chrom:
                scores_chr = _absolute_log2_window_scores(
                    query_expr, control_expr_full[sex_chrom_control_rows],
                    chr_windows, params.sex_chrom_log2_unit,
                )
                if chry_active and not male_query_mask.all():
                    # Female (and unknown-sex) query cells: chrY signal is biologically zero,
                    # any "loss" call would be spurious. Mark NaN so they're skipped downstream.
                    scores_chr[~male_query_mask] = np.nan
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
                if is_sex_chrom:
                    if chry_active:
                        # Restrict per-state control pool to males
                        control_rows = np.intersect1d(control_rows, male_control_rows, assume_unique=False)
                        if control_rows.size == 0:
                            control_rows = male_control_rows
                    state_scores = _absolute_log2_window_scores(
                        query_expr[rows], control_expr_full[control_rows],
                        chr_windows, params.sex_chrom_log2_unit,
                    )
                    if chry_active:
                        # NaN-out female / unknown query cells within this state
                        state_male_mask = male_query_mask[rows]
                        if not state_male_mask.all():
                            state_scores[~state_male_mask] = np.nan
                else:
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

    sex_chrom_set_runtime = (
        {str(c) for c in params.sex_chroms}
        if params.sex_chrom_mode == "absolute_log2"
        else set()
    )
    # Mask of autosomal windows for the burden filter. Sex-chrom scores are biological-signal-driven
    # rather than CNV-burden-like (a single donor-level chrY loss/gain produces strong signal across
    # all chrY windows), so including them in the top-5% burden statistic causes WT cells with mild
    # sex-chrom signal to spuriously cross the autosomal burden threshold.
    burden_window_mask = np.ones(len(windows), dtype=bool)
    if sex_chrom_set_runtime:
        for w_idx, w in enumerate(windows):
            if str(w.chrom) in sex_chrom_set_runtime:
                burden_window_mask[w_idx] = False

    for state in eligible_states:
        rows = state_index_map[state]
        scores_block = all_scores[rows]
        burden = _cnv_burden(scores_block[:, burden_window_mask], params.burden_quantile)
        candidate_mask = (np.abs(scores_block) >= params.low_threshold).any(axis=1)
        burden_mask = burden >= params.cnv_burden_threshold
        eligible_local = np.flatnonzero(candidate_mask & burden_mask)
        eligible_set = set(int(i) for i in eligible_local)
        for local_index, row in enumerate(rows):
            barcode = str(query.obs_names[row])
            sample = sample_series.iloc[row] if sample_series is not None else ""
            sex_chrom_only = False
            if local_index in eligible_set:
                calls = _call_intervals_for_cell(
                    scores_block[local_index], windows, slices_by_chr,
                    high_threshold=params.high_threshold,
                    low_threshold=params.low_threshold,
                    min_run_windows=params.min_run_windows,
                    min_interval_genes=params.min_interval_genes,
                    min_mean_score=params.min_mean_score,
                )
            elif sex_chrom_set_runtime and candidate_mask[local_index]:
                # Sex-chrom rescue: cell didn't pass the autosomal burden filter but may
                # have a survived sex-chrom interval call (e.g. LOY in an otherwise-quiet
                # cell). Run the interval caller and keep ONLY sex-chrom calls. This avoids
                # flooding outputs with autosomal noise calls from cells that legitimately
                # failed the burden filter.
                all_calls = _call_intervals_for_cell(
                    scores_block[local_index], windows, slices_by_chr,
                    high_threshold=params.high_threshold,
                    low_threshold=params.low_threshold,
                    min_run_windows=params.min_run_windows,
                    min_interval_genes=params.min_interval_genes,
                    min_mean_score=params.min_mean_score,
                )
                calls = [c for c in all_calls if str(c["chr"]) in sex_chrom_set_runtime]
                sex_chrom_only = bool(calls)
            else:
                calls = []
            if calls and all(str(c["chr"]) in sex_chrom_set_runtime for c in calls):
                sex_chrom_only = True

            status = "CNV" if calls else "WT"
            cell_records.append({
                "CellBarcode": barcode,
                "cell_state": state,
                "sample": sample,
                "cnv_status": status,
                "cnv_burden": float(burden[local_index]),
                "baseline_source": "control",
                "n_cnv_intervals": len(calls) if status == "CNV" else 0,
                "sex_chrom_only": sex_chrom_only,
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

    barcode_to_index = {str(b): i for i, b in enumerate(query.obs_names.astype(str))}

    zygosity_thresholds: Optional[Dict[str, float]] = None
    if cached_query_expr is not None and interval_records and params.zygosity_mode == "relative":
        cal_t0 = time.perf_counter()
        wt_indices = np.asarray(
            [barcode_to_index[b] for b, status in zip(
                (r["CellBarcode"] for r in cell_records),
                (r["cnv_status"] for r in cell_records),
            ) if status == "WT" and b in barcode_to_index],
            dtype=np.int64,
        )
        zygosity_thresholds = _calibrate_zygosity_thresholds(
            interval_records, wt_indices, cached_query_expr, control_gene_means, coords,
            seed=params.random_state,
        )
        if zygosity_thresholds is not None:
            LOGGER.info(
                "Calibrated zygosity thresholds in %s: hom_loss<=%.2f het_loss<=%.2f het_gain>=%.2f hom_gain>=%.2f",
                _format_duration(time.perf_counter() - cal_t0),
                zygosity_thresholds["homozygous_loss"], zygosity_thresholds["heterozygous_loss"],
                zygosity_thresholds["heterozygous_gain"], zygosity_thresholds["homozygous_gain"],
            )
        else:
            LOGGER.info("Could not calibrate relative zygosity thresholds; falling back to fixed.")

    sex_chrom_thresholds = (
        {
            "homozygous_loss": float(params.sex_chrom_hom_loss),
            "heterozygous_loss": float(params.sex_chrom_het_loss),
            "heterozygous_gain": float(params.sex_chrom_het_gain),
            "homozygous_gain": float(params.sex_chrom_hom_gain),
        }
        if params.sex_chrom_mode == "absolute_log2"
        else None
    )

    if cached_query_expr is not None and interval_records:
        zyg_t0 = time.perf_counter()
        interval_records = _annotate_intervals_with_zygosity(
            interval_records, cached_query_expr, control_gene_means, coords, barcode_to_index,
            thresholds=zygosity_thresholds,
            sex_chroms=params.sex_chroms if sex_chrom_thresholds is not None else (),
            sex_chrom_thresholds=sex_chrom_thresholds,
        )
        LOGGER.info("Per-cell zygosity annotation done (%s).", _format_duration(time.perf_counter() - zyg_t0))

    cell_df = pd.DataFrame(cell_records)
    if not cell_df.empty:
        cell_df["state_clone_id"] = "WT"
        cell_df["global_clone_id"] = "WT"
    helper_params = FastCNVParams(
        h5ad=params.h5ad,
        gene_coordinates=params.gene_coordinates,
        output_prefix=params.output_prefix,
        state_key=params.state_key,
        sample_key=params.sample_key,
        layer=params.layer,
        input_normalized=params.input_normalized,
        burden_quantile=params.burden_quantile,
        cnv_burden_threshold=params.cnv_burden_threshold,
        max_clones_per_state=params.max_clones_per_state,
        max_global_clones=params.max_global_clones,
        min_clone_cells=params.min_clone_cells,
        clone_similarity_threshold=params.clone_similarity_threshold,
        clone_consensus_fraction=params.clone_consensus_fraction,
        nmf_max_iter=params.nmf_max_iter,
        skip_clones=params.skip_clones,
        skip_pdf=params.skip_pdf,
        random_state=params.random_state,
    )

    state_clone_rows: List[Dict[str, object]] = []
    if not params.skip_clones and not cell_df.empty:
        clones_t0 = time.perf_counter()

        # Cells whose CNV calls are entirely on sex chromosomes (e.g. pure LOY) are excluded
        # from NMF-based autosomal clone discovery — they have no signal in the autosomal
        # feature space, and their inclusion shifts cluster boundaries and dilutes real
        # autosomal CNV clones. They are still tracked, assigned to a synthetic 'sex_chrom'
        # clone, and will be picked up by the clone-consensus interval calling pass.
        sex_chrom_only_col = cell_df["sex_chrom_only"] if "sex_chrom_only" in cell_df.columns else None
        cnv_global_indices = np.asarray(
            [barcode_to_index[b]
             for b, status, sex_only in zip(
                 cell_df["CellBarcode"], cell_df["cnv_status"],
                 sex_chrom_only_col if sex_chrom_only_col is not None else [False]*len(cell_df),
             )
             if status == "CNV" and not bool(sex_only) and b in barcode_to_index],
            dtype=np.int64,
        )
        # Build a mask excluding sex-chrom windows from clone discovery features. Sex-chrom
        # signal (LOY/LOX) is binary at the donor level and dominates NMF clustering, which
        # otherwise scrambles autosomal clone discovery. Sex-chrom CNVs are still surfaced
        # via per-cell intervals and the clone consensus pass on the full score matrix.
        autosomal_mask = np.ones(all_scores.shape[1], dtype=bool)
        if sex_chrom_set_runtime:
            for w_idx, w in enumerate(windows):
                if str(w.chrom) in sex_chrom_set_runtime:
                    autosomal_mask[w_idx] = False

        if cnv_global_indices.size >= max(params.min_clone_cells, 1):
            cnv_scores = all_scores[cnv_global_indices]
            active = (np.abs(cnv_scores) >= params.low_threshold).astype(np.float32)
            active_fraction = np.nan_to_num(active.mean(axis=0), nan=0.0)
            informative = (active_fraction >= params.clone_min_active_fraction) & autosomal_mask
            n_informative = int(informative.sum())
            if n_informative > params.clone_max_features:
                rank_score = np.nan_to_num(np.abs(cnv_scores).mean(axis=0), nan=0.0)
                rank_score[~informative] = -np.inf
                top_idx = np.argpartition(rank_score, -params.clone_max_features)[-params.clone_max_features:]
                informative = np.zeros_like(informative)
                informative[top_idx] = True
            informative_idx = np.flatnonzero(informative).astype(np.int64)
            if informative_idx.size == 0:
                informative_idx = np.flatnonzero(autosomal_mask).astype(np.int64)
            LOGGER.info(
                "Clone feature pruning: %d / %d autosomal windows retained "
                "(active >= %.0f%% across %d CNV cells; cap %d; sex-chrom windows excluded).",
                int(informative_idx.size), int(autosomal_mask.sum()),
                params.clone_min_active_fraction * 100.0, int(cnv_global_indices.size),
                params.clone_max_features,
            )
        else:
            informative_idx = np.flatnonzero(autosomal_mask).astype(np.int64)

        clone_scores = all_scores[:, informative_idx]

        group_columns = ["cell_state"]
        if params.sample_key and "sample" in cell_df.columns:
            group_columns.append("sample")

        group_specs: List[Tuple[int, str, str, np.ndarray, np.ndarray, np.ndarray]] = []
        for state_number, (group_key, state_cell_df) in enumerate(cell_df.groupby(group_columns, sort=False)):
            if isinstance(group_key, tuple):
                state = str(group_key[0])
                sample = str(group_key[1]) if len(group_key) > 1 else ""
            else:
                state = str(group_key)
                sample = ""
            valid_state_df = state_cell_df[state_cell_df["cnv_status"].isin(["CNV", "WT"])]
            if valid_state_df.empty:
                continue
            row_indices = valid_state_df.index.to_numpy()
            score_indices = np.asarray(
                [barcode_to_index[b] for b in valid_state_df["CellBarcode"] if b in barcode_to_index],
                dtype=int,
            )
            if score_indices.size == 0:
                continue
            cnv_mask = valid_state_df["cnv_status"].to_numpy(dtype=object) == "CNV"
            group_specs.append((state_number, state, sample, row_indices, score_indices, cnv_mask))

        def _run_one(spec):
            state_number, state, sample, _row_indices, score_indices, cnv_mask = spec
            labels, summary = _kmeans_state_clones(
                clone_scores[score_indices],
                cnv_mask=cnv_mask,
                max_clones=params.max_clones_per_state,
                min_cells=params.min_clone_cells,
                seed=params.random_state + state_number,
            )
            return state, sample, labels, summary

        try:
            from joblib import Parallel, delayed
            n_jobs = params.n_jobs if params.n_jobs != 0 else 1
            results = Parallel(n_jobs=n_jobs, prefer="threads", verbose=0)(
                delayed(_run_one)(spec) for spec in group_specs
            )
        except ImportError:
            results = [_run_one(spec) for spec in group_specs]

        for spec, (state, sample, labels, summary) in zip(group_specs, results):
            row_indices = spec[3]
            resolved = np.asarray(
                ["WT" if label == "WT" else f"{state}::{sample}::{label}" for label in labels],
                dtype=object,
            )
            cell_df.loc[row_indices, "state_clone_id"] = resolved
            if not summary.empty:
                summary.insert(0, "cell_state", state)
                summary.insert(1, "sample", sample)
                summary["state_clone_id"] = summary["state_clone_id"].map(
                    lambda label: f"{state}::{sample}::{label}"
                )
                state_clone_rows.extend(summary.to_dict("records"))

        _, state_profiles, state_sizes = _clone_profile_frame(cell_df, clone_scores, barcode_to_index, "state_clone_id")
        global_assignments = _merge_global_clones(state_profiles, state_sizes, helper_params)
        cell_df["global_clone_id"] = cell_df["state_clone_id"].map(lambda c: global_assignments.get(str(c), "WT"))
        for row in state_clone_rows:
            row["global_clone_id"] = global_assignments.get(str(row["state_clone_id"]), "WT")

        # Override: sex-chrom-only cells (excluded from NMF) get a synthetic clone label
        # so the consensus interval pass picks up their chrY/chrX signal as a clone-level
        # call. This produces an LOY clone for visualization without disrupting autosomal
        # NMF.
        if "sex_chrom_only" in cell_df.columns:
            sex_chrom_mask = cell_df["sex_chrom_only"].fillna(False).astype(bool)
            if sex_chrom_mask.any():
                cell_df.loc[sex_chrom_mask, "global_clone_id"] = "clone_loy"
                cell_df.loc[sex_chrom_mask, "state_clone_id"] = "loy"

        n_global = len({c for c in cell_df["global_clone_id"] if c != "WT"})
        LOGGER.info(
            "Clone discovery done (%s): %d state-local clones merged into %d global clones.",
            _format_duration(time.perf_counter() - clones_t0),
            len(state_clone_rows),
            n_global,
        )


    n_total_cells = int(cell_df.shape[0]) if not cell_df.empty else 0
    n_cnv_cells = int((cell_df["cnv_status"] == "CNV").sum()) if not cell_df.empty else 0
    cnv_fraction = n_cnv_cells / max(n_total_cells, 1)
    global_low_confidence = (
        not params.skip_clones and not cell_df.empty and cnv_fraction < params.clone_min_cnv_fraction
    )
    if global_low_confidence:
        LOGGER.warning(
            "CNV-cell fraction %.2f%% below threshold %.2f%%; flagging ALL clones as low_confidence.",
            100.0 * cnv_fraction, 100.0 * params.clone_min_cnv_fraction,
        )
    clone_size_counts = (
        cell_df["global_clone_id"].value_counts().to_dict() if "global_clone_id" in cell_df.columns else {}
    )

    def _confidence_for(clone_id: object) -> str:
        cid = str(clone_id)
        if cid == "WT" or cid == "":
            return "wt"
        if global_low_confidence:
            return "low"
        if int(clone_size_counts.get(cid, 0)) < params.clone_min_cells_confident:
            return "low"
        return "high"

    if not cell_df.empty and "global_clone_id" in cell_df.columns:
        cell_df["clone_confidence"] = cell_df["global_clone_id"].map(_confidence_for)
    for row in state_clone_rows:
        row["clone_confidence"] = _confidence_for(row.get("global_clone_id", ""))

    clone_intervals_df = (
        _call_clone_consensus_intervals(cell_df, all_scores, windows, slices_by_chr, barcode_to_index, helper_params)
        if not params.skip_clones and not cell_df.empty
        else pd.DataFrame()
    )

    if cached_query_expr is not None and not clone_intervals_df.empty:
        clone_to_query_rows_zyg: Dict[str, np.ndarray] = {}
        for clone_id in clone_intervals_df["global_clone_id"].astype(str).unique():
            barcodes = cell_df.loc[cell_df["global_clone_id"] == clone_id, "CellBarcode"].astype(str)
            indices = np.asarray([barcode_to_index[b] for b in barcodes if b in barcode_to_index], dtype=np.int64)
            if indices.size > 0:
                clone_to_query_rows_zyg[clone_id] = indices
        if clone_to_query_rows_zyg:
            clone_log2_ratios_full, clone_log2_order = _compute_clone_log2_ratios_from_cache(
                cached_query_expr, control_gene_means, clone_to_query_rows_zyg,
            )
            clone_log2_lookup = {clone_id: clone_log2_ratios_full[i] for i, clone_id in enumerate(clone_log2_order)}
            coords_chr = coords["chr"].astype(str).to_numpy()
            coords_start = coords["start"].to_numpy()
            coords_end = coords["end"].to_numpy()
            mean_log2_list: List[float] = []
            zyg_states: List[str] = []
            for _, interval_row in clone_intervals_df.iterrows():
                clone_id = str(interval_row["global_clone_id"])
                gene_log2 = clone_log2_lookup.get(clone_id)
                if gene_log2 is None:
                    mean_log2_list.append(float("nan"))
                    zyg_states.append("low_signal")
                    continue
                in_interval = (
                    (coords_chr == str(interval_row["chr"])) &
                    (coords_start >= float(interval_row["start"])) &
                    (coords_end <= float(interval_row["end"]))
                )
                if not in_interval.any():
                    mean_log2_list.append(float("nan"))
                    zyg_states.append("low_signal")
                    continue
                mean_log2 = float(np.nanmean(gene_log2[in_interval]))
                mean_log2_list.append(mean_log2)
                interval_chr = str(interval_row["chr"])
                if (
                    sex_chrom_thresholds is not None
                    and interval_chr in params.sex_chroms
                ):
                    zyg_states.append(_zygosity_state(mean_log2, sex_chrom_thresholds))
                else:
                    zyg_states.append(_zygosity_state(mean_log2, None))
            clone_intervals_df = clone_intervals_df.copy()
            clone_intervals_df["mean_log2_ratio"] = mean_log2_list
            clone_intervals_df["zygosity_state"] = zyg_states
            clone_intervals_df["clone_confidence"] = clone_intervals_df["global_clone_id"].astype(str).map(
                lambda cid: _confidence_for(cid)
            )

    cells_path = params.output_prefix.with_suffix(".cnv_cells.tsv")
    intervals_path = params.output_prefix.with_suffix(".cnv_intervals.tsv")
    windows_path = params.output_prefix.with_suffix(".cnv_windows.tsv")
    scores_path = params.output_prefix.with_suffix(".cnv_window_scores.npz")
    state_clones_path = params.output_prefix.with_suffix(".state_clones.tsv")
    clone_intervals_path = params.output_prefix.with_suffix(".clone_intervals.tsv")
    clone_pdf_path = params.output_prefix.with_suffix(".clone_genome.pdf")
    h5ad_path = params.output_prefix.with_suffix(".fastcnv.h5ad")
    sample_sex_path = params.output_prefix.with_suffix(".sample_sex.tsv")

    sex_audit_rows = []
    for sample, row in query_sex_df.iterrows():
        sex_audit_rows.append({"source": "query", "sample": sample, **row.to_dict()})
    for sample, row in control_sex_df.iterrows():
        sex_audit_rows.append({"source": "control", "sample": sample, **row.to_dict()})
    pd.DataFrame(sex_audit_rows).to_csv(sample_sex_path, sep="\t", index=False)

    cell_df.to_csv(cells_path, sep="\t", index=False)
    interval_columns = [
        "CellBarcode", "cell_state", "sample", "call", "chr", "start", "end",
        "start_gene", "end_gene", "n_windows", "n_genes", "mean_score", "max_score", "confidence",
        "mean_log2_ratio", "zygosity_state",
    ]
    pd.DataFrame(interval_records, columns=interval_columns).to_csv(intervals_path, sep="\t", index=False)
    _window_frame(windows).to_csv(windows_path, sep="\t", index=False)
    pd.DataFrame(state_clone_rows).to_csv(state_clones_path, sep="\t", index=False)
    clone_interval_columns = [
        "global_clone_id", "n_cells", "support_fraction", "call", "chr", "start", "end",
        "start_gene", "end_gene", "n_windows", "n_genes", "mean_score", "max_score", "confidence",
        "mean_log2_ratio", "zygosity_state", "clone_confidence",
    ]
    pd.DataFrame(clone_intervals_df, columns=clone_interval_columns).to_csv(clone_intervals_path, sep="\t", index=False)
    np.savez(
        scores_path,
        scores=all_scores,
        cell_barcodes=np.asarray(query.obs_names.astype(str), dtype=object),
        states=state_series.to_numpy(dtype=object),
        window_ids=np.asarray([w.window_id for w in windows], dtype=object),
    )

    pdf_path = None
    if not params.skip_pdf and not cell_df.empty and "global_clone_id" in cell_df.columns:
        pdf_t0 = time.perf_counter()
        clone_to_query_rows: Dict[str, np.ndarray] = {}
        clone_sizes_map: Dict[str, int] = {}
        clone_ids = sorted(
            (c for c in cell_df["global_clone_id"].dropna().unique() if c != "WT"),
            key=lambda c: (0, int(c[5:])) if str(c).startswith("clone") and str(c)[5:].isdigit() else (1, str(c)),
        )
        for clone_id in clone_ids:
            barcodes = cell_df.loc[cell_df["global_clone_id"] == clone_id, "CellBarcode"].astype(str)
            indices = np.asarray([barcode_to_index[b] for b in barcodes if b in barcode_to_index], dtype=np.int64)
            if indices.size == 0:
                continue
            clone_to_query_rows[str(clone_id)] = indices
            clone_sizes_map[str(clone_id)] = int(indices.size)

        if clone_to_query_rows:
            ordered_clone_ids = list(clone_to_query_rows.keys())
            clone_score_profiles = np.full(
                (len(ordered_clone_ids), all_scores.shape[1]), np.nan, dtype=np.float32,
            )
            for i, clone_id in enumerate(ordered_clone_ids):
                rows = clone_to_query_rows[clone_id]
                if rows.size == 0:
                    continue
                clone_score_profiles[i] = np.nanmean(all_scores[rows], axis=0)
            pdf_path = _write_clone_genome_pdf(
                clone_pdf_path,
                coords=coords,
                clone_score_profiles=clone_score_profiles,
                windows=windows,
                clone_order=ordered_clone_ids,
                clone_sizes=clone_sizes_map,
                clone_intervals_df=clone_intervals_df,
                score_y_clip=float(params.pdf_score_y_clip),
                score_threshold=float(params.low_threshold),
            )
        LOGGER.info("Clone PDF written in %s.", _format_duration(time.perf_counter() - pdf_t0))

    outputs = {
        "cells": cells_path,
        "intervals": intervals_path,
        "windows": windows_path,
        "scores": scores_path,
        "state_clones": state_clones_path,
        "clone_intervals": clone_intervals_path,
        "sample_sex": sample_sex_path,
    }
    if pdf_path is not None:
        outputs["clone_pdf"] = pdf_path
    if params.write_h5ad:
        obs_updates = cell_df.set_index("CellBarcode")
        for col in ("cnv_status", "cnv_burden", "state_clone_id", "global_clone_id"):
            query.obs[f"fastcnv_{col}"] = obs_updates.reindex(query.obs_names.astype(str))[col].values
        query.write_h5ad(h5ad_path)
        outputs["h5ad"] = h5ad_path

    LOGGER.info("fastCNV total runtime: %s", _format_duration(time.perf_counter() - run_start))
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
    p.add_argument("--max-clones-per-state", type=int, default=10)
    p.add_argument("--max-global-clones", type=int, default=10)
    p.add_argument("--min-clone-cells", type=int, default=5)
    p.add_argument("--clone-similarity-threshold", type=float, default=0.88)
    p.add_argument("--clone-consensus-fraction", type=float, default=0.45)
    p.add_argument("--nmf-max-iter", type=int, default=100)
    p.add_argument("--clone-min-active-fraction", type=float, default=0.05, help="Window kept for clone discovery if active in >= this fraction of CNV cells.")
    p.add_argument("--clone-max-features", type=int, default=400, help="Hard cap on windows used for clone discovery; top by mean |score| if exceeded.")
    p.add_argument("--clone-min-cnv-fraction", type=float, default=0.015, help="If CNV-cell fraction is below this, ALL clones are flagged low_confidence.")
    p.add_argument("--clone-min-cells-confident", type=int, default=30, help="Clones with fewer than this many cells are flagged low_confidence.")
    p.add_argument("--zygosity-mode", choices=["fixed", "relative"], default="relative", help="Zygosity threshold mode: 'fixed' uses literature defaults; 'relative' calibrates from WT-cell log2 distribution.")
    p.add_argument("--n-jobs", type=int, default=-1, help="Workers for parallel clone discovery; -1 = all cores, 1 = sequential.")
    p.add_argument("--pdf-smooth-genes", type=int, default=50, help="(Legacy) Rolling gene-window size for the older log2-ratio PDF; unused by the score-based PDF.")
    p.add_argument("--pdf-y-clip", type=float, default=1.0, help="(Legacy) Symmetric y-axis limit for the older log2-ratio PDF; unused by the score-based PDF.")
    p.add_argument("--pdf-label-threshold", type=float, default=0.25, help="(Legacy) Driver-gene label threshold for the older PDF; unused by the score-based PDF.")
    p.add_argument("--pdf-score-y-clip", type=float, default=4.0,
                   help="Symmetric y-axis limit for the per-clone CNV-score genome plot. Window scores past this are clipped for display only.")
    p.add_argument("--skip-clones", action="store_true", help="Skip state-local NMF clone discovery and global merge.")
    p.add_argument("--skip-pdf", action="store_true", help="Skip clone-level genome PDF export.")
    p.add_argument(
        "--sex-chrom-mode",
        choices=["off", "absolute_log2"],
        default="absolute_log2",
        help="Score chrX/chrY using absolute log2 fold-change instead of MAD-standardized "
        "residuals. The MAD denominator inflates with donor-level chrY variability and "
        "squashes real LOY signal; absolute log2 recovers it. Use 'off' for legacy behavior.",
    )
    p.add_argument(
        "--sex-chrom-log2-unit",
        type=float,
        default=0.040,
        help="Divisor mapping the mean log1p(CP10K) difference (query - control_state_mean) "
        "into standardized-score units used by the interval caller. Default 0.040 calibrated "
        "to keep healthy-male LOY false-positive rate under ~2%% per donor while retaining "
        "sensitivity to bona fide LOY clones. Lower values (e.g. 0.025) increase sensitivity "
        "at the cost of false positives; higher values (>=0.05) collapse the LOY signal entirely.",
    )
    p.add_argument("--sex-chrom-het-loss", type=float, default=-0.6,
                   help="Mean log2 fold-change (computed against the male-only control baseline "
                        "for chrY) at/below which a sex-chrom interval is called heterozygous_loss.")
    p.add_argument("--sex-chrom-hom-loss", type=float, default=-1.5,
                   help="Mean log2 fold-change at/below which a sex-chrom interval is called "
                        "homozygous_loss.")
    p.add_argument("--sex-chrom-het-gain", type=float, default=0.4)
    p.add_argument("--sex-chrom-hom-gain", type=float, default=0.7)
    p.add_argument("--control-sample-key", default=None,
                   help="Optional obs column on the control h5ad identifying donor/sample. "
                        "Used by automatic sex detection to call each control sample male/female. "
                        "If absent, the entire control is treated as a single sample.")
    p.add_argument("--sex-detection-threshold", type=float, default=SEX_DETECTION_CHRY_PCT_DEFAULT,
                   help="Percent of cells in a sample with >=1 chrY-Y-only marker UMI required "
                        "to call the sample male. Default 5%% reliably separates female samples "
                        "(typically 0.5-2%% chrY+ background dropout) from male samples even "
                        "with severe LOY (a male with ~85%% LOY still has ~15%% chrY+ cells).")
    p.add_argument("--sex-chroms", default="chrY",
                   help="Comma-separated chromosomes scored via absolute log1p-difference rather "
                        "than MAD-standardized residuals. Defaults to chrY only because chrX in "
                        "single-cell data has X-inactivation/dosage-compensation complexities and "
                        "often shows systematic batch effects vs reference. Pass 'chrX,chrY' to "
                        "include chrX (e.g. for monosomy-X detection in a sex-matched setting).")
    p.add_argument("--random-state", type=int, default=0)
    p.add_argument("--write-h5ad", action="store_true")
    p.add_argument("--per-sample", action="store_true",
                   help="Run fastCNV independently for each unique value in --sample-key. "
                        "Each sample gets its own subdirectory under the parent of --output; "
                        "the --output stem becomes the file prefix inside each subdirectory. "
                        "Example: '--output /run/fastCNV --per-sample' produces "
                        "/run/S1/fastCNV.* , /run/S2/fastCNV.* , etc. No cross-sample summary "
                        "files are written; each sample directory is self-contained.")
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
        max_clones_per_state=int(args.max_clones_per_state),
        max_global_clones=int(args.max_global_clones),
        min_clone_cells=int(args.min_clone_cells),
        clone_similarity_threshold=float(args.clone_similarity_threshold),
        clone_consensus_fraction=float(args.clone_consensus_fraction),
        nmf_max_iter=int(args.nmf_max_iter),
        clone_min_active_fraction=float(args.clone_min_active_fraction),
        clone_max_features=int(args.clone_max_features),
        clone_min_cnv_fraction=float(args.clone_min_cnv_fraction),
        clone_min_cells_confident=int(args.clone_min_cells_confident),
        zygosity_mode=str(args.zygosity_mode),
        skip_clones=bool(args.skip_clones),
        skip_pdf=bool(args.skip_pdf),
        random_state=int(args.random_state),
        write_h5ad=bool(args.write_h5ad),
        n_jobs=int(args.n_jobs),
        pdf_smooth_genes=int(args.pdf_smooth_genes),
        pdf_y_clip=float(args.pdf_y_clip),
        pdf_label_threshold=float(args.pdf_label_threshold),
        pdf_score_y_clip=float(args.pdf_score_y_clip),
        sex_chrom_mode=str(args.sex_chrom_mode),
        sex_chrom_log2_unit=float(args.sex_chrom_log2_unit),
        sex_chrom_het_loss=float(args.sex_chrom_het_loss),
        sex_chrom_hom_loss=float(args.sex_chrom_hom_loss),
        sex_chrom_het_gain=float(args.sex_chrom_het_gain),
        sex_chrom_hom_gain=float(args.sex_chrom_hom_gain),
        sex_chroms=tuple(c.strip() for c in str(args.sex_chroms).split(",") if c.strip()),
        sex_detection_threshold_pct=float(args.sex_detection_threshold),
        control_sample_key=args.control_sample_key,
    )


def _run_per_sample(args: argparse.Namespace, base_params: FastParams) -> Dict[str, Path]:
    """Split the query h5ad by --sample-key and run fastCNV independently per sample.

    Each sample gets its own subdirectory under the parent of ``--output``. The
    ``--output`` stem becomes the file prefix inside each subdirectory, e.g. ::

        --output /path/run/fastCNV --sample-key sample_id --per-sample
        ->  /path/run/S1/fastCNV.cnv_cells.tsv
            /path/run/S2/fastCNV.cnv_cells.tsv
            ...

    No cross-sample aggregates are written — each sample's directory is self-contained.
    """
    if not base_params.sample_key:
        raise ValueError("--per-sample requires --sample-key")
    LOGGER.info("Reading combined query h5ad to enumerate samples: %s", base_params.h5ad)
    full = ad.read_h5ad(base_params.h5ad)
    if base_params.sample_key not in full.obs.columns:
        raise ValueError(
            f"--sample-key '{base_params.sample_key}' not in query.obs columns: {list(full.obs.columns)}"
        )
    sample_values = full.obs[base_params.sample_key].astype(str)
    samples = [s for s in pd.unique(sample_values) if s and s.lower() not in ("nan", "none")]
    LOGGER.info("Discovered %d samples: %s", len(samples), samples)

    parent_dir = base_params.output_prefix.parent
    parent_dir.mkdir(parents=True, exist_ok=True)
    file_stem = base_params.output_prefix.name

    aggregated_outputs: Dict[str, Path] = {}

    for sample in samples:
        sample_safe = str(sample).replace("/", "_").replace(" ", "_")
        sample_dir = parent_dir / sample_safe
        sample_dir.mkdir(parents=True, exist_ok=True)
        sample_prefix = sample_dir / file_stem
        sample_h5ad = sample_prefix.with_suffix(".input.h5ad")

        sub = full[sample_values.values == sample].copy()
        LOGGER.info("[%s] %d cells -> %s/", sample, sub.n_obs, sample_dir)
        if sub.n_obs < max(base_params.min_clone_cells, base_params.min_state_cells):
            LOGGER.warning(
                "[%s] only %d cells (< min-state-cells=%d); writing inputs but expect mostly empty outputs.",
                sample, sub.n_obs, base_params.min_state_cells,
            )
        sub.write_h5ad(sample_h5ad)

        sample_params = replace(
            base_params,
            h5ad=sample_h5ad,
            output_prefix=sample_prefix,
        )
        try:
            outs = run_fast(sample_params)
        except Exception:
            LOGGER.exception("[%s] run_fast failed; continuing with remaining samples", sample)
            continue
        for k, v in outs.items():
            aggregated_outputs[f"{sample}/{k}"] = v

    return aggregated_outputs


def main(argv: Optional[Sequence[str]] = None) -> Dict[str, Path]:
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s | %(message)s")
    base_params = params_from_args(args)
    if getattr(args, "per_sample", False):
        outputs = _run_per_sample(args, base_params)
    else:
        outputs = run_fast(base_params)
    for name, path in outputs.items():
        LOGGER.info("%s: %s", name, path)
    return outputs


if __name__ == "__main__":
    main()
