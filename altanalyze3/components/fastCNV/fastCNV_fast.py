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


LOGGER = logging.getLogger("fastCNV.fast")


CANCER_DRIVER_GENES: Tuple[str, ...] = (
    "MCL1", "MYCN", "GLI2", "PDGFRA", "FGFR3", "KIT", "EGFR", "CDK6", "MET",
    "FGFR1", "MYC", "CDKN2A", "CDKN2B", "PTCH1", "JAK2", "PTEN", "FAS",
    "CCND1", "ATM", "KRAS", "CDK4", "MDM2", "RB1", "BRCA2", "FOXO1", "NF1",
    "TP53", "ERBB2", "BRCA1", "STK11", "CCNE1", "AURKA", "RUNX1", "AKT1",
    "BCL2", "SMAD4", "APC", "NOTCH1", "FBXW7", "ARID1A", "TET2", "DNMT3A",
    "ASXL1", "EZH2", "SF3B1", "U2AF1", "SRSF2", "IDH1", "IDH2", "FLT3",
    "NPM1", "WT1", "GATA2", "SETBP1", "BCOR", "STAG2", "PHF6", "RAD21",
    "CALR", "MPL",
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
    skip_clones: bool = False
    skip_pdf: bool = False
    write_h5ad: bool = False
    random_state: int = 0
    n_jobs: int = -1
    pdf_smooth_genes: int = 50
    pdf_y_clip: float = 1.0
    pdf_label_threshold: float = 0.25


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
    log2_ratios: np.ndarray,
    clone_order: Sequence[str],
    clone_sizes: Dict[str, int],
    clone_intervals_df: pd.DataFrame,
    driver_genes: Sequence[str],
    y_clip: float,
    label_threshold: float = 0.25,
) -> Optional[Path]:
    if log2_ratios.shape[0] == 0 or coords.empty:
        return None
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LinearSegmentedColormap

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    x_positions, chr_order, chr_bounds = _genome_x_positions(coords)
    chr_centers = [(s + e) / 2.0 for s, e in chr_bounds]
    chr_labels = [c.replace("chr", "") for c in chr_order]
    chr_to_bounds = dict(zip(chr_order, chr_bounds))

    cmap = LinearSegmentedColormap.from_list(
        "fastcnv_diverging",
        [(0.0, "#c8252c"), (0.5, "#cccccc"), (1.0, "#1f863d")],
    )

    driver_lookup = {g.upper(): g for g in driver_genes}
    coord_gene_upper = coords["gene"].astype(str).str.upper().to_numpy()
    coord_chr = coords["chr"].astype(str).to_numpy()
    is_driver = np.isin(coord_gene_upper, list(driver_lookup.keys()))
    driver_indices = np.flatnonzero(is_driver)
    driver_gene_names = coord_gene_upper[driver_indices]

    intervals_by_clone: Dict[str, pd.DataFrame] = {}
    if not clone_intervals_df.empty:
        intervals_by_clone = {
            str(cid): grp for cid, grp in clone_intervals_df.groupby("global_clone_id", sort=False)
        }

    with PdfPages(path) as pdf:
        for clone_idx, clone_id in enumerate(clone_order):
            clone_ratios = log2_ratios[clone_idx]
            ratios_clipped = np.clip(np.nan_to_num(clone_ratios, nan=0.0), -y_clip, y_clip)

            # Page A: genome-wide scatter
            fig, ax = plt.subplots(figsize=(14, 4.5))
            scatter = ax.scatter(
                x_positions, ratios_clipped,
                c=ratios_clipped, cmap=cmap, vmin=-y_clip, vmax=y_clip,
                s=4, linewidths=0, alpha=0.85, rasterized=True,
            )
            ax.axhline(0.0, color="#888888", linewidth=0.5)
            for chrom in chr_order:
                start_x, _ = chr_to_bounds[chrom]
                ax.axvline(start_x, color="#cccccc", linewidth=0.4, linestyle="--")
            ax.set_xticks(chr_centers)
            ax.set_xticklabels(chr_labels, fontsize=8)
            ax.set_xlim(chr_bounds[0][0], chr_bounds[-1][1])
            ax.set_ylim(-y_clip - 0.05, y_clip + 0.05)
            ax.set_ylabel("log2(clone / control)")
            ax.set_xlabel("Chromosome")
            n_cells = int(clone_sizes.get(str(clone_id), 0))
            ax.set_title(f"{clone_id}  (n={n_cells} cells)")

            # Overlay called intervals as dark blue horizontal segments at the interval mean ratio
            clone_intervals = intervals_by_clone.get(str(clone_id))
            if clone_intervals is not None and not clone_intervals.empty:
                for _, interval_row in clone_intervals.iterrows():
                    chrom = str(interval_row["chr"])
                    if chrom not in chr_to_bounds:
                        continue
                    chr_x_start, _ = chr_to_bounds[chrom]
                    chr_min_start = float(coords.loc[coords["chr"] == chrom, "start"].min())
                    seg_x0 = chr_x_start + (float(interval_row["start"]) - chr_min_start)
                    seg_x1 = chr_x_start + (float(interval_row["end"]) - chr_min_start)
                    in_chr = (coord_chr == chrom)
                    in_segment = in_chr & (
                        (coords["start"].to_numpy() >= float(interval_row["start"])) &
                        (coords["end"].to_numpy() <= float(interval_row["end"]))
                    )
                    if not in_segment.any():
                        continue
                    seg_y = float(np.nanmean(clone_ratios[in_segment]))
                    ax.plot([seg_x0, seg_x1], [seg_y, seg_y], color="#1f3a93", linewidth=2.0, solid_capstyle="butt")

            # Driver gene labels
            if driver_indices.size:
                for di, gene in zip(driver_indices, driver_gene_names):
                    y_val = ratios_clipped[di]
                    if abs(y_val) < label_threshold:
                        continue
                    ax.annotate(
                        driver_lookup.get(gene, gene),
                        xy=(x_positions[di], y_val),
                        xytext=(0, 8 if y_val >= 0 else -10),
                        textcoords="offset points",
                        fontsize=6, color="#1f3a93", ha="center",
                        va="bottom" if y_val >= 0 else "top",
                    )
                    ax.scatter([x_positions[di]], [y_val], s=18, facecolor="#1f3a93", edgecolor="white", linewidth=0.4, zorder=5)
            fig.tight_layout()
            pdf.savefig(fig, dpi=200)
            plt.close(fig)

            # Page B: driver-gene strip plot
            if driver_indices.size:
                fig, ax = plt.subplots(figsize=(14, 4.5))
                ordered_drivers = [driver_lookup.get(coord_gene_upper[i], coord_gene_upper[i]) for i in driver_indices]
                strip_x = np.arange(driver_indices.size, dtype=np.float64)
                strip_y = ratios_clipped[driver_indices]
                ax.scatter(strip_x, strip_y, c=strip_y, cmap=cmap, vmin=-y_clip, vmax=y_clip,
                           s=40, linewidths=0.4, edgecolor="white", rasterized=True)
                ax.axhline(0.0, color="#888888", linewidth=0.5)
                for x in strip_x:
                    ax.axvline(x, color="#eeeeee", linewidth=0.5, linestyle="-", zorder=0)
                ax.set_xticks(strip_x)
                ax.set_xticklabels(ordered_drivers, rotation=90, fontsize=7)
                ax.set_ylabel("log2(clone / control)")
                ax.set_ylim(-y_clip - 0.05, y_clip + 0.05)
                ax.set_xlim(-0.6, driver_indices.size - 0.4)
                ax.set_title(f"{clone_id}  driver-gene log2 ratios  (n={n_cells} cells)")
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
        ax.set_title("fastCNV (fast) clone-level genome calls")
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

    all_scores = np.full((query.n_obs, len(windows)), np.nan, dtype=np.float32)
    control_gene_means = np.full(coords.shape[0], np.nan, dtype=np.float32)
    coords_var_lookup = pd.Series(
        np.arange(coords.shape[0], dtype=np.int64), index=coords["var_index"].to_numpy()
    )
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
        coord_rows_for_chr = coords_var_lookup.reindex(query_cols).to_numpy()
        gene_chunk_means = np.nanmean(control_expr_full, axis=0).astype(np.float32)
        control_gene_means[coord_rows_for_chr] = gene_chunk_means
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

    cell_df = pd.DataFrame(cell_records)
    if not cell_df.empty:
        cell_df["state_clone_id"] = "WT"
        cell_df["global_clone_id"] = "WT"
    barcode_to_index = {str(b): i for i, b in enumerate(query.obs_names.astype(str))}
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

        cnv_global_indices = np.asarray(
            [barcode_to_index[b] for b, status in zip(cell_df["CellBarcode"], cell_df["cnv_status"])
             if status == "CNV" and b in barcode_to_index],
            dtype=np.int64,
        )
        if cnv_global_indices.size >= max(params.min_clone_cells, 1):
            cnv_scores = all_scores[cnv_global_indices]
            active = (np.abs(cnv_scores) >= params.low_threshold).astype(np.float32)
            active_fraction = np.nan_to_num(active.mean(axis=0), nan=0.0)
            informative = active_fraction >= params.clone_min_active_fraction
            n_informative = int(informative.sum())
            if n_informative > params.clone_max_features:
                rank_score = np.nan_to_num(np.abs(cnv_scores).mean(axis=0), nan=0.0)
                rank_score[~informative] = -np.inf
                top_idx = np.argpartition(rank_score, -params.clone_max_features)[-params.clone_max_features:]
                informative = np.zeros_like(informative)
                informative[top_idx] = True
            informative_idx = np.flatnonzero(informative).astype(np.int64)
            if informative_idx.size == 0:
                informative_idx = np.arange(all_scores.shape[1], dtype=np.int64)
            LOGGER.info(
                "Clone feature pruning: %d / %d windows retained (active >= %.0f%% across %d CNV cells; cap %d).",
                int(informative_idx.size), int(all_scores.shape[1]),
                params.clone_min_active_fraction * 100.0, int(cnv_global_indices.size),
                params.clone_max_features,
            )
        else:
            informative_idx = np.arange(all_scores.shape[1], dtype=np.int64)

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
            anchor_mask = np.zeros(score_indices.size, dtype=bool)
            labels, summary = _discover_state_clones(
                clone_scores[score_indices],
                cnv_mask=cnv_mask,
                anchor_mask=anchor_mask,
                params=helper_params,
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
        n_global = len({c for c in cell_df["global_clone_id"] if c != "WT"})
        LOGGER.info(
            "Clone discovery done (%s): %d state-local clones merged into %d global clones.",
            _format_duration(time.perf_counter() - clones_t0),
            len(state_clone_rows),
            n_global,
        )

    clone_intervals_df = (
        _call_clone_consensus_intervals(cell_df, all_scores, windows, slices_by_chr, barcode_to_index, helper_params)
        if not params.skip_clones and not cell_df.empty
        else pd.DataFrame()
    )

    cells_path = params.output_prefix.with_suffix(".cnv_cells.tsv")
    intervals_path = params.output_prefix.with_suffix(".cnv_intervals.tsv")
    windows_path = params.output_prefix.with_suffix(".cnv_windows.tsv")
    scores_path = params.output_prefix.with_suffix(".cnv_window_scores.npz")
    state_clones_path = params.output_prefix.with_suffix(".state_clones.tsv")
    clone_intervals_path = params.output_prefix.with_suffix(".clone_intervals.tsv")
    clone_pdf_path = params.output_prefix.with_suffix(".clone_genome.pdf")
    h5ad_path = params.output_prefix.with_suffix(".fastcnv.h5ad")

    cell_df.to_csv(cells_path, sep="\t", index=False)
    interval_columns = [
        "CellBarcode", "cell_state", "sample", "call", "chr", "start", "end",
        "start_gene", "end_gene", "n_windows", "n_genes", "mean_score", "max_score", "confidence",
    ]
    pd.DataFrame(interval_records, columns=interval_columns).to_csv(intervals_path, sep="\t", index=False)
    _window_frame(windows).to_csv(windows_path, sep="\t", index=False)
    pd.DataFrame(state_clone_rows).to_csv(state_clones_path, sep="\t", index=False)
    clone_interval_columns = [
        "global_clone_id", "n_cells", "support_fraction", "call", "chr", "start", "end",
        "start_gene", "end_gene", "n_windows", "n_genes", "mean_score", "max_score", "confidence",
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
            log2_ratios, ordered_clone_ids = _compute_clone_log2_ratios(
                query_matrix, query_lib, coords, control_gene_means,
                clone_to_query_rows, params.input_normalized,
            )
            log2_ratios = _smooth_and_center_log2_ratios(log2_ratios, coords, params.pdf_smooth_genes)
            pdf_path = _write_clone_genome_pdf(
                clone_pdf_path,
                coords=coords,
                log2_ratios=log2_ratios,
                clone_order=ordered_clone_ids,
                clone_sizes=clone_sizes_map,
                clone_intervals_df=clone_intervals_df,
                driver_genes=CANCER_DRIVER_GENES,
                y_clip=params.pdf_y_clip,
                label_threshold=params.pdf_label_threshold,
            )
        LOGGER.info("Clone PDF written in %s.", _format_duration(time.perf_counter() - pdf_t0))

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
        obs_updates = cell_df.set_index("CellBarcode")
        for col in ("cnv_status", "cnv_burden", "state_clone_id", "global_clone_id"):
            query.obs[f"fastcnv_{col}"] = obs_updates.reindex(query.obs_names.astype(str))[col].values
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
    p.add_argument("--max-clones-per-state", type=int, default=10)
    p.add_argument("--max-global-clones", type=int, default=10)
    p.add_argument("--min-clone-cells", type=int, default=5)
    p.add_argument("--clone-similarity-threshold", type=float, default=0.88)
    p.add_argument("--clone-consensus-fraction", type=float, default=0.45)
    p.add_argument("--nmf-max-iter", type=int, default=100)
    p.add_argument("--clone-min-active-fraction", type=float, default=0.05, help="Window kept for clone discovery if active in >= this fraction of CNV cells.")
    p.add_argument("--clone-max-features", type=int, default=400, help="Hard cap on windows used for clone discovery; top by mean |score| if exceeded.")
    p.add_argument("--n-jobs", type=int, default=-1, help="Workers for parallel clone discovery; -1 = all cores, 1 = sequential.")
    p.add_argument("--pdf-smooth-genes", type=int, default=50, help="Rolling gene-window size for PDF log2-ratio smoothing per chromosome.")
    p.add_argument("--pdf-y-clip", type=float, default=1.0, help="Symmetric y-axis limit for PDF panels.")
    p.add_argument("--pdf-label-threshold", type=float, default=0.25, help="Driver-gene label hidden when |smoothed log2 ratio| below this.")
    p.add_argument("--skip-clones", action="store_true", help="Skip state-local NMF clone discovery and global merge.")
    p.add_argument("--skip-pdf", action="store_true", help="Skip clone-level genome PDF export.")
    p.add_argument("--random-state", type=int, default=0)
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
        max_clones_per_state=int(args.max_clones_per_state),
        max_global_clones=int(args.max_global_clones),
        min_clone_cells=int(args.min_clone_cells),
        clone_similarity_threshold=float(args.clone_similarity_threshold),
        clone_consensus_fraction=float(args.clone_consensus_fraction),
        nmf_max_iter=int(args.nmf_max_iter),
        clone_min_active_fraction=float(args.clone_min_active_fraction),
        clone_max_features=int(args.clone_max_features),
        skip_clones=bool(args.skip_clones),
        skip_pdf=bool(args.skip_pdf),
        random_state=int(args.random_state),
        write_h5ad=bool(args.write_h5ad),
        n_jobs=int(args.n_jobs),
        pdf_smooth_genes=int(args.pdf_smooth_genes),
        pdf_y_clip=float(args.pdf_y_clip),
        pdf_label_threshold=float(args.pdf_label_threshold),
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
