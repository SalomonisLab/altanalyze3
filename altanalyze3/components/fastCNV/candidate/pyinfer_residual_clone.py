"""pyInferCNV-residual candidate clone caller.

This is a separate experimental fastCNV candidate. It deliberately does not
modify production fastCNV. The purpose is to test whether fastCNV's missed
SamplePre34 clones are recoverable when the clone geometry is computed with the
same unsupervised residual model used by the pyInferCNV positive control.

The caller remains region-generic:

* compute pyInferCNV residuals from the query count matrix;
* aggregate residuals by chromosome;
* within each cell state, identify separated altered components;
* carry validated fastCNV LOY clone labels forward from a fastCNV prefix.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


@dataclass(frozen=True)
class ResidualCloneParams:
    h5ad: Path
    fastcnv_prefix: Path
    output_prefix: Path
    state_key: str
    layer: str = "auto"
    min_state_cells: int = 60
    min_component_cells: int = 20
    min_component_fraction: float = 0.03
    min_abs_region_score: float = 0.05
    min_separation_mad: float = 1.5
    exclude_chromosomes: tuple[str, ...] = ("chrX", "chrY", "chrM")
    keep_loy: bool = True
    pyinfercnv_path: Optional[Path] = Path("/Users/saljh8/Documents/GitHub/pyInferCNV")
    write_heatmap: bool = True
    heatmap_max_cells: int = 20000
    heatmap_bin_genes: int = 50
    heatmap_min_chrom_bins: int = 25
    heatmap_filter_threshold: float = 0.03


def _matrix(adata: ad.AnnData, layer: str):
    if layer == "auto":
        return adata.layers["counts"] if "counts" in adata.layers else adata.X
    if layer == "X":
        return adata.X
    return adata.layers[layer]


def _prefix_path(prefix: Path, suffix: str) -> Path:
    return Path(f"{prefix}.{suffix}")


def _mad(values: np.ndarray) -> float:
    vals = values[np.isfinite(values)]
    if vals.size == 0:
        return 1.0
    med = float(np.median(vals))
    return max(1.4826 * float(np.median(np.abs(vals - med))), 1e-6)


def _two_component(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    finite = np.isfinite(values)
    labels = np.zeros(values.size, dtype=np.int8)
    if finite.sum() < 4:
        return labels, np.array([np.nan, np.nan])
    x = values[finite]
    centers = np.asarray(np.quantile(x, [0.25, 0.75]), dtype=float)
    if abs(centers[1] - centers[0]) < 1e-8:
        return labels, centers
    lab = np.zeros(x.size, dtype=bool)
    for _ in range(50):
        lab = np.abs(x - centers[1]) < np.abs(x - centers[0])
        new = np.array([
            float(np.mean(x[~lab])) if np.any(~lab) else centers[0],
            float(np.mean(x[lab])) if np.any(lab) else centers[1],
        ])
        if np.max(np.abs(new - centers)) < 1e-6:
            centers = new
            break
        centers = new
    if centers[0] > centers[1]:
        centers = centers[::-1]
        lab = ~lab
    labels[np.flatnonzero(finite)] = lab.astype(np.int8)
    return labels, centers


def _load_fastcnv_cells(prefix: Path, barcodes: np.ndarray) -> pd.DataFrame:
    path = _prefix_path(prefix, "cnv_cells.tsv")
    if not path.exists():
        return pd.DataFrame(index=barcodes)
    frame = pd.read_csv(path, sep="\t")
    frame["CellBarcode"] = frame["CellBarcode"].astype(str)
    return frame.set_index("CellBarcode").reindex(barcodes)


def _compute_residuals(params: ResidualCloneParams) -> tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    try:
        from pyinfercnv import GENE_POS_GENCODE_V32
        from pyinfercnv.h5ad_io import load_gene_positions, build_infercnv_object
        from pyinfercnv.run import run
    except ModuleNotFoundError:
        if not params.pyinfercnv_path or not params.pyinfercnv_path.exists():
            raise
        sys.path.insert(0, str(params.pyinfercnv_path))
        from pyinfercnv import GENE_POS_GENCODE_V32
        from pyinfercnv.h5ad_io import load_gene_positions, build_infercnv_object
        from pyinfercnv.run import run

    adata = ad.read_h5ad(params.h5ad)
    matrix = _matrix(adata, params.layer)
    matrix = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)
    gene_positions = load_gene_positions(GENE_POS_GENCODE_V32)
    obj = build_infercnv_object(
        matrix,
        list(adata.var_names),
        gene_positions,
        cell_labels=None,
        chr_exclude=params.exclude_chromosomes,
        prefilter_min_cells=3,
    )
    run(obj, cutoff=0.1, denoise=True, verbose=False)
    residual = np.nan_to_num(obj.expr_data, nan=1.0).T.astype(np.float32)
    chrom = np.asarray(obj.gene_chr).astype(str)
    genes = np.asarray(list(obj.gene_order.index)).astype(str)
    obs = adata.obs.copy()
    obs.index = obs.index.astype(str)
    return residual, chrom, genes, obs


def _clone_sort_key(value: object) -> tuple[int, str]:
    text = str(value)
    if text == "clone_loy":
        return (0, "000000")
    if text.startswith("clone") and text[5:].isdigit():
        return (1, f"{int(text[5:]):06d}")
    if text == "WT":
        return (9, text)
    return (2, text)


def _binned_residual_matrix(
    residual: np.ndarray,
    chrom: np.ndarray,
    params: ResidualCloneParams,
) -> tuple[np.ndarray, list[str], list[int]]:
    blocks: list[np.ndarray] = []
    labels: list[str] = []
    boundaries = [0]
    for chrom_name in pd.unique(chrom):
        idx = np.flatnonzero(chrom == chrom_name)
        if idx.size == 0:
            continue
        labels.append(str(chrom_name))
        bin_size = max(int(params.heatmap_bin_genes), 1)
        chrom_blocks = []
        for start in range(0, idx.size, bin_size):
            sub = idx[start:start + bin_size]
            chrom_blocks.append(np.nanmean(residual[:, sub] - 1.0, axis=1))
        block = np.vstack(chrom_blocks).T.astype(np.float32)
        if 0 < block.shape[1] < params.heatmap_min_chrom_bins:
            take = np.linspace(0, block.shape[1] - 1, int(params.heatmap_min_chrom_bins)).round().astype(np.int64)
            block = block[:, take]
        blocks.append(block)
        boundaries.append(boundaries[-1] + block.shape[1])
    if not blocks:
        return np.zeros((residual.shape[0], 0), dtype=np.float32), [], [0]
    return np.hstack(blocks).astype(np.float32), labels, boundaries


def _write_residual_heatmaps(
    params: ResidualCloneParams,
    residual: np.ndarray,
    chrom: np.ndarray,
    cell_df: pd.DataFrame,
) -> Dict[str, Path]:
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, ListedColormap

    matrix, chrom_labels, boundaries = _binned_residual_matrix(residual, chrom, params)
    if matrix.shape[1] == 0:
        return {}
    cells = cell_df.copy()
    cells["_row"] = np.arange(cells.shape[0], dtype=np.int64)
    cells["_clone_sort"] = cells["global_clone_id"].map(_clone_sort_key)
    if params.heatmap_max_cells and matrix.shape[0] > params.heatmap_max_cells:
        parts = []
        rng = np.random.default_rng(0)
        for _, group in cells.groupby(["global_clone_id", "cell_state"], sort=False):
            quota = max(1, int(np.floor(params.heatmap_max_cells * len(group) / len(cells))))
            if str(group["global_clone_id"].iloc[0]) != "WT":
                quota = max(quota, min(len(group), 50))
            if quota >= len(group):
                parts.append(group)
            else:
                parts.append(group.loc[np.sort(rng.choice(group.index.to_numpy(), size=quota, replace=False))])
        cells = pd.concat(parts, axis=0)
    cells = cells.sort_values(["_clone_sort", "cell_state", "CellBarcode"], kind="mergesort")
    row_idx = cells["_row"].to_numpy(dtype=np.int64)
    display = matrix[row_idx]
    display = display - np.nanmedian(display, axis=1, keepdims=True)
    display = np.nan_to_num(display, nan=0.0, posinf=0.0, neginf=0.0)
    vmax = max(float(np.nanquantile(np.abs(display), 0.995)), 0.05)
    display = np.clip(display, -vmax, vmax)

    clones = cells["global_clone_id"].astype(str).to_numpy()
    clone_order = sorted(pd.unique(clones), key=_clone_sort_key)
    palette = [
        "#4c78a8", "#f58518", "#54a24b", "#e45756", "#72b7b2",
        "#b279a2", "#ff9da6", "#9d755d", "#bab0ac", "#59a14f",
        "#edc948", "#af7aa1", "#76b7b2", "#d37295",
    ]
    colors = {clone: palette[i % len(palette)] for i, clone in enumerate(clone_order)}
    colors["WT"] = "#d9d9d9"
    colors["clone_loy"] = "#2458a6"
    strip = np.array([clone_order.index(c) for c in clones], dtype=np.int16)[:, None]
    strip_cmap = ListedColormap([colors[c] for c in clone_order])
    heat_cmap = LinearSegmentedColormap.from_list(
        "pyinfer_residual",
        ["#1f4fa3", "#597fc5", "#f7f7f7", "#d95f5f", "#a51f1f"],
    )
    centers = [(boundaries[i] + boundaries[i + 1]) / 2.0 for i in range(len(chrom_labels))]

    def draw(data: np.ndarray, title: str, png: Path, pdf: Path) -> None:
        fig_height = min(12.0, max(5.0, 0.00055 * data.shape[0] + 4.5))
        fig = plt.figure(figsize=(13.2, fig_height))
        gs = fig.add_gridspec(1, 3, width_ratios=[0.018, 1.0, 0.018], wspace=0.012)
        ax_strip = fig.add_subplot(gs[0, 0])
        ax = fig.add_subplot(gs[0, 1], sharey=ax_strip)
        cax = fig.add_subplot(gs[0, 2])
        ax_strip.imshow(strip, aspect="auto", interpolation="none", cmap=strip_cmap, origin="upper")
        ax_strip.set_axis_off()
        im = ax.imshow(data, aspect="auto", interpolation="none", cmap=heat_cmap, vmin=-vmax, vmax=vmax, origin="upper", rasterized=True)
        for boundary in boundaries[1:-1]:
            ax.axvline(boundary - 0.5, color="#202020", linewidth=0.35)
        ax.set_xticks(centers)
        ax.set_xticklabels([c.replace("chr", "") for c in chrom_labels], fontsize=7)
        ax.set_yticks([])
        ax.set_xlabel("chromosome")
        ax.set_title(title, fontsize=10)
        handles = [plt.Rectangle((0, 0), 1, 1, color=colors[c]) for c in clone_order]
        labels = [f"{c} ({int((clones == c).sum()):,})" for c in clone_order]
        ax.legend(handles, labels, title="candidate clone", loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=7, title_fontsize=8, frameon=False)
        fig.colorbar(im, cax=cax).set_label("pyInferCNV residual - 1", fontsize=8)
        fig.savefig(png, dpi=150, bbox_inches="tight")
        fig.savefig(pdf, bbox_inches="tight")
        plt.close(fig)

    primary_png = _prefix_path(params.output_prefix, "candidate_infercnv_heatmap.png")
    primary_pdf = _prefix_path(params.output_prefix, "candidate_infercnv_heatmap.pdf")
    filtered_png = _prefix_path(params.output_prefix, "candidate_infercnv_heatmap_filtered.png")
    filtered_pdf = _prefix_path(params.output_prefix, "candidate_infercnv_heatmap_filtered.pdf")
    draw(display, "pyInferCNV-residual candidate heatmap", primary_png, primary_pdf)
    filtered = np.where(np.abs(display) < float(params.heatmap_filter_threshold), 0.0, display).astype(np.float32)
    draw(filtered, f"pyInferCNV-residual candidate heatmap, filtered |residual|<{params.heatmap_filter_threshold:g}", filtered_png, filtered_pdf)
    return {
        "heatmap_png": primary_png,
        "heatmap_pdf": primary_pdf,
        "heatmap_filtered_png": filtered_png,
        "heatmap_filtered_pdf": filtered_pdf,
    }


def run_candidate(params: ResidualCloneParams) -> Dict[str, Path]:
    residual, chrom, genes, obs = _compute_residuals(params)
    barcodes = obs.index.astype(str).to_numpy()
    fast_cells = _load_fastcnv_cells(params.fastcnv_prefix, barcodes)
    states = obs[params.state_key].astype(str).to_numpy()
    samples = (
        fast_cells["sample"].fillna("").astype(str).to_numpy()
        if "sample" in fast_cells.columns else np.array([""] * len(barcodes), dtype=object)
    )
    calls: list[dict] = []
    chroms = [c for c in pd.unique(chrom) if c not in set(params.exclude_chromosomes)]
    for region in chroms:
        region_mask = chrom == region
        if int(region_mask.sum()) < 10:
            continue
        values = np.nanmean(residual[:, region_mask], axis=1) - 1.0
        start_gene = str(genes[np.flatnonzero(region_mask)[0]])
        end_gene = str(genes[np.flatnonzero(region_mask)[-1]])
        for state in pd.unique(states):
            idx = np.flatnonzero(states == state)
            if idx.size < params.min_state_cells:
                continue
            x = values[idx]
            labels, centers = _two_component(x)
            if not np.all(np.isfinite(centers)):
                continue
            median = float(np.nanmedian(x))
            neutral = int(np.argmin(np.abs(centers - median)))
            altered = 1 - neutral
            altered_local = np.flatnonzero(labels == altered)
            if altered_local.size < params.min_component_cells:
                continue
            if altered_local.size / idx.size < params.min_component_fraction:
                continue
            alt_center = float(centers[altered])
            neutral_center = float(centers[neutral])
            sep = abs(alt_center - neutral_center) / _mad(x)
            if abs(alt_center) < params.min_abs_region_score or sep < params.min_separation_mad:
                continue
            call = "gain" if alt_center > neutral_center else "loss"
            for row in idx[altered_local]:
                calls.append(
                    {
                        "CellBarcode": barcodes[row],
                        "cell_state": state,
                        "sample": samples[row],
                        "call": call,
                        "chr": region,
                        "start": 0,
                        "end": 0,
                        "start_gene": start_gene,
                        "end_gene": end_gene,
                        "region": region,
                        "n_windows": 1,
                        "n_genes": int(region_mask.sum()),
                        "mean_score": alt_center,
                        "max_score": abs(alt_center),
                        "confidence": sep,
                        "mean_log2": alt_center,
                        "boundary_source": "pyinfer_residual_component",
                    }
                )
    region_df = pd.DataFrame(calls)
    cell_df = pd.DataFrame({"CellBarcode": barcodes, "cell_state": states, "sample": samples})
    cell_df["cnv_status"] = "WT"
    cell_df["n_cnv_intervals"] = 0
    cell_df["state_clone_id"] = "WT"
    cell_df["global_clone_id"] = "WT"
    if not region_df.empty:
        sig = (
            region_df.assign(sig=region_df["chr"].astype(str) + ":" + region_df["call"].astype(str))
            .groupby("CellBarcode")["sig"]
            .agg(lambda s: "|".join(sorted(set(s))))
        )
        clone_names = {v: f"clone{i + 1}" for i, v in enumerate(sorted(sig.unique()))}
        nint = region_df.groupby("CellBarcode").size()
        cell_df = cell_df.set_index("CellBarcode")
        cell_df.loc[sig.index, "cnv_status"] = "CNV"
        cell_df.loc[sig.index, "n_cnv_intervals"] = nint.astype(int)
        cell_df.loc[sig.index, "state_clone_id"] = sig.map(clone_names)
        cell_df.loc[sig.index, "global_clone_id"] = sig.map(clone_names)
        cell_df = cell_df.reset_index()

    if params.keep_loy and "global_clone_id" in fast_cells.columns:
        loy_mask = fast_cells["global_clone_id"].astype(str).eq("clone_loy").fillna(False).to_numpy()
        if loy_mask.any():
            cell_df.loc[loy_mask, "cnv_status"] = "CNV"
            cell_df.loc[loy_mask, "state_clone_id"] = "clone_loy"
            cell_df.loc[loy_mask, "global_clone_id"] = "clone_loy"
            cell_df.loc[loy_mask, "n_cnv_intervals"] = np.maximum(
                cell_df.loc[loy_mask, "n_cnv_intervals"].astype(int), 1
            )
            loy_df = pd.DataFrame(
                {
                    "CellBarcode": barcodes[loy_mask],
                    "cell_state": states[loy_mask],
                    "sample": samples[loy_mask],
                    "call": "loss",
                    "chr": "chrY",
                    "start": 0,
                    "end": 0,
                    "start_gene": "",
                    "end_gene": "",
                    "region": "chrY",
                    "n_windows": 1,
                    "n_genes": 1,
                    "mean_score": np.nan,
                    "max_score": np.nan,
                    "confidence": np.nan,
                    "mean_log2": np.nan,
                    "boundary_source": "fastcnv_loy",
                }
            )
            region_df = pd.concat([region_df, loy_df], ignore_index=True)

    if region_df.empty:
        region_df = pd.DataFrame(
            columns=[
                "CellBarcode", "cell_state", "sample", "call", "chr", "start", "end",
                "start_gene", "end_gene", "region", "n_windows", "n_genes",
                "mean_score", "max_score", "confidence", "mean_log2", "boundary_source",
            ]
        )
    summary = (
        region_df.groupby(["chr", "region", "call"], dropna=False)
        .agg(n_cells=("CellBarcode", "nunique"), n_states=("cell_state", "nunique"), mean_score=("mean_score", "mean"))
        .reset_index()
        if not region_df.empty else pd.DataFrame(columns=["chr", "region", "call", "n_cells", "n_states", "mean_score"])
    )
    params.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    cells_path = _prefix_path(params.output_prefix, "candidate_cnv_cells.tsv")
    regions_path = _prefix_path(params.output_prefix, "candidate_cnv_regions.tsv")
    summary_path = _prefix_path(params.output_prefix, "candidate_global_summary.tsv")
    cell_df.to_csv(cells_path, sep="\t", index=False)
    region_df.to_csv(regions_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    outputs = {"cells": cells_path, "regions": regions_path, "summary": summary_path}
    if params.write_heatmap:
        outputs.update(_write_residual_heatmaps(params, residual, chrom, cell_df))
    return outputs


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="pyInferCNV residual candidate clone caller.")
    p.add_argument("--h5ad", required=True, type=Path)
    p.add_argument("--fastcnv-prefix", required=True, type=Path)
    p.add_argument("--output-prefix", required=True, type=Path)
    p.add_argument("--state-key", required=True)
    p.add_argument("--layer", default="auto")
    p.add_argument("--min-state-cells", type=int, default=60)
    p.add_argument("--min-component-cells", type=int, default=20)
    p.add_argument("--min-component-fraction", type=float, default=0.03)
    p.add_argument("--min-abs-region-score", type=float, default=0.05)
    p.add_argument("--min-separation-mad", type=float, default=1.5)
    p.add_argument("--exclude-chromosomes", default="chrX,chrY,chrM")
    p.add_argument("--no-keep-loy", dest="keep_loy", action="store_false")
    p.add_argument("--pyinfercnv-path", type=Path, default=Path("/Users/saljh8/Documents/GitHub/pyInferCNV"))
    p.add_argument("--no-heatmap", dest="write_heatmap", action="store_false")
    p.add_argument("--heatmap-max-cells", type=int, default=20000)
    p.add_argument("--heatmap-bin-genes", type=int, default=50)
    p.add_argument("--heatmap-min-chrom-bins", type=int, default=25)
    p.add_argument("--heatmap-filter-threshold", type=float, default=0.03)
    return p


def _split(value: str) -> tuple[str, ...]:
    return tuple(v.strip() for v in str(value).split(",") if v.strip())


def main(argv: Optional[Sequence[str]] = None) -> Dict[str, Path]:
    args = build_parser().parse_args(argv)
    params = ResidualCloneParams(
        h5ad=args.h5ad,
        fastcnv_prefix=args.fastcnv_prefix,
        output_prefix=args.output_prefix,
        state_key=args.state_key,
        layer=args.layer,
        min_state_cells=args.min_state_cells,
        min_component_cells=args.min_component_cells,
        min_component_fraction=args.min_component_fraction,
        min_abs_region_score=args.min_abs_region_score,
        min_separation_mad=args.min_separation_mad,
        exclude_chromosomes=_split(args.exclude_chromosomes),
        keep_loy=bool(args.keep_loy),
        pyinfercnv_path=args.pyinfercnv_path,
        write_heatmap=bool(args.write_heatmap),
        heatmap_max_cells=int(args.heatmap_max_cells),
        heatmap_bin_genes=int(args.heatmap_bin_genes),
        heatmap_min_chrom_bins=int(args.heatmap_min_chrom_bins),
        heatmap_filter_threshold=float(args.heatmap_filter_threshold),
    )
    outputs = run_candidate(params)
    for key, path in outputs.items():
        print(f"{key}: {path}")
    return outputs


if __name__ == "__main__":
    main()
