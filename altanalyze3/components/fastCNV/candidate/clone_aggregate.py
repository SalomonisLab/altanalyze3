"""Aggregate clone-level candidate CNV caller for fastCNV score matrices.

This candidate is intentionally separate from production fastCNV. It consumes a
standard fastCNV prefix (`*.cnv_window_scores.npz`, `*.cnv_windows.tsv`,
`*.cnv_cells.tsv`) and asks a different question from the per-cell interval
caller: is there a state-aware subpopulation with a coherent regional shift?

The method is unsupervised and region-generic:

* aggregate existing fastCNV window scores to chromosome arms;
* within each cell state, split each region into two score components;
* accept an altered component only when it is sufficiently separated from the
  state's central component and has enough cells;
* emit candidate cell and region tables in the same candidate naming convention.

Existing chrY/LOY calls from fastCNV are carried through unchanged so the LOY
positive control remains governed by the validated unbiased fastCNV chrY path.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence

import numpy as np
import pandas as pd


CENTROMERE_HG38 = {
    "chr1": 123_000_000, "chr2": 93_900_000, "chr3": 90_900_000,
    "chr4": 50_000_000, "chr5": 48_800_000, "chr6": 59_800_000,
    "chr7": 60_100_000, "chr8": 45_200_000, "chr9": 43_000_000,
    "chr10": 39_800_000, "chr11": 53_400_000, "chr12": 35_500_000,
    "chr13": 17_700_000, "chr14": 17_200_000, "chr15": 19_000_000,
    "chr16": 36_800_000, "chr17": 25_100_000, "chr18": 18_500_000,
    "chr19": 26_200_000, "chr20": 28_100_000, "chr21": 12_000_000,
    "chr22": 15_000_000,
}


@dataclass(frozen=True)
class AggregateParams:
    fastcnv_prefix: Path
    output_prefix: Path
    h5ad: Optional[Path] = None
    gene_coordinates: Optional[Path] = None
    layer: str = "auto"
    query_centered: bool = False
    window_genes: int = 101
    stride_genes: int = 10
    min_state_cells: int = 60
    min_component_cells: int = 40
    min_component_fraction: float = 0.05
    min_abs_region_score: float = 0.8
    min_separation_mad: float = 2.0
    max_neutral_abs_median: float = 0.35
    include_chromosomes: Optional[tuple[str, ...]] = None
    exclude_chromosomes: tuple[str, ...] = ("chrY",)
    keep_loy: bool = True


def _prefix_path(prefix: Path, suffix: str) -> Path:
    return Path(f"{prefix}.{suffix}")


def _load_fastcnv(prefix: Path) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame, np.ndarray]:
    scores_path = _prefix_path(prefix, "cnv_window_scores.npz")
    windows_path = _prefix_path(prefix, "cnv_windows.tsv")
    cells_path = _prefix_path(prefix, "cnv_cells.tsv")
    if not scores_path.exists():
        raise FileNotFoundError(scores_path)
    if not windows_path.exists():
        raise FileNotFoundError(windows_path)
    if not cells_path.exists():
        raise FileNotFoundError(cells_path)
    npz = np.load(scores_path, allow_pickle=True)
    scores = np.asarray(npz["scores"], dtype=np.float32)
    barcodes = np.asarray(npz["cell_barcodes"]).astype(str)
    windows = pd.read_csv(windows_path, sep="\t")
    cells = pd.read_csv(cells_path, sep="\t")
    cells["CellBarcode"] = cells["CellBarcode"].astype(str)
    cells = cells.set_index("CellBarcode").reindex(barcodes).reset_index()
    cells = cells.rename(columns={"index": "CellBarcode"})
    if "cell_state" not in cells.columns:
        states = np.asarray(npz.get("states", np.array(["unknown"] * len(barcodes)))).astype(str)
        cells["cell_state"] = states
    if "sample" not in cells.columns:
        cells["sample"] = ""
    return scores, windows, cells, barcodes


def _matrix(adata, layer: str):
    if layer == "auto":
        return adata.layers["counts"] if "counts" in adata.layers else adata.X
    if layer == "X":
        return adata.X
    return adata.layers[layer]


def _load_query_centered(params: AggregateParams) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame, np.ndarray]:
    if params.h5ad is None or params.gene_coordinates is None:
        raise ValueError("--query-centered requires --h5ad and --gene-coordinates")
    import anndata as ad
    import scipy.sparse as sp
    from altanalyze3.components.fastCNV.main import load_gene_coordinates, build_windows

    adata = ad.read_h5ad(params.h5ad)
    coords = load_gene_coordinates(params.gene_coordinates, adata.var_names).reset_index(drop=True)
    coords = coords[coords["chr"].astype(str).str.match(r"^chr([0-9]+|X|Y)$")].reset_index(drop=True)
    matrix = _matrix(adata, params.layer)
    cols = coords["var_index"].to_numpy(dtype=np.int64)
    sub = matrix[:, cols]
    if sp.issparse(sub):
        sub = sub.toarray()
    sub = np.asarray(sub, dtype=np.float32)
    if params.layer == "X" and np.nanmax(sub) < 30:
        expr = sub
    else:
        lib = np.maximum(np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32), 1.0)
        expr = np.log1p(sub * (10000.0 / lib[:, None])).astype(np.float32)
    gene_center = np.nanmedian(expr, axis=0).astype(np.float32)
    residual = expr - gene_center[None, :]
    cell_center = np.nanmedian(residual[:, coords["chr"].astype(str).str.match(r"^chr[0-9]+$").to_numpy()], axis=1)
    residual = residual - cell_center[:, None]
    wins = build_windows(coords, params.window_genes, params.stride_genes, 25)
    score = np.zeros((adata.n_obs, len(wins)), dtype=np.float32)
    rows = []
    coord_chr = coords["chr"].astype(str).to_numpy()
    coord_start = coords["start"].to_numpy()
    coord_end = coords["end"].to_numpy()
    coord_gene = np.asarray([str(adata.var_names[int(i)]) for i in coords["var_index"].to_numpy()])
    for i, win in enumerate(wins):
        vals = residual[:, win.start_offset:win.end_offset]
        score[:, i] = np.nanmean(vals, axis=1)
        rows.append(
            {
                "window_id": win.window_id,
                "chr": win.chrom,
                "start": int(coord_start[win.start_offset]),
                "end": int(coord_end[win.end_offset - 1]),
                "start_gene": str(coord_gene[win.start_offset]),
                "end_gene": str(coord_gene[win.end_offset - 1]),
            }
        )
    cells = pd.DataFrame(
        {
            "CellBarcode": adata.obs_names.astype(str),
            "cell_state": adata.obs.iloc[:, 0].astype(str).to_numpy(),
            "sample": "",
            "cnv_status": "WT",
            "global_clone_id": "WT",
        }
    )
    if params.fastcnv_prefix:
        fast_cells = pd.read_csv(_prefix_path(params.fastcnv_prefix, "cnv_cells.tsv"), sep="\t")
        fast_cells["CellBarcode"] = fast_cells["CellBarcode"].astype(str)
        for col in ("cell_state", "sample", "cnv_status", "global_clone_id"):
            if col in fast_cells.columns:
                cells[col] = fast_cells.set_index("CellBarcode").reindex(cells["CellBarcode"])[col].fillna(cells[col]).to_numpy()
    return score, pd.DataFrame(rows), cells, adata.obs_names.astype(str).to_numpy()


def _region_table(windows: pd.DataFrame, include: Optional[Iterable[str]], exclude: Iterable[str]) -> pd.DataFrame:
    include_set = set(include or [])
    exclude_set = set(exclude or [])
    rows = []
    for chrom, group in windows.reset_index().groupby("chr", sort=False):
        chrom = str(chrom)
        if include_set and chrom not in include_set:
            continue
        if chrom in exclude_set:
            continue
        if chrom in CENTROMERE_HG38:
            cen = CENTROMERE_HG38[chrom]
            parts = [(f"{chrom}p", group[group["end"] <= cen]), (f"{chrom}q", group[group["start"] >= cen])]
        else:
            parts = [(chrom, group)]
        for name, sub in parts:
            if sub.empty:
                continue
            rows.append(
                {
                    "region": name,
                    "chr": chrom,
                    "start": int(sub["start"].min()),
                    "end": int(sub["end"].max()),
                    "start_gene": str(sub["start_gene"].iloc[0]) if "start_gene" in sub.columns else "",
                    "end_gene": str(sub["end_gene"].iloc[-1]) if "end_gene" in sub.columns else "",
                    "n_windows": int(sub.shape[0]),
                    "window_indices": sub["index"].to_numpy(dtype=np.int64),
                }
            )
    return pd.DataFrame(rows)


def _two_component_split(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    finite = np.isfinite(values)
    labels = np.zeros(values.shape[0], dtype=np.int8)
    if finite.sum() < 4:
        return labels, np.array([np.nan, np.nan])
    x = values[finite]
    c = np.array([np.quantile(x, 0.25), np.quantile(x, 0.75)], dtype=float)
    if abs(c[1] - c[0]) < 1e-6:
        return labels, c
    for _ in range(30):
        d0 = np.abs(x - c[0])
        d1 = np.abs(x - c[1])
        lab = d1 < d0
        new = np.array([
            np.mean(x[~lab]) if np.any(~lab) else c[0],
            np.mean(x[lab]) if np.any(lab) else c[1],
        ])
        if np.max(np.abs(new - c)) < 1e-4:
            c = new
            break
        c = new
    if c[0] > c[1]:
        c = c[::-1]
        lab = ~lab
    labels[np.flatnonzero(finite)] = lab.astype(np.int8)
    return labels, c


def _mad(values: np.ndarray) -> float:
    vals = values[np.isfinite(values)]
    if vals.size == 0:
        return 1.0
    med = float(np.median(vals))
    scale = 1.4826 * float(np.median(np.abs(vals - med)))
    return max(scale, 1e-3)


def run_aggregate(params: AggregateParams) -> Dict[str, Path]:
    if params.query_centered:
        scores, windows, cells, barcodes = _load_query_centered(params)
    else:
        scores, windows, cells, barcodes = _load_fastcnv(params.fastcnv_prefix)
    regions = _region_table(windows, params.include_chromosomes, params.exclude_chromosomes)
    if regions.empty:
        raise ValueError("No candidate regions after chromosome include/exclude filters.")

    calls: list[dict] = []
    state_values = cells["cell_state"].astype(str).to_numpy()
    sample_values = cells["sample"].fillna("").astype(str).to_numpy()
    for state in pd.unique(state_values):
        state_rows = np.flatnonzero(state_values == state)
        if state_rows.size < params.min_state_cells:
            continue
        for reg in regions.itertuples(index=False):
            vals = np.nanmean(scores[state_rows[:, None], reg.window_indices], axis=1)
            labels, centers = _two_component_split(vals)
            if not np.all(np.isfinite(centers)):
                continue
            state_median = float(np.nanmedian(vals))
            neutral_idx = int(np.argmin(np.abs(centers - state_median)))
            altered_idx = 1 - neutral_idx
            neutral_center = float(centers[neutral_idx])
            altered_center = float(centers[altered_idx])
            if abs(neutral_center) > params.max_neutral_abs_median:
                neutral_center = state_median
            altered_local = np.flatnonzero(labels == altered_idx)
            if altered_local.size < params.min_component_cells:
                continue
            frac = altered_local.size / state_rows.size
            if frac < params.min_component_fraction:
                continue
            delta = altered_center - neutral_center
            sep = abs(delta) / _mad(vals)
            if abs(altered_center) < params.min_abs_region_score or sep < params.min_separation_mad:
                continue
            call = "gain" if altered_center > neutral_center else "loss"
            global_rows = state_rows[altered_local]
            for row in global_rows:
                calls.append(
                    {
                        "CellBarcode": barcodes[row],
                        "cell_state": state,
                        "sample": sample_values[row],
                        "call": call,
                        "chr": reg.chr,
                        "start": int(reg.start),
                        "end": int(reg.end),
                        "start_gene": reg.start_gene,
                        "end_gene": reg.end_gene,
                        "region": reg.region,
                        "n_windows": int(reg.n_windows),
                        "n_genes": int(reg.n_windows),
                        "mean_score": float(altered_center),
                        "max_score": float(abs(altered_center)),
                        "confidence": float(sep),
                        "mean_log2": float(altered_center),
                        "boundary_source": "aggregate_arm",
                    }
                )

    call_df = pd.DataFrame(calls)
    cell_out = cells[["CellBarcode", "cell_state", "sample"]].copy()
    cell_out["cnv_status"] = "WT"
    cell_out["n_cnv_intervals"] = 0
    cell_out["state_clone_id"] = "WT"
    cell_out["global_clone_id"] = "WT"
    if not call_df.empty:
        n_by_cell = call_df.groupby("CellBarcode").size()
        cell_out = cell_out.set_index("CellBarcode")
        cell_out.loc[n_by_cell.index, "cnv_status"] = "CNV"
        cell_out.loc[n_by_cell.index, "n_cnv_intervals"] = n_by_cell.astype(int)
        sig = (
            call_df.assign(sig=call_df["chr"].astype(str) + ":" + call_df["call"].astype(str))
            .groupby("CellBarcode")["sig"]
            .agg(lambda s: "|".join(sorted(set(s))))
        )
        clone_names = {sigv: f"clone{i + 1}" for i, sigv in enumerate(sorted(sig.unique()))}
        cell_out.loc[sig.index, "state_clone_id"] = sig.map(clone_names)
        cell_out.loc[sig.index, "global_clone_id"] = sig.map(clone_names)
        cell_out = cell_out.reset_index()

    if params.keep_loy and "global_clone_id" in cells.columns:
        loy_mask = cells["global_clone_id"].astype(str).eq("clone_loy")
        if loy_mask.any():
            cell_out.loc[loy_mask.to_numpy(), "cnv_status"] = "CNV"
            cell_out.loc[loy_mask.to_numpy(), "global_clone_id"] = "clone_loy"
            cell_out.loc[loy_mask.to_numpy(), "state_clone_id"] = "clone_loy"
            cell_out.loc[loy_mask.to_numpy(), "n_cnv_intervals"] = np.maximum(
                cell_out.loc[loy_mask.to_numpy(), "n_cnv_intervals"].astype(int), 1
            )
            loy_regions = pd.DataFrame(
                {
                    "CellBarcode": cell_out.loc[loy_mask.to_numpy(), "CellBarcode"],
                    "cell_state": cell_out.loc[loy_mask.to_numpy(), "cell_state"],
                    "sample": cell_out.loc[loy_mask.to_numpy(), "sample"],
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
            call_df = pd.concat([call_df, loy_regions], ignore_index=True)

    params.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    cells_path = _prefix_path(params.output_prefix, "candidate_cnv_cells.tsv")
    regions_path = _prefix_path(params.output_prefix, "candidate_cnv_regions.tsv")
    summary_path = _prefix_path(params.output_prefix, "candidate_global_summary.tsv")
    cell_out.to_csv(cells_path, sep="\t", index=False)
    if call_df.empty:
        call_df = pd.DataFrame(
            columns=[
                "CellBarcode", "cell_state", "sample", "call", "chr", "start", "end",
                "start_gene", "end_gene", "region", "n_windows", "n_genes",
                "mean_score", "max_score", "confidence", "mean_log2", "boundary_source",
            ]
        )
    call_df.to_csv(regions_path, sep="\t", index=False)
    summary = (
        call_df.groupby(["chr", "region", "call"], dropna=False)
        .agg(n_cells=("CellBarcode", "nunique"), n_states=("cell_state", "nunique"), mean_score=("mean_score", "mean"))
        .reset_index()
        if not call_df.empty else pd.DataFrame(columns=["chr", "region", "call", "n_cells", "n_states", "mean_score"])
    )
    summary.to_csv(summary_path, sep="\t", index=False)
    return {"cells": cells_path, "regions": regions_path, "summary": summary_path}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Aggregate clone-level candidate caller from a fastCNV prefix.")
    p.add_argument("--fastcnv-prefix", required=True, type=Path)
    p.add_argument("--output-prefix", required=True, type=Path)
    p.add_argument("--query-centered", action="store_true", help="Compute query-centered inferCNV-style scores from --h5ad.")
    p.add_argument("--h5ad", type=Path)
    p.add_argument("--gene-coordinates", type=Path)
    p.add_argument("--layer", default="auto")
    p.add_argument("--window-genes", type=int, default=101)
    p.add_argument("--stride-genes", type=int, default=10)
    p.add_argument("--min-state-cells", type=int, default=60)
    p.add_argument("--min-component-cells", type=int, default=40)
    p.add_argument("--min-component-fraction", type=float, default=0.05)
    p.add_argument("--min-abs-region-score", type=float, default=0.8)
    p.add_argument("--min-separation-mad", type=float, default=2.0)
    p.add_argument("--max-neutral-abs-median", type=float, default=0.35)
    p.add_argument("--include-chromosomes", default="")
    p.add_argument("--exclude-chromosomes", default="chrY")
    p.add_argument("--no-keep-loy", dest="keep_loy", action="store_false")
    return p


def _split_csv(value: str) -> Optional[tuple[str, ...]]:
    vals = tuple(v.strip() for v in str(value).split(",") if v.strip())
    return vals or None


def main(argv: Optional[Sequence[str]] = None) -> Dict[str, Path]:
    args = build_parser().parse_args(argv)
    params = AggregateParams(
        fastcnv_prefix=args.fastcnv_prefix,
        output_prefix=args.output_prefix,
        h5ad=args.h5ad,
        gene_coordinates=args.gene_coordinates,
        layer=args.layer,
        query_centered=bool(args.query_centered),
        window_genes=int(args.window_genes),
        stride_genes=int(args.stride_genes),
        min_state_cells=args.min_state_cells,
        min_component_cells=args.min_component_cells,
        min_component_fraction=args.min_component_fraction,
        min_abs_region_score=args.min_abs_region_score,
        min_separation_mad=args.min_separation_mad,
        max_neutral_abs_median=args.max_neutral_abs_median,
        include_chromosomes=_split_csv(args.include_chromosomes),
        exclude_chromosomes=_split_csv(args.exclude_chromosomes) or (),
        keep_loy=bool(args.keep_loy),
    )
    outputs = run_aggregate(params)
    for key, path in outputs.items():
        print(f"{key}: {path}")
    return outputs


if __name__ == "__main__":
    main()
