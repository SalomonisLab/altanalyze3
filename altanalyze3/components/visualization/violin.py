#!/usr/bin/env python3
"""
Lineage-ordered violin plot exports for cellHarmony-aligned AnnData outputs.

Modes
-----
1. across-lineages
   Plot expression of one or more genes across all lineage / population groups
   ordered by ``adata.uns["lineage_order"]`` when available.

2. population-obs
   Plot expression of one or more genes within a single selected population,
   split across values from another obs field.
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path
from typing import List, Optional, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from approximate_umap import (
        _build_preview_palette,
        _clean_string,
        _extract_lineage_order,
        _flatten_expr,
        _load_adata,
        _sanitize_filename_component,
        _subset_adata_for_plot,
    )
except ImportError:  # pragma: no cover
    from altanalyze3.components.visualization.approximate_umap import (
        _build_preview_palette,
        _clean_string,
        _extract_lineage_order,
        _flatten_expr,
        _load_adata,
        _sanitize_filename_component,
        _subset_adata_for_plot,
    )


def _log_timing(label: str, started_at: float) -> float:
    elapsed = max(0.0, time.perf_counter() - started_at)
    print(f"[timing] {label}={elapsed:.2f}s")
    return elapsed


def _configure_matplotlib_pdf_style() -> None:
    plt.rcParams["axes.linewidth"] = 0.5
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
    plt.rcParams["figure.facecolor"] = "white"


def _resolve_group_order(series: pd.Series, lineage_order: Optional[Sequence[str]] = None) -> List[str]:
    values = series.astype(str)
    unique_values = [str(v) for v in pd.unique(values)]
    if lineage_order:
        lineage = [str(v) for v in lineage_order]
        lineage_set = set(lineage)
        ordered = [value for value in lineage if value in unique_values]
        remainder = [value for value in unique_values if value not in lineage_set]
        return ordered + remainder
    if isinstance(series.dtype, pd.CategoricalDtype):
        categories = [str(v) for v in series.cat.categories.tolist()]
        ordered = [value for value in categories if value in unique_values]
        remainder = [value for value in unique_values if value not in set(ordered)]
        return ordered + remainder
    return unique_values


def _parse_list_arg(values: Optional[Sequence[str]]) -> List[str]:
    resolved: List[str] = []
    for entry in values or []:
        text = _clean_string(entry)
        if text:
            resolved.append(text)
    return resolved


def _plot_gene_groups(ax: matplotlib.axes.Axes, group_payload: List[dict], *, title: str) -> None:
    positions = np.arange(1, len(group_payload) + 1)
    labels = [entry["label"] for entry in group_payload]
    palette = _build_preview_palette(labels)
    parts = ax.violinplot(
        [entry["values"] for entry in group_payload],
        positions=positions,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )
    for body, entry in zip(parts["bodies"], group_payload):
        color = palette.get(entry["label"], "#64748b")
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.55)
        body.set_linewidth(0.7)
    if "cmedians" in parts:
        parts["cmedians"].set_color("#0f172a")
        parts["cmedians"].set_linewidth(0.8)
    for idx, entry in enumerate(group_payload, start=1):
        vals = np.asarray(entry["values"], dtype=float)
        if vals.size == 0:
            continue
        color = palette.get(entry["label"], "#64748b")
        jitter = np.random.default_rng(0).normal(0, 0.035, size=len(vals))
        ax.scatter(
            np.full(len(vals), idx) + jitter,
            vals,
            s=2,
            c=color,
            alpha=0.25,
            linewidths=0,
            rasterized=False,
        )
    ax.set_xticks(positions)
    ax.set_xticklabels([entry["label"] for entry in group_payload], rotation=45, ha="right")
    ax.set_title(title)
    ax.set_ylabel("Expression")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _build_across_lineages_groups(
    adata,
    *,
    gene: str,
    cluster_key: str,
    lineage_order: Optional[Sequence[str]],
) -> List[dict]:
    if cluster_key not in adata.obs:
        raise KeyError(f"obs column '{cluster_key}' was not found in the AnnData object.")
    if gene not in adata.var_names:
        raise KeyError(f"Gene '{gene}' was not found in the AnnData object.")
    clusters = adata.obs[cluster_key]
    ordered_groups = _resolve_group_order(clusters, lineage_order)
    values = _flatten_expr(adata[:, gene].X).astype(float)
    payload: List[dict] = []
    cluster_values = clusters.astype(str).to_numpy()
    for group in ordered_groups:
        mask = cluster_values == str(group)
        finite = values[mask]
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            continue
        payload.append(
            {
                "label": str(group),
                "values": finite.tolist(),
                "mean": float(np.mean(finite)),
            }
        )
    return payload


def _build_population_obs_groups(
    adata,
    *,
    gene: str,
    cluster_key: str,
    population: str,
    compare_obs_field: str,
    compare_values: Optional[Sequence[str]],
) -> List[dict]:
    if cluster_key not in adata.obs:
        raise KeyError(f"obs column '{cluster_key}' was not found in the AnnData object.")
    if compare_obs_field not in adata.obs:
        raise KeyError(f"obs column '{compare_obs_field}' was not found in the AnnData object.")
    if gene not in adata.var_names:
        raise KeyError(f"Gene '{gene}' was not found in the AnnData object.")

    population_mask = adata.obs[cluster_key].astype(str) == str(population)
    if not bool(population_mask.any()):
        raise ValueError(f"No cells matched population '{population}' in obs['{cluster_key}'].")
    subset = adata[population_mask].copy()
    compare_series = subset.obs[compare_obs_field]
    ordered_groups = _resolve_group_order(compare_series)
    if compare_values:
        selected = {str(value) for value in compare_values}
        ordered_groups = [group for group in ordered_groups if group in selected]
    values = _flatten_expr(subset[:, gene].X).astype(float)
    compare_values_array = compare_series.astype(str).to_numpy()
    payload: List[dict] = []
    for group in ordered_groups:
        mask = compare_values_array == str(group)
        finite = values[mask]
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            continue
        payload.append(
            {
                "label": str(group),
                "values": finite.tolist(),
                "mean": float(np.mean(finite)),
            }
        )
    if not payload:
        raise ValueError(
            f"No values from obs['{compare_obs_field}'] remained after filtering within population '{population}'."
        )
    return payload


def write_lineage_violin_pdf(
    h5ad,
    *,
    output_pdf: Path,
    genes: Sequence[str],
    cluster_key: str,
    mode: str = "across-lineages",
    population: Optional[str] = None,
    compare_obs_field: Optional[str] = None,
    compare_values: Optional[Sequence[str]] = None,
    restrict_obs_field: Optional[str] = None,
    restrict_obs_value: Optional[Sequence[str] | str] = None,
    restrict_obs_mode: str = "include",
) -> str:
    step_started_at = time.perf_counter()
    adata = _load_adata(h5ad)
    _log_timing("lineage_violin.load_h5ad", step_started_at)

    step_started_at = time.perf_counter()
    subset = _subset_adata_for_plot(
        adata,
        obs_field=restrict_obs_field,
        obs_value=restrict_obs_value,
        obs_mode=restrict_obs_mode,
        dataset_label="query",
    )
    _log_timing("lineage_violin.apply_restriction", step_started_at)

    genes = [str(gene) for gene in genes if _clean_string(gene)]
    if not genes:
        raise ValueError("At least one target gene must be supplied.")

    step_started_at = time.perf_counter()
    lineage_order = _extract_lineage_order(subset, cluster_key)
    if lineage_order:
        preview = ", ".join(lineage_order[:10]) + (", …" if len(lineage_order) > 10 else "")
        print(f"[info] Using lineage_order from AnnData uns (n={len(lineage_order)}): {preview}")
    else:
        print("[warn] lineage_order not found in AnnData uns; using detected group order.")
    _log_timing("lineage_violin.resolve_order", step_started_at)

    step_started_at = time.perf_counter()
    plot_payloads: List[dict] = []
    for gene in genes:
        if mode == "across-lineages":
            groups = _build_across_lineages_groups(
                subset,
                gene=gene,
                cluster_key=cluster_key,
                lineage_order=lineage_order,
            )
            title = f"{gene} across {cluster_key}"
        elif mode == "population-obs":
            if not _clean_string(population):
                raise ValueError("--population is required when mode=population-obs.")
            if not _clean_string(compare_obs_field):
                raise ValueError("--compare-obs-field is required when mode=population-obs.")
            groups = _build_population_obs_groups(
                subset,
                gene=gene,
                cluster_key=cluster_key,
                population=str(population),
                compare_obs_field=str(compare_obs_field),
                compare_values=compare_values,
            )
            title = f"{gene} in {population} across {compare_obs_field}"
        else:
            raise ValueError(f"Unsupported mode '{mode}'.")
        if not groups:
            raise ValueError(f"No finite expression values were available to plot for gene '{gene}'.")
        plot_payloads.append({"gene": gene, "groups": groups, "title": title})
    _log_timing("lineage_violin.build_groups", step_started_at)

    _configure_matplotlib_pdf_style()
    max_groups = max((len(entry["groups"]) for entry in plot_payloads), default=1)
    width = max(8.5, min(22.0, 2.0 + 0.35 * max_groups))
    height = max(4.5, 3.6 * len(plot_payloads))
    step_started_at = time.perf_counter()
    fig, axes = plt.subplots(len(plot_payloads), 1, figsize=(width, height), squeeze=False)
    try:
        for ax, payload in zip(axes.ravel(), plot_payloads):
            _plot_gene_groups(ax, payload["groups"], title=payload["title"])
        fig.tight_layout()
        output_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_pdf, format="pdf", bbox_inches="tight")
    finally:
        plt.close(fig)
    _log_timing("lineage_violin.render_pdf", step_started_at)
    return str(output_pdf)


def _parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Lineage-ordered violin PDF export from cellHarmony h5ad outputs.")
    parser.add_argument("--h5ad", required=True, help="Input cellHarmony output h5ad.")
    parser.add_argument("--cluster-key", required=True, help="obs column used for lineage / population labels.")
    parser.add_argument(
        "--mode",
        choices=("across-lineages", "population-obs"),
        default="across-lineages",
        help="Plot mode.",
    )
    parser.add_argument(
        "--gene",
        action="append",
        required=True,
        help="Target gene. Repeat for multiple genes.",
    )
    parser.add_argument("--population", help="Required for mode=population-obs.")
    parser.add_argument("--compare-obs-field", help="Required for mode=population-obs.")
    parser.add_argument(
        "--compare-obs-value",
        action="append",
        default=[],
        help="Optional obs values to retain in mode=population-obs. Repeat for multiple values.",
    )
    parser.add_argument("--restrict-obs-field", help="obs column used to restrict which cells are included.")
    parser.add_argument(
        "--restrict-obs-value",
        help="obs value(s) to include/exclude. Comma- or pipe-delimited values are accepted.",
    )
    parser.add_argument(
        "--restrict-obs-mode",
        choices=("include", "exclude"),
        default="include",
        help="Whether the supplied restrict obs values are included or excluded.",
    )
    parser.add_argument("--outdir", required=True, help="Directory where the PDF will be written.")
    parser.add_argument(
        "--output-pdf",
        help="Output PDF filename or path. Defaults to <mode>-<cluster-key>.pdf in --outdir.",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable INFO logging.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    main_started_at = time.perf_counter()
    args = _parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.output_pdf:
        output_candidate = Path(args.output_pdf)
        output_pdf = output_candidate if output_candidate.is_absolute() else outdir / output_candidate
    else:
        output_pdf = outdir / f"{args.mode}-{_sanitize_filename_component(args.cluster_key)}.pdf"

    compare_values = _parse_list_arg(args.compare_obs_value)
    write_lineage_violin_pdf(
        args.h5ad,
        output_pdf=output_pdf,
        genes=args.gene,
        cluster_key=args.cluster_key,
        mode=args.mode,
        population=args.population,
        compare_obs_field=args.compare_obs_field,
        compare_values=compare_values,
        restrict_obs_field=args.restrict_obs_field,
        restrict_obs_value=args.restrict_obs_value,
        restrict_obs_mode=args.restrict_obs_mode,
    )
    _log_timing("lineage_violin.total", main_started_at)
    print(f"[info] Wrote violin PDF to {output_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
