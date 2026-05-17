"""
Static plotting utilities for GO-Elite enrichment outputs.
"""

from __future__ import annotations

from pathlib import Path
import logging
from typing import Iterable, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import transforms as mtransforms
from matplotlib.ticker import LogFormatterMathtext
import numpy as np
import pandas as pd

GOELITE_HIGHLIGHT_KEYWORDS: tuple[str, ...] = (
    "cell cycle",
    "mitotic",
    "splicing",
    "mrna processing",
    "proliferation",
    "cytokine",
    "death",
    "chromatin",
    "lipid",
    "circadian",
    "tp53",
    "wnt",
    "tgf",
    "tnf",
    "granule",
)


def _is_finite_number(value: object) -> bool:
    try:
        return np.isfinite(float(value))
    except Exception:
        return False


def _prepare_goelite_plot_frame(frame: pd.DataFrame) -> pd.DataFrame:
    if frame is None or frame.empty:
        return pd.DataFrame()
    prepared = frame.copy()
    if "population" not in prepared.columns:
        prepared["population"] = "GO-Elite"
    if "term_name" not in prepared.columns:
        prepared["term_name"] = prepared.get("term_id", pd.Series("", index=prepared.index)).astype(str)
    prepared["population"] = prepared["population"].astype(str)
    prepared["term_name"] = prepared["term_name"].astype(str)
    prepared["selected"] = prepared.get("selected", False).fillna(False).astype(bool)
    prepared["fdr"] = pd.to_numeric(prepared.get("fdr"), errors="coerce")
    prepared["p_value"] = pd.to_numeric(prepared.get("p_value"), errors="coerce")
    prepared["z_score"] = pd.to_numeric(prepared.get("z_score"), errors="coerce")
    prepared["fdr_plot"] = prepared["fdr"].fillna(prepared["p_value"]).clip(lower=1e-300, upper=1.0)
    prepared["is_positive_sig"] = prepared["fdr_plot"].le(0.05) & prepared["z_score"].gt(2.0)
    prepared["is_selected_positive_sig"] = prepared["selected"] & prepared["is_positive_sig"]
    prepared = prepared.loc[
        prepared["z_score"].map(_is_finite_number)
        & prepared["fdr_plot"].map(_is_finite_number)
        & prepared["fdr_plot"].gt(0)
    ].copy()
    prepared.sort_values(["population", "fdr_plot", "p_value", "z_score", "term_name"], ascending=[True, True, True, False, True], inplace=True)
    return prepared


def _select_goelite_labels(
    frame: pd.DataFrame,
    keywords: Sequence[str] = GOELITE_HIGHLIGHT_KEYWORDS,
    *,
    top_label_count: int = 4,
) -> list[dict[str, object]]:
    if frame.empty:
        return []
    top_label_count = max(0, int(top_label_count))
    selected = frame.loc[frame["is_selected_positive_sig"]].copy()
    selected.sort_values(["fdr_plot", "p_value", "z_score", "term_name"], ascending=[True, True, False, True], inplace=True)
    labels: list[dict[str, object]] = []
    used_names: set[str] = set()
    for _, row in selected.head(top_label_count).iterrows():
        name_lc = str(row["term_name"]).strip().lower()
        used_names.add(name_lc)
        labels.append(
            {
                "term_name": str(row["term_name"]),
                "z_score": float(row["z_score"]),
                "fdr_plot": float(row["fdr_plot"]),
                "label_color": "#111827",
            }
        )
    keyword_term = None
    for _, row in selected.iloc[top_label_count:].iterrows():
        name_lc = str(row["term_name"]).strip().lower()
        if name_lc in used_names:
            continue
        if any(keyword in name_lc for keyword in keywords):
            keyword_term = row
            break
    if keyword_term is None:
        fallback = frame.loc[frame["z_score"].gt(0)].copy()
        fallback.sort_values(["fdr_plot", "p_value", "z_score", "term_name"], ascending=[True, True, False, True], inplace=True)
        for _, row in fallback.iterrows():
            name_lc = str(row["term_name"]).strip().lower()
            if name_lc in used_names:
                continue
            if any(keyword in name_lc for keyword in keywords):
                keyword_term = row
                break
    if keyword_term is not None:
        labels.append(
            {
                "term_name": str(keyword_term["term_name"]),
                "z_score": float(keyword_term["z_score"]),
                "fdr_plot": float(keyword_term["fdr_plot"]),
                "label_color": "#111827",
            }
        )
    return labels


def write_goelite_scatter_pdf(
    frame: pd.DataFrame,
    out_path: str | Path,
    *,
    title_prefix: str = "GO-Elite",
    keywords: Sequence[str] = GOELITE_HIGHLIGHT_KEYWORDS,
    top_label_count: int = 4,
) -> Path | None:
    prepared = _prepare_goelite_plot_frame(frame)
    if prepared.empty:
        return None

    output_path = Path(out_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
    plt.rcParams["axes.linewidth"] = 0.8
    logging.getLogger("fontTools").setLevel(logging.ERROR)
    logging.getLogger("fontTools.subset").setLevel(logging.ERROR)

    annotation_offsets = [(-6, 38), (-2, 10), (2, -10), (6, -34), (18, -60)]
    populations: Iterable[tuple[str, pd.DataFrame]] = prepared.groupby("population", sort=False)
    log_formatter = LogFormatterMathtext(base=10)

    with PdfPages(output_path) as pdf:
        for population, pop_df in populations:
            if pop_df.empty:
                continue
            labels = _select_goelite_labels(pop_df, keywords=keywords, top_label_count=top_label_count)
            positive = pop_df.loc[pop_df["is_selected_positive_sig"]]
            background = pop_df.loc[~pop_df["is_positive_sig"]]

            fig, ax = plt.subplots(figsize=(9.2, 7.8))
            if not background.empty:
                ax.scatter(
                    background["z_score"].to_numpy(dtype=float),
                    background["fdr_plot"].to_numpy(dtype=float),
                    s=42,
                    c="#d1d5db",
                    alpha=0.95,
                    linewidths=0,
                    zorder=2,
                )
            if not positive.empty:
                ax.scatter(
                    positive["z_score"].to_numpy(dtype=float),
                    positive["fdr_plot"].to_numpy(dtype=float),
                    s=48,
                    c="#1f19c7",
                    alpha=0.98,
                    linewidths=0,
                    zorder=3,
                )
            for index, label in enumerate(labels):
                dx, dy = annotation_offsets[min(index, len(annotation_offsets) - 1)]
                ax.annotate(
                    str(label["term_name"]),
                    xy=(float(label["z_score"]), float(label["fdr_plot"])),
                    xytext=(dx, dy),
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    fontsize=9,
                    color=str(label["label_color"]),
                    arrowprops={
                        "arrowstyle": "-",
                        "color": str(label["label_color"]),
                        "linewidth": 1.0,
                        "alpha": 0.9,
                        "shrinkA": 0,
                        "shrinkB": 0,
                    },
                    zorder=4,
                )

            x_values = pop_df["z_score"].to_numpy(dtype=float)
            y_values = pop_df["fdr_plot"].to_numpy(dtype=float)
            x_min = min(-10.0, float(np.floor(np.min(x_values) - 0.5)))
            x_max = max(20.0, float(np.ceil(np.max(x_values) + 2.5)))
            y_min = max(float(np.min(y_values)) * 0.5, 1e-300)
            ax.set_yscale("log")
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_min, 1.0)
            ax.set_title(f"{title_prefix}: {population}")
            ax.axvline(0.0, color="#111827", linewidth=1.2, alpha=0.95, zorder=1)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.grid(False)
            ax.tick_params(axis="y", which="both", left=False, right=False, labelleft=False)
            ax.set_xlabel("Z-Score", fontsize=12, color="#111827", labelpad=10)
            ax.set_ylabel("")
            fig.tight_layout()
            fig.canvas.draw()
            guide_transform = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
            guide_display_xy = guide_transform.transform((0.0, 0.5))
            guide_figure_xy = fig.transFigure.inverted().transform(guide_display_xy)
            y_ticks = [
                tick
                for tick in ax.get_yticks()
                if np.isfinite(tick) and tick > 0 and y_min <= tick <= 1.0
            ]
            for tick in y_ticks:
                tick_label = log_formatter(tick)
                ax.annotate(
                    "",
                    xy=(0.0, tick),
                    xycoords="data",
                    xytext=(-6, 0),
                    textcoords="offset points",
                    arrowprops={
                        "arrowstyle": "-",
                        "color": "#111827",
                        "linewidth": 1.0,
                        "alpha": 1.0,
                        "shrinkA": 0,
                        "shrinkB": 0,
                    },
                    zorder=5,
                )
                ax.annotate(
                    tick_label,
                    xy=(0.0, tick),
                    xycoords="data",
                    xytext=(-10, 0),
                    textcoords="offset points",
                    ha="right",
                    va="center",
                    fontsize=10,
                    color="#111827",
                    zorder=5,
                )
            ax.annotate(
                "Fishers FDR p",
                xy=(0.0, 0.5),
                xycoords=guide_transform,
                xytext=(-46, 0),
                textcoords="offset points",
                rotation=90,
                rotation_mode="anchor",
                va="center",
                ha="center",
                fontsize=12,
                color="#111827",
                zorder=5,
            )
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    return output_path
