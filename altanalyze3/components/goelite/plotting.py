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
from matplotlib.collections import EllipseCollection
from matplotlib.colors import to_rgba
from matplotlib import transforms as mtransforms
from matplotlib.ticker import FuncFormatter, LogLocator, NullLocator
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


def _prepare_goelite_plot_frame(frame: pd.DataFrame, significance: str = "fdr", *, max_sig: float = 0.1,
                                min_z: float = 1.96) -> pd.DataFrame:
    """The plot's significance tiers come from the run's own gate (max_sig on the plotted metric, min_z),
    never from a plotting constant: an internal FDR<=0.05/z>2 cut drew terms the run had selected at FDR 0.1
    grey and unlabelled (Nathan, 2026-10-24: "you continually fail to validate the missing more critical
    required associations in the GO-Elite plots")."""
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
    # y-axis metric: raw p-value when significance=='rawp', else FDR
    _metric = prepared["p_value"] if significance == "rawp" else prepared["fdr"]
    prepared["fdr_plot"] = _metric.fillna(prepared["p_value"]).clip(lower=1e-300, upper=1.0)
    prepared["is_positive_sig"] = prepared["fdr_plot"].le(float(max_sig)) & prepared["z_score"].ge(float(min_z))
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
    top_label_count: int | None = None,
) -> list[dict[str, object]]:
    """Every term that passes the run's gate is labelled: GO-Elite-selected terms in black and significant
    terms pruned by the hierarchy step (a child or parent of a selected term) in blue-grey, so a reader sees
    "transcription elongation by RNA polymerase II" even when its parent was elected. top_label_count caps the
    labels to the most significant ones (default: all); the keyword term is appended when it passes the gate
    and is not already labelled. A label is never put on a term that fails the gate."""
    if frame.empty:
        return []
    selected = frame.loc[frame["is_positive_sig"]].copy()
    top_label_count = GOELITE_MAX_LABELS if top_label_count is None else max(0, int(top_label_count))
    # The cap keeps the most significant gated terms whatever their tier: a pruned child such as
    # "transcription elongation by RNA polymerase II" (p 2e-4) outranks a selected term at p 0.04.
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
                "label_color": "#111827" if bool(row["selected"]) else "#5B7BB5",
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
    if keyword_term is not None:
        labels.append(
            {
                "term_name": str(keyword_term["term_name"]),
                "z_score": float(keyword_term["z_score"]),
                "fdr_plot": float(keyword_term["fdr_plot"]),
                "label_color": "#111827" if bool(keyword_term["selected"]) else "#5B7BB5",
            }
        )
    return labels


GOELITE_MAX_LABELS = 40                  # labels a reader can still read on one page; the table holds every term
GOELITE_DOT_ELITE = "#1F19C7"            # selected by GO-Elite and positively significant
GOELITE_DOT_SIGNIFICANT = "#A8C4EE"      # significant (p<=0.05, z>2) but not elected
GOELITE_DOT_NOT_SIGNIFICANT = "#D1D5DB"  # everything else


def _add_dots(ax, subset, hex_color: str, diameter_pt: float, *, zorder: int) -> None:
    """Draw dots as inline vector circles with no stroke and no transparency group.

    ax.scatter writes a PathCollection whose per-marker alpha becomes a PDF transparency
    group, and Illustrator shows each group as an empty clipping box around the dot with no
    fill and no stroke. Fully opaque colours remove the group, and EllipseCollection in point
    units keeps every dot a real circle rather than a rasterised stamp.
    """
    if subset is None or subset.empty:
        return
    offsets = np.column_stack(
        (subset["z_score"].to_numpy(dtype=float), subset["fdr_plot"].to_numpy(dtype=float))
    )
    n = offsets.shape[0]
    collection = EllipseCollection(
        widths=np.full(n, float(diameter_pt)),
        heights=np.full(n, float(diameter_pt)),
        angles=np.zeros(n),
        units="points",
        offsets=offsets,
        offset_transform=ax.transData,
        facecolors=to_rgba(hex_color, 1.0),   # opaque: no PDF transparency group
        edgecolors="none",
        linewidths=0.0,
        zorder=zorder,
    )
    collection.set_alpha(None)
    ax.add_collection(collection)


def write_goelite_scatter_pdf(
    frame: pd.DataFrame,
    out_path: str | Path,
    *,
    title_prefix: str = "GO-Elite",
    keywords: Sequence[str] = GOELITE_HIGHLIGHT_KEYWORDS,
    top_label_count: int | None = None,
    significance: str = "fdr",
    max_sig: float = 0.1,
    min_z: float = 1.96,
) -> Path | None:
    prepared = _prepare_goelite_plot_frame(frame, significance=significance, max_sig=max_sig, min_z=min_z)
    if prepared.empty:
        return None
    y_axis_label = "Fishers raw p" if significance == "rawp" else "Fishers FDR p"
    gate_text = f"{'raw p' if significance == 'rawp' else 'FDR'} <= {max_sig:g}, z >= {min_z:g}"

    output_path = Path(out_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    plt.rcParams["axes.linewidth"] = 0.8
    logging.getLogger("fontTools").setLevel(logging.ERROR)
    logging.getLogger("fontTools.subset").setLevel(logging.ERROR)

    populations: Iterable[tuple[str, pd.DataFrame]] = prepared.groupby("population", sort=False)

    def _plain_log_tick(value: float, _pos=None) -> str:
        # Plain live text ("0.01", "1e-5"), never mathtext: LogFormatterMathtext wrote "10^-1" as a base and
        # a superscript in separate glyph runs, which Illustrator shows as separate objects (Nathan, 2026-10-23).
        if value <= 0 or not np.isfinite(value):
            return ""
        exponent = int(round(np.log10(value)))
        return f"{value:g}" if exponent >= -2 else f"1e{exponent}"


    with PdfPages(output_path) as pdf:
        for population, pop_df in populations:
            if pop_df.empty:
                continue
            labels = _select_goelite_labels(pop_df, keywords=keywords, top_label_count=top_label_count)
            positive = pop_df.loc[pop_df["is_selected_positive_sig"]]
            # every non-selected term is drawn (grey); previously terms that were
            # significant-but-not-selected fell in neither set and vanished.
            background = pop_df.loc[~pop_df["is_selected_positive_sig"]]

            fig, ax = plt.subplots(figsize=(9.2, 7.8))
            # Three tiers, not two. A term that clears p<=0.05 and z>2 but that GO-Elite did not
            # elect is significant, so it must not share the grey of the non-significant cloud.
            elite = pop_df.loc[pop_df["is_selected_positive_sig"]]
            sig_not_elite = pop_df.loc[pop_df["is_positive_sig"] & ~pop_df["is_selected_positive_sig"]]
            not_sig = pop_df.loc[~pop_df["is_positive_sig"]]
            _add_dots(ax, not_sig, GOELITE_DOT_NOT_SIGNIFICANT, 6.5, zorder=2)
            _add_dots(ax, sig_not_elite, GOELITE_DOT_SIGNIFICANT, 7.0, zorder=3)
            _add_dots(ax, elite, GOELITE_DOT_ELITE, 7.4, zorder=4)
            # Term labels are drawn AFTER the axis limits/scale are final (see
            # below) so each dot->term leader can be a single straight 2-point
            # Line2D rather than a FancyArrowPatch (which exports to PDF as a
            # multi-anchor path that is awkward to edit in Illustrator).

            x_values = pop_df["z_score"].to_numpy(dtype=float)
            y_values = pop_df["fdr_plot"].to_numpy(dtype=float)
            x_min = min(-10.0, float(np.floor(np.min(x_values) - 0.5)))
            x_max = max(20.0, float(np.ceil(np.max(x_values) + 2.5)))
            y_min = max(float(np.min(y_values)) * 0.5, 1e-300)
            ax.set_yscale("log")
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_min, 1.0)
            n_gated = int(pop_df["is_positive_sig"].sum())
            label_note = (f"all {len(labels)} gated terms labelled" if len(labels) >= n_gated
                          else f"the {len(labels)} most significant of {n_gated} gated terms labelled")
            ax.set_title(f"{title_prefix}: {population}\n{len(elite)} selected of {n_gated} terms passing the gate ({gate_text}) among {len(pop_df)} tested; "
                         f"{label_note}; black = selected, blue-grey = pruned by the hierarchy step", fontsize=8)
            ax.axvline(0.0, color="#111827", linewidth=1.2, zorder=1)  # opaque: no ExtGState
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            # One real y axis on the zero line with plain-text decade labels (one formatted axis object, not a
            # text annotation per tick); minor ticks off.
            ax.spines["left"].set_position(("data", 0.0))
            ax.spines["left"].set_color("#111827")
            ax.spines["left"].set_linewidth(1.2)
            ax.grid(False)
            ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=9))
            ax.yaxis.set_minor_locator(NullLocator())
            ax.yaxis.set_major_formatter(FuncFormatter(_plain_log_tick))
            ax.tick_params(axis="y", which="major", left=True, right=False, labelleft=True, length=6, width=1.0,
                           colors="#111827", labelsize=10, pad=4)
            ax.set_xlabel("Z-Score", fontsize=12, color="#111827", labelpad=10)
            ax.set_ylabel("")
            fig.tight_layout()
            fig.canvas.draw()

            # Draw term labels + single straight 2-point leader lines now that
            # transData reflects the final log-scaled axis. The leader endpoint
            # is the label anchor (dot position offset by `annotation_offsets`
            # points); computing it in data space keeps the line aligned through
            # the tight-bbox save. Each leader is one Line2D (2 anchors only).
            # Labels never overlap: they form one column at the right of the data, ordered by the dot's
            # y position (most significant at the top), evenly spaced in axes fraction, each joined to its
            # dot by one straight 2-point Line2D. The column sits past the largest z so no label covers a dot.
            data_inv = ax.transData.inverted()
            if labels:
                # top of the column = the dot nearest the top of the axes (largest plotted p), so leaders do not cross
                labels = sorted(labels, key=lambda l: (-float(l["fdr_plot"]), float(l["z_score"])))
                x_col = float(np.max(x_values)) + 0.06 * (x_max - x_min)
                if x_col > x_max - 0.02 * (x_max - x_min):
                    x_max = x_col + 0.02 * (x_max - x_min)
                    ax.set_xlim(x_min, x_max)
                    fig.canvas.draw()
                n = len(labels)
                step = min(0.045, 0.96 / max(n, 1))
                if n * 0.0125 > 0.96:  # 7 pt text needs about 1.25% of a 7.8 in axes per line; grow the figure
                    fig.set_size_inches(9.2, 7.8 * (n * 0.0125 / 0.96))
                    fig.tight_layout(); fig.canvas.draw()
                    step = 0.96 / n
                top = 0.5 + step * (n - 1) / 2.0
                axes_to_data = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
                for index, label in enumerate(labels):
                    zx, fy = float(label["z_score"]), float(label["fdr_plot"])
                    y_frac = top - index * step
                    disp = axes_to_data.transform((x_col, y_frac))
                    x1, y1 = data_inv.transform(disp)
                    ax.plot([zx, x1], [fy, y1], color="#9CA3AF", linewidth=0.6, solid_capstyle="butt", zorder=3)
                    ax.text(x1, y1, "  " + str(label["term_name"]), ha="left", va="center", fontsize=7,
                            color=str(label["label_color"]), zorder=5, clip_on=False)

            guide_transform = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
            guide_display_xy = guide_transform.transform((0.0, 0.5))
            guide_figure_xy = fig.transFigure.inverted().transform(guide_display_xy)
            ax.annotate(
                y_axis_label,
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
            # Illustrator-ready: no clip groups, no background rectangles (figure_style.strip_masks).
            from altanalyze3.components.visualization.figure_style import strip_masks
            strip_masks(fig)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    return output_path
