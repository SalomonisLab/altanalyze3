#!/usr/bin/env python3
"""Shared matplotlib settings for Illustrator-ready vector PDFs, and the hematopoietic state order.

Every figure module in this package that writes a PDF for a manuscript uses these helpers so the
output holds live Arial text (TrueType, fonttype 42), no clipping masks, no background boxes, no
raster stamps, and every point as a real vector circle or triangle.

    apply_rcparams()            before the figure is created
    vector_circles(ax, ...)     points as an EllipseCollection in point units (no marker stamps)
    vector_triangles(ax, ...)   points as a RegularPolyCollection in point units
    strip_masks(fig)            after drawing, before fig.savefig(path, format="pdf")
    pdf_vector_report(path)     counts raster stamps and embedded TrueType fonts in a written PDF

``HEMATOPOIETIC_STATE_ORDER`` is the authoritative bone-marrow cell-state order supplied by Nathan
Salomonis on 2026-06-06. Modules that order cell states on an axis use it and drop absent states.
"""

from __future__ import annotations

from typing import Iterable, List, Sequence

import matplotlib as mpl
import numpy as np
from matplotlib.collections import EllipseCollection, RegularPolyCollection
from matplotlib.colors import to_rgba
from matplotlib.patches import Rectangle

HEMATOPOIETIC_STATE_ORDER: List[str] = [
    "HSC-1", "HSC-2", "MPP-1", "MPP-2", "MPP-MEP", "LMPP-1", "LMPP-1-cycling", "Multilin-1",
    "Multilin-2", "MEP-1", "MEP-2", "MEP-Eryth-1", "MEP-Eryth-2", "ERP-1", "ERP-2", "ERP-3",
    "ERP-4", "ERP-5", "ERP-6", "ERP-7", "ERP-8", "Erythroblast-1", "Erythroblast-2",
    "Erythroblast-3", "MKP-early", "MKP-late", "MK-Platelet", "BMCP-1", "BMCP-2", "LMPP-2", "CLP",
    "Pro-B-Early", "Pro-B-1", "Pro-B-cycling-2", "Pro-B-Early-cycling", "Pro-B-cycling-1",
    "Transitional-B-1", "Transitional-B-2", "Pro-B-2", "Pro-B-3", "pre-B", "B Memory-1",
    "B Memory-2", "Plasma Cell", "CD4 TCM", "T CD4 Activated", "T CD4 Naive-1", "T CD4 Naive-2",
    "T CD8 Naive", "CD8 TEM", "MAIT", "NK-Mature", "MultiLin-GMP-1", "MultiLin-GMP-2",
    "MultiLin-GMP-3", "preNeu", "immNeu-1", "immNeu-2", "cMOP", "Mono-1", "Mono-2",
    "Intermediate Mono-1", "Intermediate Mono-2", "Intermediate Mono-3", "Classical-Mono",
    "Non-Classical Mono-1", "Non-Classical Mono-2", "Myeloid intermediate 1",
    "Myeloid intermediate 2", "Myeloid intermediate 3", "MDP-1", "MDP-2", "MDP-3", "MDP-4",
    "MDP-5", "MDP-6", "pre-DC-1", "pre-DC-2", "pre-DC-3", "cDC1", "cDC2-1", "cDC2-2", "pDC",
    "ASDC", "Mac", "MSC Fibroblast-1", "MSC Fibroblast-2", "Osteoblast", "Stromal Vascular",
]


def order_states(states: Iterable[str], order: Sequence[str] | None = None) -> List[str]:
    """Return ``states`` in ``order`` (default: the hematopoietic order); unknown states follow, sorted."""
    ref = list(order) if order is not None else HEMATOPOIETIC_STATE_ORDER
    present = set(states)
    known = [s for s in ref if s in present]
    unknown = sorted(present - set(ref))
    return known + unknown


def apply_rcparams(base_size: float = 7) -> None:
    mpl.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": base_size, "axes.titlesize": base_size + 1, "axes.labelsize": base_size,
        "xtick.labelsize": base_size, "ytick.labelsize": base_size, "legend.fontsize": base_size,
        "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2, "ytick.major.size": 2,
        "pdf.compression": 6, "path.simplify": False,
        "svg.fonttype": "none",
        "figure.facecolor": "white",
    })


def _rgba_array(colors, n: int, alpha: float) -> np.ndarray:
    """One opaque-or-baked RGBA row per point. ``colors`` is one colour or a sequence of n colours."""
    if isinstance(colors, str) or (isinstance(colors, tuple) and len(colors) in (3, 4)
                                   and not isinstance(colors[0], (str, tuple, list))):
        rgba = np.tile(np.asarray(to_rgba(colors), dtype=float), (n, 1))
    else:
        rgba = np.asarray([to_rgba(c) for c in colors], dtype=float)
        if rgba.shape[0] != n:
            raise ValueError(f"{rgba.shape[0]} colours for {n} points")
    rgba[:, 3] = float(alpha)
    return rgba


def vector_circles(ax, x, y, diameter_pt, colors, alpha: float = 1.0, zorder: int = 3):
    """Draw points as real circles of ``diameter_pt`` points (scalar or per point)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = x.size
    if n == 0:
        return None
    d = np.broadcast_to(np.asarray(diameter_pt, dtype=float), (n,)).copy()
    coll = EllipseCollection(
        widths=d, heights=d, angles=np.zeros(n), units="points",
        offsets=np.column_stack((x, y)), offset_transform=ax.transData,
        facecolors=_rgba_array(colors, n, alpha), edgecolors="none", linewidths=0.0,
        zorder=zorder,
    )
    ax.add_collection(coll)
    return coll


def vector_triangles(ax, x, y, side_pt, colors, alpha: float = 1.0, zorder: int = 3):
    """Draw points as upward triangles whose area matches a circle of diameter ``side_pt``."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = x.size
    if n == 0:
        return None
    area = np.pi * (float(side_pt) / 2.0) ** 2
    coll = RegularPolyCollection(
        numsides=3, rotation=0.0, sizes=np.full(n, area),
        offsets=np.column_stack((x, y)), offset_transform=ax.transData,
        facecolors=_rgba_array(colors, n, alpha), edgecolors="none", linewidths=0.0,
        zorder=zorder,
    )
    ax.add_collection(coll)
    return coll


def strip_masks(fig):
    """Remove what Illustrator shows as clipping masks, empty groups and background boxes."""
    fig.patch.set_visible(False)
    for ax in fig.axes:
        ax.patch.set_visible(False)
        for art in ax.get_children():
            if hasattr(art, "set_clip_on"):
                art.set_clip_on(False)
        for line in ax.get_xgridlines() + ax.get_ygridlines():
            line.set_clip_on(False)
        for p in list(ax.patches):
            if isinstance(p, Rectangle) and (p.get_width() == 0 or p.get_height() == 0):
                p.remove()
        leg = ax.get_legend()
        if leg is not None:
            leg.get_frame().set_visible(False)
    return fig


def pdf_vector_report(path) -> dict:
    """Count raster image stamps and embedded TrueType fonts in a PDF (both must be checked)."""
    with open(path, "rb") as fh:
        b = fh.read()
    report = {
        "raster_images": b.count(b"/Subtype/Image") + b.count(b"/Subtype /Image"),
        "truetype_fonts": b.count(b"/FontFile2"),
        "bytes": len(b),
    }
    print(f"[pdf] {path}: raster stamps {report['raster_images']}, "
          f"embedded TrueType fonts {report['truetype_fonts']}, {report['bytes']} bytes")
    return report
