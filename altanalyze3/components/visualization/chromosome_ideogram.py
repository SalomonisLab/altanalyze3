#!/usr/bin/env python3
"""Banded chromosome ideograms (vector SVG) for inter-chromosomal trans-splicing junctions.

One SVG per chimeric RNA. The donor chromosome is drawn on the left and the acceptor chromosome on the
right, both to one bp scale, with GRCh38 Giemsa bands (UCSC cytoBand). A tick marks each junction
coordinate, an arrow joins donor to acceptor, and a ribbon joins the two DNA contact bins of the best
Hi-C or ChIA-PET association.

Input is the per-junction table (one row per junction), for example
``ALL_4319_JUNCTIONS_WITHOUT_own_hic.tsv``. Columns used:

    chimeric_rna, donor_gene, acceptor_gene, junction ("chrA:pos:strand > chrB:pos:strand")
    overall_bestq_* and overall_bestp_*: source, p, q, z, dataset, cell_type_or_tissue,
        donor_side_contact, acceptor_side_contact ("chr:start-end", 1-based), significant
    kinnex_junction_reads_sum_all_samples, kinnex_samples_positive, kinnex_patients_positive,
    ogawa_junction_fragments_sum_all_362_samples, ogawa_samples_positive_of_362 (caption only)

Contact row rule (Nathan Salomonis, 2026-11-10): the best-by-q row when ``overall_bestq_significant``
is ``q<0.05``; otherwise the best-by-p row; no ribbon when the table holds no best association.

The SVG holds live text (Arial), 2-point lines, filled paths and no clip paths, masks, patterns or
images. Hatching (acen, gvar) is drawn as line segments clipped analytically to each band.

    python -m altanalyze3 visualize chromosome-ideogram \
        --junctions ALL_4319_JUNCTIONS_WITHOUT_own_hic.tsv --chimera "FOSB>SF1" --outdir OUT [--layout vertical]

Layouts: ``horizontal`` (default) is a compact panel for placement under an ISV figure; ``vertical`` draws the two
chromosomes side by side with every band labelled.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import sys
import textwrap
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape

import pandas as pd

RESOURCE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")
DEFAULT_CYTOBANDS = os.path.join(RESOURCE_DIR, "cytoBand_hg38_UCSC.txt.gz")

STAIN_FILL = {
    "gneg": "#FFFFFF", "gpos25": "#D9D9D9", "gpos50": "#A6A6A6", "gpos75": "#737373",
    "gpos100": "#404040", "gvar": "#FFFFFF", "acen": "#FFFFFF", "stalk": "#FFFFFF",
}
HATCHED = {"gvar", "acen"}
DONOR_COLOR = "#1F78B4"
ACCEPTOR_COLOR = "#E31A1C"
CONTACT_COLOR = "#33A02C"
INK = "#000000"
GREY_TEXT = "#4D4D4D"
FONT = "Arial, Helvetica, sans-serif"

REQUIRED_COLUMNS = ["chimeric_rna", "donor_gene", "acceptor_gene", "junction"]
BEST_FIELDS = ["source", "p", "q", "z", "dataset", "cell_type_or_tissue", "donor_side_contact",
               "acceptor_side_contact", "tied_datasets", "significant"]


# ----------------------------------------------------------------------------------------------- data
@dataclass
class Band:
    chrom: str
    start: int  # 0-based
    end: int
    name: str
    stain: str


@dataclass
class JunctionEnd:
    gene: str
    chrom: str
    pos: int  # 1-based
    strand: str
    band: str = ""


@dataclass
class Contact:
    rule: str
    source: str
    dataset: str
    tissue: str
    p: Optional[float]
    q: Optional[float]
    z: Optional[float]
    donor_bin: Tuple[str, int, int]  # chrom, 1-based start, end
    acceptor_bin: Tuple[str, int, int]
    significant: str
    tied_datasets: str = ""


def md5sum(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_cytobands(path: str = DEFAULT_CYTOBANDS) -> Dict[str, List[Band]]:
    """Read a UCSC cytoBand table and check that each chromosome is tiled from 0 without gaps."""
    opener = gzip.open if path.endswith(".gz") else open
    out: Dict[str, List[Band]] = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 5:
                raise ValueError(f"cytoband line with {len(f)} fields: {line!r}")
            out.setdefault(f[0], []).append(Band(f[0], int(f[1]), int(f[2]), f[3], f[4]))
    for chrom, bands in out.items():
        bands.sort(key=lambda b: b.start)
        if bands[0].start != 0:
            raise ValueError(f"{chrom}: first band starts at {bands[0].start}, not 0")
        for a, b in zip(bands, bands[1:]):
            if a.end != b.start:
                raise ValueError(f"{chrom}: bands {a.name} and {b.name} do not abut ({a.end} vs {b.start})")
        for b in bands:
            if b.end <= b.start:
                raise ValueError(f"{chrom}: band {b.name} has no length")
            if b.stain not in STAIN_FILL:
                raise ValueError(f"{chrom}: band {b.name} has unknown stain {b.stain}")
    return out


def band_at(bands: Sequence[Band], pos: int) -> Band:
    """Band holding the 1-based coordinate ``pos`` (UCSC bands are 0-based, half-open)."""
    for b in bands:
        if b.start < pos <= b.end:
            return b
    raise ValueError(f"{bands[0].chrom}:{pos} lies outside the band table (length {bands[-1].end})")


def iscn(chrom: str, band: str) -> str:
    return chrom.replace("chr", "") + band


_JUNCTION = re.compile(r"^(chr[^:]+):(\d+):([+-])\s*>\s*(chr[^:]+):(\d+):([+-])$")
_REGION = re.compile(r"^(chr[^:]+):(\d+)-(\d+)$")


def parse_junction(text: str) -> Tuple[Tuple[str, int, str], Tuple[str, int, str]]:
    m = _JUNCTION.match(str(text).strip())
    if not m:
        raise ValueError(f"junction not in 'chrA:pos:strand > chrB:pos:strand' form: {text!r}")
    return (m.group(1), int(m.group(2)), m.group(3)), (m.group(4), int(m.group(5)), m.group(6))


def parse_region(text: str) -> Tuple[str, int, int]:
    m = _REGION.match(str(text).strip())
    if not m:
        raise ValueError(f"contact region not in 'chr:start-end' form: {text!r}")
    start, end = int(m.group(2)), int(m.group(3))
    if end < start:
        raise ValueError(f"contact region end before start: {text!r}")
    return m.group(1), start, end


def _num(v) -> Optional[float]:
    if v is None or (isinstance(v, float) and math.isnan(v)) or str(v).strip() == "":
        return None
    return float(v)


def select_contact(row: pd.Series) -> Optional[Contact]:
    """Best-by-q row when it passes q < 0.05, else best-by-p row, else None."""
    def has(prefix):
        return str(row.get(f"{prefix}_source", "") or "").strip() not in ("", "nan")

    if has("overall_bestq") and str(row.get("overall_bestq_significant")) == "q<0.05":
        prefix, rule = "overall_bestq", "best by q (q < 0.05)"
    elif has("overall_bestp"):
        prefix, rule = "overall_bestp", "best by p (no row with q < 0.05)"
    else:
        return None
    g = {k: row.get(f"{prefix}_{k}") for k in BEST_FIELDS}
    return Contact(rule=rule, source=str(g["source"]), dataset=str(g["dataset"]),
                   tissue=str(g["cell_type_or_tissue"]), p=_num(g["p"]), q=_num(g["q"]), z=_num(g["z"]),
                   donor_bin=parse_region(g["donor_side_contact"]),
                   acceptor_bin=parse_region(g["acceptor_side_contact"]),
                   significant=str(g["significant"]), tied_datasets=str(g["tied_datasets"]))


# ------------------------------------------------------------------------------------------- geometry
def clip_line_to_convex(poly: Sequence[Tuple[float, float]], c: float) -> Optional[Tuple[Tuple[float, float], Tuple[float, float]]]:
    """Segment of the line x + y = c inside a convex polygon, or None."""
    pts = []
    n = len(poly)
    for i in range(n):
        (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % n]
        f1, f2 = x1 + y1 - c, x2 + y2 - c
        if f1 == 0:
            pts.append((x1, y1))
        if (f1 < 0 < f2) or (f2 < 0 < f1):
            t = f1 / (f1 - f2)
            pts.append((x1 + t * (x2 - x1), y1 + t * (y2 - y1)))
    if len(pts) < 2:
        return None
    pts.sort()
    a, b = pts[0], pts[-1]
    if math.hypot(b[0] - a[0], b[1] - a[1]) < 0.2:
        return None
    return a, b


def point_in_convex(poly: Sequence[Tuple[float, float]], p: Tuple[float, float], tol: float = 1e-6) -> bool:
    sign = 0
    n = len(poly)
    for i in range(n):
        (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % n]
        cross = (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1)
        if abs(cross) <= tol * max(1.0, math.hypot(x2 - x1, y2 - y1)):
            continue
        s = 1 if cross > 0 else -1
        if sign == 0:
            sign = s
        elif s != sign:
            return False
    return True


def hatch_segments(poly, spacing=2.4):
    """45-degree hatch lines clipped to a convex polygon."""
    xs = [p[0] + p[1] for p in poly]
    lo, hi = min(xs), max(xs)
    segs = []
    c = math.floor(lo / spacing) * spacing + spacing / 2
    while c < hi:
        seg = clip_line_to_convex(poly, c)
        if seg:
            segs.append(seg)
        c += spacing
    return segs


def _arc_points(cx, cy, r, a0, a1, n=8):
    return [(cx + r * math.cos(a0 + (a1 - a0) * k / n), cy + r * math.sin(a0 + (a1 - a0) * k / n)) for k in range(n + 1)]


@dataclass
class Ideogram:
    """Band shapes of one chromosome.

    ``cx`` is the centre across the chromosome and ``y_top`` the coordinate of base 0 along it. For
    ``orientation="vertical"`` these are x and y; for ``"horizontal"`` they are y and x. ``y(pos)`` returns the
    coordinate along the chromosome in either orientation.
    """
    chrom: str
    bands: List[Band]
    cx: float
    y_top: float
    scale: float  # px per bp
    width: float = 20.0
    orientation: str = "vertical"
    shapes: List[dict] = field(default_factory=list)

    def __post_init__(self):
        if self.orientation not in ("vertical", "horizontal"):
            raise ValueError(f"orientation must be vertical or horizontal, not {self.orientation!r}")

    def P(self, across: float, along: float) -> Tuple[float, float]:
        return (across, along) if self.orientation == "vertical" else (along, across)

    def _pt(self, across: float, along: float) -> str:
        x, y = self.P(across, along)
        return f"{x:.3f},{y:.3f}"

    @property
    def _sweep(self) -> int:
        # the horizontal layout mirrors the vertical one across the diagonal, which reverses arc direction
        return 1 if self.orientation == "vertical" else 0

    def y(self, pos_bp: float) -> float:
        return self.y_top + pos_bp * self.scale

    @property
    def length(self) -> int:
        return self.bands[-1].end

    @property
    def y_bottom(self) -> float:
        return self.y(self.length)

    def half_widths(self, b: Band) -> Tuple[float, float]:
        w = self.width / 2.0
        if b.stain == "stalk":
            return 0.4 * w, 0.4 * w
        if b.stain == "acen":
            return (w, 0.45 * w) if b.name.startswith("p") else (0.45 * w, w)
        return w, w

    def build(self) -> "Ideogram":
        n = len(self.bands)
        heights = [(b.end - b.start) * self.scale for b in self.bands]
        r = min(0.35 * self.width, heights[0], heights[-1])
        self.radius = r
        self.shapes = []
        pt, sw = self._pt, self._sweep
        for i, b in enumerate(self.bands):
            ys, ye = self.y(b.start), self.y(b.end)
            ht, hb = self.half_widths(b)
            xl_t, xr_t, xl_b, xr_b = self.cx - ht, self.cx + ht, self.cx - hb, self.cx + hb
            first, last = i == 0, i == n - 1
            d = []
            poly: List[Tuple[float, float]] = []
            if first:
                d.append(f"M{pt(xl_t, ys + r)} A{r:.3f},{r:.3f} 0 0 {sw} {pt(xl_t + r, ys)} "
                         f"L{pt(xr_t - r, ys)} A{r:.3f},{r:.3f} 0 0 {sw} {pt(xr_t, ys + r)}")
                poly += [self.P(*q) for q in _arc_points(xl_t + r, ys + r, r, math.pi, 1.5 * math.pi)]
                poly += [self.P(*q) for q in _arc_points(xr_t - r, ys + r, r, 1.5 * math.pi, 2 * math.pi)]
            else:
                d.append(f"M{pt(xl_t, ys)} L{pt(xr_t, ys)}")
                poly += [self.P(xl_t, ys), self.P(xr_t, ys)]
            if last:
                d.append(f"L{pt(xr_b, ye - r)} A{r:.3f},{r:.3f} 0 0 {sw} {pt(xr_b - r, ye)} "
                         f"L{pt(xl_b + r, ye)} A{r:.3f},{r:.3f} 0 0 {sw} {pt(xl_b, ye - r)} Z")
                poly += [self.P(*q) for q in _arc_points(xr_b - r, ye - r, r, 0, 0.5 * math.pi)]
                poly += [self.P(*q) for q in _arc_points(xl_b + r, ye - r, r, 0.5 * math.pi, math.pi)]
            else:
                d.append(f"L{pt(xr_b, ye)} L{pt(xl_b, ye)} Z")
                poly += [self.P(xr_b, ye), self.P(xl_b, ye)]
            self.shapes.append({"band": b, "d": " ".join(d), "poly": poly, "ys": ys, "ye": ye})
        return self

    def outline_path(self) -> str:
        r = self.radius
        right, left = [], []
        for s in self.shapes:
            ht, hb = self.half_widths(s["band"])
            right += [(self.cx + ht, s["ys"]), (self.cx + hb, s["ye"])]
            left += [(self.cx - ht, s["ys"]), (self.cx - hb, s["ye"])]
        pt, sw = self._pt, self._sweep
        (xr0, y0), (xrn, yn) = right[0], right[-1]
        xl0, xln = left[0][0], left[-1][0]
        d = [f"M{pt(xl0, y0 + r)} A{r:.3f},{r:.3f} 0 0 {sw} {pt(xl0 + r, y0)} L{pt(xr0 - r, y0)} "
             f"A{r:.3f},{r:.3f} 0 0 {sw} {pt(xr0, y0 + r)}"]
        for x, y in right[1:-1]:
            d.append(f"L{pt(x, y)}")
        d.append(f"L{pt(xrn, yn - r)} A{r:.3f},{r:.3f} 0 0 {sw} {pt(xrn - r, yn)} L{pt(xln + r, yn)} "
                 f"A{r:.3f},{r:.3f} 0 0 {sw} {pt(xln, yn - r)}")
        for x, y in reversed(left[1:-1]):
            d.append(f"L{pt(x, y)}")
        d.append("Z")
        return " ".join(d)

    def hatch_segments(self, poly, spacing=2.4):
        return hatch_segments(poly, spacing)


def spread_labels(desired: Sequence[float], gap: float, lo: float, hi: float) -> List[float]:
    """Move label centres apart by at least ``gap``, keeping order, inside [lo, hi] where possible."""
    n = len(desired)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: desired[i])
    ys = [desired[i] for i in order]
    for _ in range(200):
        moved = False
        for k in range(1, n):
            if ys[k] - ys[k - 1] < gap - 1e-9:
                mid = (ys[k] + ys[k - 1]) / 2
                ys[k - 1], ys[k] = mid - gap / 2, mid + gap / 2
                moved = True
        if ys[0] < lo:
            shift = lo - ys[0]
            ys = [y + shift for y in ys]
        if ys[-1] > hi and ys[0] - (ys[-1] - hi) >= lo:
            shift = ys[-1] - hi
            ys = [y - shift for y in ys]
        if not moved:
            break
    for k in range(1, n):
        if ys[k] - ys[k - 1] < gap - 1e-6:
            ys[k] = ys[k - 1] + gap
    out = [0.0] * n
    for k, i in enumerate(order):
        out[i] = ys[k]
    return out


# ------------------------------------------------------------------------------------------- drawing
def _fmt_stat(v: Optional[float]) -> str:
    if v is None:
        return "n.a."
    return f"{v:.2e}" if v < 1e-3 else f"{v:.3f}"


def _int(v) -> Optional[int]:
    x = _num(v)
    return None if x is None else int(round(x))


class SVG:
    def __init__(self):
        self.items: List[str] = []

    def add(self, s: str):
        self.items.append(s)

    def line(self, x1, y1, x2, y2, stroke=INK, width=0.5, cap="butt"):
        self.add(f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" stroke="{stroke}" '
                 f'stroke-width="{width}" stroke-linecap="{cap}"/>')

    def path(self, d, fill="none", stroke="none", width=0.0, fill_opacity=None, join="miter"):
        fo = f' fill-opacity="{fill_opacity}"' if fill_opacity is not None else ""
        sw = f' stroke-width="{width}" stroke-linejoin="{join}"' if stroke != "none" else ""
        self.add(f'<path d="{d}" fill="{fill}"{fo} stroke="{stroke}"{sw}/>')

    def text(self, x, y, s, size=7.0, anchor="start", weight="normal", fill=INK, style="normal"):
        self.add(f'<text x="{x:.3f}" y="{y:.3f}" font-family="{FONT}" font-size="{size}" font-weight="{weight}" '
                 f'font-style="{style}" fill="{fill}" text-anchor="{anchor}">{escape(s)}</text>')

    def render(self, width, height) -> str:
        head = (f'<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" version="1.1" '
                f'width="{width:.0f}" height="{height:.0f}" viewBox="0 0 {width:.0f} {height:.0f}">\n')
        return head + "\n".join(self.items) + "\n</svg>\n"


def text_width(s: str, size: float) -> float:
    return 0.52 * size * len(s)


def draw_bands(svg: SVG, ideo: Ideogram) -> None:
    """Band fills, acen and gvar hatching, and the chromosome outline, in either orientation."""
    for s in ideo.shapes:
        svg.path(s["d"], fill=STAIN_FILL[s["band"].stain])
    for s in ideo.shapes:
        if s["band"].stain in HATCHED:
            for (x1, y1), (x2, y2) in ideo.hatch_segments(s["poly"]):
                svg.line(x1, y1, x2, y2, width=0.4)
    svg.path(ideo.outline_path(), stroke=INK, width=0.75, join="round")


def draw_ideogram(svg: SVG, ideo: Ideogram, label_side: str, font=6.0):
    draw_bands(svg, ideo)
    edge = ideo.cx - ideo.width / 2 if label_side == "left" else ideo.cx + ideo.width / 2
    sign = -1 if label_side == "left" else 1
    mids = [(s["ys"] + s["ye"]) / 2 for s in ideo.shapes]
    ys = spread_labels(mids, gap=font * 1.08, lo=ideo.y_top - 2, hi=ideo.y_bottom + 2)
    anchor = "end" if label_side == "left" else "start"
    for s, ym, yl in zip(ideo.shapes, mids, ys):
        name = s["band"].name[1:]
        if abs(yl - ym) > 0.6:
            svg.line(edge + sign * 1.0, ym, edge + sign * 6.0, yl, stroke=GREY_TEXT, width=0.35)
        svg.text(edge + sign * 7.5, yl + font * 0.35, name, size=font, anchor=anchor)
    p_bands = [s for s in ideo.shapes if s["band"].name.startswith("p")]
    q_bands = [s for s in ideo.shapes if s["band"].name.startswith("q")]
    arm_x = edge + sign * 32
    for arm, group in (("p", p_bands), ("q", q_bands)):
        if group:
            ya = (group[0]["ys"] + group[-1]["ye"]) / 2
            svg.text(arm_x, ya + 3, arm, size=8, anchor="middle", style="italic", fill=GREY_TEXT)
    svg.text(ideo.cx, ideo.y_bottom + 14, f"Chromosome {ideo.chrom.replace('chr', '')}", size=8, anchor="middle",
             weight="bold")


def _cubic_mid(p0, c1, c2, p3):
    return tuple(0.125 * a + 0.375 * b + 0.375 * c + 0.125 * d for a, b, c, d in zip(p0, c1, c2, p3))


def junction_ends(row: pd.Series, cytobands: Dict[str, List[Band]]):
    """Donor and acceptor ends with their bands, and the selected contact row (checked against the chromosomes)."""
    (dc, dp, ds), (ac, ap, as_) = parse_junction(row["junction"])
    if dc == ac:
        raise ValueError(f"{row['chimeric_rna']}: donor and acceptor on {dc}; this plotter draws two chromosomes")
    for c in (dc, ac):
        if c not in cytobands:
            raise ValueError(f"{c} is not in the band table")
    donor = JunctionEnd(str(row["donor_gene"]), dc, dp, ds, band_at(cytobands[dc], dp).name)
    acceptor = JunctionEnd(str(row["acceptor_gene"]), ac, ap, as_, band_at(cytobands[ac], ap).name)
    contact = select_contact(row)
    if contact is not None:
        if contact.donor_bin[0] != dc or contact.acceptor_bin[0] != ac:
            raise ValueError(f"{row['chimeric_rna']}: contact bins {contact.donor_bin[0]}/{contact.acceptor_bin[0]} "
                             f"are not on the junction chromosomes {dc}/{ac}")
    return donor, acceptor, contact


def make_values(row: pd.Series, donor: JunctionEnd, acceptor: JunctionEnd, contact: Optional[Contact],
                scale: float, layout: str) -> dict:
    dc, dp, ds = donor.chrom, donor.pos, donor.strand
    ac, ap, as_ = acceptor.chrom, acceptor.pos, acceptor.strand
    return {
        "chimeric_rna": str(row["chimeric_rna"]), "junction": str(row["junction"]), "layout": layout,
        "donor_gene": donor.gene, "donor_chrom": dc, "donor_pos": dp, "donor_strand": ds,
        "donor_band": iscn(dc, donor.band), "acceptor_gene": acceptor.gene, "acceptor_chrom": ac,
        "acceptor_pos": ap, "acceptor_strand": as_, "acceptor_band": iscn(ac, acceptor.band),
        "contact_rule": contact.rule if contact else "no best association in the table",
        "contact_source": contact.source if contact else "", "contact_dataset": contact.dataset if contact else "",
        "contact_tissue": contact.tissue if contact else "",
        "contact_p": contact.p if contact else None, "contact_q": contact.q if contact else None,
        "contact_z": contact.z if contact else None,
        "contact_donor_bin": f"{contact.donor_bin[0]}:{contact.donor_bin[1]}-{contact.donor_bin[2]}" if contact else "",
        "contact_acceptor_bin": f"{contact.acceptor_bin[0]}:{contact.acceptor_bin[1]}-{contact.acceptor_bin[2]}" if contact else "",
        "contact_significant": contact.significant if contact else "",
        "kinnex_junction_molecules": _int(row.get("kinnex_junction_reads_sum_all_samples")),
        "kinnex_samples_positive": _int(row.get("kinnex_samples_positive")),
        "kinnex_patients_positive": _int(row.get("kinnex_patients_positive")),
        "ogawa_fragments_362": _int(row.get("ogawa_junction_fragments_sum_all_362_samples")),
        "ogawa_samples_positive_of_362": _int(row.get("ogawa_samples_positive_of_362")),
        "scale_bp_per_px": 1.0 / scale,
    }


def draw_legend(svg: SVG, y: float, x0: float = 16.0) -> float:
    """Stain swatches, the junction arrow and the contact ribbon; returns the next free y."""
    x = x0
    for stain, label in (("gneg", "gneg"), ("gpos25", "gpos25"), ("gpos50", "gpos50"), ("gpos75", "gpos75"),
                         ("gpos100", "gpos100"), ("gvar", "acen, gvar")):
        svg.path(f"M{x:.3f},{y - 6:.3f} L{x + 10:.3f},{y - 6:.3f} L{x + 10:.3f},{y:.3f} L{x:.3f},{y:.3f} Z",
                 fill=STAIN_FILL[stain], stroke=INK, width=0.5)
        if stain in HATCHED:
            poly = [(x, y - 6), (x + 10, y - 6), (x + 10, y), (x, y)]
            for (a1, a2), (b1_, b2_) in hatch_segments(poly):
                svg.line(a1, a2, b1_, b2_, width=0.4)
        svg.text(x + 13, y - 0.5, label, size=6.5)
        x += 16 + text_width(label, 6.5) + 8
    y += 13
    x = x0
    svg.line(x, y - 3, x + 18, y - 3, stroke=INK, width=0.9)
    svg.path(f"M{x + 22:.3f},{y - 3:.3f} L{x + 17:.3f},{y - 5.4:.3f} L{x + 17:.3f},{y - 0.6:.3f} Z", fill=ACCEPTOR_COLOR)
    svg.text(x + 26, y - 0.5, "junction, donor to acceptor (tick at the junction coordinate)", size=6.5)
    x += 26 + text_width("junction, donor to acceptor (tick at the junction coordinate)", 6.5) + 14
    svg.path(f"M{x:.3f},{y - 6:.3f} L{x + 14:.3f},{y - 6:.3f} L{x + 14:.3f},{y:.3f} L{x:.3f},{y:.3f} Z",
             fill=CONTACT_COLOR, fill_opacity=0.35, stroke=CONTACT_COLOR, width=0.5)
    svg.text(x + 18, y - 0.5, "DNA contact bins", size=6.5)
    return y + 8


def draw_caption(svg: SVG, values: dict, y: float, width: float, x0: float = 16.0, band_labels: bool = True) -> float:
    """The written legend: claim, panel content, contact statistics, scope; returns the next free y."""
    lines = caption_text(values, band_labels=band_labels)
    y += 6
    for k, para in enumerate(lines):
        for i, ln in enumerate(textwrap.wrap(para, width=int((width - 2 * x0) / (0.50 * 6.5)))):
            svg.text(x0, y, ln, size=6.5, weight="bold" if k == 0 else "normal")
            y += 8.2
        y += 2
    return y


def render_chimera(row: pd.Series, cytobands: Dict[str, List[Band]], *, height=520.0, chrom_width=20.0,
                   stub=128.0, arc_span=110.0, caption=True) -> Tuple[str, dict]:
    """Vertical layout: donor chromosome left, acceptor chromosome right, band labels on the outer sides."""
    donor, acceptor, contact = junction_ends(row, cytobands)
    dc, dp, ds = donor.chrom, donor.pos, donor.strand
    ac, ap, as_ = acceptor.chrom, acceptor.pos, acceptor.strand

    top = 64.0
    scale = height / max(cytobands[dc][-1].end, cytobands[ac][-1].end)
    left_margin = 62.0
    gap = 2 * stub + arc_span
    cx_d = left_margin + chrom_width / 2
    cx_a = cx_d + chrom_width + gap
    width = cx_a + chrom_width / 2 + left_margin
    L = Ideogram(dc, cytobands[dc], cx_d, top, scale, chrom_width).build()
    R = Ideogram(ac, cytobands[ac], cx_a, top, scale, chrom_width).build()
    chrom_bottom = max(L.y_bottom, R.y_bottom)
    legend_top = chrom_bottom + 34

    svg = SVG()
    title = f"{donor.gene}>{acceptor.gene}"
    svg.text(width / 2, 18, title, size=10, anchor="middle", weight="bold")
    svg.text(width / 2, 30, f"Trans-splice junction ({donor.gene} {iscn(dc, donor.band)} to {acceptor.gene} "
             f"{iscn(ac, acceptor.band)}) and best DNA contact, GRCh38", size=7, anchor="middle", fill=GREY_TEXT)

    draw_ideogram(svg, L, "left")
    draw_ideogram(svg, R, "right")

    edge_d, edge_a = cx_d + chrom_width / 2, cx_a - chrom_width / 2
    yD, yA = L.y(dp), R.y(ap)

    # junction bands outlined in the gene colour
    for ideo, end, colour in ((L, donor, DONOR_COLOR), (R, acceptor, ACCEPTOR_COLOR)):
        s = next(s for s in ideo.shapes if s["band"].name == end.band)
        svg.path(s["d"], stroke=colour, width=1.1)

    # DNA contact ribbon, bowing below the junction
    contact_drawn = None
    if contact is not None:
        b1 = (L.y(contact.donor_bin[1] - 1), L.y(contact.donor_bin[2]))
        b2 = (R.y(contact.acceptor_bin[1] - 1), R.y(contact.acceptor_bin[2]))
        x1, x2 = edge_d + 3.0, edge_a - 3.0
        drop = min(70.0, legend_top - 18 - max(b1[1], b2[1]))
        drop = max(drop, 12.0)
        cA1 = (x1 + gap * 0.30, max(b1[1], b2[1]) + drop)
        cA2 = (x2 - gap * 0.30, max(b1[1], b2[1]) + drop)
        d = (f"M{x1:.3f},{b1[0]:.3f} C{cA1[0]:.3f},{cA1[1] - (b1[1] - b1[0]):.3f} {cA2[0]:.3f},{cA2[1] - (b2[1] - b2[0]):.3f} "
             f"{x2:.3f},{b2[0]:.3f} L{x2:.3f},{b2[1]:.3f} C{cA2[0]:.3f},{cA2[1]:.3f} {cA1[0]:.3f},{cA1[1]:.3f} "
             f"{x1:.3f},{b1[1]:.3f} Z")
        svg.path(d, fill=CONTACT_COLOR, fill_opacity=0.35, stroke=CONTACT_COLOR, width=0.5)
        for x0, (y0, y1) in ((edge_d, b1), (edge_a - 3.0, b2)):
            svg.path(f"M{x0:.3f},{y0:.3f} L{x0 + 3:.3f},{y0:.3f} L{x0 + 3:.3f},{y1:.3f} L{x0:.3f},{y1:.3f} Z",
                     fill=CONTACT_COLOR)
        mid = _cubic_mid((x1, b1[1]), cA1, cA2, (x2, b2[1]))
        src = "ChIA-PET" if "chia" in contact.source.lower() else "Hi-C"
        svg.text(mid[0], mid[1] + 10, f"{src} contact: {contact.tissue}", size=6.5, anchor="middle", fill=INK)
        svg.text(mid[0], mid[1] + 18, f"p = {_fmt_stat(contact.p)}, q = {_fmt_stat(contact.q)}", size=6.5,
                 anchor="middle", fill=GREY_TEXT)
        contact_drawn = {"donor_bin_y": b1, "acceptor_bin_y": b2}

    # junction stubs, gene labels and the donor-to-acceptor arrow
    sx_d, sx_a = edge_d + stub, edge_a - stub
    svg.line(edge_d, yD, sx_d, yD, stroke=DONOR_COLOR, width=1.0)
    svg.line(sx_a, yA, edge_a - 4.0, yA, stroke=ACCEPTOR_COLOR, width=1.0)
    svg.path(f"M{edge_a:.3f},{yA:.3f} L{edge_a - 5:.3f},{yA - 2.4:.3f} L{edge_a - 5:.3f},{yA + 2.4:.3f} Z",
             fill=ACCEPTOR_COLOR)
    lift = min(60.0, min(yD, yA) - 44.0)
    lift = max(lift, 8.0)
    c1 = (sx_d + arc_span * 0.35, min(yD, yA) - lift)
    c2 = (sx_a - arc_span * 0.35, min(yD, yA) - lift)
    svg.path(f"M{sx_d:.3f},{yD:.3f} C{c1[0]:.3f},{c1[1]:.3f} {c2[0]:.3f},{c2[1]:.3f} {sx_a:.3f},{yA:.3f}",
             stroke=INK, width=0.9)
    top_mid = _cubic_mid((sx_d, yD), c1, c2, (sx_a, yA))
    svg.text(top_mid[0], top_mid[1] - 4, "trans-splice junction", size=6.5, anchor="middle", fill=INK)

    for end, x, anchor, colour, y in ((donor, edge_d + 4, "start", DONOR_COLOR, yD),
                                      (acceptor, edge_a - 7, "end", ACCEPTOR_COLOR, yA)):
        role = "donor" if end is donor else "acceptor"
        svg.text(x, y - 12, f"{end.gene}  {iscn(end.chrom, end.band)}", size=7.5, anchor=anchor, weight="bold",
                 fill=colour)
        svg.text(x, y - 3, f"{role}, {end.chrom}:{end.pos:,} ({end.strand})", size=6.5, anchor=anchor, fill=GREY_TEXT)

    # legend: stains and marks
    y = draw_legend(svg, legend_top)

    values = make_values(row, donor, acceptor, contact, scale, "vertical")

    if caption:
        y = draw_caption(svg, values, y, width)
    height_total = y + 8
    return svg.render(width, height_total), values


def render_chimera_horizontal(row: pd.Series, cytobands: Dict[str, List[Band]], *, width=720.0, thickness=14.0,
                              gap=56.0, title=True, legend=True, caption=True) -> Tuple[str, dict]:
    """Compact layout for placement under an ISV panel.

    Acceptor chromosome on top and donor chromosome below, both horizontal from p-ter (left) at one bp scale,
    as in the 2026-11-10 mock-up. The acceptor gene label sits above its band and the donor gene label below.
    The junction runs from the donor tick to an arrowhead on the acceptor band; the contact ribbon joins the
    two contact bins. ``width`` sets the drawing width so it can match the ISV figure.
    """
    donor, acceptor, contact = junction_ends(row, cytobands)
    dc, ac = donor.chrom, acceptor.chrom
    name_x = 6.0 + text_width("Chr 22", 8.0)
    left = name_x + 16.0
    right = 12.0
    scale = (width - left - right) / max(cytobands[dc][-1].end, cytobands[ac][-1].end)

    svg = SVG()
    y = 4.0
    if title:
        svg.text(width / 2, y + 11, f"{donor.gene}>{acceptor.gene}", size=9, anchor="middle", weight="bold")
        svg.text(width / 2, y + 21, f"Trans-splice junction ({donor.gene} {iscn(dc, donor.band)} to {acceptor.gene} "
                 f"{iscn(ac, acceptor.band)}) and best DNA contact, GRCh38", size=6.5, anchor="middle", fill=GREY_TEXT)
        y += 26
    cy_top = y + 10 + thickness / 2
    cy_bot = cy_top + thickness + gap
    T = Ideogram(ac, cytobands[ac], cy_top, left, scale, thickness, "horizontal").build()
    B = Ideogram(dc, cytobands[dc], cy_bot, left, scale, thickness, "horizontal").build()
    for ideo in (T, B):
        draw_bands(svg, ideo)
        svg.text(left - 6, ideo.cx + 2.6, "p", size=7.5, anchor="middle", style="italic", fill=GREY_TEXT)
        svg.text(ideo.y_bottom + 6, ideo.cx + 2.6, "q", size=7.5, anchor="middle", style="italic", fill=GREY_TEXT)
        svg.text(name_x, ideo.cx + 2.9, f"Chr {ideo.chrom.replace('chr', '')}", size=8, anchor="end", weight="bold")
    for ideo, end, colour in ((T, acceptor, ACCEPTOR_COLOR), (B, donor, DONOR_COLOR)):
        s = next(s for s in ideo.shapes if s["band"].name == end.band)
        svg.path(s["d"], stroke=colour, width=1.1)

    xA, xD = T.y(acceptor.pos), B.y(donor.pos)
    yT, yB = cy_top + thickness / 2, cy_bot - thickness / 2

    def centred(x, s, size):
        half = text_width(s, size) / 2
        return min(max(x, 2 + half), width - 2 - half)

    lab_a = f"{acceptor.gene}  {iscn(ac, acceptor.band)}"
    lab_d = f"{donor.gene}  {iscn(dc, donor.band)}"
    svg.text(centred(xA, lab_a, 7), cy_top - thickness / 2 - 3, lab_a, size=7, anchor="middle", weight="bold",
             fill=ACCEPTOR_COLOR)
    svg.text(centred(xD, lab_d, 7), cy_bot + thickness / 2 + 8.5, lab_d, size=7, anchor="middle", weight="bold",
             fill=DONOR_COLOR)

    k = 0.55 * gap
    if contact is not None:
        def span(ideo, b):
            u0, u1 = ideo.y(b[1] - 1), ideo.y(b[2])
            if u1 - u0 < 1.0:
                m = (u0 + u1) / 2
                u0, u1 = m - 0.5, m + 0.5
            return u0, u1
        d0, d1 = span(B, contact.donor_bin)
        a0, a1 = span(T, contact.acceptor_bin)
        svg.path(f"M{d0:.3f},{yB:.3f} C{d0:.3f},{yB - k:.3f} {a0:.3f},{yT + k:.3f} {a0:.3f},{yT:.3f} "
                 f"L{a1:.3f},{yT:.3f} C{a1:.3f},{yT + k:.3f} {d1:.3f},{yB - k:.3f} {d1:.3f},{yB:.3f} Z",
                 fill=CONTACT_COLOR, fill_opacity=0.35, stroke=CONTACT_COLOR, width=0.5)
        for (u0, u1), yy in (((d0, d1), yB - 2.0), ((a0, a1), yT)):
            svg.path(f"M{u0:.3f},{yy:.3f} L{u1:.3f},{yy:.3f} L{u1:.3f},{yy + 2:.3f} L{u0:.3f},{yy + 2:.3f} Z",
                     fill=CONTACT_COLOR)

    base = yT + 5.0
    svg.path(f"M{xD - 0.6:.3f},{yB - 3:.3f} L{xD + 0.6:.3f},{yB - 3:.3f} L{xD + 0.6:.3f},{yB:.3f} L{xD - 0.6:.3f},{yB:.3f} Z",
             fill=DONOR_COLOR)
    # Same end heights and control points as the ribbon edges, so at every height the junction line lies between
    # the two edges whenever the contact bins hold the junction coordinates; the arrowhead covers the last 5 units.
    svg.path(f"M{xD:.3f},{yB:.3f} C{xD:.3f},{yB - k:.3f} {xA:.3f},{yT + k:.3f} {xA:.3f},{yT:.3f}",
             stroke=INK, width=0.8)
    svg.path(f"M{xA:.3f},{yT + 0.3:.3f} L{xA - 2.4:.3f},{base:.3f} L{xA + 2.4:.3f},{base:.3f} Z", fill=ACCEPTOR_COLOR)

    if contact is not None:
        mid = _cubic_mid((xD, yB), (xD, yB - k), (xA, yT + k), (xA, yT))
        src = "ChIA-PET" if "chia" in contact.source.lower() else "Hi-C"
        l1 = f"{src} contact: {contact.tissue}"
        l2 = f"p = {_fmt_stat(contact.p)}, q = {_fmt_stat(contact.q)}"
        wmax = max(text_width(l1, 6.5), text_width(l2, 6.5))
        rising = xA >= xD
        if mid[0] + 8 + wmax <= width - 4:
            x, anchor, below = mid[0] + 8, "start", rising
        else:
            x, anchor, below = mid[0] - 8, "end", not rising
        y1 = mid[1] + 6 if below else mid[1] - 9
        svg.text(x, y1, l1, size=6.5, anchor=anchor)
        svg.text(x, y1 + 8, l2, size=6.5, anchor=anchor, fill=GREY_TEXT)

    y = cy_bot + thickness / 2 + 8.5 + 6
    values = make_values(row, donor, acceptor, contact, scale, "horizontal")
    if legend:
        y = draw_legend(svg, y + 10, x0=8.0)
    if caption:
        y = draw_caption(svg, values, y, width, x0=8.0, band_labels=False)
    return svg.render(width, y + 6), values


def caption_text(v: dict, band_labels: bool = True) -> List[str]:
    dg, ag = v["donor_gene"], v["acceptor_gene"]
    if v["contact_source"]:
        if v["contact_q"] is not None and v["contact_significant"] == "q<0.05":
            claim = f"the two loci share a DNA contact at q < 0.05."
        elif v["contact_p"] is not None and v["contact_p"] < 0.05:
            claim = f"the two loci share a DNA contact at p < 0.05 (q >= 0.05)."
        else:
            claim = f"no DNA contact between the two loci reaches p < 0.05; the ribbon shows the strongest tested pair."
    else:
        claim = "the contact table holds no tested contact for this pair."
    title = f"{dg} ({v['donor_band']}) trans-splices to {ag} ({v['acceptor_band']}); {claim}"
    shows = (f"GRCh38 ideograms (UCSC cytoBand, Giemsa stain) at one scale ({v['scale_bp_per_px'] / 1e6:.3f} Mb per unit). "
             f"Ticks mark the junction coordinates, donor {v['donor_chrom']}:{v['donor_pos']:,} ({v['donor_strand']}) and "
             f"acceptor {v['acceptor_chrom']}:{v['acceptor_pos']:,} ({v['acceptor_strand']}); the outlined bands hold them. "
             f"The arrow runs from donor to acceptor.")
    if v["contact_source"]:
        stats = (f"The ribbon joins the contact bins {v['contact_donor_bin']} and {v['contact_acceptor_bin']} of the "
                 f"{v['contact_rule']} association: {v['contact_source']}, dataset {v['contact_dataset']}, "
                 f"{v['contact_tissue']}; z = {_fmt_stat(v['contact_z'])}, p = {_fmt_stat(v['contact_p'])}, "
                 f"q = {_fmt_stat(v['contact_q'])}. The best-by-q row is shown when q < 0.05, otherwise the best-by-p row.")
    else:
        stats = "No DNA contact row exists for this junction."
    support = []
    if v["kinnex_junction_molecules"] is not None:
        support.append(f"{v['kinnex_junction_molecules']} KINNEX junction molecules in {v['kinnex_samples_positive']} "
                       f"samples ({v['kinnex_patients_positive']} patients)")
    if v["ogawa_fragments_362"] is not None:
        support.append(f"{v['ogawa_fragments_362']} confirmed Ogawa bulk RNA-seq fragments in "
                       f"{v['ogawa_samples_positive_of_362']} of 362 samples")
    scope = ("Junction support: " + "; ".join(support) + ". " if support else "") + \
            ("Values come from the per-junction table; band labels give the band number within each arm."
             if band_labels else "Values come from the per-junction table.")
    return [title, shows, stats, scope]


# ------------------------------------------------------------------------------------------------ CLI
def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--junctions", required=True, help="per-junction table (TSV), one row per junction")
    parser.add_argument("--chimera", action="append", required=True,
                        help="chimeric_rna value to draw, e.g. 'FOSB>SF1'; repeat for several")
    parser.add_argument("--junction", action="append", default=None,
                        help="restrict to these junction values ('chrA:pos:strand > chrB:pos:strand'), for a gene pair with "
                             "several junctions; repeat for several")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--cytobands", default=DEFAULT_CYTOBANDS, help="UCSC cytoBand table (GRCh38)")
    parser.add_argument("--layout", choices=["horizontal", "vertical"], default="horizontal",
                        help="horizontal: compact panel for placement under an ISV (default); vertical: banded "
                             "ideograms side by side with every band labelled")
    parser.add_argument("--width", type=float, default=720.0,
                        help="horizontal layout: drawing width in units (1 unit = 1 pt in Illustrator); 720 matches "
                             "the 10-inch minimum width of isoform_structure_view figures")
    parser.add_argument("--height", type=float, default=520.0,
                        help="vertical layout: drawn length of the longer chromosome (units)")
    parser.add_argument("--no-title", action="store_true", help="horizontal layout: omit the title lines")
    parser.add_argument("--no-legend", action="store_true", help="horizontal layout: omit the stain and mark key")
    parser.add_argument("--no-caption", action="store_true", help="omit the written legend under the drawing")


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s.replace(">", "_to_"))


def run(args) -> int:
    t = pd.read_csv(args.junctions, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED_COLUMNS if c not in t.columns]
    if missing:
        raise ValueError(f"junction table lacks columns {missing}")
    bands = load_cytobands(args.cytobands)
    os.makedirs(args.outdir, exist_ok=True)
    rows_out = []
    for chim in args.chimera:
        sel = t[t.chimeric_rna == chim]
        if args.junction:
            sel = sel[sel.junction.isin(args.junction)]
        if sel.empty:
            raise ValueError(f"{chim} is not in {args.junctions}" + (" at the given --junction values" if args.junction else ""))
        for k, (_, row) in enumerate(sel.iterrows()):
            suffix = "" if len(sel) == 1 else f"_{k + 1}"
            if args.layout == "vertical":
                svg_text, values = render_chimera(row, bands, height=args.height, caption=not args.no_caption)
            else:
                svg_text, values = render_chimera_horizontal(row, bands, width=args.width, title=not args.no_title,
                                                             legend=not args.no_legend, caption=not args.no_caption)
            out = os.path.join(args.outdir, f"{safe_name(chim)}{suffix}_ideogram_{args.layout}.svg")
            if os.path.exists(out):
                raise FileExistsError(f"{out} exists; choose a new --outdir")
            with open(out, "w") as fh:
                fh.write(svg_text)
            values["svg"] = out
            rows_out.append(values)
            print(f"{chim}: {out}")
    vt = os.path.join(args.outdir, "IDEOGRAM_VALUES.tsv")
    pd.DataFrame(rows_out).to_csv(vt, sep="\t", index=False)
    params = {"module": "altanalyze3.components.visualization.chromosome_ideogram",
              "junctions": os.path.abspath(args.junctions), "junctions_md5": md5sum(args.junctions),
              "cytobands": os.path.abspath(args.cytobands), "cytobands_md5": md5sum(args.cytobands),
              "chimera": args.chimera, "junction": args.junction, "layout": args.layout, "width": args.width, "height": args.height,
              "title": not args.no_title, "legend": not args.no_legend, "caption": not args.no_caption,
              "contact_rule": "best by q when overall_bestq_significant == 'q<0.05', else best by p",
              "outputs": [r["svg"] for r in rows_out] + [vt]}
    os.makedirs(os.path.join(args.outdir, "logs"), exist_ok=True)
    with open(os.path.join(args.outdir, "logs", "chromosome_ideogram_params.json"), "w") as fh:
        json.dump(params, fh, indent=1)
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    add_arguments(p)
    return run(p.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
