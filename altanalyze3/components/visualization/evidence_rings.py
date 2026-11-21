#!/usr/bin/env python3
"""Evidence rings: three concentric rings per panel that show how a set of direct calls grows with imputed
calls, and which evidence supports each imputed call.

    Ring 1 (inner)   direct calls per item (for example cells in which a variant was read).
    Ring 2 (middle)  direct plus imputed calls per item; a pale band joins the two arcs of an item.
    Ring 3 (outer)   one thin shell per evidence type over the imputed calls only. A shell is coloured where
                     that evidence supports the call. Colour in two shells at one angle marks the same calls
                     supported by both types (the calls are laid out in a reflected order so such runs align).

Input is one table with one row per call, or one row per group of identical calls with a ``count`` column:

    panel     panel name (one ring set per panel)
    item      item name (one colour per item)
    direct    1/True for a direct call, 0/False for an imputed call
    evidence  evidence codes of an imputed call joined by "+", for example "PV+CNV"; empty for a direct call
    count     optional number of calls in the row (default 1)

Optional per-item table (``items``): panel, item, and any of label_suffix (appended to the item name),
label_detail (start of the second label line), color (#RRGGBB).

Ring modes
    ring1="aligned"   ring-1 arc at the clockwise start of its ring-2 arc; one degree = one count in all rings
    ring1="centered"  ring-1 arc centred under its ring-2 arc; one degree = one count in all rings
    ring1="full"      ring 1 is a full circle on its own scale (share of the direct calls); band = share
    band="sweep" | "curved" | "straight"   edge shape of the band between ring 1 and ring 2

The PDF is vector only (figure_style.apply_rcparams / strip_masks). Colours of the items come from
palettes.categorical_colors unless the item table gives them.

Command line:
    python3 -m altanalyze3.components.visualization.evidence_rings --calls calls.tsv --out rings.pdf \
        --shells shells.tsv --ring1 centered --top-labels 18
"""

from __future__ import annotations

import argparse
import math
import textwrap
from typing import Dict, List, Optional, Sequence

import matplotlib
import numpy as np
import pandas as pd

from .figure_style import apply_rcparams, strip_masks
from .palettes import categorical_colors

DEFAULT_SHELL_COLORS = ["#d1352b", "#2f6fb7", "#2e9a62", "#e08a1e", "#7b4ea3", "#1f9aa6", "#8c564b", "#5f6b7a"]
R2I, R2O = 0.76, 0.94                      # ring 2
R3I, SW, SG = 0.99, 0.052, 0.012           # ring 3: inner radius, shell width, white gap between shells
STEP = 0.24                                # vertical distance between two labels, ring units
UNIT_IN = 0.79                             # inches per ring unit
SMALL_FRACTION = 0.03                      # items below this share of a panel form one run at 180 degrees


def _lighten(hexc: str, f: float) -> str:
    r, g, b = (int(hexc[i:i + 2], 16) for i in (1, 3, 5))
    return "#%02x%02x%02x" % tuple(int(c + (255 - c) * f) for c in (r, g, b))


def _xy(r: float, a: float):
    return (r * math.cos(math.radians(a)), r * math.sin(math.radians(a)))


def reflected_order(combos: Sequence[tuple], n_codes: int, k: int = 0, rev: bool = True) -> List[tuple]:
    """Order membership tuples so the calls of shell 1 form one run and the calls of shell 2 form one run that
    spans the edge of the shell-1 run; coincident calls of two evidence types then sit radially above each other."""
    combos = list(combos)
    if k == n_codes or len(combos) <= 1:
        return combos
    ones = [c for c in combos if c[k]]
    zeros = [c for c in combos if not c[k]]
    if rev:
        return reflected_order(ones, n_codes, k + 1, False) + reflected_order(zeros, n_codes, k + 1, True)
    return reflected_order(zeros, n_codes, k + 1, False) + reflected_order(ones, n_codes, k + 1, True)


def _place_labels(labels, r3o: float, ymax: float):
    """Side (+1 right, -1 left) and height of each label; labels on one side keep a distance of STEP."""
    out = {}
    sides = {i: (1 if math.cos(math.radians(a)) >= 0 else -1) for i, a in enumerate(labels)}
    maxn = int((2 * ymax - 0.3) / STEP)
    for sd in (1, -1):
        while sum(1 for v in sides.values() if v == sd) > maxn:
            j = min((i for i, v in sides.items() if v == sd), key=lambda i: abs(math.cos(math.radians(labels[i]))))
            sides[j] = -sd
    for side in (1, -1):
        idx = sorted((i for i in sides if sides[i] == side), key=lambda i: math.sin(math.radians(labels[i])))
        if not idx:
            continue
        ys = [(r3o + 0.35) * math.sin(math.radians(labels[i])) for i in idx]
        for k in range(1, len(idx)):
            ys[k] = max(ys[k], ys[k - 1] + STEP)
        ys[-1] = min(ys[-1], ymax - 0.15)
        for k in range(len(idx) - 2, -1, -1):
            ys[k] = min(ys[k], ys[k + 1] - STEP)
        if ys[0] < -ymax:
            raise ValueError("too many labels on one side for the panel height; lower top_labels or raise ymax")
        for i, y in zip(idx, ys):
            out[i] = (side, y)
    return out


def prepare_calls(calls: pd.DataFrame) -> pd.DataFrame:
    """Validate the call table and return it aggregated to panel, item, direct, evidence, count."""
    need = {"panel", "item", "direct"}
    missing = need - set(calls.columns)
    if missing:
        raise ValueError("call table lacks the columns: %s" % sorted(missing))
    c = calls.copy()
    if c["direct"].dtype != bool:
        c["direct"] = c["direct"].astype(str).str.strip().str.lower().isin(["1", "true", "t", "yes"])
    c["evidence"] = c["evidence"].fillna("").astype(str) if "evidence" in c.columns else ""
    c["count"] = c["count"].astype(int) if "count" in c.columns else 1
    if (c["count"] < 0).any():
        raise ValueError("negative count in the call table")
    bad = c[(~c["direct"]) & (c["evidence"] == "")]
    if len(bad):
        raise ValueError("%d imputed rows carry no evidence code" % len(bad))
    c.loc[c["direct"], "evidence"] = ""
    return c.groupby(["panel", "item", "direct", "evidence"], sort=False)["count"].sum().reset_index()


def _draw_panel(ax, items: pd.DataFrame, calls: pd.DataFrame, shells, colour: Dict[str, str], *, xmax: float,
                ymax: float, ring1: str, band: str, top_labels: Optional[int], unit_word: str,
                center_words, label_size: float):
    from matplotlib.patches import Polygon, Rectangle, Wedge
    import matplotlib.pyplot as plt

    codes = [s[0] for s in shells]
    full = ring1 == "full"
    center = ring1 in ("centered", "full")
    r1i, r1o = (0.36, 0.52) if full else (0.40, 0.58)
    r3o = R3I + len(shells) * SW + (len(shells) - 1) * SG
    tot = int(items["total"].sum())
    if tot == 0:
        raise ValueError("a panel holds no calls")
    deg = 360.0 / tot
    ax.set_xlim(-xmax, xmax); ax.set_ylim(-ymax, ymax); ax.set_aspect("equal"); ax.axis("off")

    def wedge(r0, r1_, lo, hi, col):
        thin = (hi - lo) < 0.8
        ax.add_patch(Wedge((0, 0), r1_, lo, hi, width=r1_ - r0, facecolor=col,
                           edgecolor=col if thin else "#ffffff", linewidth=0.25))

    rows = list(items.itertuples())
    small = sum(deg * r.total for r in rows if r.total / tot < SMALL_FRACTION)
    start = 180.0 + small / 2.0 if len(rows) > 3 else 90.0
    order = [r for r in rows if r.total / tot < SMALL_FRACTION][::-1] + [r for r in rows if r.total / tot >= SMALL_FRACTION]
    ring1_arc = {}
    n_direct = float(sum(r.direct for r in order))
    if full and n_direct > 0:
        # one rotation for the whole ring: smallest direct-weighted squared angle between the two arc middles
        d1 = 360.0 / n_direct
        m2, m1, w1 = [], [], []
        a_, c_ = start, 0.0
        for r in order:
            t_, w_ = deg * r.total, d1 * r.direct
            m2.append(a_ - t_ / 2.0); m1.append(-(c_ + w_ / 2.0)); w1.append(r.direct)
            a_ -= t_; c_ += w_
        m2, m1, w1 = np.array(m2), np.array(m1), np.array(w1, dtype=float)
        offs = np.arange(0.0, 360.0, 0.5)
        cost = [float((w1 * (((m1 + o - m2 + 180.0) % 360.0) - 180.0) ** 2).sum()) for o in offs]
        off1 = float(offs[int(np.argmin(cost))])
        c_ = 0.0
        for r, mm in zip(order, m2):
            w_ = d1 * r.direct
            hi = off1 - c_
            k_ = round((mm - (hi - w_ / 2.0)) / 360.0)        # same turn as the ring-2 arc: a band never wraps
            ring1_arc[r.item] = (hi + 360.0 * k_ - w_, hi + 360.0 * k_)
            c_ += w_
    ang = start
    labels = []
    shell_calls = {c: 0 for c in codes}
    by_item = {k: g for k, g in calls[~calls["direct"]].groupby("item", sort=False)}
    for r in order:
        col = colour[r.item]
        t_arc, d_arc = deg * r.total, deg * r.direct
        a1 = ang
        flank = (t_arc - d_arc) / 2.0 if center else 0.0      # imputed calls drawn before the direct arc
        h1 = a1 - flank; l1 = h1 - d_arc
        if full and r.item in ring1_arc:
            l1, h1 = ring1_arc[r.item]
        wedge(r1i, r1o, l1, h1, col)
        wedge(R2I, R2O, a1 - t_arc, a1, col)
        ts = np.linspace(0, 1, 24)
        sm = ts if band == "sweep" else 3 * ts ** 2 - 2 * ts ** 3
        lead = [_xy(r1o + (R2I - r1o) * t, h1 + (a1 - h1) * u) for t, u in zip(ts, sm)]
        trail = [_xy(r1o + (R2I - r1o) * t, l1 + ((a1 - t_arc) - l1) * u) for t, u in zip(ts, sm)]
        if band == "straight":
            lead = [_xy(r1o, h1), _xy(R2I, a1)]; trail = [_xy(r1o, l1), _xy(R2I, a1 - t_arc)]
        inner = [_xy(r1o, a) for a in np.linspace(h1, l1, max(2, int((h1 - l1) / 2) + 2))]
        outer = [_xy(R2I, a) for a in np.linspace(a1 - t_arc, a1, max(2, int(t_arc / 2) + 2))]
        ax.add_patch(Polygon(inner + trail[1:-1] + outer + lead[::-1][1:-1], closed=True,
                             facecolor=_lighten(col, 0.84), edgecolor="none"))
        for edge, flat in ((lead, flank < 1e-9 and not full), (trail, t_arc - d_arc - flank < 1e-9 and not full)):
            if flat or len(edge) == 2:
                ax.plot(*zip(edge[0], edge[-1]), color=col, lw=0.35, solid_capstyle="butt")
            else:
                ax.plot(*zip(*edge), color=col, lw=0.35, solid_capstyle="butt")
        # ring 3: imputed calls of this item in the reflected order of their evidence combination
        inf = by_item.get(r.item)
        memb: Dict[tuple, int] = {}
        if inf is not None:
            for e, n in zip(inf["evidence"], inf["count"]):
                key = tuple(c in e.split("+") for c in codes)
                memb[key] = memb.get(key, 0) + int(n)
        if sum(memb.values()) != r.total - r.direct:
            raise AssertionError("imputed calls of %s do not sum to total minus direct" % r.item)
        off = 0.0
        runs = {c: [] for c in codes}
        for m in reflected_order(list(memb), len(codes)):
            w = deg * memb[m]
            for k, c in enumerate(codes):
                if m[k]:
                    shell_calls[c] += memb[m]
                    if runs[c] and abs(runs[c][-1][1] - off) < 1e-9:
                        runs[c][-1][1] = off + w
                    else:
                        runs[c].append([off, off + w])
            off += w
        assert abs(off - (t_arc - d_arc)) < 1e-6
        for k, (c, scol, _) in enumerate(shells):
            r0 = R3I + k * (SW + SG)
            for o0, o1 in runs[c]:
                pieces = []
                if o0 < flank:
                    pieces.append((a1 - min(o1, flank), a1 - o0))
                if o1 > flank:
                    pieces.append((a1 - d_arc - o1, a1 - d_arc - max(o0, flank)))
                for lo, hi in pieces:
                    ax.add_patch(Wedge((0, 0), r0 + SW, lo, hi, width=SW, facecolor=scol, edgecolor="none"))
        labels.append((a1 - t_arc / 2.0, r, r.total))
        ang = a1 - t_arc
    assert abs((ang + 360.0) - start) < 1e-6
    for k in range(len(shells) + 1):          # hairline circles keep an empty shell visible
        ax.add_patch(plt.Circle((0, 0), R3I + k * (SW + SG) - SG / 2.0, facecolor="none", edgecolor="#d5d8dc", linewidth=0.2))
    if top_labels is not None and len(labels) > top_labels:
        keep = sorted(range(len(labels)), key=lambda i: -labels[i][2])[:top_labels]
        labels = [labels[i] for i in sorted(keep)]
    pos = _place_labels([a for a, _, _ in labels], r3o, ymax)
    for i, (a, r, _) in enumerate(labels):
        side, y = pos[i]
        x0, y0 = _xy(R2O if center else r3o + 0.02, a)
        x1 = side * (r3o + 0.42)
        lk = dict(color=colour[r.item], lw=1.2, alpha=0.5) if center else dict(color="#9aa0a6", lw=0.3)
        ax.plot([x0, x1], [y0, y], solid_capstyle="butt", **lk)
        ax.plot([x1, side * (r3o + 0.52)], [y, y], solid_capstyle="butt", **lk)
        n = ("%s %s" % (format(r.direct, ","), unit_word)) if r.total == r.direct else (
            "%s → %s %s" % (format(r.direct, ","), format(r.total, ","), unit_word))
        line2 = (r.label_detail + ", " if r.label_detail else "") + n
        ax.add_patch(Rectangle((side * (r3o + 0.56) - (0.0 if side == 1 else 0.07), y + 0.025), 0.07, 0.07,
                               facecolor=colour[r.item], edgecolor="none"))
        ax.text(side * (r3o + 0.67), y, "%s%s\n%s" % (r.item, r.label_suffix, line2), ha="left" if side == 1 else "right",
                va="center", fontsize=label_size, color="#202124", linespacing=1.08)
    nd = int(items["direct"].sum())
    for k, (line, fs) in enumerate([(format(nd, ","), 6.2), (center_words[0], 4.8), (format(tot, ","), 6.2), (center_words[1], 4.8)]):
        ax.text(0, 0.21 - 0.14 * k, line, ha="center", va="center", fontsize=fs, color="#202124")
    return shell_calls, nd, tot


def plot_evidence_rings(calls: pd.DataFrame, out: str, *, shells: Optional[Sequence[Sequence[str]]] = None,
                        items: Optional[pd.DataFrame] = None, panels: Optional[Sequence[str]] = None,
                        panel_titles: Optional[Dict[str, str]] = None, ring1: str = "centered", band: str = "sweep",
                        top_labels: Optional[int] = None, unit_word: str = "cells",
                        center_words: Sequence[str] = ("read-called", "with imputed"), legend: Optional[str] = None,
                        xmax: Optional[Sequence[float]] = None, ymax: Optional[float] = None,
                        figsize: Optional[Sequence[float]] = None, label_size: float = 4.8) -> dict:
    """Draw the figure and write ``out`` (PDF). Returns the per-panel counts drawn.

    shells   sequence of (code, colour, name), inside to outside; default: every code in ``calls`` in order
             of first appearance, colours from DEFAULT_SHELL_COLORS.
    items    optional per-item table (see the module text). Items are drawn in decreasing order of total.
    panels   panel order; default: order of first appearance.
    xmax, ymax, figsize   layout overrides (ring units, ring units, inches); computed when omitted.
    """
    if ring1 not in ("aligned", "centered", "full"):
        raise ValueError("ring1 must be aligned, centered or full")
    if band not in ("sweep", "curved", "straight"):
        raise ValueError("band must be sweep, curved or straight")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    c = prepare_calls(calls)
    seen: List[str] = []
    for e in c.loc[~c["direct"], "evidence"].unique():
        for code in e.split("+"):
            if code not in seen:
                seen.append(code)
    if shells is None:
        if len(seen) > len(DEFAULT_SHELL_COLORS):
            raise ValueError("more evidence codes than default shell colours; pass shells")
        shells = [(code, DEFAULT_SHELL_COLORS[i], code) for i, code in enumerate(seen)]
    shells = [tuple(s) for s in shells]
    codes = [s[0] for s in shells]
    unknown = set(seen) - set(codes)
    if unknown:
        raise ValueError("evidence code without a shell: %s" % sorted(unknown))
    panels = list(panels) if panels is not None else list(pd.unique(c["panel"]))
    panel_titles = panel_titles or {}

    summ = c.assign(direct_n=c["count"].where(c["direct"], 0)).groupby(["panel", "item"], sort=False).agg(
        direct=("direct_n", "sum"), total=("count", "sum")).astype(int).reset_index()
    for col in ("label_suffix", "label_detail", "color"):
        summ[col] = ""
    if items is not None:
        it = items.copy()
        for col in ("label_suffix", "label_detail", "color"):
            if col in it.columns:
                m = dict(zip(zip(it["panel"], it["item"]), it[col].fillna("").astype(str)))
                summ[col] = [m.get(k, "") for k in zip(summ["panel"], summ["item"])]
    summ = summ[summ["panel"].isin(panels)]
    if summ.duplicated("item").any():
        raise ValueError("an item name occurs in two panels")
    ranked = summ.sort_values(["panel", "direct"], ascending=[True, False])["item"].tolist()
    colour = dict(zip(ranked, categorical_colors(len(ranked))))
    colour.update({r.item: r.color for r in summ.itertuples() if r.color})

    apply_rcparams(6)
    r3o = R3I + len(shells) * SW + (len(shells) - 1) * SG
    per_panel = {p: summ[summ["panel"] == p].sort_values("total", ascending=False) for p in panels}
    if xmax is None:
        xmax = []
        for p in panels:
            d = per_panel[p]
            shown = d if top_labels is None else d.head(top_labels)
            chars = max([max(len(r.item) + len(r.label_suffix), len(r.label_detail) + 26) for r in shown.itertuples()] + [12])
            xmax.append(r3o + 0.7 + 0.053 * chars * label_size / 4.8)
    if ymax is None:
        most = max(min(len(per_panel[p]), top_labels or 10 ** 9) for p in panels)
        ymax = max(r3o + 0.5, (math.ceil(most / 2.0) + 2) * STEP / 2.0 + 0.4)
    leg_lines = sum(len(textwrap.wrap(par, 400)) or 1 for par in legend.split("\n")) if legend else 0
    foot_in = 0.55 + 0.105 * leg_lines
    if figsize is None:
        figsize = (2 * sum(xmax) * UNIT_IN / 0.99, 2 * ymax * UNIT_IN + foot_in + 0.3)
    fig = plt.figure(figsize=tuple(figsize))
    bottom = foot_in / figsize[1]
    gs = fig.add_gridspec(1, len(panels), left=0.005, right=0.995, top=1 - 0.24 / figsize[1], bottom=bottom,
                          wspace=0.0, width_ratios=list(xmax))
    stats = {}
    for i, p in enumerate(panels):
        ax = fig.add_subplot(gs[0, i])
        d = per_panel[p]
        stats[p] = _draw_panel(ax, d, c[c["panel"] == p], shells, colour, xmax=xmax[i], ymax=ymax, ring1=ring1, band=band,
                               top_labels=top_labels, unit_word=unit_word, center_words=center_words, label_size=label_size)
        ax.set_title("%s   %s (%d)" % ("abcdefghijklmnopqrstuvwxyz"[i], panel_titles.get(p, p), len(d)), loc="left",
                     fontsize=7, fontweight="bold", pad=2)
    # key: the shells of ring 3, inside to outside
    kax = fig.add_axes([0.02, bottom - 0.42 / figsize[1], 0.96, 0.34 / figsize[1]])
    kax.set_xlim(0, 100); kax.set_ylim(0, 2); kax.axis("off")
    kax.text(0, 1.5, "Ring 1 (inner): %s.   Ring 2 (middle): %s; a band joins the two arcs of an item.   Ring 3 (outer): evidence of the "
             "imputed calls, one shell per type, white = not supported by that type." % (center_words[0], center_words[1]),
             ha="left", va="center", fontsize=5.8, color="#202124")
    kax.text(0, 0.5, "Ring 3 shells, inside to outside:", ha="left", va="center", fontsize=5.8, color="#202124")
    x = 12.0 * 17.0 / figsize[0]
    per_char = 0.262 * 17.0 / figsize[0]
    for i, (code, scol, name) in enumerate(shells):
        kax.add_patch(Rectangle((x, 0.25), 0.7, 0.5, facecolor=scol, edgecolor="none"))
        t = "%d  %s" % (i + 1, code) + ("  %s" % name if name and name != code else "")
        kax.text(x + 0.95, 0.5, t, ha="left", va="center", fontsize=5.6, color="#202124")
        x += 2.2 + per_char * len(t)
    if legend:
        fig.text(0.02, bottom - 0.47 / figsize[1], "\n".join(textwrap.fill(par, 400) for par in legend.split("\n")),
                 ha="left", va="top", fontsize=5.6, color="#202124", linespacing=1.25)
    strip_masks(fig)
    fig.savefig(out, format="pdf")
    plt.close(fig)
    return {p: dict(direct=stats[p][1], total=stats[p][2], shells=stats[p][0]) for p in panels}


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--calls", required=True, help="TSV: panel, item, direct, evidence[, count]")
    parser.add_argument("--out", required=True, help="output PDF")
    parser.add_argument("--shells", help="TSV: code, color, name; ring-3 shells inside to outside")
    parser.add_argument("--items", help="TSV: panel, item[, label_suffix, label_detail, color]")
    parser.add_argument("--panels", help="comma-separated panel order")
    parser.add_argument("--panel-titles", help="comma-separated panel=title pairs")
    parser.add_argument("--ring1", default="centered", choices=["aligned", "centered", "full"])
    parser.add_argument("--band", default="sweep", choices=["sweep", "curved", "straight"])
    parser.add_argument("--top-labels", type=int, default=None, help="label only the N largest items of each panel")
    parser.add_argument("--unit-word", default="cells")
    parser.add_argument("--center-words", default="read-called,with imputed")
    parser.add_argument("--legend", help="text file printed under the figure")
    parser.add_argument("--counts-out", help="TSV of the counts drawn per panel")


def run(args) -> dict:
    calls = pd.read_csv(args.calls, sep="\t")
    shells = None
    if args.shells:
        s = pd.read_csv(args.shells, sep="\t")
        shells = list(zip(s["code"].astype(str), s["color"].astype(str), s["name"].astype(str)))
    items = pd.read_csv(args.items, sep="\t") if args.items else None
    titles = dict(x.split("=", 1) for x in args.panel_titles.split(",")) if args.panel_titles else None
    stats = plot_evidence_rings(
        calls, args.out, shells=shells, items=items, panels=args.panels.split(",") if args.panels else None,
        panel_titles=titles, ring1=args.ring1, band=args.band, top_labels=args.top_labels, unit_word=args.unit_word,
        center_words=tuple(args.center_words.split(",", 1)), legend=open(args.legend).read().rstrip("\n") if args.legend else None)
    if args.counts_out:
        pd.DataFrame([dict(panel=p, direct=v["direct"], total=v["total"], **{"shell_" + k: n for k, n in v["shells"].items()})
                      for p, v in stats.items()]).to_csv(args.counts_out, sep="\t", index=False)
    print("wrote", args.out)
    return stats


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    add_arguments(parser)
    return run(parser.parse_args(argv))


if __name__ == "__main__":
    main()
