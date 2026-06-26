#!/usr/bin/env python3
"""Backend (headless) ISV-web read-level molecule pileup -> vector PDF. No web server, no browser.

Part of the isv_web package: drives the SAME data path the live viewer uses (``data_api``) against the
PRE-BUILT per-sample indexes, then renders the molecule pileup. Nothing is rebuilt.

  data  : data_api.RunContext(run_dir, metadata, gene_model, gene_symbol).load_all()
          -> opens <run_dir>/_isv_web_cache/mol_index/<lib>.reads.db (molecules)
             + <lib>.mol2final.db (molecule -> final collapsed isoform)
  query : data_api.query_reads(ctx, gene, cell_states, conditions, panel_by="cell_type")
          -> _query_reads_index fast path (~ms): one panel per cell state; each molecule already carries
             exon_segments, final_isoform_id, cluster_index (grouped by FINAL collapsed isoform).
  render: matplotlib reproduction of static/app.js drawReads() -> single-page vector PDF.

Run (module form):
  python -m altanalyze3.components.visualization.isv_web.isv_web_remote \\
      --run_dir /path/to/dataset \\
      --gene TNFAIP8 --cell_states HSC-1,cDC2-2 \\
      --out ./ISV_TNFAIP8_HSC-1_cDC2-2_reads.pdf
"""
import os
import sys
import argparse
from collections import OrderedDict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# ---- drawReads() layout constants (static/app.js R = {...}) ----
LEFT, RIGHT, READ_PX, ROW_MIN, ROW_CAP = 152, 24, 1.6, 0.9, 60
MAXPANEL, PANEL_GAP, PAD, HDR, RULER_H, REF_H, MAX_ROWS, RMARGIN = 340, 14, 8, 22, 30, 26, 4000, 16
KNOWN, NOVEL, REF_FILL = "#1d5fa8", "#b0306b", "#3d4a5c"
CLUSTER_COLORS = [
    "#2b6cb0", "#c0531f", "#2f8a4e", "#7a4fb5", "#c2306b", "#0e8a9c", "#b8902a",
    "#4a5bbf", "#cf5aa0", "#2f9b8a", "#9a5b2e", "#56657a", "#1f7a5a", "#a23c5e",
    "#3d7ec2", "#8a7a1f", "#6d4ca0", "#1d8f74", "#bf6a2a", "#8a2f52", "#3a6ea5",
    "#7b9b3a", "#a64d79", "#2d7d8a",
]


TWO_ISOFORM_COLORS = ["#E9B940", "#4884BC"]   # default fills when a gene shows EXACTLY two isoforms (gold, blue)


def cluster_color(i):
    n = len(CLUSTER_COLORS)
    return CLUSTER_COLORS[(int(i) % n + n) % n]


def darken(hexcol, f=0.78):
    """Edge shade derived from the exon fill: same hue, multiplicatively darkened by f (so the edge stays
    a near-relative of the fill, not a fixed dark colour)."""
    h = str(hexcol).lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return "#%02x%02x%02x" % (int(r * f), int(g * f), int(b * f))


def color_for(m):
    ci = m.get("color_index")
    return cluster_color(ci if ci is not None else m.get("cluster_index", 0))


def truncate(s, n):
    s = str(s)
    return s if len(s) <= n else s[:n - 1] + "…"


def fmt(v):
    return f"{int(round(v)):,}"


def merge_segments(segs):
    """Port of static/app.js mergeSegments(): exon (E*) blocks merged by label/overlap; intron (I*) spans."""
    segs = segs or []
    es = sorted(
        [
            {
                "a": min(s["start"], s["end"]),
                "b": max(s["start"], s["end"]),
                "label": str(s.get("exon_id") or s.get("label") or "").split(".")[0],
            }
            for s in segs
            if s.get("start") is not None and s.get("end") is not None
            and str(s.get("type") or "E").upper().startswith("E")
        ],
        key=lambda s: s["a"],
    )
    out = []
    for s in es:
        last = out[-1] if out else None
        if last and s["a"] <= last["b"] + 2 and s["label"] == last["label"]:
            last["b"] = max(last["b"], s["b"])
        elif last and s["a"] <= last["b"]:
            last["b"] = max(last["b"], s["b"])
            if s["label"] and not last["label"]:
                last["label"] = s["label"]
        else:
            out.append(dict(s))
    introns = [
        {"a": min(s["start"], s["end"]), "b": max(s["start"], s["end"])}
        for s in segs
        if s.get("start") is not None and s.get("end") is not None
        and str(s.get("type") or "").upper().startswith("I")
    ]
    return {"exons": out, "introns": introns}


def extent_of(mols, gm):
    lo, hi = np.inf, -np.inf
    if gm.get("extent"):
        lo, hi = min(lo, gm["extent"][0]), max(hi, gm["extent"][1])
    for b in gm.get("blocks", []):
        lo, hi = min(lo, b["start"]), max(hi, b["end"])
    for m in mols:
        for s in m["merged"]["exons"] + m["merged"]["introns"]:
            lo, hi = min(lo, s["a"]), max(hi, s["b"])
    if not np.isfinite(lo):
        lo, hi = 0, 1
    return lo, hi


def _merge_ivs(ivs):
    """Merge a list of (a, b) intervals into sorted, disjoint intervals."""
    out = []
    for a, b in sorted(ivs):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def build_xz(exon_ivs, lo, hi, intron_scale, x0, x1):
    """Piecewise genomic->pixel map that scales INTRON spans by `intron_scale` while exons stay 1.0 --
    mirrors isoform_structure_view.build_gene_maps so introns can be shrunk to enlarge exons. With
    intron_scale == 1.0 this is identical to a plain linear map (backward compatible)."""
    segs, d, cur = [], 0.0, lo
    for a, b in exon_ivs:
        a, b = max(a, lo), min(b, hi)
        if b <= cur:
            continue
        if a > cur:                                  # intron gap before this exon -> scaled
            w = (a - cur) * intron_scale; segs.append((cur, a, d, d + w)); d += w; cur = a
        if b > cur:                                  # exon -> true scale
            w = (b - cur); segs.append((cur, b, d, d + w)); d += w; cur = b
    if hi > cur:                                     # trailing intron
        w = (hi - cur) * intron_scale; segs.append((cur, hi, d, d + w)); d += w
    total = d or 1.0

    def xz(c):
        c = min(max(c, lo), hi)
        for gs, ge, ds, de in segs:
            if gs <= c <= ge:
                f = (c - gs) / (ge - gs) if ge > gs else 0.0
                return x0 + ((ds + f * (de - ds)) / total) * (x1 - x0)
        return x1
    return xz


ISO_ROW, ISO_EXH = 12.0, 7.0          # bottom isoform-structure track: row height / coding-exon height


def render(res, out_path, width=1100, intron_scale=1.0, isoform_track=False):
    """Render query_reads() panels as the drawReads molecule pileup to a single-page vector PDF.
    intron_scale<1.0 shrinks intron display spans (enlarging exons); 1.0 = true genomic scale.
    isoform_track=True appends a PDF-only structure track (one row per plotted isoform; UCSC CDS-tall /
    UTR-thin, introns as lines) under the gene model, aligned to the same x-axis -- NOT shown in the web app."""
    panels = [p for p in res.get("panels", []) if p.get("n_molecules", 0) > 0]
    gm = res.get("gene_model") or {}
    if not panels:
        print("No molecules for this selection -> no PDF written.")
        return False

    for p in panels:
        for m in p["molecules"]:
            m["merged"] = merge_segments(m.get("exon_segments"))
        p["mol"] = sorted(
            p["molecules"],
            key=lambda m: (m.get("cluster_index", 0), -m.get("count", 0), str(m.get("isoform_id"))),
        )
        if len(p["mol"]) > MAX_ROWS:
            p["mol"] = p["mol"][:MAX_ROWS]

    all_mol = [m for p in panels for m in p["mol"]]
    lo, hi = extent_of(all_mol, gm)
    track_x0, track_x1 = LEFT, width - RIGHT - RMARGIN

    # default two-isoform palette (gold / blue): applied ONLY when the gene shows exactly two isoforms,
    # consistently across the molecule pileup and the structure track.
    def _ckey(m):
        return m.get("color_index") if m.get("color_index") is not None else m.get("cluster_index", 0)
    _ckeys = sorted({_ckey(m) for m in all_mol})
    _two = {k: TWO_ISOFORM_COLORS[i] for i, k in enumerate(_ckeys)} if len(_ckeys) == 2 else {}

    def colof(m):
        return _two.get(_ckey(m)) or color_for(m)

    # unique plotted isoforms (one structure row each) for the optional bottom structure track
    iso_rows = []
    if isoform_track:
        seen = set()
        for m in sorted(all_mol, key=lambda mm: (mm.get("cluster_index", 0),
                                                 str(mm.get("final_isoform_id") or mm.get("isoform_id")))):
            key = m.get("final_isoform_id") or m.get("isoform_id")
            if key in seen or not (m.get("merged") or {}).get("exons"):
                continue
            seen.add(key); iso_rows.append(m)

    # intron-compressed x-map: exon intervals come from the gene-model blocks + the molecules' merged
    # exons; everything NOT covered by an exon is treated as intron and scaled by intron_scale.
    exon_ivs = []
    for b in gm.get("blocks", []):
        if b.get("start") is not None and b.get("end") is not None:
            exon_ivs.append((min(b["start"], b["end"]), max(b["start"], b["end"])))
    for m in all_mol:
        for e in m["merged"]["exons"]:
            exon_ivs.append((min(e["a"], e["b"]), max(e["a"], e["b"])))
    xz = build_xz(_merge_ivs(exon_ivs), lo, hi, max(0.01, float(intron_scale)), track_x0, track_x1)

    for p in panels:
        raw = [min(ROW_CAP, max(ROW_MIN, m.get("count", 0) * READ_PX)) for m in p["mol"]]
        raw_sum = sum(raw) or 1
        n_clusters = len({m.get("cluster_index", 0) for m in p["mol"]})
        scale = min(MAXPANEL, max(raw_sum, n_clusters * 13)) / raw_sum
        h, rows = 0.0, []
        for i, m in enumerate(p["mol"]):
            hh = max(0.55, raw[i] * scale)
            rows.append({"m": m, "dy": h, "h": hh})
            h += hh
        p["_rows"] = rows
        p["_bodyH"] = max(18, h)

    # Auto-fit (no knob): stacking every cell-state panel on one page makes it taller than a screen, so
    # viewers cut off the lower panels. If the figure would exceed one screenful, shrink ALL panel bodies
    # + their rows uniformly so the whole thing fits, keeping the shared genomic x-axis. Only shrinks;
    # figures that already fit are untouched.
    _FIT_H = 900.0
    _ov = 10 + 2 + RULER_H + REF_H + 22 + len(panels) * (HDR + PAD + PANEL_GAP)
    _budget = _FIT_H - _ov
    _body_total = sum(p["_bodyH"] for p in panels)
    if _body_total > _budget > 0:
        _f = _budget / _body_total
        for p in panels:
            p["_bodyH"] *= _f
            for r in p["_rows"]:
                r["dy"] *= _f
                r["h"] *= _f

    y = 10.0
    for p in panels:
        p["_top"] = y
        p["_y"] = y + HDR
        y += HDR + p["_bodyH"] + PAD + PANEL_GAP
    axis_top = y + 2
    height = axis_top + RULER_H + REF_H + 22
    ruler_base = axis_top + RULER_H - 4
    ref_top = axis_top + RULER_H
    struct_top = height + 6                              # bottom isoform-structure track starts here
    if iso_rows:
        height = struct_top + 16 + len(iso_rows) * ISO_ROW + 8   # +header +rows

    fig = plt.figure(figsize=(width / 72.0, height / 72.0))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.axis("off")

    def rect(x, yy, w, hh, fc, ec="none", lw=0, alpha=1.0):
        ax.add_patch(Rectangle((x, yy), w, hh, facecolor=fc, edgecolor=ec, linewidth=lw, alpha=alpha))

    def hline(x0, x1, yy, color, lw=1.0, alpha=1.0):
        ax.plot([x0, x1], [yy, yy], color=color, linewidth=lw, alpha=alpha, solid_capstyle="butt")

    for p in panels:
        top = p["_top"]
        box_h = HDR + p["_bodyH"] + PAD
        rect(6, top, width - 12, box_h, "#ffffff", ec="#d4d9e0", lw=0.8)
        rect(6, top, width - 12, HDR, "#eef1f5")
        ax.text(16, top + 15, f"{p['condition']}", fontsize=11, weight="bold", va="center", color="#1a2230")
        ax.text(16 + 8.5 * len(str(p["condition"])) + 10, top + 15,
                f"· {p['n_molecules']:,} reads", fontsize=9.5, va="center", color="#5a6675")
        cl_extent = {}
        for r in p["_rows"]:
            ci = r["m"].get("cluster_index", 0)
            if ci not in cl_extent:
                cl_extent[ci] = [r["dy"], r["dy"] + r["h"], r["m"]]
            else:
                cl_extent[ci][0] = min(cl_extent[ci][0], r["dy"])
                cl_extent[ci][1] = max(cl_extent[ci][1], r["dy"] + r["h"])
        for ci, (yt, yb, m) in cl_extent.items():
            ax.text(LEFT - 8, p["_y"] + (yt + yb) / 2,
                    truncate(m.get("final_isoform_id") or m.get("isoform_id"), 20),
                    fontsize=8, ha="right", va="center", color=(KNOWN if m.get("known") else NOVEL))
        for r in p["_rows"]:
            m = r["m"]
            y_top = p["_y"] + r["dy"]
            color = colof(m)
            merged = m["merged"]
            exs = merged["exons"]
            eh = max(0.7, r["h"] - 0.8 if r["h"] > 3 else r["h"])
            ey = y_top + (r["h"] - eh) / 2
            if exs:
                a = xz(min(e["a"] for e in exs))
                b = xz(max(e["b"] for e in exs))
                hline(a, b, y_top + r["h"] / 2, color, lw=0.7, alpha=0.55)
            for s in merged["introns"]:
                xa, xb = xz(s["a"]), xz(s["b"])
                rect(xa, y_top + r["h"] / 2 - eh * 0.3, max(0.6, xb - xa), max(0.6, eh * 0.6), color, alpha=0.45)
            for s in exs:
                xa, xb = xz(s["a"]), xz(s["b"])
                rect(xa, ey, max(1.0, xb - xa), eh, color)

    # ---- gene-model reference track (static/app.js drawAxis) ----
    hline(track_x0, track_x1, ruler_base, "#7a8696", lw=0.9)
    for tk in np.linspace(lo, hi, max(4, round((track_x1 - track_x0) / 120))):
        px = xz(tk)
        ax.plot([px, px], [ruler_base - 5, ruler_base], color="#7a8696", linewidth=0.8, solid_capstyle="butt")
        ax.text(px, ruler_base - 8, fmt(tk), fontsize=8, ha="center", va="bottom", color="#5a6675")
    cy = ref_top + REF_H / 2
    ref_lo, ref_hi = max(track_x0, xz(lo)), min(track_x1, xz(hi))
    hline(ref_lo, ref_hi, cy, "#9aa4b2", lw=1.0)
    strand = gm.get("strand")
    if strand and ref_hi - ref_lo >= 40:
        d = -1 if strand == "-" else 1
        px = ref_lo + 24
        while px < ref_hi - 6:
            ax.plot([px - 3 * d, px, px - 3 * d], [cy - 3, cy, cy + 3], color="#9aa4b2", linewidth=0.7, solid_capstyle="round")
            px += 46
    last_lbl = None
    for b in gm.get("blocks", []):
        xa, xb = xz(b["start"]), xz(b["end"])
        if xb < track_x0 or xa > track_x1:
            continue
        rect(xa, ref_top + 3, max(1.6, xb - xa), REF_H - 6, REF_FILL)
        lbl = str(b.get("exon_id") or "").split(".")[0]
        if lbl and lbl != last_lbl and (xb - xa) > 7:
            ax.text((xa + xb) / 2, ref_top + REF_H + 10, lbl, fontsize=7, ha="center", va="top", color="#5a6675")
            last_lbl = lbl
    ax.text(track_x0, ref_top + REF_H + 11, gm.get("chrom") or "", fontsize=8, ha="left", va="top", color="#5a6675")

    # ---- PDF-only bottom isoform-structure track (one row per plotted isoform; UCSC CDS-tall / UTR-thin,
    #      introns as connecting lines), aligned to the same x-axis as the molecule pileup above. ----
    if iso_rows:
        ax.text(LEFT - 8, struct_top + 9, "Isoform structures", fontsize=8.5, weight="bold",
                ha="right", va="center", color="#1a2230")
        ystruct = struct_top + 16
        for m in iso_rows:
            color = colof(m)
            merged = m["merged"]
            exs, ins = merged["exons"], merged["introns"]
            rc = ystruct + ISO_ROW / 2.0                          # row centre
            ex_y = ystruct + (ISO_ROW - ISO_EXH) / 2.0            # tall (CDS) exon top
            utr_h = max(3.0, ISO_EXH * 0.5)
            utr_y = ex_y + (ISO_EXH - utr_h) / 2.0                # thin (UTR) exon top
            if exs:                                              # backbone line across the isoform span (introns)
                hline(xz(min(e["a"] for e in exs)), xz(max(e["b"] for e in exs)), rc, color, lw=0.7, alpha=0.55)
            edge = darken(color)                                  # exon outline = same-hue darker shade of fill
            cds_lo, cds_hi = m.get("cds_min"), m.get("cds_max")
            for s in exs:
                xa, xb = xz(s["a"]), xz(s["b"])
                if cds_lo is not None and cds_hi is not None:
                    rect(xa, utr_y, max(0.8, xb - xa), utr_h, color, ec=edge, lw=0.5)     # UTR (thin) full exon
                    ca, cb = max(s["a"], cds_lo), min(s["b"], cds_hi)
                    if cb > ca:
                        cxa, cxb = xz(ca), xz(cb)
                        rect(cxa, ex_y, max(1.0, cxb - cxa), ISO_EXH, color, ec=edge, lw=0.5)  # CDS (tall) coding overlap
                else:
                    rect(xa, ex_y, max(1.0, xb - xa), ISO_EXH, color, ec=edge, lw=0.5)    # no CDS info -> full height
            ax.text(LEFT - 8, rc, truncate(m.get("final_isoform_id") or m.get("isoform_id"), 22),
                    fontsize=7, ha="right", va="center", color=(KNOWN if m.get("known") else NOVEL))
            ystruct += ISO_ROW

    fig.savefig(out_path, format="pdf")
    plt.close(fig)
    print(f"WROTE {out_path}")
    return True


def parse_group_spec(spec):
    """'HSPC=HSC-1,HSC-2,MPP-1;DC=cDC1,cDC2-1' -> OrderedDict{label: [members]}. A part with no '=' is
    self-labelled by its members joined with '+', so 'HSC-1,HSC-2;cDC1' also works."""
    groups = OrderedDict()
    for part in (spec or "").split(";"):
        part = part.strip()
        if not part:
            continue
        label, _, csv = part.partition("=") if "=" in part else ("", "", part)
        members = [x.strip() for x in csv.split(",") if x.strip()]
        if not members:
            continue
        groups[(label.strip() or "+".join(members))] = members
    return groups


_REC_FIELDS = ("isoform_id", "final_isoform_id", "known", "cluster_index", "color_index",
               "exon_segments", "protein_length", "nmd_status")


def compose_panels(ctx, res, cell_groups, cond_groups, max_isoforms):
    """Re-aggregate the SHARED, union-clustered molecules of one query_reads() call into arbitrary panels:
    panel = (cell-state group) x (condition group). Each molecule's count in a panel is the sum of its
    per-(cell_type, sample) detections that fall inside that panel's cell states AND conditions; cluster_index
    / color_index are inherited from the single union clustering, so colors stay consistent across panels."""
    master = OrderedDict()                      # uid -> record (+ merged detection list)
    for p in res.get("panels", []):
        for m in p.get("molecules", []):
            rec = master.get(m["isoform_id"])
            if rec is None:
                rec = {k: m.get(k) for k in _REC_FIELDS}
                rec["_det"] = []
                master[m["isoform_id"]] = rec
            rec["_det"].extend(m.get("detections", []))

    def cond_of(lib):
        return ctx.samples.get(lib, {}).get("group")

    show_cond = not (len(cond_groups) == 1 and next(iter(cond_groups)) == "all")
    panels = []
    for cl, states in cell_groups.items():
        sset = set(states)
        for condl, conds in cond_groups.items():
            cset = set(conds) if conds else None
            mols = []
            for rec in master.values():
                cnt, det = 0.0, []
                for d in rec["_det"]:
                    if d.get("cell_type") in sset and (cset is None or cond_of(d.get("sample")) in cset):
                        cnt += float(d.get("count") or 0.0)
                        det.append(d)
                if cnt <= 0:
                    continue
                mm = {k: rec.get(k) for k in _REC_FIELDS}
                mm["count"] = round(cnt, 3)
                mm["detections"] = det
                mols.append(mm)
            mols.sort(key=lambda m: -m["count"])           # top-N by reads (mirror query_reads' per-panel cap)
            mols = mols[:max_isoforms]
            label = cl + (f" | {condl}" if show_cond else "")
            panels.append({"condition": label,
                           "n_molecules": int(round(sum(m["count"] for m in mols))),
                           "molecules": mols})
    return {"gene": res.get("gene"), "symbol": res.get("symbol"),
            "gene_model": res.get("gene_model"), "panels": panels}


def filter_low_fraction_isoforms(res, min_pct):
    """Drop final-isoform clusters whose pooled reads are < min_pct% of the gene's TOTAL reads in the
    rendered selection (declutter rare isoforms). Applied GLOBALLY across panels so the surviving isoform
    set -- and thus the shared cluster colors/ordering -- stays consistent between panels."""
    panels = res.get("panels", [])
    by_iso, total = {}, 0.0
    for p in panels:
        for m in p.get("molecules", []):
            c = float(m.get("count", 0) or 0); total += c
            k = m.get("final_isoform_id") or m.get("isoform_id")
            by_iso[k] = by_iso.get(k, 0.0) + c
    if total <= 0:
        return res
    keep = {k for k, v in by_iso.items() if 100.0 * v / total >= min_pct}
    for p in panels:
        p["molecules"] = [m for m in p.get("molecules", [])
                          if (m.get("final_isoform_id") or m.get("isoform_id")) in keep]
        p["n_molecules"] = int(round(sum(float(m.get("count", 0) or 0) for m in p["molecules"])))
    return res


def build_arg_parser():
    ap = argparse.ArgumentParser(
        prog="python -m altanalyze3.components.visualization.isv_web.isv_web_remote",
        description="Backend ISV-web read-level molecule pileup -> vector PDF (no web server).",
    )
    ap.add_argument("--run_dir", required=True, help="dataset dir holding _isv_web_cache/mol_index/ + the metadata")
    ap.add_argument("--gene", required=True, help="ENSG id OR gene symbol (e.g. ENSG00000145779 or TNFAIP8)")
    ap.add_argument("--cell_states", default=None, help="comma-separated cell states; ONE pileup panel each (e.g. HSC-1,cDC2-2)")
    ap.add_argument("--cell_groups", default=None,
                    help="COMBINE cell states into contrasted panels: 'HSPC=HSC-1,HSC-2,MPP-1;DC=cDC1,cDC2-1,cDC2-2'. "
                         "One panel per group, molecules pooled across its member states. (label= is optional.)")
    ap.add_argument("--condition_groups", default=None,
                    help="COMBINE/contrast by condition too: 'young=young;aged=aged'. Panels become "
                         "cell_group x condition_group. Omit to pool all conditions into each panel.")
    ap.add_argument("--out", default=None, help="output PDF (default: <run_dir>/ISV_<gene>_<states>_reads.pdf)")
    ap.add_argument("--conditions", default=None, help="restrict to these sample groups (comma-separated; default: ALL groups)")
    ap.add_argument("--metadata", default=None, help="metadata TSV (default: <run_dir>/mds_metadata_bam.txt)")
    ap.add_argument("--gene_model", default=None, help="exon coords TSV (default: <run_dir>/Hs_Ensembl_exon.txt)")
    ap.add_argument("--gene_symbol", default=None, help="ENSG->symbol annotations (default: <run_dir>/cellHarmony/Hs_Ensembl-annotations.txt)")
    ap.add_argument("--max_isoforms", type=int, default=1_000_000,
                    help="max molecules shown per panel -- default is effectively UNLIMITED so the pileup "
                         "shows ALL reads and the figure auto-shrinks to fit; lower it only to thin a very "
                         "high-expression gene.")
    ap.add_argument("--min_isoform_pct", type=float, default=0.0,
                    help="exclude final isoforms whose reads are < this %% of the gene's total reads in the "
                         "rendered selection (declutter rare isoforms; default 0 = keep all)")
    ap.add_argument("--width", type=int, default=1100, help="figure width in points (default 1100)")
    ap.add_argument("--intron_scale", "--intron-scale", dest="intron_scale", type=float, default=1.0,
                    help="relative INTRON display size to enlarge exons (1.0=true genomic scale; "
                         "0.5=introns at 50%%, 0.2=20%%), like isoform_structure_view --intron-scale")
    ap.add_argument("--no_isoform_track", "--no-isoform-track", dest="isoform_track", action="store_false",
                    default=True, help="omit the PDF-only bottom isoform-structure track (CDS-aware, "
                                       "drawn by default; aligned to the molecule pileup x-axis)")
    ap.add_argument("--export_molecules", default=None,
                    help="ALSO write a TSV of (molecule_id, final_isoform_id, cell_state, sample, count) for "
                         "EVERY molecule in the selection (UNCAPPED, ignores --max_isoforms) to this path -- "
                         "for reconciliation against the raw reads.db index and the source pseudobulk counts.")
    return ap


def _export_molecule_tsv(res, path):
    """Write every molecule's (cell_state, sample) detection to a TSV -- one row per
    (molecule, cell_state, sample) -- for reconciliation against the raw reads.db index and the source
    pseudobulk counts. Callers pass an UNCAPPED query so nothing is dropped."""
    import csv
    n = 0
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["molecule_id", "final_isoform_id", "cell_state", "sample", "count"])
        for p in res.get("panels", []):
            for m in p.get("molecules", []):
                mid = m.get("isoform_id"); fin = m.get("final_isoform_id")
                for d in m.get("detections", []):
                    w.writerow([mid, fin, d.get("cell_type"), d.get("sample"),
                                round(float(d.get("count") or 0), 3)])
                    n += 1
    print(f"[isv_web_remote] exported {n} molecule rows -> {path}", flush=True)


def main(argv=None):
    a = build_arg_parser().parse_args(argv)
    from . import data_api as da

    metadata = a.metadata or os.path.join(a.run_dir, "mds_metadata_bam.txt")
    gene_model = a.gene_model or os.path.join(a.run_dir, "Hs_Ensembl_exon.txt")
    gene_symbol = a.gene_symbol or os.path.join(a.run_dir, "cellHarmony", "Hs_Ensembl-annotations.txt")
    conditions = [c.strip() for c in a.conditions.split(",") if c.strip()] if a.conditions else None
    if not (a.cell_states or a.cell_groups):
        sys.exit("ERROR: provide --cell_states (one panel per state) or --cell_groups (combined panels).")
    for f in (metadata, gene_model):
        if not os.path.exists(f):
            sys.exit(f"ERROR: not found: {f}")

    print(f"[isv_web_remote] loading run context from {a.run_dir} (pre-built indexes, lite) ...", flush=True)
    # Lightweight load for the read-level (molecule) view: run ONLY the loaders the molecule pileup needs,
    # SKIPPING the heavy/unused steps -- the 2.1M-isoform coding-regions index (CDS/UTR shading, never drawn
    # here), the protein summary/FASTA, and the interactions network. The web server (run.py) still calls
    # the full RunContext.load_all(); this path leaves it untouched.
    ctx = da.RunContext(a.run_dir, metadata, gene_model, gene_symbol)
    ctx._load_metadata_and_clusters()
    ctx._load_combined_matrix()
    ctx._load_gene_model()
    ctx._load_gene_symbols()
    ctx._index_genes_from_var_names()
    ctx._build_viewer_inputs()
    if a.isoform_track:
        # the bottom isoform-structure track needs CDS spans -> load the coding-regions index (skipped on
        # the fast molecule-only path). Heavy (~2.1M isoforms) but only for the PDF export.
        try:
            ctx._index_coding_regions()
        except Exception as _e:
            print(f"[isv_web_remote] WARN: coding regions unavailable ({_e}); structure track without CDS shading", flush=True)
    gene = ctx.resolve_gene(a.gene) or a.gene

    if a.cell_groups or a.condition_groups:
        # ---- COMBINED / contrasted panels: (cell-state group) x (condition group) ----
        cell_groups = (parse_group_spec(a.cell_groups) if a.cell_groups
                       else OrderedDict((s.strip(), [s.strip()]) for s in a.cell_states.split(",") if s.strip()))
        if a.condition_groups:
            cond_groups = parse_group_spec(a.condition_groups)
        elif conditions:
            cond_groups = OrderedDict([("all", conditions)])    # restrict (and pool) to these conditions
        else:
            cond_groups = OrderedDict([("all", None)])          # all conditions pooled into each panel
        all_states = sorted({s for v in cell_groups.values() for s in v})
        slug = "_".join(list(cell_groups) +
                        ([] if (len(cond_groups) == 1 and "all" in cond_groups) else list(cond_groups)))
        out = a.out or os.path.join(a.run_dir, f"ISV_{a.gene}_{slug}_reads.pdf")
        print(f"[isv_web_remote] gene={a.gene}->{gene} cell_groups={dict(cell_groups)} "
              f"condition_groups={ {k: (v or 'ALL') for k, v in cond_groups.items()} } "
              f"grouping={ctx.molecule_grouping}", flush=True)
        # ONE query over the union of all member states (huge cap so nothing is dropped before re-aggregation)
        # -> a single shared final-isoform clustering, then re-composed into the requested panels.
        base = da.query_reads(ctx, gene, all_states, conditions=None,
                              max_isoforms=1_000_000, panel_by="cell_type")
        res = compose_panels(ctx, base, cell_groups, cond_groups, a.max_isoforms)
    else:
        # ---- ONE panel per cell state (original behaviour) ----
        cell_states = [s.strip() for s in a.cell_states.split(",") if s.strip()]
        out = a.out or os.path.join(a.run_dir, f"ISV_{a.gene}_{'_'.join(cell_states)}_reads.pdf")
        print(f"[isv_web_remote] gene={a.gene}->{gene} states={cell_states} "
              f"conditions={conditions or 'ALL'} grouping={ctx.molecule_grouping}", flush=True)
        res = da.query_reads(ctx, gene, cell_states, conditions=conditions,
                             max_isoforms=a.max_isoforms, panel_by="cell_type")

    if a.export_molecules:
        # UNCAPPED re-query so EVERY molecule is exported (the render caps at max_isoforms; the
        # reconciliation TSV must not). One row per (molecule, cell_state, sample).
        exp_states = ([s.strip() for s in a.cell_states.split(",") if s.strip()] if a.cell_states
                      else sorted({s for v in parse_group_spec(a.cell_groups).values() for s in v}))
        full = da.query_reads(ctx, gene, exp_states, conditions=conditions,
                              max_isoforms=10**9, panel_by="cell_type")
        _export_molecule_tsv(full, a.export_molecules)

    if a.min_isoform_pct and a.min_isoform_pct > 0:
        res = filter_low_fraction_isoforms(res, a.min_isoform_pct)
    if a.isoform_track:
        # attach CDS genomic span to each plotted molecule for the bottom structure track (UCSC CDS/UTR)
        for p in res.get("panels", []):
            for m in p.get("molecules", []):
                cds = ctx.coding_region(m.get("final_isoform_id") or m.get("isoform_id"))
                m["cds_min"], m["cds_max"] = (cds[0], cds[1]) if cds else (None, None)
    print("[isv_web_remote] source={} min_isoform_pct={} panels={}".format(
        res.get("source", "index"), a.min_isoform_pct,
        ", ".join(f"{p['condition']}:{p['n_molecules']}" for p in res.get("panels", []))), flush=True)
    sys.exit(0 if render(res, out, width=a.width, intron_scale=a.intron_scale,
                         isoform_track=a.isoform_track) else 2)


if __name__ == "__main__":
    main()
