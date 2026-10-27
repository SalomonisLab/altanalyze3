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
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"]
# Arial ships here as regular and bold only, with no italic face, so italic text asks for a Helvetica-metric
# family that has one. Order keeps a real Helvetica or Arial Italic first where the system provides it.
ITALIC_FAMILY = ["Helvetica", "Nimbus Sans", "Liberation Sans", "Arial", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# ---- drawReads() layout constants (static/app.js R = {...}) ----
LEFT, RIGHT, READ_PX, ROW_MIN, ROW_CAP = 152, 24, 1.6, 0.9, 60
MAXPANEL, PANEL_GAP, PAD, HDR, RULER_H, REF_H, MAX_ROWS, RMARGIN = 340, 14, 8, 22, 30, 26, 4000, 16
ROW_VIS = 8.0          # minimum drawn height of one read row, px, enforced only where it fits in MAXPANEL
KNOWN, NOVEL, REF_FILL = "#1d5fa8", "#b0306b", "#3d4a5c"
CLUSTER_COLORS = [
    "#2b6cb0", "#c0531f", "#2f8a4e", "#7a4fb5", "#c2306b", "#0e8a9c", "#b8902a",
    "#4a5bbf", "#cf5aa0", "#2f9b8a", "#9a5b2e", "#56657a", "#1f7a5a", "#a23c5e",
    "#3d7ec2", "#8a7a1f", "#6d4ca0", "#1d8f74", "#bf6a2a", "#8a2f52", "#3a6ea5",
    "#7b9b3a", "#a64d79", "#2d7d8a",
]


TWO_ISOFORM_COLORS = ["#E9B940", "#4884BC"]   # default fills when a gene shows EXACTLY two isoforms (gold, blue)
THREE_ISOFORM_COLORS = ["#4884BC", "#E9B940", "#2f8a4e"]  # trans-splicing scheme: blue, gold, green


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


def isoform_label(m, with_length=True, id_max=None, style="paren"):
    """An isoform's name with its protein length.

    Nathan, 2026-10-15: "the AA length only goes in parentheses for the view that doesn't already
    have the AA". The read-panel group labels never carried a length, so they take
    `<id> (568aa)`, with ", NMD" inside the brackets. The isoform structure rows already carried
    one and keep their own `<id>|568` form, style="bar".
    """
    iid = str(m.get("final_isoform_id") or m.get("isoform_id") or "")
    # the id is shortened, never the length: truncating the whole label hid the "(635aa)" it was
    # added for
    if id_max and len(iid) > id_max:
        iid = iid[:max(1, id_max - 1)] + "\u2026"
    if not with_length:
        return iid
    aa = m.get("protein_length")
    nmd = str(m.get("nmd_status") or "").strip() in ("NMD", "Potential-NMD")
    if aa in (None, ""):
        return iid + (" (NMD)" if nmd else "")
    try:
        aa = int(float(aa))
    except (TypeError, ValueError):
        return iid + (" (NMD)" if nmd else "")
    if style == "bar":
        return f"{iid}|{aa}" + (" (NMD)" if nmd else "")
    return f"{iid} ({aa}aa, NMD)" if nmd else f"{iid} ({aa}aa)"


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
        # clamp BOTH ends into [lo, hi]. Clamping only the low end leaves `a` above `hi` for an exon beyond
        # the window, and the "intron gap before this exon" below then measures (a - cur), a span that is not
        # on the axis at all. On RUNX1 over chr21:35,048,300-35,049,600 the next exon block above the window
        # added 36,656 phantom units to a 1,044-unit axis and squashed the whole panel into its left 2.8%.
        # Without --region, hi is always the end of a drawn exon, so no interval could start above it.
        a, b = min(max(a, lo), hi), min(max(b, lo), hi)
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


def color_key(m):
    """Colour and row-order key: the explicit color_index when one was assigned (--isoforms), else the
    cluster index. Identical to the former cluster_index ordering whenever no color_index is set."""
    return m.get("color_index") if m.get("color_index") is not None else m.get("cluster_index", 0)


def parse_region(text):
    """'chr11:1,887,200-1,892,300', '1887200-1892300' or '1887200..1892300' -> (lo, hi); None stays None.

    A window is needed when the event sits in a small part of a long gene. The x-axis is otherwise the union
    of every plotted molecule, so a 209 bp splice inside one exon of a 261 kb gene is drawn under a point wide.
    """
    if not text:
        return None
    s = str(text).strip().replace(",", "").replace("..", "-")
    if ":" in s:
        s = s.split(":", 1)[1]
    parts = [q for q in s.split("-") if q != ""]
    if len(parts) != 2:
        raise ValueError("--region expects START-END or chr:START-END, got {!r}".format(text))
    try:
        a, b = int(parts[0]), int(parts[1])
    except ValueError:
        raise ValueError("--region bounds must be integers, got {!r}".format(text))
    if b <= a:
        raise ValueError("--region END must exceed START, got {!r}".format(text))
    return a, b


def mol_extent(merged):
    """(lo, hi, span) over EVERY segment of a molecule, exons and retained introns alike.
    Returns (0, 0, 0) for a molecule with no drawable segment."""
    segs = list((merged or {}).get("exons") or []) + list((merged or {}).get("introns") or [])
    if not segs:
        return (0, 0, 0)
    lo = min(s["a"] for s in segs)
    hi = max(s["b"] for s in segs)
    return (lo, hi, hi - lo)


def iso_struct_extent(m):
    """(lo, hi, exonic_bp) of the structure the figure DRAWS for this isoform: its registered
    combined.gff.gz exons when the isoform is registered, else the molecule's own merged segments.
    (0, 0, 0) when nothing is drawable."""
    reg = m.get("registered_exons")
    if reg:
        segs = [(min(int(a), int(b)), max(int(a), int(b))) for a, b in reg]
    else:
        mg = m.get("merged") or merge_segments(m.get("exon_segments"))
        segs = [(s["a"], s["b"]) for s in list(mg.get("exons") or []) + list(mg.get("introns") or [])]
    if not segs:
        return (0, 0, 0)
    return (min(a for a, _ in segs), max(b for _, b in segs), sum(b - a for a, b in segs))


def iso_order_map(mols):
    """final isoform id -> rank, LONGEST ISOFORM FIRST, over every molecule drawn in the figure.

    Isoform order used to be the cluster index, which the read index assigns by abundance, so the
    structure track and the colour groups above it ran most-read-first and the track read as noise:
    on FLNB the 163,830 bp ENST00000295956 was followed by the 12,979 bp ENST00000419752 and then the
    26,420 bp ENST00000481470. Length belongs to the isoform, not to the row, so ONE rank orders both
    the read rows and the structure rows and the two line up. Ties fall to exonic bp, then the 5' edge,
    then the id, so the order is deterministic. Without the structure track no registered structure is
    attached and the rank falls back to each isoform's representative read extent.
    """
    best = {}
    for m in mols:
        k = str(m.get("final_isoform_id") or m.get("isoform_id"))
        lo, hi, bp = iso_struct_extent(m)
        cur = best.get(k)
        if cur is None or (hi - lo, bp) > (cur[1] - cur[0], cur[2]):
            best[k] = (lo, hi, bp)
    ranked = sorted(best.items(), key=lambda kv: (-(kv[1][1] - kv[1][0]), -kv[1][2], kv[1][0], kv[0]))
    return {k: i for i, (k, _v) in enumerate(ranked)}


def group_key(m):
    """Which isoform group a row belongs to, and where that group sits: the length rank when render()
    stamped one, else the colour key. An explicit --isoforms order is never overridden."""
    r = m.get("_iso_rank")
    return color_key(m) if r is None else r


def row_sort_key(m):
    """Row order inside a pileup: isoform group first (by isoform LENGTH, see iso_order_map), then the
    LONGEST molecule first inside that group.

    The previous key was (colour, -count, isoform_id). Long-read molecules almost all carry count 1,
    so the count term is constant and the order collapsed to lexicographic molecule id, i.e. arbitrary.
    Ordering by genomic extent instead stacks each colour group into a staircase, so where reads start
    and stop is readable off the panel. Ties fall back to the 5' edge and then the id, so the order is
    deterministic and does not depend on dict iteration.
    """
    lo, hi, span = m.get("_span") or mol_extent(m.get("merged"))
    return (group_key(m), -span, lo, hi, str(m.get("isoform_id")))


def segs_in_window(segs, lo, hi):
    """The segments that intersect [lo, hi]. build_xz clamps coordinates INTO the window, so a segment lying
    wholly outside would otherwise be drawn as a minimum-width tick pinned to the window edge. With lo/hi
    taken from the molecules' own extent every segment qualifies, so this is a no-op without --region."""
    return [s for s in segs if not (s["b"] < lo or s["a"] > hi)]


def mol_in_window(merged, lo, hi):
    """True when any exon or intron segment of this molecule intersects [lo, hi]."""
    return bool(segs_in_window(list(merged.get("exons") or []) + list(merged.get("introns") or []), lo, hi))


ISO_ROW, ISO_EXH = 12.0, 7.0          # bottom isoform-structure track: row height / coding-exon height
# figure preset only (structure_labels_above): one label line, label-to-its-own-structure gap, and the gap
# between two adjacent structures. The label must sit CLOSER to the structure it names than the structures
# sit to each other, or the reader pairs the wrong label with the wrong structure.
ISO_LBL_H, ISO_LBL_GAP, ISO_STRUCT_GAP = 12.0, 3.0, 10.0


def struct_rows_geometry(n, FS, labels_above, y0=0.0):
    """y positions of the bottom isoform-structure track: [(row centre, label y, label va), ...], total height.

    labels_above=False keeps the original fixed ISO_ROW grid with no labels in the track (they go in the left
    margin). labels_above=True is the figure preset: every row carries its label on its own line, the last row
    of a multi-row track below its structure so the structures of one event sit next to each other.
    """
    if not n:
        return [], 0.0
    if not labels_above:
        return [(y0 + i * ISO_ROW * FS + ISO_ROW * FS / 2.0, None, None) for i in range(n)], n * ISO_ROW * FS
    exh, lbl, g_lbl, g_str = ISO_EXH * FS, ISO_LBL_H * FS, ISO_LBL_GAP * FS, ISO_STRUCT_GAP * FS
    rows, y = [], y0
    for i in range(n):
        if (i == n - 1) and n > 1:                     # label BELOW its structure
            rc = y + exh / 2.0
            rows.append((rc, rc + exh / 2.0 + g_lbl, "top"))
            y = rc + exh / 2.0 + g_lbl + lbl
        else:                                          # label ABOVE its structure
            rc = y + lbl + g_lbl + exh / 2.0
            rows.append((rc, y + lbl, "bottom"))
            y = rc + exh / 2.0
        y += g_str if i < n - 1 else 0.0
    return rows, y - y0


def read_label_positions(centers, body_height, gap):
    """Keep every ordered label within the panel with enough space for its text."""
    if not centers:
        return []
    if body_height + 1e-9 < len(centers) * gap:
        raise ValueError("Panel too short for every isoform label")
    positions = []
    for center in centers:
        positions.append(max(gap / 2, center, positions[-1] + gap if positions else 0))
    positions[-1] = min(positions[-1], body_height - gap / 2)
    for i in range(len(positions) - 2, -1, -1):
        positions[i] = min(positions[i], positions[i + 1] - gap)
    return positions


def merge_locus_segments(segments):
    """Use the established exon merging separately within each annotated gene."""
    groups = OrderedDict()
    for segment in segments or []:
        groups.setdefault(segment["gene"], []).append(segment)
    merged = {"exons": [], "introns": []}
    for gene, parts in groups.items():
        for kind, values in merge_segments(parts).items():
            merged[kind].extend(dict(s, gene=gene) for s in values)
    return merged


def prepare_locus_coordinates(response, intron_scale, trim_to_reads, isoform_track):
    """Project gene-qualified geometry onto independent sections of the native ISV axis.

    Only a drawing copy is changed. Every real read, count and final assignment is
    retained. Intron compression is the native build_xz calculation, per locus.
    Registered structures must identify their genes; untagged genomic pairs are
    ambiguous across chromosomes and cannot substitute for the registered model.
    """
    from copy import deepcopy
    from .locus_view import coding_segments_for_loci, locus_axis_groups
    res = deepcopy(response)
    records = list({id(m): m for m in
                    [r for p in res.get("panels", []) for r in p.get("molecules", [])]
                    + res.get("reference_isoforms", [])}.values())
    if not records:
        return res, []
    primary = res["gene"]
    models = {m["gene"]: m for m in res["gene_models"]}
    parts = OrderedDict()
    for m in records:
        if m.get("unplaced_genes"):
            raise ValueError("Cannot place annotated partner loci: " + str(m["unplaced_genes"]))
        if isoform_track and not m.get("registered_exon_segments"):
            raise ValueError("Gene-qualified registered structure required for multi-locus ISV: "
                             + str(m.get("final_isoform_id")))
        for s in list(m.get("exon_segments", [])) + list(m.get("registered_exon_segments", [])):
            gene = s.get("gene") or primary
            s["gene"] = gene
            parts.setdefault(gene, []).append(s)
    frames, transforms, by_gene_frame = [], {}, {}
    # Equal-width sections are explicitly separate loci, never genomic distance.
    slot, gap = 1_000_000.0, 120_000.0
    for i, genes in enumerate(locus_axis_groups([models[g] for g in parts])):
        segments = [s for g in genes for s in parts[g]]
        bounds = [(min(s["start"], s["end"]), max(s["start"], s["end"])) for s in segments]
        lo, hi = min(a for a, b in bounds), max(b for a, b in bounds)
        if not trim_to_reads:
            lo = min(lo, *(models[g]["extent"][0] for g in genes))
            hi = max(hi, *(models[g]["extent"][1] for g in genes))
        reference = [(min(b["start"], b["end"]), max(b["start"], b["end"]))
                     for g in genes for b in models[g].get("blocks", [])]
        start = i * (slot + gap)
        project = build_xz(_merge_ivs(bounds + reference), lo, hi, max(0.01, float(intron_scale)), start, start + slot)
        for gene in genes:
            model = models[gene]
            transforms[gene] = project
            projected = dict(model, blocks=[dict(b, start=project(b["start"]), end=project(b["end"]))
                                           for b in model.get("blocks", [])
                                           if b["end"] >= lo and b["start"] <= hi], extent=[start, start + slot])
            by_gene_frame[gene] = dict(model=projected, start=start, end=start + slot,
                                      genomic_extent=[lo, hi], ticks=np.linspace(lo, hi, 3), project=project)
    frames = [by_gene_frame[g] for g in parts]
    def chrom_key(value):
        key = str(value).upper().removeprefix("CHR")
        return "M" if key == "MT" else key
    # A reference row describes the surrounding chromosome, not only the span
    # of this gene. Carry it through every displayed window on that chromosome,
    # showing only real annotated exon intersections within each window.
    for frame in frames:
        model = models[frame["model"]["gene"]]
        windows = []
        seen_windows = set()
        for axis in frames:
            if chrom_key(model["chrom"]) != chrom_key(axis["model"]["chrom"]):
                continue
            window_key = axis["start"], axis["end"]
            if window_key in seen_windows:
                continue
            seen_windows.add(window_key)
            lo, hi = axis["genomic_extent"]
            project = axis["project"]
            blocks = [dict(b, start=project(max(lo, b["start"])),
                           end=project(min(hi, b["end"]))) for b in model.get("blocks", [])
                      if b["end"] >= lo and b["start"] <= hi]
            windows.append(dict(model=dict(model, blocks=blocks), start=axis["start"], end=axis["end"],
                                genomic_extent=axis["genomic_extent"], ticks=axis["ticks"],
                                project=axis["project"]))
        frame["reference_windows"] = windows
    for m in records:
        reg = m.get("registered_exon_segments") or m.get("exon_segments", [])
        if len({s["gene"] for s in reg}) > 1:
            coding = coding_segments_for_loci(reg, models, m.get("cds_min"),
                                              m.get("cds_max"), m.get("cds_chrom"))
            if coding is not None:
                m["_locus_cds_intervals"] = [dict(gene=s["gene"],
                    a=transforms[s["gene"]](s["start"]),
                    b=transforms[s["gene"]](s["end"])) for s in coding]
        for field in ("exon_segments", "registered_exon_segments"):
            runs = []
            for s in m.get(field, []):
                if not runs or runs[-1][0]["gene"] != s["gene"]:
                    runs.append([])
                runs[-1].append(s)
            joins = []
            for left, right in zip(runs, runs[1:]):
                lg, rg = left[0]["gene"], right[0]["gene"]
                if m.get("alignment_coordinates_verified"):
                    # Original alignment blocks are supplied in genomic order.
                    a, b = max(s["end"] for s in left), min(s["start"] for s in right)
                else:
                    a = min(s["start"] for s in left) if models[lg]["strand"] == "-" else max(s["end"] for s in left)
                    b = max(s["end"] for s in right) if models[rg]["strand"] == "-" else min(s["start"] for s in right)
                joins.append((transforms[lg](a), transforms[rg](b)))
            m["_" + field + "_joins"] = joins
    # Many reads share one registered list. Transform each segment once, after
    # recording every join in the original coordinates (deepcopy preserves aliases).
    segments = {id(s): s for m in records for field in ("exon_segments", "registered_exon_segments")
                for s in m.get(field, [])}
    for s in segments.values():
        s["genomic_start"], s["genomic_end"] = s["start"], s["end"]
        s["start"], s["end"] = transforms[s["gene"]](s["start"]), transforms[s["gene"]](s["end"])
    for m in records:
        reg = m.get("registered_exon_segments")
        if reg:
            m["_registered_locus_merged"] = merge_locus_segments(reg)
            # Retained intron sequence belongs to the registered RNA too. The
            # native structure track draws carried sequence as solid blocks.
            m["_registered_locus_merged"]["exons"].extend(m["_registered_locus_merged"]["introns"])
            m["_registered_locus_merged"]["introns"] = []
            m["registered_exons"] = [(s["start"], s["end"]) for s in reg]
        genes = {s["gene"] for s in reg or m.get("exon_segments", [])}
        # A single genomic CDS span cannot specify coding sequence across loci.
        if len(genes) == 1:
            project = transforms[next(iter(genes))]
            for key in ("cds_min", "cds_max"):
                if m.get(key) is not None:
                    m[key] = project(m[key])
        else:
            m["cds_min"] = m["cds_max"] = None
    res["gene_model"] = dict(blocks=[b for f in frames for b in f["model"]["blocks"]],
                             extent=[min(f["start"] for f in frames), max(f["end"] for f in frames)])
    return res, frames


def _architecture_structure_length(model, track):
    """Use audited sequence length on integrated figures, retaining source metadata."""
    observed = model.get('protein_length')
    literal = int(track['protein_aa'])
    reported = int(track['reported_aa'] or literal)
    if observed not in (None, '') and int(observed) not in (literal, reported):
        raise ValueError('ISV protein length differs from audited architecture: ' + track['isoform_id'])
    return dict(isoform_id=track['isoform_id'], source_RNA_protein_length=observed,
                source_architecture_reported_aa=track['reported_aa'],
                displayed_literal_sequence_aa=literal)


def render(res, out_path, width=1100, intron_scale=1.0, isoform_track=False,
           font_scale=1.0, show_ruler=True, show_read_labels=True, trim_to_reads=False,
           title_gene=None, title_event=None, show_read_counts=True, show_structure_header=True,
           label_protein_length=True, dense_exon_labels=False, structure_labels_above=False,
           region=None, iso_order="length", caption=None, hide_bleeding_reads=False, bleed_tolerance=0,
           protein_features=None):
    """Render query_reads() panels as the drawReads molecule pileup to a single-page vector PDF.
    intron_scale<1.0 shrinks intron display spans (enlarging exons); 1.0 = true genomic scale.
    font_scale multiplies EVERY font size and the label margins and header heights that depend on them,
    so a panel stays readable when it is shrunk into a multi-gene figure. show_ruler=False drops the
    coordinate ruler and the chromosome label, show_read_labels=False drops the per-read isoform labels
    inside the pileups (the structure track still names each isoform once), and trim_to_reads=True sets
    the x-axis from the plotted reads instead of the whole gene, so gene regions with no reads are not drawn.
    isoform_track=True appends a PDF-only structure track (one row per plotted isoform; UCSC CDS-tall /
    UTR-thin, introns as lines) under the gene model, aligned to the same x-axis -- NOT shown in the web app.
    title_gene is drawn above the panels in italic (the gene symbol) and title_event on the line under it
    (the splicing event, e.g. ENSG00000169756:E3.1-E17.1|E1.2-E17.1). show_read_counts=False drops the
    "· N reads" note from each panel header, show_structure_header=False drops the "Isoform structures"
    caption, label_protein_length=True writes each structure row as <isoform id>|<protein length> with
    (NMD) appended when the isoform is NMD, and dense_exon_labels=True labels every exon block wide enough
    to carry text, staggered over two rows, instead of only blocks wider than the font.
    structure_labels_above=True writes each structure row's label ABOVE its own structure, in that
    structure's own colour, instead of to the left in the known/novel colours; with no labels left of the
    track the left margin collapses and the panels take the full width.
    region=(lo, hi) draws ONLY that genomic window: the x-axis is the window, reads that do not overlap it
    are dropped, and segments outside it are not drawn. It overrides trim_to_reads, whose axis is the union
    of every plotted molecule and therefore hides an event that sits in a small part of a long gene."""
    from .locus_view import is_trans_splicing_view
    if protein_features and not isoform_track:
        raise ValueError('Protein architecture requires the registered isoform structure track')
    if protein_features:
        from copy import deepcopy
        res = deepcopy(res)
        if any(len(p.get('molecules', [])) > MAX_ROWS for p in res.get('panels', [])):
            raise ValueError('Protein architecture would exceed the ISV row limit; explicitly select the display read cap first')
    if hide_bleeding_reads and not res.get('read_filter', {}).get('enabled'):
        from .read_filters import filter_bleeding_reads
        res = filter_bleeding_reads(res, tolerance=bleed_tolerance)
    if res.get('read_filter', {}).get('enabled'):
        f = res['read_filter']
        note = (f"Display filter: hide reads extending into assigned-model introns; {f['tolerance_nt']} nt tolerance. "
                f"{f['eligible_RNA']:g} eligible of {f['input_RNA']:g} examined RNA records before sampling; "
                f"{f['hidden_RNA']:g} hidden. Original assignments/count matrices unchanged.")
        if res.get('display_read_cap_per_isoform'):
            note += f" Up to {res['display_read_cap_per_isoform']} eligible RNA per isoform are sampled after filtering."
        if not protein_features:
            caption = (str(caption) + ' ' if caption else '') + note
    # Overlapping annotations can share an axis without representing trans-splicing.
    # Give each reference model its own row so its title/exons remain readable.
    stack_reference_loci = is_trans_splicing_view(res) or bool(res.get('has_multiple_loci'))
    locus_frames = []
    if res.get("has_multiple_loci"):
        if region:
            raise ValueError("Multi-locus ISV requires separate locus windows; a single --region is ambiguous")
        res, locus_frames = prepare_locus_coordinates(res, intron_scale, trim_to_reads, isoform_track)
        intron_scale = 1.0  # each locus was compressed independently, once
    panels = [p for p in res.get("panels", []) if p.get("n_molecules", 0) > 0]
    gm = res.get("gene_model") or {}
    if not panels:
        print("No molecules for this selection -> no PDF written.")
        return False

    _region = parse_region(region) if isinstance(region, str) else region
    for p in panels:
        for m in p["molecules"]:
            m["merged"] = merge_locus_segments(m.get("exon_segments")) if locus_frames else merge_segments(m.get("exon_segments"))
        if _region:
            _kept = [m for m in p["molecules"] if mol_in_window(m["merged"], _region[0], _region[1])]
            print("[isv_web_remote] region {:,}-{:,}: {} keeps {} of {} molecules".format(
                _region[0], _region[1], p.get("condition"), len(_kept), len(p["molecules"])), flush=True)
            p["molecules"] = _kept
            p["n_molecules"] = len(_kept)
        for m in p["molecules"]:
            m["_span"] = mol_extent(m["merged"])
    extra_models = res.get("reference_isoforms", [])
    for m in extra_models:
        m["merged"] = merge_locus_segments(m.get("exon_segments")) if locus_frames else merge_segments(m.get("exon_segments"))
    # Isoform groups run LONGEST ISOFORM FIRST by default. One rank over the molecules of every panel,
    # so each panel groups in the same order and the structure track below them matches.
    # iso_order="given" keeps the caller's own order, which main() passes for --isoforms.
    if iso_order == "length":
        _rank = iso_order_map([m for p in panels for m in p["molecules"]])
        for p in panels:
            for m in p["molecules"]:
                m["_iso_rank"] = _rank[str(m.get("final_isoform_id") or m.get("isoform_id"))]
        print("[isv_web_remote] isoform order: longest isoform first over {} isoform(s)".format(len(_rank)),
              flush=True)
    elif iso_order != "given":
        raise ValueError("iso_order must be 'length' or 'given', got {!r}".format(iso_order))
    for p in panels:
        p["mol"] = sorted(p["molecules"], key=row_sort_key)
        if len(p["mol"]) > MAX_ROWS:
            p["mol"] = p["mol"][:MAX_ROWS]
    panels = [p for p in panels if p["mol"]]
    if not panels:
        print("No molecules overlap --region -> no PDF written.")
        return False

    all_mol = [m for p in panels for m in p["mol"]]
    FS = max(0.1, float(font_scale))
    # the left margin holds the isoform labels; it is needed only while something is labelled there
    _need_left = bool(show_read_labels or (isoform_track and not structure_labels_above))
    _chrom = str(gm.get("chrom") or "")
    # with no labels in the margin the figure preset still writes the chromosome name to the LEFT of the
    # gene model (ha="right" at track_x0 - 6, 9*FS pt). Reserve its rendered width, or it is clipped:
    # "chr2" at FS 2.6 is 56pt wide against a 39pt margin, which printed as "hr2".
    _chrom_w = (0.62 * 9.0 * FS * len(_chrom) + 12.0) if (structure_labels_above and _chrom) else 0.0
    # The margin has to hold the longest label drawn in it. Adding the protein length to every name
    # pushed "52751007.WF86-MDS-0609-6020_CD34 (635aa)" past the page edge at the fixed 152 pt.
    _lbl_w = 0.0
    if _need_left:
        for _m in all_mol:
            _lbl_w = max(_lbl_w,
                         0.60 * 8.0 * FS * len(isoform_label(_m, label_protein_length, id_max=26)))
        if isoform_track and not structure_labels_above:
            for _m in all_mol:
                _lbl_w = max(_lbl_w,
                             0.60 * 7.0 * FS * len(isoform_label(_m, label_protein_length,
                                                                  style="bar", id_max=32)))
    LEFT_ = (min(max(LEFT * FS, _lbl_w + 12.0), width * 0.45) if _need_left
             else max(28.0, 15.0 * FS, _chrom_w))
    HDR_ = HDR * FS                                   # panel header band grows with its text
    RULER_H_ = (RULER_H * FS) if show_ruler else (20 * FS if locus_frames else 0.0)
    EXON_LBL_H = 22 * FS                              # row under the gene model holding the exon labels
    # title block above the first panel: italic gene symbol, then the splicing event. Defined here because
    # the auto-fit budget below subtracts it before any panel geometry is laid out.
    _t_gene_h = (15.0 * FS) if title_gene else 0.0
    _t_event_h = (12.0 * FS) if title_event else 0.0
    TITLE_H = (_t_gene_h + _t_event_h + 6.0 * FS) if (title_gene or title_event) else 0.0
    # trim_to_reads: ignore the gene model extent so the axis spans only what the reads cover.
    # region: the caller fixed the window, so neither the reads nor the gene model may widen it.
    lo, hi = _region if _region else extent_of(all_mol, {} if trim_to_reads else gm)
    if locus_frames:
        lo, hi = min(f["start"] for f in locus_frames), max(f["end"] for f in locus_frames)
    track_x0, track_x1 = LEFT_, width - RIGHT - RMARGIN

    # Use the established two-/three-isoform schemes consistently across read and structure tracks.
    # Explicit caller colors still take precedence.
    _ckeys = sorted({color_key(m) for m in all_mol + extra_models})
    _palette = {2: TWO_ISOFORM_COLORS, 3: THREE_ISOFORM_COLORS}.get(len(_ckeys), [])
    _auto_colors = dict(zip(_ckeys, _palette))

    def colof(m):
        return m.get("display_color") or _auto_colors.get(color_key(m)) or color_for(m)

    # unique plotted isoforms (one structure row each) for the optional bottom structure track.
    #
    # The row's STRUCTURE is the isoform's REGISTERED structure from gff-output/combined.gff.gz,
    # attached as m["registered_exons"] by main() via precompute_gff_structures.StructureIndex.
    # It is NOT the representative read's structure. That read is arbitrary with respect to
    # completeness: when every read carries count 1 the `-count` tie-break below is constant and
    # the choice degenerates to lexicographic molecule-id order. On FLNB that drew
    # ENST00000295956 as 5 merged exons / 11,161 bp (58,159,623-58,170,783) in place of its
    # registered 46 exons / 163,830 bp (58,008,422-58,172,251) -- 151,201 bp of 5' structure lost.
    # A representative read is used ONLY when the isoform is not registered, or the index is absent.
    iso_rows = []
    n_registered = n_fallback = 0
    if isoform_track:
        seen = set()
        for m in sorted(all_mol + extra_models, key=lambda mm: (group_key(mm),
                                                 str(mm.get("final_isoform_id") or mm.get("isoform_id")))):
            key = m.get("final_isoform_id") or m.get("isoform_id")
            if key in seen or not (m.get("merged") or {}).get("exons"):
                continue
            seen.add(key)
            reg = m.get("registered_exons")
            if m.get("_registered_locus_merged"):
                m["struct_merged"] = m["_registered_locus_merged"]
                m["struct_source"] = "registered"
                n_registered += 1
            elif reg:
                m["struct_merged"] = {"exons": [{"a": int(a), "b": int(b), "label": ""} for a, b in reg],
                                      "introns": []}
                m["struct_source"] = "registered"
                n_registered += 1
            else:
                m["struct_merged"] = m["merged"]
                m["struct_source"] = "representative-read"
                n_fallback += 1
            iso_rows.append(m)
        print("[isv_web_remote] structure track: {} isoform(s) from the registered combined.gff.gz "
              "structures, {} fallback to a representative read".format(n_registered, n_fallback),
              flush=True)
        print("[isv_web_remote]   drawn order ({}): ".format(iso_order)
              + ", ".join("{} {:,}bp".format(str(m.get("final_isoform_id") or m.get("isoform_id")),
                                             iso_struct_extent(m)[1] - iso_struct_extent(m)[0])
                          for m in iso_rows), flush=True)
        # a registered structure can reach beyond the reads' own extent -> widen so nothing is clipped.
        # Skipped under trim_to_reads, whose whole point is an axis set by the reads, and under region,
        # whose whole point is an axis set by the caller.
        if not trim_to_reads and not _region:
            for m in iso_rows:
                for e in m["struct_merged"]["exons"]:
                    lo = min(lo, e["a"]); hi = max(hi, e["b"])

    # intron-compressed x-map: exon intervals come from the gene-model blocks + the molecules' merged
    # exons; everything NOT covered by an exon is treated as intron and scaled by intron_scale.
    exon_ivs = []
    for b in gm.get("blocks", []):
        if b.get("start") is not None and b.get("end") is not None:
            exon_ivs.append((min(b["start"], b["end"]), max(b["start"], b["end"])))
    for m in all_mol:
        for e in m["merged"]["exons"]:
            exon_ivs.append((min(e["a"], e["b"]), max(e["a"], e["b"])))
        # A RETAINED intron is sequence the read carries, so it must keep true scale like an exon.
        # Only the gaps BETWEEN a molecule's segments are spliced out, and those stay compressed.
        # Before this, intron_scale shrank retained introns too: on the RUNX1 novel bin 3,538 of 3,551
        # reads carry one (median 1,229 bp against 238 bp of exon), so at --intron_scale 0.1 each read
        # was drawn 4.1x shorter than its sequence, while the reference bin (9 of 2,143 with a retained
        # intron) was drawn at full length. One bin shrank, the other did not.
        for e in m["merged"]["introns"]:
            exon_ivs.append((min(e["a"], e["b"]), max(e["a"], e["b"])))
    for m in iso_rows:                        # registered structure exons share the same x-map
        for e in m["struct_merged"]["exons"]:
            exon_ivs.append((min(e["a"], e["b"]), max(e["a"], e["b"])))
    xz = build_xz(_merge_ivs(exon_ivs), lo, hi, max(0.01, float(intron_scale)), track_x0, track_x1)

    for p in panels:
        raw = [min(ROW_CAP, max(ROW_MIN, m.get("count", 0) * READ_PX)) for m in p["mol"]]
        raw_sum = sum(raw) or 1
        n_clusters = len({m.get("cluster_index", 0) for m in p["mol"]})
        scale = min(MAXPANEL, max(raw_sum, n_clusters * 13)) / raw_sum
        # A sparse panel used to collapse into one band. The floor above counts COLOUR GROUPS, not
        # rows, so a 5-read single-isoform panel was scaled to 13 px andeach read drew 2.6 px tall,
        # which reads as a single block. When the rows would be thinner than ROW_VIS and giving every
        # row ROW_VIS still fits inside MAXPANEL, the panel is scaled to n_rows * ROW_VIS instead.
        # Dense panels are untouched: 100 reads would need 800 px, more than MAXPANEL, so the branch
        # does not fire and their geometry is unchanged.
        if p["mol"] and min(raw) * scale < ROW_VIS and len(p["mol"]) * ROW_VIS <= MAXPANEL:
            scale = (len(p["mol"]) * ROW_VIS) / raw_sum
        h, rows = 0.0, []
        for i, m in enumerate(p["mol"]):
            hh = max(0.55, raw[i] * scale)
            rows.append({"m": m, "dy": h, "h": hh})
            h += hh
        p["_rows"] = rows
        p["_bodyH"] = max(18, h, n_clusters * 10 * FS if show_read_labels else 0)
        p["_label_minH"] = n_clusters * 10 * FS if show_read_labels else 0

    # Auto-fit (no knob): stacking every cell-state panel on one page makes it taller than a screen, so
    # viewers cut off the lower panels. If the figure would exceed one screenful, shrink ALL panel bodies
    # + their rows uniformly so the whole thing fits, keeping the shared genomic x-axis. Only shrinks;
    # figures that already fit are untouched.
    _FIT_H = 900.0
    reference_height = (RULER_H_ + REF_H + EXON_LBL_H) * (max(1, len(locus_frames)) if stack_reference_loci else 1)
    _ov = 10 + TITLE_H + 2 + reference_height + len(panels) * (HDR_ + PAD + PANEL_GAP)
    _budget = _FIT_H - _ov
    _body_total = sum(p["_bodyH"] for p in panels)
    # Do not compress a dense panel below the height needed by all labels.
    _label_total = sum(p["_label_minH"] for p in panels)
    if _body_total > _budget > 0 and _label_total <= _budget:
        _f = _budget / _body_total
        for p in panels:
            p["_bodyH"] = max(p["_bodyH"] * _f, p["_label_minH"])
            for r in p["_rows"]:
                r["dy"] *= _f
                r["h"] *= _f

    y = 10.0 + TITLE_H
    for p in panels:
        p["_top"] = y
        p["_y"] = y + HDR_
        y += HDR_ + p["_bodyH"] + PAD + PANEL_GAP
    axis_top = y + 2
    height = axis_top + reference_height
    ruler_base = axis_top + RULER_H_ - 4
    ref_top = axis_top + RULER_H_
    struct_top = height + 6                              # bottom isoform-structure track starts here
    if iso_rows:
        _, _struct_h = struct_rows_geometry(len(iso_rows), FS, structure_labels_above)
        height = struct_top + 16 * FS + _struct_h + 8

    protein_layout = None
    if protein_features:
        from ...iso2function.plotting.isv_architecture import load_architectures, select_tracks, layout_tracks
        bundle = load_architectures(protein_features)
        ids = [str(m.get('final_isoform_id') or m.get('isoform_id')) for m in iso_rows]
        tracks = select_tracks(bundle, ids, res.get('symbol'))
        length_audit = []
        for m, track in zip(iso_rows, tracks):
            if m.get('struct_source') != 'registered':
                raise ValueError('Protein architecture requires the registered native transcript model: ' + track['isoform_id'])
            length_audit.append(_architecture_structure_length(m, track))
            # res is a drawing-only deepcopy. Frozen RNA input metadata/counts stay intact.
            m['protein_length'] = int(track['protein_aa'])
            track['display_color'] = colof(m)
        protein_top = height + 8 * FS
        protein_layout = layout_tracks(tracks, track_x1 - track_x0, FS)
        protein_layout['source'] = {k: bundle[k] for k in ('root', 'manifest_sha256', 'source_manifest_sha256', 'selection')}
        protein_layout['structure_length_display_audit'] = length_audit
        protein_layout['RNA'] = dict(panels=[dict(condition=p['condition'],
            records=len(p['mol']), reads=sum(float(m.get('count', 0)) for m in p['mol'])) for p in panels],
            isoform_order=ids, read_filter=res.get('read_filter'),
            display_read_cap_per_isoform=res.get('display_read_cap_per_isoform'))
        height = protein_top + protein_layout['height']

    caption_top = height
    caption_text = None
    if caption:
        import textwrap
        caption_text = textwrap.fill(str(caption), width=max(40, int((width - 40) / (3.8 * FS))))
        height += (len(caption_text.splitlines()) * 9.1 + 20) * FS
    fig = plt.figure(figsize=(width / 72.0, height / 72.0))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.axis("off")

    def rect(x, yy, w, hh, fc, ec="none", lw=0, alpha=1.0):
        ax.add_patch(Rectangle((x, yy), w, hh, facecolor=fc, edgecolor=ec, linewidth=lw, alpha=alpha))

    def hline(x0, x1, yy, color, lw=1.0, alpha=1.0):
        ax.plot([x0, x1], [yy, yy], color=color, linewidth=lw, alpha=alpha, solid_capstyle="butt")

    def backbone(segments, yy, color, joins=None, dashed_joins=False, read_height=None):
        if not locus_frames:
            hline(xz(min(e["a"] for e in segments)), xz(max(e["b"] for e in segments)), yy, color, lw=0.7, alpha=0.55)
            return
        groups = OrderedDict()
        # Thin connectors must stay thinner than densely packed read rows;
        # otherwise neighboring hairlines overlap into apparently solid RNA.
        line_width = min(0.7, max(0.02, read_height * 0.35)) if read_height is not None else 0.7
        for e in segments:
            groups.setdefault(e["gene"], []).append(e)
        for parts in groups.values():
            hline(xz(min(e["a"] for e in parts)), xz(max(e["b"] for e in parts)), yy, color, lw=line_width, alpha=0.55)
        # Stacked, phase-aligned dashes resemble vertical exon bars. Use native
        # solid hairlines in read pileups; reserve dashes for the model row.
        for a, b in joins or []:
            ax.plot([xz(a), xz(b)], [yy, yy],
                    color=color, linewidth=line_width, alpha=0.55,
                    linestyle=(0, (3, 3)) if dashed_joins else "-", solid_capstyle="butt")

    _tx = (track_x0 + track_x1) / 2.0
    if title_gene:
        ax.text(_tx, 10.0 + _t_gene_h * 0.72, str(title_gene), fontsize=13 * FS, style="italic",
                fontfamily=ITALIC_FAMILY, ha="center", va="center", color="#1a2230")
    if title_event:
        ax.text(_tx, 10.0 + _t_gene_h + _t_event_h * 0.62, str(title_event), fontsize=8.5 * FS,
                ha="center", va="center", color="#39434f")

    for p in panels:
        top = p["_top"]
        box_h = HDR_ + p["_bodyH"] + PAD
        rect(6, top, width - 12, box_h, "#ffffff", ec="#d4d9e0", lw=0.8)
        rect(6, top, width - 12, HDR_, "#eef1f5")
        if structure_labels_above:        # figure preset: italic, centred in the header band
            ax.text((track_x0 + track_x1) / 2.0, top + 15 * FS, f"{p['condition']}", fontsize=11 * FS,
                    weight="normal" if protein_layout else "bold", style="normal" if protein_layout else "italic",
                    fontfamily=ITALIC_FAMILY, ha="center", va="center",
                    color="#1a2230")
        else:
            ax.text(16, top + 15 * FS, f"{p['condition']}", fontsize=11 * FS, weight="bold", va="center", color="#1a2230")
        if show_read_counts:
            ax.text(16 + 8.5 * FS * len(str(p["condition"])) + 10, top + 15 * FS,
                    f"· {p['n_molecules']:,} reads", fontsize=9.5 * FS, va="center", color="#5a6675")
        cl_extent = {}
        for r in p["_rows"]:
            ci = r["m"].get("cluster_index", 0)
            if ci not in cl_extent:
                cl_extent[ci] = [r["dy"], r["dy"] + r["h"], r["m"]]
            else:
                cl_extent[ci][0] = min(cl_extent[ci][0], r["dy"])
                cl_extent[ci][1] = max(cl_extent[ci][1], r["dy"] + r["h"])
        if show_read_labels:
            groups = sorted(cl_extent.values(), key=lambda item: (item[0] + item[1]) / 2)
            centers = [(yt + yb) / 2 for yt, yb, m in groups]
            label_y = read_label_positions(centers, p["_bodyH"], 10 * FS)
            for (yt, yb, m), ly in zip(groups, label_y):
                center = (yt + yb) / 2
                if abs(ly - center) > 0.1:
                    ax.plot([LEFT_ - 7, track_x0 - 2], [p["_y"] + ly, p["_y"] + center],
                            color="#b8c0ca", linewidth=0.45, zorder=1)
                ax.text(LEFT_ - 8, p["_y"] + ly,
                        isoform_label(m, label_protein_length, id_max=26),
                        fontsize=8 * FS, ha="right", va="center", color=(KNOWN if m.get("known") else NOVEL))
        for r in p["_rows"]:
            m = r["m"]
            y_top = p["_y"] + r["dy"]
            color = colof(m)
            merged = m["merged"]
            exs = segs_in_window(merged["exons"], lo, hi)        # no-op without --region
            ins = segs_in_window(merged["introns"], lo, hi)
            eh = max(0.7, r["h"] - 0.8 if r["h"] > 3 else r["h"])
            ey = y_top + (r["h"] - eh) / 2
            # The backbone must span EVERY segment of the molecule, exons and retained introns alike.
            # Taking min/max over the exons alone left a retained-intron block that lies outside the exon
            # range with no line joining it: on LSP1 the read's I block sits at 1,889,516-1,889,609 while
            # its exons start at 1,891,773, so the block floated 2.1 kb from the C-terminal exon,
            # disconnected. A read is one molecule and must read as one row.
            _span = exs + ins
            if _span:
                backbone(_span, y_top + r["h"] / 2, color, m.get("_exon_segments_joins"), read_height=r["h"])
            # A RETAINED intron is sequence the molecule carries, so it is drawn as a solid block at
            # full height, exactly like an exon. It was previously drawn at 60% height and alpha 0.45,
            # which reads as "absent" and hid intron retention: on RUNX1 the gold bin IS intron-42
            # retention (3,538 of 3,551 reads) and it appeared as a faint bar beside solid exons.
            # Only the gaps BETWEEN segments are spliced out, and those stay as the backbone line.
            for s in ins + exs:
                xa, xb = xz(s["a"]), xz(s["b"])
                rect(xa, ey, max(1.0, xb - xa), eh, color)

    # ---- Native reference track, repeated for each partner locus on the SAME pileup axis. ----
    _full_lo, _full_hi, _full_gm = lo, hi, gm
    _axis_base, _ref_base, _ruler_base = axis_top, ref_top, ruler_base
    for _frame_index, _frame in enumerate(locus_frames or [dict(model=gm, start=lo, end=hi)]):
        _offset = _frame_index * (RULER_H_ + REF_H + EXON_LBL_H) if stack_reference_loci else 0
        axis_top, ref_top, ruler_base = _axis_base + _offset, _ref_base + _offset, _ruler_base + _offset
        gm, lo, hi = _frame["model"], _frame["start"], _frame["end"]
        _row_model, _row_lo, _row_hi = gm, lo, hi
        _windows = _frame.get("reference_windows", [_frame]) if stack_reference_loci else [_frame]
        for _window in _windows:
            gm, lo, hi = _window["model"], _window["start"], _window["end"]
            if show_ruler:
                hline(xz(lo) if locus_frames else track_x0, xz(hi) if locus_frames else track_x1, ruler_base, "#7a8696", lw=0.9)
                for tk in (_window["ticks"] if locus_frames else np.linspace(lo, hi, max(4, round((track_x1 - track_x0) / (120 * FS))))):
                    px = xz(_window["project"](tk)) if locus_frames else xz(tk)
                    ax.plot([px, px], [ruler_base - 5, ruler_base], color="#7a8696", linewidth=0.8, solid_capstyle="butt")
                    ax.text(px, ruler_base - 8, fmt(tk), fontsize=8 * FS, ha="center", va="bottom", color="#5a6675")
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
            _lbl_right = -1e9                                  # right edge of the last exon label, in points
            # The default keeps the historical width threshold. dense_exon_labels instead draws one row and skips a
            # label only when its text would touch the previous one, so the exons inside the plotted region are the
            # ones named and nothing overlaps.
            _min_w = 1.0 if dense_exon_labels else 7 * FS
            _efs = (8 * FS) if dense_exon_labels else (7 * FS)
            for b in gm.get("blocks", []):
                # Only under --region. The pixel test below cannot catch a block outside the axis, because xz clamps
                # it INTO the axis, so it lands exactly on track_x0/track_x1 and is drawn as a minimum-width tick.
                # That artifact also exists under trim_to_reads (3 blocks on the ATP5J and LSP1 panels), but removing
                # it there would change every panel already rendered, so it is reported and left alone, not fixed here.
                if _region:
                    _ba, _bb = min(b["start"], b["end"]), max(b["start"], b["end"])
                    if _bb < lo or _ba > hi:
                        continue
                xa, xb = xz(b["start"]), xz(b["end"])
                if xb < track_x0 or xa > track_x1:
                    continue
                rect(xa, ref_top + 3, max(1.6, xb - xa), REF_H - 6, REF_FILL)
                lbl = str(b.get("exon_id") or "").split(".")[0]
                if lbl and lbl != last_lbl and (xb - xa) > _min_w:
                    if dense_exon_labels:
                        _half = 0.30 * _efs * len(lbl)          # half the rendered text width
                        _cx = (xa + xb) / 2
                        if _cx - _half < _lbl_right + 0.25 * _efs:
                            continue                            # would touch the previous label -> skip this one
                        _lbl_right = _cx + _half
                    ax.text((xa + xb) / 2, ref_top + REF_H + 10, lbl, fontsize=_efs, ha="center", va="top", color="#5a6675")
                    last_lbl = lbl
            if show_ruler and not locus_frames:
                ax.text(track_x0, ref_top + REF_H + 11, gm.get("chrom") or "", fontsize=8 * FS, ha="left", va="top", color="#5a6675")
            elif structure_labels_above and not locus_frames:
                ax.text(track_x0 - 6, ref_top + REF_H / 2.0, gm.get("chrom") or "", fontsize=9 * FS,
                        ha="right", va="center", color="#1a2230")
        gm, lo, hi = _row_model, _row_lo, _row_hi

        if locus_frames:
            ax.text((xz(lo) + xz(hi)) / 2, axis_top + 9 * FS,
                    "{} · {} ({})".format(gm.get("symbol") or gm.get("gene"), gm.get("chrom"), gm.get("strand")),
                    fontsize=8 * FS, ha="center", va="center", color="#1a2230")
    lo, hi, gm = _full_lo, _full_hi, _full_gm
    # ---- PDF-only bottom isoform-structure track (one row per plotted isoform; UCSC CDS-tall / UTR-thin,
    #      introns as connecting lines), aligned to the same x-axis as the molecule pileup above. ----
    if iso_rows:
        if show_structure_header:
            ax.text(LEFT_ - 8, struct_top + 9 * FS, "Isoform structures", fontsize=8.5 * FS, weight="bold",
                    ha="right", va="center", color="#1a2230")
        ystruct = struct_top + (16 if show_structure_header else 4) * FS
        ISO_EXH_ = ISO_EXH * FS
        _geo, _ = struct_rows_geometry(len(iso_rows), FS, structure_labels_above, y0=ystruct)
        for _idx, m in enumerate(iso_rows):
            color = colof(m)
            merged = m.get("struct_merged") or m["merged"]
            exs, ins = merged["exons"], merged["introns"]
            rc, _ly, _va = _geo[_idx]                             # row centre, label y, label anchor
            ex_y = rc - ISO_EXH_ / 2.0                            # tall (CDS) exon top
            utr_h = max(3.0, ISO_EXH_ * 0.5)
            utr_y = ex_y + (ISO_EXH_ - utr_h) / 2.0               # thin (UTR) exon top
            exs = [e for e in exs if not (e["b"] < lo or e["a"] > hi)]   # drop segments off a trimmed axis
            if exs:                                              # backbone line across the isoform span (introns)
                backbone(exs, rc, color, m.get("_registered_exon_segments_joins"), dashed_joins=True)
            edge = darken(color)                                  # exon outline = same-hue darker shade of fill
            cds_lo, cds_hi = m.get("cds_min"), m.get("cds_max")
            locus_cds = m.get("_locus_cds_intervals")
            for s in exs:
                xa, xb = xz(s["a"]), xz(s["b"])
                if locus_cds is not None:
                    rect(xa, utr_y, max(0.8, xb - xa), utr_h, color, ec=edge, lw=0.5)
                    for cds in locus_cds:
                        if cds["gene"] != s["gene"]:
                            continue
                        ca, cb = max(s["a"], cds["a"]), min(s["b"], cds["b"])
                        if cb > ca:
                            rect(xz(ca), ex_y, max(1.0, xz(cb) - xz(ca)), ISO_EXH_, color, ec=edge, lw=0.5)
                elif cds_lo is not None and cds_hi is not None:
                    rect(xa, utr_y, max(0.8, xb - xa), utr_h, color, ec=edge, lw=0.5)     # UTR (thin) full exon
                    ca, cb = max(s["a"], cds_lo), min(s["b"], cds_hi)
                    if cb > ca:
                        cxa, cxb = xz(ca), xz(cb)
                        rect(cxa, ex_y, max(1.0, cxb - cxa), ISO_EXH_, color, ec=edge, lw=0.5)  # CDS (tall) coding overlap
                else:
                    rect(xa, ex_y, max(1.0, xb - xa), ISO_EXH_, color, ec=edge, lw=0.5)   # no CDS info -> full height
            _iid = str(m.get("structure_label") or m.get("final_isoform_id") or m.get("isoform_id"))
            _aa = m.get("protein_length")
            _nmd = str(m.get("nmd_status") or "").strip() in ("NMD", "Potential-NMD")
            if structure_labels_above:
                # "<id>  >  <n>aa" on its own line, in this structure's colour. The last row carries its
                # label UNDER the structure so the two structures of an event sit next to each other.
                _lfs = 9.5 * FS
                ax.text(track_x0, _ly, _iid, fontsize=_lfs, ha="left", va=_va, color=color,
                        weight="normal" if protein_layout else "bold")
                if _aa not in (None, ""):
                    _w = 0.60 * _lfs * len(_iid)                      # rendered width of the id
                    _tx, _ty = track_x0 + _w + 0.45 * _lfs, _ly + (0.38 * _lfs if _va == "top" else -0.38 * _lfs)
                    _h = 0.34 * _lfs
                    ax.add_patch(Polygon([[_tx, _ty - _h], [_tx, _ty + _h], [_tx + 1.5 * _h, _ty]],
                                         closed=True, facecolor=color, edgecolor="none", zorder=45))
                    ax.text(_tx + 2.3 * _h, _ly, f"{int(_aa)}aa" + (" NMD" if _nmd else ""),
                            fontsize=_lfs, ha="left", va=_va, color=color, weight="normal" if protein_layout else "bold")
            else:
                _lab = isoform_label(m, label_protein_length, style="bar",
                                     id_max=32 if label_protein_length else 22)
                ax.text(LEFT_ - 8, rc, _lab,
                        fontsize=7 * FS, ha="right", va="center", color=(KNOWN if m.get("known") else NOVEL))

    if protein_layout:
        from ...iso2function.plotting.isv_architecture import draw_matplotlib
        protein_artists = draw_matplotlib(ax, protein_layout, track_x0, protein_top)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        labels = [artist for artist in protein_artists if hasattr(artist, 'get_text')]
        for artist in labels:
            box = artist.get_window_extent(renderer)
            if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.x1+1 or box.y1 > fig.bbox.y1+1:
                raise ValueError('Protein architecture label outside page: ' + artist.get_text())
        from itertools import combinations
        for first, second in combinations(labels, 2):
            if first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer)):
                raise ValueError('Protein architecture labels overlap: ' + first.get_text() + ' / ' + second.get_text())
    if caption_text:
        ax.text(20, caption_top + 8 * FS, caption_text, fontsize=7 * FS,
                ha="left", va="top", color="#39434f", linespacing=1.3)
    if protein_layout:
        from datetime import datetime, timezone
        epoch = datetime(2000, 1, 1, tzinfo=timezone.utc)
        fig.savefig(out_path, format="pdf", metadata=dict(CreationDate=epoch, ModDate=epoch,
                    Creator='AltAnalyze3 ISV protein architecture'))
        import json, hashlib
        from pathlib import Path
        protein_layout['pdf_sha256'] = hashlib.sha256(Path(out_path).read_bytes()).hexdigest()
        Path(out_path).with_suffix('.protein_tracks.json').write_text(json.dumps(protein_layout, indent=2, allow_nan=False)+'\n')
    else:
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
               "exon_segments", "protein_length", "nmd_status", "annotation_tokens",
               "cds_min", "cds_max", "cds_chrom",
               "gene_path", "locus_genes", "locus_status", "locus_joins", "trans_splicing_candidate")


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
            "gene_model": res.get("gene_model"), "panels": panels,
            "gene_models": res.get("gene_models"), "has_multiple_loci": res.get("has_multiple_loci", False)}


def keep_top_isoforms(res, n):
    """Keep only the n final isoforms with the most reads pooled across ALL panels (ties: isoform id), so every
    panel shows the same isoforms and colors (n=2 -> the renderer's gold/blue two-isoform palette)."""
    panels = res.get("panels", [])
    by_iso = {}
    for p in panels:
        for m in p.get("molecules", []):
            k = m.get("final_isoform_id") or m.get("isoform_id")
            by_iso[k] = by_iso.get(k, 0.0) + float(m.get("count", 0) or 0)
    keep = {k for k, _ in sorted(by_iso.items(), key=lambda kv: (-kv[1], str(kv[0])))[:int(n)]}
    for p in panels:
        p["molecules"] = [m for m in p.get("molecules", [])
                          if (m.get("final_isoform_id") or m.get("isoform_id")) in keep]
        p["n_molecules"] = int(round(sum(float(m.get("count", 0) or 0) for m in p["molecules"])))
    return res, sorted(keep, key=lambda k: -by_iso[k]), by_iso


def top_isoform_per_filter(da, ctx, gene, cell_states, conditions, patterns):
    """For EACH --isoform_filter pattern, the final isoform with the most reads pooled across the panels among
    the reads carrying that pattern (ties: isoform id). A final already picked for an earlier pattern is
    skipped, so a two-junction event yields one isoform per junction. --top_isoforms 2 can instead keep two
    isoforms that carry the SAME junction: FLI1 E3.4|E5.2 vs U1.1_128685372|E5.2 kept ENST00000527786 and
    ENST00000281428, both E3.4|E5.2. -> [(pattern, final_id, reads, n_finals_carrying_pattern)]."""
    picks, chosen = [], set()
    for pat in patterns:
        r = da.query_reads(ctx, gene, cell_states, conditions=conditions, max_isoforms=10**9,
                           panel_by="cell_type", isoform_filter=[pat])
        by_iso = {}
        for p in r.get("panels", []):
            for m in p.get("molecules", []):
                k = m.get("final_isoform_id") or m.get("isoform_id")
                by_iso[k] = by_iso.get(k, 0.0) + float(m.get("count", 0) or 0)
        ranked = [kv for kv in sorted(by_iso.items(), key=lambda kv: (-kv[1], str(kv[0]))) if kv[0] not in chosen]
        if not ranked:
            raise ValueError(f"--top_per_filter: no final isoform carries {pat} apart from {sorted(chosen)}")
        k, v = ranked[0]
        chosen.add(k)
        picks.append((pat, k, v, len(by_iso)))
    return picks


def keep_isoforms(res, keep):
    """Keep only molecules whose final isoform is in `keep`, in every panel."""
    for p in res.get("panels", []):
        p["molecules"] = [m for m in p.get("molecules", [])
                          if (m.get("final_isoform_id") or m.get("isoform_id")) in keep]
        p["n_molecules"] = int(round(sum(float(m.get("count", 0) or 0) for m in p["molecules"])))
    return res


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
    ap.add_argument("--region", default=None,
                    help="draw ONLY this genomic window: START-END or chr:START-END (commas allowed). The "
                         "x-axis becomes the window, reads that do not overlap it are dropped, and segments "
                         "outside it are not drawn. Use it when the event sits in a small part of a long "
                         "gene: the default axis is the union of every plotted molecule, so a 209 bp splice "
                         "inside one exon of a 261 kb gene is drawn under a point wide. Overrides the "
                         "read-set axis of --figure_mode.")
    ap.add_argument("--out", default=None, help="output PDF (default: <run_dir>/ISV_<gene>_<states>_reads.pdf)")
    ap.add_argument('--protein-features', default=None,
                    help='Append exon-colored protein/domain tracks from completed iso2function protein-impact-plots output with a coding-exon map')
    ap.add_argument('--render-input', default=None,
                    help='Render a saved ISV response JSON(.gz) without loading native matrices or sequence files; query/filter options cannot be combined')
    ap.add_argument('--render-input-sha256', default=None,
                    help='Require this SHA256 for the saved response before rendering')
    ap.add_argument('--caption-file', default=None,
                    help='UTF-8 scientific figure caption appended below the exported PDF')
    ap.add_argument('--save-render-input', default=None,
                    help='Save the exact selected response as JSON(.gz) for indexed, reproducible rerendering')
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
    ap.add_argument("--isoform_filter", nargs="+", default=None,
                    help="keep only reads whose structure contains ANY of these exon/junction patterns, e.g. "
                         "'E1.2|E3.1' 'E2.1|E3.1' (isoform_structure_view isoform_filter syntax: '|' joins the two exons)")
    ap.add_argument("--max_reads", default=None,
                    help="cap every panel at this many reads (after --isoform_filter) by sampling reads without "
                         "replacement (seeded), so each isoform keeps its expected share; 'auto' = the smallest "
                         "panel's total, so all panels show the same number of reads")
    ap.add_argument("--seed", type=int, default=0, help="seed for --max_reads sampling (default 0)")
    ap.add_argument("--hide-bleeding-reads", "--hide_bleeding_reads", action="store_true",
                    help="Display only: exclude reads extending across a model splice edge into its intron; filter the full pool before sampling")
    ap.add_argument("--bleed-tolerance", "--bleed_tolerance", type=int, default=0,
                    help="Allowed extension into model introns in nt (default 0: strict)")
    ap.add_argument("--max-reads-per-isoform", "--max_reads_per_isoform", type=int, default=None,
                    help="Sample this many eligible reads per isoform after filtering (default 200 with --hide-bleeding-reads)")
    ap.add_argument("--top_isoforms", type=int, default=None,
                    help="keep only the N final isoforms with the most reads across all panels (applied after "
                         "--isoform_filter, before --max_reads); N=2 draws them in the gold/blue two-isoform palette")
    ap.add_argument("--isoforms", nargs="+", default=None,
                    help="draw ONLY these final isoform ids, named exactly (e.g. ENST00000300305 "
                         "ENSG00000159216:28775276.EP_SRSF2_Ctrl). Applied after --isoform_filter and before "
                         "--max_reads; exits if any named id carries no read in the selection. Two ids draw in "
                         "the gold/blue two-isoform palette, in the order given")
    ap.add_argument("--no-protein-length", "--no_protein_length", dest="no_protein_length",
                    action="store_true",
                    help="write bare isoform ids. By default every isoform name carries its protein "
                         "length, as '<id> (568aa)', and ', NMD' inside the brackets when the "
                         "isoform is NMD.")
    ap.add_argument("--figure-mode", "--figure_mode", dest="figure_mode", action="store_true",
                    help="multi-gene figure preset: large fonts (--font_scale, default 1.9) and a panel stripped "
                         "for a small reproduction -- no coordinate ruler or chromosome label, no per-read isoform "
                         "labels inside the pileups, and an x-axis set by the plotted reads so gene regions with no "
                         "reads are not drawn. The isoform-structure track is kept and still names each isoform.")
    ap.add_argument("--title", default=None,
                    help="gene line of the title, drawn in italic above the panels. Default: the gene symbol "
                         "under --figure_mode, nothing otherwise. Pass an empty string to omit it.")
    ap.add_argument("--title_event", "--title-event", dest="title_event", default=None,
                    help="event line of the title, under the gene line. Default under --figure_mode: the "
                         "--isoform_filter junctions as <gene>:<junction A>|<junction B>.")
    ap.add_argument("--font_scale", "--font-scale", dest="font_scale", type=float, default=None,
                    help="multiply every font size, label margin and header height (default 1.0, or 1.9 with "
                         "--figure_mode). Use it alone to enlarge text without the other figure-mode changes.")
    ap.add_argument("--top_per_filter", action="store_true",
                    help="with --isoform_filter: keep the final isoform with the most reads carrying EACH filter "
                         "pattern, so a two-junction event draws one isoform per junction in the gold/blue palette. "
                         "Works with --cell_states and with --cell_groups (ranked over the union of member states); "
                         "applied before --max_reads; excludes --top_isoforms")
    ap.add_argument("--export_molecules", default=None,
                    help="ALSO write a TSV of (molecule_id, final_isoform_id, cell_state, sample, count) for "
                         "EVERY molecule in the selection (UNCAPPED, ignores --max_isoforms) to this path -- "
                         "for reconciliation against the raw reads.db index and the source pseudobulk counts.")
    return ap


def cap_panel_reads(res, max_reads, seed=0):
    """Downsample every panel to max_reads reads without replacement (multivariate hypergeometric over the
    molecules' integer read counts, seeded), keeping each molecule's expected share. 'auto' -> the smallest
    panel total. Panels at or below the cap are unchanged. Returns (res, cap)."""
    import numpy as np
    panels = res.get("panels", [])
    totals = [int(round(sum(float(m.get("count") or 0) for m in p.get("molecules", [])))) for p in panels]
    cap = min([t for t in totals if t > 0], default=0) if str(max_reads).lower() == "auto" else int(max_reads)
    if cap <= 0:
        raise ValueError(f"max_reads resolved to {cap}; panel totals {totals}")
    for k, p in enumerate(panels):
        mols = p.get("molecules", [])
        counts = np.array([int(round(float(m.get("count") or 0))) for m in mols], dtype=np.int64)
        if int(counts.sum()) <= cap:
            continue
        rng = np.random.default_rng([int(seed), k])
        drawn = rng.multivariate_hypergeometric(counts, cap)
        kept = []
        for m, n in zip(mols, drawn):
            if n > 0:
                m = dict(m); m["count"] = float(n); kept.append(m)
        if int(sum(m["count"] for m in kept)) != cap:
            raise AssertionError(f"cap_panel_reads drew {sum(m['count'] for m in kept)} reads, expected {cap}")
        p["molecules"] = kept
        p["n_molecules"] = cap
    return res, cap


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


def _render_selected(a, res, out, region):
    """One rendering path for indexed queries and saved native responses."""
    if a.save_render_input:
        import gzip, json
        opener = gzip.open if a.save_render_input.endswith('.gz') else open
        with opener(a.save_render_input, 'wt') as fh:
            json.dump(res, fh, allow_nan=False)
    fs = a.font_scale if a.font_scale is not None else (2.6 if a.figure_mode else 1.0)
    title = a.title if a.title is not None else ((res.get('symbol') or a.gene) if a.figure_mode else None)
    event = a.title_event
    if event is None and a.figure_mode and a.isoform_filter:
        event = str(res.get('gene') or a.gene) + ':' + '|'.join(str(x).replace('|', '-', 1) for x in a.isoform_filter)
    return render(res, out, width=a.width, intron_scale=a.intron_scale,
                  isoform_track=a.isoform_track, font_scale=fs,
                  show_ruler=not a.figure_mode, show_read_labels=not a.figure_mode,
                  trim_to_reads=a.figure_mode, title_gene=title or None, title_event=event or None,
                  show_read_counts=not a.figure_mode, show_structure_header=not a.figure_mode,
                  label_protein_length=not a.no_protein_length, dense_exon_labels=a.figure_mode,
                  structure_labels_above=a.figure_mode, region=region,
                  iso_order='given' if a.isoforms or a.render_input else 'length',
                  protein_features=a.protein_features,
                  caption=(Path(a.caption_file).read_text().strip()
                           if getattr(a, 'caption_file', None) else None))


def main(argv=None):
    return run_from_args(build_arg_parser().parse_args(argv))


def run_from_args(a):
    """Shared standalone and ``altanalyze3 visualize isv-domains`` CLI handler."""
    if a.render_input_sha256 and not a.render_input:
        sys.exit('--render-input-sha256 requires --render-input')
    if a.hide_bleeding_reads:
        if a.bleed_tolerance < 0:
            sys.exit('bleed_tolerance must be nonnegative')
        a.max_isoforms = 10**9  # refill from the full eligible pool, never the old display cap
    try:                                  # fail on a malformed window before the read query runs
        _region = parse_region(a.region)
    except ValueError as exc:
        sys.exit("[isv_web_remote] {}".format(exc))
    if a.render_input:
        import gzip, json
        conflicting = [name for name in ('cell_states', 'cell_groups', 'condition_groups', 'conditions',
                       'isoforms', 'isoform_filter', 'top_isoforms', 'top_per_filter', 'max_reads',
                       'max_reads_per_isoform', 'hide_bleeding_reads', 'min_isoform_pct', 'export_molecules')
                       if getattr(a, name)]
        if conflicting:
            sys.exit('--render-input is already selected; query/filter options cannot be applied: ' + ', '.join(conflicting))
        if a.max_isoforms != 1_000_000:
            sys.exit('--render-input is already selected; --max_isoforms cannot be applied')
        opener = gzip.open if a.render_input.endswith('.gz') else open
        if a.render_input_sha256:
            import hashlib
            from pathlib import Path
            if hashlib.sha256(Path(a.render_input).read_bytes()).hexdigest() != a.render_input_sha256:
                sys.exit('--render-input SHA256 mismatch; saved response retained unchanged')
        with opener(a.render_input, 'rt') as fh:
            res = json.load(fh)
        if a.gene not in (res.get('gene'), res.get('symbol')):
            sys.exit('--render-input gene does not match --gene')
        out = a.out or os.path.join(a.run_dir, 'ISV_' + a.gene + '_protein_architecture.pdf')
        sys.exit(0 if _render_selected(a, res, out, _region) else 2)
    from . import data_api as da
    # data_api imports isoform_structure_view, whose module level sets font.sans-serif to ['DejaVu Sans']
    # (isoform_structure_view.py line 25). That import runs AFTER the settings at the top of this file, so
    # every PDF written here carried DejaVu instead of Arial. Re-assert the fonts this module asks for, now
    # that the whole import chain is loaded. Nimbus Sans is the URW Helvetica clone and supplies the italic
    # face that the bundled Arial files lack.
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    metadata = a.metadata or os.path.join(a.run_dir, "mds_metadata_bam.txt")
    gene_model = a.gene_model or os.path.join(a.run_dir, "Hs_Ensembl_exon.txt")
    gene_symbol = a.gene_symbol or os.path.join(a.run_dir, "cellHarmony", "Hs_Ensembl-annotations.txt")
    conditions = [c.strip() for c in a.conditions.split(",") if c.strip()] if a.conditions else None
    if not (a.cell_states or a.cell_groups):
        sys.exit("ERROR: provide --cell_states (one panel per state) or --cell_groups (combined panels).")
    if a.top_per_filter and (not a.isoform_filter or a.top_isoforms):
        sys.exit("ERROR: --top_per_filter needs --isoform_filter and excludes --top_isoforms.")
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
    if a.figure_mode or not a.no_protein_length:
        # every isoform name carries its protein length. The lite loader above skips
        # protein_summary.txt, so without this the length is missing and the label degrades to the bare id.
        try:
            ctx._index_protein_summary()
        except Exception as _e:
            print(f"[isv_web_remote] WARN: protein summary unavailable ({_e}); labels carry the isoform id only",
                  flush=True)
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
        all_states = rank_states = sorted({s for v in cell_groups.values() for s in v})
        slug = "_".join(list(cell_groups) +
                        ([] if (len(cond_groups) == 1 and "all" in cond_groups) else list(cond_groups)))
        out = a.out or os.path.join(a.run_dir, f"ISV_{a.gene}_{slug}_reads.pdf")
        print(f"[isv_web_remote] gene={a.gene}->{gene} cell_groups={dict(cell_groups)} "
              f"condition_groups={ {k: (v or 'ALL') for k, v in cond_groups.items()} } "
              f"grouping={ctx.molecule_grouping}", flush=True)
        # ONE query over the union of all member states (huge cap so nothing is dropped before re-aggregation)
        # -> a single shared final-isoform clustering, then re-composed into the requested panels.
        # --isoform_filter must be applied HERE too: without it the composed panels kept every isoform of the
        # gene while the log still printed an "isoform_filter=..." line, so the flag was silently ignored in
        # --cell_groups mode (RUNX1 HSPC panel read 2,930 filtered-looking reads vs 1,139 actually filtered).
        base = da.query_reads(ctx, gene, all_states, conditions=None,
                              max_isoforms=(10**9 if a.hide_bleeding_reads else 1_000_000), panel_by="cell_type", isoform_filter=a.isoform_filter)
        res = compose_panels(ctx, base, cell_groups, cond_groups, a.max_isoforms)
    else:
        # ---- ONE panel per cell state (original behaviour) ----
        cell_states = rank_states = [s.strip() for s in a.cell_states.split(",") if s.strip()]
        out = a.out or os.path.join(a.run_dir, f"ISV_{a.gene}_{'_'.join(cell_states)}_reads.pdf")
        print(f"[isv_web_remote] gene={a.gene}->{gene} states={cell_states} "
              f"conditions={conditions or 'ALL'} grouping={ctx.molecule_grouping}", flush=True)
        res = da.query_reads(ctx, gene, cell_states, conditions=conditions,
                             max_isoforms=a.max_isoforms, panel_by="cell_type", isoform_filter=a.isoform_filter)

    _src = da.read_segment_stats()
    if _src:
        _tot = sum(_src.values())
        print("[isv_web_remote] read blocks: %d of %d from the sample coordinate store, %d parsed "
              "from the exon string" % (_src.get("store", 0), _tot, _src.get("parsed", 0)), flush=True)

    if a.export_molecules:
        # UNCAPPED re-query so EVERY molecule is exported (the render caps at max_isoforms; the
        # reconciliation TSV must not). One row per (molecule, cell_state, sample).
        exp_states = ([s.strip() for s in a.cell_states.split(",") if s.strip()] if a.cell_states
                      else sorted({s for v in parse_group_spec(a.cell_groups).values() for s in v}))
        full = da.query_reads(ctx, gene, exp_states, conditions=conditions,
                              max_isoforms=10**9, panel_by="cell_type", isoform_filter=a.isoform_filter)
        _export_molecule_tsv(full, a.export_molecules)

    if a.isoform_filter:
        print(f"[isv_web_remote] isoform_filter={a.isoform_filter}: panel reads "
              + ", ".join(f"{p['condition']}:{int(round(sum(m['count'] for m in p['molecules'])))}"
                          for p in res.get("panels", [])), flush=True)
    if a.top_per_filter:
        # rank over the union of the states that are drawn: the member states in --cell_groups mode, the
        # listed states otherwise. Both give one ranking shared by every panel, as the docstring states.
        picks = top_isoform_per_filter(da, ctx, gene, rank_states, conditions, a.isoform_filter)
        before = {p["condition"]: p["n_molecules"] for p in res.get("panels", [])}
        res = keep_isoforms(res, {k for _, k, _, _ in picks})
        print("[isv_web_remote] top_per_filter: "
              + "; ".join(f"{pat} -> {k} ({v:.0f} reads, top of {n} isoforms carrying it)" for pat, k, v, n in picks)
              + "; panel reads " + ", ".join(f"{p['condition']}:{p['n_molecules']} of {before[p['condition']]}"
                                             for p in res.get("panels", [])), flush=True)
    if a.top_isoforms:
        res, kept, by_iso = keep_top_isoforms(res, a.top_isoforms)
        tot = sum(by_iso.values()) or 1.0
        print(f"[isv_web_remote] top_isoforms={a.top_isoforms}: kept "
              + ", ".join(f"{k} ({by_iso[k]:.0f} reads, {100 * by_iso[k] / tot:.1f}%)" for k in kept)
              + f" of {len(by_iso)} isoforms; panel reads "
              + ", ".join(f"{p['condition']}:{p['n_molecules']}" for p in res.get("panels", [])), flush=True)
    if a.isoforms:
        # explicit isoform ids: keep exactly those final isoforms and colour them in the order given.
        # A name matches either the final isoform id as stored, or that id without its "<gene>:" prefix,
        # because the pseudobulk matrices key the same isoform as ENSG00000159216:ENST00000300305 while
        # the read index stores ENST00000300305. An ambiguous short name is an error, never a guess.
        present, alias = {}, {}
        for p in res.get("panels", []):
            for m in p.get("molecules", []):
                k = m.get("final_isoform_id") or m.get("isoform_id")
                present[k] = present.get(k, 0.0) + float(m.get("count", 0) or 0)
        for k in present:
            alias.setdefault(k.split(":", 1)[-1], set()).add(k)
        want, missing, ambiguous = [], [], []
        for name in dict.fromkeys(a.isoforms):
            if name in present:
                want.append(name)
                continue
            hits = alias.get(name.split(":", 1)[-1], set())
            if len(hits) == 1:
                want.append(next(iter(hits)))
            elif hits:
                ambiguous.append((name, sorted(hits)))
            else:
                missing.append(name)
        if missing or ambiguous:
            sys.exit(f"ERROR: --isoforms against the {len(present)} isoforms with reads in this selection: "
                     f"not found {missing}; ambiguous {ambiguous}")
        res = keep_isoforms(res, set(want))
        rank = {k: i for i, k in enumerate(want)}
        for p in res.get("panels", []):
            for m in p["molecules"]:
                m["color_index"] = rank[m.get("final_isoform_id") or m.get("isoform_id")]
        print(f"[isv_web_remote] isoforms={want}: kept "
              + ", ".join(f"{k} ({present[k]:.0f} reads)" for k in want)
              + f" of {len(present)} isoforms with reads; panel reads "
              + ", ".join(f"{p['condition']}:{p['n_molecules']}" for p in res.get("panels", [])), flush=True)
    if a.max_reads and not a.hide_bleeding_reads:
        res, cap = cap_panel_reads(res, a.max_reads, a.seed)
        print(f"[isv_web_remote] max_reads={a.max_reads} -> cap {cap}: panel reads "
              + ", ".join(f"{p['condition']}:{int(round(sum(m['count'] for m in p['molecules'])))}"
                          for p in res.get("panels", [])), flush=True)
    if a.min_isoform_pct and a.min_isoform_pct > 0:
        res = filter_low_fraction_isoforms(res, a.min_isoform_pct)
    if not a.no_protein_length:
        # the read-panel group labels carry the length too, so this runs with or without the track
        for p in res.get("panels", []):
            for m in p.get("molecules", []):
                if m.get("protein_length") in (None, ""):
                    _pm = ctx.protein_meta(m.get("final_isoform_id") or m.get("isoform_id")) or {}
                    m["protein_length"] = _pm.get("protein_length")
                    m["nmd_status"] = m.get("nmd_status") or _pm.get("nmd_status")
    if a.isoform_track or a.hide_bleeding_reads:
        # attach CDS genomic span to each plotted molecule for the bottom structure track (UCSC CDS/UTR)
        for p in res.get("panels", []):
            for m in p.get("molecules", []):
                _fin = m.get("final_isoform_id") or m.get("isoform_id")
                cds = ctx.coding_region(_fin)
                m["cds_min"], m["cds_max"] = (cds[0], cds[1]) if cds else (None, None)
                m["cds_chrom"] = ctx.coding_chromosome(_fin)
        # attach each final isoform's REGISTERED exon structure from gff-output/combined.gff.gz, so the
        # bottom structure track draws the registered isoform rather than one representative read.
        # Absent index -> registered_exons stays None and render() falls back, with a printed count.
        if res.get("has_multiple_loci"):
            final_tokens = da._load_structures(ctx, gene)
            for panel in res.get("panels", []):
                for m in panel.get("molecules", []):
                    fin = m.get("final_isoform_id") or m.get("isoform_id")
                    tokens = final_tokens.get(fin)
                    if not tokens:
                        raise ValueError("Registered transcript annotation required for multi-locus ISV: " + str(fin))
                    m["registered_exon_segments"] = [da._seg_json(s) for s in
                        da.isv.build_isoform_segments(tokens, ctx.exon_lookup, gene)]
        from .precompute_gff_structures import StructureIndex
        sidx = StructureIndex.open(a.run_dir)
        if sidx is not None:
            gm = res.get("gene_model") or {}
            locus = ((gm.get("chrom"), gm.get("strand"), gm["extent"][0], gm["extent"][1])
                     if gm.get("chrom") and gm.get("extent") else None)
            resolved = set(); unresolved = set()
            for p in res.get("panels", []):
                for m in p.get("molecules", []):
                    fin = m.get("final_isoform_id") or m.get("isoform_id")
                    ex = sidx.exons(gene, fin, locus=locus)
                    m["registered_exons"] = ex
                    if ex and res.get("has_multiple_loci"):
                        segs = m["registered_exon_segments"]
                        if len({s.get("gene") for s in segs}) == 1:
                            model_gene = segs[0]["gene"]
                            m["registered_exon_segments"] = [dict(gene=model_gene, start=a, end=b, type="E")
                                                             for a, b in ex]
                    (resolved if ex else unresolved).add(fin)
            print("[isv_web_remote] registered structures: {} of {} final isoform(s) resolved from "
                  "{}".format(len(resolved), len(resolved) + len(unresolved), sidx.path), flush=True)
            print("[isv_web_remote]   lookup outcomes (unique isoforms) locus={}: {}".format(
                locus, dict(sidx.stats)), flush=True)
            for fin, why in sidx.rejected:
                print(f"[isv_web_remote]   REJECTED registered structure for {fin}: {why}", flush=True)
            if unresolved:
                print("[isv_web_remote]   NOT registered (representative read will be drawn): {}".format(
                    ", ".join(sorted(unresolved)[:10]) + (" ..." if len(unresolved) > 10 else "")), flush=True)
    if a.hide_bleeding_reads:
        from .read_filters import cap_reads_per_isoform
        res = da.prepare_bleeding_read_filter(ctx, res, a.bleed_tolerance)
        res = cap_reads_per_isoform(res, a.max_reads_per_isoform or 200, a.seed)
        if a.max_reads:
            res, _ = cap_panel_reads(res, a.max_reads, a.seed)
    elif a.max_reads_per_isoform:
        from .read_filters import cap_reads_per_isoform
        res = cap_reads_per_isoform(res, a.max_reads_per_isoform, a.seed)
    print("[isv_web_remote] source={} min_isoform_pct={} panels={}".format(
        res.get("source", "index"), a.min_isoform_pct,
        ", ".join(f"{p['condition']}:{p['n_molecules']}" for p in res.get("panels", []))), flush=True)
    fs = a.font_scale if a.font_scale is not None else (2.6 if a.figure_mode else 1.0)
    # title: italic gene symbol, then the event written as <gene>:<junction A>|<junction B>. The filter
    # patterns use '|' between the two exons of ONE junction ("E3.1|E17.1"); the event string uses '-'.
    t_gene = a.title if a.title is not None else ((res.get("symbol") or a.gene) if a.figure_mode else None)
    t_event = a.title_event
    if t_event is None and a.figure_mode and a.isoform_filter:
        t_event = f"{gene}:" + "|".join(str(x).replace("|", "-", 1) for x in a.isoform_filter)
    t_gene = t_gene or None
    t_event = t_event or None
    if a.figure_mode or a.font_scale is not None:
        print(f"[isv_web_remote] figure_mode={a.figure_mode} font_scale={fs} ruler={not a.figure_mode} "
              f"read_labels={not a.figure_mode} trim_to_reads={a.figure_mode} region={_region}", flush=True)
    sys.exit(0 if _render_selected(a, res, out, _region) else 2)


if __name__ == "__main__":
    main()
