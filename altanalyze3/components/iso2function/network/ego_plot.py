"""Side-by-side cell-state ego-networks for a single regulator (e.g. MAX), with identical layout and
node size/color encoding the per-state abundance (isoform abundance for the regulator, gene abundance for
targets/partners) on a shared scale, so the two panels are directly comparable.

PDI edges (regulator -> DNA target) are drawn as directed arrows; PPI edges (physical partners) as lines.
"""

import os
import argparse
import logging

import numpy as np

from . import coexpr

logger = logging.getLogger(__name__)


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    return plt


def _mean(tpm, donors, state, col):
    rows = [(state, d) for d in donors.get(state, ()) if (state, d) in tpm.index]
    if not rows or col not in tpm.columns:
        return 0.0
    return float(tpm.loc[rows, col].mean())


def _draw_edge(ax, p0, p1, color, lw, directed=True, alpha=0.6, shrinkA=0.14, shrinkB=0.14):
    """Draw an edge as a SINGLE 2-point line (a 2-anchor Line2D, NEVER a multi-segment polyline) plus,
    for directed edges, one triangular arrowhead MARKER just outside the target node. shrinkA/shrinkB are
    the source/target node radii (data units) so endpoints sit on the node rim, not buried inside it."""
    x0, y0 = p0
    x1, y1 = p1
    dx, dy = x1 - x0, y1 - y0
    L = (dx * dx + dy * dy) ** 0.5
    if L < 1e-6:
        return
    ux, uy = dx / L, dy / L
    a = (x0 + ux * shrinkA, y0 + uy * shrinkA)
    tip = (x1 - ux * shrinkB, y1 - uy * shrinkB)        # just outside the target node rim
    if directed:
        from matplotlib.patches import Polygon
        hl, hw = 0.14, 0.085                           # arrowhead length / half-width in DATA units (scales)
        base = (tip[0] - ux * hl, tip[1] - uy * hl)
        ax.plot([a[0], base[0]], [a[1], base[1]], color=color, lw=lw, alpha=alpha,
                solid_capstyle="round", zorder=1)      # single 2-point shaft, ends at arrowhead base
        px, py = -uy, ux
        ax.add_patch(Polygon([tip, (base[0] + px * hw, base[1] + py * hw),
                              (base[0] - px * hw, base[1] - py * hw)], closed=True,
                             facecolor=color, edgecolor="none", alpha=alpha, zorder=2))
    else:
        ax.plot([a[0], tip[0]], [a[1], tip[1]], color=color, lw=lw, alpha=alpha,
                solid_capstyle="round", zorder=1)      # undirected (PPI): single 2-point line, no head


def plot_regulator(dataset_dir, state_a, state_b, regulator, out_path, regulator_enst=None,
                   max_neighbors=28):
    plt = _mpl()
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.patches import FancyArrowPatch
    cmap = LinearSegmentedColormap.from_list("white_red", ["#FFFFFF", "#FF0000"])

    edges = coexpr.load_edges()
    sym2ensg = coexpr.load_symbol_to_ensg(edges)
    iso_tpm, gene_tpm, donors = coexpr.load_groups(dataset_dir, [state_a, state_b])
    gene_tpm = coexpr.reference_isoform_expr(iso_tpm)   # target scored at reference-isoform resolution

    me = edges[edges["Symbol"] == regulator].copy()
    me["tgt_ensg"] = me["target"].map(sym2ensg).fillna("")
    reg_ensts = [e for e in me["best_ENST"].unique() if e in iso_tpm.columns]
    if regulator_enst:
        reg_ensts = [regulator_enst]
    if not reg_ensts:
        raise SystemExit("regulator %s has no expressed isoform in this dataset" % regulator)
    # use the most-abundant regulator isoform as THE regulator node
    reg_enst = max(reg_ensts, key=lambda e: _mean(iso_tpm, donors, state_a, e) + _mean(iso_tpm, donors, state_b, e))
    me = me[(me["best_ENST"] == reg_enst) & (me["tgt_ensg"] != "") & me["tgt_ensg"].isin(gene_tpm.columns)]

    # one node per neighbor; PDI takes precedence over PPI for typing
    nb, aff = {}, {}
    for r in me.itertuples(index=False):
        prev = nb.get(r.target)
        if prev is None or (r.interaction_type == "PDI" and prev[0] != "PDI"):
            nb[r.target] = (r.interaction_type, r.tgt_ensg)
        aff[r.target] = max(aff.get(r.target, 0.0), float(r.activity_score))
    # rank neighbors by abundance, keep the top N
    def nb_ab(ensg):
        return _mean(gene_tpm, donors, state_a, ensg) + _mean(gene_tpm, donors, state_b, ensg)
    # rank by abundance and keep only neighbours expressed in >=1 of the two states (drop all-zero genes)
    neigh = [kv for kv in sorted(nb.items(), key=lambda kv: nb_ab(kv[1][1]), reverse=True)
             if nb_ab(kv[1][1]) > 0][:max_neighbors]
    pdi = [(g, e) for g, (t, e) in neigh if t == "PDI"]
    ppi = [(g, e) for g, (t, e) in neigh if t == "PPI"]
    logger.info("[ego] %s isoform %s: %d PDI targets + %d PPI partners (top %d shown)",
                regulator, reg_enst, len(pdi), len(ppi), len(neigh))

    # ---- shared radial layout (PDI on the upper arc, PPI on the lower arc), identical in both panels
    pos = {"__REG__": (0.0, 0.0)}
    R = 1.0
    def place(items, a0, a1):
        if not items:
            return
        angs = np.linspace(a0, a1, len(items), endpoint=True) if len(items) > 1 else [ (a0 + a1) / 2 ]
        for (g, _), ang in zip(items, np.radians(angs)):
            pos[g] = (R * np.cos(ang), R * np.sin(ang))
    place(pdi, 165, 15)        # upper arc
    place(ppi, 195, 345)       # lower arc

    # ---- abundances per state (regulator isoform TPM; neighbor gene TPM), shared scale
    ab = {st: {} for st in (state_a, state_b)}
    for st in (state_a, state_b):
        ab[st]["__REG__"] = _mean(iso_tpm, donors, st, reg_enst)
        for g, (t, e) in neigh:
            ab[st][g] = _mean(gene_tpm, donors, st, e)
    allv = [v for st in ab for v in ab[st].values()]
    vmax = max(allv) if allv else 1.0
    cnorm = Normalize(0, np.log1p(vmax))

    def size(v):
        return 120 + 1400 * np.sqrt(max(v, 0) / (vmax + 1e-9))

    # per-state edge activity EA = w * sqrt(P_regulator * P_target) -> edge width (shared scale)
    def edge_EA(st, g):
        return aff.get(g, 0.0) * np.sqrt(max(ab[st]["__REG__"], 0) * max(ab[st].get(g, 0.0), 0))
    ea_max = max([edge_EA(st, g) for st in (state_a, state_b) for g, _ in neigh] + [1e-9])

    def ew(st, g):
        return 0.4 + 5.5 * (edge_EA(st, g) / ea_max)

    fig, axes = plt.subplots(1, 2, figsize=(13, 6.5))
    for ax, st in zip(axes, (state_a, state_b)):
        # edges first (single 2-point paths; arrows for PDI, lines for PPI); width = interaction activity
        for g, _ in pdi:
            _draw_edge(ax, pos["__REG__"], pos[g], "#7A4FA3", ew(st, g), directed=True)
        for g, _ in ppi:
            _draw_edge(ax, pos["__REG__"], pos[g], "#888888", ew(st, g), directed=False)
        # nodes
        for name, (x, y) in pos.items():
            v = ab[st].get(name, 0.0)
            is_reg = name == "__REG__"
            ax.scatter([x], [y], s=size(v) * (1.6 if is_reg else 1.0),
                       c=[cmap(cnorm(np.log1p(v)))], marker="D" if is_reg else "o",
                       edgecolors="#222222", linewidths=1.3 if is_reg else 0.6, zorder=3)
            label = "%s\n(%s)" % (regulator, reg_enst.split(".")[0]) if is_reg else name
            ax.annotate(label, (x, y), fontsize=7 if is_reg else 6, ha="center",
                        va="center" if is_reg else "bottom",
                        xytext=(0, 0 if is_reg else 11), textcoords="offset points", zorder=4)
        ax.set_title("%s\nnode = abundance (CPM); edge width = interaction activity" % st, fontsize=10)
        ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.55); ax.set_aspect("equal"); ax.axis("off")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=cnorm); sm.set_array([])
    cb = fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.02)
    cb.set_label("log1p CPM (shared scale)")
    fig.suptitle("%s regulator network — %s vs %s  (purple=PDI target, grey=PPI partner)"
                 % (regulator, state_a, state_b), fontsize=11)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    logger.info("[ego] wrote %s", out_path)


def plot_regulator_isoforms(dataset_dir, state_a, state_b, regulator, out_path, group="young",
                            tgt_min=1.0, iso_min=1.0):
    """Side-by-side networks of ALL catalogued isoforms of `regulator` with functional interactions.
    Each regulator isoform is a node whose abundance is STRUCTURE-RESOLVED (sum of the expressed
    transcripts the file catalogues for that isoform's structure); its targets are drawn around it with
    abundance = the MAX isoform CPM of the target gene. Nodes/edges below `iso_min`/`tgt_min` CPM are
    hidden; donors restricted to `group`. Node size/colour (grey->red) = abundance (CPM)."""
    plt = _mpl()
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.patches import FancyArrowPatch
    from matplotlib.lines import Line2D
    cmap = LinearSegmentedColormap.from_list("white_red", ["#FFFFFF", "#FF0000"])
    ISO_COLORS = ["#D1495B", "#2E86AB", "#3C6E47", "#B5651D", "#7A4FA3", "#8A8D00"]

    edges = coexpr.load_edges()
    sym2ensg = coexpr.load_symbol_to_ensg(edges)
    iso_tpm, _, donors = coexpr.load_groups(dataset_dir, [state_a, state_b])
    gene_tpm = coexpr.load_max_isoform_expr(dataset_dir, [state_a, state_b])     # target = max-isoform CPM
    if group:
        from .cross_state import young_libraries
        keep = young_libraries(dataset_dir, group)
        donors = {s: [d for d in ds if d in keep] for s, ds in donors.items()}
    me = edges[edges["Symbol"] == regulator].copy()
    me["tgt_ensg"] = me["target"].map(sym2ensg).fillna("")

    def iso_ab(st, ensts):                # regulator isoform abundance = sum of its catalogued ENSTs (CPM)
        rows = [(st, d) for d in donors[st] if (st, d) in iso_tpm.index]
        return float(iso_tpm.loc[rows, ensts].sum(axis=1).mean()) if rows and ensts else 0.0

    # catalogued isoforms of the regulator -> expressed ENSTs (structure-resolved); apply CPM thresholds
    isos = []
    for sid, grp in me.groupby("source_isoform_id"):
        ensts = sorted({e.split(".")[0] for e in grp["best_ENST"] if e and e.split(".")[0] in iso_tpm.columns})
        if not ensts or max(iso_ab(state_a, ensts), iso_ab(state_b, ensts)) < iso_min:
            continue                                       # hide regulator isoforms below threshold
        tgts = {}
        for r in grp.itertuples(index=False):
            if r.tgt_ensg and r.tgt_ensg in gene_tpm.columns:
                tgts.setdefault(r.target, (r.interaction_type, r.tgt_ensg))
        # target abundance = MAX isoform CPM of the gene; keep targets >= tgt_min in >=1 state
        tgts = {g: v for g, v in tgts.items()
                if max(_mean(gene_tpm, donors, state_a, v[1]), _mean(gene_tpm, donors, state_b, v[1])) >= tgt_min}
        if tgts:
            isos.append((sid, ensts, tgts))
    isos.sort(key=lambda x: -len(x[2]))
    logger.info("[ego] %s: %d isoforms with functional interactions above threshold (iso>=%g,tgt>=%g CPM)",
                regulator, len(isos), iso_min, tgt_min)

    all_targets = sorted({g for _, _, t in isos for g in t})
    # shared layout: regulator isoforms across the bottom, target genes on the top semicircle
    pos, R = {}, 1.0
    nri = len(isos)
    for i, (sid, _, _) in enumerate(isos):
        x = 0.0 if nri == 1 else -0.75 + 1.5 * i / (nri - 1)
        pos[("ISO", sid)] = (x, -1.05)
    for j, g in enumerate(all_targets):
        ang = np.radians(np.linspace(160, 20, len(all_targets))[j] if len(all_targets) > 1 else 90)
        pos[("TGT", g)] = (1.35 * np.cos(ang), 0.15 + 1.05 * np.sin(ang))

    ab = {st: {} for st in (state_a, state_b)}
    for st in (state_a, state_b):
        for sid, ensts, _ in isos:
            ab[st][("ISO", sid)] = iso_ab(st, ensts)
        for g in all_targets:
            e = next(v[1] for _, _, t in isos for k, v in t.items() if k == g)
            ab[st][("TGT", g)] = _mean(gene_tpm, donors, st, e)
    vmax = max([v for st in ab for v in ab[st].values()] + [1.0])
    cnorm = Normalize(0, np.log1p(vmax))
    sz = lambda v: 120 + 1300 * np.sqrt(max(v, 0) / (vmax + 1e-9))

    def expressed(key):                  # is this node above threshold IN this state's panel?
        return ab_st.get(key, 0.0) >= (iso_min if key[0] == "ISO" else tgt_min)

    fig, axes = plt.subplots(1, 2, figsize=(15, 7.5))
    for ax, st in zip(axes, (state_a, state_b)):
        ab_st = ab[st]
        # edges only when BOTH endpoints are expressed in this state (no edges where not expressed)
        for i, (sid, ensts, tgts) in enumerate(isos):
            col = ISO_COLORS[i % len(ISO_COLORS)]
            for g, (typ, _) in tgts.items():
                if not (expressed(("ISO", sid)) and expressed(("TGT", g))):
                    continue
                _draw_edge(ax, pos[("ISO", sid)], pos[("TGT", g)], col, 1.2, directed=(typ == "PDI"))
        for key, (x, y) in pos.items():
            v = ab_st.get(key, 0.0)
            is_iso = key[0] == "ISO"
            alpha = 1.0 if expressed(key) else 0.2          # below-threshold nodes are near-transparent
            ax.scatter([x], [y], s=sz(v) * (1.5 if is_iso else 1.0), c=[cmap(cnorm(np.log1p(v)))],
                       marker="D" if is_iso else "o", edgecolors="#222222",
                       linewidths=1.3 if is_iso else 0.6, alpha=alpha, zorder=3)
            ax.annotate(key[1], (x, y), fontsize=6.5 if is_iso else 6, ha="center",
                        va="top" if is_iso else "bottom", alpha=alpha,
                        xytext=(0, -12 if is_iso else 10), textcoords="offset points", zorder=4)
        ax.set_title("%s\nnode size/colour = abundance (CPM)" % st, fontsize=10)
        ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.55, 1.45); ax.set_aspect("equal"); ax.axis("off")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=cnorm); sm.set_array([])
    cb = fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.02); cb.set_label("log1p CPM (shared scale)")
    handles = [Line2D([0], [0], color=ISO_COLORS[i % len(ISO_COLORS)], lw=2, label=sid)
               for i, (sid, _, _) in enumerate(isos)]
    axes[0].legend(handles=handles, loc="upper left", fontsize=6, title="%s isoform" % regulator, frameon=False)
    fig.suptitle("%s — all catalogued isoforms with functional interactions — %s vs %s"
                 % (regulator, state_a, state_b), fontsize=11)
    fig.savefig(out_path, bbox_inches="tight"); plt.close(fig)
    logger.info("[ego] wrote %s", out_path)


def plot_regulator_states(dataset_dir, regulator, states, out_path, group="young",
                          tgt_min=1.0, iso_min=1.0, iso_tpm=None, gene_tpm=None, donors=None,
                          only_connected=False):
    """Single figure, one panel per cell state (same node layout), showing `regulator`'s isoform network
    across the given states — to read network REWIRING/divergence at a glance. Regulator isoform
    abundance is structure-resolved; target abundance = max-isoform CPM; below-threshold nodes are 20%
    opacity and carry no edges (white->red, shared CPM scale)."""
    import math
    plt = _mpl()
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.patches import FancyArrowPatch
    from matplotlib.lines import Line2D
    cmap = LinearSegmentedColormap.from_list("white_red", ["#FFFFFF", "#FF0000"])
    ISO_COLORS = ["#D1495B", "#2E86AB", "#3C6E47", "#B5651D", "#7A4FA3", "#8A8D00"]

    edges = coexpr.load_edges()
    sym2ensg = coexpr.load_symbol_to_ensg(edges)
    if iso_tpm is None:                       # load if not supplied by the caller (driver pre-loads once)
        iso_tpm, _, donors = coexpr.load_groups(dataset_dir, states)
        gene_tpm = coexpr.load_max_isoform_expr(dataset_dir, states)
        if group:
            from .cross_state import young_libraries
            keep = young_libraries(dataset_dir, group)
            donors = {s: [d for d in donors.get(s, []) if d in keep] for s in states}
    me = edges[edges["Symbol"] == regulator].copy()
    me["tgt_ensg"] = me["target"].map(sym2ensg).fillna("")

    def iso_ab(st, ensts):
        rows = [(st, d) for d in donors.get(st, []) if (st, d) in iso_tpm.index]
        return float(iso_tpm.loc[rows, ensts].sum(axis=1).mean()) if rows and ensts else 0.0

    isos = []
    for sid, grp in me.groupby("source_isoform_id"):
        ensts = sorted({e.split(".")[0] for e in grp["best_ENST"] if e and e.split(".")[0] in iso_tpm.columns})
        if not ensts or max((iso_ab(s, ensts) for s in states), default=0) < iso_min:
            continue
        tgts = {}
        for r in grp.itertuples(index=False):
            if r.tgt_ensg and r.tgt_ensg in gene_tpm.columns:
                tgts.setdefault(r.target, (r.interaction_type, r.tgt_ensg))
        tgts = {g: v for g, v in tgts.items()
                if max((_mean(gene_tpm, donors, s, v[1]) for s in states), default=0) >= tgt_min}
        if tgts:
            isos.append((sid, ensts, tgts))
    isos.sort(key=lambda x: -len(x[2]))
    # order targets by isoform-connectivity SIGNATURE (which isoforms regulate each target):
    # iso1-unique, {1,2}, iso2-unique, {1,3}, {2,3}, iso3-unique, ... ; targets shared by >2 isoforms
    # go to the far right.
    iso_of = {}
    for i, (sid, ensts, tgts) in enumerate(isos):
        for g in tgts:
            iso_of.setdefault(g, set()).add(i + 1)

    def _sig_rank(g):
        s = tuple(sorted(iso_of.get(g, ())))
        if len(s) >= 3:
            return (10 ** 9, 0, g)                 # shared by >2 isoforms -> far right
        k = s[-1] if s else 10 ** 9
        return (k, (s[0] if len(s) == 2 else k), g)
    all_targets = sorted({g for _, _, t in isos for g in t}, key=_sig_rank)
    logger.info("[ego] %s across %d states: %d isoforms, %d targets", regulator, len(states), len(isos), len(all_targets))

    # half-dome ELLIPSE of targets just over the TFs: horizontal radius grows with target count (space),
    # vertical radius stays small so the dome doesn't float high above the TFs; base ~ TF level.
    pos = {}
    nri = len(isos)
    N = len(all_targets)
    Rx = max(1.6, 0.085 * N)
    Ry = 1.5
    y0 = -1.0
    for i, (sid, _, _) in enumerate(isos):
        pos[("ISO", sid)] = (0.0 if nri == 1 else -0.9 + 1.8 * i / (nri - 1), y0)
    for j, g in enumerate(all_targets):
        ang = np.radians(np.linspace(178, 2, N)[j] if N > 1 else 90)
        pos[("TGT", g)] = (Rx * np.cos(ang), y0 + Ry * np.sin(ang))

    ab = {st: {} for st in states}
    for st in states:
        for sid, ensts, _ in isos:
            ab[st][("ISO", sid)] = iso_ab(st, ensts)
        for g in all_targets:
            e = next(v[1] for _, _, t in isos for k, v in t.items() if k == g)
            ab[st][("TGT", g)] = _mean(gene_tpm, donors, st, e)
    from matplotlib.patches import Circle, RegularPolygon
    vmax = max([v for st in ab for v in ab[st].values()] + [1.0])
    cnorm = Normalize(0, np.log1p(vmax))
    # node RADIUS in data units (so edges can stop on the rim); capped to avoid overlap in the dome
    def radius(v, is_iso):
        r = 0.05 + 0.10 * np.sqrt(max(v, 0) / (vmax + 1e-9))
        return r * 1.5 if is_iso else r

    ncol = min(3, len(states))
    nrow = math.ceil(len(states) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.1 * (Rx + 0.5) * ncol, (Ry + 1.9) * nrow), squeeze=False)
    axflat = axes.flat
    for ax, st in zip(axflat, states):
        ab_st = ab[st]
        exp = lambda key: ab_st.get(key, 0.0) >= (iso_min if key[0] == "ISO" else tgt_min)
        rad = {key: radius(ab_st.get(key, 0.0), key[0] == "ISO") for key in pos}
        connected = set()                      # nodes participating in >=1 drawn edge this panel
        for i, (sid, ensts, tgts) in enumerate(isos):
            col = ISO_COLORS[i % len(ISO_COLORS)]
            for g, (typ, _) in tgts.items():
                if not (exp(("ISO", sid)) and exp(("TGT", g))):
                    continue
                _draw_edge(ax, pos[("ISO", sid)], pos[("TGT", g)], col, 1.1, directed=(typ == "PDI"),
                           shrinkA=rad[("ISO", sid)] + 0.01, shrinkB=rad[("TGT", g)] + 0.025)
                connected.add(("ISO", sid)); connected.add(("TGT", g))
        for key, (x, y) in pos.items():
            if only_connected and key not in connected:
                continue                       # hide edgeless/hanging nodes for this cell-state/condition
            v = ab_st.get(key, 0.0)
            is_iso = key[0] == "ISO"
            alpha = 1.0 if exp(key) else 0.2
            fc = cmap(cnorm(np.log1p(v)))
            r = rad[key]
            if is_iso:
                ax.add_patch(RegularPolygon((x, y), 4, radius=r, facecolor=fc, edgecolor="#222222",
                                            lw=1.1, alpha=alpha, zorder=3))
                ax.annotate(key[1].replace("|", "|\n"), (x, y - r - 0.04), fontsize=4.5, ha="center",
                            va="top", alpha=alpha, zorder=4)
            else:
                ax.add_patch(Circle((x, y), radius=r, facecolor=fc, edgecolor="#222222",
                                    lw=0.5, alpha=alpha, zorder=3))
                ax.annotate(key[1], (x, y), fontsize=4, ha="center", va="center", alpha=alpha,
                            color="#000000", zorder=4)         # gene label centred inside the circle
        ax.set_title(st, fontsize=9)
        ax.set_xlim(-Rx - 0.45, Rx + 0.45); ax.set_ylim(y0 - 0.35, y0 + Ry + 0.35)
        ax.set_aspect("equal"); ax.axis("off")
    for ax in list(axflat):
        ax.axis("off")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=cnorm); sm.set_array([])
    cb = fig.colorbar(sm, ax=axes, fraction=0.018, pad=0.01); cb.set_label("log1p CPM (max-isoform)")
    handles = [Line2D([0], [0], color=ISO_COLORS[i % len(ISO_COLORS)], lw=2, label=sid)
               for i, (sid, _, _) in enumerate(isos)]
    fig.legend(handles=handles, loc="lower center", ncol=min(len(isos), 6), fontsize=6,
               title="%s isoform (edge colour)" % regulator, frameon=False)
    sub = ("%s vs %s" % (states[0], states[1])) if len(states) == 2 \
        else "the %d most-divergent cell states" % len(states)
    fig.suptitle("%s network — %s  (edge = supported interaction; node = max-isoform CPM)"
                 % (regulator, sub), fontsize=11)
    fig.savefig(out_path, bbox_inches="tight"); plt.close(fig)
    logger.info("[ego] wrote %s", out_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--state-a", required=True)
    ap.add_argument("--state-b", required=True)
    ap.add_argument("--regulator", required=True)
    ap.add_argument("--enst", default=None)
    ap.add_argument("--all-isoforms", action="store_true",
                    help="show ALL catalogued isoforms of the regulator with functional interactions")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    if a.all_isoforms:
        plot_regulator_isoforms(a.dataset, a.state_a, a.state_b, a.regulator, a.out)
    else:
        plot_regulator(a.dataset, a.state_a, a.state_b, a.regulator, a.out, a.enst)


if __name__ == "__main__":
    main()
