"""Switch-pair visualizations: for a pair of isoforms of the same TF that swap dominance across cell
states, render (a) TF->DNA-bait interaction networks in the states with the largest isoform-ratio shift,
(b) the per-state isoform-usage line plot, and (c) a paired UMAP.

These are library functions: they take ALREADY-COMPUTED tidy inputs (a `pair` descriptor + a per-state
usage table) so the same drawing works for any cohort / dataset. The PDI baits and M1H activation values
are pulled from the iso2function data tables via ``config._load`` (``pdi_ey1h`` and ``activation_m1h``),
so this stays anchored to the atlas the rest of the component uses.

Drawing conventions match the rest of network/: edges are single 2-point paths + a triangular arrowhead
(``ego_plot._draw_edge``), isoform nodes are diamonds, targets are circles, colours are explicit RGB hex.

Input contract
--------------
pair : Mapping with keys
    gene, ensg, isoA, isoB      -- gene symbol, ENSG, and the two isoform ids (ENST or clone-mapped)
    A_clone, B_clone            -- TFIsoDB clone ids (used to look up baits + M1H)
    A_m1h_call, B_m1h_call      -- optional activator/repressor/neither labels for titles
state_usage : pandas.DataFrame indexed by cell-state name, columns
    usage_A, usage_B            -- fraction of the gene's expression carried by isoform A / B in that state
    ratio_A                     -- A / (A + B) in that state
    n_donors                    -- donors contributing (annotated; states with <4 can never reach p<0.05)
    lineage                     -- optional lineage label, for the lineage-level line plot
"""
import os
import logging

import numpy as np
import pandas as pd

from .. import config
from . import ego_plot

logger = logging.getLogger("altanalyze3")

CA, CB, CBOTH = "#F2B705", "#2E86C1", "#BDBDBD"   # isoform A gold / B blue / shared grey


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    matplotlib.rcParams["font.family"] = "Arial"
    import matplotlib.pyplot as plt
    return plt


def load_pdi_m1h(data_dir=None):
    """Return (detected_baits_per_clone, tested_baits_per_clone, m1h_frame) from the iso2function tables.

    detected/tested are {clone_id -> set(bait_id)}; m1h_frame has clone_id, M1H_rep1..3, M1H_mean, m1h_call.
    Reuses the same ``_load`` + ``pdi_ey1h`` / ``activation_m1h`` tables as ``associate.py``.
    """
    from ..associate import _load, _detected_mask
    data_dir = data_dir or config.DATA_DIR
    ey1h = _load(data_dir, "pdi_ey1h")
    det, tes = {}, {}
    if ey1h is not None and len(ey1h):
        d = ey1h[_detected_mask(ey1h)]
        det = d.groupby("clone_id")["bait_id"].apply(lambda s: set(x for x in s if x)).to_dict()
        tes = ey1h.groupby("clone_id")["bait_id"].apply(lambda s: set(x for x in s if x)).to_dict()
    m1h = _load(data_dir, "activation_m1h")
    if m1h is not None:
        for c in ("M1H_rep1", "M1H_rep2", "M1H_rep3", "M1H_mean"):
            if c in m1h.columns:
                m1h[c] = pd.to_numeric(m1h[c], errors="coerce")
    return det, tes, m1h


def _clone_set(clone_field):
    """A clone field may hold >1 collapsed clone id (pipe-joined); return the list."""
    return [c.strip() for c in str(clone_field).split("|") if c.strip()]


def _shared_baits(det, tes, cloneA, cloneB):
    """Baits BOTH clones were tested on (so 'A binds / B does not' is a real difference, not a panel gap).
    Returns (shared_detected_by_A, shared_detected_by_B, all_shared_or_detected)."""
    tA = set().union(*[tes.get(c, set()) for c in _clone_set(cloneA)]) if cloneA else set()
    tB = set().union(*[tes.get(c, set()) for c in _clone_set(cloneB)]) if cloneB else set()
    dA = set().union(*[det.get(c, set()) for c in _clone_set(cloneA)]) if cloneA else set()
    dB = set().union(*[det.get(c, set()) for c in _clone_set(cloneB)]) if cloneB else set()
    shared = tA & tB
    if not shared:
        return dA, dB, sorted(dA | dB)
    return dA & shared, dB & shared, sorted((dA | dB) & shared)


ISO_PALETTE = ["#D1495B", "#2E86AB", "#3C6E47", "#B5651D", "#7A4FA3", "#8A8D00"]


def isoform_set_pdi_m1h(gene, clones, out_path, data_dir=None, iso_colors=None):
    """M1H activation + PDI (Y1H DNA-bait) display for a SET of TFIsoDB clones (not just a pair).

    Top panel: a bait network — each clone is a diamond (left) linked by a directed edge to every DNA bait
    it binds (Y1H-detected; circles, right); a clone with no Y1H assay is marked '(not Y1H-tested)'.
    Bottom panel: the M1H activation bar (M1H_mean bar + 3 replicate points, with activator/repressor
    thresholds) for the same clones. Data via load_pdi_m1h(); thresholds via config. Writes ``out_path``
    (.pdf) and returns (path, {clone: (m1h_mean, n_tested, n_detected)})."""
    plt = _mpl()
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import RegularPolygon, Circle
    det, tes, m1h = load_pdi_m1h(data_dir)
    cols = iso_colors or {c: ISO_PALETTE[i % len(ISO_PALETTE)] for i, c in enumerate(clones)}
    g = m1h[m1h["gene_symbol"] == gene].copy() if m1h is not None else pd.DataFrame()
    all_baits = sorted(set().union(*[det.get(c, set()) for c in clones])) if clones else []
    summary = {}

    fig = plt.figure(figsize=(max(8.0, 0.30 * len(all_baits) + 6.0), 8.6))
    gs = GridSpec(2, 1, height_ratios=[2.4, 1.0], hspace=0.32)

    # ---- top: bait network (clones left -> detected DNA baits right) ----
    ax = fig.add_subplot(gs[0])
    ncl = len(clones)
    cl_y = np.linspace(1.0, -1.0, ncl) if ncl > 1 else np.array([0.0])
    bx = 2.2
    by_ = np.linspace(1.25, -1.25, max(len(all_baits), 1))
    bait_pos = {b: (bx, y) for b, y in zip(all_baits, by_)}
    for b, (x, y) in bait_pos.items():
        ax.add_patch(Circle((x, y), 0.045, facecolor="#f0f0f0", edgecolor="#555555", lw=0.7, zorder=3))
        ax.text(x + 0.08, y, b, fontsize=4.6, ha="left", va="center", zorder=4)
    for cy, c in zip(cl_y, clones):
        col = cols[c]
        dcount, tcount = len(det.get(c, set())), len(tes.get(c, set()))
        ax.add_patch(RegularPolygon((0.0, cy), 4, radius=0.075, facecolor=col, edgecolor="#222222",
                                    lw=1.1, zorder=3))
        tag = "(not Y1H-tested)" if tcount == 0 else "%d baits" % dcount
        ax.text(-0.12, cy, "%s\n%s" % (c, tag), fontsize=6.2, ha="right", va="center", zorder=4)
        for b in det.get(c, set()):
            if b in bait_pos:
                ego_plot._draw_edge(ax, (0.075, cy), bait_pos[b], col, 0.9, directed=True,
                                    shrinkA=0.08, shrinkB=0.06)
    ax.set_xlim(-1.15, bx + 0.9); ax.set_ylim(-1.5, 1.5); ax.axis("off")
    ax.set_title("%s PDI — Y1H DNA-bait interactions per isoform (edge = detected binding)" % gene, fontsize=9)

    # ---- bottom: M1H activation bar ----
    axb = fig.add_subplot(gs[1])
    order = list(clones)[::-1]
    for k, cid in enumerate(order):
        row = g[g["clone_id"] == cid]
        mm = float(row["M1H_mean"].iloc[0]) if len(row) else np.nan
        summary[cid] = (mm, len(tes.get(cid, set())), len(det.get(cid, set())))
        axb.barh(k, mm, color=cols[cid], edgecolor="#222222", lw=0.7, height=0.66, zorder=2)
        if len(row):
            for rep in ("M1H_rep1", "M1H_rep2", "M1H_rep3"):
                if rep in row and pd.notna(row[rep].iloc[0]):
                    axb.scatter([float(row[rep].iloc[0])], [k], facecolors="none", edgecolors="#222222",
                                s=24, lw=0.7, zorder=3)
    axb.set_yticks(range(len(order))); axb.set_yticklabels(order, fontsize=7)
    axb.axvline(0, color="#000000", lw=1.0)
    axb.axvline(config.M1H_ACTIVATOR_MIN, color="#000000", lw=0.9, ls="--")
    axb.axvline(config.M1H_REPRESSOR_MAX, color="#000000", lw=0.9, ls="--")
    axb.set_xlabel("M1H log2(activation fold change)  —  dashed = activator/repressor thresholds", fontsize=8)
    axb.set_title("%s M1H activation (mean bar + replicate points)" % gene, fontsize=8)
    for s_ in ("top", "right"):
        axb.spines[s_].set_visible(False)

    fig.suptitle("%s isoforms %s — M1H activation + PDI DNA baits (iso2function)"
                 % (gene, ", ".join(clones)), fontsize=10)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight"); plt.close(fig)
    logger.info("[switch_pair_plots] wrote %s", out_path)
    return out_path, summary


def _net_panel(ax, state, ratio_A, cloneA, cloneB, baitsA, baitsB, usage_A, usage_B):
    baits = sorted(set(baitsA) | set(baitsB))
    ax.scatter([0.0], [0.6], s=300 * (0.3 + usage_A), marker="D", c=CA, edgecolors="#222222", lw=0.6, zorder=3)
    ax.scatter([0.0], [-0.6], s=300 * (0.3 + usage_B), marker="D", c=CB, edgecolors="#222222", lw=0.6, zorder=3)
    ax.text(-0.05, 0.6, f"{cloneA}\n{usage_A*100:.0f}%", fontsize=6, ha="right", va="center")
    ax.text(-0.05, -0.6, f"{cloneB}\n{usage_B*100:.0f}%", fontsize=6, ha="right", va="center")
    ys = np.linspace(1.15, -1.15, max(len(baits), 1))
    for bait, by in zip(baits, ys):
        inA, inB = bait in baitsA, bait in baitsB
        col = CBOTH if (inA and inB) else (CA if inA else CB)
        ax.scatter([1.5], [by], s=42, c="#f0f0f0", edgecolors=col, lw=1.1, zorder=3)
        ax.text(1.56, by, bait, fontsize=4.3, ha="left", va="center")
        if inA:
            ego_plot._draw_edge(ax, (0.06, 0.6), (1.5, by), CA, 0.8, directed=True, shrinkA=0.06, shrinkB=0.05)
        if inB:
            ego_plot._draw_edge(ax, (0.06, -0.6), (1.5, by), CB, 0.8, directed=True, shrinkA=0.06, shrinkB=0.05)
    ax.set_title(f"{state}\nA/(A+B)={ratio_A*100:.0f}%", fontsize=7)
    ax.set_xlim(-0.9, 2.4)
    ax.set_ylim(-1.35, 1.35)
    ax.axis("off")


def bait_network_with_m1h(pair, state_usage, out_path, data_dir=None, top_states=6):
    """Combined figure: top-N network panels (states with the largest |ratio_A - median| shift) + an M1H
    activation bar chart over all clones of the gene. Writes ``out_path`` (.pdf); returns the path."""
    plt = _mpl()
    from matplotlib.gridspec import GridSpec
    det, tes, m1h = load_pdi_m1h(data_dir)
    baitsA, baitsB, _ = _shared_baits(det, tes, pair.get("A_clone"), pair.get("B_clone"))

    su = state_usage.dropna(subset=["ratio_A"]).copy()
    if len(su) < 2:
        logger.warning("[switch_pair_plots] %s: <2 usable states, skipping network", pair.get("gene"))
        return None
    med = np.median(su["ratio_A"].to_numpy())
    su["_shift"] = (su["ratio_A"] - med).abs()
    top = su.sort_values("_shift", ascending=False).head(top_states).sort_values("ratio_A")

    fig = plt.figure(figsize=(2.5 * len(top), 7.4))
    gs = GridSpec(2, len(top), height_ratios=[2.1, 1.0], hspace=0.45, wspace=0.25)
    for c, (st, r) in enumerate(top.iterrows()):
        ax = fig.add_subplot(gs[0, c])
        _net_panel(ax, st, r["ratio_A"], pair["A_clone"], pair["B_clone"], baitsA, baitsB,
                   r["usage_A"], r["usage_B"])

    axb = fig.add_subplot(gs[1, :])
    g = m1h[m1h["gene_symbol"] == pair["gene"]].copy() if m1h is not None else pd.DataFrame()
    order = g["clone_id"].tolist()[::-1]
    for k, cid in enumerate(order):
        row = g[g["clone_id"] == cid].iloc[0]
        bar_c = CA if cid == pair["A_clone"] else (CB if cid == pair["B_clone"] else "#c9c9c9")
        hl = cid in (pair["A_clone"], pair["B_clone"])
        axb.barh(k, row["M1H_mean"], color=bar_c, edgecolor="#222222", lw=0.7 if hl else 0.4, height=0.7, zorder=2)
        for rep in ("M1H_rep1", "M1H_rep2", "M1H_rep3"):
            if rep in row and pd.notna(row[rep]):
                axb.scatter([row[rep]], [k], facecolors="none", edgecolors="#222222", s=22, lw=0.7, zorder=3)
    axb.set_yticks(range(len(order)))
    axb.set_yticklabels(order, fontsize=7)
    axb.axvline(0, color="#000000", lw=1.0)
    axb.axvline(config.M1H_ACTIVATOR_MIN, color="#000000", lw=0.9, ls="--")
    axb.axvline(config.M1H_REPRESSOR_MAX, color="#000000", lw=0.9, ls="--")
    axb.set_xlabel("log2(activation fold change)", fontsize=8)
    axb.set_title(f"{pair['gene']} M1H activation (all clones; {pair['A_clone']}=gold, {pair['B_clone']}=blue)",
                  fontsize=8)
    for s_ in ("top", "right"):
        axb.spines[s_].set_visible(False)

    fig.suptitle(f"{pair['gene']} switch: {pair['A_clone']} [{pair.get('A_m1h_call','')}] vs "
                 f"{pair['B_clone']} [{pair.get('B_m1h_call','')}] — TF->bait interactions in top-{len(top)} "
                 f"isoform-ratio-shift states (edge gold=A binds / blue=B binds / grey=both)", fontsize=9)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    logger.info("[switch_pair_plots] wrote %s (baits A=%d B=%d)", out_path, len(baitsA), len(baitsB))
    return out_path


def usage_line(pair, state_usage, out_path, level="lineage", lineage_order=None):
    """Per-state isoform-usage line plot (isoform A gold, B blue). ``level`` = 'lineage' groups by the
    ``lineage`` column; 'cellstate' uses the cell-state index. Empty groups are dropped so the line is
    continuous. Writes ``out_path``; returns the path."""
    plt = _mpl()
    su = state_usage.copy()
    if level == "lineage":
        if "lineage" not in su.columns:
            raise ValueError("usage_line(level='lineage') needs a 'lineage' column in state_usage")
        agg = su.groupby("lineage").agg(usage_A=("usage_A", "mean"), usage_B=("usage_B", "mean"),
                                        n_donors=("n_donors", "max"))
        if lineage_order:
            agg = agg.reindex([L for L in lineage_order if L in agg.index])
    else:
        agg = su[["usage_A", "usage_B", "n_donors"]]
    agg = agg.dropna(subset=["usage_A", "usage_B"])
    x = np.arange(len(agg))
    fig, ax = plt.subplots(figsize=(max(4.0, 0.5 * len(agg) + 3), 3.6))
    ax.plot(x, agg["usage_A"] * 100, color=CA, marker="o", ms=4, lw=1.6, label=f"{pair['A_clone']} ({pair['isoA']})")
    ax.plot(x, agg["usage_B"] * 100, color=CB, marker="o", ms=4, lw=1.6, label=f"{pair['B_clone']} ({pair['isoB']})")
    for k, nd in enumerate(agg["n_donors"]):
        ax.text(k, -8, f"n={int(nd)}", ha="center", fontsize=5.5, color="#888888")
    ax.set_xticks(x)
    ax.set_xticklabels(agg.index, rotation=45, ha="right", fontsize=(7 if level == "lineage" else 5))
    ax.set_ylabel("% isoform usage", fontsize=8)
    ax.set_title(f"{pair['gene']} — {pair['A_clone']} vs {pair['B_clone']} | {level}", fontsize=8)
    ax.legend(fontsize=6, frameon=False, loc="upper right")
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def state_usage_from_pseudobulk(adata, ensg, iso_a, iso_b, donors=None, min_reads=20,
                                lineage_map=None):
    """Compute the per-cell-state ``state_usage`` frame for one isoform pair from an AltAnalyze pseudobulk
    AnnData (obs_names ``<cellstate>.<sample>-isoform``; var_names ``<ENSG>:<isoform>``).

    Returns a DataFrame indexed by cell state with usage_A, usage_B, ratio_A, n_donors (+ lineage if a
    ``lineage_map`` {cell_state -> lineage} is given). ``donors`` (iterable of sample tokens) restricts to
    a cohort; None uses all samples.
    """
    import scipy.sparse as sp
    obs = np.asarray(adata.obs_names, dtype=str)
    cell = np.array([o.split(".", 1)[0] for o in obs])
    samp = np.array([o.split(".", 1)[1].rsplit("-isoform", 1)[0] if "." in o else "" for o in obs])
    keep = np.ones(len(obs), bool) if donors is None else np.array([s in set(donors) for s in samp])
    var = np.asarray(adata.var_names, dtype=str)
    gof = np.array([v.split(":", 1)[0] for v in var])
    vidx = {v: i for i, v in enumerate(var)}
    fa, fb = f"{ensg}:{iso_a}", f"{ensg}:{iso_b}"
    if fa not in vidx or fb not in vidx:
        return pd.DataFrame(columns=["usage_A", "usage_B", "ratio_A", "n_donors"])
    cols = np.where(gof == ensg)[0]
    X = adata.X
    X = (np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X))[keep]
    cell = cell[keep]
    iA, iB = vidx[fa], vidx[fb]
    rows = {}
    for st in sorted(set(cell)):
        m = (cell == st)
        tot = X[np.ix_(m, cols)].sum(1)
        ok = tot >= min_reads
        if ok.sum() < 2:
            continue
        A = X[m][ok, iA].mean()
        B = X[m][ok, iB].mean()
        if (A + B) <= 0:
            continue
        gt = X[np.ix_(m, cols)][ok].sum(1).mean()
        rows[st] = dict(usage_A=A / gt if gt else 0.0, usage_B=B / gt if gt else 0.0,
                        ratio_A=A / (A + B), n_donors=int(ok.sum()))
    df = pd.DataFrame.from_dict(rows, orient="index")
    if lineage_map is not None and len(df):
        df["lineage"] = [lineage_map.get(s, "") for s in df.index]
    return df


def render_switch_pairs(pairs, pseudobulk, out_dir, data_dir=None, donors=None,
                        lineage_map=None, lineage_order=None, top_states=6):
    """CLI-facing driver: for every pair in ``pairs`` render the bait network + M1H bar chart and the
    per-cell-state (and, if a lineage_map is given, per-lineage) usage line plots.

    pairs : DataFrame or path to a TSV with columns gene, ensg, isoA, isoB, A_clone, B_clone
            (+ optional A_m1h_call, B_m1h_call).
    pseudobulk : AnnData or path to an isoform pseudobulk h5ad.
    Returns the list of written figure paths.
    """
    import anndata as ad
    if isinstance(pairs, str):
        pairs = pd.read_csv(pairs, sep="\t", low_memory=False)
    if isinstance(pseudobulk, str):
        pseudobulk = ad.read_h5ad(pseudobulk)
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for _, row in pairs.iterrows():
        pair = {k: row[k] for k in ("gene", "ensg", "isoA", "isoB", "A_clone", "B_clone") if k in row}
        for k in ("A_m1h_call", "B_m1h_call"):
            pair[k] = row[k] if k in row and pd.notna(row[k]) else ""
        su = state_usage_from_pseudobulk(pseudobulk, pair["ensg"], pair["isoA"], pair["isoB"],
                                         donors=donors, lineage_map=lineage_map)
        if len(su) < 2:
            logger.warning("[switch_pair_plots] %s %s vs %s: <2 usable states, skipping",
                           pair["gene"], pair["A_clone"], pair["B_clone"])
            continue
        tag = f"{pair['gene']}_{pair['A_clone']}_vs_{pair['B_clone']}"
        p = bait_network_with_m1h(pair, su, os.path.join(out_dir, f"NET_{tag}.pdf"),
                                  data_dir=data_dir, top_states=top_states)
        if p:
            written.append(p)
        written.append(usage_line(pair, su, os.path.join(out_dir, f"LINE_{tag}_cellstate.pdf"),
                                  level="cellstate"))
        if lineage_map is not None:
            written.append(usage_line(pair, su, os.path.join(out_dir, f"LINE_{tag}_lineage.pdf"),
                                      level="lineage", lineage_order=lineage_order))
    logger.info("[switch_pair_plots] rendered %d figures for %d pairs -> %s",
                len(written), len(pairs), out_dir)
    return written


def pair_umap(pair, coords, expr_A, expr_B, out_path, dot=1.6, vmax=2.0):
    """Paired UMAP: isoform A and B side by side on a shared embedding, fixed 0-``vmax`` grey->red scale.
    ``coords`` is (n_cells, 2); ``expr_A``/``expr_B`` are per-cell expression vectors. Writes ``out_path``."""
    plt = _mpl()
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("gr", ["#d9d9d9", "#fee0d2", "#fc9272", "#ef3b2c", "#a50f15"])
    coords = np.asarray(coords, float)
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.4), squeeze=False)
    for ax, (e, clone, iso, call) in zip(axes.ravel(),
                                         [(np.asarray(expr_A), pair["A_clone"], pair["isoA"], pair.get("A_m1h_call", "")),
                                          (np.asarray(expr_B), pair["B_clone"], pair["isoB"], pair.get("B_m1h_call", ""))]):
        o = np.argsort(e)
        sc = ax.scatter(coords[o, 0], coords[o, 1], c=e[o], cmap=cmap, vmin=0, vmax=vmax, s=dot,
                        linewidths=0, clip_on=False)
        ax.set_title(f"{clone} [{call}]\n{iso}\n{int(e.sum()):,} reads, {int((e > 0).sum())} cells", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        for s_ in ax.spines.values():
            s_.set_visible(False)
        cb = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(f"expr (fixed 0-{vmax:g})", fontsize=7)
        cb.outline.set_visible(False)
    fig.suptitle(f"{pair['gene']} switch pair: {pair['A_clone']} vs {pair['B_clone']}", fontsize=11)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out_path
