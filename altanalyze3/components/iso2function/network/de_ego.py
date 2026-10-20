"""Per-cell-state ego-networks for a regulator, coloured by DIFFERENTIAL EXPRESSION (state vs all other
states) instead of raw abundance. Reuses the fixed layout and edge logic of :mod:`ego_plot` so panels are
directly comparable, but a node is coloured only when it is significantly differential in that state
(FDR < threshold): red = up in the state, blue = down; grey = not differential / not measured. An edge is
drawn only when BOTH of its endpoint nodes are differential in that state (a co-differential interaction).

DE values come from precomputed AltAnalyze cell-state stats files (``diff-cluster-gene`` for targets,
``diff-cluster-isoform`` for the regulator's isoforms). The caller supplies them as nested dicts keyed by
state so this module performs no file IO of its own.
"""

import math
import logging

import numpy as np

from . import coexpr
from .ego_plot import _draw_edge

logger = logging.getLogger(__name__)

UP_COLOR = "#B2182B"        # up in this cell state vs all others
DOWN_COLOR = "#2166AC"      # down in this cell state vs all others
NEUTRAL_COLOR = "#D9D9D9"   # not differential (FDR >= thresh) or not measured
ISO_COLORS = ["#D1495B", "#2E86AB", "#3C6E47", "#B5651D", "#7A4FA3", "#8A8D00"]


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    return plt


def load_final_isoform_map(regulator, data_dir=None):
    """Return {gff_clone_label -> (clone_id, final_isoform_id)} from clone_to_structure.tsv — the curated
    1-to-1 TFIsoDB-clone -> KINNEX-isoform identity (final_isoform_id is a known ENST or a novel
    <num>.<sample> token, and is UNIQUE per clone). Clones with no assigned final_isoform_id are omitted."""
    import os
    import pandas as pd
    if data_dir is None:
        data_dir = os.path.join(os.path.dirname(__file__), "..", "data")
    ct = pd.read_csv(os.path.join(data_dir, "clone_to_structure.tsv"), sep="\t", dtype=str)
    ct = ct[ct["gene_symbol"] == regulator]
    out = {}
    for r in ct.itertuples(index=False):
        fid = r.final_isoform_id
        if isinstance(fid, str) and fid and fid.lower() != "nan":
            out[r.gff_clone_label] = (r.clone_id, fid)     # fid used verbatim (versionless ENST or novel id)
    return out


def build_regulator_structure(regulator, iso_universe, edges=None, sym2ensg=None,
                              collapse_redundant_clones=True, key_by_final_isoform=False,
                              data_dir=None, isoform_links=None):
    """Return (isos, all_targets, sym2ensg) for `regulator`: the same isoform/target structure and target
    ordering used by ego_plot.plot_regulator_states, but built from `iso_universe` (the set of isoform
    tokens actually measured) with NO expression threshold. `iso_universe` is a set of versionless ENST
    ids. `isos` is a list of (node_label, [ensts], {target_symbol: (interaction_type, tgt_ensg)}).

    `key_by_final_isoform` (RECOMMENDED for the DE network): key each node on the clone's curated
    `final_isoform_id` (clone_to_structure.tsv) — a strict 1-to-1 TFIso<->KINNEX mapping, so every KINNEX
    isoform appears exactly once and `best_ENST` many-to-one collisions (which showed one transcript as
    several diamonds) are impossible. Clones whose final_isoform_id is not in `iso_universe` (novel IDs
    absent from this catalog, or unassigned) are DROPPED as unmappable ("disagreements"). Node label is the
    clone's `gff_clone_label`.

    ``isoform_links`` optionally supplies explicit clone/transcript links with ``gff_clone_label``,
    ``clone_id`` and ``final_isoform_id`` columns. Each linked transcript has its own node and DE test;
    unmeasured links remain unavailable. The caller must verify biological mapping evidence before use.

    When `collapse_redundant_clones` is True (best_ENST mode only), clones with an IDENTICAL best_ENST set
    are merged; this is discouraged because such clones can have DISTINCT targets."""
    if edges is None:
        edges = coexpr.load_edges()
    if sym2ensg is None:
        sym2ensg = coexpr.load_symbol_to_ensg(edges)
    me = edges[edges["Symbol"] == regulator].copy()
    me["tgt_ensg"] = me["target"].map(sym2ensg).fillna("")
    isos = []
    if isoform_links is not None:
        # Explicit, source-verified clone -> current transcript links. Each transcript remains
        # a separate expression test; never pool same-protein transcripts or choose one for them.
        required = {"gff_clone_label", "clone_id", "final_isoform_id"}
        if not required.issubset(isoform_links.columns):
            raise ValueError("isoform_links requires columns: " + ", ".join(sorted(required)))
        links = isoform_links.copy()
        if links[list(required)].isna().any().any() or links[list(required)].eq("").any().any():
            raise ValueError("Explicit clone links must have nonempty identities")
        if links.duplicated(["gff_clone_label", "final_isoform_id"]).any():
            raise ValueError("Duplicate explicit clone/transcript link")
        if not set(links.gff_clone_label).issubset(set(me.source_isoform_id)):
            raise ValueError("Explicit clone link has no source interaction record")
        for sid, grp in me.groupby("source_isoform_id"):
            tgts = {}
            for r in grp.itertuples(index=False):
                if r.tgt_ensg:
                    tgts.setdefault(r.target, (r.interaction_type, r.tgt_ensg))
            for r in links[links.gff_clone_label == sid].itertuples(index=False):
                # Keep catalog links without expression rows as unavailable, not negative tests.
                if tgts:
                    isos.append((r.clone_id + "|" + r.final_isoform_id,
                                 [r.final_isoform_id], dict(tgts)))
        isos.sort(key=lambda x: (-len(x[2]), x[0]))
        return isos, _order_targets(isos), sym2ensg
    if key_by_final_isoform:
        fmap = load_final_isoform_map(regulator, data_dir=data_dir)
        for sid, grp in me.groupby("source_isoform_id"):
            if sid not in fmap:
                continue
            _clone_id, fid = fmap[sid]
            if fid not in iso_universe:               # unmappable to this KINNEX catalog -> drop (disagreement)
                logger.info("[de-ego] %s drop unmapped clone %s -> %s (not in iso_universe)", regulator, sid, fid)
                continue
            tgts = {}
            for r in grp.itertuples(index=False):
                if r.tgt_ensg:
                    tgts.setdefault(r.target, (r.interaction_type, r.tgt_ensg))
            if tgts:
                isos.append((sid, [fid], tgts))       # 1-to-1: exactly one KINNEX isoform per node
        isos.sort(key=lambda x: -len(x[2]))
        all_targets = _order_targets(isos)
        return isos, all_targets, sym2ensg
    for sid, grp in me.groupby("source_isoform_id"):
        ensts = sorted({e.split(".")[0] for e in grp["best_ENST"]
                        if e and e.split(".")[0] in iso_universe})
        if not ensts:
            continue
        tgts = {}
        for r in grp.itertuples(index=False):
            if r.tgt_ensg:
                tgts.setdefault(r.target, (r.interaction_type, r.tgt_ensg))
        if tgts:
            isos.append((sid, ensts, tgts))
    if collapse_redundant_clones:
        groups = {}
        for sid, ensts, tgts in isos:
            groups.setdefault(frozenset(ensts), []).append((sid, ensts, tgts))
        merged = []
        for _, members in groups.items():
            if len(members) == 1:
                merged.append(members[0]); continue
            sids = [m[0] for m in members]
            ensts = members[0][1]
            tg = {}
            for _, _, t in members:
                for g, v in t.items():
                    tg.setdefault(g, v)
            gene = sids[0].split("|", 1)[0]
            wells = [s.split("|")[-1] for s in sids]
            rep = "%s|%s" % (gene, "/".join(wells))
            merged.append((rep, ensts, tg))
            logger.info("[de-ego] %s collapsed redundant clones %s -> %s", regulator, sids, rep)
        isos = merged
    isos.sort(key=lambda x: -len(x[2]))
    all_targets = _order_targets(isos)
    return isos, all_targets, sym2ensg


def _order_targets(isos):
    """Target left-to-right ordering by isoform-connectivity signature (iso1-unique, {1,2}, iso2-unique,
    ...; targets shared by >2 isoforms go far right). Recomputed after any pruning so the layout is tight."""
    iso_of = {}
    for i, (sid, ensts, tgts) in enumerate(isos):
        for g in tgts:
            iso_of.setdefault(g, set()).add(i + 1)

    def _sig_rank(g):
        s = tuple(sorted(iso_of.get(g, ())))
        if len(s) >= 3:
            return (10 ** 9, 0, g)
        k = s[-1] if s else 10 ** 9
        return (k, (s[0] if len(s) == 2 else k), g)

    return sorted({g for _, _, t in isos for g in t}, key=_sig_rank)


def _iso_de(ide_state, ensts, thresh):
    """Aggregate isoform-level DE for a clone's constituent ENSTs in one state: pick the most-significant
    (smallest adjp) row present. Returns (signed_diff, adjp) or None if none of the ENSTs was measured."""
    best = None
    for e in ensts:
        v = ide_state.get(e)
        if v is None:
            continue
        if best is None or v[1] < best[1]:
            best = v
    return best


def plot_regulator_states_de(regulator, states, out_path, gde, ide, iso_universe,
                             edges=None, sym2ensg=None, sig_thresh=0.05, min_log2fc=0.0,
                             direction="up", title_suffix=None, restrict_to_differential=False,
                             collapse_redundant_clones=False, key_by_final_isoform=True, iso_order=None,
                             data_dir=None, unique_connections_only=False, figure_note=None,
                             isoform_links=None, expressed_isoforms=None, effect_label="log2FC",
                             eligible_isoforms=None, isoform_annotations=None):
    """One figure, one panel per state, fixed layout. `gde[state][ensg] = (signed_diff, adjp)` (targets),
    `ide[state][enst] = (signed_diff, adjp)` (regulator isoforms). A node is DIFFERENTIAL in a state when
    ``adjp < sig_thresh`` AND its signed log2 fold-change passes ``min_log2fc`` in ``direction``
    ('up' => diff > min_log2fc; 'down' => diff < -min_log2fc; 'both' => |diff| > min_log2fc). Differential
    nodes are coloured by direction (red up / blue down); others grey. Edges are drawn only when both
    endpoints are differential in that panel. When `restrict_to_differential` is True, prune the network to
    only the isoforms/targets that form at least one differential edge in >=1 shown state before layout.
    ``unique_connections_only`` shows a qualifying partner only when exactly one mapped isoform
    has a qualifying edge to it IN THAT STATE. Competitors are evaluated before display pruning;
    statistical inputs are unchanged. ``data_dir`` selects the curated clone map, and
    ``figure_note`` adds a provenance/method legend. ``isoform_links`` passes independently verified
    explicit clone/transcript links to the structure builder without pooling their DE tests.
    ``expressed_isoforms[state]`` optionally supplies independently verified expression eligibility
    instead of requiring regulator DE; gene DE gates remain unchanged. These TF nodes are gold,
    explicitly distinguished from gene DE colours. No DE values are replaced or manufactured.
    ``effect_label`` names the supplied effect scale in the legend; it does not transform values.
    ``eligible_isoforms`` optionally intersects either regulator rule with a verified per-state RNA
    eligibility set. ``isoform_annotations[state][transcript]`` adds measured support to node labels.
    Returns per-state DE counts."""
    plt = _mpl()
    from matplotlib.patches import Circle, RegularPolygon
    from matplotlib.lines import Line2D
    if edges is None:
        edges = coexpr.load_edges()
    isos, all_targets, sym2ensg = build_regulator_structure(
        regulator, iso_universe, edges, sym2ensg, collapse_redundant_clones=collapse_redundant_clones,
        key_by_final_isoform=key_by_final_isoform, data_dir=data_dir, isoform_links=isoform_links)

    def _sig(de):
        if de is None or de[1] is None or np.isnan(de[1]) or de[1] >= sig_thresh:
            return False
        d = de[0]
        if d is None or (isinstance(d, float) and np.isnan(d)):
            return False
        if direction == "up":
            return d > min_log2fc
        if direction == "down":
            return d < -min_log2fc
        return abs(d) > min_log2fc

    def _iso_pick(ide_state, ensts):
        """Isoform-NODE evaluation: a node is differential iff ANY of its constituent transcripts passes
        _sig (not merely the most-significant one, which could fail the fold gate while another passes).
        Returns (representative_de, sig_bool); representative = passing transcript with the largest
        |log2FC| (for colour), else the most-significant transcript for a grey node."""
        des = [ide_state.get(e) for e in ensts]
        des = [d for d in des if d is not None]
        passers = [d for d in des if _sig(d)]
        if passers:
            return max(passers, key=lambda d: abs(d[0])), True
        if des:
            return min(des, key=lambda d: (d[1] if (d[1] is not None and not np.isnan(d[1])) else 1.0)), False
        return None, False

    def _iso_eligible(st, ensts):
        if eligible_isoforms is not None and not any(e in eligible_isoforms.get(st, set()) for e in ensts):
            return False
        if expressed_isoforms is None:
            return _iso_pick(ide.get(st, {}), ensts)[1]
        return any(e in expressed_isoforms.get(st, set()) for e in ensts)

    # Evaluate the complete mapped network before any display pruning, so pruning cannot
    # manufacture uniqueness. The default reproduces the original endpoint significance rule.
    eligible_pairs = {}
    for st in states:
        passing = {(sid, g) for sid, ensts, tgts in isos
                   if _iso_eligible(st, ensts)
                   for g, (_, ge) in tgts.items() if _sig(gde.get(st, {}).get(ge))}
        if unique_connections_only:
            from collections import Counter
            competitors = Counter(g for _, g in passing)
            passing = {(sid, g) for sid, g in passing if competitors[g] == 1}
        eligible_pairs[st] = passing

    if restrict_to_differential:
        n0i, n0t = len(isos), len(all_targets)
        keep_pairs = set().union(*eligible_pairs.values())
        isos = [(sid, ensts, {g: v for g, v in tgts.items() if (sid, g) in keep_pairs})
                for sid, ensts, tgts in isos]
        isos = [x for x in isos if x[2]]
        all_targets = _order_targets(isos)
        logger.info("[de-ego] %s restrict: isoforms %d->%d, targets %d->%d (kept only nodes in a "
                    "differential edge)", regulator, n0i, len(isos), n0t, len(all_targets))

    if iso_order:                       # user-specified left-to-right isoform-node order (by node label)
        rank = {s: i for i, s in enumerate(iso_order)}
        isos.sort(key=lambda x: (rank.get(x[0], len(iso_order)), x[0]))
        all_targets = _order_targets(isos)
        logger.info("[de-ego] %s isoform node order: %s", regulator, [s for s, _, _ in isos])
    tgt_ensg = {g: next(v[1] for _, _, t in isos for k, v in t.items() if k == g) for g in all_targets}
    logger.info("[de-ego] %s: %d isoforms, %d targets, %d states", regulator, len(isos), len(all_targets), len(states))

    # ---- fixed layout (identical to ego_plot.plot_regulator_states) ----
    pos = {}
    nri = len(isos)
    N = len(all_targets)
    Rx = max(1.6, 0.085 * N)
    if isoform_annotations is not None:
        # Full transcript IDs and measured support need physical room between RNA nodes.
        Rx = max(Rx, 0.65 * nri)
    Ry = 1.5
    y0 = -1.0
    for i, (sid, _, _) in enumerate(isos):
        span = max(0.9, Rx * 0.85) if isoform_links is not None else 0.9
        pos[("ISO", sid)] = (0.0 if nri == 1 else -span + 2 * span * i / (nri - 1), y0)
    for j, g in enumerate(all_targets):
        ang = np.radians(np.linspace(178, 2, N)[j] if N > 1 else 90)
        pos[("TGT", g)] = (Rx * np.cos(ang), y0 + Ry * np.sin(ang))

    def node_eval(st, key):
        """Return (representative_de, sig_bool) for a node in a state."""
        if key[0] == "ISO":
            ensts = next(e for s, e, _ in isos if s == key[1])
            de, sig = _iso_pick(ide.get(st, {}), ensts)
            return de, _iso_eligible(st, ensts)
        de = gde.get(st, {}).get(tgt_ensg.get(key[1]))
        return de, _sig(de)

    def facecolor(de, sig):
        if not sig:
            return NEUTRAL_COLOR
        return UP_COLOR if de[0] > 0 else DOWN_COLOR

    r_tgt, r_iso = 0.075, 0.115
    ncol = min(3, len(states))
    nrow = math.ceil(len(states) / ncol)
    fig_width, fig_height = 2.1 * (Rx + 0.5) * ncol, (Ry + 1.9) * nrow
    if isoform_links is not None:
        fig_width, fig_height = max(9.0, fig_width), max(4.8, fig_height)
    fig, axes = plt.subplots(nrow, ncol, figsize=(fig_width, fig_height), squeeze=False)
    clone_colors = {}
    if isoform_links is not None:
        clone_colors = {c: ISO_COLORS[i % len(ISO_COLORS)]
                        for i, c in enumerate(sorted(set(isoform_links.clone_id)))}
    axflat = list(axes.flat)
    counts = {}
    for ax, st in zip(axflat, states):
        ev = {key: node_eval(st, key) for key in pos}
        de_of = {k: v[0] for k, v in ev.items()}
        sig_of = {k: v[1] for k, v in ev.items()}
        nup = sum(1 for k, d in de_of.items() if sig_of[k] and d is not None and d[0] > 0 and (expressed_isoforms is None or k[0]!='ISO'))
        ndn = sum(1 for k, d in de_of.items() if sig_of[k] and d is not None and d[0] < 0 and (expressed_isoforms is None or k[0]!='ISO'))
        counts[st] = (nup, ndn)
        rad = {key: (r_iso if key[0] == "ISO" else r_tgt) for key in pos}
        # edges: only when BOTH endpoints differential in this state
        for i, (sid, ensts, tgts) in enumerate(isos):
            col = (clone_colors[sid.split("|", 1)[0]] if isoform_links is not None
                   else ISO_COLORS[i % len(ISO_COLORS)])
            for g, (typ, _) in tgts.items():
                if (sid, g) not in eligible_pairs[st]:
                    continue
                _draw_edge(ax, pos[("ISO", sid)], pos[("TGT", g)], col, 1.1, directed=(typ == "PDI"),
                           shrinkA=rad[("ISO", sid)] + 0.01, shrinkB=rad[("TGT", g)] + 0.025)
        for key, (x, y) in pos.items():
            de = de_of[key]
            sig = sig_of[key]
            fc = "#E3C16F" if expressed_isoforms is not None and key[0]=="ISO" and sig else facecolor(de, sig)
            alpha = 1.0 if sig else 0.35
            lab_alpha = 0.95 if sig else 0.3
            r = rad[key]
            if key[0] == "ISO":
                ax.add_patch(RegularPolygon((x, y), 4, radius=r, facecolor=fc,
                                            edgecolor="#222222" if sig else "#BBBBBB",
                                            lw=1.1 if sig else 0.5, alpha=alpha, zorder=3))
                label = key[1].replace("|", "|\n")
                if isoform_links is not None:
                    import textwrap
                    clone, transcript = key[1].split("|", 1)
                    label = clone + "\n" + "\n".join(textwrap.wrap(transcript, width=23))
                    if isoform_annotations is not None:
                        extra = isoform_annotations.get(st, {}).get(transcript, "")
                        if extra:
                            label += "\n" + extra
                ax.annotate(label, (x, y - r - 0.04), fontsize=6 if isoform_annotations is not None else 4.5, ha="center",
                            va="top", alpha=lab_alpha, zorder=4)
            else:
                ax.add_patch(Circle((x, y), radius=r, facecolor=fc,
                                    edgecolor="#222222" if sig else "#BBBBBB",
                                    lw=0.5 if sig else 0.3, alpha=alpha, zorder=3))
                ax.annotate(key[1], (x, y), fontsize=6 if isoform_annotations is not None else 4, ha="center", va="center", alpha=lab_alpha,
                            color="#000000", zorder=4)
        ax.set_title("%s  (↑%d ↓%d)" % (st, nup, ndn), fontsize=8)
        if not pos:
            ax.text(0.5, 0.5, "No PPI meets both expression thresholds\nin the available tests for this state",
                    transform=ax.transAxes, ha="center", va="center", fontsize=9, color="#666666")
        ax.set_xlim(-Rx - 0.45, Rx + 0.45); ax.set_ylim(y0 - (0.65 if isoform_annotations is not None else 0.35), y0 + Ry + 0.35)
        ax.set_aspect("equal"); ax.axis("off")
    for ax in axflat:
        ax.axis("off")
    gate_txt = "FDR<%.2g & |%s|>%.2g" % (sig_thresh, effect_label, min_log2fc) if direction == "both" \
        else "FDR<%.2g & %s%s%.2g" % (sig_thresh, effect_label, ">" if direction == "up" else "<-", min_log2fc)
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=UP_COLOR, markersize=8,
                      label="up in state (%s)" % gate_txt),
               Line2D([0], [0], marker="o", color="none", markerfacecolor=DOWN_COLOR, markersize=8,
                      label="down in state (%s)" % gate_txt),
               Line2D([0], [0], marker="o", color="none", markerfacecolor=NEUTRAL_COLOR, markersize=8,
                      label="not differential" if expressed_isoforms is None else "does not meet panel rule / unavailable")]
    if expressed_isoforms is not None:
        handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=UP_COLOR,
                          markersize=8, label="partner passes upregulation rule"),
                   Line2D([0], [0], marker="o", color="none", markerfacecolor=NEUTRAL_COLOR,
                          markersize=8, label="ineligible / unavailable"),
                   Line2D([0], [0], marker="D", color="none", markerfacecolor="#E3C16F",
                          markersize=8, label="TF meets expression rule (DE not required)")]
    elif isoform_annotations is not None:
        handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=UP_COLOR,
                          markersize=8, label="passes upregulation rule"),
                   Line2D([0], [0], marker="o", color="none", markerfacecolor=NEUTRAL_COLOR,
                          markersize=8, label="ineligible / unavailable")]
    if isoform_links is not None:
        handles += [Line2D([0], [0], color=clone_colors[c], lw=2, label=c)
                    for c in sorted({sid.split("|", 1)[0] for sid, _, _ in isos})]
    else:
        handles += [Line2D([0], [0], color=ISO_COLORS[i % len(ISO_COLORS)], lw=2, label=sid)
                    for i, (sid, _, _) in enumerate(isos)]
    legend_title = ("node = DE direction; edge colour = %s isoform (drawn only if both endpoints DE)" % regulator
                    if expressed_isoforms is None else
                    "gold TF = expressed; red partner = upregulated; edge colour = assayed %s clone" % regulator)
    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 3 if expressed_isoforms is not None or isoform_annotations is not None else 6), fontsize=6,
               title=legend_title,
               frameon=False)
    suff = title_suffix or ("%d cell states" % len(states))
    if expressed_isoforms is None:
        fig.suptitle("%s differential-expression network — %s  "
                     "(node coloured only if %s in state vs all others; edges require both endpoints DE)"
                     % (regulator, suff, gate_txt), fontsize=11)
    else:
        fig.suptitle("%s — %s (TF expressed; partner gene %s vs other states)" % (regulator, suff, gate_txt), fontsize=11)
    if isoform_links is not None:
        import textwrap
        fig._suptitle.set_text(textwrap.fill(fig._suptitle.get_text(), width=105))
    if figure_note:
        import textwrap
        width, height = fig.get_size_inches()
        if unique_connections_only:
            figure_note = "Display: partner qualifies for exactly one mapped isoform within that state. " + figure_note
        wrapped_note = textwrap.fill(figure_note, width=max(90, int(width * 18)))
        note_height = len(wrapped_note.splitlines()) * 6.5 * 1.4 / 72
        total_height = height + note_height + 1.1
        fig.set_size_inches(width, total_height)
        fig.subplots_adjust(bottom=(note_height + 1.04) / total_height, top=0.84)
        for legend in fig.legends:
            legend.set_bbox_to_anchor((0.5, (note_height + 0.32) / total_height))
        fig.text(0.03, 0.14 / total_height, wrapped_note,
                 ha="left", va="bottom", fontsize=6.5, linespacing=1.4)
    fig.savefig(out_path, bbox_inches="tight"); plt.close(fig)
    logger.info("[de-ego] wrote %s", out_path)
    return counts
