"""General cross-cell-state TF-isoform interaction-activity scoring.

For every catalogued TF isoform and every cell state, computes an interaction-activity score
  A(isoform, state) = sum over the isoform's edges of  w_ij * sqrt(P_i * P_j)
where P_i = the isoform's expression (best_ENST TPM) and P_j = the target gene's expression, averaged
over the donors of a chosen group ("young" by default), with a tunable expression threshold below which a
molecule does not contribute. It then ranks, for each isoform, the cell-state PAIR with the largest
activity difference (where its interaction network changes most), and ranks all TFs genome-wide by that
spread. Used to find where MAX activity (and any TF's) differs most across cell states.

Outputs a (TF isoform x cell state) activity matrix, a ranked differential table, and a grey->red heatmap.
"""

import os
import argparse
import logging

import numpy as np
import pandas as pd

from . import coexpr

logger = logging.getLogger(__name__)


def young_libraries(dataset_dir, group="young", meta=None):
    """Donor keys whose covariate group == `group`. Returns BOTH the `uid` and `library` values, since
    pseudobulk obs donor keys may use either form across samples."""
    m = pd.read_csv(os.path.join(dataset_dir, meta or coexpr.META), sep="\t", dtype=str).fillna("")
    sel = m[m["groups"].str.lower() == group.lower()]
    libs = set(sel["library"]) | (set(sel["uid"]) if "uid" in sel.columns else set())
    logger.info("[cross_state] %s donor keys (uid|library): %s", group, sorted(libs))
    return libs


def state_activity(dataset_dir, group="young", min_expr=1.0, min_state_donors=1):
    """(TF isoform x cell state) interaction-activity matrix, averaged over `group` donors, with an
    expression threshold (molecules below `min_expr` TPM do not contribute to an edge)."""
    edges = coexpr.load_edges()
    sym2ensg = coexpr.load_symbol_to_ensg(edges)
    keep = young_libraries(dataset_dir, group)
    iso_tpm, _, _ = coexpr.load_groups(dataset_dir, states=None, keep_donors=keep)      # every cell state
    max_iso = coexpr.load_max_isoform_expr(dataset_dir, states=None, keep_donors=keep)  # max-isoform CPM

    def state_mean(tpm):
        sub = tpm[tpm.index.get_level_values("donor").isin(keep)]
        cnt = sub.groupby(level="state").size()
        good = cnt[cnt >= min_state_donors].index
        return sub.groupby(level="state").mean().loc[good]

    iso_s, gene_s = state_mean(iso_tpm), state_mean(max_iso)
    states = list(iso_s.index)
    logger.info("[cross_state] %d cell states with >=%d %s donor(s)", len(states), min_state_donors, group)

    e = edges.copy()
    e["tgt_ensg"] = e["target"].map(sym2ensg).fillna("")
    e = e[(e["best_ENST"] != "") & e["best_ENST"].isin(iso_s.columns)
          & e["tgt_ensg"].isin(gene_s.columns)].drop_duplicates(["best_ENST", "tgt_ensg", "interaction_type"])

    Pi = iso_s[e["best_ENST"].values].to_numpy().T            # (edges x states)
    Pj = gene_s[e["tgt_ensg"].values].to_numpy().T
    Pi = np.where(Pi >= min_expr, Pi, 0.0)
    Pj = np.where(Pj >= min_expr, Pj, 0.0)
    EA = e["activity_score"].to_numpy()[:, None] * np.sqrt(Pi * Pj)

    a = pd.DataFrame(EA, columns=states)
    a["Symbol"], a["isoform"] = e["Symbol"].values, e["source_isoform_id"].values
    activity = a.groupby(["Symbol", "isoform"])[states].sum()
    logger.info("[cross_state] activity matrix: %d TF isoforms x %d states (%d edges, min_expr=%g)",
                activity.shape[0], activity.shape[1], len(e), min_expr)
    return activity


def aged_vs_young(dataset_dir, group_a="aged", group_b="young", min_expr=1.0, min_edges=3,
                  iso_tpm=None, max_iso=None):
    """Compare each TF's interaction network between two covariate groups WITHIN every matched cell state
    (states present in both groups). Loads both groups' donors once. Returns two ranked tables:
    (1) activity change  A_aged - A_young per (TF isoform, cell state); (2) rewiring (jaccard divergence
    of supported edge sets) per (gene, cell state) with isoform-switch flag. Both contrast aging while
    holding cell state fixed."""
    import itertools
    edges = coexpr.load_edges()
    sym2ensg = coexpr.load_symbol_to_ensg(edges)
    da = young_libraries(dataset_dir, group_a)
    db = young_libraries(dataset_dir, group_b)
    if iso_tpm is None:                       # else caller pre-loaded a matrix covering both groups
        keep = da | db
        iso_tpm, _, _ = coexpr.load_groups(dataset_dir, states=None, keep_donors=keep)
        max_iso = coexpr.load_max_isoform_expr(dataset_dir, states=None, keep_donors=keep)

    def state_mean(tpm, donors):
        sub = tpm[tpm.index.get_level_values("donor").isin(donors)]
        return sub.groupby(level="state").mean()

    iso_a, iso_b = state_mean(iso_tpm, da), state_mean(iso_tpm, db)
    gene_a, gene_b = state_mean(max_iso, da), state_mean(max_iso, db)
    matched = sorted(set(iso_a.index) & set(iso_b.index) & set(gene_a.index) & set(gene_b.index))
    e = edges.copy()
    e["tgt_ensg"] = e["target"].map(sym2ensg).fillna("")
    e = e[(e["best_ENST"] != "") & e["best_ENST"].isin(iso_a.columns) & e["tgt_ensg"].isin(gene_a.columns)] \
        .drop_duplicates(["Symbol", "source_isoform_id", "best_ENST", "target"]).reset_index(drop=True)
    logger.info("[cross_state] aged-vs-young: %d matched cell states, %d edges", len(matched), len(e))

    def parts(iso_s, gene_s):
        Pi = iso_s.reindex(index=matched)[e["best_ENST"].values].to_numpy().T
        Pj = gene_s.reindex(index=matched)[e["tgt_ensg"].values].to_numpy().T
        EA = e["activity_score"].to_numpy()[:, None] * np.sqrt(np.where(Pi >= min_expr, Pi, 0)
                                                               * np.where(Pj >= min_expr, Pj, 0))
        sup = (Pi >= min_expr) & (Pj >= min_expr)
        return EA, sup

    EAa, supa = parts(iso_a, gene_a)
    EAb, supb = parts(iso_b, gene_b)

    def activity(EA):
        a = pd.DataFrame(EA, columns=matched)
        a["Symbol"], a["isoform"] = e["Symbol"].values, e["source_isoform_id"].values
        return a.groupby(["Symbol", "isoform"])[matched].sum()
    Aa, Ab = activity(EAa), activity(EAb)
    diff = (Aa - Ab).stack().rename("activity_change").reset_index()
    diff.columns = ["Symbol", "isoform", "cell_state", "activity_change"]
    diff["activity_%s" % group_a] = Aa.stack().values
    diff["activity_%s" % group_b] = Ab.stack().values
    diff = diff[diff[["activity_%s" % group_a, "activity_%s" % group_b]].max(1) >= 1.0]
    diff = diff.reindex(diff["activity_change"].abs().sort_values(ascending=False).index).reset_index(drop=True)

    # rewiring per (gene, matched state): edge sets a vs b
    rows = []
    for si, st in enumerate(matched):
        ga, gb = {}, {}
        for i in range(len(e)):
            g = e.at[i, "Symbol"]
            it = (e.at[i, "source_isoform_id"], e.at[i, "target"])
            if supa[i, si]:
                ga.setdefault(g, set()).add(it)
            if supb[i, si]:
                gb.setdefault(g, set()).add(it)
        for g in set(ga) | set(gb):
            A, B = ga.get(g, set()), gb.get(g, set())
            if len(A) < min_edges and len(B) < min_edges:
                continue
            u = A | B
            jac = 1 - len(A & B) / len(u) if u else 0.0
            doma = _dominant_iso(A); domb = _dominant_iso(B)
            rows.append((g, st, len(B), len(A), len(A & B), round(jac, 3), domb, doma, doma != domb))
    rew = pd.DataFrame(rows, columns=["gene", "cell_state", "edges_%s" % group_b, "edges_%s" % group_a,
                                      "shared", "jaccard_divergence", "dominant_%s" % group_b,
                                      "dominant_%s" % group_a, "isoform_switch"])
    rew = rew.sort_values("jaccard_divergence", ascending=False).reset_index(drop=True)
    return diff, rew


def _dominant_iso(edge_set):
    if not edge_set:
        return ""
    c = {}
    for iso, _ in edge_set:
        c[iso] = c.get(iso, 0) + 1
    return max(c, key=c.get)


def rank_differential(activity, min_activity=1.0):
    """Per TF isoform: the cell-state pair with the largest activity difference, ranked by that spread."""
    vals = activity.to_numpy()
    amax, amin = vals.max(1), vals.min(1)
    smax = activity.columns[vals.argmax(1)]
    smin = activity.columns[vals.argmin(1)]
    out = activity.reset_index()[["Symbol", "isoform"]].copy()
    out["state_high"], out["state_low"] = smax, smin
    out["activity_high"], out["activity_low"] = amax, amin
    out["activity_range"] = amax - amin
    out = out[out["activity_high"] >= min_activity].sort_values("activity_range", ascending=False)
    return out.reset_index(drop=True)


def heatmap(activity, ranking, out_path, top_tfs=30, top_states=24):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    plt.rcParams["font.family"] = "sans-serif"; plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42
    grey_red = LinearSegmentedColormap.from_list("white_red", ["#FFFFFF", "#FF0000"])

    rows = ranking.head(top_tfs)
    sub = activity.loc[[(r.Symbol, r.isoform) for r in rows.itertuples()]]
    # show the most-involved states (highest variance across the selected isoforms)
    sc = sub.var(0).sort_values(ascending=False).head(top_states).index
    sub = sub[sorted(sc)]
    M = np.log1p(sub.to_numpy())
    fig, ax = plt.subplots(figsize=(max(7, 0.42 * sub.shape[1]), max(5, 0.3 * sub.shape[0])))
    im = ax.imshow(M, cmap=grey_red, aspect="auto", norm=Normalize(0, M.max()))
    ax.set_xticks(range(sub.shape[1])); ax.set_xticklabels(sub.columns, rotation=90, fontsize=6)
    ax.set_yticks(range(sub.shape[0]))
    ax.set_yticklabels(["%s %s" % (i[0], i[1]) for i in sub.index], fontsize=6)
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01); cb.set_label("log1p interaction activity")
    ax.set_title("Top cross-cell-state differential TF-isoform activity (%s donors)" % "young", fontsize=10)
    fig.tight_layout(); fig.savefig(out_path); plt.close(fig)
    logger.info("[cross_state] wrote %s", out_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--group", default="young")
    ap.add_argument("--min-expr", type=float, default=1.0)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    os.makedirs(a.out, exist_ok=True)
    activity = state_activity(a.dataset, a.group, a.min_expr)
    ranking = rank_differential(activity)
    activity.to_csv(os.path.join(a.out, "tf_isoform_activity_by_state.csv"))
    ranking.to_csv(os.path.join(a.out, "tf_differential_ranking.csv"), index=False)
    heatmap(activity, ranking, os.path.join(a.out, "fig_cross_state_tf_activity.pdf"))
    mx = ranking[ranking["Symbol"] == "MAX"]
    logger.info("[cross_state] MAX cell-state activity differences:\n%s",
                mx[["isoform", "state_high", "state_low", "activity_high", "activity_low", "activity_range"]]
                .head(6).to_string(index=False))
    logger.info("[cross_state] top differential TFs:\n%s",
                ranking[["Symbol", "isoform", "state_high", "state_low", "activity_range"]].head(15)
                .to_string(index=False))


if __name__ == "__main__":
    main()
