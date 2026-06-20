"""Topological network DIVERGENCE (rewiring) between cell states — distinct from presence/absence.

For each gene's interaction network across cell-state pairs, this scores how much the *thresholded* edge
set is REWIRED between the two states, restricted to pairs where the network is present in BOTH states
(so "the network is just gone in one state" is NOT rewarded). Rewiring is driven especially by a switch
in the predominant TF isoform (different isoform -> different edges). To keep a gene's graph a single
connected component for topology, the gene's isoforms are joined by INVISIBLE structural edges, so
topological change is measured on a connected graph rather than disjoint per-isoform stars.

Metrics per (gene, state_a, state_b) with both networks present:
  - jaccard_divergence  = 1 - |E_a & E_b| / |E_a | E_b|   (edge-set rewiring; the headline score)
  - isoform_switch      = the predominant isoform (most supported edges) differs between the states
  - degree_divergence   = mean |deg_a(node) - deg_b(node)| over the union node set of the connected graph
An edge (isoform -> target) is "supported" in a state when BOTH endpoints are expressed >= min_expr CPM.
"""

import os
import argparse
import logging
import itertools

import numpy as np
import pandas as pd

from . import coexpr

logger = logging.getLogger(__name__)


def edge_support(dataset_dir, group="young", min_expr=1.0, min_state_donors=1):
    """Per-edge x per-state support boolean (both endpoints >= min_expr CPM), with edge metadata
    (gene, isoform, target). Source = isoform CPM (best_ENST); target = gene max-isoform CPM."""
    from .cross_state import young_libraries
    edges = coexpr.load_edges()
    sym2ensg = coexpr.load_symbol_to_ensg(edges)
    keep = young_libraries(dataset_dir, group)
    iso_tpm, _, _ = coexpr.load_groups(dataset_dir, states=None, keep_donors=keep)
    max_iso = coexpr.load_max_isoform_expr(dataset_dir, states=None, keep_donors=keep)

    def state_mean(tpm):
        sub = tpm[tpm.index.get_level_values("donor").isin(keep)]
        cnt = sub.groupby(level="state").size()
        return sub.groupby(level="state").mean().loc[cnt[cnt >= min_state_donors].index]

    iso_s, gene_s = state_mean(iso_tpm), state_mean(max_iso)
    states = list(iso_s.index)
    e = edges.copy()
    e["tgt_ensg"] = e["target"].map(sym2ensg).fillna("")
    e = e[(e["best_ENST"] != "") & e["best_ENST"].isin(iso_s.columns)
          & e["tgt_ensg"].isin(gene_s.columns)].drop_duplicates(
        ["Symbol", "source_isoform_id", "best_ENST", "target"]).reset_index(drop=True)
    Pi = iso_s[e["best_ENST"].values].to_numpy().T
    Pj = gene_s[e["tgt_ensg"].values].to_numpy().T
    support = (Pi >= min_expr) & (Pj >= min_expr)        # (edges x states)
    logger.info("[rewiring] %d edges x %d states; supported edge-state cells: %d",
                len(e), len(states), int(support.sum()))
    return e, support, states


def gene_state_edge_sets(dataset_dir, group="young", min_expr=1.0):
    """{gene -> {state -> set of (isoform, target) supported edges}} and the state list."""
    e, support, states = edge_support(dataset_dir, group, min_expr)
    sidx = {s: j for j, s in enumerate(states)}
    by_gene = {}
    for i, row in e.iterrows():
        gs = by_gene.setdefault(row["Symbol"], {s: set() for s in states})
        sup = support[i]
        for s in states:
            if sup[sidx[s]]:
                gs[s].add((row["source_isoform_id"], row["target"]))
    return by_gene, states


def select_divergent_states(by_gene, gene, k=6, min_edges=3):
    """The k cell states whose networks for `gene` are maximally divergent — greedy farthest-point on
    edge-set Jaccard distance among states where the network is present (>= min_edges supported)."""
    sd = by_gene.get(gene, {})
    present = [s for s in sd if len(sd[s]) >= min_edges]
    if len(present) <= k:
        return present

    def jd(a, b):
        u = sd[a] | sd[b]
        return 1 - len(sd[a] & sd[b]) / len(u) if u else 0.0

    start = max(itertools.combinations(present, 2), key=lambda p: jd(*p))
    sel = list(start)
    while len(sel) < k:
        sel.append(max([s for s in present if s not in sel], key=lambda s: min(jd(s, x) for x in sel)))
    return sel


def rewiring_scores(dataset_dir, group="young", min_expr=1.0, min_edges=3):
    """Rank (gene, state_a, state_b) by edge-set rewiring among pairs PRESENT IN BOTH states."""
    e, support, states = edge_support(dataset_dir, group, min_expr)
    sidx = {s: j for j, s in enumerate(states)}
    # gene -> state -> set of (isoform, target) supported edges  (+ per-isoform counts)
    by_gene = {}
    for i, row in e.iterrows():
        g = row["Symbol"]
        gs = by_gene.setdefault(g, {s: set() for s in states})
        sup = support[i]
        for s in states:
            if sup[sidx[s]]:
                gs[s].add((row["source_isoform_id"], row["target"]))

    rows = []
    for g, sd in by_gene.items():
        present = [s for s in states if len(sd[s]) >= min_edges]
        if len(present) < 2:
            continue
        for a, b in itertools.combinations(present, 2):
            Ea, Eb = sd[a], sd[b]
            union = Ea | Eb
            inter = Ea & Eb
            jac = 1 - len(inter) / len(union) if union else 0.0
            if jac == 0:
                continue                                  # identical wiring -> not divergent
            doma = _dominant_isoform(Ea); domb = _dominant_isoform(Eb)
            # degree divergence on the connected graph (nodes = isoforms+targets; isoforms joined)
            deg = _degree_divergence(Ea, Eb)
            rows.append((g, a, b, len(Ea), len(Eb), len(inter), round(jac, 3),
                         doma, domb, doma != domb, round(deg, 3)))
    out = pd.DataFrame(rows, columns=[
        "gene", "state_a", "state_b", "edges_a", "edges_b", "shared_edges", "jaccard_divergence",
        "dominant_iso_a", "dominant_iso_b", "isoform_switch", "degree_divergence"])
    out = out.sort_values(["jaccard_divergence", "degree_divergence"], ascending=False).reset_index(drop=True)
    logger.info("[rewiring] %d divergent (gene,state,state) pairs present-in-both | isoform-switch: %d",
                len(out), int(out["isoform_switch"].sum()))
    return out


def _dominant_isoform(edge_set):
    if not edge_set:
        return ""
    c = {}
    for iso, _ in edge_set:
        c[iso] = c.get(iso, 0) + 1
    return max(c, key=c.get)


def _degree_divergence(Ea, Eb):
    """mean |deg_a - deg_b| over the union of nodes; the same-gene isoforms are linked by an invisible
    edge so the graph is connected (each isoform gets +1 structural degree in both states -> cancels)."""
    def deg(E):
        d = {}
        for iso, tgt in E:
            d[iso] = d.get(iso, 0) + 1
            d[tgt] = d.get(tgt, 0) + 1
        return d
    da, db = deg(Ea), deg(Eb)
    nodes = set(da) | set(db)
    return float(np.mean([abs(da.get(n, 0) - db.get(n, 0)) for n in nodes])) if nodes else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--group", default="young")
    ap.add_argument("--min-expr", type=float, default=1.0)
    ap.add_argument("--min-edges", type=int, default=3)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    os.makedirs(a.out, exist_ok=True)
    out = rewiring_scores(a.dataset, a.group, a.min_expr, a.min_edges)
    out.to_csv(os.path.join(a.out, "network_rewiring_ranking.csv"), index=False)
    sw = out[out["isoform_switch"]]
    logger.info("[rewiring] TOP isoform-switch-driven rewiring (network present in both states):\n%s",
                sw[["gene", "state_a", "state_b", "edges_a", "edges_b", "shared_edges",
                    "jaccard_divergence", "dominant_iso_a", "dominant_iso_b"]].head(15).to_string(index=False))
    logger.info("[rewiring] MAX rewiring (present in both):\n%s",
                out[out["gene"] == "MAX"][["state_a", "state_b", "edges_a", "edges_b", "shared_edges",
                    "jaccard_divergence", "dominant_iso_a", "dominant_iso_b", "isoform_switch"]]
                .head(8).to_string(index=False))


if __name__ == "__main__":
    main()
