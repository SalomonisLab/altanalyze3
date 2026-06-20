"""Deploy the expression+assay condition-network contrast between two groups (cell-states/covariates).

Pipeline: load donor-paired pseudobulks -> per-edge activity -> donor-paired differential edges ->
gene DE -> PDI regulon coherence + direction -> driver TF isoforms -> WikiPathways pathway contrast ->
network topology summary. Writes CSV tables and editable PDF figures.

Usage:
  python -m altanalyze3.components.iso2function.network.run_contrast \
      --dataset /path/to/Isoform-Workflow-BAM --state-a "HSC-1" --state-b "Myeloid intermediate 1" \
      --out /path/to/output
"""

import os
import argparse
import logging

import numpy as np
import pandas as pd
from scipy import stats

from . import coexpr, wikipathways

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


# --------------------------------------------------------------------------- pathway contrast
def pathway_contrast(edge_diff, gene_de, paths, min_overlap=3):
    """Two complementary views per WikiPathway: (1) enrichment of differentially-active-edge genes in
    the pathway (hypergeometric), (2) sub-network differential activity (mean edge log2FC of edges with
    >=1 endpoint in the pathway, donor-agnostic effect size)."""
    # gene universe = all genes participating in the tested network (edge targets + measured genes)
    universe = set(edge_diff["target_ENSG"]) | set(gene_de.index)
    sig = edge_diff[edge_diff["differential"]]
    sig_genes = set(sig["target_ENSG"])
    rows = []
    for wp, info in paths.items():
        pg = info["genes"] & universe
        if len(pg) < min_overlap:
            continue
        hit = sig_genes & pg
        # hypergeometric enrichment of significant-edge genes in the pathway
        M, n, N = len(universe), len(sig_genes & universe), len(pg)
        k = len(hit)
        p_enrich = stats.hypergeom.sf(k - 1, M, n, N) if (n and N) else 1.0
        # sub-network activity differential: edges touching pathway genes
        sub = edge_diff[edge_diff["target_ENSG"].isin(pg)]
        d_act = float(np.nanmean(sub["log2FC"])) if len(sub) else 0.0
        rows.append((wp, info["name"], len(pg), k, p_enrich, len(sub), d_act))
    out = pd.DataFrame(rows, columns=["WPID", "pathway", "n_genes_in_net", "n_sig_genes",
                                      "p_enrichment", "n_edges", "mean_edge_log2FC"])
    if not out.empty:
        out["FDR_enrichment"] = coexpr._bh(out["p_enrichment"].to_numpy())
        out = out.sort_values("p_enrichment").reset_index(drop=True)
    return out


# --------------------------------------------------------------------------- topology
def topology_summary(edge_diff, state_a, state_b):
    """Node-level differential strength + a global summary of how different the two networks are."""
    ca, cb = "mean_EA_%s" % state_a, "mean_EA_%s" % state_b
    df = edge_diff[["Symbol", "source_isoform_id", ca, cb, "log2FC", "FDR", "differential"]].rename(
        columns={"Symbol": "node", "source_isoform_id": "isoform", ca: "EA_a", cb: "EA_b"})
    node = df.groupby(["node", "isoform"]).agg(
        strength_a=("EA_a", "sum"), strength_b=("EA_b", "sum"),
        n_edges=("log2FC", "size"), n_differential=("differential", "sum")).reset_index()
    node["delta_strength"] = np.log2(node["strength_a"] + coexpr.EPS) - np.log2(node["strength_b"] + coexpr.EPS)
    node = node.reindex(node["delta_strength"].abs().sort_values(ascending=False).index).reset_index(drop=True)
    diff = edge_diff[edge_diff["differential"]]
    glob = {
        "edges_tested": len(edge_diff),
        "edges_differential": len(diff),
        "edges_significant_rawp": int(edge_diff["significant"].sum()),
        "edges_FDR_sig": int(edge_diff["significant_FDR"].sum()),
        "gained_or_lost": int(diff["class"].str.startswith(("gained", "lost")).sum()),
        "rewired": int((diff["class"] == "rewired").sum()),
        "quantitative": int((diff["class"] == "quantitative").sum()),
        "network_distance_sumabsLFC": float(np.nansum(np.abs(diff["log2FC"]))),
    }
    return node, glob


# --------------------------------------------------------------------------- figures
def make_figures(out, state_a, state_b, edge_diff, drivers, pdi_coh, paths_tbl):
    plt = _mpl()
    C_A, C_B, C_GREY = "#D1495B", "#2E86AB", "#B0B0B0"

    # 1. volcano of differential edges (raw p; colored by effect-size differential + direction)
    fig, ax = plt.subplots(figsize=(5, 4.5))
    x = edge_diff["log2FC"].to_numpy()
    y = -np.log10(np.clip(edge_diff["p"].to_numpy(), 1e-300, 1))
    d = edge_diff["differential"].to_numpy()
    ax.scatter(x[~d], y[~d], s=6, c=C_GREY, linewidths=0, alpha=0.5)
    ax.scatter(x[d & (x > 0)], y[d & (x > 0)], s=12, c=C_A, linewidths=0)
    ax.scatter(x[d & (x < 0)], y[d & (x < 0)], s=12, c=C_B, linewidths=0)
    ax.axvline(1, color="#444444", lw=0.8); ax.axvline(-1, color="#444444", lw=0.8)
    ax.set_xlabel("edge activity log2FC (%s vs %s)" % (state_a, state_b))
    ax.set_ylabel("-log10 p (moderated paired)")
    ax.set_title("Differential interaction edges  (red=up in %s)" % state_a)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig_volcano_edges.pdf")); plt.close(fig)

    # 2. top driver isoforms
    top = drivers.head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.barh(range(len(top)), top["network_impact"], color="#3C6E47", height=0.7)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["Symbol"].fillna("") + " (" + top["source_isoform_id"].fillna("").astype(str) + ")",
                       fontsize=6)
    ax.set_xlabel("network impact score")
    ax.set_title("Driver TF isoforms (%s vs %s)" % (state_a, state_b))
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig_driver_isoforms.pdf")); plt.close(fig)

    # 3. PDI regulon direction: isoform log2FC vs target median log2FC
    if pdi_coh is not None and not pdi_coh.empty:
        fig, ax = plt.subplots(figsize=(5, 4.5))
        sig = pdi_coh["FDR_combined"] < 0.1
        ax.scatter(pdi_coh["isoform_edge_log2FC"][~sig], pdi_coh["target_median_log2FC"][~sig],
                   s=8, c=C_GREY, linewidths=0, alpha=0.6)
        ax.scatter(pdi_coh["isoform_edge_log2FC"][sig], pdi_coh["target_median_log2FC"][sig],
                   s=18, c="#7A4FA3", linewidths=0)
        ax.axhline(0, color="#444444", lw=0.8); ax.axvline(0, color="#444444", lw=0.8)
        ax.set_xlabel("TF isoform edge log2FC")
        ax.set_ylabel("target-gene median log2FC")
        ax.set_title("PDI regulatory coherence")
        for _, r in pdi_coh[sig].head(12).iterrows():
            ax.annotate(r["Symbol"], (r["isoform_edge_log2FC"], r["target_median_log2FC"]), fontsize=5)
        fig.tight_layout(); fig.savefig(os.path.join(out, "fig_pdi_coherence.pdf")); plt.close(fig)

    # 4. top pathways by enrichment
    if paths_tbl is not None and not paths_tbl.empty:
        tp = paths_tbl.head(20).iloc[::-1]
        fig, ax = plt.subplots(figsize=(6.5, 5))
        ax.barh(range(len(tp)), -np.log10(np.clip(tp["p_enrichment"], 1e-300, 1)), color="#B5651D", height=0.7)
        ax.set_yticks(range(len(tp))); ax.set_yticklabels(tp["pathway"].str.slice(0, 48), fontsize=6)
        ax.set_xlabel("-log10 enrichment p")
        ax.set_title("WikiPathways enriched for rewired edges")
        fig.tight_layout(); fig.savefig(os.path.join(out, "fig_pathways.pdf")); plt.close(fig)


# --------------------------------------------------------------------------- orchestration
def run(dataset_dir, state_a, state_b, out_dir, min_detect=2):
    os.makedirs(out_dir, exist_ok=True)
    edges = coexpr.load_edges()
    sym2ensg = coexpr.load_symbol_to_ensg(edges)
    iso_tpm, gene_tpm, donors = coexpr.load_groups(dataset_dir, [state_a, state_b])
    # isoform-resolved on BOTH sides: source = the TF's best_ENST; target = the gene's reference
    # (MANE-Select, productive/Not-NMD) transcript. Gene-level total TPM is NOT used.
    ref_iso = coexpr.reference_isoform_expr(iso_tpm)
    common = sorted(set(donors[state_a]) & set(donors[state_b]))
    logger.info("[run] donor-paired contrast %s vs %s on %d donors: %s",
                state_a, state_b, len(common), common)

    ea = coexpr.edge_activity(edges, iso_tpm, ref_iso, sym2ensg, [state_a, state_b], donors)
    edge_diff = coexpr.differential_edges(ea, state_a, state_b, common, min_detect=min_detect)
    gene_de = coexpr.gene_differential(ref_iso, state_a, state_b, common)
    m1h = coexpr.load_m1h()
    pdi_coh = coexpr.pdi_coherence(edge_diff, gene_de, state_a, state_b, m1h)
    drivers = coexpr.driver_isoforms(edge_diff, pdi_coh)
    paths = wikipathways.load_pathways()
    paths_tbl = pathway_contrast(edge_diff, gene_de, paths)
    node, glob = topology_summary(edge_diff, state_a, state_b)

    edge_diff.to_csv(os.path.join(out_dir, "differential_edges.csv"), index=False)
    gene_de.reset_index().to_csv(os.path.join(out_dir, "gene_differential.csv"), index=False)
    pdi_coh.to_csv(os.path.join(out_dir, "pdi_regulon_coherence.csv"), index=False)
    drivers.to_csv(os.path.join(out_dir, "driver_isoforms.csv"), index=False)
    paths_tbl.to_csv(os.path.join(out_dir, "pathway_contrast.csv"), index=False)
    node.to_csv(os.path.join(out_dir, "node_topology.csv"), index=False)
    pd.Series(glob).to_csv(os.path.join(out_dir, "network_summary.csv"))

    make_figures(out_dir, state_a, state_b, edge_diff, drivers, pdi_coh, paths_tbl)
    logger.info("[run] DONE. global summary: %s", glob)
    logger.info("[run] outputs in %s", out_dir)
    return {"edge_diff": edge_diff, "pdi": pdi_coh, "drivers": drivers, "paths": paths_tbl, "global": glob}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--state-a", required=True)
    ap.add_argument("--state-b", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-detect", type=int, default=2, help="min donors detecting an edge to test it")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    run(a.dataset, a.state_a, a.state_b, a.out, min_detect=a.min_detect)


if __name__ == "__main__":
    main()
