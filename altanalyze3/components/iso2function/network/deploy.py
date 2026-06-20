"""Condition interaction-network visualization step for the AltAnalyze3 long-read workflow.

Runs immediately after differential isoform analysis (``sclr-diff``) on the pseudobulk h5ads it produced,
and writes the cross-cell-state + differential interaction-network figures. Exposed three ways:
  * ``run_network_visualization(...)`` — the callable API (used by the workflow hook and remote calls),
  * ``altanalyze3 sclr-iso2func-network ...`` — a standalone CLI subcommand (submittable as a cluster job),
  * automatically inside ``sclr-diff`` when ``--network-viz`` is passed.

For each requested interaction type (all / PDI):
  1. cross-cell-state TF-isoform interaction activity + ranking + grey/white->red heatmap (reference group),
  2. network REWIRING (topological divergence) ranking + per-TF "top-6 divergent cell states" figures,
  3. for each requested group contrast (e.g. disease vs control): per-matched-cell-state activity change
     + rewiring tables and young-vs-aged-style isoform-switch figures.
All figures use the standard style (lower half-dome, signature-ordered targets, centered gene labels,
data-unit PDI arrowheads, only-connected option).
"""

import os
import logging

import pandas as pd

from . import coexpr, cross_state, rewiring, ego_plot

logger = logging.getLogger(__name__)

# default pseudobulk file names produced by the long-read workflow (sclr-diff)
DEFAULT_ISO_H5AD = "isoform_combined_pseudo_cluster_cpm-filtered.h5ad"
DEFAULT_GENE_H5AD = "gene_combined_pseudo_cluster_cpm-filtered.h5ad"


def _load_expr(dataset_dir, keep):
    iso, _, don = coexpr.load_groups(dataset_dir, None, keep_donors=keep)
    mx = coexpr.load_max_isoform_expr(dataset_dir, None, keep_donors=keep)
    return iso, mx, don


def _relabel(tpm, donor_group):
    tpm = tpm.copy()
    tpm.index = pd.MultiIndex.from_tuples(
        [(f"{s} ({donor_group.get(d, '?')})", d) for (s, d) in tpm.index], names=["state", "donor"])
    return tpm


def cross_state_block(dataset_dir, group, out, top_tfs=10, min_expr=1.0, min_edges=3, k_states=6,
                      extra_tfs=None):
    """Cross-state activity + rewiring among one reference group, with per-TF top-divergent figures.

    extra_tfs : always plot these TFs in addition to the top-ranked, even if they lack a dominant-isoform
    switch (e.g. multi-isoform regulators whose isoforms hit different targets, or single-isoform TFs with
    expression divergence). Useful when an interaction type is sparse (PDI)."""
    os.makedirs(out, exist_ok=True)
    act = cross_state.state_activity(dataset_dir, group, min_expr)
    rank = cross_state.rank_differential(act)
    act.to_csv(os.path.join(out, f"{group}_tf_activity_by_state.csv"))
    rank.to_csv(os.path.join(out, f"{group}_tf_differential_ranking.csv"), index=False)
    cross_state.heatmap(act, rank, os.path.join(out, f"fig_{group}_cross_state_tf_activity.pdf"))
    rw = rewiring.rewiring_scores(dataset_dir, group, min_expr, min_edges)
    rw.to_csv(os.path.join(out, f"{group}_network_rewiring_ranking.csv"), index=False)
    sw = rw[rw["isoform_switch"]]
    gtf = (sw.sort_values("jaccard_divergence", ascending=False).groupby("gene")
           .agg(best_div=("jaccard_divergence", "max"), n=("gene", "size")).reset_index()
           .sort_values(["best_div", "n"], ascending=False))
    gtf.to_csv(os.path.join(out, f"{group}_rewiring_ranked_TFs.csv"), index=False)

    by_gene, _ = rewiring.gene_state_edge_sets(dataset_dir, group, min_expr)
    iso, mx, don = _load_expr(dataset_dir, cross_state.young_libraries(dataset_dir, group))
    fa = os.path.join(out, f"divergent_TFs_{group}")
    fc = os.path.join(fa, "only_connected")
    os.makedirs(fc, exist_ok=True)
    plot_tfs = list(gtf["gene"].head(top_tfs))
    plot_tfs += [t for t in (extra_tfs or []) if t not in set(plot_tfs)]
    for tf in plot_tfs:
        sel = rewiring.select_divergent_states(by_gene, tf, k=k_states)
        if len(sel) < 2:
            logger.info("[iso2func-netviz] %s: <2 divergent states in %s, skipped", tf, group)
            continue
        for oc, o in ((False, fa), (True, fc)):
            ego_plot.plot_regulator_states(dataset_dir, tf, sel, os.path.join(o, f"fig_{tf}_top{k_states}divergent.pdf"),
                                           group=None, iso_tpm=iso, gene_tpm=mx, donors=don, only_connected=oc)
    logger.info("[iso2func-netviz] cross-state(%s): %d states, %d rewiring TFs, top %d plotted",
                group, act.shape[1], len(gtf), top_tfs)
    return gtf


def comparison_block(dataset_dir, group_a, group_b, out, min_expr=1.0, min_edges=3, max_cases=12):
    """Per-matched-cell-state group contrast (group_a vs group_b): activity change + rewiring + figures."""
    cname = f"{group_a}-{group_b}"
    da = cross_state.young_libraries(dataset_dir, group_a)
    db = cross_state.young_libraries(dataset_dir, group_b)
    iso, mx, _ = _load_expr(dataset_dir, da | db)
    dg = {d: group_a for d in da}
    dg.update({d: group_b for d in db})
    iso_r, mx_r = _relabel(iso, dg), _relabel(mx, dg)
    don_r = {}
    for (ps, d) in iso_r.index:
        don_r.setdefault(ps, set()).add(d)
    don_r = {k: sorted(v) for k, v in don_r.items()}

    diff, rew = cross_state.aged_vs_young(dataset_dir, group_a, group_b, min_expr, min_edges,
                                          iso_tpm=iso, max_iso=mx)
    diff.to_csv(os.path.join(out, f"{cname}_activity_change_per_state.csv"), index=False)
    rew.to_csv(os.path.join(out, f"{cname}_rewiring_per_state.csv"), index=False)
    fa = os.path.join(out, f"divergent_TFs_{cname}")
    fc = os.path.join(fa, "only_connected")
    os.makedirs(fc, exist_ok=True)
    ec = [c for c in rew.columns if c.startswith("edges_")]
    sw = rew[rew["isoform_switch"] & (rew[ec].min(1) >= min_edges)]
    cases = sw.sort_values("jaccard_divergence", ascending=False).head(max_cases).drop_duplicates(["gene", "cell_state"])
    for r in cases.itertuples(index=False):
        st = r.cell_state
        for oc, o in ((False, fa), (True, fc)):
            ego_plot.plot_regulator_states(dataset_dir, r.gene, [f"{st} ({group_a})", f"{st} ({group_b})"],
                                           os.path.join(o, f"fig_{r.gene}_{st}_{cname}.pdf".replace(" ", "_")),
                                           group=None, iso_tpm=iso_r, gene_tpm=mx_r, donors=don_r, only_connected=oc)
    logger.info("[iso2func-netviz] contrast %s: %d rewiring pairs, %d figures", cname, len(rew), len(cases))


def run_network_visualization(dataset_dir, out_dir=None, metadata=None, iso_h5ad=None, gene_h5ad=None,
                              reference_group="young", comparisons=None, edge_types=("all", "PDI"),
                              top_tfs=10, min_expr=1.0, min_edges=3, extra_tfs=None):
    """Generate all condition interaction-network figures for a long-read analysis result directory.

    dataset_dir : directory holding the pseudobulk h5ads + metadata (the sclr-diff work dir).
    comparisons : list of (group_a, group_b) covariate contrasts (e.g. [("AML-SRSF2","aged")]); optional.
    edge_types  : any of "all"/"PDI"/"PPI"; each gets its own sub-output (PDI -> PDI_only/).
    """
    out_dir = out_dir or os.path.join(dataset_dir, "iso2function_network")
    coexpr.ISO_H5AD = iso_h5ad or DEFAULT_ISO_H5AD
    coexpr.GENE_H5AD = gene_h5ad or DEFAULT_GENE_H5AD
    if metadata:
        coexpr.META = metadata
    for et in edge_types:
        coexpr.EDGE_TYPE = None if et in (None, "all") else et
        out = out_dir if coexpr.EDGE_TYPE is None else os.path.join(out_dir, f"{et}_only")
        os.makedirs(out, exist_ok=True)
        logger.info("[iso2func-netviz] interaction=%s -> %s", et, out)
        cross_state_block(dataset_dir, reference_group, out, top_tfs, min_expr, min_edges, extra_tfs=extra_tfs)
        for (ga, gb) in (comparisons or []):
            comparison_block(dataset_dir, ga, gb, out, min_expr, min_edges)
    coexpr.EDGE_TYPE = None
    logger.info("[iso2func-netviz] DONE -> %s", out_dir)
    return out_dir
