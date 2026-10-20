"""Combined query network export using ISV's shared Fruchterman–Reingold layout.

Rendering is a presentation of explicit assay calls, not an inferred activity model.
Both sides of every qualifying difference survive: positive and tested negative.
"""

import json
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

from . import cytoscape


def query_graph(divergent, pairs):
    nodes, edges = {}, {}
    by_clone = {}
    for r in pairs:
        for side in ("A", "B"):
            by_clone[r["clone_" + side]] = dict(gene=r["gene"], aa=r["aa_" + side],
                native_aliases=r["native_aliases_" + side], reads=r["max_single_native_reads_" + side],
                native_evidence=r.get("native_evidence_"+side, "exact"))
    for r in divergent:
        target = r["modality"] + ":" + r["target"]
        nodes[target] = dict(id=target, label=r["target"], node_type="protein" if r["modality"] == "ppi" else "dna_target",
                             gene_symbol=r["target"] if r["modality"] == "ppi" else "", aa="", native_aliases="", reads="",
                             max_state_mean_cpm=r["max_state_mean_cpm"])
        for side in ("A", "B"):
            clone = r["clone_" + side]
            nodes[clone] = dict(id=clone, label=clone + ("*" if by_clone[clone]["native_evidence"] == "AA_substitutions" else ""),
                                node_type="isoform", gene_symbol=by_clone[clone]["gene"],
                                aa=by_clone[clone]["aa"], native_aliases=by_clone[clone]["native_aliases"],
                                reads=by_clone[clone]["reads"], native_evidence=by_clone[clone]["native_evidence"], max_state_mean_cpm=None)
            key = clone, target, r["modality"]
            rec = edges.setdefault(key, dict(source=clone, target=target, edge_type=r["modality"],
                assay="y2h" if r["modality"] == "ppi" else "ey1h", detected=r["call_" + side], tested=True,
                score="", comparison_pairs=[]))
            if rec["detected"] != r["call_" + side]:
                raise ValueError("Inconsistent combined-network assay call")
            pair = r["clone_A"] + " / " + r["clone_B"]
            if pair not in rec["comparison_pairs"]:
                rec["comparison_pairs"].append(pair)
    for r in edges.values():
        r["comparison_pairs"] = json.dumps(sorted(r["comparison_pairs"]))
    return pd.DataFrame(nodes.values()), pd.DataFrame(edges.values())


def _layout(nodes, edges):
    from altanalyze3.components.visualization.isv_web.network import _combined_layout
    # Reuse ISV's union layout, including its same-gene isoform glue. No gene-CPM
    # aggregation or NMD filtering from the separate ISV expression model is invoked.
    panel = dict(nodes=[dict(id=r.id, kind="isoform" if r.node_type == "isoform" else "gene",
                            gene=r.gene_symbol or r.id) for r in nodes.itertuples()],
                 edges=[dict(src_iso=r.source, tgt_gene=r.target) for r in edges.itertuples()])
    positions = _combined_layout(panel, dict(nodes=[], edges=[]))
    if not positions:
        raise ImportError("Combined network requires python-igraph for the existing ISV FR layout")
    return positions


def _legend(manifest):
    p = manifest["parameters"]
    expression = manifest["coverage"].get("expression", {})
    states = expression.get("states", p.get("state") or [])
    n = expression.get("donors", "not requested")
    threshold = p.get("min_partner_cpm")
    substitutions = p.get("max_aa_substitutions", 0)
    native_match = (f"same-gene, same-length protein match with <= {substitutions} AA substitutions; "
                    "asterisks mark proteins supported through substitutions rather than exact identity; "
                    "their assay results predict possible native function. " if substitutions else
                    "literal same-gene protein match; ")
    protein = (f"Protein partners: mean native gene CPM > {threshold:g} in at least one state; zeros retained. "
               if threshold is not None else "Protein partner expression was not thresholded. ")
    return ("Discovery screen of isoform-dependent measured interactions. Diamonds: assayed TF proteins (labels include AA length); "
            "circles: PPI partner genes; squares: DNA baits. Blue solid edges: assay positive; gray dashed edges: explicitly tested negative. "
            "PPI lines are undirected; PDI arrows denote binding to DNA, not activation or repression. "
            "Only targets with an explicit True/False difference in at least one retained pair are drawn; both calls are shown. "
            "Native support, when requested: " + native_match +
            f">= {p.get('min_isoform_reads', 0):g} assigned reads per individual isoform; aliases are not pooled. "
            + protein + f"States: {', '.join(states) or 'not restricted'}; conditions: {', '.join(p.get('condition') or []) or 'all'}. "
            f"Gene-expression biological units: {n} donor libraries per state. "
            "DNA-bait expression is undefined and unfiltered; no bait is relabeled as a gene. "
            "Uniform node/edge sizes have no quantitative activity meaning. No statistical test, confidence interval, or inferred activity score is used. "
            "Atlas assays were not performed in these HSPCs; RNA presence supports relevance, not in-cell interaction validation. "
            "Unknown and unchanged comparisons remain in all_pair_target_calls.tsv and all_pair_decisions.tsv; absent assay data never imply a loss.")


def _draw(ax, nodes, edges, positions, title, overview=False):
    from .ego_plot import _draw_edge
    ax.set_title(title, fontsize=11)
    if edges.empty:
        ax.text(.5, .5, "No qualifying differences", ha="center", transform=ax.transAxes)
        ax.axis("off")
        return
    ids = set(edges.source) | set(edges.target)
    nd = nodes[nodes.id.isin(ids)]
    xy = np.array([positions[n] for n in ids])
    span = max(float(np.ptp(xy[:, 0])), float(np.ptp(xy[:, 1])), 1)
    scale = 12 / span
    pos = {n: (positions[n][0] * scale, positions[n][1] * scale) for n in ids}
    for r in edges.itertuples():
        color = "#2878A8" if r.detected else "#A9ADB2"
        if r.detected:
            _draw_edge(ax, pos[r.source], pos[r.target], color, .75, directed=r.edge_type == "pdi",
                       alpha=.55, shrinkA=.08, shrinkB=.07)
        else:
            a, b = pos[r.source], pos[r.target]
            ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=.55, alpha=.6, ls="--", zorder=0)
            if r.edge_type == "pdi":
                # A negative PDI is an assay comparison, not a directed regulatory effect.
                ax.annotate("", xy=b, xytext=(a[0]+.92*(b[0]-a[0]), a[1]+.92*(b[1]-a[1])),
                            arrowprops=dict(arrowstyle="-|>", color=color, lw=.55))
    texts = []
    for r in nd.itertuples():
        marker, color = {"isoform": ("D", "#F2CA73"), "protein": ("o", "#D9E9F1"),
                         "dna_target": ("s", "#DED3E9")}[r.node_type]
        x, y = pos[r.id]
        ax.scatter([x], [y], marker=marker, c=[color], s=55 if r.node_type == "isoform" else 22,
                   edgecolors="#4E545A", linewidths=.5, zorder=3)
        if not overview or r.node_type == "isoform":
            label = f"{r.label} ({r.aa} aa)" if r.node_type == "isoform" else r.label
            texts.append(ax.text(x, y + .08, label, ha="center", fontsize=8,
                         bbox=dict(facecolor="white", alpha=.7, edgecolor="none", pad=.15), zorder=4))
    ax.margins(.10)
    ax.set_aspect("equal")
    ax.axis("off")
    if texts:
        import inspect
        from adjustText import adjust_text
        # adjustText 0.8 and 1.x use different iteration-limit names.
        options = {"lim": 100} if "lim" in inspect.signature(adjust_text).parameters else {"iter_lim": 100}
        adjust_text(texts, x=[pos[n][0] for n in ids], y=[pos[n][1] for n in ids], ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#777777", lw=.35), **options)


def export_query_network(divergent, pairs, out_dir, manifest):
    """Write graph tables/JSON, combined PDF and complete per-gene detail pages."""
    from matplotlib.backends.backend_pdf import PdfPages
    from .ego_plot import _mpl
    import adjustText  # fail before exports if label separation is unavailable
    out = Path(out_dir)
    if not divergent:
        return dict(status="no_qualifying_interactions")
    nodes, edges = query_graph(divergent, pairs)
    positions = _layout(nodes, edges)
    nodes["x"] = nodes.id.map(lambda n: positions[n][0])
    nodes["y"] = nodes.id.map(lambda n: positions[n][1])
    nodes.to_csv(out / "combined_network.nodes.tsv", sep="\t", index=False)
    edges.to_csv(out / "combined_network.edges.tsv", sep="\t", index=False)
    cytoscape.write_cytoscape_json(nodes, edges, str(out / "combined_network.cyjs.json"), include_undetected=True)
    plt = _mpl()
    pdf = out / "combined_divergent_interactions.pdf"
    legend = _legend(manifest)
    genes = sorted({r["gene"] for r in pairs})
    with PdfPages(pdf) as pages:
        for gene in [None] + genes:
            sub = edges if gene is None else edges[edges.source.isin(nodes.loc[nodes.gene_symbol.eq(gene) & nodes.node_type.eq("isoform"), "id"])]
            modalities = ["ppi", "pdi"] if gene is None else [m for m in ["ppi", "pdi"] if sub.edge_type.eq(m).any()]
            size = (28, 16) if gene is None else ((20, 13) if len(modalities) == 2 else (15, 12))
            fig, axes = plt.subplots(1, len(modalities), figsize=size, squeeze=False)
            fig.subplots_adjust(left=.025, right=.98, bottom=.17, top=.91, wspace=.08)
            title = "All supported TF pairs" if gene is None else gene
            detail_positions = positions if gene is None else _layout(nodes[nodes.id.isin(set(sub.source) | set(sub.target))], sub)
            for ax, modality in zip(axes.flat, modalities):
                es = sub[sub.edge_type.eq(modality)]
                _draw(ax, nodes, es, detail_positions, f"{modality.upper()}: {len(es)} explicit assay calls", overview=gene is None)
            fig.suptitle(f"{title}: isoform-dependent PPI and DNA-bait binding", fontsize=17)
            page_legend = legend + (" Overview labels identify TF proteins; all target labels appear on the following complete gene-detail pages." if gene is None else "")
            fig.text(.03, .025, textwrap.fill(page_legend, width=250 if gene is None else (180 if len(modalities) == 2 else 145)),
                     fontsize=9, va="bottom", linespacing=1.45)
            pages.savefig(fig)
            plt.close(fig)
    return dict(pdf=str(pdf), nodes=len(nodes), edges=len(edges), genes=len(genes),
                pages=1+len(genes), layout="ISV _combined_layout / igraph Fruchterman-Reingold",
                label_separation="adjustText " + getattr(adjustText, "__version__", "unknown"),
                edge_policy="Unique clone/modality/target calls; positive and negative comparisons both shown")
