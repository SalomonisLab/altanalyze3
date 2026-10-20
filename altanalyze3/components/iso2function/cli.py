"""Command-line entry point for iso2function.

Exposes the staged workflow as one subcommand with actions: ``ingest`` | ``crosswalk`` | ``network`` |
``enrich`` | ``switch-plots`` | ``all``. Designed to plug into ``altanalyze3/utilities/parser.py`` as the
``sclr-iso2func`` subcommand: register :func:`add_arguments` on a subparser and set
``func=run_iso2function``. Heavy imports are deferred into the handlers so importing this module (at
parser-build time) stays cheap and cannot break the global CLI.

Standalone use:
    python -m altanalyze3.components.iso2function.cli all
    python -m altanalyze3.components.iso2function.cli ingest --out-dir <dir>
"""

import argparse
import logging
import sys

logger = logging.getLogger(__name__)


def add_arguments(parser, include_loglevel=True):
    """Attach iso2function arguments to an existing argparse parser/subparser. Returns the parser so a
    caller in ``utilities/parser.py`` can do ``add_arguments(sub.add_parser('sclr-iso2func'))``."""
    parser.add_argument("action",
                        choices=["ingest", "crosswalk", "associate", "uniprot", "network", "enrich",
                                 "export", "switch-plots", "stacked-bar", "query", "protein-impact", "protein-impact-plots", "protein-exon-map", "assay-match", "report", "all"],
                        help="pipeline stage to run")
    parser.add_argument("--out-dir", default=None,
                        help="output dir for parsed/crosswalk tables (default: component data/ dir)")
    parser.add_argument("--artifacts-dir", default=None,
                        help="output dir for network/figure artifacts (default: component artifacts/)")
    parser.add_argument("--clone", action="append", default=None,
                        help="restrict the network to these paper2 clone_ids (repeatable)")
    parser.add_argument("--detected-only", action="store_true",
                        help="drop assayed-negative edges from the exported network")
    # ---- switch-plots action ----
    parser.add_argument("--pairs", default=None,
                        help="switch-plots/query: TSV of isoform-switch pairs "
                             "(cols gene,ensg,isoA,isoB,A_clone,B_clone[,A_m1h_call,B_m1h_call])")
    parser.add_argument("--pseudobulk", default=None,
                        help="switch-plots/query: isoform pseudobulk h5ad (query requires counts) "
                             "(obs '<cellstate>.<sample>-isoform', var '<ENSG>:<isoform>')")
    parser.add_argument("--donor", action="append", default=None,
                        help="switch-plots/query: restrict samples (query: original_library; repeatable); default all")
    parser.add_argument("--lineage-map", default=None,
                        help="switch-plots: optional 2-col TSV (cell_state<TAB>lineage) to add lineage plots")
    parser.add_argument("--top-states", type=int, default=6,
                        help="switch-plots: number of largest-ratio-shift cell states in the network (default 6)")
    # ---- stacked-bar action (plot_stacked_isoform_bar) ----
    parser.add_argument("--gene", default=None, help="stacked-bar/query: gene symbol")
    parser.add_argument("--ensg", default=None, help="stacked-bar: ENSG id")
    parser.add_argument("--isoforms", default=None,
                        help="stacked-bar: '|'- or ','-separated ordered isoform ids to show")
    parser.add_argument("--gene-model", default=None, help="stacked-bar: Ensembl exon coords TSV (Hs_Ensembl_exon.txt)")
    parser.add_argument("--gff", default=None, help="stacked-bar: dataset gff-output dir")
    parser.add_argument("--level", default="cellstate", choices=["cellstate", "lineage"],
                        help="stacked-bar: per cell state or aggregate by lineage (default cellstate)")
    parser.add_argument("--plottype", default="bar", choices=["bar", "line"], help="stacked-bar: bar|line")
    parser.add_argument("--cohort-sample", action="append", default=None,
                        help="stacked-bar: keep only these sample tokens (repeatable); default all")
    parser.add_argument("--cluster-order", default=None,
                        help="stacked-bar: file with one cell-state/lineage per line (x-axis order)")
    parser.add_argument("--clone-map", default=None,
                        help="stacked-bar/switch-plots: 2-col TSV (isoform<TAB>clone_label) for legend labels")
    parser.add_argument("--out", default=None, help="stacked-bar: output base path; export --query-results: text file")
    parser.add_argument("--query-results", help="export/report: completed query directory")
    # Shared by the top-level and standalone CLI; no dataset-specific query wrappers.
    parser.add_argument("--interaction-type", choices=["all", "ppi", "pdi"], default="all")
    parser.add_argument("--metadata", help="query: library metadata TSV (library, groups, original_library)")
    parser.add_argument("--state", action="append", help="query: cell state to retain (repeatable)")
    parser.add_argument("--condition", action="append", help="query: metadata groups to retain (repeatable)")
    parser.add_argument("--protein-fasta", help="query: native protein FASTA, optionally gzip compressed")
    parser.add_argument("--min-isoform-reads", type=float, default=30,
                        help="query: minimum assigned reads for each native isoform (inclusive; default 30)")
    parser.add_argument("--max-aa-substitutions", type=int, default=0,
                        help="query: permit this many AA substitutions for same-gene, same-length proteins "
                             "without an exact assay match (default 0); report transferred assay evidence as putative")
    parser.add_argument("--gene-cpm", help="query: native gene-level CPM h5ad")
    parser.add_argument("--gene-counts", help="query: gene counts h5ad to verify native CPM against")
    parser.add_argument("--min-partner-cpm", type=float,
                        help="query: strict threshold on mean donor gene CPM in at least one selected state")
    parser.add_argument("--min-differential-partners", type=int, default=1,
                        help="query: minimum qualifying explicit True/False differences (default 1)")
    parser.add_argument("--clonelist", help="query: assayed clone catalogue TSV (default: configured CloneList)")
    parser.add_argument("--gene-annotation", help="query: gene ID/symbol TSV (default: configured Ensembl91)")
    parser.add_argument("--mane-summary", help="query: MANE summary TSV, optionally gzip compressed")
    parser.add_argument("--symbol-map", help="query: additive identifier TSV (symbol, ensg, source)")
    parser.add_argument("--combined-network", action="store_true",
                        help="query: export combined divergent-interaction network using ISV layout")
    parser.add_argument("--native-models", help="protein-impact: native model TSV with protein/CDS/cDNA sequences")
    parser.add_argument("--required-roster", help="protein-impact: independent ordered roster JSON with source hash")
    parser.add_argument("--reference-bundle", help="protein-impact: pinned AltAnalyze protein reference bundle JSON")
    parser.add_argument("--coding-regions", help="protein-impact: native coding_regions.txt (optionally gzipped)")
    parser.add_argument("--prior-feature-audit", help="protein-impact: earlier reference_mapped_features.tsv for regression")
    parser.add_argument("--impact-results", help="report/protein-impact-plots: completed protein-impact result directory")
    parser.add_argument("--genes", help="protein-impact-plots: comma-separated genes; default all genes with audited pairs")
    parser.add_argument("--protein-plot-style", choices=["compact", "detailed"], default="compact",
                        help="protein-impact-plots: compact domain architecture (default), or detailed audit figure")
    parser.add_argument("--protein-figure-layout", choices=["pages", "stacked"], default="pages",
                        help="protein-impact-plots: separate PDF pages (default), or stack all compact panels on one vector PDF page")
    parser.add_argument("--protein-figure-width", type=float, default=12,
                        help="protein-impact-plots: compact figure width in inches (default 12); increase for narrow exon labels")
    parser.add_argument("--isv-structure-index", help="protein-exon-map: existing registered ISV final_structures.db")
    parser.add_argument("--gene-model-index", help="protein-exon-map: existing ISV gene_model.db for AltAnalyze exon IDs")
    parser.add_argument("--coding-exon-map", help="completed protein-exon-map directory; plots read it, protein-exon-map reuses it without rescanning GFF/CDS")
    parser.add_argument("--feature-type", action="append", default=None,
                        help="protein-impact-plots: annotation type to display (repeatable); default functional domains, motifs and regions")
    parser.add_argument("--annotation-source", action="append", default=None,
                        help="protein-impact-plots: source provider to display (repeatable); default reviewed UniProt and literature annotations")
    parser.add_argument("--feature-reference", action="append", default=None,
                        help="protein-impact-plots: display a specified gene's reference protein as GENE=REFERENCE_PROTEIN (repeatable); all other source rows stay in the audit")
    parser.add_argument("--assay-results", help="report: completed assay-match result directory")
    if include_loglevel:
        parser.add_argument("--loglevel", default="INFO",
                            help="logging level (DEBUG/INFO/WARNING/ERROR)")
    return parser


def run_iso2function(args):
    """Dispatch the requested action. Imports each layer lazily."""
    logging.basicConfig(level=getattr(logging, str(getattr(args, "loglevel", "info")).upper(),
                                       logging.INFO), format="%(message)s")
    from . import config
    config.ensure_dirs()
    action = args.action

    if action == "protein-impact":
        from .protein_impact import run_protein_impact
        return run_protein_impact(args)
    if action == "protein-impact-plots":
        from .plotting.protein_features import run_protein_impact_plots
        return run_protein_impact_plots(args)
    if action == "protein-exon-map":
        from .coding_exons import run_protein_exon_map
        return run_protein_exon_map(args)
    if action == "assay-match":
        from .assay_match import run_assay_match
        return run_assay_match(args)
    if action == "report":
        from .evidence_report import run_report
        return run_report(args)
    if action == "query":
        from .query import run_query
        return run_query(args)
    if action == "export" and getattr(args, "query_results", None):
        from .export_interactions import export_query_results
        return export_query_results(args.query_results, args.out, args.interaction_type)

    if action in ("ingest", "all"):
        from .ingest import paper2, paper1
        manifest = paper2.parse_all(out_dir=args.out_dir)
        bad = manifest[~manifest["qc_ok"]]
        if len(bad):
            logger.error("[iso2function] paper2 ingest QC failures:\n%s", bad.to_string(index=False))
        p1 = paper1.parse_all(out_dir=args.out_dir)
        bad1 = p1[~p1["qc_ok"]]
        if len(bad1):
            logger.error("[iso2function] paper1 ingest QC failures:\n%s", bad1.to_string(index=False))

    if action in ("crosswalk", "all"):
        from .crosswalk import clone_map, paper1_map
        clone_map.build_crosswalk(data_dir=args.out_dir)
        try:
            paper1_map.build_paper1_crosswalk(data_dir=args.out_dir)
        except FileNotFoundError as e:
            logger.warning("[iso2function] skipping paper1 crosswalk: %s", e)

    if action in ("associate", "all"):
        from . import associate
        associate.build_all(data_dir=args.out_dir)

    if action in ("uniprot", "all"):
        from . import uniprot
        try:
            uniprot.build_all(data_dir=args.out_dir)
        except FileNotFoundError as e:
            logger.warning("[iso2function] skipping UniProt layer: %s", e)

    if action in ("network", "all"):
        from .network import build, cytoscape
        nodes, edges = build.build_graph(clone_ids=args.clone, data_dir=args.out_dir,
                                         detected_only=args.detected_only)
        build.write_graph(nodes, edges, out_dir=args.artifacts_dir)
        cytoscape.write_cytoscape_json(
            nodes, edges,
            include_undetected=not args.detected_only,
        )
        cytoscape.write_sif(edges)

    if action in ("switch-plots",):
        from .network import switch_pair_plots
        import pandas as pd
        if not args.pairs or not args.pseudobulk:
            logger.error("[iso2function] switch-plots requires --pairs <tsv> and --pseudobulk <h5ad>")
            return 2
        lineage_map = None
        if args.lineage_map:
            lm = pd.read_csv(args.lineage_map, sep="\t", header=None, dtype=str)
            lineage_map = dict(zip(lm.iloc[:, 0], lm.iloc[:, 1]))
        out = args.artifacts_dir or "switch_pair_plots"
        switch_pair_plots.render_switch_pairs(
            args.pairs, args.pseudobulk, out, data_dir=args.out_dir,
            donors=args.donor, lineage_map=lineage_map, top_states=args.top_states)

    if action in ("stacked-bar",):
        from .plotting import plot_stacked_isoform_bar
        import pandas as pd
        need = [args.gene, args.ensg, args.isoforms, args.pseudobulk, args.gene_model, args.gff, args.out]
        if not all(need):
            logger.error("[iso2function] stacked-bar requires --gene --ensg --isoforms --pseudobulk "
                         "--gene-model --gff --out")
            return 2
        isoforms = [x for x in args.isoforms.replace(",", "|").split("|") if x]
        lineage_map = None
        if args.lineage_map:
            lm = pd.read_csv(args.lineage_map, sep="\t", header=None, dtype=str)
            lineage_map = dict(zip(lm.iloc[:, 0], lm.iloc[:, 1]))
        clone_map = None
        if args.clone_map:
            cm = pd.read_csv(args.clone_map, sep="\t", dtype=str)
            c0, c1 = (("isoform", "clone") if {"isoform", "clone"}.issubset(cm.columns)
                      else (cm.columns[0], cm.columns[1]))
            clone_map = dict(zip(cm[c0], cm[c1]))
        cluster_order = None
        if args.cluster_order:
            cluster_order = [ln.strip() for ln in open(args.cluster_order) if ln.strip()]
        plot_stacked_isoform_bar(
            args.gene, args.ensg, isoforms, args.pseudobulk, args.out,
            gene_model=args.gene_model, gff_dir=args.gff, level=args.level, plottype=args.plottype,
            cohort_samples=args.cohort_sample, lineage_map=lineage_map, cluster_order=cluster_order,
            clone_map=clone_map)

    if action in ("enrich", "all"):
        from .enrichment import switch_enrichment
        switch_enrichment.build_all(data_dir=args.out_dir)

    if action in ("export", "all"):
        from . import export_interactions
        export_interactions.build_interactions_txt(data_dir=args.out_dir)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(prog="iso2function", description=__doc__.splitlines()[0])
    add_arguments(parser)
    args = parser.parse_args(argv)
    return run_iso2function(args)


if __name__ == "__main__":
    sys.exit(main())
