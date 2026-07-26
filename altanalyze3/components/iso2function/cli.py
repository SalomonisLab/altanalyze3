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


def add_arguments(parser):
    """Attach iso2function arguments to an existing argparse parser/subparser. Returns the parser so a
    caller in ``utilities/parser.py`` can do ``add_arguments(sub.add_parser('sclr-iso2func'))``."""
    parser.add_argument("action",
                        choices=["ingest", "crosswalk", "associate", "uniprot", "network", "enrich",
                                 "export", "switch-plots", "stacked-bar", "all"],
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
                        help="switch-plots: TSV of isoform-switch pairs "
                             "(cols gene,ensg,isoA,isoB,A_clone,B_clone[,A_m1h_call,B_m1h_call])")
    parser.add_argument("--pseudobulk", default=None,
                        help="switch-plots: isoform pseudobulk h5ad "
                             "(obs '<cellstate>.<sample>-isoform', var '<ENSG>:<isoform>')")
    parser.add_argument("--donor", action="append", default=None,
                        help="switch-plots: restrict to these sample tokens (repeatable); default all")
    parser.add_argument("--lineage-map", default=None,
                        help="switch-plots: optional 2-col TSV (cell_state<TAB>lineage) to add lineage plots")
    parser.add_argument("--top-states", type=int, default=6,
                        help="switch-plots: number of largest-ratio-shift cell states in the network (default 6)")
    # ---- stacked-bar action (plot_stacked_isoform_bar) ----
    parser.add_argument("--gene", default=None, help="stacked-bar: gene symbol")
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
    parser.add_argument("--out", default=None, help="stacked-bar: output base path (.pdf/.png appended)")
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
