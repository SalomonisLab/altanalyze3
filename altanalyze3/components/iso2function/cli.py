"""Command-line entry point for iso2function.

Exposes the staged workflow as one subcommand with actions: ``ingest`` | ``crosswalk`` | ``network`` |
``enrich`` | ``all``. Designed to plug into ``altanalyze3/utilities/parser.py`` as the
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
                                 "export", "all"],
                        help="pipeline stage to run")
    parser.add_argument("--out-dir", default=None,
                        help="output dir for parsed/crosswalk tables (default: component data/ dir)")
    parser.add_argument("--artifacts-dir", default=None,
                        help="output dir for network/figure artifacts (default: component artifacts/)")
    parser.add_argument("--clone", action="append", default=None,
                        help="restrict the network to these paper2 clone_ids (repeatable)")
    parser.add_argument("--detected-only", action="store_true",
                        help="drop assayed-negative edges from the exported network")
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
