#!/usr/bin/env python3
"""
CLI entry point for GO-Elite component.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import List

import pandas as pd

import logging

if __package__ in (None, ""):
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    from altanalyze3.components.goelite.parser import ParsedGO, build_tree_with_annotations
    from altanalyze3.components.goelite.resources import AVAILABLE_SPECIES, prepare_species_resources
    from altanalyze3.components.goelite.runner import EnrichmentSettings, GOEliteRunner
    from altanalyze3.components.goelite.plotting import write_goelite_scatter_pdf
else:
    from .parser import ParsedGO, build_tree_with_annotations
    from .resources import AVAILABLE_SPECIES, prepare_species_resources
    from .runner import EnrichmentSettings, GOEliteRunner
    from .plotting import write_goelite_scatter_pdf


def _read_gene_list(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def main(argv: List[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Run GO-Elite enrichment inside AltAnalyze3.")
    ap.add_argument("--obo", help="GO ontology OBO file (optionally gzipped).")
    ap.add_argument("--gaf", help="GO annotation GAF file.")
    ap.add_argument(
        "--species",
        choices=AVAILABLE_SPECIES,
        help="Preferred species to load from the GO-Elite cache (e.g. human, mouse).",
    )
    ap.add_argument("--cache-dir", help="Override the GO-Elite cache directory.")
    ap.add_argument("--version", help="Specific cached version to use/build (default: latest available).")
    ap.add_argument("--force-refresh", action="store_true", help="Force rebuilding cached resources for the species.")
    ap.add_argument("--query", required=True, help="Newline-delimited file of query genes.")
    ap.add_argument("--background", help="Optional newline-delimited file of background genes.")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--min-term-size", type=int, default=5)
    ap.add_argument("--max-term-size", type=int, default=2000)
    ap.add_argument("--min-z", type=float, default=1.96)
    ap.add_argument("--max-fdr", type=float, default=0.1)
    ap.add_argument("--min-overlap", type=int, default=2)
    ap.add_argument("--delta-z", type=float, default=0.5)
    ap.add_argument("--top-label-count", type=int, default=4, help="Number of top enriched terms to label in the PDF.")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("go_elite_cli")
    os.makedirs(args.outdir, exist_ok=True)

    parsed: ParsedGO
    if args.obo and args.gaf:
        logger.info("Loading GO ontology and annotations from explicit paths...")
        parsed = build_tree_with_annotations(args.obo, args.gaf)
    elif args.species:
        logger.info("Preparing GO-Elite cache for species '%s'...", args.species)
        parsed = prepare_species_resources(
            args.species,
            cache_dir=args.cache_dir,
            version=args.version,
            force=args.force_refresh,
        )
    else:
        ap.error("Either specify --obo and --gaf, or provide --species to use cached resources.")

    settings = EnrichmentSettings(
        min_term_size=args.min_term_size,
        max_term_size=args.max_term_size,
    )
    settings.prioritization.min_z = args.min_z
    settings.prioritization.max_fdr = args.max_fdr
    settings.prioritization.min_overlap = args.min_overlap
    settings.prioritization.delta_z = args.delta_z

    runner = GOEliteRunner(parsed, logger=logger, settings=settings)

    query_genes = _read_gene_list(args.query)
    if args.background:
        background_genes = _read_gene_list(args.background)
    else:
        background_genes = sorted({g for genes in parsed.term_to_genes.values() for g in genes})

    logger.info("Running GO-Elite enrichment for %d query genes.", len(query_genes))
    results = runner.run(query_genes, background_genes)

    if not results:
        logger.warning("No enriched GO terms detected.")
        return 0
    prepared_genes = runner.prepare_background(background_genes)

    query_map = {}
    for gene in query_genes:
        key = str(gene).upper()
        if key and key not in query_map:
            query_map[key] = str(gene)

    records = [
        {
            "term_id": res.term_id,
            "term_name": (runner.go_tree.get(res.term_id).name if runner.go_tree.get(res.term_id) is not None else ""),
            "z_score": res.z_score,
            "p_value": res.p_value,
            "fdr": res.fdr,
            "overlap": res.overlap,
            "term_total": res.total_genes,
            "background_total": res.background,
            "selected": res.selected,
            "blocked_by": res.blocked_by,
            "overlap_genes": ",".join(
                sorted(
                    query_map[key]
                    for key in prepared_genes.term_genes.get(res.term_id, set()) & set(query_map.keys())
                    if key in query_map
                )
            ),
        }
        for res in results
    ]
    df = pd.DataFrame.from_records(records)

    tsv_path = os.path.join(args.outdir, "goelite_results.tsv")
    df.to_csv(tsv_path, sep="\t", index=False)
    logger.info("Wrote results table: %s", tsv_path)

    pdf_path = os.path.join(args.outdir, "goelite_results.pdf")
    written_pdf = write_goelite_scatter_pdf(
        df,
        pdf_path,
        title_prefix="GO-Elite",
        top_label_count=args.top_label_count,
    )
    if written_pdf is not None:
        logger.info("Wrote results PDF: %s", written_pdf)

    json_path = os.path.join(args.outdir, "goelite_results.json")
    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(records, handle, indent=2)
    logger.info("Wrote results JSON: %s", json_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
