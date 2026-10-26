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
    from altanalyze3.components.goelite.parser import build_flat_gene_sets
    from altanalyze3.components.goelite.runner import EnrichmentSettings, GOEliteRunner
    from altanalyze3.components.goelite.plotting import write_goelite_scatter_pdf
else:
    from .parser import ParsedGO, build_tree_with_annotations
    from .resources import AVAILABLE_SPECIES, prepare_species_resources
    from .parser import build_flat_gene_sets
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
    ap.add_argument("--gene-sets", help="Flat gene-set collection instead of GO: an AltAnalyze gene-mapp file "
                                        "(UID, SystemCode, MAPP; e.g. EnsMart100 Ensembl-MAPP.txt = WikiPathways) "
                                        "or a gene/set[/name] table. Query and background must use the file's gene ids.")
    ap.add_argument("--gene-sets-label", default="WikiPathways", help="Namespace label written for --gene-sets terms.")
    ap.add_argument("--cache-dir", help="Override the GO-Elite cache directory.")
    ap.add_argument("--version", help="Specific cached version to use/build (default: latest available).")
    ap.add_argument("--force-refresh", action="store_true", help="Force rebuilding cached resources for the species.")
    ap.add_argument("--gene2go-path", help="Local NCBI gene2go.gz for a rebuild (default: download once to the cache root).")
    ap.add_argument("--gene-info-path", help="Local NCBI gene_info.gz for a rebuild (default: download the species file).")
    ap.add_argument("--no-gene2go", action="store_true",
                    help="Rebuild from the GOA GAF alone. The default protocol is the UNION of GOA and NCBI gene2go.")
    ap.add_argument("--no-ortholog-transfer", action="store_true",
                    help="Rebuild without the ortholog species' annotations (default: human gets mouse, mouse gets "
                         "human, through NCBI gene_orthologs; ToppFun does the same).")
    ap.add_argument("--no-regulates-propagation", action="store_true",
                    help="Rebuild propagating over is_a/part_of only (default also propagates over regulates, as ToppFun).")
    ap.add_argument("--gene-orthologs-path", help="Local NCBI gene_orthologs.gz for a rebuild (default: download once to the cache root).")
    ap.add_argument("--ortholog-gaf-path", help="Local GOA GAF of the ortholog species for a rebuild.")
    ap.add_argument("--ortholog-gene-info-path", help="Local NCBI gene_info.gz of the ortholog species for a rebuild.")
    ap.add_argument("--query", required=True, help="Newline-delimited file of query genes.")
    # Nathan Salomonis, 2026-10-24: "The denominator should always be all genes ... ensure the default is always
    # all genes." The default background is every annotated gene of the resource. A custom background changes the
    # question and is refused unless --allow-custom-background states that this is intended.
    ap.add_argument("--background", help="Custom newline-delimited background genes. REFUSED unless --allow-custom-background is "
                                         "also given: the standing rule is that the denominator is every annotated gene.")
    ap.add_argument("--allow-custom-background", action="store_true",
                    help="Confirms that a custom --background is intended (standing rule: all annotated genes).")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--min-term-size", type=int, default=5)
    ap.add_argument("--max-term-size", type=int, default=2000)
    ap.add_argument("--min-z", type=float, default=1.96)
    ap.add_argument("--max-fdr", type=float, default=0.1)
    ap.add_argument("--max-rawp", type=float, default=None,
                    help="Gate term significance on raw p-value (not FDR); also plots raw p on the y-axis.")
    ap.add_argument("--min-overlap", type=int, default=2)
    ap.add_argument("--delta-z", type=float, default=0.5)
    ap.add_argument("--top-label-count", type=int, default=None,
                    help="Number of PDF labels: the selected terms first, then the most significant pruned terms (default 40; "
                         "the table holds every term).")
    ap.add_argument("--name", type=str, default=None,
                    help="Comparison name; prepended to output filenames and used as the plot title.")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("go_elite_cli")
    os.makedirs(args.outdir, exist_ok=True)

    parsed: ParsedGO
    if args.gene_sets:
        logger.info("Loading flat gene sets (%s) from %s", args.gene_sets_label, args.gene_sets)
        parsed = build_flat_gene_sets(args.gene_sets, label=args.gene_sets_label)
    elif args.obo and args.gaf:
        logger.info("Loading GO ontology and annotations from explicit paths...")
        parsed = build_tree_with_annotations(args.obo, args.gaf)
    elif args.species:
        logger.info("Preparing GO-Elite cache for species '%s'...", args.species)
        parsed = prepare_species_resources(
            args.species,
            cache_dir=args.cache_dir,
            version=args.version,
            force=args.force_refresh,
            gene2go_path=args.gene2go_path,
            gene_info_path=args.gene_info_path,
            use_gene2go=not args.no_gene2go,
            use_orthologs=not args.no_ortholog_transfer,
            use_regulates=not args.no_regulates_propagation,
            gene_orthologs_path=args.gene_orthologs_path,
            ortholog_gaf_path=args.ortholog_gaf_path,
            ortholog_gene_info_path=args.ortholog_gene_info_path,
        )
    else:
        ap.error("Either specify --obo and --gaf, --gene-sets, or provide --species to use cached resources.")

    settings = EnrichmentSettings(
        min_term_size=args.min_term_size,
        max_term_size=args.max_term_size,
    )
    settings.prioritization.min_z = args.min_z
    settings.prioritization.max_fdr = args.max_fdr
    settings.prioritization.max_rawp = args.max_rawp
    settings.prioritization.min_overlap = args.min_overlap
    settings.prioritization.delta_z = args.delta_z

    runner = GOEliteRunner(parsed, logger=logger, settings=settings)

    query_genes = _read_gene_list(args.query)
    # Map outdated symbols to the annotation's current symbols (GAF synonyms); every decision is written.
    annotated = {g.upper() for genes in parsed.term_to_genes.values() for g in genes}
    synonyms = {}
    if args.species and getattr(parsed, "version_dir", None):
        from altanalyze3.components.goelite.resources import load_symbol_synonyms
        synonyms = load_symbol_synonyms(parsed.version_dir)
    mapping_rows = []
    mapped_query = []
    for gene in query_genes:
        key = str(gene).upper()
        if key in annotated:
            mapped_query.append(gene); mapping_rows.append(dict(input=gene, used=gene, status="direct"))
        elif key in synonyms:
            mapped_query.append(synonyms[key]); mapping_rows.append(dict(input=gene, used=synonyms[key], status="synonym"))
        else:
            mapped_query.append(gene); mapping_rows.append(dict(input=gene, used=gene, status="unannotated"))
    n_syn = sum(r["status"] == "synonym" for r in mapping_rows)
    n_un = sum(r["status"] == "unannotated" for r in mapping_rows)
    logger.info("Query genes: %d; annotated directly %d; mapped through a synonym %d; unannotated %d.",
                len(query_genes), len(query_genes) - n_syn - n_un, n_syn, n_un)
    for r in mapping_rows:
        if r["status"] == "synonym":
            logger.info("  symbol mapped: %s -> %s", r["input"], r["used"])
    query_genes = mapped_query
    if args.background and not args.allow_custom_background:
        ap.error("--background is refused: the GO-Elite denominator is every annotated gene (Nathan, 2026-10-24). "
                 "Omit --background, or add --allow-custom-background if a custom denominator is really intended.")
    if args.background:
        logger.warning("Custom background in use (%s); the standing rule is all annotated genes.", args.background)
        background_genes = _read_gene_list(args.background)
    else:
        background_genes = sorted({g for genes in parsed.term_to_genes.values() for g in genes})
        logger.info("Background: all %d annotated genes of the resource (standing rule).", len(background_genes))

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

    prefix_map = f"{args.name}_" if args.name else ""
    pd.DataFrame(mapping_rows).to_csv(os.path.join(args.outdir, f"{prefix_map}gene_mapping.tsv"), sep="\t", index=False)
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

    # optional comparison name -> prepended to output filenames and used as the plot title
    prefix = f"{args.name}_" if args.name else ""
    title_prefix = args.name if args.name else "GO-Elite"

    tsv_path = os.path.join(args.outdir, f"{prefix}goelite_results.tsv")
    df.to_csv(tsv_path, sep="\t", index=False)
    logger.info("Wrote results table: %s", tsv_path)

    pdf_path = os.path.join(args.outdir, f"{prefix}goelite_results.pdf")
    written_pdf = write_goelite_scatter_pdf(
        df,
        pdf_path,
        title_prefix=title_prefix,
        top_label_count=args.top_label_count,
        significance=("rawp" if args.max_rawp is not None else "fdr"),
        max_sig=(args.max_rawp if args.max_rawp is not None else args.max_fdr),
        min_z=args.min_z,
    )
    if written_pdf is not None:
        logger.info("Wrote results PDF: %s", written_pdf)

    json_path = os.path.join(args.outdir, f"{prefix}goelite_results.json")
    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(records, handle, indent=2)
    logger.info("Wrote results JSON: %s", json_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
