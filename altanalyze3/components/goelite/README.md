# GO-Elite Component (AltAnalyze3)

This document summarises the modernised GO-Elite implementation that now lives
under `altanalyze3/components/goelite/`. It details the main modules, public
interfaces, and how ontology / gene relationships are downloaded and cached for
reuse.

---

## Overview

The component exposes three primary layers:

| Layer | Location | Responsibility |
| --- | --- | --- |
| Parsing & Structures | `parser.py`, `structures.py` | Read GO OBO / GOA GAF files and represent the DAG + term -> gene relationships in memory. |
| Enrichment & Prioritisation | `runner.py`, `prio.py` | Compute per-term enrichment statistics (hypergeometric + Z-score) and prune redundant terms using the GO DAG. |
| Resource Management | `resources.py` | Automate download, caching, and versioning of ontology and annotation resources for supported species. |

`goelite/main.py` wires everything together into a CLI entry point, while
`__init__.py` re-exports the public API for other components.

---

## Key Data Structures

Defined in `structures.py`:

* `GOTermNode` — immutable representation of a GO term, including parents, children, depth, and the set of mapped genes.
* `GOTree` — container for all GO nodes with convenience methods (`ancestors`, `descendants`, `ensure_depths`, etc).
* `EnrichmentResult` — per-term statistics (z-score, overlap, adjusted p-value, selection flags).
* `compute_z_score` — helper for normal-approximation fold enrichment.

---

## Parsing Utilities

Defined in `parser.py`:

* `parse_obo(path)` — minimal OBO parser that builds the DAG as a `GOTree`.
* `parse_gaf(path, evidence=None)` — read GOA GAF gene-term annotations, filtering by evidence codes if provided.
* `build_tree_with_annotations(obo_path, gaf_path, acceptable_evidence=None)` — convenience wrapper returning a `ParsedGO` dataclass containing `(GOTree, term_to_genes)`. Gene sets are attached to the tree nodes.

Use these when you already have ontology files locally and do not need caching.

---

## Enrichment & Prioritisation

* `GOEliteRunner` (`runner.py`) — orchestrates the enrichment workflow.
  * `run(query_genes, background_genes, apply_prioritization=True)` — returns a list of `EnrichmentResult` objects. z-scores, raw p-values, and BH-FDR are computed automatically.
  * `EnrichmentSettings` — configure minimum/maximum term size and prioritisation thresholds.
* `prioritize_terms(tree, results, settings=None)` (`prio.py`) — implements GO-Elite’s DAG pruning logic given pre-computed enrichment results.
  * `PrioritizationSettings` allows tuning minimum z-score, FDR, overlap, and the z-score delta required for a child term to survive when a parent is selected.

---

## Automated Resource Management

The `resources.py` module adds a caching layer around GO downloads so that
pipelines no longer need to ship ontology files manually.

### Annotation protocol (default for every species, set 2026-10-24, revised the same day)

A species cache is built from the UNION of three annotation sources:

* the UniProt-GOA GAF (`goa_<species>.gaf.gz`);
* NCBI `gene2go.gz` rows of the species taxon, with GeneIDs mapped to symbols through the species
  `gene_info.gz` (gene2go covers every taxon, about 1.4 GB; it is downloaded once to
  `<cache root>/downloads/gene2go.gz` and reused);
* the same GAF + gene2go union of the ortholog species (`ortholog_species` in `SPECIES_CONFIG`: human
  gets mouse, mouse gets human), carried over through NCBI `gene_orthologs.gz` (about 130 MB, downloaded
  once to `<cache root>/downloads/gene_orthologs.gz`). ToppFun does this: IKZF1, IGHM and LYST reach
  "regulation of leukocyte activation" only through their mouse orthologs (2026-10-24 comparison).

NOT-qualified rows are excluded from every source, and every gene is propagated to the parent terms of
its annotations over is_a and part_of (the GO true-path rule) AND over regulates, positively_regulates
and negatively_regulates (as ToppFun does; "leukocyte activation" holds 512 human genes with is_a/part_of
alone and 1,157 with the full protocol). When a source cannot be obtained (download failure, missing
configuration) the build proceeds without it and records it: `metadata.json` carries
`annotation_protocol` (`union(goa_gaf,ncbi_gene2go)+ortholog(mouse:union(goa_gaf,ncbi_gene2go))`,
`union(...)`, or `single(<source>)`), `propagation` (`is_a,part_of,regulates` or `is_a,part_of`),
`ortholog_transfer` (pairs, transferred rows and genes), the per-source and union row counts
(`annotation_rows`) and `build_problems`. The version name ends in `_union` or the single source used,
then `_orth` when orthologs were transferred and `_reg` when regulates propagation was used.
`data/synonyms.tsv.gz` (GAF synonyms plus gene_info synonyms) lets the CLI map outdated query symbols to
current ones; every mapping is written to `<name>_gene_mapping.tsv`.
A cache built before this protocol is still readable: propagation is applied at load with the relations
recorded at build time, and the richest cache (ortholog transfer, then union, then single source) is
preferred across every version and cache root. A rebuild never removes an existing cache version; it is
moved under `<species>/delete/`.

Current built-in species (extensible through `SPECIES_CONFIG`: GAF URL, taxon id, gene_info URL,
ortholog species):
* `human` — `goa_human.gaf.gz`, taxon 9606, `Homo_sapiens.gene_info.gz`, orthologs from `mouse`
* `mouse` — `goa_mouse.gaf.gz`, taxon 10090, `Mus_musculus.gene_info.gz`, orthologs from `human`

The URL targets point to the current releases of the GO Consortium and NCBI. Cached entries
record the OBO `data-version` and the GAF `!date-generated` (or `!date` / `!Generated`) release.
Rebuild with `--force-refresh` (optionally `--obo`, `--gaf`, `--gene-info-path`, `--gene2go-path`,
`--gene-orthologs-path`, `--ortholog-gaf-path`, `--ortholog-gene-info-path` to reuse downloaded files;
`--no-gene2go` for a GAF-only build, `--no-ortholog-transfer` and `--no-regulates-propagation` to switch
the two ToppFun-matching steps off).

### Denominator

The enrichment denominator is every annotated gene of the resource. `--background` is refused
unless `--allow-custom-background` is also given (rule of record, 2026-10-24).

### Cache Location

By default resources are cached at `~/.altanalyze3/goelite`. The location can be
customised by:

* Passing `cache_dir` into the resource helpers.
* Setting environment variable `ALTA_GOELITE_DATA`.

### API Surface

* `resolve_cache_dir(cache_dir=None)` — resolve/create the root cache directory.
* `prepare_species_resources(species, cache_dir=None, version=None, obo_path=None, gaf_path=None, force=False)` — ensure the cache has the requested species/version. Downloads and builds the DAG + annotations if missing (or `force=True`).
* `load_cached_resources(species, cache_dir=None, version=None)` — load cached `ParsedGO` without attempting to refresh.
* `AVAILABLE_SPECIES` — tuple of species keys recognised by the downloader.
* `DEFAULT_CACHE_ENV` — name of the environment variable controlling the cache root.

### Storage Layout

```
<cache_root>/
  human/
    latest.json
    versions.json
    <version>/
      metadata.json
      downloads/
        go-basic.obo.gz
        goa_human.gaf.gz
      data/
        go_tree.json
        term_gene.parquet (or term_gene.tsv.gz fallback)
```

`metadata.json` captures source URLs, release identifiers, and counts; this is
useful for auditing or rebuilding downstream artefacts.

---

## CLI Usage

`goelite/main.py` exposes a command-line interface. Example invocations:

```bash
# Run with cached human resources (auto-download if not present)
python3 -m altanalyze3.components.goelite.main \
  --species human \
  --query genes_of_interest.txt \
  --background background_genes.txt \
  --outdir results/goelite-human

# Run using explicit ontology files (no caching)
python3 -m altanalyze3.components.goelite.main \
  --obo path/to/go-basic.obo \
  --gaf path/to/goa_human.gaf \
  --query genes_of_interest.txt \
  --background background_genes.txt \
  --outdir results/goelite-custom
```

Additional options:

| Flag | Description |
| --- | --- |
| `--cache-dir` | Override cache root (default `~/.altanalyze3/goelite` or env). |
| `--version` | Select or build a specific cached version. |
| `--force-refresh` | Rebuild the cached resources even if already present. |
| `--min-term-size`, `--max-term-size` | Filter GO terms by size. |
| `--min-z`, `--max-fdr`, `--min-overlap`, `--delta-z` | Configure prioritisation thresholds. |

---

## Programmatic Usage Snippets

```python
from altanalyze3.components.goelite import (
    GOEliteRunner,
    EnrichmentSettings,
    prepare_species_resources,
)

# Ensure cached resources exist (downloads if missing)
parsed = prepare_species_resources("human")

runner = GOEliteRunner(parsed, settings=EnrichmentSettings(min_term_size=10))

results = runner.run(
    query_genes={"STAT1", "IFI16", "IRF7"},
    background_genes={"STAT1", "IFI16", "IRF7", "ACTB", "GAPDH", "MYC"},
)

selected_terms = [res for res in results if res.selected]
```

When you need full control, bypass the cache:

```python
from altanalyze3.components.goelite import build_tree_with_annotations, GOEliteRunner

parsed = build_tree_with_annotations("go-basic.obo", "goa_human.gaf")
runner = GOEliteRunner(parsed)
```

---

## Testing

Unit tests live in `altanalyze3/components/tests/test_goelite_basic.py`:

* Validates prioritisation blocking semantics.
* Checks enrichment scoring.
* Exercises the cache round-trip workflow using synthetic OBO/GAF content.

Run with:

```bash
python3 -m pytest altanalyze3/components/tests/test_goelite_basic.py
```

---

## Extending / Customising

* To add a new species, extend `SPECIES_CONFIG` in `resources.py` with the GAF
  URL and filename. All caching logic reuses this mapping automatically.
* For other ontologies or custom gene sets, reuse `build_tree_with_annotations`
  to create `ParsedGO` instances and feed them into the existing runner.
* The cache layout is versioned; you can safely maintain multiple historical
  builds and switch between them via the `--version` flag or the `version`
  parameter to `prepare_species_resources`.

---

Feel free to adapt this module-specific documentation into the broader project
docs or developer guides as needed. Contributions and improvements are welcome!
