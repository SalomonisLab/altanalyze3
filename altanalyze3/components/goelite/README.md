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

### Supported Species

Current built-in species (extensible):
* `human` — downloads `goa_human.gaf.gz`
* `mouse` — downloads `goa_mouse.gaf.gz`

The URL targets point to the “current” release distributed by the GO
Consortium. Cached entries include the exact `data-version` (from the OBO) and
`!date` or `!Generated` metadata (from the GAF).

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
