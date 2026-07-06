# SNAF-B (surface / B-antigen) workflow

Predicts tumor-specific **cell-surface** antigens (altered extracellular topology) from splicing.
Reuses the SNAF-T GTEx tumor-specificity check, then recovers ORFs and predicts transmembrane
topology. Cross-platform pure-Python; **no TMHMM binary required**.

## CLI

```bash
altanalyze3 snaf-b \
  --juncounts counts.txt --db_dir <SNAF_reference> \
  --freq_path frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt \
  --mode short_read --validation_gtf ref.gtf --genome_fasta genome.fa \
  --min_samples 1 --cpus 8 --output out/
```

`--freq_path` is a SNAF-T frequency table (SNAF-B depends on the T pipeline for the GTEx check).
`--mode` = `short_read` or `long_read`.

## Pipeline stages

1. **get_membrane_tuples** — filter to membrane/surface splicing events that are tumor-specific.
2. **surface.run** — ORF recovery + transmembrane topology per candidate, at stringency 3/4/5.
3. **generate_full_results** — stringency × style × overlap result tables.

## Dependency modernization (no external binaries; graceful degradation)

| Legacy dependency | Replacement | Behavior if unavailable |
|---|---|---|
| TMHMM 2.0 binary | pure-Python `tmhmm.py` parsing the bundled `TMHMM2.0.model` | logs "no TMHMM library"; **continues** with remaining steps (surface results still produced at lower topology confidence) |
| EMBOSS `needle` | Biopython `PairwiseAligner` (BLOSUM62) | — |
| Ensembl REST / UCSC | offline resolver via `--genome_fasta` + local GTF | — |
| mygene (ENSG→symbol) | batched + cached; graceful fallback to ENSG | keeps ENSG if network down |

The legacy TMHMM/EMBOSS binaries are still used if a `software_path` is supplied.

## Speed optimizations (validated)

| Optimization | Speedup | Identity |
|---|---|---|
| `is_support` memoization by (uid, ORF, strict) + GTF-parse cache + region-table cache | >7× (>30 min timeout → 258 s) | stringency-3 candidate set identical |
| scratch-cleanup guards (only run when `software_path` given / file exists) | removes spurious errors | — |

All SNAF-B features are exposed via the CLI (`snaf-b`) **and** the nf-core workflow (`--mode surface`).

## Notes
- If a pure-Python library for a step is missing, SNAF-B prints that it lacks the library and proceeds with the remaining steps (never silently fails).
- Optimized copies of external tools live under `snaf/` (or in altanalyze3 if already present), per the project convention.
