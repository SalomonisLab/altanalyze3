# SNAF-B (surface / B-antigen) workflow

Predicts tumor-specific **cell-surface** antigens (altered extracellular topology) from splicing.
Reuses the SNAF-T GTEx tumor-specificity check, then recovers ORFs and predicts transmembrane
topology. Cross-platform pure-Python; **no TMHMM binary required**.

## CLI

```bash
altanalyze3 snaf-b \
  --juncounts counts.txt --db_dir <SNAF_reference> \
  --mode short_read --isoform_method learned --genome_fasta genome.fa \
  --min_samples 1 --cpus 8 --output out/
```

SNAF-B performs its own GTEx tumor-specificity check. `--freq_path` is optional; when omitted,
SNAF-B derives the frequency table, with no prior SNAF-T run or HLA types required.
`--mode` supports `short_read`, `long_read`, and the legacy `find_full_length` workflow.
Long-read validation is optional: add `--validation_gtf catalog.gtf` for stringency 4/5.

## Pipeline stages

1. **get_membrane_tuples** — filter to membrane/surface splicing events that are tumor-specific.
2. **surface.run** — reconstruct and rank synthetic isoforms, then compare proteins and predict
   transmembrane topology. In short-read CLI mode the learned predictor is the default.
3. **generate_full_results** — stringency × style × overlap result tables.

## Synthetic isoforms without long reads

Both revised methods are integrated into `altanalyze3 snaf-b --mode short_read`:

| `--isoform_method` | Behavior | Use |
|---|---|---|
| **`learned` (default)** | Bundled linear ranker; union of searches through four edits | Best mean protein similarity in the internal evaluation |
| `evidence` | Explicit sample-support rule; two-edit search | Best exact splice-chain recovery in the internal evaluation |
| `legacy` | Previous reference transcript editing and ORF selection | Reproduce earlier SNAF-B runs |

The complete `--juncounts` matrix supplies joint sample evidence, including background
junctions that did not pass tumor-specificity sifting. Reconstruction uses reference ENST
chains from Alt91; long-read validation catalogs never enter candidate construction.
The bundled ranker needs no external model file. Override it with `--isoform_ranker model.json`
or a fitted pickle bundle. `--isoform_sample SAMPLE` restricts combinations to that sample;
`--isoform_min_reads`, `--isoform_max_edits`, and `--isoform_top_k` control detection and search.
Defaults are three reads and five ranked alternatives passed to surface checks.

For selected first exons, supply `--first_exon_junctions first_exons.txt`, one junction UID
per line. If its donor has no annotated first-exon 5′ boundary, the exon is assumed to be
**250 nt long** in transcript orientation. Known boundaries take precedence. This flag
selects a first-exon hypothesis; an intronic splice site alone does not establish a new start.

`--genome_fasta` supplies genomic sequence, including novel first exons. If omitted, the
predictor uses Alt91's gene FASTA and its 2-kb flanks, with explicit coverage checks.
Candidates extending beyond available sequence remain unresolved; provide the genome
FASTA when the flanked gene sequence is insufficient.
Use references and counts from the same genome assembly.

Outputs under `synthetic_isoforms/` include ranked `predictions.tsv`, top-ranked GTF,
mRNA/protein FASTAs, and a configuration manifest. Predictions retain NMD annotations;
the existing stringency-3 surface gate excludes predicted NMD candidates. Stringency 4/5
add independent long-read support when a validation catalog is provided. Reports label
these candidates `synthetic_learned` or `synthetic_evidence`, with artifact IDs linking
to the detailed predictions. Synthetic GTF files are predictions, not validation catalogs.

See [Synthetic isoform inference](SYNTHETIC_ISOFORM_INFERENCE.md) for the algorithm,
evaluation, three exemplar structures, Python API, standalone CLI, and limitations.

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

The Nextflow workflow (`--mode surface`) also defaults to `--isoform_method learned` and
accepts `--isoform_method evidence`, `--isoform_ranker`, `--first_exon_junctions`, and
`--genome_fasta`. Additional SNAF-B flags can be passed through `SNAF_B`'s `task.ext.args`.
It publishes `synthetic_isoforms/` alongside the surface results.

## Notes
- If a pure-Python library for a step is missing, SNAF-B prints that it lacks the library and proceeds with the remaining steps (never silently fails).
- Optimized copies of external tools live under `snaf/` (or in altanalyze3 if already present), per the project convention.
- Evidence reconstruction currently runs serially; `--cpus` still controls existing parallel
  pipeline stages. Legacy reference-indexed visualization methods do not display synthetic
  candidate geometry; use the exported GTF and prediction table for structure review.
