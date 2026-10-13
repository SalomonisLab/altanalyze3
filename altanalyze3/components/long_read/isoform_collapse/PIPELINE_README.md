# Long-read Isoform Collapse Pipeline

The benchmarks and scored-collapse description below document an earlier iteration.
The current default is `scored_collapse.collapse_gene_jc`, called by `pipeline.py`,
with junction-core WTA, reference preference and native known/novel support floors.
See the collapse-only option and its validation status at the end of this document.

Scored, expression-and-structure-aware collapse of single-cell long-read isoforms into a
standardized cross-dataset catalog with per-sample re-keyed h5ads. Designed for millions of
reads/sample × up to ~200 samples. All numbers below are **measured** on the iPSC Ctrl + Aza
KINNEX samples (2 samples, ~16M pooled reads), not estimated.

## Inputs
Per sample, produced upstream by `gff_collapse` / `bam/isoform_structure_extract`:
- `gff-output/transcript_associations.txt` — TSV: `gene, strand, structure, molecule_id, source`.
  `structure` = AltAnalyze exon string (`E1.1|E2.1|...`), terminal 5'/3' coords already stripped.
- `<sample>.h5ad` — raw per-read matrix: `cells × molecules`, `var_names = "GENE:molecule_id"`, X int32.

## The algorithm (per gene)
Confirmed on the ITGA2B exemplar. Implemented in `scored_collapse.collapse_gene`.

1. **Exon blocks** for length: `E1.1|E1.2 -> E1` (introns excluded). Used ONLY to size bins /
   score. Collapse identity always uses the FULL structure string.
2. **Long bin** = structures with `exon_blocks > 0.8 * max_blocks`. Every long-bin structure is a
   DISTINCT representative isoform (full-length forms never collapse into each other).
3. **Score** `= (blocks / max_blocks) * reads`. A non-long structure collapses into the
   HIGHEST-SCORING long isoform it is a perfect CONTIGUOUS SUBSTRING of — each child to exactly
   ONE parent.
4. **Other bin** = structures not consumed by the long bin. Same scoring (no 0.8 gate): a structure
   is a representative unless it is a perfect substring of a higher-scoring other-bin structure.
5. A structure that is a substring of nothing remains its own isoform.

Hard invariants (enforced in `_check_invariants`): one-parent per structure; perfect-substring-only
collapse edges; representatives never collapse into each other; partition (every structure is a rep
xor a child).

Why the score: a short-but-abundant structure (e.g. 12-block / 5000 reads) must NOT become a sink —
the 0.8 length gate keeps it out of the long bin so it can't swallow full-length fragments. Among
the longest, expression breaks ties. This is the "expression AND structure" requirement.

No Ensembl needed at this stage; representative selection is purely score-based.

## Stages
- **Stage 1 — integrated collapse, chromosome-parallel** (`run_pipeline.stage1`):
  pool reads across ALL samples per gene (identical structures unify), collapse per gene, partition
  genes by chromosome, run in parallel. Bit-identical to single-process (genes are chromosome-local;
  verified 0/68,008 reference genes span >1 chromosome). Keeps ALL structures (no per-sample filter).
- **Stage 2 — cross-sample outlier removal** (`run_pipeline.stage2`):
  drop final isoforms with total reads across all samples `< min_total` (default 2 — "observed only
  once across all samples is ignored"). Writes `FINAL_isoform_catalog.tsv` +
  `FINAL_structure_to_exemplar.tsv`.
- **Stage 3 — per-sample h5ad re-key** (`_stage3.py`):
  build `var(GENE:mol) -> final_exemplar` (only surviving molecules), sparse grouping matrix `G`,
  `X_final = X_raw @ G`. Output `<sample>-final_isoform.h5ad` keyed on standardized isoform ids.

## Final isoform naming
`final_isoform_id` = the exemplar molecule id (sample-tagged first molecule of the representative
structure, e.g. `Ctrl:6020262`). The SAME id is used across all samples for the same structure, so
per-sample h5ads are directly comparable for differential abundance.

## Measured benchmarks (2 samples, 16.1M pooled reads)
| stage | result | time | peak RSS |
|---|---|---|---|
| 1 (parallel collapse, 8 workers) | 16,908 genes -> 709,103 isoforms | 19.8s | 171 MB / worker |
| 2 (outlier removal, >=2 reads) | 709,103 -> 233,591 catalog | 0.1s | — |
| 3 Ctrl (sparse re-key) | 8.9M molecules -> 212,779 finals | 23.0s | 3,523 MB |
| 3 Aza | 7.1M molecules -> 199,496 finals | 19.3s | 4,691 MB |

Read conservation: 97% of reads preserved in the final h5ads (3% dropped = single-read outliers +
unmapped molecules). Stage-1 chromosome partitioning gave a 10.7x memory reduction vs single-process
(1,424 MB -> 133 MB/worker). Stage-3 `var2final` optimization cut peak RSS ~36% (5,528 -> 3,523 MB).

## ITGA2B assessment (cross-sample standardization works)
Same final ids in both samples, per-sample counts ready for differential abundance:
```
Ctrl:6020262   (23-block full-length)   Ctrl 13,926   Aza 11,140
Ctrl:115902251 (1-block E1|I1, real)    Ctrl  3,064   Aza  2,382
Ctrl:7346294   (23-block)               Ctrl  1,670   Aza  1,314
```

## Known limitations / next steps
- **Stage 3 memory floor ~3.5-4.7 GB/sample**: dominated by the AnnData object over an 8.9M-string
  `var` index plus the `var2final` map (almost all molecules survive). To go lower needs gene/chrom
  -chunked h5ad slicing (the viewer's `gene_indexes_v2` already builds a gene-column index that could
  drive this). Stage 1 memory is already solved by chromosome parallelism; Stage 3 chunking is the
  remaining scale lever for 200 samples.
- **Packaging**: this dir (`isoform-collapse/`) is hyphenated and not an importable package, which
  breaks `multiprocessing` pickling of worker functions (worked around via disk hand-off of the
  structure->exemplar map). For production, move to an underscore package (e.g. `isoform_collapse/`).
- Per-molecule h5ad re-key currently re-reads `transcript_associations` per sample; could be cached.

## Files
- `scored_collapse.py` — core per-gene algorithm + invariants (importable).
- `isoform_collapse_utils.py` (parent dir) — `is_contiguous_subsequence`, `exon_block_count`, token helpers.
- `run_pipeline.py` — Stage 1 (parallel) + Stage 2; writes catalog + structure->exemplar map.
- `_stage3.py` — Stage 3 h5ad re-key.
- `_bench_*.py` — benchmark scripts (per-sample, integrated, parallel-by-chromosome).

## End-to-end validation through existing quantification + differential code (measured)
The collapse pipeline output ties into the EXISTING downstream code unchanged:
`run_pipeline` -> final isoform h5ads (var = `gene:exemplar_molecule`, single colon so existing
`split(':')` / `transcript_dict[gene][mol]` lookups work) -> existing `isoform_matrix.pseudo_cluster_counts`
-> `isoform_combined_pseudo_cluster_{tpm,ratio}-filtered.txt` -> existing `comparisons.compute_differentials`.

Measured (Ctrl+Aza, random cluster assignment for evaluation, 3 pseudo-replicates/condition so the
Mann-Whitney test has n=3 vs n=3):
- Quantification: 233,591 final isoforms x (cluster.sample) pseudobulk TPM/ratio matrices; 325 ITGA2B
  isoforms present.
- Differential: `compute_differentials` ran `mwuCompute` per cluster and wrote full stats files
  (`diff-cluster-isoform/*-stats.txt`, 233,591 isoforms x {N, median, mean, DiffMeans, mwuStat,
  mwuSign, mwuPval, mwuAdjPval}). 325 ITGA2B isoforms tested. p-values non-significant as expected
  for RANDOM clusters (validates the machinery, not biology).

### Cross-platform / robustness (verified)
Stage-1 multiprocessing is cross-platform: prefers `fork` on Unix, `spawn` on Mac/Windows, and
AUTO-FALLS-BACK to single-processor mode on ANY multiprocessing failure (verified: parallel, serial
nproc=1, and a forced-failure run all produce BIT-IDENTICAL output, 709,103 isoforms, ITGA2B identical).

## Existing-code finding (NOT modified, per "respect existing code")
`annotation/junction_isoform.annotate_iso_stats_file` reads `row['Feature']`, but
`oncosplice/metadataAnalysis.mwuCompute` writes the isoform ids as the pandas INDEX (unnamed first
column), so `pd.read_csv` yields an `Unnamed: 0` column, not `Feature` -> `KeyError: 'Feature'` at the
ANNOTATION step (after the statistics are already computed and written). This is pre-existing and would
affect any isoform differential run; the differential STATISTICS themselves are complete and correct.
Suggested one-line fix (for the user to apply): in `annotate_iso_stats_file`, read with
`index_col=0` and use the index as the feature, or rename `Unnamed: 0` -> `Feature`. Left for the user.

## Optional terminal intronic evidence protection

`--protect-intronic-evidence` is available for `sclr-isoforms` and `blr-isoforms`;
it is disabled by default. The Python entry point is
`pipeline.run_pipeline(..., protect_intronic_evidence=True)`. Use a separate
output directory with fresh Tier-1 junction-core WTA aggregation. The option
requires the original raw annotation table, primary source GFF exon coordinates,
and the exact reference GFF used for the reference cache. It rejects EM,
sample-level consensus folding and reused Tier-1 tables.

The option preserves terminal intronic and external U-exon splice junctions supported by at least
three pooled source molecules, with a bounded ±3-nt tolerance. Continuous terminal
extension of at least 100 nt into an annotated intron or an actual reference-model
gap is protected independently of that recurrence requirement. Both transcript
ends and strands are supported. An RNA endpoint alone is not a splice junction.
Variable endpoints are pooled rather than creating a class for every coordinate.
For full-cohort runs, the reference driver copies every complete source SQLite
database to node-local storage and verifies its SHA **before** pooled evidence
preparation. Both the initial SQL count checks and gene workers must use those
copies. Copying only before WTA collapse leaves preparation doing slow concurrent
reads on shared storage. Keep the original shared source manifest for cross-node
audits and reverify the complete local cache before reusing it for collapse.
External U protection requires an original own-gene U coordinate at the actual
terminal splice boundary and a terminal exon outside the annotated E/I gene body.
The native preferred Ensembl target is retained whenever the physical evidence
allows that association; evidence splitting does not change ordinary known priority.

Evidence classes form before identical legacy strings are aggregated. The existing
WTA still selects parents, bins and reference priorities. An optional compatibility
check rejects parents that lack the selected splice junction or splice out the
covered intronic sequence. A shorter retention model can remain compatible when
it retains part of the measured continuation. The reciprocal check prevents a
trusted observed splice from being folded into contradictory selected retention.
Upstream annotation and native translation/NMD classification are unchanged.

Original source identities and counts are verified. Repeated library labels
receive distinct physical-source identities and explicit original-file bindings.
Corrected references cannot overwrite observed source annotations when they share
an old structure-cache key. The exact known reference IDs and primary exon
coordinates are retained in emitted GFFs. P4 loads the persisted molecule-specific
classes so count matrices and viewer molecule mappings use the same assignments.

Validation records are in [the review directory](validation/splice_evidence_20261011/).
The default replay exactly reproduces the current saved four-source baseline.
The latest LSP1/RUNX1 local correction adds 36 models to 308,857 (0.01165%).
All 8,878 source reads and 47 references are retained; all 4,834 protected
assignments pass geometry checks. All four matrices preserve ordered cells and
unmodified counts; all 8,186 frozen junction controls are accounted for.
Historical HOPX, exemplar geometry, viewer and complete junction/PSI controls are
documented there. The 101-test suite includes actual full-cohort failure fixtures.

The original native cache strings remain matching/ranking cores. Recovered
coordinates supply physical evidence and export, avoiding artificial E/I matching
transitions. A shared source/reference cache key always receives separate internal
identities. Only exact aggregates that contradict their own selected physical
evidence are refined by observed terminal splice sites; variable RNA endpoints
and unsupported rare sites do not create additional classes.

Earlier full diagnostic candidates were rejected by their audits and are
superseded. The final all-50 KINNEX-5 run and complete downstream audits are active;
use its versioned reports for final counts and acceptance. Two existing terminal-U
splice mismatches remain recorded for discussion before a U-specific extension.
The option does not claim to protect every novel U splice.
