# Noncanonical long-read integration

This branch adds a complete-workflow entry point with standard analyses, same-chromosome noncolinear/circRNA discovery, the validated AML/MDS fusion–KMT2A PTD–FLT3 ITD profile, and inter-chromosomal extraction/classification. The October 16 deployment is preserved as a separate comparison baseline.

```bash
python -m altanalyze3 lr-integrate --config /absolute/path/noncanonical.json
# Equivalent API:
# isoform_automate.run_integrated_workflow('/absolute/path/noncanonical.json')

# Source-local execution and explicit barriers:
python -m altanalyze3 lr-integrate --config /absolute/path/noncanonical.json --stage extract --source source006
python -m altanalyze3 lr-integrate --config /absolute/path/noncanonical.json --stage aggregate
python -m altanalyze3 lr-integrate --config /absolute/path/noncanonical.json --stage standard

# Submit the same stages to LSF, from an environment available on compute hosts:
python -m altanalyze3.components.long_read.noncanonical.submit --config /absolute/path/noncanonical.json --queue rhel9
```

`modalities` defaults to `standard,circrna,oncofusion,trans-splicing`. Explicit selections are recorded in the manifest. Individual existing `sclr*` phase commands retain their original behavior; `sclr --noncanonical-config ...` adds source evidence extraction. The full cluster launcher defaults to the integrated graph through `NONCANONICAL_CONFIG`; `LEGACY_ONLY=1` explicitly selects its established standard-only graph. Missing inputs never silently disable a requested complete-workflow branch.

## Manifest

Each ordered `sources` entry requires a unique stable `source_id`, the original `uid`, `library`, `group`, barcode `reverse` orientation, BAM path/SHA256 and index path/SHA256. Repeated UIDs and library labels are legitimate. Physical sources remain distinct even when query names or BAM content match.

Reference fields are `exon_model`, `genome`, `reference_gff`, `gene_symbols`, `reference` (the validated trans-splicing resource directory), `hgnc`, `barcode_clusters`, `patients`, and `expression_dir`. Classification uses the explicitly pinned compatibility expression view, so classification does not depend on its own count reallocation. `reference_pins` maps resource paths to full byte SHA256 values; `code_sha256` is `manifest.code_digest()`. `runtime_versions` can bind the native runtime versions. Input/source hashes are verified before execution and completion receipts bind outputs to the configuration and code.

The oncofusion branch also requires `samtools` on `PATH`, with its executable SHA256 and first `--version` line in `runtime_executables.samtools`. The validated cluster profile loads `python3/3.12.3 samtools/1.21.0`; the module reports `samtools 1.21`. Missing or changed executables fail preflight. The LSF submitter pins child processes to one CPU and a deterministic Python hash seed.

`standard_parameters` explicitly records `collapse_method`, `min_total`, `protect_intronic_evidence`, `protect_novel_junctions`, `analyses`, `method`, and `conditions`. The deployed profile uses native WTA, reviewed intronic/novel-junction protection, and the native Python `limma` implementation. Test subsets must declare their source roster and test contrasts independently. This manifest does not infer donor identities or authorize new cohort decisions.

The validated oncofusion profile remains the 68-gene AML/MDS panel and FLT3 exon 14–15 scan. It is not a genome-wide pan-cancer PTD/ITD detector. `mappy==2.30` is a core dependency; the separately scoped exact-sequence utilities require the `noncanonical-sequence` extra and their own admitted inputs.

## Outputs and count semantics

- `compatibility/`: native standard extraction, catalog, quantification, PSI and differential outputs.
- `circrna/`: all raw noncolinear pairs, barcode coverage denominators and lossless event aggregation.
- `oncofusion/samples/`: source-qualified raw scanner output. `oncofusion/native_profile/` restores the original `uid__library` grouping for the unchanged native reporting rules and writes a reversible source map.
- `trans_splicing/`: all raw path joints, native WTA catalog/ordered-segment translation/quantification, complete source statuses, reporting and classification. The native catalog threshold is 3; sources with no accepted catalog or annotated cells retain their raw matrices and explicit status.
- `reconciliation/`: SQLite evidence, original read/segment identities, all event support/statuses, assignments, exact count deltas, translation caveats and conservation summary.
- `integrated/standard/`: reconciled native standard analyses on the complete compatibility catalog. Raw matrices retain `legacy_counts`; isoform ratios retain compatibility gene denominators and the original filtered ratio feature universe.
- `integrated/events/`: complete declared annotated cells plus observed source/barcode identities, UID/cell-state event-count matrices, status layers and unresolved cell mappings. Annotated cells with no observed reads remain explicit zero rows. Features and cells are retained when their integrated standard count is zero.

Each event counts distinct `(source_id, query_name)` support. Repeated crossings are retained as evidence rows but count once per event. Matrices overlap intentionally; their totals must not be summed as unique-read totals. Original native oncofusion reports retain their historical counting semantics separately.

The first routing policy subtracts a conventional whole-read contribution once when kept cross-chromosome evidence receives a native fusion or trans-splicing candidate classification. Compatible variant identities and unresolved PTD/circle evidence retain their conventional contribution. Ordinary local splice-junction and gene-expression products retain the compatibility counting view. Every subtraction must match the exact original barcode/feature count; otherwise execution stops with an identity diagnostic.

KMT2A alternate exon 2 is retained. Native PTD/circle evidence classes and FLT3 cluster membership are explicit; an individual tandem score cannot replace the native cluster median and modal insertion length. Candidate annotations do not establish a genomic rearrangement or physical RNA circularity. Reference exon-chain proteins do not resolve inserted sequences or ambiguous circular/PTD topology; affected reads receive translation caveats.

## Validation boundaries

The official 10x barcode recognizer is corrected to accept an existing numeric lane suffix such as `-1`, so the requested reverse complement reaches the standard junction path. Regression checks preserve all counts and leave bulk identifiers unchanged. The raw annotation writer is also closed before compression, preventing buffered source rows from being lost. The independent standard comparison records these two one-line compatibility fixes separately from the unmodified official source, while retaining the complete-roster gate.

`COMPLETE.json` is a completion receipt for the requested detection, classification, counting and configured standard-analysis stages. It is not an assertion of biological validation. The package retains the accepted full-query/semiglobal/admission utilities with their existing gates, and the native circRNA statistical adapter retains its explicit decision gates. These deeper sequence and cohort-statistical utilities are separately scoped; discovery does not imply every event was assessed by them. The current integrated event outputs are descriptive counts, not new differential circRNA or host-gene-ratio P-values.

Tests cover source collisions, ambiguous/coexisting identities, repeated crossings, exact BAM-to-matrix reconciliation, one subtraction for multiple accepted joints, unchanged gene denominators, stale-code/output rejection, alternate KMT2A exon 2 and native renderer compatibility. Native real-genome circRNA and cross-chromosome rule tests are executable separately. Cluster regression scripts additionally rerun the original implementations on identical source BAMs and compare outputs, then compare the native standard and cross-chromosome catalog workflows.

The November 11 cluster fixture contains six physical sources representing five analysis UIDs, all records overlapping 13 declared loci, and recursively recovered exact SA-linked alignments. It preserves 2,128,598 source-qualified queries and 2,209,567 alignment records. Its configured standard contrasts exercise within-aged cell-state comparisons; the full 50-source deployment and complete cohort contrasts are outside this regression. The pinned archived compatibility expression matrices remain the classifier inputs. A full-cohort run requires its own complete source roster, references and established contrasts in the manifest.

`SOURCE_PROVENANCE.json` records the accepted original files and explicit adapter/bug-fix changes. A package change invalidates old receipts; a failed or incomplete stage cannot count as a successful empty result.
