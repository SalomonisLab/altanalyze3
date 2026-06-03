# 4-phase parallel refactor — design

Principle (the thing that should have driven this from the start): **anything computable on ONE
sample in isolation is a per-sample parallel job; only genuine cross-sample reductions are serial.**

Evidence for the phase placement (from this repo's own run logs, sclr_isoforms_test.log):
- per-sample junction quant (`exportJunctionMatrix`) — per-sample inputs only → parallel.
- per-sample junction/isoform pseudobulk (`pseudo_cluster_counts_optimized`) — writes one
  `<sample>.txt`, no cross-sample data → parallel.
- isoform per-sample re-key (stage3) — up to 94.5s, ~19GB RSS/sample → heavy → parallel (NOT bundled
  into the 1-job collapse).
- the collapse CATALOG build (stage1+tier2) — ~92s, 305MB → cheap → one job.
- combine-pseudobulk + PSI + differentials — need all samples (union of features, group sizes) →
  genuinely cross-sample → one job.

## Phases

| Phase | Jobs | Work | Depends on |
|------|------|------|-----------|
| **1 `sclr`** | 36 (per uid) | BAM→molecule h5ad/gff → junction h5ad → **junction pseudobulk** (`<lib>-junction.txt/_tpm/_ratio`) → per-sample collapse inputs (`transcript_associations.txt`). (+optional cellHarmony.) Multi-BAM uid (e.g. AML12) = one job, BAMs merged. | — |
| **2 `sclr-junctions`** | 1 | combine 36 junction pseudobulks → filter → **PSI** → **splicing differentials** | all P1 |
| **3 `sclr-isoforms`** | 1 | combine per-sample isoform models (P1) → two-tier **collapse** (catalog) → **translation + FASTA** | all P1 |
| **4 `sclr-isoquant`** | 36 (per uid) | per-sample **isoform re-key** against P3 catalog → **isoform pseudobulk** (`<lib>-isoform.txt/_tpm/_ratio`) | P3 |
| **4b (in same job or 1 job)** | 1 | combine isoform pseudobulks → **isoform differentials** | all P4 |

(P2 and P3 both depend only on P1, so they can run concurrently. P4 depends on P3. Final isoform
differential needs all P4.)

## Function-level changes (isoform_automate.py)

Split the monoliths into a PER-SAMPLE half and a COMBINE half. Both the CLI and the single-process
driver call the SAME functions — the driver just calls them in a loop in one process.

NEW per-sample functions:
- `process_sample_junctions(sample_entry, ensembl_exon_dir, barcode_sample_dict)` —
  `export_junction_matrix` (writes `<lib>-junction.h5ad`) + per-sample junction pseudobulk
  (`pseudo_cluster_counts_optimized` → `<lib>-junction.txt`). Handles multi-BAM uid by looping the
  uid's entries and concatenating to `<uid>-junction.h5ad`, exactly as `export_junction_h5ad` does
  today (lift that logic here, per-uid).
- `process_sample_isoquant(sample_entry, catalog, ...)` — per-sample isoform re-key (stage3 against
  the final catalog) + per-sample isoform pseudobulk (`<lib>-isoform.txt`).

NEW combine-only functions (cross-sample):
- `combine_junction_pseudobulks(metadata, barcode_cluster_dirs)` — `export_pseudo_counts('junction')`
  CONCAT half only (per-sample `.txt` already exist) + filter.
- `run_splicing(...)` — PSI + splice differentials.
- `combine_isoform_pseudobulks(...)` + isoform differentials.

REFACTOR existing (keep working, same outputs):
- `export_pseudo_counts` already loops per-sample then concats. Add a `skip_existing=True` so when the
  per-sample `<sample>.txt` already exist (built in P1/P4), it goes straight to concat. Single-process
  driver leaves them to be built inline (skip_existing still works — builds then concats).
- `pre_process_samples` (single-process driver entry) = call `process_sample_junctions` for every
  sample, then `combine_junction_pseudobulks`, then `run_splicing`. Identical outputs to today.
- `combine_processed_samples` = collapse (catalog) [P3] + per-sample `process_sample_isoquant` loop
  [P4] + combine + isoform diff. Identical outputs to today.

So: single-process driver = the same per-sample functions called in a serial loop; CLI = the same
functions fanned out across bsub jobs. One code path, two schedulers.

## CLI (cli.py)

- `run_sclr` (P1): after extraction (+cellHarmony), call `process_sample_junctions` for the sample.
- `run_sclr_junctions` (P2): `combine_junction_pseudobulks` + `run_splicing` only. No per-sample work.
- `run_sclr_isoforms` (P3): collapse catalog + translation/FASTA only (no per-sample pseudobulk).
- NEW `run_sclr_isoquant` (P4): `process_sample_isoquant` for `--sample` + (final job) combine +
  isoform diff. Same `--sample`-or-all pattern as `run_sclr`.

## Submit script

P1: 36 bsub `sclr --sample <uid>`  -> DEP_P1
P2: 1 bsub `sclr-junctions`  -w DEP_P1
P3: 1 bsub `sclr-isoforms`   -w DEP_P1
P4: 36 bsub `sclr-isoquant --sample <uid>`  -w done(sclr_isoforms)  -> DEP_P4
P4-diff: 1 bsub `sclr-diff`  -w DEP_P4   (isoform differentials; splice diff already in P2)

## Validations to keep
- read-conservation check in stage3 (already present).
- per-sample pseudobulk file existence gate before each combine (loud if a sample's `.txt` missing).
- completeness gate: every uid has its expected per-sample outputs before a combine job.
