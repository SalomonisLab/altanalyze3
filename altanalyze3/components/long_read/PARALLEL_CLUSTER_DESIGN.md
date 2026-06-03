# Long-read single-cell workflow: parallel / cluster CLI design

Status: design (no implementation yet). The existing single-process driver scripts
(`run_full_workflow.py`, `run_workflow_simple.py`) are unchanged and remain the recommended
path for small cohorts (~10 samples). This adds an opt-in CLI for large cohorts (~200 samples)
where per-sample BAM processing dominates wall time.

## Principle

The per-sample work (phase 1) is embarrassingly parallel; the cross-sample work (phases 2-4) is
not. Industry standard (STAR/Salmon per-sample -> tximport; GATK HaplotypeCaller GVCF -> joint
GenotypeGVCFs; cellranger count -> aggr) is: expose a clean **per-sample command** and a clean
**aggregate command**, and let a shell loop / scheduler fan out the per-sample jobs. We follow
that: one per-sample subcommand, three integration subcommands, all idempotent so jobs are
re-submittable after a node failure with no recompute.

The two paths produce identical outputs. Parallel mode is a scheduling wrapper around the same
functions; no algorithm changes.

## The 4 phases -> CLI subcommands

| Phase | Subcommand | Parallel? | Backing function |
|------|------------|-----------|------------------|
| 1. BAM -> per-sample junction h5ad (+ optional cluster labels) | `sclr` | **YES, one job/sample** | `import_metadata(extract_from_bams=True)` + `export_junction_matrix` (+ optional cellHarmony) |
| 2. junction aggregation + PSI + differential splicing | `sclr-junctions` | no (single job) | `pre_process_samples` |
| 3. isoform collapse + per-sample isoform h5ads | `sclr-isoforms` | no (single job) | `combine_processed_samples` |
| 4. differential isoform analysis | `sclr-diff` | no (single job) | `compute_differentials` |

Why this split is mostly scheduling, not a rewrite (facts already in the code):
- `import_metadata(extract_from_bams=False)` already assumes BAMs were processed and rebuilds each
  sample's gff/matrix path from `<bam_parent_dir>/<library>` -> **BAMs in different directories
  already work**, outputs land beside each BAM.
- `export_junction_matrix` short-circuits if `<library>-junction.h5ad` exists -> integration
  **reloads** per-sample caches instead of recomputing.
- `reverse` is already per-metadata-row, consumed as `rev=reverse` -> **never goes on the command
  line**; it travels in the metadata every job reads.

## Phase 1: `altanalyze3 sclr`

```
altanalyze3 sclr \
    --metadata samples.txt          # required; master (or any) metadata file
    --sample 251H                   # optional; OMIT -> serial loop over ALL uids in the metadata
    --ref_gff gencode.v49.gtf.gz    # required; GENCODE GFF (run through gff_process)
    --species human                 # human|mouse, default human; selects bundled defaults + ref registry species id
    --exon_annot <file>             # optional; default = bundled gzipped Hs_Ensembl_exon.txt (see Defaults)
    --gene_symbol <file>            # optional; default = bundled gzipped Hs_Ensembl-annotations.txt
    --cellHarmony_ref <id|path>     # optional; ALIGN to this reference (registry id or centroid .txt path)
    --cell_annot <file>             # optional; USE existing barcode->cluster annotations (cellHarmony format)
    --threads 2 --cpus 12           # existing common args
```

Per-sample steps: BAM extract -> gene aggregate -> (optional cluster labels) -> per-sample
junction h5ad + pseudobulk-ready matrix. Idempotent: a sample whose `-junction.h5ad` already
exists is skipped.

- `--sample` omitted => **serial loop** over every uid in the metadata (the convenience / small-N
  fallback). Parallelism comes from the submit script running one `sclr --sample <uid>` per bsub.
  Each sample still uses `--cpus` internally for its own extraction.

### Cluster labeling: two mutually exclusive inputs

1. **`--cellHarmony_ref <id|path>`** — *run* cellHarmony-lite to produce barcode->cluster.
   Resolution (reuse the existing web resolver, do not invent a parallel registry):
   - If the value matches a reference **id** in `cellHarmony/flask/reference_config.json` for the
     selected `--species` (e.g. `hs_bm_reference`, `hs_lung_hlca_reference`, `mm_bm_reference`),
     resolve via `flask.pipeline._lookup_reference(species, id, registry_path)` -> its
     `states_tsv` (centroid file). This exposes ALL of cellHarmony's stored references.
   - Otherwise treat the value as a **file path** to a user-provided centroid `.txt`.
   Then call `run_cellHarmony_lite.run_cellharmony_lite(gene_h5ad, ..., cellharmony_ref=<resolved>,
   gene_translation_file=<--gene_symbol>)` -> writes `<library>_barcode_clusters.txt`.
2. **`--cell_annot <file>`** — *use* a pre-existing barcode->cluster file already in cellHarmony
   format (need not have come from cellHarmony). Skips alignment; this file becomes the sample's
   `barcode_cluster_dir` directly.

`--cellHarmony_ref` and `--cell_annot` are mutually exclusive. If neither is given, phase 1
attaches no cluster labels (they can be supplied to the integration phase).

## Phases 2-4: integration subcommands (each its own CLI call / bsub job)

```
altanalyze3 sclr-junctions --metadata samples.txt --exon_annot <file> [--barcode_clusters <dir>] ...
altanalyze3 sclr-isoforms  --metadata samples.txt --ref_gff <gencode> --genome_fasta <fa> ...
altanalyze3 sclr-diff      --metadata samples.txt [--method mwu|limma] ...
```

- `sclr-junctions` -> `pre_process_samples(metadata, barcode_cluster_dirs, exon_annot)`
  (junction aggregation + pseudobulk + PSI + differential splicing). Reloads phase-1 caches.
- `sclr-isoforms` -> `combine_processed_samples(...)` (two-tier collapse + isoform h5ads).
- `sclr-diff` -> `compute_differentials(..., method=...)` (default mwu; optional limma/eBayes).

## Defaults bundled in the repo (gzipped)

- `--exon_annot` default: a gzipped `Hs_Ensembl_exon.txt` (and `Mm_Ensembl_exon.txt` for
  `--species mouse` when added) packaged in a repo folder. Equivalent content is produced by
  `altanalyze3 index --gtf <gencode> --output <dir>` (writes `gene_model_all.tsv`, same E#.#/I#.#
  exon-segment model); a user who wants a different build passes `--exon_annot <their file>`.
  NOTE: the bundled `Hs_Ensembl_exon.txt` is headered (cols `gene, exon-id, chromosome, strand,
  exon-region-start(s), exon-region-stop(s)`, read by `parse_exon_file` via DictReader); the
  `index`-built `gene_model_all.tsv` is headerless with the same column order -- if a user points
  `--exon_annot` at an `index` output, sclr reads it positionally.
- `--gene_symbol` default: a gzipped `Hs_Ensembl-annotations.txt` (Ensembl-id -> symbol), used to
  translate BAM Ensembl gene ids to the symbol-keyed cellHarmony references.
- `--species` selects which bundled defaults + which registry `species` id to use. Human bundled
  now; mouse files added when available.

## Orchestration (submit script, user-run)

`submit_longread_cluster.sh` (written into the project/output dir, not /tmp):

```sh
#!/bin/sh
META=$1; REF_GFF=$2          # master metadata, gencode gff
# Phase 1: one bsub per sample. reverse + bam-dir come from the metadata row.
awk -F'\t' 'NR>1{print $1}' "$META" | sort -u | while read uid; do
  bsub -J "sclr_$uid" -n 12 -M 32000 \
    altanalyze3 sclr --metadata "$META" --sample "$uid" --ref_gff "$REF_GFF" \
      --species human --cellHarmony_ref hs_bm_reference --cpus 12
done
# Phases 2-4: run AFTER all sclr_* finish (manual two-phase, or add bsub -w "done(sclr_*)").
altanalyze3 sclr-junctions --metadata "$META"
altanalyze3 sclr-isoforms  --metadata "$META" --ref_gff "$REF_GFF"
altanalyze3 sclr-diff      --metadata "$META"
```

Manual two-phase per the user's preference: run the loop, wait for `sclr_*` jobs to clear, then
run the three integration commands. (LSF `-w "done(sclr_*)"` can gate them automatically if
wanted.) BAMs in different directories and per-sample reverse-complement are handled entirely by
the metadata rows, so no per-sample paths/flags appear on the command line.

## Validations built into the code (standing rule)

- **Phase-1 completeness gate** before any integration command: assert every uid in the metadata
  has its `<library>-junction.h5ad` (and `.h5ad`) on disk; fail loudly listing missing samples
  rather than silently integrating a partial cohort.
- **Mutual exclusivity**: error if both `--cellHarmony_ref` and `--cell_annot` are given.
- **Reference resolution**: if `--cellHarmony_ref` is neither a known registry id (for the chosen
  species) nor an existing file, error with the list of valid ids.
- **Idempotency**: phase-1 skips samples whose outputs already exist.

## Parser wiring (utilities/parser.py)

Add four subparsers mirroring the existing pattern (`set_defaults(func=...)`, `add_common_arguments`,
per-command assert in `assert_args`): `sclr`, `sclr-junctions`, `sclr-isoforms`, `sclr-diff`.
Resolve `--ref_gff/--exon_annot/--gene_symbol/--cell_annot/--metadata` paths in `resolve_path`.
`func` targets thin wrappers in a new `components/long_read/cli.py` that call the existing
`isoform_automate` / `comparisons` functions; `--cellHarmony_ref` is resolved through the existing
`cellHarmony/flask/reference_config.json` registry.
```
