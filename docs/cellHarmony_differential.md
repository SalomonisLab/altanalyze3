# cellHarmony_differential.py

## Overview
`cellHarmony_differential.py` tests reference-aligned single-cell datasets for differential expression between conditions. The script accepts a cellHarmony-aligned `.h5ad`, merges external covariates (e.g. sample or disease labels), performs per-population comparisons, and summarises global, local, and co-regulated gene programs. Optional pseudobulk aggregation, co-regulation heatmaps, and interaction-network exports are included to streamline downstream interpretation.

## Required Inputs
- **Aligned `.h5ad`**: Output from `cellHarmony_lite.py` (or another AnnData object) with projected populations in `.obs`.
- **Covariate table (TSV)**: Contains library IDs, sample IDs, and the covariate/condition label used in comparisons.
- **Reference population column**: Name of the `.obs` column holding cellHarmony population assignments.
- **Comparisons list**: Semi-colon separated `CASE|CONTROL` pairs defining the contrasts to evaluate (e.g. `AML|Healthy;Relapse|Diagnosis`).

The covariate table must include the columns referenced by `--library_col`, `--sample_col`, and `--covariate_col`. Library IDs must match the values stored in `.obs[Library]` unless overridden.

---

## Core Command-Line Parameters
These options cover the typical differential workflow. Defaults appear in parentheses.

| Parameter | Description |
|-----------|-------------|
| `--h5ad` | Input AnnData with cellHarmony-assigned populations. |
| `--covariates` | TSV mapping library IDs to sample IDs and covariate/condition labels. |
| `--population_col` | `.obs` column containing cellHarmony populations used to stratify DE (e.g. reference name). |
| `--comparisons` | Semi-colon separated list of `CASE|CONTROL` pairs to test. |
| `--covariate_col` (`Condition`) | Column in the covariate table (and merged `.obs`) representing the condition label. |
| `--library_col` (`Library`) / `--sample_col` (`Sample`) | Column names used to join covariates to the AnnData object. |
| `--method` (`wilcoxon`) | Scanpy `rank_genes_groups` method (`wilcoxon`, `t-test`, `t-test_overestim_var`, `logreg`). |
| `--alpha` (`0.05`) | Significance threshold (FDR or raw p-value depending on `--use_rawp`). |
| `--fc` (`1.2`) | Absolute fold-change threshold applied alongside `--alpha`. |
| `--min_cells_per_group` (`20`) | Minimum cells per condition within each population (auto-relaxed for small datasets). |
| `--outdir` (`cellHarmony_DE_out`) | Destination for reports, heatmaps, AnnData exports, and interaction plots. |

---

## Accessory Parameters
Additional controls are grouped by theme to keep the core workflow concise.

### Pseudobulk Configuration
| Parameter | Description |
|-----------|-------------|
| `--make_pseudobulk` | Aggregate cells into (population × sample) pseudobulks before testing (limma-style moderated t-test). |
| `--pseudobulk_min_cells` (`10`) | Minimum cells required to form each pseudobulk entry. |

### Statistical Behaviour
| Parameter | Description |
|-----------|-------------|
| `--use_rawp` | Use raw p-values (instead of FDR) when applying `--alpha` for significance and downstream filtering. |

### Output Controls & Extras
| Parameter | Description |
|-----------|-------------|
| `--skip_grn` | Disable interaction-network (GRN) exports based on differential gene sets. |

---

## Generated Outputs
All files are written into `--outdir`. For each comparison `CASE_vs_CONTROL` the script produces:

| File | Description |
|------|-------------|
| `DEG_detailed_CASE_vs_CONTROL.tsv` | Per-population gene statistics (log2FC, p-values/FDR, counts). |
| `DEG_summary_CASE_vs_CONTROL.tsv` | Summary of cell counts, DEG counts, and tested genes per population. |
| `DEG_assigned_groups_CASE_vs_CONTROL.tsv` | Gene classifications (global, local, co-regulated) with key populations/patterns. |
| `DEG_pooled_overall_CASE_vs_CONTROL.tsv` | Optional pooled “global” differential results (if sufficient cells). |
| `DEG_coreg_pooled_CASE_vs_CONTROL.tsv` | Pooled results for top co-regulation patterns. |
| `heatmap_CASE_vs_CONTROL_by_POPULATIONCOL.pdf/tsv` | Fixed-order co-regulation heatmap plus underlying matrix. |
| `differentials_only_CASE_vs_CONTROL.h5ad` | AnnData containing only the cells used in DE with DEG annotations in `.uns['cellHarmony_DE']`. |
| `interaction-plots/CASE_vs_CONTROL/<population>.*` | Per-population gene regulatory networks (PDF, PNG, TSV) generated via NetPerspective, unless `--skip_grn` is set. |

Additional artefacts include:
- `unsupervised_leiden_clusters.tsv` and `metacells.h5ad` when those features are triggered upstream in the `.h5ad` supplied.
- `pseudobulk` AnnData (`*_pseudobulk.h5ad`) when pseudobulk mode is enabled.

---

## Interaction Network Generation
By default, significant up/down genes per population feed into the NetPerspective workflow, leveraging curated TF–target, protein–protein, and signalling interactions from `components/visualization/interactions`. PDF/PNG plots and a TSV summary of aggregated targets are emitted for each population with ≥2 significant genes. Use `--skip_grn` to disable this step or install `python-igraph` to enable plotting locally.

---

## Example Workflow
```bash
python cellHarmony_differential.py \
  --h5ad results/cellHarmony_aligned.h5ad \
  --covariates metadata/patient_annotations.tsv \
  --population_col AnnotatedCellTypes \
  --comparisons "AML|Healthy;Relapse|Diagnosis" \
  --method wilcoxon \
  --alpha 0.05 --fc 1.3 \
  --make_pseudobulk --pseudobulk_min_cells 15
```

1. Cells and covariates are merged by library ID, ensuring population assignments and sample labels are valid.
2. Optional pseudobulk profiles are created per (population × sample) when requested.
3. Differential expression is computed per population for each comparison, with global and co-regulated summaries.
4. Heatmaps, DEG tables, AnnData subsets, and interaction plots (if enabled) are written to `cellHarmony_DE_out`.

Adjust the accessory switches to refine statistical behaviour, skip network exports, or tailor pseudobulk requirements for sparse datasets.
