# fastComm Evaluation

This document records the current prototype evaluation commands, observed
runtime, and benchmark interpretation for `fastComm`.

## Current Prototype Dataset

Local benchmark dataset:

```text
/Users/saljh8/Dropbox/Manuscripts/InProgress/Human-Titration/ShinyCell/adata_combined_rna_adt_annotated-titrated.subset50.male.h5ad
```

Observed input shape:

```text
7,919 cells x 36,743 features
```

The prototype loads only genes required for ligand-receptor and receiver
response scoring. For the current seed LR/response resources, this was `74`
genes.

## Level 2 Scoring Run

Command:

```bash
PYTHONPATH=.. /opt/homebrew/opt/python@3.11/bin/python3.11 \
  -m altanalyze3.components.fastComm.cli score \
  --h5ad /Users/saljh8/Dropbox/Manuscripts/InProgress/Human-Titration/ShinyCell/adata_combined_rna_adt_annotated-titrated.subset50.male.h5ad \
  --state-key 'Level 2' \
  --min-cells 30 \
  --exclude-self-edges \
  --output components/fastComm/artifacts/prototype_male_level2/fastcomm_scores.tsv \
  --state-pair-output components/fastComm/artifacts/prototype_male_level2/state_pair_summary.tsv \
  --state-expression-output components/fastComm/artifacts/prototype_male_level2/state_expression.tsv \
  --summary-json components/fastComm/artifacts/prototype_male_level2/summary.json
```

Observed summary:

```text
n_cells: 7,919
n_loaded_genes: 74
n_states: 27
n_lr_edges: 27
n_scored_edges: 3,391
```

Top interaction:

```text
Stroma -> Ba/Ma/Eo : CXCL12 -> CXCR4
prototype_probability: 0.590
```

## Exemplar Report

Command:

```bash
PYTHONPATH=.. /opt/homebrew/opt/python@3.11/bin/python3.11 \
  -m altanalyze3.components.fastComm.cli exemplar-report \
  --scores components/fastComm/artifacts/prototype_male_level2/fastcomm_scores.tsv \
  --output-tsv components/fastComm/artifacts/prototype_male_level2/exemplar_interactions.tsv \
  --output-md components/fastComm/artifacts/prototype_male_level2/exemplar_interactions.md \
  --split-stability components/fastComm/benchmarks/prototype_male_level2_donor/split_stability.tsv \
  --top-n 20 \
  --min-score 0.25 \
  --title 'fastComm Exemplar Interactions: Human Marrow Male Subset Level 2'
```

Top exemplar interactions:

```text
Stroma -> Ba/Ma/Eo          CXCL12 -> CXCR4        0.590
Stroma -> Transitional-B    CXCL12 -> CXCR4        0.550
Stroma -> B cell            CXCL12 -> CXCR4        0.527
Stroma -> pre-B             CXCL12 -> CXCR4        0.501
Mac -> Ba/Ma/Eo             CXCL12 -> CXCR4        0.470
Mac -> GMP                  VCAM1 -> ITGA4+ITGB1   0.433
Stroma -> GMP               VCAM1 -> ITGA4+ITGB1   0.427
T/NK -> Stroma              TNF -> TNFRSF1A        0.372
B cell -> GMP               IL6 -> IL6ST           0.349
```

Interpretation:

- `prototype_probability` is a 0-1 confidence-like score, not yet calibrated
  against perturbation or spatial benchmarks.
- `pathway_or_factor` and `supporting_response_genes` come from the active
  receiver-response matrix.
- Current factors are seed prototype signatures, not yet learned from CytoSig,
  NicheNet, Immune Dictionary, or perturbation atlases.

## Donor Split Stability Benchmark

Command:

```bash
PYTHONPATH=.. /opt/homebrew/opt/python@3.11/bin/python3.11 \
  -m altanalyze3.components.fastComm.cli benchmark-h5ad \
  --h5ad /Users/saljh8/Dropbox/Manuscripts/InProgress/Human-Titration/ShinyCell/adata_combined_rna_adt_annotated-titrated.subset50.male.h5ad \
  --state-key 'Level 2' \
  --split-key Donor \
  --min-cells 30 \
  --output-dir components/fastComm/benchmarks/prototype_male_level2_donor
```

Observed stability:

| split | n_cells | n_states | n_edges | elapsed_seconds | top_interaction | top_100_jaccard_vs_full | score_corr_vs_full |
| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| BM27 | 3,893 | 24 | 2,768 | 0.3796 | Stroma->Transitional-B:CXCL12->CXCR4 | 0.481 | 0.947 |
| WM34 | 4,026 | 27 | 3,250 | 0.4403 | Stroma->Ba/Ma/Eo:CXCL12->CXCR4 | 0.600 | 0.975 |

Interpretation:

- Rank correlation against the full run is strong for both donors.
- Top-edge Jaccard is moderate, as expected for donor-specific expression and
  state-composition differences.
- Runtime shown here is scoring time after the reduced expression matrix is
  loaded; full CLI wall time on this machine was a few seconds.

## Level 3 Multimodal Stress Run

Command:

```bash
PYTHONPATH=.. /opt/homebrew/opt/python@3.11/bin/python3.11 \
  -m altanalyze3.components.fastComm.cli score \
  --h5ad /Users/saljh8/Dropbox/Manuscripts/InProgress/Human-Titration/ShinyCell/adata_combined_rna_adt_annotated-titrated.subset50.male.h5ad \
  --state-key 'Level 3 Multimodal' \
  --min-cells 40 \
  --exclude-self-edges \
  --top-n 5000 \
  --output components/fastComm/artifacts/prototype_male_l3m/fastcomm_top5000.tsv \
  --state-pair-output components/fastComm/artifacts/prototype_male_l3m/state_pair_summary.tsv \
  --summary-json components/fastComm/artifacts/prototype_male_l3m/summary.json
```

Observed summary:

```text
n_cells: 7,919
n_loaded_genes: 74
n_states: 80
n_lr_edges: 27
n_scored_edges: 5,000 top retained
```

Top interaction:

```text
Stromal Vascular -> BMCP-2 : CXCL12 -> CXCR4
prototype_probability: 0.671
```

## Required Next Benchmarks

The current evaluation is a prototype smoke/stability benchmark. Before
claiming calibrated biological probability, `fastComm` needs:

- perturbation recovery: rank true ligand/cytokine in receiver response data
- receptor-dependence: score drop after receptor or receptor-subunit loss
- spatial plausibility: contact LR enrichment among adjacent cells
- protein support: improved confidence when receptor CITE-seq signal agrees
- runtime scaling: 10k, 100k, and 1M cell datasets with production LR resources
- cross-dataset generalization: hold out tissues, donors, and ligand families

## Test Command

```bash
PYTHONPATH=.. pytest -q components/tests/test_fastcomm_basic.py
```

Current result:

```text
6 passed
```

## Mouse Compatibility Check

The unit suite includes a synthetic mouse `.h5ad` regression test using
title-case mouse symbols:

```text
Tgfb1, Tgfbr1, Tgfbr2, Serpine1, Smad7
```

The test confirms these map to the canonical LR/response symbols:

```text
TGFB1, TGFBR1, TGFBR2, SERPINE1, SMAD7
```

Expected result:

```text
n_matched_required_genes: 5
n_missing_required_genes: 0
top interaction: TGFB1 -> TGFBR1+TGFBR2
receiver_response_score > 0.9
```

This verifies mouse title-case datasets can run through the current workflow.
For production mouse analyses, species-specific LR resources remain preferable
for non-conserved interactions and mouse-specific response programs.
