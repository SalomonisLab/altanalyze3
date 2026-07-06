# fastCNV candidate workflow

This directory is intentionally isolated from the production `fastCNV.py`.
Nothing here should be imported by the main CLI until it passes the benchmark
gates in `/Users/saljh8/Dropbox/Revio/fastCNV/SUCCESS_CRITERIA.md`.

## Current candidate: `per_gene.py`

`per_gene.py` is the only candidate workflow. LOY, AML −7, +8, del(5q),
balanced rearrangements, and normal controls must all be evaluated from this
same output. There should not be a separate LOY caller, separate LOY threshold,
or benchmark-specific decision path.

Design:

- Python-only, vectorized per-gene fold-change workflow.
- Cell-state-specific healthy reference centroids.
- Genes with low reference detection or low reference expression receive zero
  weight.
- Regional statistics are baseline-expression-weighted aggregates of per-gene
  log2 fold-change in CP10k space.
- The primary CNV unit is a genomic interval. Sliding windows are used only as
  fast/noise-robust detection scaffolds; called runs are refined back to the
  first/last altered informative gene before writing `candidate_cnv_regions`.
  Global/shared CNVs should therefore be evaluated as overlapping intervals
  across cell states, not as requiring identical arm/chromosome boundaries.
- Calls are gated by a healthy-reference null distribution and simulated copy
  states. The same logic is used for haploid and diploid chromosomes; chrY is
  treated as haploid only when selected as a haploid chromosome, not through a
  separate LOY caller.
- Sparse chromosome outputs can optionally be augmented by unsupervised
  state/chromosome two-component mixture calibration on the same per-cell
  scores. This is implemented, but current LOY tests show it is not yet robust
  enough to claim success.
- Count layers are normalized by full-cell library size before gene subsetting,
  so chromosome-restricted runs do not inflate sparse chromosomes.
- Outputs:
  `*.candidate_cnv_cells.tsv`,
  `*.candidate_cnv_regions.tsv`,
  `*.candidate_global_summary.tsv`,
  `*.candidate_reference_summary.tsv`.
- `candidate_cnv_regions` includes `boundary_source`,
  `n_informative_genes`, `n_evidence_genes`, and
  `evidence_gene_fraction` to make each interval auditable.
- Use `interval_consensus.py` to collapse per-cell intervals into sample-level
  consensus intervals. This reports overlapping support from cell states while
  preserving start/end boundary variability.

## Current status

Normal-karyotype smoke test:

```bash
python -m altanalyze3.components.fastCNV.candidate.per_gene \
  --h5ad /Users/saljh8/Dropbox/Revio/fastCNV/aml_benchmark/samples/AML0160.h5ad \
  --control-h5ad /Users/saljh8/Dropbox/Revio/fastCNV/aml_benchmark/control_reference.h5ad \
  --output-prefix /tmp/fastcnv_candidate_smoke/AML0160 \
  --state-key Hs-BM-titrated-reference-centroid \
  --control-state-key Hs-BM-titrated-reference-centroid \
  --sample-key Sample \
  --max-cells 80
```

Result after adding reference-null and copy-state gates:

- Runtime: ~4.5 s for 80 cells.
- Negative-control smoke status: **PASS for this subset**.
- Current overcall: 0/80 AML0160 normal-karyotype cells called CNV.

LOY chrY benchmark command:

```bash
python -m altanalyze3.components.fastCNV.candidate.per_gene \
  --h5ad /Users/saljh8/Dropbox/Collaborations/Grimes/MDS-AK-CITE-Seq-LOY/MDS/MDS-samples-only/pooled/fastCNV.input.h5ad \
  --control-h5ad /Users/saljh8/Dropbox/Revio/fastCNV/reference/healthy_marrow_reference.corrected.aug100.h5ad \
  --output-prefix /tmp/fastcnv_candidate_v2/LOY_chrY_z4_m3 \
  --state-key Hs-MarrowAtlas-L3M \
  --control-state-key Hs-MarrowAtlas-L3M \
  --sample-key Library \
  --layer counts \
  --counts-input \
  --control-filter-field sex \
  --control-filter-value Male \
  --chromosomes chrY \
  --min-chr-genes 10 \
  --window-genes 41 \
  --stride-genes 7 \
  --chunk-cells 2048
```

Current LOY result using the high-specificity defaults:

- Runtime: ~5.3 s for 80,283 cells restricted to chrY.
- Calls: ~9,790 chrY-loss cells.
- Compared with the prior heuristic labels: conservative, not passing the
  requested 13k-14k global LOY range or the requested per-state
  sensitivity/specificity gate.

Lowering `--min-copy-z` increases LOY recall but produces discontinuous
state-level behavior: `3.78` calls 15,238 cells, while `3.795` calls 9,912
cells. This indicates a single global scalar threshold is not a defensible final
solution.

Mixture calibration status:

- `--mixture-min-separation-z 2.5`: conservative, ~9,790 chrY-loss cells.
- `--mixture-min-separation-z 2.0`: conservative, 10,421 chrY-loss cells.
- `--mixture-min-separation-z 1.8`: overcalls, 21,455 chrY-loss cells.
- `--mixture-min-separation-z 1.5`: severe overcall, 35,582 chrY-loss cells.

## Required next fixes before any success claim

- Improve state/chromosome mixture calibration so it estimates retained-vs-
  altered components per state without using LOY labels and without a separate
  LOY path. The current two-means implementation is too discontinuous.
- Evaluate held-out controls beyond the 80-cell AML0160 smoke subset.
- Add simulation tests that scale reference cells by deletion/gain factors and
  report sensitivity/specificity against known synthetic truth.
- Only then benchmark LOY and known AML cytogenetic lesions from the same
  workflow outputs.
