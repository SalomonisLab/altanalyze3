# fastCNV Candidate Workflow

This directory contains experimental fastCNV callers and benchmark helpers.
Production `fastCNV.py` remains the primary workflow, but July 5, 2026 work
promoted the pyInferCNV-residual candidate as an opt-in primary-script export:

```bash
python -m altanalyze3.components.fastCNV.fastCNV \
  ...standard fastCNV arguments... \
  --residual-candidate
```

The primary script writes the residual candidate with the prefix
`<output>.residual_candidate.*` so it does not overwrite standard fastCNV
tables.

## July 5, 2026 Update

Problem investigated:

- Standard fastCNV retained the validated unbiased LOY behavior, but did not
  recover the pyInferCNV SamplePre34 positive-control clone structure.
- A whole-genome fastCNV score heatmap looked blank for SamplePre34 because it
  rendered the old fastCNV window-score matrix and then filtered values below a
  z-score threshold. That plot was not showing the residual signal that passed
  the SamplePre34 benchmark.

Resolution implemented:

- Added `pyinfer_residual_clone.py`, an unsupervised residual candidate that:
  - computes pyInferCNV residuals from the query count matrix;
  - scans all autosomes generically, not chromosome-specifically;
  - splits each chromosome within each cell state into neutral/altered
    two-component residual groups;
  - requires minimum component size, component fraction, absolute residual
    shift, and MAD separation;
  - carries validated fastCNV `clone_loy` labels forward rather than replacing
    the unbiased LOY path.
- Added residual-scale inferCNV-style heatmap exports to the candidate:
  - `*.candidate_infercnv_heatmap.png`
  - `*.candidate_infercnv_heatmap.pdf`
  - `*.candidate_infercnv_heatmap_filtered.png`
  - `*.candidate_infercnv_heatmap_filtered.pdf`
- Added default whole-genome fastCNV score heatmap exports to primary
  `fastCNV.py`:
  - `*.infercnv_heatmap.png`
  - `*.infercnv_heatmap.pdf`
  - `*.infercnv_heatmap_filtered.png`
  - `*.infercnv_heatmap_filtered.pdf`
- Added `--residual-candidate` to primary `fastCNV.py` so the July 5 candidate
  can be produced from the main workflow.

## July 5 Benchmark Settings

Residual candidate settings that produced the improved SamplePre34 result:

```text
min_abs_region_score = 0.05
min_separation_mad = 1.5
min_state_cells = 60
min_component_cells = 20
min_component_fraction = 0.03
exclude_chromosomes = chrX,chrY,chrM
heatmap_filter_threshold = 0.03
```

Primary-script flags:

```bash
--residual-candidate
--residual-candidate-min-abs-region-score 0.05
--residual-candidate-min-separation-mad 1.5
--residual-candidate-heatmap-filter-threshold 0.03
```

## July 5 QC Results

SamplePre34 positive-control benchmark directory:

```text
/Users/saljh8/Dropbox/Revio/PyInferCNV/pyInferCNV-inferCNVpy-fastCNV-comp
```

Latest regenerated residual-candidate prefix:

```text
/Users/saljh8/Dropbox/Revio/PyInferCNV/pyInferCNV-inferCNVpy-fastCNV-comp/pyinfer_residual_candidate_fixed_heatmap
```

SamplePre34 results:

- Standard `fastcnv_codex_default`: `3 / 853` chr8 positive-control cells
  recovered.
- fastCNV score-matrix aggregate candidate: `98 / 853` recovered.
- July 5 residual candidate: `852 / 853` recovered.
- Residual candidate total: `3,172` CNV cells and `5,355` WT cells.
- Residual candidate chr8 gain summary: `1,985` cells across `8` states.
- Carried-forward SamplePre34 LOY: `113` `clone_loy` cells.

Corrected SamplePre34 residual-candidate heatmaps:

```text
/Users/saljh8/Dropbox/Revio/fastCNV/samplepre34/fastCNV.residual_candidate_infercnv_heatmap.png
/Users/saljh8/Dropbox/Revio/fastCNV/samplepre34/fastCNV.residual_candidate_infercnv_heatmap.pdf
/Users/saljh8/Dropbox/Revio/fastCNV/samplepre34/fastCNV.residual_candidate_infercnv_heatmap_filtered.png
/Users/saljh8/Dropbox/Revio/fastCNV/samplepre34/fastCNV.residual_candidate_infercnv_heatmap_filtered.pdf
```

LOY positive-control benchmark:

```text
/Users/saljh8/Dropbox/Revio/fastCNV/mds_unbiased_run
```

Validated LOY results retained:

- `14,706` `clone_loy` cells in `fastCNV.cnv_cells.tsv`.
- `14,513` known LOY cells in `classification_accuracy_per_state.csv`.
- Weighted LOY sensitivity: `98.9%`.
- `43 / 46` LOY-bearing states with sensitivity `>= 0.95`.
- `46 / 46` LOY-bearing states with `specificity_vs_supervised >= 0.95`.
- `clone_loy` consensus interval: chrY loss, homozygous-loss, high confidence.

## Candidate Files

- `pyinfer_residual_clone.py`: July 5 residual candidate and residual-scale
  heatmap export. This is the current SamplePre34-positive benchmark path.
- `clone_aggregate.py`: exploratory aggregate caller over standard fastCNV
  score matrices. It improved over the old fastCNV result but did not recover
  the SamplePre34 positive control sufficiently.
- `per_gene.py`: older per-gene candidate retained for reference and negative
  control experiments; not the July 5 positive-control solution.

## Validation Commands

The July 5 code path was validated with:

```bash
python3.11 -m py_compile \
  fastCNV/fastCNV.py \
  fastCNV/candidate/pyinfer_residual_clone.py \
  fastCNV/candidate/clone_aggregate.py

python3.11 -m pytest tests/test_fastcnv_basic.py -q
```

Expected pytest result:

```text
1 passed
```

## Potential false-positive predictions (residual-candidate arm)

The `--residual-candidate` arm is tuned for **sensitivity** (to recover subclones the production
simulation classifier misses). That sensitivity comes with a **recurrent false-positive floor** that
callers should be aware of and filter.

**Observation (13-donor MDS/AML cohort, 2026-07-19).** The residual candidate emits
**`chr18:gain`, `chr13:gain`, and `chr21:gain` in ~10–18% of cells in *every* donor** — including
mutation-free / near-normal marrows. A "CNV" that recurs at the same rate in essentially every sample,
and is not a classic myeloid lesion (−5q / −7 / +8 / −20q / −Y), is almost certainly artifactual.

**Root cause.** For each (cell state × chromosome) the caller fits a two-component (neutral vs altered)
split and labels the minority "altered" component as a gain/loss. With many states × chromosomes this
split will separate a small minority component **even when there is no true copy-number change** — a
calling-sensitivity artifact, not a signal in the underlying residual.

**Evidence it is a calling artifact, not real and not shared** (per-chromosome mean signal, z-scored
across the 13 donors):
- The underlying **pyInferCNV residual mean does not elevate chr18/13/21** (pyInferCNV uses the same
  residual engine but calls clones differently and shows no systematic chr18/13/21 gain).
- fastCNV's own **default simulation-classifier arm suppresses chr18/21** (consistently negative); its
  mild recurrent tendency is chr6/chr11/chr19 instead.
- **inferCNVpy**'s recurrent tendency is chr14/chr16/chr8; chr21 is suppressed.
- The three methods' weak recurrent tendencies **do not agree** (chr18/13/21 vs chr6/11/19 vs
  chr14/16) — exactly what method-specific normalization / gene-density artifacts look like. A *real*
  recurrent CNV would be concordant across methods.

**Options in effect when this occurs** (July 5 benchmark defaults):
```bash
python -m altanalyze3.components.fastCNV.fastCNV ...standard fastCNV args... --residual-candidate
# residual-candidate defaults:
#   --residual-candidate-min-abs-region-score 0.05
#   --residual-candidate-min-separation-mad   1.5
#   min_state_cells 60   min_component_cells 20   min_component_fraction 0.03
#   heatmap_bin_genes 50   heatmap_filter_threshold 0.03
```

**Trustworthy vs suspect calls.** Trustworthy: **LOY** (`chrY:loss`, sex-gated to predicted-male
cells), del5q/+8, complex karyotypes, and donor-specific events that rise clearly above the recurrent
floor. Suspect: the pan-donor `chr18/13/21:gain` layer — treat as artifact unless independently
confirmed.

**Mitigation.**
1. For a **specificity-first CNV landscape**, use the production simulation-classifier arm (calls
   ~95–100% of cells WT and flags only high-confidence events), and reserve `--residual-candidate` for
   subclone recovery.
2. **Require cross-method concordance** — a real CNV should appear in ≥2 of {simulation classifier,
   residual candidate, inferCNVpy, pyInferCNV}.
3. **Tighten the residual-candidate thresholds** to suppress the minority-component false calls:
   raise `--residual-candidate-min-abs-region-score` (e.g. 0.05→0.10),
   `--residual-candidate-min-separation-mad` (1.5→2.0+), and/or `min_component_fraction` (0.03→higher).
