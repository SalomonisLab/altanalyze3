# Variant-Impact: Marker-Based Imputation of Mutation Status

Given the read-level genotypes from **[VARIANT_DISCOVERY.md](VARIANT_DISCOVERY.md)** (a minority
of cells whose mutation status is directly observed), this workflow learns a **gene/isoform
expression signature that separates MUT from WT cells within each cell state**, then uses that
supervised classifier to **impute the genotype of the many undetected (0-read) cells** and expand
the mutation map. It is the `altanalyze3 variant-impact` subcommand.

Code: `altanalyze3/components/bam/variant_impact.py` (single module; orchestration +
`genotype_from_bam` + `impute_unknown`; markers/heatmaps come from the benchmarked
`cellHarmony` MarkerFinder). CLI wiring: `altanalyze3/utilities/parser.py` (subcommand
`variant-impact`, `func=run_variant_impact`).

Production run dir (all numbers below are read from it):
`/Volumes/salomonis-archive/LabFiles/Nathan/Revio/MDS-AML-KINNEX-4/variant_impact_RUNX1_prototype/`
(cluster: `/data/salomonis-archive/...`).

---

## 1. What the pipeline does (per sample × variant × level)

1. **Genotype** each cell from the BAM `=`/`X` CIGAR (`genotype_from_bam`, no FASTA).
2. For each `--level` (`gene` | `isoform` | `both`) and **each cell state** with ≥ `--min-cells`
   confident MUT **and** WT cells: run the production **MarkerFinder**
   (`cellHarmony.generate_marker_heatmap_from_adata`, `marker_method="markerfinder"`) MUT-vs-WT,
   take the top `--top-n` markers per group as the **signature**, render a labeled MUT/WT heatmap PDF.
3. **Impute** the undetected (`UNK`) cells of that state from the signature (`impute_unknown`) and
   **expand** the genotype table.
4. A **combined** (all states pooled) MarkerFinder run on the confident calls.
5. A cross-cell-state **marker consistency** table (`marker_consistency.tsv`).

A per-state expansion is written to the final table **only if its classifier re-predicts its own
training calls at ≥ 85%** (`ACCURACY_MIN = 0.85`).

---

## 2. The classifier (`impute_unknown`) — interface and method

```python
impute_unknown(expr_all, geno_called, geno_unknown, signature, min_called=10)
    -> DataFrame[barcode, imputed(MUT/WT/UNK), score, threshold]
       .attrs: concordance, mut_recovery, wt_recovery, n_called_mut, n_called_wt
```

**Method = median-centered nearest-centroid (MUT vs WT):**
- `mut_c`, `wt_c` = per-feature **median** expression of the called MUT / WT cells over the
  signature features (normalized CP10k+log1p).
- `ref` = median-of-the-two-centroids (validated "method B"); center every cell by `ref`.
- Score = `‖v − wt_c‖ − ‖v − mut_c‖` (Euclidean); `>0` → nearer MUT.
- **Threshold** = 10th percentile of the correctly-scored *called* cells' margins; an `UNK` cell is
  imputed MUT if `score ≥ thr`, WT if `score ≤ −thr`, else left `UNK`.

**Built-in verification (resubstitution).** The classifier re-predicts the original called cells it
was built from; `concordance` = fraction whose `sign(score)` matches the read genotype, with
per-class `mut_recovery` / `wt_recovery`. Written to `<tag>_classifier_verification.txt`
(`status = PASS` / `FAIL (<85%; not deployed)`).

> **This is a resubstitution (train == test) metric, by design** — the code comment states the
> centroid classifier makes train==test acceptable. It is an **optimistic ceiling**, not held-out
> accuracy. See §5 for how the independent DEG-concordance test qualifies it.

---

## 3. CLI interface (verified against `parser.py`)

```
altanalyze3 variant-impact \
  --metadata <sclr metadata.tsv>     # required; columns uid, bam, library (+ reverse, groups)
  --variants <chr:pos:label;... | mutations file>   # required; genotyped from =/X CIGAR (no FASTA)
  --level    gene|isoform|both       # default gene
  --cell_annot <barcode\tcluster>    # else uses obs['cluster'] in the h5ad
  --mut-min  1                       # min X reads to call MUT (combined confident set)
  --wt-min   2                       # min = reads to call WT (combined confident set)
  --min-cells 20                     # min MUT and WT per state to run MarkerFinder
  --top-n    50                      # top markers per group
  --min-mapq 20
  --indel-window 5
  --impute                           # impute UNK cells from the per-state signature
  --gene_symbol <ensembl->symbol>    # optional gene-level heatmap labels
  --samples  <uid ...>               # optional subset
  --output   <dir>                   # default ./variant_impact
```

Per-sample matrices are read from `<bam_dir>/<library>-<level>.h5ad` (i.e. beside the BAM). BAM
paths and libraries come from the `--metadata` TSV.

---

## 4. Exact production commands and inputs (sample `5801M_pre`)

Both runs call the **unchanged** `run_variant_impact` (repo injected via
`sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")`), interpreter
`/opt/homebrew/opt/python@3.11/bin/python3.11`.

**RUNX1** — `run_runx1_wt1.py`:
```python
args = SimpleNamespace(
    metadata   = "metadata_local.txt",
    output     = "variant_impact_RUNX1_wt1",
    level      = "both",
    variants   = "RUNX1_mutation.txt",         # chr21 34799432 34799432 RUNX1p.W279* SNV
    cell_annot = cell_annot_5801M_pre.txt,     # built from the isoform h5ad obs['cluster']
    mut_min=1, wt_min=1, min_cells=15, top_n=50,
    min_mapq=20, indel_window=5,
    gene_symbol=None, samples=["5801M_pre"], impute=True)
run_variant_impact(args)
```

**SRSF2** — `run_srsf2_wt1.py`: identical except `variants=SRSF2_mutation.txt`,
`output=variant_impact_SRSF2_wt1`, `min_cells=20`.

Equivalent CLI form:
```bash
altanalyze3 variant-impact --metadata metadata_local.txt --variants RUNX1_mutation.txt \
  --level both --cell_annot cell_annot_5801M_pre.txt \
  --mut-min 1 --wt-min 1 --min-cells 15 --top-n 50 --impute \
  --samples 5801M_pre --output variant_impact_RUNX1_wt1
```

Inputs:
- `metadata_local.txt` — cols `uid  bam  library  reverse  groups`; `5801M_pre` →
  `.../SamplePre34_outputs/BAM-clusters/scisoseq.mapped.bam`, library `WM71-0121-5801_CD34`.
- `cell_annot_5801M_pre.txt` — 7,978 rows `barcode<TAB>cell_state` (e.g. `...AA-1  MEP-2`), built
  from `WM71-0121-5801_CD34-isoform.h5ad` `obs['cluster']` (the gene h5ad has no `cluster`).

**Why `wt_min=1`:** at `wt_min=2`, long-read ~1-read/cell coverage collapsed RUNX1 WT to ~20 cells
total → **no state testable**. `RUNX1_testability_by_state.tsv` shows every state `testable_wtmin2=False`;
only MPP-1 (26/51) and MPP-MEP (21/44) become testable at `wt_min=1`. Trade-off: `≥1`-read WT admits
allelic-dropout false-WT.

---

## 5. Performance (read from the run's summary tables)

### 5a. Per-cell-state classifier accuracy — RUNX1 W279* (`variant_impact_RUNX1_wt1/variant_impact_summary.tsv`)

| level | cell state | n_MUT | n_WT | accuracy % | pass ≥85 | imputed MUT | imputed WT | deployed |
|---|---|---|---|---|---|---|---|---|
| gene | **combined (pooled)** | 256 | 435 | **59.8** | ✗ | 3386 | 3233 | ✗ |
| gene | LMPP-1-cycling | 15 | 19 | 100.0 | ✓ | 13 | 55 | ✓ |
| gene | MEP-2 | 15 | 30 | 97.8 | ✓ | 39 | 174 | ✓ |
| gene | MPP-MEP | 21 | 44 | 90.8 | ✓ | 199 | 350 | ✓ |
| gene | MPP-1 | 26 | 51 | 76.6 | ✗ | 390 | 355 | ✗ |
| isoform | **combined (pooled)** | 256 | 435 | **68.7** | ✗ | 1935 | 4338 | ✗ |
| isoform | LMPP-1-cycling | 15 | 19 | 97.1 | ✓ | 18 | 116 | ✓ |
| isoform | MEP-2 | 15 | 30 | 97.8 | ✓ | 44 | 245 | ✓ |
| isoform | MPP-MEP | 21 | 44 | 90.8 | ✓ | 168 | 401 | ✓ |
| isoform | MPP-1 | 26 | 51 | 81.8 | ✗ | 198 | 620 | ✗ |

**Central finding: the pooled/combined classifier FAILS (gene 59.8%, isoform 68.7%); per-cell-state
classifiers PASS (90.8–100%).** Mixing states buries the within-state MUT/WT signal under the much
larger between-state expression variance — **per-state training is required**. Only states passing
the 85% gate contribute imputed cells to `sample_variant_calls.tsv`.

### 5b. SRSF2 P95R (`variant_impact_SRSF2_wt1` / `_wt2/variant_impact_summary.tsv`)

Pooled again fails (gene 58.9% / 59.3% at wt1; 60.4% / 55.7% at wt2). Many **per-state** classifiers
pass at 85–100% (e.g. gene MultiLin-GMP-2 100.0%, MEP-Eryth-2 96.2%, MEP-1 92.8%; isoform MEP-Eryth-2
94.9%, MultiLin-GMP-2 93.0%, ERP-5 92.8%). At wt2 a RUNX1 combined row appears at gene 86.2%
(256 MUT / **20** WT) — the WT-collapse that motivated the wt1 re-run.

### 5c. Independent validation — DEG concordance (does the imputation recover real biology?)

`deg_concordance.py` runs limma (`oncosplice.metadataAnalysis.limmaCompute`) MUT-vs-WT on the
**original called** cells vs the **newly imputed** cells (mutually exclusive) and correlates
per-feature log-fold-changes. Summary (`imputation_concordance_summary.tsv` /
`deg_concordance/deg_concordance_summary.tsv`):

| variant | level | state | orig cells | imputed gained | expansion× | concordant frac | Pearson (DEGs-in-both) |
|---|---|---|---|---|---|---|---|
| RUNX1 wt1 | isoform | MPP-MEP | 65 | 569 | 8.75 | 0.91 | 0.88 |
| RUNX1 wt1 | isoform | MEP-2 | 45 | 289 | 6.42 | 0.92 | 0.84 |
| RUNX1 wt1 | gene | MEP-2 | 45 | 213 | 4.73 | 0.86 | 0.87 |
| RUNX1 wt1 | gene | MPP-MEP | 65 | 549 | 8.45 | 0.78 | 0.75 |
| SRSF2 wt2 | gene | MEP-2 | 135 | 269 | 1.99 | 0.82 | 0.85 |
| SRSF2 wt1 | gene | MultiLin-GMP-1 | 77 | 73 | 0.95 | 0.76 | 0.64 |

**RUNX1 (a transcription factor) is the strong positive**: 4/5 states pass the rho ≥ 0.6 bar on
co-detected DEGs (0.75–0.88), directional concordance 76–92%.

**SRSF2 P95R (a splicing factor) is negative/inapplicable at the expression level.** Documented in
`variant_impact_SRSF2_wt2/deg_concordance/README.md`: *"the imputation is neither confirmed nor
refuted by this test, because there is no within-state DE signal to be concordant about."* Three
diagnostics support **absence-of-signal, not imputation failure**:
- Null p-value distribution (MEP-2 gene: 231 observed vs 223 expected DEGs at p<0.05; at p<0.01,
  38 observed vs 45 expected — *below* chance).
- **Split-half reproducibility ceiling ≈ 0.00** on the original within-genotype contrast
  (−0.02 … +0.03).
- **Positive control** — the *between-cell-state* contrast reaches split-half rho **0.55–0.61**
  (`positive_control_between_state.tsv`: MEP-2 vs LMPP-1-cycling 0.612, vs ERP-1 0.563), proving the
  pipeline detects real DE.

Biology: within one cell state, SRSF2-MUT and -WT cells are near-identical in gene/isoform
*expression magnitude*; the mutation acts on **splicing (junction/isoform-ratio usage)**, which
per-feature expression DE does not capture. **Caveat:** even the real between-state control caps at
~0.6, so single-cell sparsity puts the 0.6 pass-bar essentially at the ceiling — the DEG test is
most useful as a **negative screen**.

> **FDR caveat (important for honesty).** In `deg_concordance_fdr/`, once DEGs are FDR-gated,
> `n_DEG_orig = 0` in every RUNX1 state, so the strong "DEGs-in-both" Pearson in the non-FDR file is
> carried by genes selected in the *imputed* set. The within-state genotype contrast has ~zero
> reproducible signal by the split-half ceiling even for RUNX1. So: **the ≥85% resubstitution
> accuracy is the primary support for the imputed calls; the DEG test confirms RUNX1 and is
> inapplicable to SRSF2 — it does not independently prove either mutation map.**

---

## 6. Output files

Under `<output>/`:
- **`sample_variant_calls.tsv`** — final per-cell map:
  `sample  variant  level  barcode  cell_state  genotype(MUT/WT)  source(called/imputed)`.
  Contains every read-called cell + imputed cells from states that passed 85%.
- **`variant_impact_summary.tsv`** — the per-unit table in §5a
  (`unit, n_MUT, n_WT, pct_correct_MUT, pct_correct_WT, accuracy_pct, pass_85, n_imputed_MUT,
  n_imputed_WT, deployed`).
- `<uid>/<label>/<level>/combined/` and `.../per_celltype/<state>/`:
  `<tag>_marker_heatmap.pdf`, `<tag>_marker_heatmap_markers.tsv` (the signature),
  `<tag>_classifier_verification.txt`, `<tag>_imputed.csv` (`barcode,imputed,score,threshold`).
- `marker_consistency.tsv` — features that recur as MUT/WT markers across states.

---

## 7. Limitations (summary)

1. **Deployment accuracy is resubstitution (train==test)**, an optimistic ceiling — not held-out CV.
2. **Pooled classifier does not work** — you must train per cell state (see §5a); low-cell states
   (e.g. MPP-1) fall below the 85% gate and are not deployed.
3. **WT scarcity** for clonal mutations at ~1 read/cell forces `wt_min=1`, which admits
   allelic-dropout false-WT contaminating the training labels.
4. **Expression-DEG validation is blind to splicing-factor variants** (SRSF2): it neither confirms
   nor refutes them. A junction-PSI / isoform-ratio (oncosplice) readout, or a TF variant with
   ≥20 cells/state, would be needed to validate those.
5. **Barcode orientation** (isoform BAM = reverse-complement of junction/short-read barcodes) is
   handled by `_pick_orientation`; a 0-cell overlap means the orientation guess or the `cell_annot`
   barcodes are wrong.
6. **Reproducibility gaps in the run dir**: no `run_parameters.json`, no pinned altanalyze3 commit,
   reference genome is a symlink. Parameters live only inline in `run_runx1_wt1.py` /
   `run_srsf2_wt1.py`. Pin these before publishing any number.

---

## 8. Reproduce

```bash
PY=/opt/homebrew/opt/python@3.11/bin/python3.11
cd /Volumes/salomonis-archive/LabFiles/Nathan/Revio/MDS-AML-KINNEX-4/variant_impact_RUNX1_prototype
$PY run_runx1_wt1.py          # or run_srsf2_wt1.py; or the `altanalyze3 variant-impact` CLI in §4
$PY deg_concordance.py        # independent DEG-concordance validation (no BAM access needed)
```

*Cross-reference: memory `variant-impact-deg-concordance`, `rigorous-validation`,
`imputation-extrapolate-not-retrieve`; upstream genotyping →
[VARIANT_DISCOVERY.md](VARIANT_DISCOVERY.md).*
