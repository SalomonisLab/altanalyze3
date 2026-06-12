# pyudon ↔ altanalyze3/udon integration

Reconciliation of the original developer package (`pyudon`, "former student") with the current
altanalyze3 UDON code. Goal: best of both worlds — the original's validated SATAY-UDON API + the
current code's matched-control normalization, study-aware integration, gene filter, GO-Elite, and
richer figures — as a **reproducible, command-driven** workflow.

## Conventions chosen (pyudon vs current)
| topic | decision |
|---|---|
| SATAY-UDON API | **pyudon API is canonical** (`make_mean_based_binary`, `fishers_clinical_feats`, `cmh_clinical_feats`, `fdr_correction`, `satay_udon`) — see `satay_udon_core.py`, importable from `satay_udon`. |
| SATAY donor floor | **current**: collapse to ≥3 distinct **donors** (avoids pseudoreplication). pyudon used ≥3 pseudobulks. `donor_of=None` reproduces pyudon's per-pseudobulk counting. |
| Numeric covariates | **strict** by default (error/pre-bin); `--mean-binarize cols` opt-in (pyudon's `make_mean_based_binary`). |
| NMF `n_run` | **default 5** (pyudon stability); `n_run=1` opt-in (validated deterministic NNDSVD). |
| Enrichment | **GO-Elite** (`goelite_enrichment.py`) over pyudon's PathwayCommons (both retained). |
| Gene filter | **current**: protein-coding only + drop RPL/RPS + XIST/TSIX, keep Y (loss-of-Y). |
| Figures | **current**: donor/celltype composition + final_program + per-study + `satay_byfield` + `satay_heatmap` (Adobe-safe, embedded fonts). |

## Adopted FROM the developer (restorations / ports)
- `satay_udon_core.py` — pyudon API; **Fisher numerically identical to pyudon** (`validate_satay_core.py`, max|Δ|=0).
- `cmh_clinical_feats` — batch-stratified Cochran-Mantel-Haenszel; **ported and FIXED** (pyudon indexed an int array with a string batch and had an operator-precedence bug). Drives the `--batch-key` step.
- `make_mean_based_binary` — verbatim (numeric→binary at the column mean).
- `classificationFunctions.py` — restored `gene_pearson_corr`; restored the `min(shape)` PCA `n_comps` bound.
- `nmf.run_nmf` — `n_run` parameter, default 5 (the developer's stability default).
- `determine_controls.py` — ported `determine_cell_type_deg` (per-cell-type Wilcoxon DEGs) and `check_similarity_distribution` (adaptive HVG threshold) from `sashimi_udon.py`.
- `run_udon_pseudobulk.py` — `--min-cells` QC (drops low-cell pseudobulks; the developer's QC, adapted to the pre-made-pseudobulk path).

## Standard metadata input (the key fix)
The original hardcoded a 2-tab AML xlsx. SATAY-UDON now takes a **standard tab-delimited covariate
table** via `--metadata` (`SATAY_METADATA_FORMAT.md`): `Sample` required; `Donor_ID`/`Study` optional;
other columns binary (0/1) or categorical. Dataset-specific cleaning lives **with the data**, not the
software — e.g. `UDON/prepare_satay_metadata.py` converts the AML xlsx → standard tsv.

## Module / package completeness (every pyudon module now present)
The udon directory is now an importable **package and a strict superset of pyudon** (`__init__.py`
exposes the full pyudon API + the altanalyze3 additions; verified all 18 pyudon API functions import).
Missing modules were copied in and their interfaces wired:
- `sashimi_udon.py` — copied (control-cell determination: `determine_control_cells`,
  `determine_cell_type_hvg/deg`, `determine_control_cells_cosine_similarity/one_class_svm`,
  `check_similarity_distribution`); `determine_controls.py` is now a thin alias to it.
- `pre_process.py` — replaced my reduced copy with the **full** pyudon version (restores
  `pseudobulk_wrapper`, `matched_controls`, `collective_controls`, and the pseudobulk QC).
- `datasets/` (`loaders.py` + `__init__.py`) — pyudon demo-data fetchers (`fetch_geo_raw`,
  `fetch_processed_h5ad`); optional subpackage (needs `pooch`), not imported by the main package.
- `__init__.py` — pyudon-compatible package API (so `run_example_dataset.py`-style code works) plus
  `pseudobulk_protocol`, `goelite_enrichment`, `markerFinder`.
Standalone CLI execution is unaffected (scripts still run directly; the package uses a sys.path shim).

## Canonical / legacy
- **Canonical:** `satay_udon_core.py` (engine + CLI), `run_workflow.py` (master, logged+timed).
- **Legacy (retained):** `satay_metadata_enrichment.py` (multi-run cross-comparison orchestration).
- **Redundancy noted:** `fast_feature_selection.py` is a validated `--fast` opt-in (not the default path); `visualizations.plot_metadata_heatmap` is superseded by `udon_binary_heatmaps.py` for donor/celltype.

## Validation
- `validate_satay_core.py` — engine vs original pyudon (Fisher identical; CMH runs; metadata loader strict).
- `validate_workflow_outputs.py <outdir>` — files present, PDFs Adobe-safe + fonts embedded, enrichments non-trivial, logs/timeline complete.

## Reproduce
```
# dataset prep (in the output/data dir, not the software):
python UDON/prepare_satay_metadata.py            # AML xlsx -> satay_metadata.tsv

# full analysis into a fresh, tracked directory with logs + timing:
python run_workflow.py --outdir RUN_DIR \
  --pseudobulk PB.h5ad --counts COUNTS.h5ad --sample-metadata SAMPLE.tsv \
  --metadata satay_metadata.tsv --mean-binarize age,blast_pct --batch-key Study --fast
python validate_workflow_outputs.py RUN_DIR
```
