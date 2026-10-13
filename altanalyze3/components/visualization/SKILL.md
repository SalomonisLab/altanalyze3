---
name: run-markerfinder
description: Run MarkerFinder marker discovery through the established AltAnalyze3 cellHarmony/scALABLE hook, producing its standard marker tables, heatmap, caches and companion outputs. Use for requests to run or rerun MarkerFinder on RNA or other expression modalities.
---

# Run standard MarkerFinder

## Required route

Use the existing public hook in `altanalyze3.components.visualization.marker_heatmap_h5ad`:

- API: `generate_marker_heatmap_from_adata(..., marker_method="markerfinder")`.
- CLI: `python -m altanalyze3.components.visualization.marker_heatmap_h5ad --marker-method markerfinder ...`.

This is the hook invoked by `components/cellHarmony/flask/pipeline.py`; scALABLE consumes its standard outputs in `components/visualization/scalable_viewer/prepare_assets.py`. The local checkout is `/Users/saljh8/Documents/GitHub/altanalyze3`, with code under `altanalyze3/components`. Check the current hook and harness arguments when using another checkout or version.

Never write a replacement MarkerFinder runner, marker-selection algorithm, correlation/FDR calculation, heatmap renderer, or output formatter. Do not bypass the public hook by directly calling `components.cellHarmony.markerFinder.marker_finder`, its private selection functions, or its standalone alternate CLI. Even when correlations agree, that bypass produces different selection criteria and omits standard artifacts. Do not substitute scanpy ranking. Orchestrating existing commands and using the public API are permitted; bespoke analytical wrappers are not.

## Filtered AnnData objects: preferred invocation

Pass the filtered in-memory AnnData object directly to the public hook. It already supports this; a file is not required. Create the output directory before the API call (the CLI does this automatically). Preserve X, requested layers, obs, var and `uns['lineage_order']` when subsetting. Ordinary AnnData filtering and a direct hook invocation are permitted input preparation, not a replacement analytical workflow.

```python
from altanalyze3.components.visualization.marker_heatmap_h5ad import generate_marker_heatmap_from_adata
from pathlib import Path

filtered = adata[verified_cell_mask].copy()
Path("/absolute/path/MarkerFinder").mkdir(parents=True, exist_ok=True)
outputs = generate_marker_heatmap_from_adata(
    filtered,
    cluster_key="Mm-MarrowAtlas-L4",
    out="/absolute/path/MarkerFinder/cell_state_marker_heatmap.pdf",
    marker_method="markerfinder",
    top_n=50, cells_per_cluster=100, seed=0,
    export_networks=True, network_top_n=1000, network_jobs=4,
    species="mouse",
    write_heatmap_tsv=False, write_expression_tsv=False,
    write_heatmap_cache=True, render_heatmap=True, write_svg=True,
)
```

Use supported API arguments for requested changes. Do not serialize filtered objects solely because the CLI uses files. If input preparation or a supported argument is insufficient, diagnose the existing hook; do not fall back to custom MarkerFinder code. Running outside this standard workflow is an unacceptable failure mode.

## File inputs and invocation

Inspect the source modalities, `obs` grouping labels, expression scale and `uns['lineage_order']`. Use full-precision source expression rather than rounded exports. Preserve the original lineage order (including qHSC then aHSC when present). Resolve known tissue/library inconsistencies before selecting cells, without altering the original H5AD.

For subsets, use the existing `altanalyze3.components.aggregate.h5ad_subset` CLI/API. It preserves expression, metadata and lineage order. Select verified library values if a tissue column is wrong; do not invent preprocessing scripts. `--drop-layers` is appropriate only when the requested analysis uses normalized X and no raw/counts-dependent output is needed.

The standard RNA web-harness invocation is:

```bash
PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 \
/opt/homebrew/opt/python@3.11/bin/python3.11 \
  -m altanalyze3.components.visualization.marker_heatmap_h5ad \
  --h5ad /absolute/path/input.h5ad \
  --cluster-key 'Mm-MarrowAtlas-L4' \
  --marker-method markerfinder \
  --top-n 50 --cells-per-cluster 100 --seed 0 \
  --export-networks --network-top-n 1000 --network-jobs 4 \
  --species mouse --skip-expression-tsv \
  --out /absolute/path/MarkerFinder/cell_state_marker_heatmap.pdf
```

Use the requested cluster key, species and destination; do not assume mouse or marrow labels for other projects. `--top-n` controls the displayed heatmap subset, not the complete assigned-marker table. `--cells-per-cluster` samples the display only; discovery uses all input cells. Keep the hook's selection/statistical defaults unless the user requests supported changes. Do not carry over the alternate low-level runner's rho=0.2 threshold.

RNA scaling validation stays enabled. If raw counts need normalization, use the hook's supported `--scale-data` option. For signed or imputed modalities such as DSB ADT, inspect the harness's modality route: it uses `--no-scaling-check --centroid-method mean` and typically `--top-n 5`. Do not apply this RNA normalization bypass to ordinary RNA. Do not modify or reimplement the hook to handle a failed input.

`--skip-expression-tsv` matches the web harness and omits its large expression/fold TSVs while retaining markers, centroids, cache and static heatmaps. Omit that flag when those matrix TSVs are requested. Do not pass `--markers-only`, `--skip-heatmap-render`, `--skip-svg` or `--skip-heatmap-cache` for a normal MarkerFinder request: the standard well-established heatmap and companion artifacts are required. Network exports use the existing NetPerspective route; report missing dependencies rather than fabricating substitutes.

## Splicing events (PSI) and other inputs with missing values

The same hook runs PSI. NaN is a missing value, never a measurement. MarkerFinder scores each event
over the samples that measured it, with that event's own n - 2 degrees of freedom, as AltAnalyze2
`markerFinder.advancedPearsonCorrelationAnalysis` did. Means (`Query Exp`, `Ref Exp`, centroids)
divide by the measured samples. Input without NaN takes the unchanged code path.

The PSI options are PSI-specific. None changes gene-expression output.
- `--max-missing-fraction` acts only on a matrix that holds NaN.
- `--fdr-threshold` defaults to the API's `pval_threshold` of 0.001.
- `--max-sample-missing-above-group-median`, `--label-events` and `--collapse-shared-exons` are off unless given. `--psi` turns on the first and third.
- The `SYMBOL:junction|junction` label applies only to splicing-event IDs.

Validated 2026-10-11 on the scVelo pancreas data (3,696 cells, raw counts) against commit 3e979dd:
- **Default scanpy route:** every TSV is byte-identical, and the PDF text and images are identical.
- **RNA MarkerFinder route** (`--marker-method markerfinder --scale-data --top-n 50 --cells-per-cluster 100`): the same holds. The heatmap cache h5ad has identical contents.

Log: `logs/validate_gene_defaults_unchanged.log` in the approved PSI run folder.

PSI invocation: `--psi` applies the PSI defaults (Nathan, 2026-10-11: "appropriate for PSI analysis in
general"). They are the settings of the KINNEX-5 healthy run he approved:
`/Volumes/salomonis-archive/LabFiles/Nathan/Revio/MDS-AML-KINNEX-5/analysis/PSI/healthy8_markerfinder_min5_depthz3_state70_pbmiss10_present45_top100_rho025_fdr01_noMultiLinGMP3_MEP1_20261011/MarkerFinder/cell_state_marker_heatmap.pdf`.
A flag given explicitly overrides its preset value. API: `generate_marker_heatmap_from_adata(..., **PSI_DEFAULTS)`.
Pass `--top-n` explicitly so the heatmap cut is a stated choice (100 in the approved run; 50 also made).

```bash
python -m altanalyze3.components.visualization.marker_heatmap_h5ad --h5ad <psi.h5ad> \
  --cluster-key <state column> --psi --top-n 100 --seed 0 --species human \
  --out <dir>/cell_state_marker_heatmap.pdf
```

| Option | Hook default | `--psi` value | Effect |
|---|---|---|---|
| `--marker-method` | scanpy | markerfinder | MarkerFinder statistics |
| `--no-scaling-check` | off | on | PSI is not a count matrix |
| `--centroid-method` | log2_cp10k | mean | Centroids are mean PSI over measured samples |
| `--max-sample-missing-above-group-median` | off | 0.10 | Remove a sample (pseudobulk) whose NaN fraction exceeds its state's median by more than 10 points: a low-coverage outlier. Every decision goes to `*_sample_missing_filter.tsv` |
| `--max-missing-fraction` | 0.15 | 0.55 | An event missing in this fraction of samples or more is not considered |
| `--min-rho` | none (rho > 0) | 0.25 | Minimum MarkerFinder score for assignment |
| `--fdr-threshold` | 0.001 | 0.01 | Maximum Benjamini-Hochberg FDR (per population) for assignment; the API's `pval_threshold` |
| `--collapse-shared-exons` | off | on | Of the assigned events that share an exon block, keep the most significant (FDR, then rho); the rest go to `*_shared_exon_excluded.tsv` |
| `--top-n` | 5 | 100 | Events per state on the HEATMAP only; the marker table keeps every assigned event |

Upstream of the hook, the approved run prepared its input with analysis scripts, not the hook: healthy donors only;
pseudobulks with |robust z| > 3 of log10 junction reads within their state removed; states with >= 5 pseudobulks and
<= 70% missing PSI kept; then one HSC, one LMPP and one MEP state, and other states with >= 3 events
(`scripts/` in that run folder).

`--psi` validated 2026-10-11: on the unfiltered 318-pseudobulk input restricted to the approved run's 32 states, it
reproduced that run's marker table byte for byte (3,567 events), the same 199 retained pseudobulks and the same 2,955
shared-exon exclusions. Without `--psi`, output is unchanged (marker table byte-identical to a prior run). Log:
`logs/validate_psi_preset.log` in the approved run folder.

Keep the hook's rasterized heatmap PDF for PSI. Nathan rejected a `--no-rasterize` (one pixel per cell) render on
2026-10-11.

Row labels for splicing events (Nathan, 2026-10-11). The heatmap prints each event as
`SYMBOL:junction|junction`, for example `ATP5J:E1.7-E3.1|E1.11-E3.1`. The gene ID becomes its symbol
from the bundled `long_read/resources/<Hs|Mm>_Ensembl-annotations.txt.gz` (`--species`). The second
junction drops the repeated gene. Other labels, such as gene symbols, are unchanged. Tables keep the full UIDs.

`--label-events <file>` labels only the listed events and their related rows, as callouts with leader
lines to the right of the heatmap. It changes no heatmap row, column or value.
- **File:** text or TSV. The hook reads its `UID` column, else its first column.
- **Which rows get a label:** a listed event in either junction order (A|B or B|A). A row that shares an
  exon block with a listed event is labelled too; that is the shared-exon collapse's redundancy rule
  (Nathan, 2026-10-11). `--label-events-exact` turns the shared-exon rows off.
- **Label text:** each label shows the heatmap row's own junction order.
- **Record:** `*_label_events.tsv` lists every labelled row, its state, its listed event and the shared exons.
- **Fit:** the callout font shrinks so the longest label fits on the page; the page and heatmap keep their size.
- **Where it works:** on a full run and with `--render-from-cache`.

The approved run's 16-event list is `MarkerFinder_label_whitelist/label_events.tsv` in that run folder.
Its top-100 heatmap labels 14 of its 2,160 rows: 3 listed events (ERICH1, ERP29 and MBNL1, all in
reciprocal order) and 11 rows that share an exon with a listed event.

```bash
python -m altanalyze3.components.visualization.marker_heatmap_h5ad \
  --render-from-cache <run>/MarkerFinder/cell_state_marker_heatmap_fold_matrix.h5ad \
  --out <run>/MarkerFinder_label_whitelist/cell_state_marker_heatmap.pdf \
  --label-events <events.tsv> --species human
```

The cache fixes the heatmap rows, so `--render-from-cache` takes no `--top-n`.

Never fill NaN before the hook: a PSI of 0 means full exclusion. Never pass `--matrix-tsv` PSI:
that loader still fills NaN with 0 (`_load_text_expression_inputs`). Build an h5ad with NaN kept.

Validated 2026-10-11 against the AltAnalyze2 healthy PSI run
(`/Users/saljh8/Dropbox/Revio/BoneMarrow/BAMs/AltAnalyze-Healthy/Restricted`): rho for 13,382 events x 19
states and all 13,382 p-values reproduced (max |diff| 5.5e-15; relative p 9.8e-13). Scripts and logs:
`/Volumes/salomonis-archive/LabFiles/Nathan/Revio/MDS-AML-KINNEX-5/analysis/PSI/healthy8_markerfinder_test_20261010/`.

## Completion checks

Confirm success in the hook's log and inspect the standard outputs:

- `*_markers.tsv`: complete unique assigned-marker statistics, not a manually capped list.
- `*_redundant_markers.tsv`: population-specific companion statistics.
- Standard heatmap PDF and SVG from the hook's renderer.
- `*_fold_matrix.h5ad` cache (older versions may use NPZ).
- `*_fold_matrix.centroids.tsv`, logs, and networks when enabled.
- Expression/fold TSVs when requested through the supported switches.

Check input counts, tissue scope, marker/state coverage and reference ordering using existing tools. Inspect the rendered standard PDF for legibility. If redraw is needed, use this module's `--render-from-cache`, never a custom plotting script. scALABLE ingests these artifacts; do not rebuild a heatmap for it.

Preserve failed or superseded runs separately and identify the authoritative standard-hook result. Report actual output paths and relevant limitations: display sampling versus full-data scoring, sparse populations, RNA-derived labels, and cell-level statistics versus biological replication. Do not claim marker agreement independently validates labels derived from the same RNA data.
