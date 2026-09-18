# Marker output benchmark — 2026-09-18

Optional output controls select which marker files a run writes. On 2026-09-18 the user made **Skip static heatmap files** the scALABLE default, so a default run now takes the second column below. Marker discovery still uses all cells. Interactive heatmaps and Download PDF remain available. The `marker_heatmap_h5ad` package API and CLI keep their own defaults; only scALABLE changed.

The tables below label the pre-2026-09-18 web setting "Default". That column is now the opt-in **Generate static PDF/SVG** setting.

## Study and method

Loaded job: `1b0faa928a9848fbb3cd225ca26a0e42`, marrow study, 147,376 cells and 86 cell states. Saved expression matrices contain 49,056 RNA genes, 129 ADTs, 2,533 metabolites, 1,009 lipids and 217 TF activities.

Each case reran MarkerFinder and its output generation on the same saved matrix in a separate output folder. Input loading, imputation, alignment, downstream differential analysis and network export are excluded from the reported marker-stage times. No source H5AD or saved analysis was modified. Network export was disabled equally in every benchmark case; its web default is unchanged.

Baseline matches the pre-2026-09-18 web marker output defaults: 100 displayed cells per state, 2400-DPI PDF plus companion SVG, compact H5AD heatmap cache, and no expression/fold-matrix TSVs. RNA selects 50 unique markers per state; imputed modalities select 5. Both unique and redundant marker tables and state centroids are always retained.

Two passes tested the three principal modes; the confirmation pass reversed their order (no static files, markers only, default). The first pass also tested three RNA rendering alternatives. These are local wall-clock measurements, not a fresh end-to-end imputation run. The first pass overlapped some regression testing; the confirmation pass ran without that test workload.

## Confirmation pass: seconds per modality

| Modality | Default | Skip static PDF/SVG; keep interactive cache | Markers only; no plotting matrix/cache |
|---|---:|---:|---:|
| RNA | 58.50 | 18.18 | 13.47 |
| ADT | 28.97 | 1.19 | 0.71 |
| Metabolites | 40.20 | 12.24 | 11.74 |
| Lipids | 32.49 | 4.96 | 4.53 |
| TF activity | 27.94 | 1.48 | 1.03 |
| **Total** | **188.10** | **38.05** | **31.48** |

Skipping static figures cut these five marker stages by **79.8%**, a **4.94×** speedup in the confirmation pass.

| Mode | First pass (s) | Confirmation (s) | Generated artifacts (MB, first pass) |
|---|---:|---:|---:|
| Default | 216.86 | 188.10 | 123.46 |
| Skip static figures | 36.28 | 38.05 | 47.96 |
| Markers only | 33.61 | 31.48 | 7.74 |

Artifact sizes exclude logs and the primary expression H5ADs. The latter are retained for Explore, Chat, differentials and session reload. The web path already skips the selected-marker expression/fold TSV exports; it does not export a full expression TSV through MarkerFinder.

## RNA rendering alternatives (first pass)

| Option | Marker stage (s) | Rendering only (s) |
|---|---:|---:|
| Default: 100 cells/state, 2400 DPI, PDF + SVG | 62.16 | 43.04 |
| 25 displayed cells/state | 50.41 | 36.40 |
| 600 DPI | 27.42 | 12.70 |
| PDF only, omit SVG | 44.07 | 28.88 |
| No static figures; keep interactive cache | 15.22 | 0.00 |
| No static figures or plotting matrices/cache | 13.50 | 0.00 |

The 100-cell cap displayed 6,852 cells and the 25-cell cap displayed 1,869, both with the same 4,205 selected RNA markers. Reducing display cells helps less than skipping static rendering because rasterization still draws a large fixed-DPI image. Lower DPI reduces image resolution while preserving embedded vector text (verified with `pdffonts`: embedded CID TrueType with Unicode mapping).

## Validation and availability

- All 33 benchmark cases produced byte-identical unique-marker tables, redundant-marker tables and state-centroid tables within each modality, including across both passes. SHA-256 digests and detailed timings are recorded in the release validation JSON.
- 37 relevant tests passed: output controls, CLI markers-only, cache re-render, unchanged scoring, optional modality outputs, real differential tests on synthetic fixtures, and upload pipeline integration.
- Upload API tests verified interactive heatmap data and on-demand PDF downloads when no static heatmap was generated. A skipped old PDF is excluded from a rerun archive.
- Chrome verified control defaults, disabled irrelevant SVG/DPI menus, and the submitted QC payload in the loaded session. Submission was intercepted so the saved study was not rerun or changed. The COPD viewer loaded successfully with the shared UI changes.
- The local upload server on port 8000 includes the new controls.
- API: `render_heatmap=False`; optional `write_svg=False`, `heatmap_dpi=600`. Disable `write_heatmap_cache`, `write_heatmap_tsv`, and `write_expression_tsv` too to avoid constructing a plotting matrix.
- CLI: `--skip-heatmap-render`, `--skip-svg`, `--dpi 600`, `--markers-only`. The standalone markers-only mode has no interactive cache. scALABLE retains that cache for interactive exploration.

## Reproduction

From the repository root:

```bash
.venv/bin/python tests/benchmark_marker_outputs.py PATH_TO_JOB_JSON NEW_OUTPUT_DIRECTORY
.venv/bin/python tests/benchmark_marker_outputs.py PATH_TO_JOB_JSON ANOTHER_NEW_DIRECTORY --variants no_static markers_only default
```

Platform: macOS-13.7.8-arm64-arm-64bit. Versions: numpy 2.4.6, pandas 2.2.3, scipy 1.17.1, anndata 0.10.9, matplotlib 3.10.9, scanpy 1.10.3.
