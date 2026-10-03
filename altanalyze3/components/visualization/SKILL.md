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
