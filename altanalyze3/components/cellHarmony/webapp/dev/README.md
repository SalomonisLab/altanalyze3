# Explore and Chat repairs, 2026-08-25

## Explore serving updates, 2026-10-06

Cold CombPlot requests no longer need to occupy an Apache proxy connection while
the plot is calculated. The web/discover route returns `202` with `Retry-After: 1`
for background preparation and `200` with the original payload when ready. Inputs
with at least 10,000 cells use this automatically; the updated Explore and Chat
interfaces explicitly request `deferred=true`. Small legacy requests remain
synchronous. API clients must poll the same URL on `202`; `deferred=false` retains
the synchronous API for existing integrations that explicitly need it.

The background queue has one builder per web process, at most eight distinct
pending/completed entries, and five-minute completed-result retention. Completed
JSON lives in temporary files rather than another expression dictionary. Its
temporary directory is removed at application shutdown; expired files are removed
on subsequent requests. Errors retain their HTTP status and detail. Use the
existing single ASGI web worker: additional web processes have separate queues
and caches and can duplicate builds. Analysis workers retain their existing
configuration and admission policy.

Additional changes preserve the analysis and figure contents:

- The serving bundle computes the existing float64 group contrast and sequential
  float32 totals in one feature pass, with temporary batches bounded by one million
  stored entries. Gene choices are cached per expression entry/grouping (64 KiB,
  eight selections, ten-minute TTL). The formulas and all input features are retained.
- CSR CombPlot requests extract up to 64 requested features in one sparse slice
  instead of rescanning the matrix for each feature. Sampling, five-decimal values,
  gene order, missing-gene reporting, filters and cell identities are unchanged.
- Compact UMAP/expression payloads are constructed directly from columns, avoiding
  hundreds of thousands of temporary Python row dictionaries. Legacy row payloads
  remain available. Compact and legacy builders have exact parity checks.
- The shared browser caches recent encoded responses for five minutes, at most
  eight entries and 384 MiB estimated parsed storage (four times JSON text length).
  This is an estimate for cached responses, not a bound on the browser's total heap
  or Plotly buffers. Full query URLs and the result identity govern reuse; resets,
  job changes and state-layer changes invalidate it. Pending requests are deduplicated,
  and older responses cannot overwrite a newly selected view.
- Large CombPlots retain nonzero positive/negative SVG bars and every zero-valued
  observation's hover point on transparent SVG lines. This removes zero-height
  rectangle nodes without dropping observations, changing values or rasterizing
  PDF exports. Dense modalities may still require many SVG bars.

Validation is recorded in `plot_performance_20261006.json`. A synthetic 600,000-cell,
1,024-feature, 24,576,000-entry serving benchmark ran in the 30 GiB/4 CPU container.
Across the second fresh-process pair, compact UMAP construction fell from 1.019 s
to 0.139 s; JSON serialization still took about 0.61 s. Peak process RSS fell from
1.525 to 1.096 GiB. Entire UMAP and sampled CombPlot payload hashes matched exactly.
The original full analytical workflow was not rerun for this serving change.

On the isolated local web server, cold CombPlot returned `202` in 0.0075 s;
a subsequent cached 7.6 MB result returned `200` in 0.041 s. The browser rendered
50,800 sampled cells and 12 genes using 24,562 bar nodes rather than 609,600 gene
bars. A zero-valued cell's barcode, gene, state, sample and zero value were verified
in its tooltip. The downloaded PDF contained all 12 genes and no raster images.
These are synthetic serving checks, not timings for the reported online visitor
job. Complete browser draw time was not measured. Dense imputed-modality CombPlots,
all-point violin rendering and validation against the visitor's saved job remain.

Deployment: rebuild/restart with `app.py`, `job_bundle.py`, `plot_payload.py`, the
new `plot_build.py`, and the shared `static/app.js` together. Refresh browser pages
so older JavaScript is not interpreting a `202` as a completed plot. No new Python
dependency or Apache-wide `ProxyTimeout` increase is required. Shared browser
transport/caching and SVG improvements apply to web, discover and viewer; the
standalone viewer keeps its specialized immediate CombPlot endpoint.

Reproduce in a scratch directory with `benchmark_plot_builds.py prepare ROOT`,
then `baseline ROOT` and `optimized ROOT` in separate fresh Python processes.
`serve ROOT --port 8012` starts the isolated synthetic UI. These inputs are
computational test data, not biological replicates or inference outputs.

I repaired three faults in the scALABLE webapp and I added two controls. This
folder holds the scripts that prove each change.

## What I changed

All edits sit in
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/`.

1. `app.py` `_state_colors` and `_group_axis`. Both wrote `numpy_array or []`.
   Python calls `bool()` on the array and raises. Every DotPlot and CombPlot
   request returned HTTP 500.
2. `app.py` `_state_colors` also read one label per cell, not one per category.
   The length test never matched, so every cell state came back grey. The
   function now reads the categories. A dataset without stored colours gets the
   12 Paired colours that `static/app.js` draws cell states with.
3. `app.py` `_default_marker_genes`. The old loop sliced one CSR column per
   gene per group. An indicator matrix now gives every group sum in one pass.
4. `app.py` 26 read-only GET handlers changed from `async def` to `def`.
   FastAPI runs a plain `def` handler in its threadpool, so one slow figure no
   longer blocks the other panels.
5. `app.py`, `static/app.js`, `templates/index.html`. The UMAP cell-type view
   takes a colour column and a coordinate set. The app hides the reference
   atlas while either choice stays off-default.
6. `app.py`, `static/app.js`, `templates/index.html`. The Chat tab gained the
   `/api/jobs/{id}/chat` route and example questions built from the job.
7. `app.py` `_get_marker_heatmap_cache_entry`. The reader called `np.load` on a
   cache that `visualization/marker_heatmap_h5ad.py:73` now writes as an h5ad.
   Every MarkerHeatmap request raised. The webapp now calls
   `_read_heatmap_cache`, the writer's own reader, which opens both formats.
8. `app.py`, `static/app.js`, `templates/index.html`. Any two numeric obs
   columns can serve as the axes, the way ShinyCell plots one cell annotation
   against another. Pick "obs columns" in the Coordinates list and the X and Y
   lists appear.

## Measurements

| Check | Before | After |
|---|---|---|
| DotPlot and CombPlot requests | HTTP 500, every request | HTTP 200 |
| DotPlot with an empty gene set | 7.3 min, then HTTP 500 | 0.13 s |
| CSR column slices per figure | 57,168 | 0 |
| Light request under one heavy request | timeout at 5 s | 0.024 s under six |
| Chat question | HTTP 404, no route | HTTP 200 |
| MarkerHeatmap matrix | HTTP 500, every job | HTTP 200 on 3 of 3 jobs |

## Scripts

Run each from this folder with
`/opt/homebrew/opt/python@3.11/bin/python3.11`.

- `validate_default_markers.py` runs the old nested loop and the new one-pass
  code over the same candidate genes. Both must pick the same gene for every
  group. Result: 12 of 12 groups identical, and the largest gap difference over
  480 pairs is 2.4e-15.
- `validate_umap_options.py` builds an AnnData with three embeddings, two
  annotation columns and four numeric columns. It proves the panel reads `obsm`
  and `obs`. It proves an unknown key, a single axis and a boolean axis each
  fall back to the cellHarmony projection and restore the reference. It also
  proves the count of cells the panel cannot place: an axis recorded for 20 of
  40 cells draws 20 points and reports 20 as dropped.
- `validate_explore_endpoints.sh <job id>` sweeps every Explore endpoint of a
  running server and prints each status code and time.
- `restart_app.sh` restarts the server on port 8000. It waits for the old
  process to release the port, so a probe cannot reach the dying process.

## Limits

- Every h5ad on this machine carries one embedding, `X_umap`, and the three
  uploads carry none. The three job h5ads also hold one numeric obs column
  each, `pct_counts_mt`. So `validate_umap_options.py` proves both switches on a
  synthetic AnnData. On the live server I proved the obs axis against a real
  job: the served x values equal `obs['pct_counts_mt']` for all 2,797 cells.
- The MarkerHeatmap repair reproduces the cache exactly. The served TSV names
  the same 943 rows and 500 columns as the h5ad, and 30,000 compared values
  match to the last digit the TSV prints.
- The Chat route asks the LungMAP assistant at
  `http://127.0.0.1:8001/api/assistant/viewer-intent` to read each question.
  The route answers HTTP 503 when that server is down.

## Record of the edits

`applied_patches/` holds the ten scripts that made the edits, in order. Each
script asserts that it finds its target text exactly once, so the set reads as a
precise log of every line I changed. Do not run them again: the source now holds
the new text, and each assertion fails.
