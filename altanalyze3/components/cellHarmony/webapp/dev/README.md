# Explore and Chat repairs, 2026-08-25

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
