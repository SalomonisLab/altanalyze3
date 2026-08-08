# scalable_viewer (scALABLE)

A multi-dataset single-cell browser in pure Python. It replaces the ShinyCell viewer.
The stack is FastAPI plus vanilla JavaScript. No R and no Shiny run anywhere in it.

Location:
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/scalable_viewer/`

## Current architecture: the viewer IS scALABLE

Start the server through `run.py`, which calls `scalable_app.create_scalable_app`. That
function builds scALABLE's own application with
`altanalyze3.components.cellHarmony.webapp.app.create_app()`, then replaces
`app.state.job_store` with `bundle_meta.BundleJobStore`.

Every scALABLE endpoint opens with `store.job_exists(job_id)` and `store.get_job(job_id)`
(app.py:4499 onward) and then calls a builder that takes only `(app, meta, ...)`. The store
swap therefore serves a precomputed bundle through scALABLE's whole application: its routes,
its Jinja template, its stylesheet, its 5,350-line front end, its payload builders and its
matplotlib PDF renderers. No plot, colour ramp or payload is rewritten here.

Two functions would otherwise open an h5ad. `bundle_meta.seed_expression_cache` and
`bundle_meta.seed_marker_heatmap_cache` pre-fill `app.state.expression_cache` and
`app.state.marker_heatmap_cache`, so `_get_expression_cache` (app.py:1353) and
`_get_marker_heatmap_cache_entry` (app.py:1260) both return on a cache hit. A gene column
then comes from one slice of the memory-mapped store.

`prepare_assets.py` runs the scALABLE backends that make the artifacts a bundle lacks:
fastComm scores, per-cell-state marker networks, per-contrast differential networks and
GO-Elite terms.

`server.py` below still exists. `scalable_app` mounts it at `/fast`, so the binary and
memory-mapped endpoints stay available and lose no speed.

## Design of the bundle

The viewer follows the precompute-then-serve pattern of
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/isv_web/`.
`precompute.py` scans the h5ad once and writes a bundle. `server.py` reads only the bundle.
**No endpoint opens an h5ad.** One server process serves many bundles through `/api/catalog`,
and each dataset loads on its first request.

The bundle extends the cellHarmony reference bundle at
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/flask/references/Human/Lung/HLCA/`.
Every legacy file keeps its name and its columns: `<prefix>.txt`, `<prefix>_clusters.tsv`,
`<prefix>_umap.tsv`, `<prefix>_metadata.json`, `<prefix>_config_snippet.json`. The viewer adds
binary sidecars and a `scalable_viewer` block inside the metadata JSON. A legacy reader ignores
the new block; `bundle.discover()` uses it to tell a viewer bundle from a plain reference.

The Morpheus heatmap reuses the embed of
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/app.py`
(lines 349, 377-378 and 573): the Broad CDN scripts, a `fetch` of a TSV, a `Blob`, and
`new morpheus.HeatMap({dataset: blobUrl, ...})`. **Morpheus loads from
`software.broadinstitute.org`.** No Morpheus asset ships in this repository, so the heatmap tab
needs network access. Every other tab works offline.

The job-based `/api/jobs/{job_id}/` design of the webapp is deliberately not copied. That design
runs analyses. This viewer browses published ones.

## Bundle layout

| File | Content |
| --- | --- |
| `<prefix>.txt` | centroid matrix, `UID` + cell-state columns in **canonical order** |
| `<prefix>_clusters.tsv` | `barcode`, cluster key, `Population` |
| `<prefix>_umap.tsv` | `barcode`, `UMAP1`, `UMAP2` |
| `<prefix>_metadata.json` | legacy keys plus the `scalable_viewer` block |
| `<prefix>_config_snippet.json` | cellHarmony registration snippet |
| `<prefix>_cells.npz` | embedding `(N,2)` f32, `state_code` `(N,)` i16, one array per covariate |
| `<prefix>_genes.tsv` | row index to gene id and symbol |
| `<prefix>_stats_mean.npy` | `(G,S)` f32 mean of the layer per gene per cell state |
| `<prefix>_stats_frac.npy` | `(G,S)` f32 fraction of cells above zero |
| `<prefix>_stats_n.npy` | `(S,)` i64 cells per state, the denominator |
| `<prefix>_expr_indptr.npy` | `(G+1,)` i64 gene-major index |
| `<prefix>_expr_indices.npy` | `(nnz,)` u32 cell index |
| `<prefix>_expr_data.npy` | `(nnz,)` f32 layer value |
| `<prefix>_markers.tsv` | normalised marker table |
| `<prefix>_deg/*.tsv` | every `DEG_detailed_*` / `DEG_pooled_overall_*` copied verbatim |
| `<prefix>_deg_manifest.json` | comparison id to file, columns, row count |
| `<prefix>_ccc.tsv` | cell-cell communication edges, optional |

The gene-major store makes a gene overlay instant. The h5ad holds a cell-major CSR, so
reading one gene means scanning every cell. `precompute.py` transposes it once with a two-pass
counting sort, so the server reads one contiguous memory-mapped slice.

## Methods, stated plainly

* **Statistics.** `mean` is the unweighted mean of the layer over the cells of a state. The code
  does not weight a metacell by the cell count behind it. `frac` counts the cells above zero.
  Zero values are never dropped: a violin fills each state to its known size with zeros.
* **Highly variable genes.** `dispersion = variance / mean` on the layer values, 20 equal-count
  bins of mean, dispersion z-scored inside each bin, top *n* kept. I describe the method here in
  full. It does **not** call `scanpy.pp.highly_variable_genes`.
* **Embedding.** z-score the HVG matrix, clip at plus or minus 10, `randomized_svd` PCA, then
  UMAP. `--embedding-from obsm:<key>` reuses an embedding the h5ad already carries.
* **Expression store precision.** `--expr-dtype float32` (the default) copies the source float32
  values unchanged. `--expr-dtype float16` halves the size and rounds the values.
* **Cell-state order and colour** come from `--order` or from `uns['lineage_order']` and
  `uns['cluster_colors_json']`. The client never re-sorts. Alphabetical order is never used.

## Run it

Precompute a bundle (once per dataset):

```bash
PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 \
/opt/homebrew/opt/python@3.11/bin/python3.11 \
  -m altanalyze3.components.visualization.scalable_viewer.precompute \
  --h5ad    /Users/saljh8/Dropbox/LungMAP/LungMAP.net/Datasets/COPD-atlas/results/COPD_metacells.deid.h5ad \
  --out     /Users/saljh8/Dropbox/LungMAP/LungMAP.net/Datasets/COPD-atlas/scalable_viewer/bundles/COPD-metacells \
  --prefix  Hs-Lung-COPD-metacells \
  --dataset-id COPD-metacells \
  --label   "COPD lung atlas (metacells)" \
  --cluster-key cell_state \
  --layer   lognorm \
  --markers /Users/saljh8/Dropbox/LungMAP/LungMAP.net/Datasets/COPD-atlas/analysis/markers/COPD_metacell_markers.tsv \
  --order   /Users/saljh8/Dropbox/LungMAP/LungMAP.net/Datasets/COPD-atlas/analysis/canonical_HLCA_cell_state_order.tsv \
  --deg     /Users/saljh8/Dropbox/LungMAP/LungMAP.net/Datasets/COPD-atlas/analysis/differential \
  --n-hvg 2000 --n-pcs 50 --expr-dtype float32
```

Serve every bundle under one root:

```bash
PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 \
/opt/homebrew/opt/python@3.11/bin/python3.11 \
  -m altanalyze3.components.visualization.scalable_viewer.run \
  --root /Users/saljh8/Dropbox/LungMAP/LungMAP.net/Datasets/COPD-atlas/scalable_viewer/bundles \
  --host 127.0.0.1 --port 8060
```

Then open `http://127.0.0.1:8060`. Add a dataset by precomputing another bundle under the same
root and restarting the server.

## Endpoints

| Method and path | Returns |
| --- | --- |
| `GET /` | the viewer shell |
| `GET /healthz` | dataset count, which datasets are loaded, load errors |
| `GET /api/catalog` | every dataset, from metadata JSON only |
| `GET /api/dataset/{ds}/meta` | states, colours, covariates, methods, warnings, paths |
| `GET /api/dataset/{ds}/genes?q=&limit=` | gene search |
| `GET /api/dataset/{ds}/embedding.bin` | raw float32 `[x0,y0,x1,y1,...]` |
| `GET /api/dataset/{ds}/colorby?key=` | legend for a colouring |
| `GET /api/dataset/{ds}/colorby.bin?key=` | raw codes (int16, int32) or values (float32) |
| `GET /api/dataset/{ds}/expression.bin?gene=` | `[n:int32][idx:u32*n][val:f32*n]` |
| `GET /api/dataset/{ds}/violin?gene=&bins=` | per-state histogram, quantiles, mean, frac |
| `GET /api/dataset/{ds}/dotplot?genes=` | mean and frac per gene per state |
| `GET /api/dataset/{ds}/markers?state=&top=` | marker rows |
| `GET /api/dataset/{ds}/heatmap.tsv?states=&top=&row_scale=` | Morpheus TSV |
| `GET /api/dataset/{ds}/morpheus?top=&row_scale=&scheme=` | the Morpheus iframe page |
| `GET /api/dataset/{ds}/deg` | the comparisons in the bundle |
| `GET /api/dataset/{ds}/deg/table?id=&fdr_max=&state=&max_rows=` | table plus volcano points |
| `GET /api/dataset/{ds}/ccc?limit=` | communication edges, or `available:false` |

`colorby` accepts `cell_state` or any covariate name from `meta.covariates`.

## Deep links

Every view has a URL, so a colleague can open the exact plot you saw:

```
http://127.0.0.1:8060/?ds=COPD-metacells&tab=embed&gene=SFTPC
http://127.0.0.1:8060/?tab=embed&colorby=FEV1
http://127.0.0.1:8060/?tab=expr&genes=SFTPC,SCGB1A1,COL1A1,MARCO
http://127.0.0.1:8060/?tab=deg&deg=COPD_vs_non-COPD::per_cell_state&fdr=0.05
http://127.0.0.1:8060/?tab=heat&top=5
```

## Validate a bundle

Run this after every precompute. It re-reads the source h5ad and fails loudly on any
mismatch.

```bash
PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 \
/opt/homebrew/opt/python@3.11/bin/python3.11 \
  -m altanalyze3.components.visualization.scalable_viewer.validate \
  --bundle /Users/saljh8/Dropbox/LungMAP/LungMAP.net/Datasets/COPD-atlas/scalable_viewer/bundles/COPD-metacells \
  --prefix Hs-Lung-COPD-metacells --n-cells 150 --n-genes 40
```

## Tabs

1. **Embedding** - UMAP on a canvas, coloured by cell state or by any covariate, with a gene
   overlay on a cyan-to-yellow ramp.
2. **Expression** - violin per cell state and a dot plot, both in canonical order.
3. **Marker heatmap** - Morpheus, top *n* markers per state, optional row z-score.
4. **Differential** - a precomputed comparison as a volcano plot and a table.
5. **Cell-cell communication** - the edge table when the bundle has one.
6. **Provenance** - counts, methods, source paths, covariates and every precompute warning.

## Measured on the COPD dataset

Source: `/Users/saljh8/Dropbox/LungMAP/LungMAP.net/Datasets/COPD-atlas/results/COPD_metacells.deid.h5ad`,
161,432 metacells x 36,249 genes, 975,729,760 non-zeros in `layers['lognorm']`.
Bundle: `/Users/saljh8/Dropbox/LungMAP/LungMAP.net/Datasets/COPD-atlas/scalable_viewer/bundles/COPD-metacells`,
7.3 GB, built in 543 s on a 64 GB M-series Mac. Log: `precompute.log` in that directory.

| Step | Time |
| --- | --- |
| count non-zeros per gene (pass 1) | 34.9 s |
| transpose to gene-major + statistics (pass 2) | 309.8 s |
| pick 2,000 HVGs and materialise the matrix | 0.6 s |
| PCA, 50 components | 67.7 s |
| UMAP, 161,432 x 50 | 160.2 s |
| write every bundle file | 0.6 s |

Server: startup 0.6 s, resident memory 106 MB while serving the 7.3 GB bundle, because the
arrays are memory-mapped. Response times, warm:

| Endpoint | Time |
| --- | --- |
| `/api/catalog`, `/meta`, `/deg`, `/markers` | 3 to 4 ms |
| `/embedding.bin` (1.29 MB) | 4 ms |
| `/expression.bin` for 9 different genes | 3 to 29 ms |
| `/violin` | 12 to 17 ms |
| `/dotplot`, 6 genes | 14 ms |
| `/heatmap.tsv`, 117 rows x 39 columns | 8 ms |
| `/deg/table`, 314 rows plus volcano | 9 ms |

Nothing exceeded 30 ms.

## Limits

* The heatmap tab needs the Broad CDN. Everything else runs offline.
* `--skip-expr` cannot build an embedding, because the per-gene variance is not stored. Pass
  `--embedding-from obsm:<key>` with it, or rerun the full precompute.
* The server reads the catalog at startup. A bundle added later needs a server restart.
* The DEG tab reads the tables the differential workflow wrote. It recomputes no statistic.
