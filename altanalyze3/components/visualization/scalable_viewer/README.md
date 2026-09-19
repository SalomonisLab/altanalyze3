# scALABLE-viewer (scalable_viewer)

scALABLE-viewer serves a published single-cell atlas through the scALABLE application without a
Run tab. `precompute.py` reads an h5ad once and writes a bundle of memory-mapped arrays;
the server reads only the bundle. One process serves many bundles, and no endpoint opens an
h5ad. The viewer replaces the ShinyCell viewer. The stack is FastAPI plus the web app's own
JavaScript. No R runs anywhere in it.

Location:
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/scalable_viewer/`

This file covers architecture, the command-line tools that build and serve a bundle, the bundle
format, and the computations the viewer adds. Every default below comes from the source as of
2026-09-14.

| Question | File |
| --- | --- |
| How do I use the viewer? | `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/scalable_viewer/HOW_TO_USE.md` |
| What did the chat validation measure? | `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/scalable_viewer/VALIDATION.md` |
| What do the statistics, imputation models and differential thresholds mean? | `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/README.md` |
| Which controls does the shared interface offer? | `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/HOW_TO_USE.md` |

## Architecture

`run.py` calls `scalable_app.create_scalable_app`, which builds the scALABLE application with
`altanalyze3.components.cellHarmony.webapp.app.create_app()` and then replaces
`app.state.job_store` with `bundle_meta.BundleJobStore`. Every scALABLE endpoint opens with
`store.get_job(job_id)` and passes the job metadata to a builder, so a bundle flows through the
web app's routes, template, stylesheet, front end, payload builders and PDF renderers unchanged.
The viewer serves `index.html`, `app.js` and `styles.css` straight out of the webapp directory.

`bundle_meta.seed_expression_cache` and `seed_marker_heatmap_cache` pre-fill the two caches
that would otherwise call `read_h5ad`, one entry per dataset and modality, with a `BundleAnnData`
wrapper over the memory-mapped store. A gene column is then one contiguous slice.

Where a web app handler would need `adata.X`, the viewer removes the route and registers its own.

| Route | Viewer action | Reason |
| --- | --- | --- |
| `POST /api/jobs/{id}/chat`, `GET .../chat-examples` | replaced | the web app version reads `adata.X` |
| `GET /api/jobs/{id}/dotplot`, `.../combplot` | replaced | read `stats_mean`, `stats_frac` and the store directly |
| `GET /api/jobs/{id}/genes`, `.../differential/interactive/gene` | replaced | resolve a feature's display name to its store key |
| `POST /api/jobs/{id}/differential` | blocked | a bundle runs nothing; the route returns the current payload |
| `GET .../expression`, `.../expression/pdf` | wrapped | adds a `covariate` parameter for the violin |
| `GET .../differential/interactive/go` | wrapped | adds three significance tiers and a legend |
| `GET /api/jobs/{id}/state-colors` | added | the bundle's own cell-state colours |
| `GET /api/catalog`, `GET /api/study` | added | dataset list and the Study tab record |
| `GET .../differential/contrasts`, `POST .../differential/select` | added | list and switch precomputed contrasts |
| `GET .../grn/regulator-network`, `.../grn/tf-activity` | added | GRN views over a bundle, see below |

A middleware rewrites the served page: the title becomes `scALABLE-viewer`, the two documentation
links point at this directory, and `viewer.css` plus `viewer_bootstrap.js` load before `</body>`.
A second middleware sets `Cache-Control: no-cache, must-revalidate` on every static file, because
an ETag without it let browsers replay stale JavaScript. `server.py`, the earlier binary API,
stays mounted at `/fast`.

## Run the server

```bash
PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 \
/opt/homebrew/opt/python@3.11/bin/python3.11 \
  -m altanalyze3.components.visualization.scalable_viewer.run \
  --root /path/to/bundles --assets /path/to/assets --state-dir /path/to/runtime \
  --host 127.0.0.1 --port 8062
```

| Flag | Default | Meaning |
| --- | --- | --- |
| `--root` | none | directory tree holding precomputed bundles |
| `--catalog` | none | JSON `{"datasets":[{"bundle_dir":..,"prefix":..}]}`; relative `bundle_dir` resolves against the file |
| `--assets` | none | asset root written by `prepare_assets.py`; without it every asset-driven view is empty |
| `--state-dir` | `<root>/../viewer_runtime` | per-dataset logs; nothing is written into a bundle |
| `--host` | `127.0.0.1` | |
| `--port` | `8062` | |

Give `--root` or `--catalog`. The server exits with code 2 when it finds no servable
bundle. It prints one line per dataset, `id: N cells x M genes, S states`, then
`serving scALABLE on http://host:port`. The catalog loads at start; a bundle added later needs a
restart.

| Environment variable | Default | Effect |
| --- | --- | --- |
| `SCALABLE_ASSISTANT_URL` | `http://127.0.0.1:8001/api/assistant/viewer-intent` | the chat intent router |
| `LUNGMAP_SITE_DB` | `/Users/saljh8/Dropbox/LungMAP/refactored_website/lungmap-data/site/breath.sqlite` | the study records, opened read-only |
| `LUNGMAP_SOURCE_TABLES` | `/Users/saljh8/Dropbox/LungMAP/refactored_website/build/lungmap-data-new/data/metadata` | fallback study tables |
| `LUNGMAP_STUDY_IDS` | empty | deployment-wide study id list |
| `LUNGMAP_SITE_BASE` | `http://127.0.0.1:8001` | base URL for links on the Study tab |
| `CELLHARMONY_ASSISTANT_URL`, `CELLHARMONY_MAX_FILES`, `CELLHARMONY_JOB_WORKERS` | as in the web app | still read; the viewer overrides `JOB_STORAGE`, `ROOT_PATH` and the template paths |

The viewer needs the same dependencies as the web app (`requirements.docker.txt` in the webapp
directory) plus `h5py`, `scikit-learn` and `umap-learn` for `precompute.py`. Its own
`requirements.txt` lists only the serving subset and omits `pandas`, which `bundle_meta.py` and
`prepare_assets.py` import.

## Run the server in Docker

```bash
./build.sh                                              # image: scalable-viewer:latest
./run.sh --bundles /path/to/bundles --assets /path/to/assets
./deploy.sh --tag 2026-09-16                            # build, then push to ECR
```

`run.sh` publishes port 8003 and replaces any container of the same name, so it is also the
restart. It prints the URL, the health URL and the two commands for logs and stop.

| Script | What it does |
| --- | --- |
| `build.sh` | builds with `altanalyze3/components` as the context, because the viewer imports the cellHarmony web app, then imports `create_scalable_app` inside the image so a missing component fails the build rather than a request |
| `run.sh` | mounts the bundles and assets read-only, mounts a writable state directory, and starts the container detached with `--restart unless-stopped` |
| `deploy.sh` | builds, creates the ECR repository when absent, logs in and pushes; it restarts nothing on the target host |

| Path in the container | Source |
| --- | --- |
| `/data/bundles` | `--bundles`, read-only; `VIEWER_ROOT` inside |
| `/data/assets` | `--assets`, read-only; `VIEWER_ASSETS` inside |
| `/srv/scalable_viewer/state` | `--state`, writable, defaults to `viewer_runtime` beside the scripts |

The image holds the code only. `PORT`, `BIND`, `IMAGE`, `NAME` and `STATE` override the
defaults; `SCALABLE_ASSISTANT_URL`, `LUNGMAP_SITE_DB`, `LUNGMAP_SOURCE_TABLES`,
`LUNGMAP_STUDY_IDS` and `LUNGMAP_SITE_BASE` pass through when set, and the two that name a
host path are mounted at that same path inside.

A bundle records two absolute paths, and `run.sh` mounts both read-only at those same paths
when the host has them. A catalog's `bundle_dir` is the dataset itself, since `--catalog` mounts
the file at its own path and relative entries resolve against it. `source_h5ad` is the h5ad
precompute read, which the viewer does not need: `2e98942` made the expression cache bundle-owned,
so a bundle serves wherever it is copied. Mounting the file where it exists still gives the Gene
Detail fallback for a gene the comparison does not carry.

`/fast/healthz` is the health check, on the binary API mounted at `/fast`. It answers
`{"ok": true, "n_datasets": N, "loaded": [...], "load_errors": [...]}` without rendering a page;
the container's own HEALTHCHECK calls it.

The entrypoint builds run.py's flags from the environment. Arguments passed to `docker run`
replace them, so `docker run ... scalable-viewer --port 9000` re-flags the server and
`docker run ... scalable-viewer python -m ...precompute --h5ad ...` builds a bundle with the
image's own interpreter.

## Build a bundle: precompute.py

```bash
PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 \
/opt/homebrew/opt/python@3.11/bin/python3.11 \
  -m altanalyze3.components.visualization.scalable_viewer.precompute \
  --h5ad /path/to/data.h5ad --out /path/to/bundles/MyDataset --prefix My-Dataset \
  --label "My dataset" --cluster-key cell_state --layer lognorm \
  --markers /path/to/markers.tsv --order /path/to/canonical_order.tsv \
  --deg /path/to/differential
```

| Flag | Default | Meaning |
| --- | --- | --- |
| `--h5ad` | none | source h5ad; required unless `--modalities-only` |
| `--out` | required | bundle directory, created if absent |
| `--prefix` | required | file-name prefix inside the bundle |
| `--label` | none | human label for the catalog |
| `--dataset-id` | `--prefix` | catalog id |
| `--study-id` | none | LungMAP study id, for example `lmdata:LMEX0000009416`; without it the Study tab shows no record |
| `--cluster-key` | `cell_state` | obs column holding the cell state; must be categorical or string |
| `--layer` | `lognorm` | layer to serve; `X` for `adata.X` |
| `--markers` | none | marker table |
| `--order` | none | canonical order TSV with `order`, `cell_state`, `color` columns |
| `--deg` | none | directory tree holding the RNA `DEG_*.tsv` tables |
| `--deg-modality ID=DIR` | none | one modality's own `DEG_*.tsv` tree; repeatable |
| `--ccc` | none | cell-cell communication TSV |
| `--n-hvg` | `2000` | highly variable genes for the embedding |
| `--n-pcs` | `50` | principal components |
| `--n-neighbors` | `15` | UMAP neighbours |
| `--min-dist` | `0.3` | UMAP min_dist |
| `--seed` | `0` | |
| `--row-block` | `8192` | rows per transpose block |
| `--expr-dtype` | `float32` | `float32` copies the values; `float16` halves the store and rounds |
| `--embedding-from` | none | `obsm:<key>` to reuse an embedding instead of computing one |
| `--max-centroid-genes` | `4000` | rows in `<prefix>.txt`, marker genes first |
| `--skip-expr` | off | reuse an existing store and statistics; needs `--embedding-from` |
| `--modality ID=PATH` | none | add a modality store; repeatable; see below |
| `--modality-label ID=LABEL` | none | menu label for a modality |
| `--modality-feature-label ID=NOUN` | none | what one feature is called, for example `lipids=lipid` |
| `--modality-display-names ID=TSV` | none | two-column `feature<TAB>display`; the key stays searchable |
| `--modalities-only` | off | write only the modality sidecars into an existing bundle |

The script stops with an error when the cluster key is absent or numeric, when a cell has no
state, or when no canonical order exists in `--order` or `uns['lineage_order']`. A
`--modality-label` for an undeclared id, or `--modalities-only` without a `--modality` or
`--deg-modality`, is also an error. It returns 0 on success and raises when the bundle is
incomplete at the end.

### Methods

| Quantity | Definition |
| --- | --- |
| `stats_mean` | unweighted mean of the layer over the cells of a state; a metacell counts once |
| `stats_frac` | fraction of the state's cells with a value above zero |
| `stats_n` | cells per state, the denominator |
| highly variable genes | `dispersion = variance / mean` on the layer, 20 equal-count bins of mean, z-score inside each bin, top n; not `scanpy.pp.highly_variable_genes` |
| embedding | z-score the HVG matrix, clip at plus or minus 10, randomised-SVD PCA, UMAP |
| expression store | one gene-major CSC store, built by a two-pass counting sort, so one gene is one contiguous slice |
| cell-state order and colour | `--order`, else `uns['lineage_order']` and `uns['cluster_colors_json']`; the client never re-sorts |

## Add modalities to a bundle

`--modality adt=/path/to/predictions.csv` adds a second feature matrix over the same cells.
`precompute.py` sniffs the file's bytes: HDF5 opens as h5ad, anything else as CSV or TSV, so an
`rna2*` prediction written as CSV under an `.h5ad` name still reads correctly. The row labels
decide the store kind.

| Row labels | Store kind | Files | Explore |
| --- | --- | --- | --- |
| cell barcodes | `per_cell` | genes, stats, and the three CSC arrays | offered |
| cell-state names | `per_state` | genes and stats only; each cell takes its state's value | not offered |

A `per_state` store gives every metacell of a population one value, so a UMAP would show flat
patches and invite a within-state reading it cannot support. On the COPD-metacells bundle the GRN
edge `ARNT|GAL3ST4` took 11 distinct values over 83,416 metacells, uniform within all 50 of 50
populations. The rule keys on store granularity, not on the modality name, and it leaves the
differential list untouched. The log reports per-cell retention as matched of total cells and
prints a warning below 90 percent; unmatched cells keep zeros.

`--modalities-only` re-derives barcodes, states and counts from the bundle's own files, fails on
any disagreement, writes the sidecars, merges modality DEG roots into the manifest, and stamps
`modalities_updated_utc`. Adding a modality this way takes minutes; a full rebuild takes hours.

## Prepare the assets: prepare_assets.py

The bundle holds the matrices. The assets hold the artifacts a view needs beyond them: the
MarkerFinder heatmap cache, marker networks, fastComm scores, differential networks, GO-Elite
terms and fold matrices. The script writes one file, `<out>/<prefix>_assets.json`, plus the
network TSVs it exports.

```bash
PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 \
/opt/homebrew/opt/python@3.11/bin/python3.11 \
  -m altanalyze3.components.visualization.scalable_viewer.prepare_assets \
  --bundle-dir <dir> --prefix <prefix> --out <asset dir> --markers-tsv <marker table> \
  --goelite-from <differential root> --diff-networks-from <differential root> \
  --fastcomm-dir <fastComm run dir> --study-id lmdata:LMEX0000009416
```

| Flag | Default | Meaning |
| --- | --- | --- |
| `--bundle-dir`, `--prefix`, `--out` | required | `--out` is never the bundle directory |
| `--markers-tsv` | none | marker table with `Gene`, `Fold`, `cluster` |
| `--marker-heatmap-npz` | `*_fold_matrix.npz` beside `--markers-tsv` | the MarkerFinder heatmap cache |
| `--marker-heatmap-npz-full` | none | a second, all-column MarkerFinder run offered as the `All cells` density |
| `--marker-gct` | none | convert a per-cell GCT when the npz is gone |
| `--no-marker-heatmap` | off | build without the marker heatmap, on purpose |
| `--network-markers-tsv` | `*_redundant_markers.tsv` beside the marker table | 250 genes per state for the marker networks |
| `--skip-networks` | off | no marker networks and no differential networks |
| `--order-tsv` | none | canonical cell-state order |
| `--study-id` | none | the Study tab record; overrides the bundle's own |
| `--fastcomm-dir` | `<scalable_viewer>/fastComm/<dataset id>/` | finished fastComm run |
| `--no-fastcomm` | off | build without cell communication, on purpose |
| `--allow-stale-fastcomm` | off | ship a run whose cell count, states or state key disagree with the bundle |
| `--fastcomm-sample-key` | empty | obs column of the per-sample split |
| `--species` | `human` | |
| `--goelite-from` | none | ingest the differential's own `GeneSetEnrichment/` tables |
| `--goelite` | off | recompute GO-Elite from the DEG table; enriches every row instead of the assigned-group lists and omits `coreg_*`; prefer `--goelite-from` |
| `--no-goelite` | off | build without GO terms, on purpose |
| `--diff-networks-from` | none | ingest the differential's own `interaction-plots/` |
| `--folds-from` | `--diff-networks-from`, then `--goelite-from` | each contrast's full feature by cell-state fold matrix |
| `--differential-from-modality ID=DIR` | none | one modality's own differential root; repeatable; never read from the RNA root |
| `--no-folds` | off | draw an uncalled cell state as 0 in the heatmap |

`--goelite` with `--goelite-from`, `--no-marker-heatmap` with a heatmap source, and
`--no-fastcomm` with `--fastcomm-dir` are errors.

The script ends with a completeness report. Each input is present, or a named flag waived it.
A build that satisfies neither prints `MISSING` and returns exit code 2.

| Viewer input | Source | Waiver |
| --- | --- | --- |
| Explore / MarkerHeatmap | the npz beside `--markers-tsv` | `--no-marker-heatmap` |
| Explore / MarkerNetwork | the redundant marker table | `--skip-networks` |
| Explore / Cell communication | `--fastcomm-dir` | `--no-fastcomm` |
| Differential / Network | `--diff-networks-from` | `--skip-networks` |
| Differential / GO Terms | `--goelite-from` | `--no-goelite` |
| Differential / folds | `--folds-from` | `--no-folds`, or no fold source given |

The gate exists because a build before 2026-08-12 returned 0 with no heatmap and no cell
communication, and the viewer served an Explore tab missing both without any log line.

A fastComm run must match the bundle. `verify_fastcomm_matches_bundle` compares the run's
`summary.json` with the bundle on state key, cell count, and cell states. When fastComm skipped
a state through `--min-cells`, the check reports it and continues. The check exists
because one run scored 161,432 cells of one h5ad while the bundle served 123,076 cells of
another, both with 39 states, so no number on the page looked wrong.

The manifest also records `study_id`, `dotplot_default`, `marker_fold_lookup_tsv`, and per-contrast
entries with their modality. The server accepts the manifest at `<assets>/<prefix>_assets.json`
or `<assets>/<id>/<prefix>_assets.json`, and resolves relative paths against the bundle root's
third parent.

## Validate a bundle: validate.py

```bash
PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 \
/opt/homebrew/opt/python@3.11/bin/python3.11 \
  -m altanalyze3.components.visualization.scalable_viewer.validate \
  --bundle /path/to/bundles/MyDataset --prefix My-Dataset --n-cells 200 --n-genes 40 --seed 7
```

`--h5ad` defaults to `source_h5ad` from the metadata. The script re-reads the source and runs
eight checks, each reported with its denominator:

| Check | Passes when |
| --- | --- |
| A | non-zero count of the source layer equals the gene-major store |
| B | every non-zero of N random cells sits in the store with the same value |
| C | cell indices ascend strictly inside each gene |
| D | `stats_mean` and `stats_frac` recomputed from the store match the stored matrices within 1e-4 |
| E | `stats_n` matches the cell-state counts in obs |
| F | the embedding is finite with one row per cell |
| G | the legacy TSVs hold exactly one row per cell or per gene |
| H | the centroid matrix columns follow the canonical cell-state order |

Exit code 0 when every check passes, 1 otherwise. Run it after every precompute.

## Bundle layout

Every legacy cellHarmony reference file keeps its name and columns, so a legacy reader ignores
the additions. `bundle.discover()` treats a directory as a bundle when its `*_metadata.json`
holds a `scalable_viewer` block, so pointing the server at a shared reference tree is safe.

| File | Content |
| --- | --- |
| `<prefix>.txt` | centroid matrix, `UID` plus cell-state columns in canonical order, up to `--max-centroid-genes` rows |
| `<prefix>_clusters.tsv` | `barcode`, cluster key, `Population` |
| `<prefix>_umap.tsv` | `barcode`, `UMAP1`, `UMAP2` |
| `<prefix>_metadata.json` | legacy keys plus the `scalable_viewer` block |
| `<prefix>_config_snippet.json` | cellHarmony registration snippet |
| `<prefix>_cells.npz` | `embedding` (N,2) f32, `state_code` (N,) i16, `cov_num_<name>` and `cov_cat_<name>` arrays |
| `<prefix>_genes.tsv` | row index, gene id, symbol, and a `display` column when display names were given |
| `<prefix>_stats_mean.npy` | (G,S) f32 mean per gene per state |
| `<prefix>_stats_frac.npy` | (G,S) f32 fraction above zero |
| `<prefix>_stats_n.npy` | (S,) i64 cells per state |
| `<prefix>_expr_indptr.npy`, `_expr_indices.npy`, `_expr_data.npy` | gene-major CSC store: (G+1,) i64, (nnz,) u32, (nnz,) f32 |
| `<prefix>_markers.tsv` | normalised marker table |
| `<prefix>_deg/<file>.tsv`, `<prefix>_deg_manifest.json` | every `DEG_detailed_*` and `DEG_pooled_overall_*` table copied verbatim, plus the manifest |
| `<prefix>_ccc.tsv` | cell-cell communication edges, optional |
| `<prefix>_<modality>_genes.tsv`, `_stats_mean.npy`, `_stats_frac.npy`, `_expr_*.npy` | one sidecar set per modality; `per_state` stores omit the `_expr_*` files |

The server refuses a bundle without `metadata`, `cells`, `genes`, the three `stats` files and
the three `expr` files. The `.txt`, `_clusters.tsv`, `_umap.tsv`, `_markers.tsv`, `_ccc.tsv`
and DEG files are optional.

The `scalable_viewer` metadata block records `bundle_version` (1), `id`, `label`, `study_id`,
`prefix`, `built_utc`, `layer`, `expr_dtype`, `nnz`, `n_states`, `states`, `state_n`,
`canonical_order_source`, `covariates`, `embedding_method`, `stats_method`, `hvg_method`,
`markers`, `deg`, `ccc`, `modalities`, `centroid_genes`, and `warnings`. Each modality entry
records `id`, `label`, `feature_label`, `kind`, `source`, `source_mtime`, `n_features`,
`n_matched`, `n_expected`, `retention`, `nnz` and `expr_dtype`.

DEG manifest ids follow one convention. An RNA comparison is `<contrast>::<comparison>::<kind>`
and any other modality is `<modality>::<contrast>::<comparison>::<kind>`, where `kind` is
`per_cell_state` or `pooled_overall`. Non-RNA tables sit in a per-modality subdirectory because
two modalities share a file name. The front end never parses the `::` string; it reads the
manifest's own `modality` field.

## Catalog and Study tab

`data_api.build_catalog` merges `--root` discovery with the `--catalog` file and de-duplicates on
`(bundle_dir, prefix)`. The catalog reads only metadata, so start-up costs under a second, and
each dataset loads on its first request. A duplicate id or an incomplete bundle lands in
`load_errors`, which `/api/catalog` returns beside the dataset rows. Each row carries `id`,
`label`, `n_cells`, `n_genes`, `n_states`, `contrasts`, `has_markers`, `fastcomm` and
`default_gene`.

`/api/study` resolves one study id per dataset, most specific source first:

1. `study_id` in the asset manifest, from `prepare_assets.py --study-id`
2. `scalable_viewer.study_id` in the bundle metadata, from `precompute.py --study-id`
3. `LUNGMAP_STUDY_IDS`
4. nothing, and the Study tab reports that no source names an id

The manifest wins because a manifest rebuild costs minutes and a bundle rebuild costs hours.
Step 4 never guesses. Until 2026-08-12 the module carried a default id, and the COPD bundle
served the Study tab of a different study. The record comes from the LungMAP site database,
opened read-only, with the metadata tables as a fallback.

## Computations the viewer adds

| Feature | Definition |
| --- | --- |
| DotPlot fast path | with no subset filter and the full state list, `mean` and `frac` come straight from `stats_mean` and `stats_frac`; otherwise the store is scanned, and `frac` counts stored non-zero entries |
| DotPlot default genes | the marker table's top gene per state; an imputed modality without markers uses its first 12 features |
| CombPlot | individual cells by default, in cell-state order; `Display: Donor means` explicitly enables averaging (`unit=donor`, `min_cells` default 5). Covariate bands use each cell's annotations or the donor group's modal annotation. Column details appear on hover only. |
| plot variables | categorical covariates with 2 to 60 levels; numeric axes exclude `n_counts`, `n_cells`, `n_cells_total`, `n_genes_detected`, `metacell`, `n_donors`, `meta_sample_n_donors`, any `*__n_obs`, and any field with a missing value |
| state colours | the bundle's `cluster_colors` when every requested label is a cell state; otherwise the web app's ramp |
| violin covariate | the violin groups by any categorical covariate; the default path goes through the web app's renderer so the default PDF is unchanged |
| GO tiers | `Representative` in `#1f19c7` when GO-Elite selected the term, `Significant` in `#60a5fa` when `FDR <= 0.05` without selection, `Other` in `#d1d5db` |
| feature display names | `/genes` returns the display column and a `keys` list, so a reader's name such as `18:2 Cholesterol ester` resolves to the store key |
| contrast select | rebuilds the differential block for the chosen id, invalidates the cache, refreshes display names; it computes nothing |
| marker heatmap density | `heatmap_cache` (10 metacells per population) and `heatmap_cache_full`; on the COPD v7 build the two runs share 760 of 1,248 marker rows |

### GRN regulator network and TF activity

`grn_network.regulator_network` builds a network for one cell state and contrast in four steps.

| Step | Rule | Default |
| --- | --- | --- |
| targets | the top rows of the contrast in the state at `FDR <= 0.05`, re-ranked by absolute log2 fold; `features=` overrides | `limit` 200, maximum 2,000 |
| edges | edges from the `grn` store whose score in the state clears a percentile of that state's own edge scores | `edge_percentile` 50 |
| factors | a factor must clear a percentile of all modelled factors' RNA mean in the state; factors below are reported as silent | `expression_percentile` 50 |
| per-factor statistics | activity fold and FDR from the sibling `grn_tf` contrast, expression fold and FDR from the RNA contrast, kept apart on each node | |

A bundle averages over every metacell of the state, not over the contrast's donors; the payload
states this in `arm_scope`. `tf_activity_profile` returns each factor's activity in one state,
or the mean over states, with the contrast's fold and FDR; a factor with no row in the contrast
did not pass that run's fold and FDR gates. Both routes read a `grn_tf` store beside the `grn`
edge store, which the COPD bundle carries as `['rna','adt','grn','grn_tf','lipid']`.

## Chat in the viewer

The viewer sends each question to the router at `SCALABLE_ASSISTANT_URL` with the bundle's
states, contrasts, covariates and modalities, then computes the answer from the bundle. The web
app README describes the router. The viewer answers more protocols than the web app because a
bundle carries per-state statistics for every cell.

| Protocol | Executor | Plot kind |
| --- | --- | --- |
| `cell_identity`, `state_comparison`, `expression_lookup`, `state_contrast` and its aliases | as in the web app, from `stats_mean` weighted by `state_n` and from the DEG tables at `FDR < 0.05` | `dotplot`, `volcano` |
| `severity_gradient` | per-donor pseudobulk against a numeric covariate | `gradient` |
| `coexpression_module` | per-donor co-expression around a seed gene, top 5 partners | `combplot` |
| `composition_shift` | per-donor cell-state composition between two groups | `frequency` |
| `annotation_concordance` | cross-tabulation of two annotations | `heatmap` |
| `dose_response` | expression across ordered stages | `combplot` |
| `donor_heterogeneity` | a signature score per donor | `signature` |
| `most_affected_state` | states ranked by significant genes, with their abundance per group | `frequency` or `barchart` |
| `pathway_program` | GO-Elite terms of the contrast in the state | `barchart` |
| `regulatory_driver` with the `grn_tf` modality | `tf_activity_profile` | `barchart` |
| `regulatory_driver` with a state | `regulator_network` | `network` |

The viewer returns up to 4 follow-up questions built from the vocabulary of the answer, and it
trusts the router's slots without the web app's sentence repair. When the router is down the
route answers HTTP 503. `VALIDATION.md` records 55 of 55, 68 of 68 and 15 of 15 passing
questions on the COPD-metacells bundle on 2026-08-26, with the router calling the model 0 times
in the 68-question sweep, and lists the scope limits of that measurement.

## Measured on the COPD dataset

These timings describe the 161,432-metacell build of
`/Users/saljh8/Dropbox/LungMAP/LungMAP.net/Datasets/COPD-atlas/results/COPD_metacells.deid.h5ad`
(36,249 genes, 975,729,760 non-zeros in `layers['lognorm']`) into
`/Users/saljh8/Dropbox/LungMAP/LungMAP.net/Datasets/COPD-atlas/scalable_viewer/bundles/COPD-metacells`,
7.3 GB, on a 64 GB M-series Mac. The log is `precompute.log` in that directory. The bundle now
served has 123,076 cells, so these numbers are an upper bound for it.

| Step | Time |
| --- | --- |
| count non-zeros per gene (pass 1) | 34.9 s |
| transpose to gene-major and statistics (pass 2) | 309.8 s |
| pick 2,000 HVGs and materialise the matrix | 0.6 s |
| PCA, 50 components | 67.7 s |
| UMAP, 161,432 x 50 | 160.2 s |
| write every bundle file | 0.6 s |

The server started in 0.6 s and held 106 MB resident while serving the 7.3 GB bundle, because
the arrays are memory-mapped. Warm responses on the `/fast` endpoints ranged from 3 ms for the
catalog to 29 ms for a gene overlay; nothing exceeded 30 ms.

## Limits

- The MarkerHeatmap loads Morpheus from `software.broadinstitute.org`; every other view runs offline.
- `--skip-expr` cannot build an embedding, because the per-gene variance is not stored; pass `--embedding-from`.
- A bundle added after start needs a server restart.
- The Differential tab reads the tables the differential workflow wrote; it recomputes no statistic.
- A `per_state` modality store never appears in Explore.
- Four helpers take `id.split("::")[0]` as the contrast name. For a four-part non-RNA id that yields the modality. The calls sit in `scalable_app.py` near lines 1413, 1473, 3634 and 3983 as of 2026-09-14.
- The chat validation covers one RNA-only bundle and no browser rendering.


### Shared GRN support for uploaded jobs

Regulatory networks, factor profiles, comparison matching, and regulatory question
routing are implemented in `cellHarmony/grn_analysis.py`. The viewer's
`grn_network.py` is a compatibility import. Uploaded jobs adapt their own h5ad and
completed differential files to the same functions; neither deployment needs the
LungMAP database to draw regulatory results. See
[GRN release validation](../../cellHarmony/webapp/GRN_RELEASE.md).

CombPlot preserves the stored observations without further aggregation by default. For a bundle explicitly containing metacells, `sv.observation_unit = "metacells"` labels those observations accurately; the general default is `cells`. The sampling selector is shared by individual-cell CombPlot and MarkerHeatmap.

CombPlot and MarkerHeatmap share **Cells per sample per cell type** controls in both panels: **5, 10 (default), 20, 50, All cells**. Sampling preserves individual observations and uses a reproducible cell-ID ranking within each sample × cell-type group. Groups smaller than the cap keep all available cells. The caption names the sample annotation; datasets without one are treated as one sample. Display filters apply before sampling. Donor means remain an explicit separate CombPlot mode.

MarkerHeatmap uses the same fixed marker rows for every sample size and reads the selected cells from the underlying expression store. Per-gene standardization uses the full source dataset, so values stay comparable when changing sample size or display filters. Screen, TSV, and PDF receive the same sample-size parameter and cell selection. Original analysis matrices are preserved.
# Unidentified metabolite descriptors

The shared scALABLE presentation layer supports audited CPTAC AML unidentified
metabolite descriptors (original ID, m/z, HILIC+/RP+ assay and retention time),
including source-table provenance on hover and labels in PDF exports. Declare
`annotation_source: "PDC000561"` in the metabolite modality manifest when using
that source model. Existing `Metabolite (AML)` modality labels are recognized.
Study-local Unknown IDs from unrelated assays are never mapped by ID alone.
