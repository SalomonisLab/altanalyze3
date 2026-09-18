## Pathway color scales and window switching (2026-09-18)

Cross-modality pathway diagrams now have a numeric color bar per modality:
MarkerFinder r for cell-state markers, or a separate symmetric log2FC scale for
each modality in contrasts. Domains use the mapped features of the selected
pathway. A node stripe uses its modality's strongest absolute retained score;
all original feature scores remain in hover text. Modalities without mapped
values are explicitly labelled, without an invented numeric range. The same
vector color bars and ranges are included in PDF exports.

Switching between one and two windows resizes integrated views without drawing
cached expression data over them. Pathway, state, zoom, and the other panel's
view are preserved. SVG gradient IDs are unique to each mounted panel.

## Named TF activity by cell type (2026-09-18)

“Where is SPI1 TF-activity most enriched?” ranks the named TF's stored mean
activity across cell states, names the leading states in the answer, and supplies
all states in a sortable/filterable table and matching bar chart. The default
no longer averages away the cell-state axis. An explicit top-N limit is honored.
Questions naming a cell state retain the existing within-state factor ranking.
This is an activity-level ranking, without a new differential or enrichment test.
The shared analysis and chart apply to both uploaded sessions and the viewer.

## Cross-modality pathway representation (2026-09-18)

Chat accepts “What pathways have the best cross-modality HSC-1 representation?”
for retained positive MarkerFinder scores, without a differential run. Add
“for contrasts” to use the selected completed comparison. Explicit named
comparisons must match a saved comparison label. All joined modality runs must
share the comparison, grouping field, cell-state field, and comparison type.

The sortable, searchable table and stacked bar chart report unique original
feature IDs per pathway per modality. Repeated diagram nodes do not increase
counts. Checked modality boxes are AND requirements, each requiring at least
one mapped hit. Filters and pagination share one result set and do not change
normalization. Chart heights represent raw counts, not ranking scores.

Ranking uses **modality count + balanced coverage**. For each modality, divide
the pathway count by its maximum count across all evaluated diagrams; balanced
coverage is the mean of those fractions over modalities with mapped hits.
Modality breadth therefore ranks first, and each mapped modality contributes
equally to the secondary score regardless of assay size. This is descriptive
representation across the 68 bundled WikiPathways diagrams, not an enrichment
p-value or an exhaustive pathway database search.

RNA/TF activity map by exact gene or Ensembl ID; ADTs use curated gene mappings;
lipid species map through Discover lipid classes; named metabolites match exact
names. Unknown metabolites do not acquire inferred identities. GRN edges count
once only when all regulators and the target belong to the pathway. Missing
marker tables, unmatched comparisons, and unmapped modalities are reported.
No additional significance filter is applied to retained differential calls.

Selecting a table pathway or bar opens its diagram in Explore (markers) or
Differential (contrasts), preserving state, comparison, and pathway ID. Node
stripes show the contributing modalities, with individual feature IDs and
retained statistics on hover. Effect sizes from different modalities are never
averaged into one color. Both the stacked chart and diagram offer vector PDF
export through the shared exporter. Implemented in the shared web app and used
by the precomputed viewer as well.

# scALABLE (cellHarmony web)

scALABLE is the analysis web application of AltAnalyze3. It aligns single-cell RNA data to a
reference atlas, then runs marker discovery, approximate UMAP placement, cell-cell
communication scoring, optional multimodal imputation, and group differential analysis on the
aligned cells. A Chat tab reads a question, routes it to one of the analyses, and answers from
the job's own data. The stack is FastAPI plus vanilla JavaScript. No R runs anywhere in it.

Location:
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/`

This file is the technical reference: how to run the server, what each pipeline step computes,
what each statistic means, and where each output lands. The user guide is
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/HOW_TO_USE.md`.

scALABLE-viewer is a second program that serves a precomputed bundle through this same
application, without a Run tab. Its documentation lives in
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/scalable_viewer/`.

Every line count and default in this file comes from the source as of 2026-09-14
(`app.py` 6,904 lines, `static/app.js` 7,593 lines, `flask/pipeline.py`, `cellHarmony_differential.py`).

## Documentation map

| Question | File |
| --- | --- |
| How do I use the interface? | `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/HOW_TO_USE.md` |
| How does the alignment CLI work? | `/Users/saljh8/Documents/GitHub/altanalyze3/docs/cellHarmony.md` |
| How does the differential CLI work? | `/Users/saljh8/Documents/GitHub/altanalyze3/docs/cellHarmony_differential.md` |
| What did the 2026-08-25 Explore and Chat repairs change? | `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/dev/README.md` |
| How does the approximate UMAP endpoint work? | `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/approximate_umap_api.md` |
| How do I serve a precomputed bundle? | `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/scalable_viewer/README.md` |
| Which imputation model produces which modality? | the `README.md` of each `rna2*` component, listed under Modalities below |

## Run the server

From the repository root:

```bash
cd /Users/saljh8/Documents/GitHub/altanalyze3
/opt/homebrew/opt/python@3.11/bin/python3.11 -m uvicorn \
  altanalyze3.components.cellHarmony.webapp.app:app --host 127.0.0.1 --port 8000
```

Open `http://127.0.0.1:8000`. The ASGI target exists because `app.py` ends with
`app = create_app()`. `dev/restart_app.sh` restarts this exact command and logs to
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/server_8000.log`.

The Chat tab needs a second process, the LungMAP assistant, on port 8001. Without it the chat
route answers HTTP 503. See Chat below.

### Environment variables

`config.py` reads these at start. The code holds every other parameter as a constant.

| Variable | Default | Effect |
| --- | --- | --- |
| `CELLHARMONY_ROOT_PATH` | empty | FastAPI `root_path` for a reverse proxy |
| `CELLHARMONY_JOB_STORAGE` | `<webapp>/jobs` | where job folders live |
| `CELLHARMONY_REFERENCE_REGISTRY` | `<webapp>/../flask/reference_config.json` | the reference registry |
| `CELLHARMONY_MAX_FILES` | `7` | files per job |
| `CELLHARMONY_JOB_WORKERS` | `2` | pipeline threads |
| `CELLHARMONY_EXPORT_APPROX_PDFS` | `false` | write approximate UMAP comparison PDFs |
| `CELLHARMONY_H5AD_COMPRESSION` | `lzf` | `lzf`, `gzip` or `none` for every h5ad the pipeline writes |
| `CELLHARMONY_ASSISTANT_URL` | `http://127.0.0.1:8001/api/assistant/viewer-intent` | the chat intent router |
| `CELLHARMONY_ENABLE_FASTCNV` | `false` | run fastCNV clone analysis after fastComm |
| `CELLHARMONY_FASTCNV_EXPORT_PDF` | `true` | write the fastCNV clone PDF |

Boolean variables accept `1 true t yes y on` and `0 false f no n off`. Any other value keeps
the default. The server accepts at most 1 GiB per request and only the extensions `.h5` and
`.h5ad`.

### Docker

`Dockerfile` starts from `python:3.11-slim`, installs `requirements.docker.txt`, copies the
`altanalyze3/components` tree to `/app/altanalyze3/components`, exposes port 8000 and runs
uvicorn on `0.0.0.0:8000`. The build context is `altanalyze3/components`:

```bash
cd /Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components
docker build -f cellHarmony/webapp/Dockerfile -t cellharmony-web:latest .
docker compose -f cellHarmony/webapp/docker-compose.yml up
```

`docker-compose.yml` defines one service, `cellharmony-web`, maps port 8000, bind-mounts
`./jobs` to `/srv/cellharmony/jobs`, and sets the registry path inside the image.

### Job storage and retention

Each job is one folder under the job storage directory, named by a 32-character hex id, holding
`job.json`, `uploads/`, `outputs/` and `logs/pipeline.log`. Two mechanisms remove old jobs:

| Mechanism | Rule |
| --- | --- |
| on every new upload | `job_manager.py` deletes jobs older than 8.0 hours whose status is `completed`, `failed` or `cancelled` |
| `cleanup_jobs.py` | a scheduled purge, see the flags below |

```bash
/opt/homebrew/opt/python@3.11/bin/python3.11 \
  /Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/cleanup_jobs.py \
  --job-root <dir> --retain-days 7 --keep-latest 0 --dry-run
```

| Flag | Default | Meaning |
| --- | --- | --- |
| `--job-root` | `CELLHARMONY_JOB_STORAGE` | directory of job folders |
| `--retain-days` | `7` | keep jobs updated within N days |
| `--keep-latest` | `0` | always keep the newest N purge-eligible jobs |
| `--dry-run` | off | report without deleting |

The script never removes a job whose status is not terminal or whose differential status is
`queued` or `processing`. Exit code 1 means a missing job root or a negative flag. `launchd/org.cellharmony.jobs-retention.plist`
schedules the script at 03:15 daily with `--retain-days 7 --keep-latest 5`. Its three paths
contain the placeholder `your-user` and need editing before `launchctl load`.

## The reference registry

`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/flask/reference_config.json`
lists every species and reference the Run tab offers. Relative paths resolve against the
registry's own directory.

| Key | Required | Meaning |
| --- | --- | --- |
| `id`, `label` | yes | the menu value and the menu text |
| `states_tsv` | yes | cellHarmony centroid matrix; its file stem becomes the query cluster key |
| `reference_clusters_tsv` | yes | barcode to cell state, for the reference UMAP |
| `reference_coords_tsv` | yes | barcode to UMAP1, UMAP2 |
| `cluster_key` | yes | obs column of the reference labels |
| `impute_modalities` | no | list of modality ids the Impute modality menu may offer |
| `impute_config.<id>.bundle_path` | no | model bundle for that modality |
| `impute_config.<id>.expression_scale` | no | `linear`, `log1p` or `log2`; the scale of the imputed matrix |
| `impute_config.<id>.log_base` | no | `e` or `2`, used to derive a linear `counts` layer |
| `impute_config.<id>.pseudobulk_statistic` | no | GRN only, see Modalities |
| `ambient_options` | no | present in the file but no Run control sends it |

References registered on 2026-09-14, with the modalities each one can impute:

| Species | Reference id | Label | Imputable modalities |
| --- | --- | --- | --- |
| Human | `hs_lung_cellref2_reference` | CellRef2.0 | `adt`, `lipids`, `grn` |
| Human | `hs_lung_cellref_reference` | LungMAP CellRef v1.1 (Guo 2023) | `adt`, `lipids`, `grn` |
| Human | `hs_bm_reference` | Bone marrow CITE-Seq (Zhang 2024) | `adt`, `metabolite`, `lipid`, `grn` |
| Human | `hs_lung_hlca_reference` | Lung HLCA (Sikkema 2023) | `adt`, `lipids`, `grn` |
| Human | `hs_lung_natri_reference` | Lung ILD Atlas (Natri 2024) | `adt`, `lipids`, `grn` |
| Human | `hs_lung_bpd_sun_reference` | BPD Atlas (Sun) | `adt`, `lipids`, `grn` |
| Mouse | `mm_lung_cellref_reference` | LungMAP CellRef v1.0 | none |
| Mouse | `mm_bm_reference` | Bone marrow CITE-Seq (Ferchen 2025) | `adt` |
| Mouse | `mm_lung_adult_airway_reference` | LungMAP Adult Lung/Airway v1.0 | none |
| Mouse | `mm_lung_flu_reference` | Lung Regeneration-Infection (Niethamer 2025) | none |

To add a reference, export its three TSVs from an annotated h5ad, then add an entry:

```bash
PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 \
/opt/homebrew/opt/python@3.11/bin/python3.11 \
  /Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/export_reference_umap_bundle.py \
  --h5ad <reference.h5ad> --cluster-key <obs column> --outdir <dir> \
  [--output-prefix <stem>] [--umap-key X_umap] [--downsample <n cells>]
```

The script writes `<prefix>_reference_umap.tsv`, `<prefix>_reference_clusters.tsv`,
`<prefix>_reference_metadata.json` and a `<prefix>_reference_config_snippet.json` stub. The
stub carries a placeholder `states_tsv` path and a `soupx_options` key that the registry does
not read. Edit both before pasting the entry.

## The pipeline

`flask/pipeline.py` `run_cellharmony_pipeline` runs once per job on a worker thread. The status
endpoint reports the stage from the log tail. The table gives every step with its compiled
parameters.

| Progress | Step | Function and parameters | Output under `outputs/` |
| --- | --- | --- | --- |
| 10 | load reference | `_lookup_reference`; default gene = first row of `states_tsv` | log line `Reference metadata loaded.` |
| 10 | QC and alignment | `cellHarmony_lite.combine_and_align_h5`: `alignment_mode="cosine"`, `metacell_align=False`, `unsupervised_cluster=False`, QC values from the form, `ambient_correct_cutoff="auto"` when the user chose Yes | `cellHarmony_lite_assignments.txt` |
| 72 | markers | `generate_marker_heatmap_from_adata`: `marker_method="markerfinder"`, `top_n=50`, `cells_per_cluster=100`, `seed=0`, `export_networks=True`, `network_top_n=1000`, `write_heatmap_cache=True` | `marker_heatmap/cell_state_marker_heatmap.pdf`, heatmap cache, marker TSVs |
| 78 | marker networks | NetPerspective per cell state; zip of `.pdf .tsv .png .npz` | `cell_state_marker_genes.zip` |
| 82 | approximate UMAP | `approximate_umap`: `umap_key="X_umap"`, `jitter=0.05`, `num_reference_cells=1` | `approximate/`, UMAP columns appended to the assignments |
| 88 | imputation | one call per selected modality, see Modalities | `combined_with_umap_and_markers_<id>.h5ad` and companions |
| 88 | write h5ad | `combined_with_umap_and_markers.h5ad` with `CELLHARMONY_H5AD_COMPRESSION` | the Explore data source for RNA |
| 90 | cell communication | fastComm, see below | `fastComm/` |
| 92 | fastCNV | only when `CELLHARMONY_ENABLE_FASTCNV` is on and species is human or mouse | `fastCNV/fastcnv*` |
| 100 | finish | differential options, modality list, message `Approximate UMAP completed.` | `job.json` |

Alignment scores each cell against every reference centroid by cosine similarity and keeps the
best state. Cells below the cosine cutoff leave the dataset before every later step.

The QC form defaults are `min_genes 500`, `min_counts 1000`, `min_cells 0`, `mit_percent 15`,
`align_cutoff 0.4`. The pipeline's own fallbacks differ (`200`, `500`, `3`, `10`, `0.1`), but
the form always sends its values, so the fallbacks only apply to a job created outside the
form.

### Optional marker outputs

The Run tab's **Optional marker heatmap exports** controls apply to RNA and every
imputed modality. scALABLE skips the static PDF and SVG render by default, changed on
2026-09-18. The default keeps marker statistics, state centroids, the compact
selected-marker cache, interactive heatmaps, and on-demand Download PDF. Only the
packaged static files change, so the browser loses no function. On a 147,376-cell study
the five marker stages took 188.10 seconds with the static render and 38.05 seconds
without it, a 4.94x speedup (`MARKER_OUTPUT_BENCHMARK.md`).

To write the static files, set **Static heatmap files → Generate static PDF/SVG**. That
option writes a PDF plus SVG at 2400 DPI (or `MARKER_HEATMAP_DPI`) and displays up to
100 cells per state. `marker_write_svg` and `marker_heatmap_dpi` apply only when the
static render runs.

scALABLE already disables both marker expression-matrix TSV exports. The pipeline keeps
the primary expression H5AD for Explore, Chat, differential analysis and session reload.
The `marker_heatmap_h5ad` package API and CLI keep their own defaults; only scALABLE
changed.

The QC API accepts `marker_render_heatmap` (default `false`), `marker_write_svg`
(`true`), `marker_heatmap_dpi` (`null`, use renderer default) and
`marker_cells_per_cluster` (`100`, `0` displays all). Values are saved with the job.
Each modality records output options and stage timings. Sampling affects
only display; marker scoring and cell-state averages still use all input cells.

For standalone Python callers, `generate_marker_heatmap_from_adata` accepts
`render_heatmap=False`, `write_svg=False`, and `heatmap_dpi=600`. To avoid constructing
any plotting matrix, also set `write_heatmap_cache=False`, `write_heatmap_tsv=False`,
and `write_expression_tsv=False`. Marker tables and marker-by-state centroids remain.
The CLI equivalents are `--skip-heatmap-render`, `--skip-svg`, `--dpi 600`, and
`--markers-only` (the last combines all heatmap-output skips). Existing CLI/API defaults
remain unchanged. The standalone markers-only mode does not supply an interactive
heatmap cache; use the web option when interactive exploration is required.

`tests/benchmark_marker_outputs.py JOB_JSON NEW_OUTPUT_DIRECTORY` benchmarks these
options against saved matrices, checks identical marker tables and centroids, and
writes stage timings and artifact sizes without changing the saved job or imputations.
Measured results for the 147,376-cell marrow study are in [MARKER_OUTPUT_BENCHMARK.md](MARKER_OUTPUT_BENCHMARK.md).

### fastComm parameters

`FastCommParams`: `lr_sources=("CellChatDB",)`, `min_cells=5`, `min_lr_expression_score=0.2`,
`max_lr_candidates_per_state_pair=5`, `include_self_edges=False`, species passed only for human
or mouse. The exemplar report lists interactions with `fastcomm_score >= 0.25`, and the job
records `significance_threshold = 0.25`. When obs holds one of `Library`, `group`, `sample`
or `Donor`, the pipeline also scores each sample separately into `fastComm/per_sample/`, which
the per-sample plot and the cell-communication differential read. A fastComm failure never fails
the job; the job records `enabled: false, status: "failed"`.

Files: `fastComm/fastcomm_scores.tsv`, `state_pair_summary.tsv`, `state_expression.tsv`,
`significant_interactions.tsv`, `significant_interactions.md`, `per_sample/split_scores_long.tsv`,
`cell_communication_fastcomm.zip`.

## Modalities

Each modality holds one feature matrix over the same cells. RNA is the measured one. The others are
imputed from the aligned RNA by models that ship with AltAnalyze3. `_DEFAULT_MODALITY_DEFINITIONS`
in `app.py` and `_MODALITY_DEFINITIONS` in `flask/pipeline.py` hold the same table.

| id | Menu label | One feature is a | Marker heatmap | Marker network | Differential network and GO |
| --- | --- | --- | --- | --- | --- |
| `rna` | RNA | gene | yes | yes | yes |
| `lipids` | Lipids | lipid | yes | no | no |
| `adt` | ADT (CITE-seq) | ADT | yes | no | no |
| `metabolite` | Metabolite (AML) | metabolite | yes | no | no |
| `lipid` | Lipid (AML) | lipid | yes | no | no |
| `grn` | GRN (edges) | edge | no | no | no |
| `grn_tf` | TF activity (imputed) | factor | yes | no | no |
| `cell_communication` | Cell communication | ligand-receptor interaction | no | no | no |

`cell_communication` is not imputed. It appears only as a Differential modality, once a
fastComm run exists. Aliases such as `cite`, `regulon`, `rna2grn` and `fastcomm` normalise to
these ids in `_normalize_modality_id`.

### How a job selects modalities

The Run form field `impute_modality` offers `none`, `all`, and the ids the reference declares.
The browser sends `impute_modalities` as a list of zero or one id. The `/qc` route rejects an
id the reference does not declare with HTTP 400. `_selected_impute_modalities` expands `all`
to the reference's full list and intersects any other request with it.

### The models

| id | Model | Training data | Input | Output scale | Bundle |
| --- | --- | --- | --- | --- | --- |
| `adt` | one ElasticNet per ADT (`alpha=0.01`, `l1_ratio=0.5`) over a curated gene panel, z-scored in and out | human bone marrow CITE-seq: 3,132 genes in, 129 ADTs out; human lung: 81 genes in, 56 ADTs out, 40,000 training cells; mouse bone marrow: 172 genes in, 103 ADTs out | per cell | `log1p`, base e; negative predictions clipped to 0 and counted | `rna2adt/rna2adt_bm_bundle.pkl`, `rna2adt/lung/rna2adt_hs_lung_bundle.pkl`, `rna2adt/mouse/rna2adt_mm_bundle.pkl` |
| `lipids` | one ElasticNetCV per lipid over that lipid's top correlated genes; 202 lipids, 1,303 genes | 50 sorted cell-type profiles from 10 human lung donors, 0 bulk profiles | per cell | `log2` | `rna2lipid/rna2lipid_hs_lung_lipidwise_bundle.pkl` |
| `metabolite` | one ridge per metabolite (`l1_ratio=0`, 1,000 genes per target, alpha 100); 2,533 metabolites | CPTAC AML, 84 cases with RNA and metabolomics | pseudobulk = sample x cell state | `log2` | `rna2metabolite/artifacts/rna2metabolite_aml_bundle.pkl.gz` |
| `lipid` | one ridge per lipid (120 genes per target, alpha 100); 1,009 lipids | CPTAC AML, 87 cases with RNA and lipidomics | pseudobulk = sample x cell state | `log2` | `rna2lipid/aml/artifacts/rna2lipid_aml_bundle.pkl.gz` |
| `grn` | one ridge per TF-to-target edge on three standardised features: target expression, TF expression, regulon mean | leukemia reference 7,486 edges; lung hybrid reference 63,647 edges, 284 TFs, 39 pseudobulks | pseudobulk = sample x cell state | edge score, then a per-cell TF activity | `rna2grn/rna2grn_bundle.pkl.gz`, `rna2grn/rna2grn_lung_hybrid_bundle.pkl.gz` |

Each model's own `README.md` and `VALIDATION.md` under
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2adt/lung/`,
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2lipid/`,
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2lipid/aml/`,
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2metabolite/` and
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2grn/` give the held-out
performance. Two limits from those files matter for interpretation. The lung lipid model
reports Pearson 0.871 on a holdout that shares donors with training, so that number is an
upper bound. The AML metabolite model has a held-out median Spearman of 0.267 over 2,533
metabolites; 1,084 of 2,533 exceed 0.3, and the imputable subset is the intended deliverable.

The pseudobulk models group cells by `sample|cell_state`, using the first obs column found among
`sample`, `Sample`, `library`, `Library`, `orig.ident`, `dataset`, `Dataset`. GRN reads the
reference's `pseudobulk_statistic`: `mean_over_cells_of_log1p_cp10k` averages the normalised
matrix, which the lung references declare; otherwise the model sums `layers['counts']`.

### GRN: edges and TF activity

Selecting GRN imputation produces separate edge and factor outputs. Factor activity is
the **sum of predicted outgoing edge scores**, matching the precomputed viewer. Cell
predictions are computed in bounded chunks and discarded after summing by factor.
Differentials use independent sample x cell-state predictions; cells are never counted
as independent biological replicates for either GRN modality.

| Object | Features | Rows | Used by |
| --- | --- | --- | --- |
| `combined_with_umap_and_markers_grn_tf.h5ad` | factors | cells | TF Explore plots and differential detail violins |
| `combined_with_umap_and_markers_grn_tf_pseudobulk.h5ad` | factors | sample x cell state | TF activity differential |
| `combined_with_umap_and_markers_grn_edges.h5ad` | `TF|target` edges | sample x cell state | edge differential and GRN network |

Older jobs used standardized target-set enrichment for per-cell TF exploration. They
are labelled **TF enrichment (legacy)**, retain their edge differentials, and must be
reprocessed to obtain predicted TF activity comparisons. No existing scientific output
is rewritten automatically.

MarkerFinder operates directly on the continuous predicted TF scores; RNA depth-scaling
validation does not apply to them. When no TF markers pass selection, the job completes
with the activity and differential outputs and explains why no marker heatmap exists.

### Output files per modality

| id | Explore h5ad | Differential h5ad | Marker outputs | Zip |
| --- | --- | --- | --- | --- |
| `lipids` | `combined_with_umap_and_markers_lipids.h5ad` | none | `marker_heatmap_lipids/`, `cell_state_lipid_markers.zip` | `lipids_results.zip` |
| `adt` | `combined_with_umap_and_markers_adt.h5ad` | none | `marker_heatmap_adt/`, `cell_state_adt_markers.zip` | `adt_results.zip` |
| `metabolite` | `..._metabolite.h5ad` | `..._metabolite_pseudobulk.h5ad` | `marker_heatmap_metabolite/` | `metabolite_results.zip` |
| `lipid` | `..._lipid.h5ad` | `..._lipid_pseudobulk.h5ad` | `marker_heatmap_lipid/` | `lipid_results.zip` |
| `grn` | network only | `..._grn_edges.h5ad` | none | `grn_results.zip` |
| `grn_tf` | `..._grn_tf.h5ad` | `..._grn_tf_pseudobulk.h5ad` | `marker_heatmap_grn_tf/` | `grn_tf_results.zip` |

Explore opens one h5ad per modality; it never reads layers of the RNA object. Each imputed
h5ad also carries `layers['counts']`, the linear-scale values derived from `expression_scale`
and `log_base`. Marker heatmaps of imputed modalities use `top_n=5` per cell state. The pipeline also
writes `obsm['X_<id>']` and `uns['imputed_modalities']` into the RNA h5ad for bookkeeping.

Two defects in that bookkeeping exist on 2026-09-14. The `lipids` and `adt` branches assign
`uns['imputed_modalities']` instead of adding to it, so a job with both keeps only the `adt`
entry. The `lipids` and `lipid` branches both write `uns['lipid_feature_names']`. No shipped
reference declares both, so the second collision cannot occur yet.

## The differential engine

### Request

`POST /api/jobs/{job_id}/differential` takes `modality` (default `rna`), `population_col`,
`sample_field`, `group1_samples`, `group2_samples`, and `comparison_type` (`cells` or
`pseudobulk`). The request carries no threshold: the server sets thresholds from the modality.
The engine tests every cell state in `population_col`. The UI offers `pseudobulk` only for a single
h5ad upload or four or more files; otherwise the server downgrades it to `cells`.

### Grouping and aggregation

The pipeline subsets cells to the two groups and writes a synthetic covariate column with a case
label and a control label. Under `pseudobulk`, `compute_pseudobulk_per_population` sums
`layers['counts']` over the cells of each `population|sample` group with at least 10 cells,
normalises to counts per 10,000, and takes `log2(x + 1)`. Any cell state with fewer than 2
pseudobulks in either condition leaves the object before testing. The pseudobulk modalities
(`metabolite`, `lipid`, `grn`, `grn_tf`) arrive already aggregated and skip this step, because re-summing
predicted values would distort them.

### Tests

| Comparison | Test | Where |
| --- | --- | --- |
| `pseudobulk` | moderated t-test with empirical Bayes variance shrinkage: pooled variance `s2`, prior `median(s2)`, `s2_shrunk = 0.8*s2 + 0.2*median(s2)`, two-sided t with `n_case + n_ctrl - 2` degrees of freedom | `cellHarmony_differential._moderated_t_test` |
| `cells` | `scanpy.tl.rank_genes_groups` with `method="wilcoxon"` on CP10k, log1p values | `cellHarmony_differential._rank_genes_scanpy` |

The `method` parameter applies only to the `cells` branch. The moderated t-test returns an empty
table when either arm has fewer than 2 replicates.

Independent filtering restricts the Benjamini-Hochberg correction to genes with a value above
0.1 in at least 2 rows. A gene outside that filter keeps its fold and raw p and receives
FDR 1.0. The reported `log2fc` is not the test statistic: the engine de-logs the matrix, takes
each arm's mean in linear space, and reports `log2((mean_case + 1) / (mean_ctrl + 1))`. The
engine calls a gene when `|log2fc| > log2(fc_thresh)` and the chosen p falls below `alpha`.

### Thresholds per modality

`_differential_runtime_params` in `flask/pipeline.py`:

| Modality | alpha | fold threshold | min cells per group | p used for the call |
| --- | --- | --- | --- | --- |
| `rna` and default | 0.05 | 1.2 | 20 for `cells`, 1 for `pseudobulk` | raw p for `pseudobulk`, FDR for `cells` |
| `lipids`, `adt` | 0.05 | 1.0 | 10 for `cells`, 1 for `pseudobulk` | raw p for `pseudobulk`, FDR for `cells` |
| `metabolite`, `lipid` | 0.05 | 1.2 | 2 | raw p |
| `grn`, `grn_tf` | 0.05 | 1.1 | 2 | raw p |

Under `pseudobulk` the gate on a cell state is `min_replicates_per_group = 2`, not the table's
min-cells column. Under `cells` with fewer than 200 cells in the contrast, the gate rises to
`max(4, min_cells)`. A cell state that fails the gate appears in `DEG_summary_*` with 0 DEGs
and 0 tested genes and contributes no rows to `DEG_detailed_*`, so the Differential tab does
not list it. When no cell state yields a DEG, the pipeline promotes the pooled result into the
detailed table under the label `Pooled overall`.

### GRN edge differentials

The `grn` differential tests edges, not TFs: its input is the `TF|target` pseudobulk object,
with fold threshold 1.1 and raw p. The Explore mode `GRN edges` reads the same object and draws
edges touching a requested gene whose mean score over the selected pseudobulks exceeds a
threshold, up to `max_edges` (default 300).
It defaults to the first cell state in lineage order and a TF with the largest summed
absolute outgoing score in that state (or retains the panel's selected TF). Sample
`(any)` averages samples within that state. Neutral nodes and edges do not imply
up/down regulation; TFs are diamonds, targets are circles, and width encodes score.

Explore's `Regulatory network` uses stored RNA markers for the selected state and
state-specific predicted GRN edges. It does not require differential runs. Its
filters are marker fold versus other states, absolute edge score, and mean stored
RNA expression of the TF within the state. Unmeasured or below-threshold TFs are
omitted. Node colours use actual marker log2 folds; missing marker folds are gray.
The Differentials regulatory view retains its matched RNA/edge differential gates.
Cell communication is available in either Explore panel regardless of its modality.

### Cell communication differentials

`_run_cell_communication_differential` is a separate engine. It pivots
`per_sample/split_scores_long.tsv` into an interaction by sample matrix keyed
`sender|receiver|ligand|receptor`, taking the maximum score per cell. With 2 or more samples
per arm it runs a two-sided Mann-Whitney U per interaction and a local Benjamini-Hochberg. With
one sample in an arm it sets p to 1.0 and reports an effect-rank pseudo-p, `rank(|delta|)/N`;
the manifest records this. The fold is `log2((case_mean + 1e-4) / (control_mean + 1e-4))`. The
summary marks an interaction significant at `fdr <= 0.10` and `|delta| >= 0.05`. The receiver
state fills the `population` column. No network and no GO run for this modality.

### GO-Elite

GO-Elite runs for `rna` only, into `GeneSetEnrichment/`, three times: the assigned-group gene
lists, the up-regulated genes per cell state, and the down-regulated genes per cell state. It
matches gene symbols, never Ensembl ids, so the h5ad needs `var['gene_symbols']`, `var['features']`
or symbol `var_names`. The background is the tested gene list; terms keep 5 to 2,000 genes. The
test is hypergeometric with a z-score and Benjamini-Hochberg per query set. Prioritisation keeps
a term as `selected` when `z >= 1.96`, `FDR <= 0.1`, overlap `>= 3`, and no parent or child term
represents it better. The Differential tab reads only the assigned-group file `goelite_tsv`.

### Differential networks

NetPerspective builds one network per cell state with 2 or more distinct DEGs, from
`Hs_Ensembl-TF-BioGRID-Pathway.txt` or the mouse equivalent under
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/interactions/`.
It keeps an edge only when both genes are DEGs, drops self edges, and excludes BioGRID edges by
default. Node colour follows the fold sign. The PDF and the summary TSV colour an exact zero
grey; the interactive view colours it as up. The step needs `python-igraph`.

### Output files per contrast

`DEG_detailed_<tag>.tsv`, `DEG_summary_<tag>.tsv`, `DEG_assigned_groups_<tag>.tsv`,
`DEG_pooled_overall_<tag>.tsv`, `DEG_coreg_pooled_<tag>.tsv`, `DEG_fold_matrix_<tag>.tsv`
(the complete gene by cell-state fold matrix the heatmap reads), `heatmap_<tag>_by_<col>.pdf`
and `.tsv`, `differentials_only_<tag>.h5ad`, `GeneSetEnrichment/GOElite_<tag>.tsv` and `.pdf`,
`interaction-plots/<tag>/<state>.pdf .png _interactions.tsv`, cell-frequency plots,
`manifest.json`, and a zip of the run. The tag is `<group1>_vs_<group2>` with sample names
joined by `__`.

## Chat

### Architecture

Neither the web app nor the viewer contains a language model. `POST /api/jobs/{job_id}/chat`
sends the question, the cell-state names, the completed contrast label, the groupable covariate
names, and the modality ids to the intent router at `CELLHARMONY_ASSISTANT_URL`. The router
runs inside the LungMAP site process at
`/Users/saljh8/Dropbox/LungMAP/refactored_website/app/lungmap/assistant/viewer.py`.
Three passes read each question:

| Pass | Method | Model called |
| --- | --- | --- |
| 1 | deterministic slot extraction of genes, cell states, contrasts, covariates | no |
| 2 | scoring against a weighted cue table of 17 protocols | no |
| 3 | only when passes 1 and 2 return `clarify` or `unsupported` and the sentence names something in the dataset | yes |

The model in pass 3 is Qwen2.5-1.5B-Instruct quantised to Q4_K_M, 1.07 GB, run on CPU through
`llama_cpp` in a worker subprocess, with `temperature 0.0`, `max_tokens 160`, a JSON schema, and
a 6-second cap. The prompt contains the protocol menu and the dataset's names only. The model
never sees a value and never states a number. The web app then computes every figure and
statistic from the job's h5ad. When the router does not answer after two attempts, the route
returns HTTP 503.

### Protocols and executors

The router returns one of 17 protocol ids. `_CHAT_PROTOCOL_ALIAS` maps them onto the web app's
executors.

| Protocol | Executor | Answer |
| --- | --- | --- |
| `cell_identity` | markers | 25 marker rows of one state (`gene, cluster, fold, p, source`) and a DotPlot of the top 12 |
| `state_comparison` | compare | 30 marker rows for two states and a DotPlot |
| `expression_lookup` | expression | top 5 states by mean for each gene (`gene, cell state, mean, fraction`) and a DotPlot |
| `state_contrast`, `contrast_specificity`, `shared_vs_state_specific`, `patient_stratification`, `donor_heterogeneity`, `most_affected_state` | differential | top 25 rows of the completed contrast in the named state, sorted by FDR, and a volcano |
| `regulatory_driver` | shared regulatory network | inline network with factor activity and expression statistics from matching completed comparisons |
| `tf_activity`, `regulator_activity` | shared factor profile | inline ranked activity bars and a table of reported fold changes/FDRs |
| `pathway_program`, `communication_rewiring` | existing view | names the Differential or Explore tab that holds the analysis |
| `severity_gradient`, `dose_response`, `composition_shift`, `coexpression_module`, `annotation_concordance` | not implemented here | names the missing statistic instead of answering with a neighbouring analysis |

The web app repairs the router's reading against its own names before dispatch: it matches
cell states longest-first and keeps a token as a gene only when the dataset holds that gene.
Statuses in the response: `read_from_question`, `clarify`, `not_implemented`, `use_existing_view`,
`not_found`, `not_run`, `not_available`, `not_covered`. A `not_covered` answer names the states
the contrast does cover and states that missing data is not absence of change.

The `plot` field is a specification object such as `{"kind": "dotplot", "genes": [...]}`. The
browser re-calls the same `/dotplot` or `/combplot` endpoint the Explore tab uses, or draws a
volcano from the returned table. The chat cannot produce a figure the Explore tab cannot.

`GET /api/jobs/{job_id}/chat-examples` builds up to 6 example questions from the job: the two
largest states, two marker genes of the largest state, and two contrast questions when a
differential has completed.

The viewer's chat answers more protocols from a bundle's precomputed statistics; see
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/scalable_viewer/README.md`.
The measured routing accuracy, 55 of 55 and 68 of 68 questions on the COPD bundle, is in
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/visualization/scalable_viewer/VALIDATION.md`.
No equivalent suite exists for the web app's own chat route.

## Explore payload definitions

| Payload | Definition |
| --- | --- |
| DotPlot `mean` | arithmetic mean of the modality's matrix over the cells of a group |
| DotPlot `frac` | fraction of those cells with a value strictly greater than 0 |
| DotPlot default genes | one gene per group with the largest (group mean minus mean of all other cells), up to 12 groups |
| CombPlot value | for each gene, the mean of the matrix over one donor's cells within one group; a donor with fewer than `min_cells` (default 5) cells in a group is dropped and counted in `n_groups_dropped` |
| CombPlot donor column | `sample_field`, else the first of `meta_sample`, `donor`, `Donor`, `Library`, `sample`, `pool` |
| Group by / Show | any obs column with 2 to 60 categorical levels |
| Colour by | the same rule |
| X and Y axes | the two cellHarmony UMAP axes, every 2-D `obsm` entry, and every float obs column with at least one non-integral value; whole-number columns are counts or indices and are excluded |
| Display filter fields | non-numeric obs columns with 2 to 120 levels, excluding UMAP columns; the browser sends one value per row |
| Cell frequency | fraction of each sample's cells per state, over the filtered cells |
| Violin | per-state values over the filtered cells, top states by mean, 10 by default and 30 with one panel |
| State colours | `uns['<cluster_key>_colors']` when its length matches, else the 12 Paired colours by position |

The DotPlot draws marker size `4 + 18*frac` px and colour by mean on a white to red ramp. The
CombPlot draws a colour strip of states, one bar row per gene, and, in the viewer only,
covariate annotation bands. Explore PDFs come from the browser through svg2pdf, so they hold
vector shapes and editable text. Server PDFs use matplotlib with `pdf.fonttype 42` and
DejaVu Sans.

## API endpoints

| Method and path | Purpose |
| --- | --- |
| `GET /api/meta/species` | the registry |
| `GET /api/meta/reference-preview?species=&reference=` | reference UMAP points and median label positions |
| `POST /api/jobs` | upload `species`, `reference`, `sample_names[]`, `files[]` |
| `POST /api/jobs/{id}/qc` | QC settings and `impute_modalities` |
| `POST /api/jobs/{id}/configure` | change species or reference; clears outputs |
| `POST /api/jobs/{id}/run` | queue the pipeline |
| `GET /api/jobs/{id}/status` | job.json plus log head and tail and `differential_ui` |
| `GET /api/jobs/{id}/umap` | cells with coordinates for the UMAP modes |
| `GET /api/jobs/{id}/expression` | one feature over the cells, UMAP and violin payloads |
| `GET /api/jobs/{id}/genes?modality=` | feature names of one modality |
| `GET /api/jobs/{id}/dotplot` | mean and fraction per gene per group |
| `GET /api/jobs/{id}/combplot` | per-donor means per gene per group |
| `GET /api/jobs/{id}/plot-variables?modality=` | groupable columns, colour columns, coordinate sets, numeric axes |
| `GET /api/jobs/{id}/display-filters?modality=` | filter fields and their levels |
| `GET /api/jobs/{id}/marker/network?population=` | NetPerspective marker network of one state |
| `GET /api/jobs/{id}/grn/network?genes=&cell_state=&sample=&threshold=&max_edges=` | GRN edges touching the genes |
| `GET /api/jobs/{id}/fastcomm/plot?population=&plot_type=&limit=` | seven cell-communication plot types |
| `GET /api/jobs/{id}/fastcomm/network` | incoming or outgoing pairs; no browser caller |
| `GET /api/jobs/{id}/marker/heatmap.tsv`, `.pdf`, `/jobs/{id}/marker/heatmap/viewer` | MarkerFinder matrix, its PDF, the Morpheus page |
| `GET /api/jobs/{id}/umap/pdf`, `/expression/pdf` | server-rendered PDFs; the Explore buttons use the browser path instead |
| `GET /api/jobs/{id}/download/{artifact}`, `/log` | artifacts and the pipeline log |
| `POST /api/jobs/{id}/differential`, `GET .../differential/status` | run and poll a contrast |
| `GET /api/jobs/{id}/differential/interactive/{summary,heatmap,volcano,go,network,table,gene}` | the Differential views |
| `GET /api/jobs/{id}/differential/interactive/pdf?mode=`, `.../gene/pdf` | server PDFs of the views |
| `GET /api/jobs/{id}/differential/archive`, `/artifact/{key}`, `/heatmap?format=`, `/network/{id}?format=` | downloads |
| `GET /api/jobs/{id}/chat-examples`, `POST /api/jobs/{id}/chat` | Chat |
| `POST /api/jobs/{id}/client-log` | append a browser message to the job log |
| `POST /api/tools/approximate-umap` | run approximate UMAP on files already on disk |

The Morpheus heatmap page loads its scripts from `software.broadinstitute.org`, so the
MarkerHeatmap mode needs network access. Every other mode works offline.

## Known gaps on 2026-09-14

- The Run form never sends `ambient_option`, so the registry's `ambient_options` lists have no effect.
- The browser sends `subset2_by` for the second annotation row. The web app's `/dotplot` and `/combplot` accept only `subset_by`, so only the viewer applies that row to those plots.
- The CombPlot annotation bands render from `tracks` fields that only the viewer's builder emits.
- `/api/jobs/{id}/fastcomm/network`, `/umap/pdf` and `/expression/pdf` have no browser caller.
- The two `uns['imputed_modalities']` defects listed under Modalities.
- No automated test covers the web app's chat route.


### GRN release validation

See [GRN_RELEASE.md](GRN_RELEASE.md) for the shared implementation, release checks,
compatibility notes, and the evaluated boundary of upload/viewer parity. Completed
uploaded comparisons are retained and selectable in **Differential → Completed
comparison**. Regulatory Chat questions route locally in both deployments; other
questions still use the configured intent service.

CombPlot defaults to individual cells in both visualization panels. **Display → Donor means** opts into per-donor, per-group averaging and reveals **Min cells**. Both annotation filters apply to the cell view. Bar and annotation-band details appear only on hover, including in the shared viewer; PDFs omit those hover labels. **Cells per sample per cell type** controls individual-cell CombPlot and MarkerHeatmap.

CombPlot and MarkerHeatmap share **Cells per sample per cell type** controls in both panels: **5, 10 (default), 20, 50, All cells**. Sampling preserves individual observations and uses a reproducible cell-ID ranking within each sample × cell-type group. Groups smaller than the cap keep all available cells. The caption names the sample annotation; datasets without one are treated as one sample. Display filters apply before sampling. Donor means remain an explicit separate CombPlot mode.

MarkerHeatmap uses the same fixed marker rows for every sample size and reads the selected cells from the underlying expression store. Per-gene standardization uses the full source dataset, so values stay comparable when changing sample size or display filters. Screen, TSV, and PDF receive the same sample-size parameter and cell selection. Original analysis matrices are preserved.
# Unidentified AML metabolite labels

For the CPTAC AML imputation model, unidentified features retain their source ID
and display m/z, assay/polarity and retention time, for example
`Unknown 0235 · m/z 234.09706 · RP+ · RT 1.221 min`. Plot hovers also identify
PDC000561, the original supplementary-table row, DOI and any recorded formula.
The same labels are used in heatmaps and PDF exports; queries and plot clicks
continue to use the original IDs. These descriptors support tracing and candidate
comparison, not chemical identification: retention times depend on the assay.
Other metabolomics studies do not inherit this study's Unknown-ID annotations.
Viewer manifests can explicitly identify this catalog with
`annotation_source: "PDC000561"` on their metabolite modality.
