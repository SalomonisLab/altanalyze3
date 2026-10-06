# scALABLE-discover

scALABLE-discover is scALABLE-web with ICGS3 unsupervised clustering in place of reference
alignment. It keeps scALABLE-web's upload, QC, ambient RNA correction, Explore views and
Chat. It runs no differential expression and no modality imputation. Species: human or mouse.

Controls are described in [HOW_TO_USE.md](HOW_TO_USE.md).

## Method

| Step | Module | Parameters |
| --- | --- | --- |
| Load, ambient RNA, QC, normalize | `cellHarmony.cellHarmony_lite.combine_and_align_h5` with `cellharmony_ref=None` | the user's QC values; scALABLE-web's other arguments |
| Clustering | `clustering.ICGS.run_icgs3` | `ICGS3Config` defaults, except the five below |
| Marker networks | ICGS3's MarkerFinder heatmap call, `export_networks=True` | top 1000 markers per cluster |
| Cell communication | `fastComm.api.run_fastcomm` | scALABLE-web's values: CellChatDB, `min_cells` 5, score 0.2, 5 pairs |
| Serving bundle | `scalable_viewer.precompute`, at 10,000 cells or more | scALABLE-web's |

ICGS3 parameters that differ from `ICGS3Config`:

| Parameter | Value | Reason |
| --- | --- | --- |
| `min_genes`, `min_cells`, `min_counts` | 0 | QC ran once, upstream, with the user's values |
| `mito_percent` | None | same |
| `species` | `Hs` or `Mm` | the species the user chose |
| `louvain_downsample_cutoff` | 10,000 | step 1, Louvain community sampling, runs above 10,000 cells |
| `pre_pagerank_cells` | 10,000 | step 1 keeps 10,000 cells |
| `pagerank_cells` | 5,000 | step 2, PageRank, keeps 5,000 cells; these train NMF |
| `export_marker_networks` | True | the MarkerNetwork view reads these networks |
| `normalization` | `auto` (`none` when the upload holds no counts) | ICGS3 reads counts |
| `minimal_outputs`, `write_h5ad` | True, False under the default `minimal` exports | see Outputs |

ICGS3's own downsampling defaults are 30,000 (cutoff), 120,000 (step 1) and 30,000 (step 2).
The SVM still assigns every QC-retained cell, and UMAP and the heatmap cache hold every
clustered cell.

The QC-only call stops `combine_and_align_h5` before alignment. ICGS3 reads the
QC-retained counts. The combined h5ad holds every QC-retained gene, normalized by
`cellHarmony_lite.normalize_adata` as scALABLE-web does, for the cells ICGS3 placed in a
cluster. ICGS3's own `icgs3_result.h5ad` holds only its protein-coding gene set.

### Cell annotations and layers

The combined h5ad carries ICGS3's columns without the `ICGS3_` prefix: `cell_state_predicted`,
`cluster`, `original_NMF_cluster`, `SVM_score`, `SVM_margin`. `cell_state_predicted` is ICGS3's
GO-Elite BioMarkers label of each cluster, `<label>_c<n>`, one label per cluster; the pipeline
refuses a run where two clusters share a label or one cluster has two.

Two cell-state layers drive every view. A viewer picks one under `Cell-state layer`; the
choice rides in the `discover_layer` cookie and the server reads the job through it.

| Layer | Default | Cell states | MarkerHeatmap, MarkerNetwork, Chat markers | Cell communication |
| --- | --- | --- | --- | --- |
| Predicted cell states | yes | `cell_state_predicted` | copies of ICGS3's MarkerFinder tables under `outputs/layers/cell_state_predicted/`, cluster ids replaced by labels | fastComm run on `cell_state_predicted`, `outputs/fastComm/` |
| Clusters | no | `cluster` (C1..Cn) | ICGS3's own `ICGS3/MarkerFinder/` tables | fastComm run on `cluster`, `outputs/layers/cluster/fastComm/` |

The relabelled copies change only the cluster field and the centroid header; every marker,
fold, cell and order is ICGS3's. A viewer who switches layer rebuilds the serving caches.

### GO-Elite BioMarkers plot

`GO-Elite BioMarkers` draws, for one cell state, every BioMarkers term that overlaps its
markers, in the form of scALABLE's GO Terms view: GO-Elite z-score against FDR.

| Value | Source |
| --- | --- |
| p, FDR, overlap, query and term size | ICGS3's `GO-Elite/icgs3_biomarker_enrichment.tsv` (hypergeometric p; BH FDR over every cluster and term) |
| z-score | GO-Elite's `compute_z_score` (`goelite/structures.py`) on ICGS3's counts and the gene universe ICGS3 tested; ICGS3 stores no z-score |
| Highlight | FDR <= 0.05 and z > 2, the rule scALABLE's GO Terms view uses; BioMarkers sets have no ontology, so GO-Elite's pruning does not apply |
| Labelled terms | the term behind the cell-state label (blue) and the top highlighted terms |

`Download PDF` draws the same plot with scALABLE's GO Terms PDF renderer.

### Progress

While a job runs, the Run tab shows the percentage and the current stage: QC (ambient
correction, cell filters, normalization), then ICGS3 steps 1 to 9 (loading, normalizing,
downsampling, guide genes and NMF rank, NMF, MarkerFinder and SVM, BioMarkers, UMAP,
MarkerFinder set and networks), then the combined h5ad and fastComm. The runner reads the
line each step writes to the ICGS3 log.

## Changes to shared altanalyze3 code

| File | Change | Default behaviour |
| --- | --- | --- |
| `cellHarmony/cellHarmony_lite.py` | `cellharmony_ref=None` runs QC only and returns before alignment | unchanged: validation V1 |
| `clustering/ICGS.py` | `export_marker_networks`, `marker_network_top_n`, `marker_network_jobs`; CLI `--export-marker-networks` | off: validation V2 |
| `clustering/ICGS.py` | `minimal_outputs`; CLI `--minimal-outputs`: write only the final h5ad, the MarkerFinder marker set, the GO-Elite BioMarkers tables, `icgs3_config.json` and the log | off; not validated |
| `clustering/ICGS.py` | cell-state labels read `<label>_c18` for cluster C18, not `_cC18` | changes the label text of every ICGS3 run; not validated |
| `cellHarmony/cellHarmony_lite.py` | `unaligned_h5ad=<path>` writes the QC-passed cells below `min_alignment_score` | off unless a path is given; not validated |
| `cellHarmony/flask/pipeline.py`, `webapp/static/app.js` | scALABLE-web passes that path and offers "Download unaligned QC-passed cells (h5ad)" when any cell falls below the cutoff | new; not validated |
| `clustering/ICGS.py` | BioMarkers gene sets come from `clustering/biomarkers/`, or from `$ICGS3_BIOMARKER_DIR` | identical tables: validation V2 |
| `clustering/biomarkers/{Hs,Mm}/Ensembl-BioMarkers.txt.gz` | AltDatabase EnsMart72 BioMarkers files, gzipped unchanged | new |
| `cellHarmony/flask/tasks.py` | `JobRunner.WORKER_MODULE` names the isolated worker module | same module as before |
| `cellHarmony/webapp/requirements.docker.txt` | `python-louvain==0.16`, which ICGS3 downsampling imports | new |

Before this change ICGS3 read the BioMarkers files from one workstation's absolute path.
Any other host, a container included, produced `UNK-c<cluster>` labels without an error.
ICGS3 now logs the BioMarkers file it reads, or that it found none.

## Running it

Locally, from the repository root (the directory that holds `altanalyze3/`):

```bash
python3.11 -m uvicorn altanalyze3.components.cellHarmony.scalable_discover.app:app \
    --host 127.0.0.1 --port 8010
```

### Online deployment

scALABLE-discover runs in its own image and container, separate from scALABLE-web. Its
Dockerfile follows scALABLE-web's: the same base image, system libraries and Python
requirements list (`cellHarmony/webapp/requirements.docker.txt`), the same `/docs` health
check, and the same job runner with isolated workers.

| | scALABLE-web | scALABLE-discover |
| --- | --- | --- |
| Dockerfile | `cellHarmony/webapp/Dockerfile` | `cellHarmony/scalable_discover/Dockerfile` |
| compose file | `cellHarmony/webapp/docker-compose.scalable.yml` | `cellHarmony/scalable_discover/docker-compose.discover.yml` |
| image, container | `scalable-web` | `scalable-discover` |
| entry point | `cellHarmony.webapp.app:app` | `cellHarmony.scalable_discover.app:app` |
| host port | `127.0.0.1:8006` | `127.0.0.1:8007` |
| path prefix | `/scalable` (`CELLHARMONY_ROOT_PATH`) | `/scalable-discover` (`SCALABLE_DISCOVER_ROOT_PATH`) |
| job volume | `cellHarmony/webapp/jobs` -> `/srv/cellharmony/jobs` | `cellHarmony/scalable_discover/jobs` -> `/srv/scalable-discover/jobs` |
| API paths in `openapi.json` | see scALABLE-web's deploy guide | 41 |
| reference files | the registry, baked into the image | none; BioMarkers gene sets in `clustering/biomarkers/` |

On the host:

```bash
git clone --depth 1 https://github.com/SalomonisLab/altanalyze3.git
cd altanalyze3/altanalyze3/components/cellHarmony/scalable_discover
docker network create lungmap_default 2>/dev/null || true   # the compose file joins it, as scALABLE-web's does
docker compose -f docker-compose.discover.yml up -d --build
curl -s http://127.0.0.1:8007/openapi.json | jq '.paths | length'   # 41
```

Then route the proxy, as `/scalable/` is routed, passing the prefix through:

```apache
ProxyPass        /scalable-discover/ http://127.0.0.1:8007/scalable-discover/
ProxyPassReverse /scalable-discover/ http://127.0.0.1:8007/scalable-discover/
```

and check `curl -s https://<site origin>/scalable-discover/openapi.json | jq '.paths | length'`.

Bind `./jobs` to real disk; uploads and results land there. `altanalyze3/components/.dockerignore`
excludes that directory, so jobs never enter a rebuilt image. Job retention is scALABLE-web's:
each new upload deletes finished jobs older than 8 hours, and `cellHarmony/webapp/cleanup_jobs.py
--job-root <jobs dir>` purges on a schedule. Chat needs `CELLHARMONY_ASSISTANT_URL` to reach the
LungMAP.net assistant, which the compose file sets for the `lungmap_default` network.

| Variable | Default | Effect |
| --- | --- | --- |
| `SCALABLE_DISCOVER_ROOT_PATH` | empty | FastAPI `root_path` behind a proxy |
| `SCALABLE_DISCOVER_JOB_STORAGE` | `scalable_discover/jobs` | job folders; kept apart from scALABLE-web's |
| `SCALABLE_DISCOVER_EXPORTS` | `minimal` | `minimal` or `full`; see Outputs |
| `ICGS3_BIOMARKER_DIR` | unset | a directory of `<Hs|Mm>/Ensembl-BioMarkers.txt[.gz]` that replaces the bundled sets |
| `CELLHARMONY_*` | as scALABLE-web | workers, memory limits, caches, bundle threshold, chat assistant URL |

## Outputs

One folder per job, `<job storage>/<job id>/`, with `uploads/`, `logs/pipeline.log` and
`outputs/`:

`SCALABLE_DISCOVER_EXPORTS` selects the export set. `minimal`, the default, exports only the
final h5ad and scALABLE's MarkerFinder output set, and skips intermediates. Every value is
still computed; only the files differ.

| Path under `outputs/` | `minimal` | `full` | Content |
| --- | --- | --- | --- |
| `combined_with_umap_and_markers.h5ad` | yes, download | yes, download | the final h5ad: every QC gene, normalized, clustered cells, ICGS3 labels, UMAP |
| `ICGS3/MarkerFinder/` markers, redundant markers, centroids, fold-matrix h5ad, networks | yes | yes | scALABLE's MarkerFinder set; `icgs3_marker_genes.zip` is the download |
| `ICGS3/logs/`, `ICGS3/icgs3_config.json`, `scalable_discover_parameters.json` | yes | yes, parameters also a download | provenance |
| `icgs3_umap_coordinates.tsv`, `fastComm/` tables | yes | yes | read by the Explore and Cell communication views; not downloads |
| `bundle/` | at 10,000 cells or more | same | disk-backed serving store |
| ICGS3 `GO-Elite/` BioMarkers tables | yes | yes | read by the GO-Elite BioMarkers plot |
| `layers/` relabelled MarkerFinder copies and the second fastComm run | yes | yes | read by the alternative cell-state layer |
| ICGS3 `sNMF/` and `UMAPs/` tables, heatmap PDF/SVG, UMAP PDFs, `icgs3_result.h5ad` | no | yes | intermediates (ICGS3 `--minimal-outputs` skips them) |
| ICGS3 input counts and ambient raw counts | temporary, removed after use | `ICGS3_input/`, `ambient/` | staged inputs |
| `icgs3_cluster_assignments.txt`, `icgs3_results.zip`, `cell_communication_fastcomm.zip` | no | yes, downloads | exports |

`ambient/soupx_summary.tsv` comes from the shared QC code and appears in both modes when
ambient correction runs. Under `minimal` the ICGS3 folders `sNMF/` and `UMAPs/` exist but
stay empty.

## Validation (2026-10-05)

Scripts and results are in `validation/`. Each script takes its input paths as arguments
and records inputs by file name and SHA-256.

| Check | Input | Result |
| --- | --- | --- |
| V1: `cellHarmony_lite` reference path unchanged | 2 human lung 10x files, 10,107 cells, ambient on | assignments file and every matrix hash identical, 7,445 aligned cells |
| V1: QC-only mode equals scALABLE QC | same | 8,485 cells x 32,738 genes; X, `counts`, `soupx_raw` identical on all 7,445 aligned cells |
| V2: ICGS3 defaults unchanged | 8,485 QC cells | 9 of 9 output files byte-identical, GO-Elite tables included |
| V2: network option changes no result | same | 9 of 9 files identical; 24 network tables for 25 clusters |
| V3: human, HTTP end to end | lung, ambient on | 50 of 50 checks; 8,099 of 8,485 cells (95.5%) in 25 clusters; 140 s from run request to completion |
| V3: mouse, HTTP end to end | 2 mouse 10x files, 26,726 cells, ambient off | 50 of 50 checks; 25,261 of 26,048 cells (97.0%) in 11 clusters; bundle built; 300 s |
| V4: Docker requirement set | clean Python 3.11 env from `requirements.docker.txt` (scanpy 1.11.5, numpy 2.4.6) | 50 of 50 checks; 8,098 of 8,099 cells share the V3 cluster label (ARI 0.99967) |

UI check (`validation/ui_window_switch_check.py`, headless Chrome, 2026-10-06, one human job
of 2,999 clustered cells): 10 of 10 Explore plot types draw the same plot after Windows 2 -> 1,
with no new error. A pathway opened from Chat leaves no markup after the panel switches to
Cell frequency, and the frequency axis gives 237 px to a 38-character state name. Result:
`validation/ui_window_switch_check.json`.

Existing tests: 66 of 66 cellHarmony and runner tests pass. `clustering/test_ICGS.py` passes
32 of 37; the same 5 tests fail on the unmodified ICGS.py.

V1 to V4 ran before the 10,000/5,000 downsampling values, the `minimal` export mode, the
cell-state layers, the GO-Elite BioMarkers plot, the stage progress, the `_c18` label fix and
the unaligned-cell export. None of these has been run; all are not verified.

Also not verified: a Docker build and run (the local container VM did not start).

## Where the time goes

Measured before the downsampling and export changes, on the first V3 human job and the V3
mouse job, macOS, one isolated worker each, by
`validation/stage_timings.py <job dir>`. Each stage runs from one log marker to the next,
so the stages sum to the total.

| Stage | Human, 8,485 cells (s) | Mouse, 26,048 cells (s) |
| --- | ---: | ---: |
| scALABLE load, ambient RNA, QC, ICGS3 input write | 4.5 | 8.9 |
| ICGS3 load, normalize, gene filter | 1.6 | 5.2 |
| ICGS3 feature selection, NMF rank estimate | 14.2 | 51.0 |
| ICGS3 NMF | 8.0 | 63.0 |
| ICGS3 MarkerFinder, SVM, MarkerFinder | 2.0 | 8.0 |
| ICGS3 GO-Elite BioMarkers | 4.9 | 1.8 |
| ICGS3 UMAP fit and plots | 43.9 | 70.8 |
| ICGS3 marker heatmap; static render of every cell alone | 47.1; 44.3 | 56.0; 50.4 |
| ICGS3 h5ad write (gzip) | 3.8 | 17.5 |
| Combined h5ad, fastComm, downloads, bundle | 4.2 | 12.3 |
| Total | 134.2 | 294.5 |

Neither job exceeded 30,000 cells, so PageRank downsampling selected every cell (0.0 s and
0.1 s) and NMF trained on all of them.

## Limits

- ICGS3 leaves cells unassigned when their best SVM score is 0 or lower: 386 of 8,485 and
  787 of 26,048 in the V3 runs. The status line reports the count; those cells appear in no view.
- UMAP coordinates differ between library versions; cluster labels agreed on 8,098 of 8,099 cells.
- The Morpheus MarkerHeatmap uses scALABLE-web's colour scheme, not the Spectral standing rule.
