# scALABLE web: serving large datasets without loading them into RAM

Status update, 2026-10-02: steps A/B are implemented; step C now runs each heavy
analysis in a disposable subprocess supervised by the shared thread queue.
Serving caches have byte/entry/TTL limits, and the default bundle threshold is
10,000 cells. Full all-modality validation at 128,388 aligned cells, repeated
saved-job opens, concurrent visitors and Chrome view checks are documented in
`README.md` and `dev/memory_workflow_results.json`. The original design below is
retained as the proposal and records its original threshold and validation goals.

Design, 2026-09-28. Implemented: the fast expression-store builder of step A (section 3,
"Build time, measured"). Not implemented: the step A call in `flask/pipeline.py`, step B
and step C. Covers the scALABLE upload app
(`cellHarmony/webapp`, `cellHarmony/flask`) and reuses the scALABLE viewer
(`visualization/scalable_viewer`). Target: a job of more than 100,000 cells finishes and
serves Explore, Differential and Chat without RAM growing with cell count, and the
extra build step takes no more than 2 minutes at 100,000 cells.

## 1. Where the RAM goes today

Measured on job `cfad37d18ece4729bed20e2f80bda9b7`: 7 human marrow h5 files, 42,796
cells after QC, every marrow imputation model.

| # | Cause | Code | Effect |
|---|---|---|---|
| 1 | Jobs run as threads inside the web server | `flask/tasks.py:54` `ThreadPoolExecutor` | The pipeline peak (10.6 GB) lives in the server process. Python keeps it: the server held 8.5 GB after the job ended. Two workers put two jobs in one process. |
| 2 | Explore reads each result h5ad whole | `webapp/app.py:1536` `ad.read_h5ad(h5ad_path)` in `_get_expression_cache`, once per modality | The RNA file holds `X` and `layers['counts']`, 124,430,079 non-zeros each, plus four `obsm` modality copies. |
| 3 | The GRN edge view reads a dense per-cell matrix whole | `webapp/app.py:2433` in `_grn_edges_adata` | 42,796 × 7,486 float32 = 1.28 GB. |
| 4 | Differential gene detail reads its h5ad whole | `webapp/app.py:1295`, `:1306` in `_open_gene_detail_adata` | One more full copy per open. |

Causes 2 to 4 grow linearly with cells. By extrapolation from the 42,796-cell job (not
measured at these sizes), the first Explore request needs about 4.7 GB for RNA at
100,000 cells and about 47 GB at 1,000,000 cells, before GRN edges add 3 GB and 30 GB.
Cause 1 adds the pipeline peak on top and never returns it.

## 2. What the viewer already solves

The viewer serves an 83,416-cell, 5-modality atlas at 2.1 to 2.3 GB RSS (measured). It
never opens an h5ad while serving:

- `data_api.Dataset` and `ModalityStore` read a gene-major CSC store
  (`<prefix>_expr_data.npy`, `_expr_indices.npy`, `_expr_indptr.npy`) with
  `np.load(..., mmap_mode="r")`. One gene is one contiguous slice
  (`data_api.py:359`). The same loader memory-maps per-state `stats_mean` and `stats_frac`.
- `bundle_meta.BundleJobStore` (`bundle_meta.py:810`) and `seed_expression_cache`
  (`bundle_meta.py:709`) hand the unchanged web routes a `BundleAnnData` over those
  memmaps. The same `app.py` code answers both apps.

The viewer does **not** keep expression in SQLite. SQLite holds only tabular results:
differential calls (`integrated_pseudobulk/differentials.sqlite`) and the study record.
Expression in SQLite would need one row per non-zero (290 million at 100,000 cells) and
an index build far beyond 2 minutes, or one blob per gene, which is the CSC store again.

## 3. Design: every large job becomes a bundle

| Step | Where | What |
|---|---|---|
| A. Build | the job, at the end of `flask/pipeline.py`, while the matrix is still in memory | Write `outputs/bundle/` in the viewer's `BundlePaths` format: the RNA CSC store, per-state `stats_mean`/`stats_frac`/`stats_n`, `cells.npz` (barcodes, cell state, sample, covariates, UMAP), and one feature-major store per imputed modality. |
| B. Serve | `webapp/app.py` | For a bundled job, `_get_expression_cache`, `_grn_edges_adata` and `_open_gene_detail_adata` answer from the bundle through the existing `BundleAnnData` path. The server never reads a result h5ad whole. |
| C. Isolate | `flask/tasks.py` | Run each job in a spawned child process (`ProcessPoolExecutor(mp_context=spawn, max_tasks_per_child=1)`), so the pipeline's memory returns to the operating system when the job ends. `CELLHARMONY_JOB_WORKERS` keeps its meaning. |

The h5ad outputs stay on disk as downloads and as input to later differential runs.

### Why the build fits in 2 minutes

The viewer's builder re-reads a compressed h5ad twice and sorts each block in Python.
Its own log for 123,076 cells, 749,284,851 non-zeros:

| Stage | Seconds |
|---|---:|
| pass 1, count per gene | 25.4 |
| pass 2, 16 blocks, sort and scatter | 175.9 |
| expression store total | 201.1 |
| HVG, PCA, UMAP | about 104 |
| whole build | 305.2 |

A web job needs neither re-read nor embedding. It already has the matrix in memory and a
UMAP from alignment. So step A:

1. Transposes the in-memory CSR with SciPy's C routine `csr_tocsc`, writing straight into
   `np.memmap` output files. One O(non-zeros) pass, no Python sort, no second copy in
   anonymous RAM.
2. Computes per-state sums and non-zero counts as sparse products with a one-hot state
   matrix.
3. Writes modality stores with the viewer's existing `ingest_modality` layout.

Projected cost at 100,000 cells (about 290 million non-zeros at this job's density of
0.089), not measured: under 60 seconds. Even the viewer's slower two-pass algorithm,
scaled linearly by non-zeros, gives about 78 seconds for the expression store.

### Build time, measured

`visualization/scalable_viewer/fast_store.py` builds the store, and `precompute.py` uses it
by default. It reads each 8,192-row block once, transposes it with `csr_tocsc` in 8 worker
processes, spills it, then fills each range of genes in one contiguous write. Measured on
2026-09-28 through `precompute.py`, every bundle file byte-identical to the original
builder's:

| Dataset | Cells | Non-zeros | Original builder (s) | New builder (s) |
|---|---:|---:|---:|---:|
| human marrow job | 42,796 | 124,430,079 | 21.6 | 3.7 |
| COPD metacells | 83,416 | 666,852,999 | 146.9 | 13.4 |
| COVID lung | 302,922 | 1,216,355,183 | 247.5 | 24.3 |

`fast_store.build_expression_store_from_csr` is the step A entry point for a matrix
already in memory: 2.2 s on the marrow job, byte-identical output. Logs:
`/Users/saljh8/Dropbox/LungMAP/LungMAP.net/Datasets/COPD-atlas/scalable_viewer/checks/fast_builder_20260928/`.

### Threshold

Apply the bundle path to jobs above `CELLHARMONY_BUNDLE_MIN_CELLS`, default 100,000.
Two serving paths need two sets of tests. Once the parity gate in section 5 passes,
setting the threshold to 0 leaves one path.

## 4. Local and hosted

- One format and one reader on both. The bundle is plain `.npy` and `.npz` files in the job
  folder: no database server and no build at server start. A restarted server opens the
  memmaps and seeds its caches in seconds, as the viewer does today.
- The job volume must be a local block device (local SSD or EBS). Memory-mapping over NFS
  or EFS turns every gene lookup into network reads.
- In a container, file-backed pages count toward the cgroup limit but are reclaimable.
  Anonymous RSS stays near the viewer's 2 to 3 GB. `docker/compose.yaml` needs no change
  beyond the job volume.

## 5. Validation gates

| Gate | Pass condition |
|---|---|
| Parity | A job under 100,000 cells served both ways returns identical values from every Explore, Differential, marker and Chat endpoint, compared on at least 50 genes per modality, sample order and identity included. |
| Build time | 100,000 cells: step A at most 120 s. Measured on a real dataset, not a subset. |
| Scale | `/Users/saljh8/Dropbox/Transfer/COVID_combined_with_umap_and_markers.h5ad` (302,922 cells) completes. Server RSS after the job returns to within 1 GB of its idle value. |
| Memory bound | Server RSS while serving the 302,922-cell job stays below 4 GB. |

## 6. Not covered

- The pipeline's own peak during alignment and imputation (10.6 GB at 42,796 cells).
  Step C returns it after each job but does not lower it. Lowering it means chunked
  imputation, which needs its own measurement first.
- cellHarmony-differential reads the h5ad when a comparison runs. It runs in the job
  process, so step C bounds it; moving it onto the bundle is a later step.


Analysis retention now uses streamed prediction H5ADs, lossless sample/state
profile indices for broadcast modalities, row-block MarkerFinder statistics and
bounded differential reads. GRN X/counts hard links avoid duplicate buffers;
small enough GRN inputs use shared RAM and larger inputs use an anonymous
uncompressed memory map. The map is temporary and is released by the OS when a
worker exits. A job's combined RNA file omits duplicate wide imputation `obsm`
arrays above 20,000 cells; separate modality artifacts remain portable H5ADs.

The 128k-cell full pipeline and million-row statistics measurements are in the
README and `dev/memory_workflow_results.json`. These establish bounded prediction
and statistical processing, not a guarantee for a full million-cell upload.
Initial sparse RNA import and normalization still materialize the RNA matrix.
