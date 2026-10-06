# Discover PCA UMAP comparison — 2026-10-06

**Historical comparison:** the PCA choice below was superseded after the user reported
overlapping populations and required at least 200 landmarks per cluster. See
`UMAP_SEPARATION_BENCHMARK.md` for the subsequent identical-input, eight-setting
comparison and current feature-based Discover default. These original measurements
remain preserved; they were taken under the earlier sampling policy.

The user authorized a PCA-based Accelerated embedding. This comparison changes only
the final UMAP representation and graph; it preserves the original ICGS3 clustering,
MarkerFinder feature panel, cell order, samples and expression values. ICGS3's default
and Discover's **All cells** option retain the original full feature UMAP.

## Paired embedding benchmark

Actual input: Discover job `2c4231425f4846b99c1ec752c2e41fbf`, uploaded as 65,662 cells
with 24,857 genes. All three methods consume the same 62,498 assigned cells and the
exact recorded 1,260 final MarkerFinder features. The original ICGS3 procedure leaves
3,164 cells unassigned; this comparison makes no additional exclusions.

Each method runs in a fresh Python process, with four OpenBLAS/OMP/Numba threads.
Elapsed time includes imports and input conversion inside the embedding call, but
excludes loading the verified input files and saving results. Peak RSS covers the
entire isolated process. These are local Mac measurements, not container cgroup peaks.

| Method | Seconds | Peak RSS (GiB) | 15-neighbor state consistency |
| --- | ---: | ---: | ---: |
| Previous Accelerated: original features, correlation, 50 neighbors, 30k landmarks | 204.48 | 2.395 | 0.8625 |
| Centered Scanpy PCA, 50 PCs, Euclidean, 15 neighbors, full fit | 59.34 | 2.057 | 0.8718 |
| **Selected: same PCA, 30k landmarks and remaining-cell transform** | **45.24** | **1.988** | **0.8679** |

The selected method is 4.52 times faster in this embedding comparison, with 17% lower
isolated-process peak RSS. Every cell receives finite coordinates, and all 47 existing
states are represented in the fitting set. Neighbor consistency measures agreement
with existing ICGS3 clusters in 2D; it does not establish biological accuracy or prove
that either embedding preserves every expression-space relationship.

PCA is centered, unscaled Scanpy ARPACK PCA on the complete ordered marker panel.
Both candidates use `min_dist=0.75`, seed 0 and up to 50 PCs. Landmark selection retains
up to 50 cells from each state, then fills the 30,000 budget randomly; it expands the
budget when needed to cover every state. Remaining cells map in balanced blocks of
at most 50,000. Jobs below the fitting budget fit all cells using the same PCA method.

## Provenance and reproduction

`benchmark_pca_umap.py` requires a verified baseline manifest, full feature/cell/state
identity hashes and explicit PCA benchmark authorization. It calls the production
helper `accelerated_umap.pca_umap`. The input was prepared from the saved expression,
the recorded ordered UMAP feature panel and original clusters; no panel was reconstructed.
The three results have identical input, feature, cell and cluster hashes.

Results, versions, identity hashes and per-state diagnostics:
[`discover_pca_umap_results.json`](discover_pca_umap_results.json).

```bash
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4 \
python -m altanalyze3.components.clustering.benchmarking.benchmark_pca_umap \
  /path/to/verified/input /path/to/new/output --method landmark
```

Use `--method legacy` and `--method scanpy` for the paired alternatives. Each output
directory must be new. The input contract is documented and enforced by
`benchmark_landmark_umap.load_verified_input`.

## Complete Discover rerun

A separate verification job, `481c86ad85984ed2af792de660fc42d7`, reused the original
uploaded H5AD and QC settings through the local server and its isolated worker.
It completed in 304.65 seconds, compared with 321.73 seconds for the earlier run.
UMAP output fell from 189.51 to 75.67 seconds in these full runs. Other steps varied:
NMF/MarkerFinder/SVM took 98.9 seconds in the rerun versus 44.5 seconds earlier.
These two full runs are not a controlled throughput comparison; the isolated paired
embedding measurements above provide the stronger timing comparison.

The rerun's Ensembl lookup produces BioMarkers predictions for 41 of 47 clusters.
Six clusters without marker evidence use `c7`, `c11`, `c22`, `c34`, `c35`, `c44`.
All 24,857 supplied symbol aliases resolve from the uploaded `feature_name` column,
while native Ensembl IDs remain unchanged. The existing fastComm procedure now runs
using those aliases. Its output still reports 15 unmatched catalog genes; this is
recorded coverage, not a claim of complete ligand/receptor coverage.

Exact comparison passed for all ordered cells and genes, source obs/var fields,
every expression value, cluster assignments, original NMF clusters, SVM scores and
margins. The ordered UMAP feature list, canonical MarkerFinder marker/redundant
tables and centroid matrix are byte-identical. Only the authorized embedding,
predictions and symbol alias field change.

The saved-job browser test shows lowercase labels and confirms that filtering `c18`
still requests the original category. Automated checks passed: 61 in the test
environment and 14 identifier/UMAP-input checks in the live server's Python runtime.
The shared serving alias fix additionally passed 20 bundle/lookup/plot checks and
five targeted checks in the live Python runtime. On the real job, `SFTPC` and
`ENSG00000168484` return identical violin values. Every primary ID remains accepted,
with 24,857 readable gene suggestions rather than twice as many browser options.

The full-run memory sampler stopped partway through because of a macOS psutil
`proc_cmdline` error. Consequently, no complete full-run memory peak or reduction
is claimed. Isolated embedding peak RSS above was measured independently with
`getrusage`. Full checks and coverage are retained in
[`discover_pca_workflow_validation.json`](discover_pca_workflow_validation.json).
