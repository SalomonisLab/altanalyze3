# Discover UMAP separation comparison — 2026-10-06

The user reported overlapping populations in the PCA-based Accelerated layout, requested
centroid separation as a selection criterion, and explicitly required at least 200 fitting
cells per cluster, every cell from smaller clusters, and sampling within larger clusters.
Only the final embedding and its landmark selection change. The upstream ICGS3 procedure,
expression, annotations, samples, feature roster, MarkerFinder and cluster assignments remain
unchanged. No labels supervise the UMAP objective and no clusters are moved after fitting.

## Identical-input comparison

Source: job `3989b35d2f9547deb5d208895df0f67f`, corresponding to the reported screenshot.
The complete final dataset contains 60,164 assigned cells, 24,857 expression features and
69 clusters. Every candidate uses all 1,273 recorded final MarkerFinder UMAP features and
returns coordinates for every assigned cell. The original workflow retained 65,662 cells
through QC and assigned 60,164; this comparison makes no additional cell exclusions.

The earlier feature-based job `6b5c827ff6234dd5b53890753e339e74` has identical ordered cells,
genes, clusters, UMAP feature panel, and every stored expression value and sparse index.
Inputs and complete ordered rosters are hash-gated before fitting. All landmark candidates
use the same within-cluster selection: minimum 200 where available, every smaller-cluster
cell, and proportional allocation of remaining places to reach 30,000. The budget expands
if minimum cluster coverage requires more cells. Random seed is 0, transformation blocks
are at most 50,000 cells, and each fit runs in a separate local Python process with four
OpenBLAS/OMP/Numba threads. Full fits include every cell without landmark selection.

| Representation / fit | Neighbors | min_dist | Seconds | Peak RSS GiB | Median separation | Lower-tail separation | Macro neighbor consistency | Input-neighbor recall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **Marker features / landmarks — selected** | **15** | **0.75** | **87.63** | **2.646** | **0.4328** | **0.1558** | **0.7407** | **0.2044** |
| Marker features / landmarks | 15 | 0.15 | 88.29 | 2.689 | 0.3925 | 0.1303 | 0.7594 | 0.2110 |
| Marker features / landmarks | 50 | 0.75 | 184.47 | 2.749 | 0.3489 | 0.1535 | 0.7313 | 0.1921 |
| Marker features / full fit | 50 | 0.75 | 254.08 | 3.402 | 0.3883 | 0.0979 | 0.7166 | 0.2000 |
| 50 PCs / landmarks | 15 | 0.75 | 45.47 | 2.075 | 0.3899 | 0.1371 | 0.6499 | 0.1540 |
| 50 PCs / landmarks | 15 | 0.15 | 48.89 | 2.033 | 0.2839 | 0.0788 | 0.6861 | 0.1741 |
| 50 PCs / full fit | 15 | 0.75 | 50.87 | 2.234 | 0.3333 | 0.1472 | 0.6323 | 0.1353 |
| 50 PCs / landmarks, repulsion 2 | 15 | 0.75 | 42.87 | 2.054 | 0.3225 | 0.0934 | 0.6528 | 0.1557 |

Marker-feature UMAP uses correlation distance. PCA is centered, unscaled Scanpy ARPACK
PCA with Euclidean UMAP. Repulsion is 1 except where explicitly stated. These are single
paired runs on a local Mac, not a replicated throughput or container-memory study.
Elapsed time covers embedding input conversion and fitting/transformation, excluding
verified-input loading and output saving. Peak RSS covers the isolated embedding process;
it is not the memory peak of the complete multimodal workflow.

## Selection criteria and limits

For each cluster, calculate its centroid and RMS cell-to-centroid radius. For every other
cluster calculate centroid distance divided by the sum of both radii; take the minimum
ratio for that cluster. Report the median and tenth percentile across **all 69 clusters**.
This dimensionless separation cannot improve merely by stretching coordinates. The selected
candidate maximizes both statistics among these eight settings. The screenshot's saved PCA
layout scores 0.3242 and 0.1454 respectively; the selected candidate's median is 33.5% higher.
The saved screenshot used the previous sampling policy, so the table's 200-cell PCA control
provides the stronger comparison for representation changes.

Neighbor consistency uses up to 50 fixed evaluation cells per cluster and 15 two-dimensional
neighbors, averaged within each cluster and then across clusters. Input-neighbor recall uses
345 fixed evaluation cells (up to five per cluster) against exact 15-neighbor correlation
distances in the original marker space. Rare and overlapping states are retained in every
diagnostic. The compact feature variant has slightly higher neighbor consistency and recall;
the selected variant has higher centroid separation and lower-tail separation at similar speed.

The selected feature method is 2.90 times faster and uses 22% less isolated-process peak RSS
than the full feature fit. It uses more memory and time than the PCA candidates. These
descriptive metrics measure layout and consistency with existing clusters; they do not prove
biological accuracy, disease specificity, or calibrated distances between biological states.
Related populations can remain connected.

## Implementation and verification

Discover's saved `landmark` choice now passes feature landmark mode and 15 UMAP neighbors
to the original ICGS3 hook. Its **All cells** choice and ICGS3's default remain the original
full correlation-distance feature fit with 50 neighbors and min_dist 0.75. Landmark selection
and provenance live in `clustering/umap_fit.py`; the minimum is 200 and the strategy is
`cluster_stratified_proportional`. Core clustering parameters are unchanged.

`benchmark_umap_separation.py` enforces verified input hashes and recorded authorization.
`compile_comparison` checks every fitting roster against the stated selection rule and
evaluates all candidates against the same expression-neighbor reference. Full reports,
versions, per-cluster values and original source hashes are preserved in
[`discover_umap_separation_results.json`](discover_umap_separation_results.json).

```bash
PYTHONPATH=/path/to/altanalyze3 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4 \
python -m altanalyze3.components.clustering.benchmarking.benchmark_umap_separation \
  /path/to/verified/input /path/to/new/output --neighbors 15
```

Rebuild and restart the Discover image to apply this profile. No dependency, environment
variable or data migration is added. Existing completed jobs retain their saved coordinates.

The actual `ICGS.compute_umap_outputs` hook was rerun on a separate copy of the reported
job, with the same profile, and completed in 88.19 seconds. Its coordinates are bit-identical
to the selected benchmark. Every stored X value and sparse index, ordered cell and gene,
obs/var field, cluster assignment and full marker panel passed comparison against source;
the original job's file hash is unchanged. The preview is
`964772417ed44d5084c7483cbb5b13f1`. It reruns only UMAP and retains the source analysis's
other artifacts and original full-analysis timing. The live Explore endpoint returns all
60,164 barcodes and exact revised coordinates at stored float32 precision without omissions
(CSV parsing introduces less than 4e-15 float64 roundoff), and SFTPC
violin values equal those of the original job. The browser rendered the revised UMAP.
The focused suite passed 31 checks; an additional 150 mixed-size allocations checked
minimum coverage, exact expanded budgets, no duplicates and complete rare-state coverage.
