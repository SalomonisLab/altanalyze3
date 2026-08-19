# scTriangulate in AltAnalyze3

scTriangulate reconciles conflicting single cell cluster annotations. It scores every
candidate cluster for biological stability. Then it runs a Shapley-style game per cell,
so the annotations compete and one wins.

This directory holds the upstream package, rewritten for speed and for macOS.

## Provenance

| Item | Value |
| --- | --- |
| Upstream repository | https://github.com/frankligy/scTriangulate |
| Upstream version | 0.13.0 |
| Upstream commit | `8b9598cfcbf269856f41740514cadcd75f1ee2c6` (2026-03-30) |
| Integrated on | 2026-08-06 |
| Component path | `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/sctriangulate` |

`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/sctriangulate/_reference/upstream_shapley.py`
holds the untouched upstream Shapley code. The equivalence tests compare against that file.

## How to run

Command line:

```
cd /Users/saljh8/Documents/GitHub/altanalyze3
python3.11 -m altanalyze3.components.sctriangulate.cli \
    --h5ad   /path/to/input.h5ad \
    --query  ann1,ann2,ann3 \
    --outdir /path/to/output
```

Python:

```python
import scanpy as sc
from altanalyze3.components.sctriangulate import ScTriangulate

adata = sc.read('/path/to/input.h5ad')
sctri = ScTriangulate(dir='/path/to/output', adata=adata, query=['ann1', 'ann2', 'ann3'])
sctri.lazy_run()
```

The CLI writes `run_parameters.json` beside the results. That file records every parameter,
the input path and the version of each library.

Use `/opt/homebrew/opt/python@3.11/bin/python3.11` on this Mac. The default `python3` lacks numpy.

### Every flag

Required: `--h5ad`, `--query`, `--outdir`.

**These four change the answer.** Each one belongs in a methods section, and each records
what it did in `run_parameters.json`.

| Flag | Default | What it does |
| --- | --- | --- |
| `--downsample N` | off | Cap cells per cluster with one whitelist shared across annotations. Cell counts are the denominator of every stability metric. |
| `--prefilter-markers N` | off | Keep the union of the top N MarkerFinder genes per cluster. Circular by construction; read the table below before using it. |
| `--scale-sccaf` | off | Mean-centre before the SCCAF fit. Also forces the legacy dense SCCAF path, because centring destroys sparsity. |
| `--criterion 1..6` | 2 | Which gene classes count as artifacts. See `metrics.read_artifact_genes`. |

**These change how the code reaches the same answer.**

| Flag | Default | What it does |
| --- | --- | --- |
| `--sequential` | off | One annotation at a time. About half the peak memory for about 9 s more on an 8,000-cell run. |
| `--cores N` | `min(n_annotations, n_cpus)` | Worker processes for the metric step. |
| `--sccaf-mode` | `optimized` | `optimized` keeps the matrix sparse; `legacy` is upstream's dense procedure. Same accuracies. |
| `--enrichment` / `--no-enrichment` | off | The enrichr and GSEA artifact annotation. Feeds no metric and no Shapley decision. `--viewer-cluster` turns it on for you. |
| `--shapley-mode` | `shapley_all_or_none` under 16 annotations, else `rank` | Which decision rule picks the winning annotation. |
| `--shapley-bonus` | 0.01 | Margin a challenger must beat. |
| `--win-fraction-cutoff` | 0.25 | Minimum share of its own cells a cluster must keep to survive pruning. |
| `--reassign-abs-thresh` | 10 | Clusters below this many cells go to the invalid list. |
| `--layer NAME` | none | Layer holding raw counts, for tf-idf when X is skewed. |
| `--prefilter-cells N` | 1000 | Cells per cluster used to rank markers for `--prefilter-markers`. |
| `--subsample-seed N` | 0 | Seeds `--downsample` and `--prefilter-markers`. |
| `--species` | human | `human` or `mouse`. Picks the artifact gene database. |
| `--reference NAME` | first `--query` entry | Annotation used as the reference for pruning. |
| `--predict-doublet` | off | Run scrublet. Off means a constant doublet score of 0.5, as upstream. |
| `--assess-raw`, `--assess-pruned` | off | Score the raw and pruned labels as if they were another annotation. |
| `--viewer-cluster`, `--viewer-heterogeneity` | off | Write the HTML viewer pages. Untested here. |

**This one changes the output file, not the result.**

| Flag | Default | What it does |
| --- | --- | --- |
| `--keep-layers` | off | Keep every layer. Without it, layers that `--layer` does not name are freed from memory and are absent from `sctriangulate.h5ad`. |

### Defaults that differ from upstream scTriangulate

| Behaviour | Upstream | Here | Why |
| --- | --- | --- | --- |
| Enrichment columns | always on | off | Feed no metric and no Shapley decision; 3.7 s of a 21 s demo run |
| SCCAF | dense, unseeded | sparse, seeded | Same accuracies; upstream matched its own first run in 1 to 3 of 50 repeats |
| Compute on X | densified | sparse | 142.7 GB dense at 984,119 by 36,249; identical labels, 2.5x less memory |
| Reassign PCA | randomized SVD | exact SVD | Randomized missed the true subspace by 0.8% of variance and drove the instability |
| Unused layers | kept | freed | Nothing reads them; `--keep-layers` restores |

The tables below measure every one of these.

## Steps 4 to 6 after pruning, on by default

`lazy_run` stops at `obs['pruned']`. Those labels name a source cluster, such as
`ICGS3_RNA@C3`, and no biology. Every reader then repeated the same three steps by hand.
`annotate.py` runs them, and the CLI runs it by default.

| Step | What runs | Which altanalyze3 function |
| --- | --- | --- |
| 4 | MarkerFinder on the pruned clusters, 60 markers and 100 displayed cells per cluster | `visualization.marker_heatmap_h5ad.generate_marker_heatmap_from_adata` |
| 5 | Cell-state naming by GO-Elite BioMarkers gene-set enrichment | `clustering.ICGS.biomarker_enrichment`, `clustering.ICGS.clean_biomarker_prediction_labels` |
| 6 | HOPACH over the cluster centroids, then MarkerFinder again in that order | `clustering.hopach.hopach` |

Nothing in `annotate.py` reimplements a method. Each step calls the function the matching CLI
calls, with the same arguments.

New outputs under `--outdir`:

```
MarkerFinder/pruned/               step 4: markers, redundant markers, fold matrix, centroids, heatmap
GO-Elite/                          step 5: enrichment table and the cell-state prediction per cluster
hopach/hopach_centroid_clusters.tsv    step 6: order, group and name source per cluster
hopach/hopach_ordered_centroids.tsv    centroids, columns in HOPACH order, columns renamed
MarkerFinder/hopach_ordered/       step 6: the heatmap redrawn in HOPACH order
MarkerFinder/hopach_ordered_ADT/   the same over AB_ features alone, when the object carries them
```

New `obs` and `uns`: `obs['cluster_name']`, `obs['hopach_cluster']`, `uns['lineage_order']`
(which `marker_heatmap_h5ad._resolve_cluster_order` reads) and `uns['sctriangulate_annotation']`.

Turn the whole stage off with `--no-annotate`, which restores the earlier behaviour exactly:
the run then writes no `cluster_name` and stops at `obs['pruned']`.

### Naming, and the two gates that stop it inventing a cell type

`biomarker_enrichment` returns the single best term per cluster whether or not the term means
anything. Two gates drop a weak term, and the cluster then takes its strongest marker as a name,
such as `Zim3-high (streets_celltype@Mature_cycling)`.

| Gate | Default | Why that value |
| --- | --- | --- |
| `--annotate-max-fdr` | `1e-5` | Measured over the 79 clusters of three runs on one mouse thymus dataset. Best-term FDR ran from 2.1e-71 to 5.3e-02. Every term that FDR 0.05 accepted and 1e-5 rejected, 18 of them, overlapped by only 2 to 8 of about 36 to 60 markers. FDR 0.05 keeps 78 of 79 clusters; 1e-5 keeps 60 of 79. |
| `--annotate-min-overlap` | `5` | An overlap of 2 genes names nothing, however small its FDR. |

Without the gates the same dataset produced `Retina Rheaume et al.` at FDR 0.027 and
`Marrow Oligodendrocyte` at FDR 0.049 as cluster names.

`--annotate-lead OBS_COLUMN` lets a curated annotation lead instead, for example a published cell
type that also competes in `--query`. Its dominant label then names the cluster and the enriched
term becomes the fallback. A dominant label that means "no call", such as `unassigned`, never
names a cluster; see `UNINFORMATIVE_LEAD_LABELS`.

### Validation

Measured 2026-08-17 against a hand-run of the same three steps on a mouse thymus CITE-seq
dataset, 13,210 cells by 15,495 features, 30 pruned clusters from four annotations.

| Check | Result |
| --- | --- |
| Step-4 marker table against the hand run | 1,749 of 1,749 rows, identical `(cluster, gene)` index, `Fold`, `Query Exp`, `Ref Exp` and `FDR p-value` differ by 0.0 |
| Step-4 centroid matrix | shape and columns identical, max absolute difference 0.0 |
| Step-6 HOPACH order and group membership | identical, k = 11, MSS = 0.3949 reproduced |
| End to end through the CLI | `pruned` labels identical on 13,210 of 13,210 cells |
| `--no-annotate` | writes no `cluster_name`, same 30 pruned clusters |
| Tests | 9 new in `tests/test_annotate.py`; 69 of 69 sctriangulate tests pass |

Scripts: `27_validate_annotate_default.py`, `28_validate_annotate_default_paths.py` and
`29_e2e_scTriangulate_annotate_default.sh` under
`/Users/saljh8/Dropbox/Collaborations/Grimes/Thymus/Streets-CITE/scripts`.

### One dependency carries a known defect

`clustering/hopach.py` fails 2 of its own 3 tests at commit `2f93c2b`, before and after this
change. On the toy two-cluster dataset in `test_hopach.py` it returns k = 4 with sizes 1, 2, 1, 2
where the test expects k = 2. The ORDER stays correct, and step 6 uses the order, so the heatmap
is unaffected. `obs['hopach_cluster']` may over-split. I did not touch `hopach.py`.

## Headline result

I ran both versions through `lazy_run`, the documented top-level entry point, on the demo
dataset upstream ships. Every number below comes from
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/sctriangulate/benchmarks/end_to_end.py`.

| Run | Wall clock | Speedup |
| --- | ---: | ---: |
| upstream, one process | 130.76 s | — |
| upstream, parallel (its default) | 85.41 s | — |
| this version, one process, defaults | 17.8 s | **7.35x** vs upstream one process |
| this version, parallel, defaults | 16.7 s | **5.11x** vs upstream parallel |
| this version, parallel, vs upstream one process | 16.7 s | **7.83x** |
| this version, one process, `--enrichment` | 21.4 s | 6.11x vs upstream one process |

The upstream rows and the `--enrichment` row come from `benchmarks/end_to_end.py`. The two
default rows come from the CLI, which now leaves the artifact enrichment off. Both harnesses
time the same span: object construction plus `lazy_run`.

The metric and Shapley stages alone, timed one at a time with enrichment on, went from
135.14 s to 19.07 s, which is 7.09x. The Shapley engine on its own reaches 205x to 138,563x,
and the gap widens with the number of annotations. See the tables below.

**I did not reach 10x.** The demo lands at 7.35x like for like, and a larger simulated dataset
at 4.56x on the staged metrics. One step sets the ceiling: the SCCAF L1 logistic regression,
a single-threaded C solver. I tried two ways to cut it and measured both; each changed the
metric for a gain of 1.03x to 1.56x, so I shipped neither. Section "Where the remaining time
goes" has the numbers.

## Fidelity: this version reproduces upstream better than upstream reproduces itself

### There is no answer key

The upstream repository ships no expected results. `test/mini_test.py` contains no assert
statement, so it checks only that the code runs without raising. The demo h5ad carries no
ground-truth cell type column, only the three leiden annotations that scTriangulate must
reconcile. So "correct" can only mean "what upstream's own code produces".

Upstream also leaves two sklearn estimators unseeded, so two upstream runs on identical input
disagree with each other. That variability is the yardstick, which is why the harness runs
each implementation twice, under `PYTHONHASHSEED` 1000 and 2000.

### Exact label agreement, ARI and NMI

`benchmarks/label_agreement.py` scores every pair. ARI ignores label names and compares the
partitions; exact agreement demands the identical label string.

| Comparison | Column | Exact | ARI | NMI |
| --- | --- | ---: | ---: | ---: |
| upstream reproduces itself | raw | 2687/2700 | 0.9956 | 0.9951 |
| | pruned | 2567/2700 | 0.9640 | 0.9236 |
| **this version reproduces itself** | raw | **2700/2700** | **1.0000** | **1.0000** |
| | pruned | **2700/2700** | **1.0000** | **1.0000** |
| this version, one process vs parallel | raw | 2700/2700 | 1.0000 | 1.0000 |
| | pruned | 2700/2700 | 1.0000 | 1.0000 |
| this version vs upstream | raw | 2692/2700 | 0.9979 | 0.9972 |
| | pruned | 2586/2700 | 0.9679 | 0.9311 |
| this version vs upstream, other seed | raw | 2693/2700 | 0.9976 | 0.9973 |
| | pruned | 2582/2700 | 0.9680 | 0.9305 |

Read the first and fourth blocks together. This version matches upstream at ARI 0.9979 on the
raw labels and 0.9679 on the pruned labels. Two upstream runs match each other at 0.9956 and
0.9640. So this version agrees with upstream more closely than upstream agrees with itself,
and it returns the identical answer on every repeat.

Upstream produced 46 raw clusters in one run and 45 in the other. This version produced 45 in
all three of its runs.

## Stage timings on the demo dataset

The demo is PBMC3k, 2,700 cells and 32,738 genes, with 9, 19 and 33 clusters in the three
annotations. Both runs used one process and ran the stages in the same order.

Input file: `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/sctriangulate/benchmarks/data/demo_pbmc3k.h5ad`
(sha256 `91ede794e28352731186b8b2dd917ea1d6a573e5e73a3adc559bc746d08ddca2`)

| Stage | Upstream 0.13.0 | This version | Speedup |
| --- | ---: | ---: | ---: |
| marker_gene, 3 annotations | 96.31 s | 8.14 s | 11.8x |
| tf-idf 10 and tf-idf 5 | 24.32 s | 1.49 s | 16.3x |
| SCCAF | 6.48 s | 6.41 s | 1.01x |
| reassign | 3.18 s | 2.87 s | 1.11x |
| Shapley, 4 modes | 4.61 s | 0.10 s | 46.1x |
| run_assign | 0.20 s | 0.00 s | — |
| doublet | 0.04 s | 0.00 s | — |
| **Total** | **135.14 s** | **19.07 s** | **7.09x** |

Repeated runs of the same command vary by a few percent on this machine. An earlier run of
the same code measured 17.94 s.

### A larger, simulated dataset

`benchmarks/simulate_dataset.py` writes 8,000 cells by 20,000 genes, 6.1% nonzero, with
10, 20 and 30 clusters in three nested annotations. It borrows real gene symbols from the demo
file, so 1,053 of the 20,000 genes are human artifact genes and the enrichment step costs what
it would cost on real data.

| Stage | Upstream | This version | Speedup |
| --- | ---: | ---: | ---: |
| marker_gene | 103.85 s | 11.07 s | 9.4x |
| tf-idf 10 and tf-idf 5 | 41.94 s | 1.83 s | 22.9x |
| SCCAF | 19.32 s | 19.11 s | 1.01x |
| reassign | 9.32 s | 9.05 s | 1.03x |
| Shapley, 4 modes | 12.79 s | 0.05 s | 256x |
| run_assign | 0.52 s | 0.00 s | — |
| **Total** | **187.81 s** | **41.15 s** | **4.56x** |

The overall factor drops from 7.09x to 4.56x because the two steps I could not change,
SCCAF and reassign, grow with cell count and now take 28.16 s of the 41.15 s. Amdahl's law
caps the achievable factor at 187.81 / 28.16 = 6.7x for this shape.

The same equivalence check ran on this dataset. tf-idf 10, tf-idf 5, doublet, the marker gene
lists, the purified lists, the top-50 exclusive genes, enrichr and GSEA all matched exactly,
across all 60 clusters and 600 enrichment values. That repeats the demo result on data with a
different size, density, cluster count and gene composition.

Inside `marker_gene`, measured separately:

| Step | Upstream | This version |
| --- | ---: | ---: |
| scanpy `rank_genes_groups` | 11.55 s | 3.49 s |
| gene-to-cluster assignment | 15.20 s | 0.01 s |
| enrichr over artifact classes | about 30 s | 0.04 s |
| gseapy GSEA prerank | about 40 s | 3.68 s |

Logs:
`benchmarks/out_baseline/baseline_profile.log` and `benchmarks/out_optimized/optimized_profile.log`,
both under the component path given above.

## Shapley against the number of annotations

The demo holds 3 annotations, where each player faces only 4 coalitions and the upstream loop
looks cheap. The loop costs `O(2**(n-1))` per player per cell, so I measured the scaling on
simulated score matrices with 20,000 cells and 4 metrics. `benchmarks/shapley_scaling.py`
builds them, and it checks every configuration against the verbatim upstream functions.

| Annotations | Coalitions per cell | Upstream | Batched | Speedup | Distinct matrices | Verified |
| ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| 3 | 9 | 6.69 s | 0.033 s | 205x | 1,728 of 20,000 | yes |
| 4 | 28 | 20.69 s | 0.046 s | 447x | 12,859 of 20,000 | yes |
| 5 | 75 | 54.64 s | 0.055 s | 992x | 19,204 of 20,000 | yes |
| 6 | 186 | 136.20 s | 0.064 s | 2,124x | 19,932 of 20,000 | yes |
| 8 | 1,016 | 742.36 s | 0.082 s | 9,076x | 20,000 of 20,000 | yes |
| 10 | 5,110 | 3,718.60 s | 0.121 s | 30,806x | 20,000 of 20,000 | yes |
| 12 | 24,564 | 17,988.35 s | 0.130 s | 138,563x | 20,000 of 20,000 | yes |

At 12 annotations upstream needs 5.0 hours for 20,000 cells. The batched engine needs 0.13 s.

I measured upstream on the first 300 cells of each configuration and scaled the time to 20,000.
The upstream loop touches one cell at a time and shares nothing between cells, so the cost
grows linearly in cells. Those same 300 cells also carry the equivalence check.

The "distinct matrices" column shows the two wins separating. Below 5 annotations the
deduplication does much of the work. Above 8, every cell holds a distinct matrix and the
closed form does all of it.

## What I changed, and why the numbers still match

### 1. Shapley: closed form instead of coalition enumeration

Upstream scores one player of one cell by looping over every non-empty coalition of the
other `n-1` players. The loop costs `O(2**(n-1))` per player per cell.

I derived a closed form. For each score column, three counts decide the whole sum:
how many other players the cell already beats, how many the bonus rescues, and how many it
loses to. A small integer lookup table, indexed by those counts, replaces the loop.
`shapley.py` carries the derivation.

The tables hold integers scaled by `n!`, and the code divides once at the end. Upstream adds
thousands of floating point products, so upstream carries rounding noise that the table
avoids. Two players with mathematically equal values come out bitwise equal here, which makes
the tie-break deterministic.

Second win: a cell's score matrix depends only on which cluster it falls in per annotation.
Cells that share a cluster tuple share a bitwise identical matrix. The engine scores each
distinct matrix once.

### 2. marker_gene: the gene assignment loop

Upstream assigned each gene to its best cluster with `np.nonzero(df['names'].values == gene)`
inside a gene loop inside a cluster loop, so it cost `O(G**2 * C)` comparisons. On the
33-cluster annotation that reaches 3.5e10 element comparisons.

I scatter each cluster's ranks into a `(genes, clusters)` matrix with one vectorised
`get_indexer` per cluster, then take the argmin. The tie-break and the stable sort match
upstream, so the output lists match exactly.

### 3. scanpy group statistics

scanpy computes the mean and variance of each group AND of everything outside each group by
slicing X, so it makes `2C` passes. Group sums and group sums-of-squares determine all of it,
and one chunked matmul produces both.

I patch only `_RankGenes._basic_stats`, and scanpy does the rest. scanpy squares in the input
dtype, so I square in the input dtype too. On the demo data the aligned statistics came out
bitwise identical, and all 61 marker gene lists matched.

The patch refuses the job when the group masks do not cover every cell, because the "rest"
statistics come from a column total minus one group. If a scanpy upgrade moves the method,
the patch logs a warning and steps aside. A scanpy upgrade can make the run slower. It cannot
make the run wrong.

### 4. tf-idf

Upstream built a dense `cells x genes` DataFrame, then recounted the nonzeros of the WHOLE
matrix once per cluster. It also read the 18,000-row artifact table once per cluster, and it
recomputed the same vectors three times for tf-idf 1, 5 and 10.

One indicator matmul in float64 now counts every cluster at once. Sums of 0 and 1 below 2**53
are exact in float64, so the counts equal `np.count_nonzero` exactly. The artifact table and
the code caches the sorted vectors.

### 5. enrichr

Upstream called `gseapy.enrichr` per cluster, which built an object, ran a hypergeometric test
and wrote a directory of output. I inlined the test. The term order, the `x < 1` skip, the
hypergeometric call and the Benjamini-Hochberg code come from gseapy 0.10.4 line for line.
`tests/test_metrics_equivalence.py` checks the two against each other.

### 6. GSEA

The enrichr and GSEA columns are now **off by default**. Neither feeds a stability metric or
a Shapley decision. Only `plot_cluster_feature(feature='enrichment')`,
`plot_multi_modal_feature_rank`, `penalize_artifact(mode='cellcycle')` and the viewer read
them, and those four raise a message naming `--enrichment` when the columns are absent.
`--viewer-cluster` turns them on for you. Measured: 21.4 s with, 17.8 s without, and the two
runs agree on 2700 of 2700 cells for every label column.

When they are on, I kept gseapy, because the NES depends on a seeded permutation stream. I
passed `outdir=None, no_plot=True, verbose=False`. Upstream wrote a figure and a table per
cluster. One call dropped from 0.714 s to 0.119 s with identical NES.

Upstream wrapped the call in a bare `except:` that returned zeros without a word. The call
does fail on real input: 6 of the 33 clusters in `sctri_rna_leiden_3` raise
`No gene sets passed through filtering condition`, because their marker list overlaps no
artifact class. This version logs each one and then fills the same zeros, so the numbers stay
comparable and the user learns it happened.

### 7. Exclusive gene storage

`uns['exclusive_genes'][annotation][cluster]` was a plain dict with one entry per gene.
The 61 clusters held 1,901,980 entries between them, and `lazy_run` pickles the object three
times. Pickling the dict form took 2.28 s and 46.2 MB per checkpoint.

`ExclusiveGeneScores` stores the same ordered mapping as a gene index plus a float64 array.
Pickling takes 0.09 s and 25.3 MB, a 25x drop. Every caller keeps working, because the class
is a `Mapping` and supports `keys()`, `values()`, `items()`, `x['GENE']`, `len(x)`, iteration
and `to_dict()`. It also supports `x[:10]`, which one plotting function already assumed and
which a dict never allowed.

### 8. Per-cell pandas loops

`run_assign`, `_prefixing`, `which_to_take` and two pruning loops walked cells one at a time
with `.iloc`. Each now runs as a vectorised pass. Two of the pruning loops also wrote through
a chained indexer, which pandas warns about.

### 9. penalize_artifact

Upstream rebuilt the label matrix and rewrote the metric block once per stamp. A cell matches
at most one stamp, so one membership test gives the same answer in a single pass.

## macOS fixes

| Problem | Effect on this Mac | Fix |
| --- | --- | --- |
| `import umap` at module scope pulls in TensorFlow | `import sctriangulate` raised `AttributeError: module 'numpy' has no attribute 'dtypes'` under numpy 1.23.5, so nothing ran at all | `_compat.py` loads umap on first use. Upstream touches it in one function. |
| `import scrublet` at module scope | Import cost, and a hard dependency the default path never uses | Deferred the same way |
| `import squidpy` at module scope in `spatial.py` | Same | Deferred the same way |
| Upstream parallel mode deadlocks | macOS defaults multiprocessing to `spawn`. Each worker re-imports the package, dies on the umap import, and the parent waits at `pool.join()` forever. I killed one run after 14 minutes with zero output. | The component defers the failing import and uses a `fork` context on Darwin |
| `pool.apply_async(each_key_run, (self, ...))` | `apply_async` pickles its arguments whatever the start method, so upstream shipped a full AnnData per worker | Under `fork`, the object is published as a module global before forking, so children inherit it copy-on-write |
| Worker processes oversubscribe BLAS | Each worker's Accelerate pool also claims every core | `_single_threaded_blas` pins BLAS to one thread inside workers |
| `subprocess.run(['rm','-r', path])` in 7 places | Deleted a directory through a formatted path, and printed an error when the directory was absent | `_cleanup_enrichr_tmp` deletes nothing. The new enrichr writes no scratch directory. A leftover directory gets reported, not removed. |
| `_to_dense()` called `.toarray()` unconditionally | Raised `AttributeError` when X was already dense | Guarded with `issparse` |
| `tf_idf*` used a bare `.cat.categories` | Raised unless scanpy had sanitized the column first | `astype('category')` first, which is a no-op when the column already is one |
| `regress_size=True` | The boolean parameter shadowed the function of the same name, so any caller passing True hit `TypeError: 'bool' object is not callable` | Callers use `_regress_size_impl`. The public name stays. |

To reproduce the deadlock claim, `benchmarks/_shimpath/sitecustomize.py` exists only so the
upstream timing runs reach their spawned children. Without it upstream's parallel mode never
returns on this machine.

## Equivalence

`benchmarks/compare_outputs.py` compares 46 fields between the two runs on the demo data.

Exact match, 0.000e+00 maximum difference, for every one of these:

- tf-idf 10, tf-idf 5 and doublet scores, for all 61 clusters
- marker gene lists, purified gene lists and top-50 exclusive genes, for all 61 clusters
- enrichr `-log10(adjusted p)`, 305 of 305 values
- GSEA NES and hit counts, 305 of 305 values
- Shapley scores and winners in `rank_all_or_none` mode, 8,100 of 8,100 values

`reassign` and `SCCAF` do not match, because upstream does not return the same number twice.
The next section measures that.

Unit tests:

- `tests/test_shapley_equivalence.py`, 29 tests. It checks the batched engine against the
  verbatim upstream functions across 4 modes and 2 to 7 players. The cases include exact ties,
  gaps inside the bonus window, all-zero matrices and coarse grids.
- `tests/test_metrics_equivalence.py`, 14 tests. It checks the gene assignment, the tf-idf
  matrix for dense and sparse input, `doublet_compute`, `run_enrichr` against gseapy,
  `penalize_artifact_void`, the artifact table cache and the cluster sizes.
- `tests/test_ttest_equivalence.py`, 6 tests. It checks the group statistics patch.

Run them:

```
cd /Users/saljh8/Documents/GitHub/altanalyze3
python3.11 -m pytest altanalyze3/components/sctriangulate/tests -q
```

All 49 pass.

## Two metrics do not match, and upstream is the reason

`reassign_score` and `SCCAF_score` are the two stability metrics whose values differ.
Neither upstream function returns the same number twice on identical input.

I ran each upstream function 50 times per annotation, across 5 processes with different
`PYTHONHASHSEED` values. The script is `benchmarks/determinism_audit.py`.

| Annotation | Metric | Upstream runs equal to the first | Upstream spread, one process | Upstream spread, across processes |
| --- | --- | ---: | ---: | ---: |
| leiden_1 (9 clusters) | reassign | 3 of 50 | 0.006 | 0.006 |
| leiden_1 | SCCAF | 1 of 50 | 0.012 | 0.200 |
| leiden_2 (19 clusters) | reassign | 1 of 50 | 0.200 | 0.243 |
| leiden_2 | SCCAF | 1 of 50 | 0.056 | 0.056 |
| leiden_3 (33 clusters) | reassign | 1 of 50 | 0.258 | 0.290 |
| leiden_3 | SCCAF | 1 of 50 | 0.167 | 0.167 |

Both metrics are accuracies between 0 and 1. A spread of 0.29 on repeated identical input
makes the raw value hard to trust, and the Shapley game reads both values.

Three unseeded sources cause the spread:

1. `PCA(n_components=30)` in `reassign_score` leaves `random_state=None`. At this shape
   sklearn picks the randomized solver, so repeated calls differ inside one process.
2. `pool = list(set(pool))` in `reassign_score` orders genes by string hash. Python randomizes
   string hashing per interpreter, so the gene column order differed between processes,
   including between the worker processes that upstream spawns.
3. `LogisticRegression(solver='liblinear')` in `SCCAF_score` leaves `random_state=None`.
   sklearn uses that argument to shuffle the data for the liblinear solver.

I sorted the pool and seeded both estimators, in `metrics.py` and in `prune.py`. This version
returned one value across 5 processes and 3 repeats each, for all 6 metric-annotation pairs.

For 121 of 122 cluster-metric values, this version sits inside the range upstream produced.
One cluster differs, `reassign` on leiden_2 cluster 9: this version gives 0.3036, upstream's
50 runs gave 0.1518 to 0.2857. The gap of 0.018 is smaller than that annotation's own mean
upstream spread of 0.089.

Seeding fixes reproducibility. It does not make either metric more accurate. A user who cares
about the absolute value of `reassign` on a small cluster should treat it as uncertain to
about plus or minus 0.1 on this dataset, whichever version produced it.

## Memory, and running this on a laptop

`benchmarks/memory_profile.py` measures peak resident memory two ways: the kernel's
high-water mark from `/usr/bin/time -l`, and the summed resident size of the whole process
tree sampled every 0.25 s. The second one catches concurrent workers, which the first does not
sum. Sampling can miss a spike shorter than the interval, so the sampled figure is a lower
bound.

Simulated dataset, 8,000 cells by 20,000 genes, 3 annotations:

| Version | Mode | Kernel peak | Process tree peak | Seconds |
| --- | --- | ---: | ---: | ---: |
| upstream | one process | 5,182 MB | 4,956 MB | 197.0 |
| this version | one process | 4,523 MB | 4,219 MB | 41.3 |
| upstream | parallel | 4,413 MB | 10,205 MB | 115.3 |
| this version | parallel | 3,364 MB | 9,418 MB | 32.1 |

Demo dataset, 2,700 cells by 32,738 genes, 3 annotations:

| Version | Mode | Kernel peak | Process tree peak | Seconds |
| --- | --- | ---: | ---: | ---: |
| upstream | one process | 5,836 MB | 5,582 MB | 134.1 |
| this version | one process | 3,285 MB | 3,140 MB | 17.6 |
| upstream | parallel | 8,381 MB | 14,452 MB | 89.1 |
| this version | parallel | 2,065 MB | 6,583 MB | 16.6 |

**Read this honestly: memory is better, but not by the margin that the runtime improved.**
This version cuts memory 1.15x to 1.31x on the larger dataset and 1.8x to 4.1x on the demo.
Two things I did not change hold the peak up, because changing them would change the numbers:
the dense expression matrix that scanpy and sklearn both need, and the float64 copy that
liblinear makes of the SCCAF training half.

Practical guidance for a laptop:

- **Use `--sequential`.** On the simulated dataset it needs 4.2 GB against 9.4 GB for the
  parallel mode, and costs 9 seconds more (41.3 s against 32.1 s). Three workers each need
  their own working memory; only the shared expression matrix is free.
- **Budget roughly 8 times the size of the dense float32 matrix** for a one-process run.
  Measured: 3.3 GB against a 354 MB matrix on the demo, 4.5 GB against a 640 MB matrix on the
  simulation. Multiply cells by genes by 4 bytes, then by 8.
- A 50,000 cell by 25,000 gene dataset implies a 5 GB dense matrix, so about 40 GB peak in one
  process. That does not fit a 16 GB laptop, in either version.

Two memory faults found and fixed while measuring this, both mine, not upstream's:

1. The tf-idf rewrite built a float64 nonzero indicator for the whole matrix, an extra
   `cells * genes * 8` bytes: 707 MB on the demo and 10 GB at 50,000 by 25,000. It now counts
   in 128 MB chunks, so that allocation is flat regardless of dataset size.
2. Under `fork`, each worker called `.toarray()` and allocated its own dense copy, 640 MB
   each. The parent now densifies once before forking and the workers share those pages
   copy-on-write. Parallel kernel peak on the simulation fell from 3,949 MB to 3,364 MB.

`SCCAF_score` also now frees the full matrix before fitting, since the two halves are already
independent copies. Neither fix changes any output; all 49 tests still pass.

## Where the remaining time goes

Of the 21.07 s that `lazy_run` took in one process on the demo WITH the enrichment columns
(17.8 s without them, which is now the default):

| Part | Time | Can it shrink? |
| --- | ---: | --- |
| SCCAF liblinear L1 fit | 6.4 s | Not without changing the numbers. The solver is C, single threaded, and the cost grows with cells and with cluster count. It sets the ceiling. |
| gseapy GSEA prerank | 3.7 s | Already off by default. Pass `--enrichment` to get it back. |
| scanpy `rank_genes_groups`, remainder | 3.5 s | The group statistics are already one pass. What is left is scanpy's own sorting, Benjamini-Hochberg and result assembly. |
| reassign PCA and KNN | 2.9 s | Small. |
| tf-idf, pruning, plots, h5ad write, 3 pickled checkpoints | about 4 s | The checkpoints already dropped from 2.35 s to 0.09 s each. |

Turning the enrichment columns off by default removes 3.7 s. The CLI measures it:
`--sequential` takes 21.4 s with `--enrichment` and 17.8 s without, and the two runs agree on
2700 of 2700 cells for `final_annotation`, `raw` and `pruned`. So upstream sequential to this
version sequential is 130.76 / 17.8 = 7.35x.

### SCCAF: an optimized default, with the legacy path kept

`SCCAF_score` now has two modes. `--sccaf-mode optimized` is the default;
`--sccaf-mode legacy` runs upstream's procedure. A seed fixes both, so both
reproduce; upstream did neither.

The optimized mode keeps the matrix sparse. `compute_metrics` parks the original
CSR in a layer before densifying for the other metrics, so SCCAF reads it directly
and never converts. liblinear keeps its own one-versus-rest loop, so it sees the
same numbers in a smaller container.

| Dataset | Annotation | Clusters | Legacy | Optimized | Speedup | Legacy matrix | Optimized matrix | max abs delta |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| demo | leiden_1 | 9 | 1.01 s | 0.54 s | 1.86x | 337 MB | 16 MB | 0.0000 |
| demo | leiden_2 | 19 | 1.69 s | 1.23 s | 1.37x | 337 MB | 16 MB | 0.0000 |
| demo | leiden_3 | 33 | 2.67 s | 2.24 s | 1.20x | 337 MB | 16 MB | 0.0000 |
| sim 8k | leiden_10 | 10 | 3.45 s | 2.59 s | 1.33x | 606 MB | 75 MB | 0.0000 |
| sim 8k | leiden_20 | 20 | 5.37 s | 4.54 s | 1.18x | 606 MB | 75 MB | 0.0000 |
| sim 8k | leiden_30 | 30 | 7.45 s | 6.63 s | 1.12x | 606 MB | 75 MB | 0.0000 |

The matrix handed to the solver shrinks 8x to 21x. Every one of the 121 cluster
scores across both datasets is identical to the legacy value, to the last bit.

End to end through `lazy_run`, the labels are the same partition:

| Dataset | Column | Exact | ARI | NMI |
| --- | --- | ---: | ---: | ---: |
| demo | final_annotation | 2700/2700 | 1.0000 | 1.0000 |
| demo | raw | 2700/2700 | 1.0000 | 1.0000 |
| demo | pruned | 2700/2700 | 1.0000 | 1.0000 |
| sim 8k | final_annotation | 8000/8000 | 1.0000 | 1.0000 |
| sim 8k | raw | 8000/8000 | 1.0000 | 1.0000 |
| sim 8k | pruned | 8000/8000 | 1.0000 | 1.0000 |

`benchmarks/sccaf_compare.py` produces both tables.

**One thing that did not work, and why the default is not more aggressive.**
Splitting liblinear's one-versus-rest loop into separate per-class fits allows
real parallelism. On threads it reached 2.5x to 4.0x, and it was wrong:
liblinear seeds and draws from a process-global C `rand()`, so concurrent
threads interleave on one stream and repeated runs disagreed. The measurement
caught it. Processes reproduce, and the per-class values are bitwise identical
whether run serially, on 2 processes or on 10, but they differ from the legacy
values by up to 0.167 and the pool startup outweighs the gain below about 20
clusters. That route is available as `metrics.SCCAF_PER_CLASS = True`, off by
default.

`tests/test_sccaf_modes.py` holds the guards: optimized equals legacy exactly,
both modes reproduce, the per-class mode does not depend on the worker count,
and the sparse view equals the dense subset.

### The two routes I tried before this, which both failed

I tested the two routes that could have reached 10x. `benchmarks/sccaf_parallel_probe.py`
runs the first. Both changed the metric, so I shipped neither.

**Route 1, fit each one-versus-rest sub-problem in a separate process.** sklearn hands the
whole multiclass problem to liblinear, which loops the classes inside C. Splitting that loop
across cores should have been a straight win. Neither axis delivered:

| Annotation | Classes | sklearn | Per-class, 8 cores | Speedup | Same accuracies? |
| --- | ---: | ---: | ---: | ---: | :--- |
| leiden_1 | 9 | 0.65 s | 4.55 s | 0.14x | no, 2 of 9 clusters differ, by up to 0.006 |
| leiden_2 | 19 | 1.35 s | 0.98 s | 1.38x | no, 5 of 19 clusters differ, by up to 0.056 |
| leiden_3 | 33 | 2.42 s | 1.55 s | 1.56x | no, 3 of 33 clusters differ, by up to 0.167 |

The coefficients differ by up to 0.26 in absolute value. sklearn draws one seed and
liblinear's coordinate descent consumes that random stream across the sub-problems in
sequence, so sub-problem k depends on 0 through k-1. Fitting them separately restarts the
stream, and each solver stops at `tol=1e-4` on a different iterate. At 9 classes the process
startup also costs more than the work saved.

**Route 2, drop features that are all-zero in the training split.** Half the matrix qualifies,
16,349 of 31,180 columns, so this looked free.

| Annotation | All-zero features | Full | Pruned | Speedup | Same accuracies? |
| --- | ---: | ---: | ---: | ---: | :--- |
| leiden_1 | 16,349 of 31,180 | 0.65 s | 0.63 s | 1.03x | no, 3 of 9 clusters differ |
| leiden_2 | 16,310 of 31,180 | 1.38 s | 1.35 s | 1.02x | no, 2 of 19 clusters differ |
| leiden_3 | 16,326 of 31,180 | 2.44 s | 2.24 s | 1.09x | no, 2 of 33 clusters differ |

liblinear already skips zero columns almost for free, so removing them buys 1.03x to 1.09x.
It still changes the answer, because the feature count changes the shuffled coordinate order.

Both differences are the same size as upstream's own run-to-run spread, so neither is more
wrong than upstream. Neither is worth trading a reproducible number for at 1.1x to 1.6x.
The SCCAF fit stays as it is.

## Three memory changes that keep every number

All three came from a review that proposed ten. I dropped seven: measurement
refuted two, three pointed the right way but named the wrong size or the wrong
place, and two need the real object before anyone can decide them. I measured
the three below, and each one produces labels identical to the labels before it.

**A. Share the index arrays instead of copying the matrix.** Two places needed a
matrix with the same sparsity pattern and different values: the t-test squares
the data, the tf-idf indicator replaces it with ones. Both called `src.copy()`,
which duplicates `indices` and `indptr` and then discards the data it just
copied. `metrics._shared_index_csr` builds the new matrix over the original
index arrays. It copies instead when the indices arrive unsorted, because scipy
would sort them in place and permute them out of step with `src.data`.

Measured at 100,000,000 nonzeros by `benchmarks/copy_avoidance_probe.py`:

| Site | Copy | Shared indices | `data`, `indices`, `indptr` |
| --- | ---: | ---: | --- |
| t-test squaring | 1200.4 MB, 0.11 s | 400.0 MB, 0.04 s | byte-equal |
| tf-idf indicator | 1200.4 MB, 0.11 s | 400.0 MB, 0.03 s | byte-equal |

**B. Stop returning a view when there is nothing to drop.**
`check_filter_single_cluster` removes clusters holding one cell. When none
exist, upstream still returned an AnnData view, and a view rebuilds `X` on every
attribute access: measured 81.1 MB on the first `.X` and 81.0 MB on the second,
against an 80.2 MB matrix, with no caching. `each_key_run` reaches `.X` from
`marker_gene`, `reassign_score`, `tf_idf10_for_cluster` and `SCCAF_score`, so
the view cost one full rebuild per reader.

The view was also providing accidental write-protection. `tests/
test_no_write_to_adata.py` replaces that guarantee: it checksums `X`, `obs`,
`var` and every layer across a full `each_key_run` and fails if one byte moves.
It runs both with and without singleton clusters, so it covers both branches.

**C. Free layers that `--layer` does not name.** `sc.read` materialises every
layer, and nothing reads the ones `--layer` did not name. Deleting them does not
lower the read peak; it lowers everything after it, which is where the real peak
is. It also removes them from the written h5ad, so `--keep-layers` opts out.
Verified on a 400 by 3,000 file carrying `counts` and `spliced`: the default
freed both and the output h5ad had none, `--keep-layers` kept both, and
`--layer counts` freed only `spliced`. `run_parameters.json` records which.

### Do the labels move? No.

`benchmarks/three_changes_validate.py` runs the same pipeline twice in one
interpreter, restoring the old code by monkeypatch for the first run:

| Dataset | Peak, old | Peak, new | `raw` | `pruned` | `final_annotation` | `confidence` |
| --- | ---: | ---: | :-: | :-: | :-: | :-: |
| demo, 2,700 cells | 269.3 MB | 192.8 MB | 2700/2700 | 2700/2700 | 2700/2700 | 2700/2700 |
| simulated, 8,000 cells | 520.5 MB | 193.3 MB | 8000/8000 | 8000/8000 | 8000/8000 | 8000/8000 |

1.40x less on the demo and 2.69x less on the simulation, with every cell
identical on all four label columns. Peak here is tracemalloc's python
allocation total, so read it as the difference between the two runs, not as the
process high-water mark.

**On time: no measurable change, and my first measurement was wrong.** That run
reported the new code slower, 18.3 s against 25.8 s. tracemalloc was on, and its
overhead, not the change, produced the gap.
`benchmarks/three_changes_ablation.py` runs all four combinations of A and B,
twice each, in both orders, with tracemalloc off:

| A | B | Run 1 | Run 2 | Median |
| --- | --- | ---: | ---: | ---: |
| old | old | 11.6 s | 10.6 s | 11.1 s |
| old | new | 14.3 s | 11.8 s | 13.0 s |
| new | old | 11.2 s | 12.3 s | 11.8 s |
| new | new | 11.5 s | 12.4 s | 11.9 s |

Every combination sits inside the 2.5 s spread that a single combination shows
between its own two runs, so the demo cannot resolve a difference here. All 8
runs returned identical `final_annotation` on 2,700 of 2,700 cells. The copy
sites themselves are about 3x faster in isolation; they are too small a share of
a demo run to move the total.

### The seven I did not do

| Proposal | Why not |
| --- | --- |
| Skip the sparse sidecar | Refuted. `layers[LAYER] is X` and shares memory. Pickle memoises it: 16,082,033 bytes became 16,082,043, +0.0%. A real copy in that slot costs +99.6%. |
| Sparse one-hot as a memory fix | Backwards. It is 4.8x faster and bit-identical, but matmul peak rose 819.2 MB to 1629.0 MB, because scipy pre-allocates the sparse product from an upper bound. The one-hot array shrinks; the operation does not. |
| Chunk the two copy sites | Right problem, wrong fix. Chunking changes float64 addition order and risks the bit-identity against scanpy. Sharing the indices gets the same saving with no arithmetic change. |
| SCCAF: rows before columns | Real but small. `_sccaf_optimized` holds one intermediate, so about one matrix, not the 4x implied. |
| Subset cells at read time | The most valuable idea on the list for a 984k object, but it cites `scripts/01_build_annotated_h5ad.py`, which does not exist in this repository. There is no block reader to reuse. |
| int32 indices | scipy already picks int32 when it fits. The proposal's own figure implies about 2.25 billion nonzeros, above the 2,147,483,648 that int32 can address, so int32 cannot apply there. |
| Cap workers by memory | Sensible, changes no result. `--cores` is already the lever; a memory-aware default needs a per-worker footprint that has to be measured. |

Also found while measuring, not fixed: `reassign_score` asks PCA for 30
components without checking that the marker pool holds 30 genes, and raises
`n_components=30 must be between 0 and min(n_samples, n_features)=18` when it
does not. `reassign_pruning` guards the same situation for its KNN neighbour
count. Small structured datasets hit this.

## Two optional reductions for very large inputs

Both are off. Both change what the metrics see. Neither one earns its cost at the sizes
measured here, because the sparse path already removed the dense matrix that set the limit.
Turn one on when a run does not fit, not to make a run faster.

### `--downsample N`: at most N cells per cluster, one shared whitelist

The naive way to cap cells is to draw N per cluster for each annotation separately. That
double-pays, because annotations that agree name the same cells. This draws once into a
shared whitelist:

1. The first `--query` annotation draws up to N cells per cluster.
2. The second annotation counts, for each of its own clusters, how many cells the whitelist
   already holds. A cluster holding 980 of them draws 20 more, and it draws them from the
   cells the whitelist does not hold. A cluster already at N draws nothing.
3. Every later annotation does the same.

Every cluster of every annotation still reaches `min(cluster size, N)` cells, so no
annotation loses a cluster and the Shapley step still sees the same competitors. Verified as
`min_cluster_coverage = 1.000` at every point in the table below.

What the sharing saves, measured by `benchmarks/downsample_validate.py`:

| Dataset | Cap | Shared whitelist | Independent per annotation | Naive quota sum |
| --- | ---: | ---: | ---: | ---: |
| demo, 2,700 cells, 3 annotations, 9/19/33 clusters | 50 | 1,437 | 1,806 | 2,617 |
| | 100 | 2,010 | 2,367 | 4,156 |
| | 250 | 2,610 | 2,693 | 6,322 |
| | 500 | 2,700 | 2,700 | 7,426 |
| simulated, 100,000 cells, 3 annotations, 120 clusters | 100 | 12,000 | 22,332 | 24,000 |
| | 250 | 30,000 | 49,540 | 60,000 |
| | 500 | 60,000 | 80,808 | 120,000 |
| | 1,000 | 100,000 | 100,000 | 213,222 |

On the simulation the shared whitelist is 1.86x smaller than the independent draw at cap 100
and 1.35x smaller at cap 500. The last row is the honest warning: at cap 1,000 every cluster
of that simulation is already under the cap, so the flag keeps all 100,000 cells and saves
nothing. **Whether `--downsample 1000` reduces anything depends on your cluster sizes, not on
your cell count.** Check `n_cells_analysed` in `run_parameters.json`.

Does it change the answer? The demo, full pipeline, labels compared on the cells both runs
analysed:

| Cap | Cells kept | ARI raw | ARI pruned | ARI final | Exact final | Pruned clusters |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| none | 2,700 | — | — | — | — | 27 |
| 1,000 | 2,700 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 27 |
| 500 | 2,700 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 27 |
| 250 | 2,610 | 0.9749 | 0.9536 | 0.7176 | 0.8747 | 26 |
| 100 | 2,010 | 0.9326 | 0.8265 | 0.7564 | 0.9214 | 28 |

The two rows that keep every cell reproduce the full run exactly, on all 2,700 cells, for
`raw`, `pruned` and `final_annotation`. Those two rows are the control: the machinery adds
no drift of its own, so everything below them measures the cost of the missing cells. Dropping 90 of 2,700 cells
(3.3%) already moved `final_annotation` to ARI 0.7176. Cell counts are the denominator of
every stability metric, so the drop follows from the method and is not a defect.

I did not report seconds for these runs. Two runs on the identical 2,700 cells differed by
2x (71.2 s and 35.0 s), so the timings are warm-up noise, not measurement.

### `--prefilter-markers N`: keep the top N MarkerFinder genes per cluster

Scores genes against cluster membership with
`altanalyze3.components.cellHarmony.markerFinder.marker_finder`, which is already sparse end
to end, on at most `--prefilter-cells` cells per cluster, then keeps the union of the top N
per cluster over every annotation.

On the demo, measured by `benchmarks/prefilter_validate.py`:

| N | Genes kept of 32,738 | ARI raw | ARI pruned | ARI final |
| ---: | ---: | ---: | ---: | ---: |
| 500 | 11,938 | 0.9237 | 0.8870 | 0.6123 |
| 200 | 6,380 | 0.7888 | 0.7712 | 0.1190 |
| 50 | 1,969 | 0.8865 | 0.8485 | 0.3531 |

**Do not turn this on without reading that table.** Two things are wrong with it. The
selection is circular: it picks genes using the same labels it then scores. And ARI is not
monotonic in N, which reads as instability rather than graceful degradation. The sparse path
delivered the memory saving without either problem, so this flag exists as a last resort for
an object that still does not fit.

### Which to reach for

The sparse compute path is the default and costs nothing: on the 100,000 by 20,000 simulation
with 120 clusters it took peak memory from 28.3 GB to 10.7 GB and the labels were identical.
Use `--sequential` next, which trades about 9 s for roughly half the peak. Only then reach for
`--downsample`, and only then for `--prefilter-markers`. Whichever you use, say so in the
methods. `run_parameters.json` records both, under `downsample_info` and
`prefilter_info`.

**I have not run either flag at 984,119 cells by 36,249 genes.** The largest input measured
here is 100,000 by 20,000.

## Reproducing every number

```
cd /Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/sctriangulate/benchmarks
PY=/opt/homebrew/opt/python@3.11/bin/python3.11

$PY baseline_profile.py  data/demo_pbmc3k.h5ad out_baseline    # upstream, staged
$PY optimized_profile.py data/demo_pbmc3k.h5ad out_optimized   # this version, staged
$PY compare_outputs.py out_baseline/baseline_reference.json \
                       out_optimized/optimized_reference.json
$PY determinism_audit.py data/demo_pbmc3k.h5ad                 # the noise floor
$PY shapley_scaling.py --cells 20000                           # Shapley vs annotation count
$PY end_to_end.py data/demo_pbmc3k.h5ad                        # lazy_run, both versions
$PY label_agreement.py                                         # exact, ARI and NMI
$PY sccaf_parallel_probe.py                                    # the failed SCCAF attempt
$PY simulate_dataset.py data/sim_8k_20k.h5ad --cells 8000 --genes 20000 \
       --gene-names-from data/demo_pbmc3k.h5ad                 # larger simulated input
$PY sparse_validate.py                                         # sparse against dense
$PY prefilter_validate.py                                      # --prefilter-markers effect
$PY downsample_validate.py                                     # --downsample effect
$PY copy_avoidance_probe.py                                    # the three proposed copy fixes
$PY three_changes_validate.py                                  # old against new: labels, peak
$PY three_changes_ablation.py                                  # which change moved the time
```

Which script produced which claim:

| Claim in this file | Script | Result file |
| --- | --- | --- |
| 7.35x end to end, upstream against this version | `end_to_end.py` | printed |
| Staged stage timings, 135.14 s to 19.07 s | `baseline_profile.py`, `optimized_profile.py`, `compare_outputs.py` | `out_baseline/`, `out_optimized/` |
| Upstream's run-to-run spread, 1 to 3 of 50 | `determinism_audit.py` | `out_determinism_run.log` |
| Shapley 205x to 138,563x | `shapley_scaling.py` | printed |
| Exact agreement, ARI, NMI against upstream | `label_agreement.py` | printed |
| Peak memory, 4 configurations | `memory_profile.py` | `out_mem_*/` |
| Sparse against dense, ARI 1.0000 | `sparse_validate.py` | printed |
| `--prefilter-markers` gene counts and ARI | `prefilter_validate.py` | `prefilter_validate.json` |
| `--downsample` whitelist sizes and ARI | `downsample_validate.py` | `downsample_validate.json` |
| Copy against shared indices, 1200.4 to 400.0 MB | `copy_avoidance_probe.py` | `copy_avoidance_probe.json` |
| Three memory changes: labels identical, peak | `three_changes_validate.py` | `three_changes_validate.json` |
| No resolvable time change from those | `three_changes_ablation.py` | `three_changes_ablation.json` |
| The two SCCAF attempts that failed | `sccaf_parallel_probe.py`, `sccaf_compare.py` | printed |

`benchmarks/.gitignore` keeps the 23 MB demo file and every output directory out of git.

## Tests

```
cd /Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/sctriangulate
/opt/homebrew/opt/python@3.11/bin/python3.11 -m pytest tests/ -q
```

73 pass, in about 75 s.

| File | Tests | What it holds to account |
| --- | ---: | --- |
| `test_shapley_equivalence.py` | 29 | The closed-form Shapley engine against upstream's coalition enumeration, in `_reference/upstream_shapley.py` |
| `test_metrics_equivalence.py` | 14 | Marker assignment, tf-idf, enrichr and GSEA against `_reference/upstream_metrics.py` |
| `test_sccaf_modes.py` | 8 | `optimized` against `legacy`, and that each repeats itself |
| `test_ttest_equivalence.py` | 6 | The fast group statistics against scanpy, bit for bit |
| `test_no_write_to_adata.py` | 6 | That `each_key_run` writes nothing to adata, and that `_shared_index_csr` is byte-equal to the copy it replaces |
| `test_downsample.py` | 10 | The shared whitelist against a plain reference implementation, and the shortfall arithmetic |

## What I did not do

- I did not speed up `SCCAF_score`. I tried two routes and measured both; each changed the
  metric for a gain of 1.03x to 1.56x. See "I tried to break the SCCAF ceiling, and failed".
  I removed one dense copy.
- I did not touch `reference_pruning` in `prune.py`. `lazy_run` never calls it. It still starts
  one process per reference cluster and passes the whole obs table to each.
- I did not rewrite gseapy's GSEA permutation. Reproducing its RNG stream would be fragile.
- I did not test the spatial module. It needs squidpy 1.2.0.
- I did not test the viewer HTML output.
- **I did not run anything at 984,119 cells by 36,249 genes.** The largest input measured
  here is 100,000 by 20,000. Every statement about the full object is arithmetic from that,
  and says so where it appears.
- I did not remove every dense component. Still dense, and all outside the metric path:
  `_sccaf_legacy`, which is dense by definition and opt-in; the heatmap helpers in
  `main_class.py` (`plot_heterogeneity` near lines 2104 and 2273, `plot_long_heatmap` near
  2864 and 2905), all calling `make_sure_mat_dense`; the utilities in `preprocessing.py`
  at lines 192, 357, 548, 681, 804, 836 and 1311; and `spatial.py:418`.
- I did not fix a latent crash I found while testing: `reassign_score` asks PCA for 30
  components without checking that the marker pool holds 30 genes, and raises
  `n_components=30 must be between 0 and min(n_samples, n_features)=18` when it does not.
  `reassign_pruning` guards the same situation for its KNN neighbour count. Small or
  narrowly structured datasets hit this.
- I did not add scTriangulate to the Sphinx site. `docs/index.rst` has no entry for it, and
  `docs/` documents no other component this way except cellHarmony. This README is the
  documentation of record.
- I did not commit anything to git.

## Change log

Each entry names the script that measured it. Dates give the day the work landed.

| Date | Change | Effect | Evidence |
| --- | --- | --- | --- |
| 2026-08-06 | Integrated upstream 0.13.0; closed-form Shapley; vectorised marker assignment, tf-idf, enrichr, GSEA; macOS import and fork fixes | 7.35x end to end; Shapley 205x to 138,563x | `end_to_end.py`, `shapley_scaling.py`, 49 tests |
| 2026-08-06 | Seeded the three unseeded sources; replaced the randomized reassign SVD with an exact one | Upstream matched its own first run in 1 to 3 of 50 repeats; this version 50 of 50 | `determinism_audit.py` |
| 2026-08-07 | SCCAF split into `optimized` (default) and `legacy` | Same accuracies, sparse and seeded | `sccaf_compare.py`, `test_sccaf_modes.py` |
| 2026-08-07 | Sparse compute made the default | Demo 2.71 to 1.11 GB, sim 100k 28.3 to 10.7 GB, ARI 1.0000 | `sparse_validate.py` |
| 2026-08-07 | `--prefilter-markers`, as a Python method only | 5.1x fewer genes, but ARI 0.79 on raw labels | `prefilter_validate.py` |
| 2026-08-08 | `--prefilter-markers` wired to the CLI; `--downsample`, `--prefilter-cells`, `--subsample-seed` added | Shared whitelist 1.86x smaller than independent draws at 100k cells; cap above every cluster reproduces the full run exactly | `downsample_validate.py`, `test_downsample.py` |
| 2026-08-08 | Reviewed ten proposed memory changes; applied three, refuted two, deferred two | Demo peak 269.3 to 192.8 MB, sim 8k 520.5 to 193.3 MB, every label identical | `copy_avoidance_probe.py`, `three_changes_validate.py`, `three_changes_ablation.py` |
| 2026-08-08 | `--keep-layers`; unused layers freed by default | Verified three ways on a file carrying `counts` and `spliced` | printed CLI output, `run_parameters.json` |
| 2026-08-17 | Steps 4 to 6 (MarkerFinder, cell-state enrichment naming, HOPACH order) added as a post-pruning default; `--no-annotate` restores the old behaviour | Reproduces a hand run exactly: 1,749 of 1,749 markers identical, HOPACH k=11 MSS 0.3949, 13,210 of 13,210 pruned labels identical | `27_validate_annotate_default.py`, `28_validate_annotate_default_paths.py`, `29_e2e_scTriangulate_annotate_default.sh`, `test_annotate.py` |
