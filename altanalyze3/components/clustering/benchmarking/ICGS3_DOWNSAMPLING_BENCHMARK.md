# ICGS3 Downsampling Benchmark: Adams HLCA Rare Cell Preservation

Date: 2026-08-15

## Goal

Evaluate whether 5,000-cell downsampling can preserve rare and transitional cell
states from the Adams HLCA lung dataset while minimizing memory and runtime.

Dataset:

`/Users/saljh8/Dropbox/Transfer/Adams_HLCA_with_umap_and_markers.h5ad`

Benchmark script:

`altanalyze3/components/clustering/benchmark_downsampling.py`

Primary output directories:

- `/private/tmp/icgs3_Adams_downsampling_benchmark_20260815_v2`
- `/private/tmp/icgs3_Adams_downsampling_benchmark_fullgene_20260815`
- `/private/tmp/icgs3_Adams_downsampling_benchmark_aggressive_20260815`

## Shared Preprocessing

All methods used the same ICGS3 preprocessing:

- Input cells: 280,446
- Input genes: 33,234
- QC:
  - `min_genes=500`
  - `min_cells=5`
  - `min_counts=1000`
  - `mito_percent=30`
- Post-QC cells: 280,446
- Post-QC genes: 32,343
- RNA unsupervised gene filter:
  - protein-coding lookup when available
  - removed RPL/RPS, MT-, dotted, GM, XIS/TSI, RSP, HLA, and `*Y`
- Filtered unsupervised genes: 18,100

Preprocessing time was about 101-106 seconds.

Memory note: peak RSS is process-level peak including the loaded/preprocessed
AnnData object. Method timings below exclude shared preprocessing. The graph
methods materially increase runtime and memory pressure beyond preprocessing;
the sparse full-gene methods complete in seconds after preprocessing.

## Methods Tested

| Method | Concept |
|---|---|
| `icgs2_louvain_pagerank` | Current ICGS2-like protocol: top 500 dispersion genes, Annoy k=10, Louvain pre-reduction to `target*4`, then PageRank to target. |
| `louvain_medoid_direct` | Same top-500 graph and Louvain communities, but directly select medoid-proximal community representatives to 5k; skips PageRank. |
| `rarity_quantile` | Top-500 dispersion genes, sparse inverse-detection cell rarity score, quantile-balanced sampling. |
| `all_gene_rarity_quantile` | Full 18,100 filtered genes, sparse inverse-detection cell rarity score, quantile-balanced sampling. |
| `random_projection_kmeans` | Top-500 genes, sparse random projection to 32D, MiniBatchKMeans landmarks. |
| `sparse_gene_coverage` | Full filtered genes; each gene nominates top expressing cells weighted by inverse detection frequency; 80% coverage-selected, 20% random reserve. |
| `sparse_gene_coverage_aggressive` | Full filtered genes; same as above with more nominations per gene and 90% coverage-selected cells. |
| `random` | Uniform random baseline. |

## Summary Results

Target cells: 5,000.

Rare population definition: HLCA categories with `20 < cells < 1000`.

| Method | Sampled cells | Time sec | Peak RSS MB | Rare states | States <20 retained | Rare-cell capture | Mean rare fraction | Median rare fraction | Ionocytes retained |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `sparse_gene_coverage` | 5,000 | 5.7 | 30,850 | 30 | 16 | 0.0556 | 0.0913 | 0.0468 | 9 / 22 |
| `sparse_gene_coverage_aggressive` | 5,000 | 7.1 | 32,694 | 30 | 17 | 0.0682 | 0.1226 | 0.0570 | 15 / 22 |
| `random_projection_kmeans` | 5,000 | 62.9 | 30,777 | 30 | 28 | 0.0215 | 0.0252 | 0.0228 | 2 / 22 |
| `louvain_medoid_direct` | 4,114 | 633.6 | 30,777 | 30 | 28 | 0.0115 | 0.0172 | 0.0123 | 2 / 22 |
| `random` | 5,000 | 0.3 | 30,777 | 30 | 28 | 0.0196 | 0.0187 | 0.0203 | 0 / 22 |
| `all_gene_rarity_quantile` | 5,000 | 2.8 | 30,850 | 30 | 28 | 0.0175 | 0.0172 | 0.0168 | 0 / 22 |
| `icgs2_louvain_pagerank` | 5,000 | 739.2 | 30,777 | 30 | 28 | 0.0152 | 0.0130 | 0.0060 | 0 / 22 |
| `rarity_quantile` | 5,000 | 42.2 | 30,777 | 30 | 30 | 0.0187 | 0.0196 | 0.0184 | 0 / 22 |

## Interpretation

The current ICGS2-like 5k downsampling is not sufficient for this Adams dataset.
It retained 0/22 ionocytes and retained fewer than 20 cells for 28/30 rare HLCA
states. The failure is consistent with the earlier 30k diagnostic: the top-500
dispersion feature set contains only one of the checked ionocyte markers
(`ASCL3`), so the graph is weak for ionocyte preservation before NMF/SVM.

The best tested non-label-aware alternative was
`sparse_gene_coverage_aggressive`. It retained 15/22 ionocytes and had the
highest total rare-cell capture fraction, while requiring only 7.1 seconds after
preprocessing. This method is mathematically distinct from graph sampling: it is
a sparse set-coverage sketch where rare/high-expression features nominate cells.
It therefore directly addresses the failure mode caused by top-500 dispersion
features excluding rare-state markers.

The less aggressive `sparse_gene_coverage` had slightly fewer rare states below
20 cells, but retained only 9/22 ionocytes. It may be less biased, but for the
ionocyte use case the aggressive variant is better.

Random projection k-means was faster than graph methods but did not preserve
rare states well enough. Random and rarity quantile methods were not adequate.

## Recommendation

Do not use 5k ICGS2-like Louvain/PageRank alone as the default for very large
RNA datasets where rare-state discovery is required.

The most defensible next candidate is a hybrid:

1. Apply official RNA filtering.
2. Build a 5k sparse gene-coverage sketch using the full filtered sparse matrix.
3. Optionally union this with a smaller ICGS2 PageRank/Louvain sketch when the
   user requests graph/geometric representation.
4. Deduplicate and cap to the target.
5. Run UDON feature selection and NMF on that retained set.

This keeps the expensive graph optional and uses sparse operations for the
rare-state guard. It is not yet an official replacement for ICGS2 PageRank; it
should next be tested through NMF/MarkerFinder/SVM to confirm that retained rare
cells become stable marker-defined clusters rather than only retained barcodes.

## Reproduction

Current ICGS2-like and initial alternatives:

```bash
python3.11 components/clustering/benchmark_downsampling.py \
  --input /Users/saljh8/Dropbox/Transfer/Adams_HLCA_with_umap_and_markers.h5ad \
  --output-dir /private/tmp/icgs3_Adams_downsampling_benchmark_20260815_v2 \
  --input-normalized \
  --target-cells 5000 \
  --species Hs \
  --modality rna \
  --audit-obs HLCA \
  --methods icgs2_louvain_pagerank,louvain_medoid_direct,rarity_quantile,random_projection_kmeans,random
```

Full-gene sparse alternatives:

```bash
python3.11 components/clustering/benchmark_downsampling.py \
  --input /Users/saljh8/Dropbox/Transfer/Adams_HLCA_with_umap_and_markers.h5ad \
  --output-dir /private/tmp/icgs3_Adams_downsampling_benchmark_fullgene_20260815 \
  --input-normalized \
  --target-cells 5000 \
  --species Hs \
  --modality rna \
  --audit-obs HLCA \
  --methods all_gene_rarity_quantile,sparse_gene_coverage
```

Aggressive sparse gene-coverage:

```bash
python3.11 components/clustering/benchmark_downsampling.py \
  --input /Users/saljh8/Dropbox/Transfer/Adams_HLCA_with_umap_and_markers.h5ad \
  --output-dir /private/tmp/icgs3_Adams_downsampling_benchmark_aggressive_20260815 \
  --input-normalized \
  --target-cells 5000 \
  --species Hs \
  --modality rna \
  --audit-obs HLCA \
  --methods sparse_gene_coverage_aggressive
```
