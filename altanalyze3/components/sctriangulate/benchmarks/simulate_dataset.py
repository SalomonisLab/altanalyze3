"""Write a simulated h5ad for scaling tests.

The demo holds 2,700 cells. Several of the rewritten steps carry costs that grow
faster than linearly in cells or genes, so a single small dataset understates
them. Two examples: the upstream gene-to-cluster loop costs O(genes**2 * clusters),
and the upstream tf-idf builds a dense cells-by-genes DataFrame.

The simulation aims to be realistic where the runtime depends on it, not to model
biology:
  * counts are negative binomial, then CP10k-normalised and log1p-transformed,
    which is the state the demo file is already in;
  * roughly 90 percent of entries are zero, close to the demo's 97.4 percent;
  * every cluster gets a set of up-regulated marker genes, so the t-test finds
    real signal and the marker lists are not empty;
  * the three annotations are nested refinements of one another, which is what
    scTriangulate is built to reconcile.

Usage:
  python3.11 simulate_dataset.py <out.h5ad> --cells 10000 --genes 15000 --seed 0
"""

import argparse

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix, vstack


def simulate(n_cells, n_genes, cluster_counts=(10, 20, 30), seed=0,
             n_markers=40, marker_boost=6.0, zero_fraction=0.90,
             gene_names=None):
    rng = np.random.default_rng(seed)

    finest = max(cluster_counts)
    base_labels = rng.integers(0, finest, n_cells)

    # baseline per-gene expression, heavy tailed like real counts
    base_rate = rng.lognormal(mean=0.0, sigma=1.2, size=n_genes)

    rate = np.tile(base_rate, (finest, 1))                      # (clusters, genes)
    for c in range(finest):
        markers = rng.choice(n_genes, size=n_markers, replace=False)
        rate[c, markers] *= marker_boost

    # Generate in row blocks and keep only the sparse result. A 100,000 by 20,000
    # Poisson draw in float64 is 16 GB; the CSR of the same data is about 1.4 GB.
    block_rows = max(1, int(64 * 1024 * 1024 / (n_genes * 8)))
    blocks = []
    for start in range(0, n_cells, block_rows):
        stop = min(start + block_rows, n_cells)
        counts = rng.poisson(rate[base_labels[start:stop]]).astype(np.float32)
        counts *= (rng.random(counts.shape) > zero_fraction)
        total = counts.sum(axis=1, keepdims=True)
        total[total == 0] = 1.0
        blocks.append(csr_matrix(np.log1p(counts / total * 1e4).astype(np.float32)))
        del counts, total
    X = vstack(blocks, format='csr') if len(blocks) > 1 else blocks[0]
    del blocks

    obs = pd.DataFrame(index=['cell{}'.format(i) for i in range(n_cells)])
    for n_clusters in cluster_counts:
        # nested refinement: coarser annotations merge neighbouring fine clusters
        mapped = (base_labels * n_clusters) // finest
        obs['sctri_sim_leiden_{}'.format(n_clusters)] = \
            pd.Categorical(mapped.astype(str))

    if gene_names is not None:
        names = list(gene_names)[:n_genes]
        if len(names) < n_genes:
            names += ['SIMGENE{:05d}'.format(i) for i in range(n_genes - len(names))]
    else:
        names = ['SIMGENE{:05d}'.format(i) for i in range(n_genes)]
    var = pd.DataFrame(index=names)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm['X_umap'] = rng.normal(size=(n_cells, 2)).astype(np.float32)
    return adata


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--cells', type=int, default=10000)
    ap.add_argument('--genes', type=int, default=15000)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--clusters', default='10,20,30')
    ap.add_argument('--gene-names-from', default=None,
                    help='h5ad to borrow real gene symbols from, so the artifact gene '
                         'overlap (and therefore the enrichment cost) is realistic')
    args = ap.parse_args()

    gene_names = None
    if args.gene_names_from:
        import scanpy as sc
        gene_names = list(sc.read(args.gene_names_from).var_names)

    cluster_counts = tuple(int(x) for x in args.clusters.split(','))
    adata = simulate(args.cells, args.genes, cluster_counts, args.seed,
                     gene_names=gene_names)
    adata.write(args.out)

    density = adata.X.nnz / (adata.shape[0] * adata.shape[1])
    print('wrote {}'.format(args.out))
    print('  {} cells x {} genes, {:.1f}% nonzero'.format(
        adata.shape[0], adata.shape[1], 100 * density))
    for c in adata.obs.columns:
        print('  {}: {} clusters'.format(c, adata.obs[c].nunique()))
    import sys, os
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                    '..', '..', '..', '..')))
    from altanalyze3.components.sctriangulate.metrics import _artifact_gene_set
    art = _artifact_gene_set('human', 2)
    n_art = sum(1 for g in adata.var_names if g in art)
    print('  {} of {} genes are human artifact genes (criterion 2)'.format(
        n_art, adata.shape[1]))


if __name__ == '__main__':
    main()
