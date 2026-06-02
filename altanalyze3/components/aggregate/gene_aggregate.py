#!/usr/bin/env python3
"""Fast molecule-level -> gene-level aggregation for BAM-derived isoform h5ad files.

The BAM extraction (``isoform_structure_extract``) produces a per-sample h5ad whose columns are
individual molecules/isoform-structures (``var_names`` like ``ENSG00000026025:344``) and whose
``var['gene']`` column already carries the Ensembl gene id for every feature. There can be millions
of such columns, so this aggregation is done with a single sparse matrix multiply rather than any
per-gene Python loop:

    X_gene[cell, gene] = sum over molecules m with var['gene'][m] == gene  of  X[cell, m]

i.e. ``X_gene = X @ G`` where ``G`` is a (n_molecules x n_genes) 0/1 grouping matrix. This is
O(nnz) and keeps only the (cells x genes) result in memory.

The result (cells x genes raw counts) is written as ``<library>-gene.h5ad`` next to the input and
the **path** is returned, so a master script can hand it straight to cellHarmony without holding the
AnnData in memory.
"""

from __future__ import annotations

import os
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def aggregate_to_gene_level(h5ad_path, output_path=None, gene_col="gene", log=print):
    """Aggregate a molecule-level h5ad to gene level. Returns the output h5ad path.

    Parameters
    ----------
    h5ad_path : str
        BAM-derived molecule-level h5ad (var_names ``gene:molecule``, ``var[gene_col]`` = gene id).
    output_path : str, optional
        Where to write the gene-level h5ad. Defaults to ``<input_dir>/<library>-gene.h5ad`` where
        ``<library>`` is the input basename with the ``.h5ad`` stripped.
    gene_col : str
        ``var`` column holding the gene id (default ``'gene'``). If absent, the gene id is parsed
        from ``var_names`` as the substring before the first ``':'``.

    Returns
    -------
    str
        Path to the written gene-level h5ad (cells x genes, raw summed counts in ``X``).
    """
    t0 = time.time()
    a = ad.read_h5ad(h5ad_path)
    n_cells, n_mol = a.shape
    log(f"[gene-agg] loaded {os.path.basename(h5ad_path)}: {n_cells:,} cells x {n_mol:,} molecules "
        f"({time.time()-t0:.1f}s)")

    # gene id per molecule: prefer the var column written by the BAM extractor, else parse var_names.
    if gene_col in a.var.columns:
        genes = a.var[gene_col].astype(str).to_numpy()
    else:
        genes = np.array([str(v).split(":", 1)[0] for v in a.var_names])

    # map each molecule -> its gene's column index (stable, first-seen order)
    gene_ids, inv = np.unique(genes, return_inverse=True)   # gene_ids sorted; inv[m] = gene index
    n_genes = gene_ids.shape[0]

    # grouping matrix G (n_mol x n_genes); X_gene = X @ G
    G = sp.csr_matrix(
        (np.ones(n_mol, dtype=np.int64), (np.arange(n_mol, dtype=np.int64), inv.astype(np.int64))),
        shape=(n_mol, n_genes),
    )
    X = a.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X_gene = (X.astype(np.int64) @ G).tocsr()
    log(f"[gene-agg] aggregated to {n_genes:,} genes; total counts "
        f"{int(X.sum()):,} -> {int(X_gene.sum()):,} ({time.time()-t0:.1f}s)")

    out = ad.AnnData(
        X=X_gene,
        obs=a.obs.copy(),
        var=pd.DataFrame(index=pd.Index(gene_ids, name=None)),
    )

    if output_path is None:
        base = os.path.basename(h5ad_path).split(".h5ad")[0]
        output_path = os.path.join(os.path.dirname(h5ad_path) or ".", f"{base}-gene.h5ad")
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    out.write_h5ad(output_path, compression="gzip")
    log(f"[gene-agg] wrote {output_path} ({n_cells:,} cells x {n_genes:,} genes) "
        f"in {time.time()-t0:.1f}s")
    return output_path


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Aggregate a molecule-level h5ad to gene level.")
    p.add_argument("--h5ad", required=True, help="Input molecule-level h5ad.")
    p.add_argument("--output-h5ad", default=None, help="Output gene-level h5ad (default <library>-gene.h5ad).")
    p.add_argument("--gene-col", default="gene", help="var column with the gene id (default 'gene').")
    args = p.parse_args()
    print(aggregate_to_gene_level(args.h5ad, args.output_h5ad, args.gene_col))
