"""Reconciled isoform ratios use the pinned compatibility gene denominator."""
import numpy as np
import pandas as pd
from scipy import sparse


def anchored_ratio(integrated, baseline, features, cpm_threshold=1):
    from altanalyze3.components.long_read.isoform_matrix import cpm_and_ratio_from_counts
    x=sparse.csr_matrix(integrated,dtype=np.float64)
    b=sparse.csr_matrix(baseline,dtype=np.float64)
    if x.shape!=b.shape or x.shape[1]!=len(features):
        raise ValueError('Complete aligned feature and pseudobulk universes required')
    if not np.isfinite(x.data).all() or (x.data<0).any():
        raise ValueError('Invalid integrated counts')
    difference=b-x
    if difference.nnz and (difference.data<0).any():
        raise ValueError('Integrated additions lack an authorized additions ledger')
    _,_,expressed,genes=cpm_and_ratio_from_counts(b,features,cpm_threshold=cpm_threshold)
    codes=pd.Index(genes).get_indexer([str(f).split(':')[0] for f in features])
    mapping=sparse.csr_matrix((np.ones(len(features)),(np.arange(len(features)),codes)),shape=(len(features),len(genes)))
    denominators=(b@mapping).tocsr()
    coo=x.tocoo();on=np.asarray(expressed[coo.row,codes[coo.col]]).ravel()>0
    d=np.asarray(denominators[coo.row,codes[coo.col]]).ravel()
    keep=on&(d>0)&(coo.data>0)
    ratio=sparse.csr_matrix(((coo.data[keep]/d[keep]).astype(np.float32),(coo.row[keep],coo.col[keep])),shape=x.shape)
    return ratio,expressed,genes


def write_anchored_ratios(out,baseline):
    import anndata as ad
    from .manifest import digest
    name='isoform_combined_pseudo_cluster_counts.h5ad'
    x=ad.read_h5ad(out/name);b=ad.read_h5ad(baseline/name)
    if not x.obs_names.equals(b.obs_names) or not x.var_names.equals(b.var_names):
        raise ValueError('Integrated pseudobulk universe differs from complete baseline; no intersection allowed')
    ratio,mask,genes=anchored_ratio(x.X,b.X,list(x.var_names))
    result=ad.AnnData(ratio,obs=x.obs.copy(),var=x.var.copy(),obsm={'ratio_gene_expressed':mask},
        uns={'ratio_genes':genes,'ratio_cpm_threshold':1.,'ratio_denominator_view':'compatibility gene totals',
             'ratio_denominator_counts_sha256':digest(baseline/name),
             'ratio_encoding':'numerator from integrated standard assignments; gene denominator and expression mask frozen to compatibility'})
    result.write_h5ad(out/'isoform_combined_pseudo_cluster_ratio.h5ad',compression='gzip')
    # Preserve the pre-existing filtered comparison universe, including newly zero features.
    prior=ad.read_h5ad(baseline/'isoform_combined_pseudo_cluster_ratio-filtered.h5ad')
    result[prior.obs_names,prior.var_names].copy().write_h5ad(out/'isoform_combined_pseudo_cluster_ratio-filtered.h5ad',compression='gzip')
