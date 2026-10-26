"""Bulk deconvolution scores equal the single-cell alignment metric, state by state."""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from altanalyze3.components.cellHarmony.cellHarmony_lite import _cosine_best_matches, deconvolve_bulk


def _write_inputs(tmp_path, counts):
    rng=np.random.default_rng(7)
    genes=[f"G{i}" for i in range(60)]
    states=['HSC','MEP','Mono','B']
    ref=pd.DataFrame(rng.gamma(1.0,1.0,(60,4)),index=genes,columns=states)
    ref.index.name='UID'
    ref.to_csv(tmp_path/'ref.txt',sep='\t')
    # 10 bulk genes are absent from the reference and 5 reference genes absent from the bulk table
    bulk_genes=genes[5:]+[f"X{i}" for i in range(10)]
    table=pd.DataFrame(counts(rng,len(bulk_genes)),index=bulk_genes,columns=['S1','S2','S3'])
    table.index.name='UID'
    table.to_csv(tmp_path/'bulk.txt',sep='\t')
    return ref,table


def _expected_log(table, ref):
    obj=ad.AnnData(table.T.to_numpy(dtype=np.float64))
    sc.pp.normalize_total(obj,target_sum=1e4)
    sc.pp.log1p(obj)
    log=pd.DataFrame(obj.X,index=table.columns,columns=table.index)
    shared=[g for g in ref.index if g in log.columns]
    return log[shared],ref.loc[shared].T


def test_cosine_matches_alignment(tmp_path):
    ref,table=_write_inputs(tmp_path,lambda rng,n:rng.poisson(rng.gamma(2.0,50.0,(n,1)),(n,3)))
    scores=deconvolve_bulk(str(tmp_path/'bulk.txt'),str(tmp_path/'ref.txt'),output_dir=str(tmp_path/'out'))
    query,reference=_expected_log(table,ref)
    best,best_scores=_cosine_best_matches(query.values,reference.values)
    assert list(scores.index)==list(ref.columns) and list(scores.columns)==['S1','S2','S3']
    assert list(scores.idxmax(axis=0))==list(reference.index[best])
    np.testing.assert_allclose(scores.max(axis=0).values,best_scores,rtol=1e-10)
    written=pd.read_csv(tmp_path/'out'/'cellHarmony_deconvolution_scores.txt',sep='\t',index_col=0)
    np.testing.assert_allclose(written.values,scores.values,rtol=1e-12)


def test_classic_is_pearson(tmp_path):
    ref,table=_write_inputs(tmp_path,lambda rng,n:rng.poisson(rng.gamma(2.0,50.0,(n,1)),(n,3)))
    scores=deconvolve_bulk(str(tmp_path/'bulk.txt'),str(tmp_path/'ref.txt'),output_dir=str(tmp_path/'out'),
                           alignment_mode='classic')
    query,reference=_expected_log(table,ref)
    expected=np.corrcoef(query.values,reference.values)[:3,3:]
    np.testing.assert_allclose(scores.values,expected.T,rtol=1e-10)


def test_sample_without_signal_is_nan(tmp_path):
    def counts(rng,n):
        x=rng.poisson(50,(n,3))
        x[:-10,1]=0  # S2 keeps counts only on genes outside the reference
        return x
    _write_inputs(tmp_path,counts)
    scores=deconvolve_bulk(str(tmp_path/'bulk.txt'),str(tmp_path/'ref.txt'),output_dir=str(tmp_path/'out'))
    assert scores['S2'].isna().all()
    assert scores[['S1','S3']].notna().all().all()


def test_refuses_non_counts(tmp_path):
    _write_inputs(tmp_path,lambda rng,n:1e6*rng.dirichlet(np.ones(n),3).T)  # TPM-like
    with pytest.raises(ValueError,match='needs bulk counts'):
        deconvolve_bulk(str(tmp_path/'bulk.txt'),str(tmp_path/'ref.txt'),output_dir=str(tmp_path/'out'))


def test_refuses_annotation_columns_and_community(tmp_path):
    _,table=_write_inputs(tmp_path,lambda rng,n:rng.poisson(50,(n,3)))
    table.insert(0,'Chr','chr1')
    table.to_csv(tmp_path/'annotated.txt',sep='\t')
    with pytest.raises(ValueError,match='not numeric'):
        deconvolve_bulk(str(tmp_path/'annotated.txt'),str(tmp_path/'ref.txt'),output_dir=str(tmp_path/'out'))
    with pytest.raises(ValueError,match="'cosine' or 'classic'"):
        deconvolve_bulk(str(tmp_path/'bulk.txt'),str(tmp_path/'ref.txt'),output_dir=str(tmp_path/'out'),
                        alignment_mode='community')
