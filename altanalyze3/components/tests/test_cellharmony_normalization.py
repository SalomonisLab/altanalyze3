"""Chunked normalization preserves Scanpy values, counts and fallback behavior."""
import anndata as ad
import numpy as np
import pytest
import scanpy as sc
import scipy.sparse as sp
from altanalyze3.components.cellHarmony.cellHarmony_lite import normalize_adata


@pytest.mark.parametrize('dtype',[np.float32,np.float64,np.int32])
@pytest.mark.parametrize('format',['csr','csc','dense'])
def test_matches_scanpy(dtype,format):
    x=np.random.default_rng(4).integers(0,30,(25,40)).astype(dtype)
    x[3]=0
    original=x.copy()
    matrix={'csr':sp.csr_matrix,'csc':sp.csc_matrix,'dense':lambda v:v}[format](x)
    obj=ad.AnnData(matrix)
    obj.layers['counts']=matrix.copy()
    expected=obj.copy()
    sc.pp.normalize_total(expected,target_sum=1e4)
    sc.pp.log1p(expected)
    normalize_adata(obj)
    convert=lambda m:m.toarray() if sp.issparse(m) else m
    np.testing.assert_array_equal(convert(obj.X),convert(expected.X))
    np.testing.assert_array_equal(convert(obj.layers['counts']),original)
    assert obj.uns==expected.uns


def test_crosses_block_boundaries_without_copying_matrix():
    # Include zero-length rows around the 4M-nonzero boundary.
    n=4_100_000
    obj=ad.AnnData(sp.csr_matrix((np.ones(n,dtype=np.float32),np.tile(np.arange(1000,dtype=np.int32),4100),np.r_[0,0,np.arange(1000,n+1,1000)]),shape=(4101,1000)))
    expected=obj.copy()
    data=obj.X.data
    sc.pp.normalize_total(expected,target_sum=1e4)
    sc.pp.log1p(expected)
    normalize_adata(obj)
    assert np.shares_memory(data,obj.X.data)
    np.testing.assert_array_equal(obj.X.data,expected.X.data)


def test_view_preserves_parent():
    obj=ad.AnnData(sp.csr_matrix(np.arange(24,dtype=np.float32).reshape(6,4)))
    original=obj.X.copy()
    normalize_adata(obj[:3])
    np.testing.assert_array_equal(obj.X.toarray(),original.toarray())
