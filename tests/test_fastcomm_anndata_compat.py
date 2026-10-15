"""I/O compatibility only: identical matrices and split scores across API locations."""
import builtins
import sys
import types
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from altanalyze3.components.fastComm import api
from altanalyze3.components.fastComm.benchmark import FastCommBenchmarkParams,run_benchmark


def legacy_imports(monkeypatch):
    read_elem,sparse_dataset=api._anndata_readers()
    experimental=types.ModuleType('anndata.experimental')
    experimental.read_elem=read_elem;experimental.sparse_dataset=sparse_dataset
    monkeypatch.setitem(sys.modules,'anndata.experimental',experimental)
    original=builtins.__import__
    def import_legacy(name,*args,**kwargs):
        if name=='anndata.io':raise ModuleNotFoundError("No module named 'anndata.io'",name='anndata.io')
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',import_legacy)


def inputs(tmp_path,encoding):
    x=np.array([[5,0,0],[4,0,0],[0,3,2],[0,2,4]]*2,dtype=np.float32)
    matrix=x if encoding=='dense' else getattr(sparse,encoding+'_matrix')(x)
    obj=ad.AnnData(matrix,obs=pd.DataFrame({'state':['sender','sender','receiver','receiver']*2,
        'Library':['sample1']*4+['sample2']*4},index=[f'cell{i}' for i in range(8)]),
        var=pd.DataFrame({'symbol':['TGFB1','TGFBR1','TGFBR2']},index=['g1','g2','g3']))
    path=tmp_path/'counts.h5ad';obj.write_h5ad(path)
    lr=tmp_path/'lr.tsv';pd.DataFrame({'ligand':['TGFB1'],'receptor':['TGFBR1+TGFBR2']}).to_csv(lr,sep='\t',index=False)
    return obj,path,lr


@pytest.mark.parametrize('encoding',['dense','csr','csc'])
def test_legacy_readers_preserve_all_input_rows_features_and_values(tmp_path,monkeypatch,encoding):
    obj,path,_=inputs(tmp_path,encoding)
    params=api.FastCommParams(h5ad=path,gene_symbol_col='symbol')
    expected=api._matrix_inputs_from_h5ad(params)
    legacy_imports(monkeypatch)
    actual=api._matrix_inputs_from_h5ad(params)
    for a,b in [(expected[0],actual[0])]:
        np.testing.assert_array_equal(a.toarray() if sparse.issparse(a) else a,b.toarray() if sparse.issparse(b) else b)
    assert actual[1]==expected[1]==list(obj.obs_names)
    assert actual[2]==expected[2]==['TGFB1','TGFBR1','TGFBR2']
    pd.testing.assert_frame_equal(actual[3],obj.obs)
    assert actual[4]==expected[4]
    assert api._peek_input_genes(params)==['TGFB1','TGFBR1','TGFBR2']


def test_split_scoring_works_on_legacy_api_with_identical_native_results(tmp_path,monkeypatch):
    obj,path,lr=inputs(tmp_path,'csr')
    options=dict(h5ad=path,lr_table=lr,response_matrix=None,species='human',gene_symbol_col='symbol',
        state_key='state',split_key='Library',min_cells=2,include_self_edges=False)
    expected=run_benchmark(FastCommBenchmarkParams(**options,output_dir=tmp_path/'modern'))
    legacy_imports(monkeypatch)
    actual=run_benchmark(FastCommBenchmarkParams(**options,output_dir=tmp_path/'legacy'))
    for filename in ['full_scores.tsv','split_scores_long.tsv']:
        a=pd.read_csv(tmp_path/'modern'/filename,sep='\t');b=pd.read_csv(tmp_path/'legacy'/filename,sep='\t')
        pd.testing.assert_frame_equal(a,b,check_exact=True)
    assert set(pd.read_csv(tmp_path/'legacy'/'split_scores_long.tsv',sep='\t')['split'])=={'sample1','sample2'}


def test_missing_package_has_accurate_diagnostic(monkeypatch):
    original=builtins.__import__
    def missing(name,*args,**kwargs):
        if name=='anndata':raise ModuleNotFoundError("No module named anndata",name='anndata')
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',missing)
    with pytest.raises(ImportError,match="requires the 'anndata' package"):api._anndata_readers()
