"""scALABLE's cap applies to selected samples and both DE comparison modes."""
import numpy as np
import pandas as pd
import anndata as ad
import pytest
from altanalyze3.components.cellHarmony.flask.pipeline import _sample_differential_cells, _read_differential_h5ad
from altanalyze3.components.cellHarmony.cellHarmony_differential import compute_pseudobulk_per_population
import importlib
web = importlib.import_module("altanalyze3.components.cellHarmony.webapp.app")


def annotations():
    return pd.DataFrame([(state, sample) for state in ['A', 'B'] for sample in ['s1','s2','s3','s4','s5','s6','excluded'] for _ in range(620)], columns=['state','donor']).set_axis([f'c{i}' for i in range(8680)])


def test_random_cap_is_per_state_and_selected_sample():
    obs = annotations()
    samples = [f's{i}' for i in range(1,7)]
    rows = _sample_differential_cells(obs, 'state', 'donor', samples)
    assert len(rows) == 6000
    assert obs.iloc[rows].groupby(['state','donor']).size().eq(500).all()
    assert 'excluded' not in set(obs.iloc[rows].donor)
    np.testing.assert_array_equal(rows, _sample_differential_cells(obs, 'state', 'donor', samples))
    assert not np.array_equal(rows[:500], np.arange(500))
    small = obs.iloc[:12]
    np.testing.assert_array_equal(_sample_differential_cells(small,'state','donor',['s1']),np.arange(12))


@pytest.mark.parametrize('dense', [True,False])
def test_capped_backed_pseudobulk_matches_capped_memory(tmp_path, dense):
    import scipy.sparse as sp
    obs = annotations()
    values = np.random.default_rng(5).integers(1,20,(len(obs),8)).astype(np.float32)
    obj = ad.AnnData(values if dense else sp.csr_matrix(values), obs=obs)
    obj.layers['counts'] = obj.X.copy()
    obj.obs['condition'] = np.where(obj.obs.donor.isin(['s1','s2','s3']), 'case','control')
    path=tmp_path/'data.h5ad'
    obj.write_h5ad(path)
    rows=_sample_differential_cells(obs,'state','donor',[f's{i}' for i in range(1,7)])
    disk=_read_differential_h5ad(path,disk_backed=True)
    try:
        _, expected=compute_pseudobulk_per_population(obj[rows].copy(),'state','donor','condition',10,str(tmp_path/'ram'))
        _, actual=compute_pseudobulk_per_population(disk[rows],'state','donor','condition',10,str(tmp_path/'disk'))
        a,b=ad.read_h5ad(expected),ad.read_h5ad(actual)
        np.testing.assert_allclose(a.X,b.X,rtol=1e-6)
        assert a.obs.equals(b.obs)
        assert len(b)==12
    finally:
        disk._analysis_h5_handle.close()


def test_web_default_and_explicit_mode(monkeypatch):
    samples=[f's{i}' for i in range(1,7)]
    options=dict(enabled=True,modalities=[{'id':'rna'}],population_columns=[{'value':'state'}],sample_fields=[{'value':'donor'}],default_sample_field='donor',sample_values={'donor':samples},sample_names=samples)
    monkeypatch.setattr(web,'_differential_options',lambda meta:options)
    meta={'status':'completed','differential_options':{'upload_profile':{'single_h5ad':True}}}
    args=dict(population_col='state',sample_field='donor',group1_samples=samples[:3],group2_samples=samples[3:])
    config=web._validate_differential_request(meta,web.DifferentialSettings(**args))
    assert config['comparison_type']=='pseudobulk'
    assert config['max_cells_per_state_sample']==500
    assert web._validate_differential_request(meta,web.DifferentialSettings(**args,comparison_type='cells'))['comparison_type']=='cells'
    args['group2_samples']=samples[3:5]
    assert web._validate_differential_request(meta,web.DifferentialSettings(**args))['comparison_type']=='cells'
    with pytest.raises(ValueError):
        web.DifferentialSettings(**args,max_cells_per_state_sample=501)


def test_selected_sparse_expression_is_loaded_once_with_bounded_budget(tmp_path,monkeypatch):
    import scipy.sparse as sp
    from altanalyze3.components.cellHarmony.flask import pipeline
    from altanalyze3.components.cellHarmony import disk_differential as disk
    matrix=sp.csr_matrix(np.arange(240,dtype=np.float32).reshape(12,20))
    obj=ad.AnnData(matrix)
    obj.layers['counts']=matrix.copy()
    path=tmp_path/'sparse.h5ad'
    obj.write_h5ad(path)
    source=_read_differential_h5ad(path,disk_backed=True)
    reads=[]
    original=disk.read_rows
    def tracked(matrix,rows,columns=slice(None)):
        reads.append(rows)
        return original(matrix,rows,columns)
    monkeypatch.setattr(disk,'read_rows',tracked)
    try:
        chosen=np.array([1,4,8])
        got=pipeline._materialize_differential_selection(source[chosen],'rna')
        np.testing.assert_array_equal(got.X.toarray(),matrix[chosen].toarray())
        np.testing.assert_array_equal(got.layers['counts'].toarray(),matrix[chosen].toarray())
        assert disk.root_and_rows(got) is None
        assert len(reads)==2
        for rows in reads:
            np.testing.assert_array_equal(rows,chosen)
    finally:
        source._analysis_h5_handle.close()
