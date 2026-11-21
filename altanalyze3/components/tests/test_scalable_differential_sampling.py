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
    frame=pd.DataFrame([(state, sample) for state in ['A', 'B'] for sample in ['s1','s2','s3','s4','s5','s6','excluded'] for _ in range(620)], columns=['state','donor']).set_axis([f'c{i}' for i in range(8680)])
    frame['AlignmentScore']=np.tile(np.arange(620,dtype=float),14)
    return frame


def test_highest_score_cap_is_per_state_and_selected_sample():
    obs = annotations()
    samples = [f's{i}' for i in range(1,7)]
    rows = _sample_differential_cells(obs, 'state', 'donor', samples)
    assert len(rows) == 6000
    assert obs.iloc[rows].groupby(['state','donor']).size().eq(500).all()
    assert 'excluded' not in set(obs.iloc[rows].donor)
    np.testing.assert_array_equal(rows, _sample_differential_cells(obs, 'state', 'donor', samples))
    assert obs.iloc[rows].AlignmentScore.min()==120
    shuffled=obs.sample(frac=1,random_state=9)
    assert obs.index[rows].equals(shuffled.index[_sample_differential_cells(shuffled,'state','donor',samples)])
    small = obs.iloc[:12]
    assert small.index[_sample_differential_cells(small,'state','donor',['s1'])].tolist()==sorted(small.index)


@pytest.mark.parametrize('dense', [True,False])
@pytest.mark.parametrize('score_column', ['AlignmentScore', 'SVM_score'])
def test_capped_backed_pseudobulk_matches_capped_memory(tmp_path, dense, score_column):
    import scipy.sparse as sp
    obs = annotations()
    obs=obs.rename(columns={'AlignmentScore':score_column})
    values = np.random.default_rng(5).integers(1,20,(len(obs),8)).astype(np.float32)
    obj = ad.AnnData(values if dense else sp.csr_matrix(values), obs=obs)
    obj.layers['counts'] = obj.X.copy()
    obj.obs['condition'] = np.where(obj.obs.donor.isin(['s1','s2','s3']), 'case','control')
    path=tmp_path/'data.h5ad'
    obj.write_h5ad(path)
    rows=_sample_differential_cells(obs,'state','donor',[f's{i}' for i in range(1,7)],scores=obs[score_column])
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


@pytest.mark.parametrize('score_column', ['AlignmentScore', 'SVM_score'])
def test_native_deg_statistics_repeat_after_row_permutation(tmp_path, score_column):
    from altanalyze3.components.cellHarmony import cellHarmony_differential as native
    obs=pd.DataFrame([(state,sample,'case' if sample=='s1' else 'control',float(i))
        for state in ['A','B'] for sample in ['s1','s2'] for i in range(80)],
        columns=['state','donor','condition',score_column],index=[f'cell{i:04d}' for i in range(320)])
    x=np.random.default_rng(137).uniform(.2,2.,(320,30)).astype(np.float32)
    x[obs.condition.eq('case').to_numpy(),:5]+=2
    obj=ad.AnnData(x,obs=obs,var=pd.DataFrame(index=[f'g{i}' for i in range(30)]))
    obj.uns['log1p']={'base':None}
    if score_column=='SVM_score':
        from altanalyze3.components.cellHarmony.differential_sampling import web_assignment_scores
        path=tmp_path/'native_icgs3.h5ad'
        obj.write_h5ad(path)
        scores,_,_,_=web_assignment_scores({'analysis_mode':'unsupervised','artifacts':{'combined_h5ad':str(path)}},obj.obs_names)
    else:
        scores=obs[score_column]
    permuted=obj[np.random.default_rng(53).permutation(len(obj))].copy()
    outputs=[]
    for data in (obj,permuted):
        rows=_sample_differential_cells(data.obs,'state','donor',['s1','s2'],limit=30,scores=scores)
        selected=data[rows].copy()
        ranks=[]
        for state in ['A','B']:
            names,fdr,fc,p=native._rank_genes_scanpy(selected[selected.obs.state.eq(state)].copy(),
                                                   'condition','case','control','wilcoxon')
            ranks.append(pd.DataFrame({'fdr':fdr.values,'fc':fc.values,'p':p.values},index=names).sort_index())
        outputs.append((selected.obs_names,ranks))
    assert outputs[0][0].equals(outputs[1][0])
    for first,second in zip(outputs[0][1],outputs[1][1]):
        pd.testing.assert_frame_equal(first,second,check_exact=True)
        assert (first.fdr<.05).any()  # A real called-gene result, not an empty-output equality.
