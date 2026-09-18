"""Chat parity, sample-level statistics, conservative routing and cache invalidation."""
import json
import numpy as np
import pandas as pd
import anndata as ad
import pytest
from scipy.stats import spearmanr
from test_grn_web_release import uploaded,seed_comparisons,web
from altanalyze3.components.cellHarmony.webapp.chat_service import read_question
from altanalyze3.components.cellHarmony.webapp.chat_data import data_for_job
from altanalyze3.components.cellHarmony import chat_protocols


def prepare(u):
    seed_comparisons(u)
    a=ad.read_h5ad(u.paths['rna']);donor=a.obs.Library.str.replace('donor','').astype(int)
    a.obs['age']=40+donor*4
    a.obs['stage']=pd.Categorical(np.where(donor<3,'early',np.where(donor<6,'middle','late')),categories=['early','middle','late'],ordered=True)
    a.obs['annotation']=a.obs.cell_type.astype(str).replace({'AT1':'alveolar 1','AT2':'alveolar 2'})
    a.write_h5ad(u.paths['rna']);u.app.state.expression_cache.clear()
    return f'/api/jobs/{u.job}'


def test_local_router_does_not_select_bookkeeping_or_ambiguous_request():
    r=read_question('Which genes in AT2 track FEV1?', ['AT1','AT2'],['SFTPC'],['fev1_percent_predicted__n_obs','wmean_fev1_percent_predicted','fev1_percent_predicted'])
    assert r['covariate']=='fev1_percent_predicted'
    assert read_question('What treatment should I use?', ['AT2'],[],[]) is None


def test_bundle_without_real_samples_retains_precomputed_queries(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from altanalyze3.components.cellHarmony.webapp.chat_data import protocol_adapter
    from altanalyze3.components.visualization.scalable_viewer import bundle_meta
    original=SimpleNamespace(sv={'integrated_root':str(tmp_path)},cluster_key='cell_state',symbols=['TF1'])
    obs=pd.DataFrame({'cell_state':['AT1'],'meta_sample':['imputed1']})
    monkeypatch.setattr(bundle_meta,'build_obs',lambda ds:(obs,None))
    ds=protocol_adapter(original)
    assert ds.states==['AT1'] and ds.symbols==['TF1']
    assert ds.donor_axis()==(None,None)
    assert protocol_adapter(original) is ds
    with pytest.raises(ValueError,match='Real sample pseudobulks'):
        ds.donor_pseudobulk([0],'AT1')


def test_uploaded_protocols_and_cached_queries(uploaded,monkeypatch):
    u=uploaded;base=prepare(u)
    monkeypatch.setattr(web.urllib.request,'urlopen',lambda *a,**k:pytest.fail('These queries must be local'))
    examples=[('Which genes in AT1 track age?','severity_gradient','gradient'),
              ('Which genes coexpress with TF1 in AT1?','coexpression_module','gradient'),
              ('Which cell states differ in abundance by condition?','composition_shift','frequency'),
              ('Compare annotation agreement in AT1','annotation_concordance','heatmap'),
              ('Show TF1 across stage levels in AT1','dose_response','gradient'),
              ('Which cell type is most affected?','most_affected_state','frequency')]
    for question,intent,kind in examples:
        r=u.client.post(base+'/chat',json={'question':question});assert r.status_code==200,r.text
        d=r.json();assert d['intent']==intent,d
        assert d['status']=='answered',d
        if kind:assert d['plot']['kind']==kind,d
        if intent=='most_affected_state':assert d['plot']['n_donors']==[4,4]
        else:assert not d['reading']['contrast']
        repeat=u.client.post(base+'/chat',json={'question':question}).json()
        assert repeat['performance']['cache_hit']
        assert repeat['table']==d['table']
    ds=data_for_job(u.app,u.store.get_job(u.job));matrix,donors,n=ds.donor_pseudobulk([0,1,2,3,4],'AT1')
    assert len(donors)==8 and set(n)=={12}
    r=chat_protocols._run_severity_gradient(ds,'AT1','age')
    ages=np.array([40+int(d.removeprefix('donor'))*4 for d in donors])
    for row in r['table']['rows']:
        i=ds.resolve_gene(row[0]);expected=spearmanr(ages,matrix[i])
        assert row[1]==pytest.approx(expected.statistic,abs=.0001)
    # A new generation cannot reuse the prior data/statistics.
    u.store.update_job(u.job,default_gene='TF2')
    assert data_for_job(u.app,u.store.get_job(u.job)) is not ds


def test_raw_counts_and_missing_fields_are_explicit(uploaded):
    u=uploaded;base=prepare(u)
    assert u.client.post(base+'/chat',json={'question':'Which genes in AT1 track age?'}).json()['status']=='answered'
    r=u.client.post(base+'/chat',json={'question':'Which genes track age?'}).json()
    assert r['status']=='clarify' and 'cell state' in r['answer']
    a=ad.read_h5ad(u.paths['rna']);del a.layers['counts'];a.write_h5ad(u.paths['rna'])
    r=u.client.post(base+'/chat',json={'question':'Which genes in AT1 track age?'}).json()
    assert r['status']=='not_covered' and 'Raw RNA counts' in r['answer']


def test_constant_and_missing_metadata_do_not_invent_correlations(uploaded):
    u=uploaded;prepare(u);ds=data_for_job(u.app,u.store.get_job(u.job))
    a=ds.rna['adata'];a.layers['counts'][:]=1.
    result=chat_protocols._run_coexpression(ds,'AT1','TF1')
    assert result['status']=='not_covered' and 'undefined' in result['answer']
    ds.obs.loc[ds.obs.Library=='donor0','condition']=None
    ds.__dict__.pop('_covariates',None)
    result=chat_protocols._run_composition(ds,'condition')
    assert result['plot']['n_donors']==[4,3] or result['plot']['n_donors']==[3,4]


def test_stage_order_missing_stage_and_all_negative_signature(uploaded):
    u=uploaded;prepare(u);ds=data_for_job(u.app,u.store.get_job(u.job))
    ds.obs['gold']=np.where(ds.obs.Library=='donor0','not applicable',np.where(ds.obs.Library.isin(['donor1','donor2']),'GOLD I','GOLD IV'))
    ds.__dict__.pop('_covariates',None)
    _,codes,levels=ds.ordered_covariate_values('gold')
    assert levels==['GOLD I','GOLD IV']
    assert np.all(codes[ds.obs.Library=='donor0']==-1)
    ds.obs['unordered']=['high','low']*(len(ds.obs)//2);ds.__dict__.pop('_covariates',None)
    with pytest.raises(ValueError,match='explicit ordered'):ds.ordered_covariate_values('unordered')
    original=ds.deg_table
    ds.deg_table=lambda *a,**k:{'rows':[dict(gene=g,log2fc=-1.,fdr=.01,population='AT1') for g in ds.symbols]}
    result=chat_protocols._run_donor_heterogeneity(ds,'AT1','run-rna')
    assert result['status']=='answered'
    assert all(np.isfinite(row[1]) for row in result['table']['rows'])
    assert len(result['plot']['donors'])==8
