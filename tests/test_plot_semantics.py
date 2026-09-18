"""Missing values, signed measurements and display filters must stay truthful."""
from importlib import import_module
from types import SimpleNamespace
import numpy as np
import pandas as pd
from fastapi import FastAPI
from fastapi.testclient import TestClient
from test_grn_web_release import uploaded

W = import_module('altanalyze3.components.cellHarmony.webapp.app')
V = import_module('altanalyze3.components.visualization.scalable_viewer.scalable_app')


def test_uploaded_dotplot_uses_both_filters_and_empty_stays_empty(uploaded):
    u = uploaded
    base = f'/api/jobs/{u.job}/dotplot'
    response = u.client.get(base, params={'genes':'TF1','subset_by':'condition','subset_values':'case',
                                          'subset2_by':'cell_type','subset2_values':'AT2'})
    assert response.status_code == 200, response.text
    p = response.json()
    assert p['state_n'] == [0, 48]
    mask = (u.a.obs.condition == 'case') & (u.a.obs.cell_type == 'AT2')
    assert np.isclose(p['mean'][0][1], u.a.X[mask, 0].mean())
    for params in ({'subset_by':'condition','subset_values':'missing'}, {'groups':'missing'}):
        p = u.client.get(base, params={'genes':'TF1', **params}).json()
        assert sum(p['state_n']) == 0
    r = u.client.get(f'/api/jobs/{u.job}/combplot', params={'genes':'TF1','unit':'donor',
                     'subset_by':'condition','subset_values':'missing'})
    assert r.status_code == 404


def test_expression_never_substitutes_another_gene(uploaded, monkeypatch):
    u = uploaded
    base = f'/api/jobs/{u.job}/expression'
    r = u.client.get(base, params={'gene':'NOT_A_GENE','modality':'grn_tf'})
    assert r.status_code == 404
    monkeypatch.setattr(W, '_load_reference_states_table', lambda meta: pd.DataFrame())
    r = u.client.get(base, params={'gene':'NOT_A_GENE'})
    assert r.status_code == 200
    assert r.json()['gene'] == 'NOT_A_GENE'
    assert r.json()['source'] == 'missing' and not r.json()['umap']


def test_expression_filter_applies_to_every_payload_and_empty_pdf(uploaded):
    u = uploaded
    base = f'/api/jobs/{u.job}/expression'
    params = {'gene':'TF1','filter1_field':'condition','filter1_values':'case',
              'filter2_field':'cell_type','filter2_values':'AT2'}
    p = u.client.get(base, params=params).json()
    assert len(p['umap']) == len(p['scatter']) == 48
    assert [v['population'] for v in p['violin']] == ['AT2']
    params['filter2_values'] = 'missing'
    p = u.client.get(base, params=params).json()
    assert not p['umap'] and not p['scatter'] and not p['violin']
    r = u.client.get(base+'/pdf', params={**params,'mode':'violin'})
    assert r.status_code == 200 and r.content.startswith(b'%PDF')


def test_network_missing_folds_stay_missing(tmp_path, monkeypatch):
    markers = tmp_path/'markers.tsv'
    pd.DataFrame({'Gene':['UP','DOWN'], 'Fold':[2.,-2.], 'cluster':['A','A']}).to_csv(markers,sep='\t',index=False)
    edges = tmp_path/'edges.tsv'
    pd.DataFrame({'Symbol1':['UP','UP'], 'Symbol2':['UNKNOWN','DOWN'], 'Direction':['neutral','neutral']}).to_csv(edges,sep='\t',index=False)
    entry = {'id':'A', 'population':'A', 'tsv':str(edges)}
    meta = {'marker_analysis':{'markers_tsv':str(markers),'networks':[entry]}}
    p = W._build_marker_network_payload(meta, 'A')
    values = {e['data']['id']:e['data']['log2fc'] for e in p['elements'] if 'source' not in e['data']}
    assert values == {'UP':2.,'UNKNOWN':None,'DOWN':-2.}
    monkeypatch.setattr(W, '_differential_network_entry',lambda *a:entry)
    monkeypatch.setattr(W, '_get_differential_detail_table',lambda *a:pd.DataFrame([
        {'population':'A','gene':'UP','log2fc':2.,'fdr':.01,'pval':.001},
        {'population':'A','gene':'DOWN','log2fc':-2.,'fdr':.01,'pval':.001}]))
    monkeypatch.setattr(W, '_get_differential_cache_entry',lambda *a:{'network_tables':{}})
    p = W._build_differential_network_payload(None,{'job_id':'demo'},'A')
    values = {e['data']['id']:e['data']['log2fc'] for e in p['elements'] if 'source' not in e['data']}
    assert values == {'UP':2.,'UNKNOWN':None,'DOWN':-2.}
    assert W._render_network_pdf(p,'network').read(4) == b'%PDF'


def test_viewer_dotplot_positive_fraction_consistent_with_and_without_filters():
    ds = SimpleNamespace(sv={'cluster_key':'state'},states=['A','B'],state_code=np.array([0,0,1,1]),
        state_n=[2,2], colors={'A':'red','B':'blue'}, names_have_pipe=False, markers=lambda:[],
        resolve_gene=lambda g:0 if g=='G' else None, stats_mean=np.array([[1.,0.]]),
        stats_frac=np.array([[1.,1.]]), gene_column=lambda r:(np.arange(4),np.array([-2.,4.,0.,0.])),
        covariate_values=lambda field:('categorical', np.zeros(4,dtype=int),['all']))
    app = FastAPI(); V._install_dotplot_routes(app,SimpleNamespace(dataset=lambda _:ds),{})
    client = TestClient(app); url = '/api/jobs/demo/dotplot'
    for params in ({}, {'subset_by':'condition','subset_values':'all'}):
        r = client.get(url,params={'genes':'G',**params})
        assert r.status_code == 200,r.text
        assert r.json()['frac'] == [[.5,0.]]
    for params in ({'subset_by':'state','subset_values':'missing'},
                   {'subset_by':'condition','subset_values':'missing'}, {'groups':'missing'}):
        r = client.get(url,params={'genes':'G',**params})
        assert r.status_code == 200,r.text
        assert sum(r.json()['state_n']) == 0


def test_volcano_never_labels_raw_p_as_fdr(monkeypatch):
    frame = pd.DataFrame([
        {'population':'A','gene':'FDR','log2fc':1.,'fdr':.05,'pval':.001},
        {'population':'A','gene':'RAW','log2fc':-2.,'fdr':np.nan,'pval':.0001},
        {'population':'A','gene':'MISSING','log2fc':3.,'fdr':np.nan,'pval':np.nan},
        {'population':'A','gene':'INVALID','log2fc':4.,'fdr':-1.,'pval':-1.},
        {'population':'A','gene':'ZERO','log2fc':1.,'fdr':0.,'pval':0.},
    ])
    monkeypatch.setattr(W, '_get_differential_detail_table',lambda *args:frame)
    p = W._build_differential_volcano_payload(None,{},'A')
    assert p['statistic_label'] == 'FDR'
    assert {r['gene'] for r in p['points']} == {'FDR','ZERO'}
    assert p['n_missing_statistic'] == 3
    assert next(r['score'] for r in p['points'] if r['gene']=='ZERO') == 300
    frame['fdr'] = np.nan
    p = W._build_differential_volcano_payload(None,{},'A')
    assert p['statistic_label'] == 'p-value'
    assert {r['gene'] for r in p['points']} == {'FDR','RAW','ZERO'}
    assert next(r['score'] for r in p['points'] if r['gene']=='RAW') == 4
