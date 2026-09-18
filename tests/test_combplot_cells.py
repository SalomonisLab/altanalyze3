"""Individual observations are the CombPlot default in both applications."""
from importlib import import_module
from types import SimpleNamespace
import anndata as ad
import numpy as np
import pandas as pd
from fastapi import FastAPI
from fastapi.testclient import TestClient
from scipy import sparse
W = import_module('altanalyze3.components.cellHarmony.webapp.app')
V = import_module('altanalyze3.components.visualization.scalable_viewer.scalable_app')


def test_uploaded_cells_and_donor_means():
    obs = pd.DataFrame({'state': pd.Categorical(['B','A','A','B'], categories=['A','B']),
                        'Library': ['d1','d1','d1','d2'], 'condition': ['case','case','control','control']},
                       index=['c0','c1','c2','c3'])
    a = ad.AnnData(sparse.csr_matrix([[0.],[1.],[9.],[3.]]), obs=obs, var=pd.DataFrame(index=['G1']))
    cache = {'adata':a, 'cluster_key':'state', 'sample_field':'Library', 'var_names':a.var_names.to_numpy(),
             'populations':a.obs.state.astype(str).to_numpy()}
    p = W._gene_cell_values(cache, ['G1'])
    assert p['unit'] == 'cells' and p['values'] == [[1.,9.,0.,3.]]
    assert [c['cell'] for c in p['columns']] == ['c1','c2','c0','c3']
    p = W._gene_cell_values(cache, ['G1'], subset_by='Library', subset_values=['d1'],
                            subset2_by='condition', subset2_values=['control'])
    assert p['values'] == [[9.]]
    assert W._gene_cell_values(cache, ['G1'], subset_by='Library', subset_values=['missing'])['columns'] == []
    p = W._gene_donor_state_means(cache, ['G1'], 1)
    assert p['unit'] == 'donor' and p['values'] == [[5.,0.,3.]]


def test_viewer_cells_filters_tracks_and_donor_option(tmp_path):
    clusters = tmp_path/'clusters.tsv'
    clusters.write_text('barcode\tcluster\nc0\tB\nc1\tA\nc2\tA\nc3\tB\n')
    cats = {'Library':(np.array([0,0,0,1]), ['d1','d2']), 'condition':(np.array([0,0,1,1]), ['case','control'])}
    ds = SimpleNamespace(n_cells=4, sv={'cluster_key':'state'}, states=['A','B'], state_code=np.array([1,0,0,1]),
                         colors={'A':'#f00','B':'#00f'}, paths=SimpleNamespace(clusters=str(clusters)), names_have_pipe=False,
                         resolve_gene=lambda gene: 0 if gene=='G1' else None,
                         gene_column=lambda row:(np.array([1,2,3]),np.array([1.,9.,3.])),
                         covariate_names=lambda:{}, covariate_values=lambda key:('categorical', *cats[key]))
    app=FastAPI(); V._install_combplot_routes(app, SimpleNamespace(dataset=lambda job:ds), {})
    client=TestClient(app); url='/api/jobs/demo/combplot'
    r=client.get(url, params={'genes':'G1'})
    assert r.status_code==200, r.text
    assert r.json()['unit']=='cells' and r.json()['values']==[[1.,9.,0.,3.]]
    assert [c['cell'] for c in r.json()['columns']]==['c1','c2','c0','c3']
    ds.covariate_names=lambda:{k:{'kind':'categorical'} for k in cats}
    r=client.get(url,params={'genes':'G1','subset_by':'Library','subset_values':'d1',
                           'subset2_by':'condition','subset2_values':'control','tracks':'condition'})
    assert r.status_code==200,r.text
    assert r.json()['values']==[[9.]] and r.json()['tracks']=={'condition':['control']}
    r=client.get(url,params={'genes':'G1','unit':'donor','donor_key':'Library','min_cells':1})
    assert r.status_code==200,r.text
    assert r.json()['unit']=='donor' and r.json()['values']==[[5.,0.,3.]]
    ds.sv['observation_unit']='metacells'
    assert client.get(url,params={'genes':'G1'}).json()['observation_unit']=='metacells'
    assert client.get(url,params={'genes':'G1','unit':'invalid'}).status_code==422
