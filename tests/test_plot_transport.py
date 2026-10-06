"""Compact plots preserve every row, value, identity and filter result."""
import importlib
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from fastapi.testclient import TestClient

from altanalyze3.components.cellHarmony.webapp.plot_payload import compact_plot_payload

W = importlib.import_module('altanalyze3.components.cellHarmony.webapp.app')


def unpack(payload):
    payload = dict(payload)
    for name in ('query', 'reference', 'umap', 'scatter'):
        table = payload.get(name)
        if isinstance(table, dict) and table.get('encoding') == 'columns-v1':
            payload[name] = [{key: table['dictionaries'][key][values[i]]
                             if key in table['dictionaries'] else values[i]
                             for key, values in table['columns'].items()}
                            for i in range(table['length'])]
    return payload


def test_transport_exact_round_trip():
    rows = [dict(barcode=f'cell-{i}α', population='state'*10 + str(i % 3),
                 sample=str(i % 2), x=i/7, y=-i/3, value=i/13) for i in range(1000)]
    payload = dict(query=rows, reference=[], message=None, n_cells_selected=len(rows))
    encoded = compact_plot_payload(payload)
    assert unpack(json.loads(json.dumps(encoded))) == payload
    assert len(json.dumps(encoded)) < len(json.dumps(payload)) * .6
    assert payload['query'] is rows


@pytest.mark.parametrize("viewer_wrapper", [False, True])
def test_view_specific_and_compact_api_match_legacy(tmp_path, monkeypatch, viewer_wrapper):
    app = W.create_app({'JOB_STORAGE': str(tmp_path)})
    store = app.state.job_store
    meta = store.create_job('human', 'ref', None, files=[])
    job = meta['job_id']
    obs = pd.DataFrame({'state':['A','B','A','B'], 'Library':['S1','S1','S2','S2']},
                       index=['c1','c2','c3','c4'])
    a = ad.AnnData(X=np.array([[0],[3],[2],[1]], dtype=np.float32), obs=obs,
                   var=pd.DataFrame(index=['G']))
    cache = dict(adata=a, obs_names=a.obs_names.to_numpy(), var_names=a.var_names.to_numpy(),
                 populations=obs.state.to_numpy(), umap_x=np.arange(4.), umap_y=np.arange(4.),
                 obs_filter_values={'Library':obs.Library.to_numpy()})
    monkeypatch.setattr(W, '_get_expression_cache', lambda *args, **kwargs: cache)
    if viewer_wrapper:
        V = importlib.import_module('altanalyze3.components.visualization.scalable_viewer.scalable_app')
        monkeypatch.setattr(V, '_to_feature_key', lambda app, job, modality, gene: gene)
        monkeypatch.setattr(V, '_retitle_with_display_name', lambda app, job, modality, response: response)
        V._install_expression_covariate_routes(app)
    client = TestClient(app)
    try:
        for filters in ({}, {'filter1_field':'Library', 'filter1_values':'S2'}):
            params = dict(gene='G', **filters)
            legacy = client.get(f'/api/jobs/{job}/expression', params=params).json()
            for view in ('all','umap','violin'):
                r = client.get(f'/api/jobs/{job}/expression', params=dict(params,view=view,compact='true'))
                assert r.status_code == 200, r.text
                decoded = unpack(r.json())
                assert decoded['umap'] == (legacy['umap'] if view != 'violin' else [])
                assert decoded['violin'] == (legacy['violin'] if view != 'umap' else [])
                assert decoded['scatter'] == (legacy['scatter'] if view == 'all' else [])
                assert decoded['global_min'] == legacy['global_min']
                assert decoded['global_max'] == legacy['global_max']
    finally:
        app.state.job_runner.executor.shutdown(wait=True)
