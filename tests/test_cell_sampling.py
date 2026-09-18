from collections import Counter
import numpy as np
import pytest
from altanalyze3.components.visualization.cell_sampling import sample_cell_indices


def test_sampling_is_balanced_nested_and_independent_of_input_order():
    ids = np.array([f'cell{i}' for i in range(240)])
    samples = np.repeat(['one', 'two'], 120)
    states = np.tile(np.repeat(['A', 'B'], 60), 2)
    previous = set()
    for limit in (5, 10, 20, 50):
        indices = sample_cell_indices(ids, samples, limit=limit, group_labels=states)
        assert set(Counter(zip(samples[indices], states[indices])).values()) == {limit}
        assert previous <= set(ids[indices])
        previous = set(ids[indices])
        shuffled = np.random.default_rng(7).permutation(len(ids))
        again = sample_cell_indices(ids[shuffled], samples[shuffled], limit=limit, group_labels=states[shuffled])
        assert set(ids[shuffled][again]) == set(ids[indices])
    assert len(sample_cell_indices(ids, samples, limit=10)) == 20
    assert len(sample_cell_indices(ids, samples, limit=0)) == len(ids)


def test_small_groups_are_retained_without_duplication():
    assert sample_cell_indices(['a','b'], ['sample','sample'], limit=50).tolist() == [0,1]
    assert sample_cell_indices([], [], limit=10).tolist() == []
    with pytest.raises(ValueError):
        sample_cell_indices(['a'], [], limit=10)
    with pytest.raises(ValueError):
        sample_cell_indices(['a'], ['sample'], limit=15)


def test_upload_sampling_matches_heatmap_and_pdf(uploaded, monkeypatch):
    from importlib import import_module
    import io
    import pandas as pd
    W = import_module('altanalyze3.components.cellHarmony.webapp.app')
    u = uploaded
    # The old plotted matrix has only one cell; the new display must use source cells.
    base = {'signature':{'source':'test'}, 'row_ids':np.array(['AT1:TF1','AT2:TF2']),
            'matrix':np.zeros((2,1)), 'col_ids':np.array(['AT1:cell0']), 'col_barcodes':np.array(['cell0'])}
    monkeypatch.setattr(W,'_get_marker_heatmap_cache_entry',lambda *a,**kw:base)
    url = f'/api/jobs/{u.job}'
    previous = set()
    previous_frame = None
    for limit in (5,10,20,50):
        response = u.client.get(url+'/combplot',params={'genes':'TF1','cells_per_sample':limit})
        assert response.status_code == 200,response.text
        data = response.json()
        field = data['sampling']['sample_field']
        selected = [c['cell'] for c in data['columns']]
        counts = Counter((str(u.a.obs.loc[c,field]),str(u.a.obs.loc[c,'cell_type'])) for c in selected)
        assert max(counts.values()) <= limit
        assert previous <= set(selected)
        previous = set(selected)
        response = u.client.get(url+'/marker/heatmap.tsv',params={'cells_per_sample':limit})
        assert response.status_code == 200,response.text
        frame = pd.read_csv(io.StringIO(response.text),sep='\t',index_col=0)
        assert {c.split(':',1)[1] for c in frame.columns} == set(selected)
        assert frame.index.tolist() == ['AT1:TF1','AT2:TF2']
        if previous_frame is not None:
            assert np.allclose(frame[previous_frame.columns], previous_frame)
        previous_frame = frame
    head = u.client.head(url+'/marker/heatmap.tsv',params={'cells_per_sample':10})
    frame = pd.read_csv(io.StringIO(u.client.get(url+'/marker/heatmap.tsv',params={'cells_per_sample':10}).text),sep='\t',index_col=0)
    assert int(head.headers['X-Marker-Columns']) == frame.shape[1]
    captured = []
    monkeypatch.setattr(W,'_render_marker_heatmap_pdf',lambda f,**kw:(captured.append(f) or io.BytesIO(b'%PDF-test')))
    assert u.client.get(url+'/marker/heatmap.pdf',params={'cells_per_sample':10}).status_code == 200
    assert captured[0].columns.tolist() == frame.columns.tolist()
    assert u.client.get(url+'/marker/heatmap.tsv',params={'cells_per_sample':15}).status_code == 400
    all_cells = u.client.get(url+'/combplot',params={'genes':'TF1','cells_per_sample':0}).json()
    assert len(all_cells['columns']) == u.a.n_obs
    assert u.client.get(url+'/combplot',params={'genes':'TF1'}).json()['sampling']['cells_per_sample'] == 10


from test_grn_web_release import uploaded
