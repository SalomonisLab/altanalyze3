"""Generic CITE-only bundle and serial virtual-gate contract (no project-specific data)."""
import json
import numpy as np
import pytest
from altanalyze3.components.visualization.flow_viewer.app import create_app


@pytest.fixture
def client(tmp_path):
    rng=np.random.default_rng(15);X=rng.uniform(0,1,(4000,4)).astype('float32')
    target=(X[:,0]>.8)&(X[:,1]<.4)
    # Deliberately perfect RNA annotation signal must not become a sort channel.
    X[:,2]=target;features=['mouse_CD117_c_kit','mouse_CD25','RNA:Spi1','DEAD']
    (tmp_path/'arrays').mkdir();X.tofile(tmp_path/'arrays/features.bin')
    target.astype('int16').tofile(tmp_path/'arrays/labels.bin')
    (X[:,3]<.6).astype('uint8').tofile(tmp_path/'arrays/live.bin')
    spec=dict(n=len(X),features=features,features_file='features.bin',embeddings={},
              labels=dict(state=dict(file='labels.bin',levels=['other','rare'])))
    manifest=dict(spaces=dict(cite_test=spec),gatesets=dict(capture=dict(space='cite_test',nodes=[dict(path='Live',mask_file='live.bin')])))
    (tmp_path/'manifest.json').write_text(json.dumps(manifest))
    return create_app(str(tmp_path)).test_client()


def test_cite_only_discovery_uses_adts_and_serial_counts_match(client):
    q=dict(space='cite_test',label_set='state',populations=['rare'],max_markers=2)
    r=client.post('/api/optimize',json=q);assert r.status_code==200,r.json
    result=r.json
    assert set(result['best']['channels'])=={'mouse_CD117_c_kit','mouse_CD25'}
    assert result['best']['test']['f1']>.9
    steps=[dict(conditions=[c]) for c in result['best']['conditions']]
    audit=client.post('/api/isolation/audit',json=dict(q,steps=steps)).json
    assert audit['stages'][-1]['n_gated']==result['best']['all_events']['n_gated']
    assert sum(s['target_lost'] for s in audit['stages'])==audit['stages'][0]['n_target']-audit['stages'][-1]['true_positive']


def test_optional_capture_and_manual_serial_rectangle(client):
    q=dict(space='cite_test',label_set='state',populations=['rare'],parent=dict(gateset='capture',path='Live'),
           steps=[dict(name='Manual rectangle',conditions=[dict(channel='mouse_CD117_c_kit',operator='>',threshold=.8),
                                                          dict(channel='mouse_CD25',operator='<=',threshold=.4)])])
    r=client.post('/api/isolation/audit',json=q);assert r.status_code==200,r.json
    stages=r.json['stages'];assert stages[1]['stage'].startswith('Capture:')
    assert stages[-1]['precision']==1.
    assert .5<stages[-1]['recall']<.7
    assert stages[-1]['stage_recovery']==1.
    assert stages[1]['target_lost']>0


def test_antibody_aliases_resolve_and_absent_capture_markers_fail(client):
    q=dict(space='cite_test',label_set='state',populations=['rare'],signs={'CD117':'>','CD25':'<='})
    r=client.post('/api/optimize',json=q);assert r.status_code==200,r.json
    assert {c['channel'] for c in r.json['best']['conditions']}=={'mouse_CD117_c_kit','mouse_CD25'}
    q.pop('signs');q['parent_conditions']=[dict(channel='CD34',operator='>',threshold=.4)]
    r=client.post('/api/optimize',json=q);assert r.status_code==400 and 'CD34' in r.json['error']
    templates=client.get('/api/capture_templates?space=cite_test').json['templates']
    assert templates[0]['conditions'][1]['channel']=='mouse_CD117_c_kit'
    assert 'CD34' in templates[0]['missing']


def test_generic_adt_bundle_preserves_cell_order_and_labels(tmp_path):
    import anndata as ad
    import pandas as pd
    from scipy import sparse
    from altanalyze3.components.visualization.flow_viewer.precompute_cite import export_cite
    from altanalyze3.components.visualization.flow_viewer.data_api import FlowBundle
    X=np.array([[1.,2.],[3.,4.],[5.,6.]],np.float32)
    data=ad.AnnData(sparse.csr_matrix(X),obs=pd.DataFrame({'state':['CLP','HSC','CLP']},index=['c3','c1','c2']),var=pd.DataFrame(index=['CD117','CD25']))
    data.obsm['X_umap']=np.array([[0.,1.],[1.,1.],[1.,0.]])
    export_cite(data,tmp_path,name='experiment',labels=['state'],scale='DSB')
    bundle=FlowBundle(str(tmp_path));assert bundle.default_space=='experiment'
    np.testing.assert_array_equal(bundle.feature_matrix(),X)
    assert (tmp_path/'arrays/cite_experiment_cell_ids.txt').read_text().splitlines()==['c3','c1','c2']
    decoded=[bundle.levels('state')[v] for v in bundle.labels('state')]
    assert decoded==['CLP','HSC','CLP']
    with pytest.raises(ValueError,match='already exists'):export_cite(data,tmp_path,name='experiment',labels=['state'],scale='DSB')
