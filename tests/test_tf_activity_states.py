from types import SimpleNamespace
import numpy as np
from altanalyze3.components.cellHarmony import grn_analysis as G
from test_grn_web_release import uploaded


def dataset():
    store=SimpleNamespace(features=['SPI1','OTHER'],stats_mean=np.array([[1.,8.,np.nan,3.],[99.,9.,7.,6.]]),kind='per_cell',label='TF activity',statistic='mean stored TF activity')
    return SimpleNamespace(states=['HSC-1','Mono','missing','DC'],modality=lambda _:store)


def test_named_factor_ranks_states_without_other_factors_or_nonfinite_values():
    r=G.tf_activity_profile(dataset(),factors=['spi1'],contrast='unused')
    assert r['by_cell_state'] and not r['contrast']
    assert [v['cell_state'] for v in r['rows']]==['Mono','DC','HSC-1']
    assert [v['activity'] for v in r['rows']]==[8.,3.,1.]
    assert r['n_states']==3
    view=G.tf_activity_state_chat(r)
    assert 'Mono (8)' in view['answer'] and 'DC (3)' in view['answer']
    assert view['plot']['label_column']=='cell_state'
    assert view['table']['columns']==['cell_state','factor','activity']
    assert G.tf_activity_profile(dataset(),factors=['SPI1'],limit=1)['rows'][0]['cell_state']=='Mono'


def test_missing_factor_does_not_fall_back_to_other_factors():
    r=G.tf_activity_profile(dataset(),factors=['ABSENT'])
    assert not r['rows'] and r['missing_factors']==['ABSENT']
    assert G.tf_activity_state_chat(r)['status']=='not_found'


def test_with_state_preserves_factor_ranking():
    ds=dataset()
    r=G.tf_activity_profile(ds,cell_state='Mono')
    assert not r.get('by_cell_state')
    assert [v['factor'] for v in r['rows']]==['OTHER','SPI1']


def test_question_routes_where_without_limiting_default_to_25_states():
    for q in ['Where is SPI1 TF-activity most enriched.....','Which cell types have the highest SPI1 TF activity?']:
        r=G.read_regulatory_question(q,['HSC-1','Mono'],['SPI1'])
        assert r['intent']=='tf_activity' and r['genes']==['SPI1'] and not r['cell_state']
        assert r['limit']==G.MAX_LIMIT
    assert G.read_regulatory_question('Show top 2 cell types for SPI1 TF activity',['Mono'],['SPI1'])['limit']==2


def test_uploaded_chat_names_states(uploaded):
    u=uploaded
    r=u.client.post(f'/api/jobs/{u.job}/chat',json={'question':'Where is TF1 TF-activity most enriched?'})
    assert r.status_code==200,r.text
    d=r.json()
    assert d['status']=='ok' and d['by_cell_state']
    assert set(v['cell_state'] for v in d['table']['rows'])=={'AT1','AT2'}
    assert 'AT1' in d['answer'] and 'AT2' in d['answer']
