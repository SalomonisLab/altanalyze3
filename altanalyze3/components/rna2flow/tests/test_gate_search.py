import numpy as np
import pytest
from altanalyze3.components.rna2flow.gate_search import optimize_population, membership
from altanalyze3.components.rna2flow.crosswalk import normalize_marker


def test_search_recovers_negative_marker_and_two_sided_interval():
    rng=np.random.default_rng(1);X=rng.uniform(-1,1,(5000,4))
    truth=(X[:,0]>.1)&(X[:,0]<.8)&(X[:,1]<-.2)
    r=optimize_population(X,np.where(truth,'state','other'),['A','B','noise1','noise2'],['state'])
    assert r['best']['test']['f1']>.85
    assert r['best']['shuffled_test_labels']['f1']<.3
    assert r['best']['test']['n_events']==1000
    assert r['best']['validation']['n_events']==1000
    assert set(c['channel'] for c in r['best']['conditions'])>={'A','B'}


def test_missing_population_is_not_silently_changed():
    with pytest.raises(ValueError,match='No events'):
        optimize_population(np.ones((100,2)),np.repeat('ML-1a',100),['A','B'],['ML-1b'])


def test_tied_threshold_membership_is_exact():
    X=np.array([[0.,1.],[1.,2.],[2.,3.]])
    assert membership(X,['A','B'],[dict(channel='A',operator='>',threshold=1.)]).tolist()==[False,False,True]
    with pytest.raises(ValueError):membership(X,['A','B'],[dict(channel='A',operator='bad',threshold=1.)])


@pytest.mark.parametrize('a,b', [('CD62L','CD62P-P-selectin-'),('CD11b','CD11c-Integrin-alpha-X'),('CD49d','CD49f-Integrin')])
def test_described_antibody_suffixes_stay_distinct(a,b):
    assert normalize_marker(a)!=normalize_marker(b)


def test_descriptions_do_not_become_cd_chains():
    assert normalize_marker('mouse_CD117_c_kit')=='cd117'
    assert normalize_marker('CD127-IL-7R')==normalize_marker('IL-7R')


def test_signed_gate_preserves_polarity_on_held_out_events():
    from altanalyze3.components.rna2flow.gate_search import optimize_signed_gate
    rng=np.random.default_rng(4);X=rng.normal(size=(3000,3));truth=(X[:,0]>.2)&(X[:,1]<-.1)
    r=optimize_signed_gate(X,np.where(truth,'CLP','other'),['Kit','CD11b','noise'],['CLP'],{'Kit':'>','CD11b':'<='})
    assert r['best']['test']['f1']>.9
    assert [(c['channel'],c['operator']) for c in r['best']['conditions']]==[('Kit','>'),('CD11b','<=')]
    assert r['best']['test']['n_events']==600
