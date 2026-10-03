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


def test_parent_losses_are_counted_in_full_starting_population():
    rng=np.random.default_rng(81);X=rng.uniform(-1,1,(6000,3))
    labels=np.where((X[:,0]>.1)&(X[:,1]<.4),'state','other')
    parent=[dict(channel='CD34',operator='>',threshold=0.)]
    r=optimize_population(X,labels,['Kit','CD25','CD34'],['state'],parent_conditions=parent)
    assert r['best']['test']['n_events']==1200
    assert r['best']['test_within_parent']['recall']>.85
    assert r['best']['test']['recall']<.6
    assert .4<r['parent_summary']['recall']<.6
    assert r['parent_warning']
    assert r['best']['conditions'][0]['fixed_parent']
    mask=membership(X,['Kit','CD25','CD34'],r['best']['conditions'])
    assert not mask[X[:,2]<=0].any()


def test_missing_parent_marker_is_not_imputed_or_ignored():
    rng=np.random.default_rng(2);X=rng.normal(size=(500,2))
    with pytest.raises(ValueError,match='Parent marker unavailable'):
        optimize_population(X,np.where(X[:,0]>0,'state','other'),['Kit','CD25'],['state'],
                            parent_conditions=[dict(channel='CD34',operator='>',threshold=0.)])


def test_greedy_strategies_preserve_training_retention_and_marker_budget():
    from altanalyze3.components.rna2flow.gate_search import _greedy_candidates
    rng=np.random.default_rng(31);X=rng.normal(size=(3000,5));target=(X[:,:3]>.0).all(1)
    candidates=_greedy_candidates(X,target,['A','B','C','D','E'],1.,3)
    assert candidates
    last=target.sum()
    for c in candidates:
        mask=membership(X,['A','B','C','D','E'],c['conditions']);n=target[mask].sum()
        assert n>=.95*last and len(c['channels'])<=3
        last=n


def test_coordinate_refinement_recovers_rare_target_without_extra_markers():
    from altanalyze3.components.rna2flow.gate_search import _refine_candidates, metrics
    rng=np.random.default_rng(101);X=rng.uniform(0,1,(12000,3));target=(X[:,0]>.8)&(X[:,1]>.9)
    seed=dict(conditions=[dict(channel='A',operator='>',threshold=.9),dict(channel='B',operator='>',threshold=.85)])
    before=metrics(membership(X,['A','B','noise'],seed['conditions']),target)['f1']
    cs=_refine_candidates(X,target,['A','B','noise'],[seed],1.,2)
    best=max(metrics(membership(X,['A','B','noise'],c['conditions']),target)['f1'] for c in cs)
    assert best>.95 and best>before+.2
    assert all(len(c['channels'])<=2 for c in cs)


def test_serial_gating_counts_losses_and_final_intersection():
    from altanalyze3.components.rna2flow.gate_search import serial_audit
    X=np.array([[0.,0.],[1.,0.],[1.,1.],[2.,2.]])
    target=np.array([True,True,True,False]);cs=[dict(channel='A',operator='>',threshold=.5),dict(channel='B',operator='<=',threshold=.5)]
    rows=serial_audit(X,['A','B'],target,cs)
    assert rows[1]['target_lost']==1 and rows[2]['target_lost']==1
    assert rows[2]['recall']==pytest.approx(1/3)
    assert rows[2]['stage_recovery']==.5
    assert rows[-1]['n_gated']==membership(X,['A','B'],cs).sum()


def test_impossible_isolation_constraints_are_explicit():
    rng=np.random.default_rng(93);X=rng.normal(size=(3000,3));labels=np.where(rng.random(3000)<.1,'rare','other')
    r=optimize_population(X,labels,['A','B','C'],['rare'],min_purity=.99,min_recovery=.99)
    assert not r['constraints']['met']
    assert r['best']['test']['n_events']==600


def test_single_marker_isolation_budget():
    rng=np.random.default_rng(41);X=rng.uniform(0,1,(3000,3));labels=np.where(X[:,0]>.85,'rare','other')
    r=optimize_population(X,labels,['Kit','noise1','noise2'],['rare'],max_markers=1)
    assert r['best']['test']['f1']>.95
    assert all(len(set(c['channels']))<=1 for c in r['candidates'])
