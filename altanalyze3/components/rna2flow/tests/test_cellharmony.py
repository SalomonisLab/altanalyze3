import numpy as np
from altanalyze3.components.rna2flow.cellharmony import nearest_pearson,assign_communities,unit_rows


def test_nearest_pearson_matches_notebook_dense_argmax():
    rng=np.random.default_rng(4);C=rng.normal(size=(300,16));F=rng.normal(size=(150,16))
    score=unit_rows(F)@unit_rows(C).T;idx,rho=nearest_pearson(C,F,batch_size=17)
    assert np.array_equal(idx,score.argmax(1));assert np.allclose(rho,score.max(1))


def test_five_stage_assignment_matches_dense_notebook_with_frozen_communities():
    rng=np.random.default_rng(3);C=rng.normal(size=(80,9));F=rng.normal(size=(60,9));labels=np.array(['ML','CLP','other','MPP']*20);rg=np.arange(80)%5;qg=np.arange(60)%4
    actual,detail=assign_communities(C,F,labels,rg,qg)
    z=lambda x:(x-x.mean(0))/x.std(0);c,f=z(C),z(F)
    corr=unit_rows(np.vstack([f[qg==i].mean(0) for i in range(4)]))@unit_rows(np.vstack([c[rg==i].mean(0) for i in range(5)])).T
    maps={i:{int(corr[i].argmax())} for i in range(4)}
    for j in range(5):maps[int(corr[:,j].argmax())].add(j)
    initial=np.empty(60,object)
    for i in range(4):
        qi=np.flatnonzero(qg==i);ri=np.flatnonzero(np.isin(rg,list(maps[i])));scores=unit_rows(f[qi])@unit_rows(c[ri]).T;initial[qi]=labels[ri[scores.argmax(1)]]
    lev=np.unique(initial);cent=np.vstack([f[initial==l].mean(0) for l in lev]);expected=lev[(unit_rows(f)@unit_rows(cent).T).argmax(1)]
    assert np.array_equal(actual,expected)
