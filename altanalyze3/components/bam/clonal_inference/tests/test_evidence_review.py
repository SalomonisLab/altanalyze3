import numpy as np
import pandas as pd
import pytest
from altanalyze3.components.bam.clonal_inference.evidence_review import evidence_consistency_audit


def fixture():
    X=np.array([[1,-1],[1,1],[1,0],[-1,-1]])
    cells=pd.DataFrame({'cell_id':['a','b','c','d'],'sample':['pre','post','post','post']})
    events=pd.DataFrame({'feature':['f','branch'],'event_label':['Founder','Branch']})
    model=dict(G=np.array([[0,0],[1,0],[1,1]]),active=np.arange(3),
        cell_ids=cells.cell_id.to_numpy(),features=events.feature.to_numpy(),
        log_evidence=np.array([[-5,0,0],[-5,-3,0],[-5,0,-1],[0,0,0]]))
    assignments=cells.assign(node=[2,1,1,0])
    return X,cells,events,model,assignments


def test_unknown_branch_and_omitted_positive_are_distinct():
    args=fixture();before=args[0].copy();assign_before=args[4].copy()
    e,p,s=evidence_consistency_audit(*args)
    pre=p[p['sample']=='pre'].iloc[0]
    assert pre.distinguishing_loci_all_unknown==1 and pre.equal_data_likelihood_cells==1
    post=e[(e['sample']=='post')&(e.node==1)&(e.event=='Branch')].iloc[0]
    assert post.positive_calls_omitted_from_genotype==1
    assert post.accepted_reference==1 and post.unknown==0
    assert s['omitted_positive_cell_event_calls']==1
    assert s['cells']==4 and s['events']==2 and s['excluded_cells']==s['excluded_events']==0
    np.testing.assert_array_equal(args[0],before);pd.testing.assert_frame_equal(args[4],assign_before)


def test_identity_swap_rejected():
    X,c,e,m,a=fixture();m['cell_ids']=m['cell_ids'][::-1]
    with pytest.raises(Exception,match='identity/order'):
        evidence_consistency_audit(X,c,e,m,a)


def test_feature_permutation_preserves_diagnostics():
    X,c,e,m,a=fixture();r,_,s=evidence_consistency_audit(X,c,e,m,a)
    m=dict(m,G=m['G'][:,::-1],features=m['features'][::-1])
    t,_,u=evidence_consistency_audit(X[:,::-1],c,e.iloc[::-1].reset_index(drop=True),m,a)
    cols=['sample','node','feature']
    pd.testing.assert_frame_equal(r.sort_values(cols).reset_index(drop=True),t.sort_values(cols).reset_index(drop=True))
    assert s==u


def test_unknowns_do_not_automatically_reject_or_reassign_clone():
    X,c,e,m,a=fixture();X[1,1]=-1
    _,p,s=evidence_consistency_audit(X,c,e,m,a)
    assert s['omitted_positive_cell_event_calls']==0 and not s['changed_assignments']
    assert len(p)==2
