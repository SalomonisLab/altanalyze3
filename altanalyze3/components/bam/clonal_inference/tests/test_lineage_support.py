import numpy as np
from altanalyze3.components.bam.clonal_inference.lineage_support import supported_lineage

def test_shared_lineage_unknown_is_not_a_descendant():
 G=np.array([[0,0,0],[1,0,0],[1,1,0],[1,0,1]])
 X=np.array([[1,-1,-1],[1,1,-1],[1,-1,1],[1,1,1],[-1,-1,-1],[0,0,-1]])
 r=supported_lineage(X,G,np.arange(4));assert r['node'].tolist()==[1,2,3,-1,0,0]
 assert r['incompatible'].tolist()==[False,False,False,True,False,False]

def test_negative_rna_is_flagged_not_forced_absence():
 G=np.array([[0,0],[1,0],[1,1]]);X=np.array([[0,1]])
 r=supported_lineage(X,G,np.arange(3));assert r['node'][0]==2 and r['negative_evidence_conflicts'][0]==1

def test_cell_time_order_does_not_change_lineage():
 G=np.array([[0,0],[1,0],[1,1]]);X=np.array([[1,-1],[1,1],[1,-1]])
 r=supported_lineage(X,G,np.arange(3));perm=np.array([2,0,1]);q=supported_lineage(X[perm],G,np.arange(3));assert np.array_equal(r['node'][perm],q['node'])

def test_absent_sampled_parent_does_not_invent_branch():
 G=np.array([[0,0,0],[1,0,0],[1,1,0],[1,0,1]]);X=np.array([[1,-1,-1]])
 r=supported_lineage(X,G,np.array([0,2,3]));assert r['node'][0]==-1 and r['ambiguous_minimum'][0]

def test_unequal_incomparable_branches_remain_unresolved():
 G=np.array([[0,0,0,0],[1,1,0,0],[1,0,1,1]])
 r=supported_lineage(np.array([[1,-1,-1,-1]]),G,np.arange(3))
 assert r['node'][0]==-1 and r['ambiguous_minimum'][0]

def test_no_active_malignant_nodes_retains_positive_cells_unresolved():
 r=supported_lineage(np.array([[1],[-1]]),np.array([[0]]),np.array([0]))
 assert r['node'].tolist()==[-1,0]

def test_invalid_node_rejected():
 import pytest
 with pytest.raises(ValueError):supported_lineage(np.array([[1]]),np.array([[0],[1]]),np.array([0,5]))
