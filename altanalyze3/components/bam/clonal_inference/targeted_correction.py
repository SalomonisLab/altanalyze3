"""Explicitly authorized likelihood sensitivity test; native core is unchanged."""
import numpy as np
from .detection_tree import evidence

def correct_evidence(X,A,T,mu,sample,l0,l1,records,concentration=2.,fp=.01,fn=.30):
    b0,b1=evidence(X,A,T,mu,'counts',concentration,fp,fn)
    c0,c1=l0.copy(),l1.copy();changes=[]
    for r in records:
        r=dict(r);s,j=r['sample_index'],r['event_index']
        reversed_=r['TPR']<=r['FPR']
        empty=min(r['effective_mutant_validation_cells'],r['effective_WT_validation_cells'])<=1e-8
        mask=(sample==s)&(T[:,j]==0)
        if reversed_ or empty:
            c0[mask,j]=b0[mask,j];c1[mask,j]=b1[mask,j]
        r.update(fallback=bool(reversed_ or empty),reversed_calibration=bool(reversed_),empty_validation_class=bool(empty),fallback_entries=int(mask.sum()) if reversed_ or empty else 0)
        changes.append(r)
    unknown=(T==0)&(X==-1)
    c0[unknown]=0.;c1[unknown]=0.
    assert np.array_equal(c0[T>0],l0[T>0]) and np.array_equal(c1[T>0],l1[T>0])
    assert np.isfinite(c0).all() and np.isfinite(c1).all()
    return c0,c1,changes
