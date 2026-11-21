import numpy as np
from altanalyze3.components.bam.clonal_inference.detection_tree import evidence
from altanalyze3.components.bam.clonal_inference.targeted_correction import correct_evidence

def test_targeted_behavior():
 X=np.array([[-1,1,0],[1,0,-1],[0,1,1]],dtype=int);A=np.zeros_like(X);T=np.zeros_like(X);T[2]=3;A[2]=[0,1,2]
 mu=np.full(3,.5);s=np.zeros(3,int);b0,b1=evidence(X,A,T,mu)
 l0=b0.copy();l1=b1.copy();l0[:2]=-.8;l1[:2]=-.2
 records=[dict(sample_index=0,event_index=j,TPR=t,FPR=f,effective_mutant_validation_cells=n,effective_WT_validation_cells=5) for j,t,f,n in [(0,.1,.2,5),(1,.7,.1,0),(2,.7,.1,5)]]
 c0,c1,r=correct_evidence(X,A,T,mu,s,l0,l1,records)
 assert (c0[(T==0)&(X==-1)]==0).all() and (c1[(T==0)&(X==-1)]==0).all()
 assert np.array_equal(c0[T>0],l0[T>0]) and np.array_equal(c1[T>0],l1[T>0])
 assert np.array_equal(c0[:2,:2],b0[:2,:2]) and np.array_equal(c1[:2,:2],b1[:2,:2])
 assert c0[0,2]==l0[0,2] and c1[0,2]==l1[0,2]
 assert [v['fallback'] for v in r]==[True,True,False]
 assert X.tolist()==[[-1,1,0],[1,0,-1],[0,1,1]]


def test_missing_distinguishing_event_has_no_likelihood_preference():
 X=np.array([[1,-1],[1,-1]]);A=np.zeros_like(X);T=np.zeros_like(X);mu=np.full(2,.5);s=np.zeros(2,int)
 b0,b1=evidence(X,A,T,mu);l0=b0.copy();l1=b1.copy();l0[:,1]=-.2;l1[:,1]=-.1
 c0,c1,_=correct_evidence(X,A,T,mu,s,l0,l1,[])
 G=np.array([[1,0],[1,1]])
 ll=c0.sum(1)[:,None]+(c1-c0)@G.T
 assert np.array_equal(ll[:,0],ll[:,1])


def test_direct_coverage_overrides_unknown_accepted_call():
 X=np.array([[-1],[-1]]);A=np.array([[2],[0]]);T=np.array([[2],[2]]);mu=np.array([.5]);s=np.zeros(2,int)
 l0,l1=evidence(X,A,T,mu)
 r=[dict(sample_index=0,event_index=0,TPR=.1,FPR=.2,effective_mutant_validation_cells=0,effective_WT_validation_cells=0)]
 c0,c1,_=correct_evidence(X,A,T,mu,s,l0,l1,r)
 assert np.array_equal(c0,l0) and np.array_equal(c1,l1)
 assert c1[0,0]>c0[0,0] and c0[1,0]>c1[1,0]
