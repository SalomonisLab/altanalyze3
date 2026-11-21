"""Shared finite mixtures with specimen-specific weights and explicit missing calls.

binary: discrete complete genotypes, fixed observation FP/FN, coordinate EM.
soft: empirical MUT-call propensities, Beta(1.5,1.5) regularization; not calibrated
genetic probabilities. Both use Dirichlet(1.5) specimen-weight regularization.
All input rows/events enter. Compression only groups identical sample/call rows.
"""
import numpy as np
from scipy.special import logsumexp

def compress(X, sample):
    X=np.asarray(X); sample=np.asarray(sample)
    assert X.ndim==2 and len(X)==len(sample) and np.isin(X,[-1,0,1]).all()
    u, inv, n=np.unique(np.column_stack([sample,X]),axis=0,return_inverse=True,return_counts=True)
    return u[:,1:],u[:,0],inv,n.astype(float)

def expectation(X, sample, p, weights):
    # einsum avoids the previously observed architecture-specific BLAS issue.
    ll=np.einsum('ij,kj->ik',(X==1).astype(float),np.log(p))
    ll+=np.einsum('ij,kj->ik',(X==0).astype(float),np.log1p(-p))
    logjoint=ll+np.log(weights[sample]); norm=logsumexp(logjoint,axis=1)
    return np.exp(logjoint-norm[:,None]),norm

def objective(norm,n,p,w,mode,errors=None):
    v=float(n@norm+0.5*np.log(w).sum())
    if mode=='soft':v+=float(0.5*(np.log(p)+np.log1p(-p)).sum())
    if mode=='learned':v+=float(.5*(np.log(errors)+np.log1p(-errors)).sum())
    return v

def fit(X,sample,K,mode,seed,fp=.1,fn=.1,maxiter=1000,tol=1e-7,fp_cap=.49,fn_cap=.49):
    assert mode in ('binary','soft','learned') and 0<fp<.5 and 0<fn<.5
    assert .001 <= fp_cap < .5 and .001 <= fn_cap < .5
    x,s,inv,n=compress(X,sample); S=int(np.max(sample))+1; J=X.shape[1]
    rng=np.random.default_rng(seed)
    mut=(x==1).astype(float); wt=(x==0).astype(float)
    rate=(np.einsum('i,ij->j',n,mut)+.5)/(np.einsum('i,ij->j',n,mut+wt)+1)
    chosen=rng.choice(len(x),K,replace=len(x)<K,p=n/n.sum())
    p=np.where(x[chosen]>=0,x[chosen],rate[None,:]).astype(float)
    p=np.clip(.8*p+.2*rng.uniform(.1,.9,size=p.shape),.02,.98)
    w=np.stack([rng.dirichlet(np.ones(K)*2) for _ in range(S)])
    # Soft warm start is initialization only; no inferred call replaces an input.
    phases=['soft',mode] if mode!='soft' else ['soft']
    errors=np.stack([np.full(J,min(fp,fp_cap)),np.full(J,max(1-fn,1-fn_cap))])
    trace=[]; converged=False
    for phase_i,phase in enumerate(phases):
        limit=25 if mode!='soft' and phase_i==0 else maxiter
        if phase!='soft':
            g=p>(fp+1-fn)/2
            p=np.where(g,errors[1],errors[0])
        prev=None
        for it in range(limit):
            R,norm=expectation(x,s,p,w); val=objective(norm,n,p,w,phase,errors)
            if prev is not None:
                assert val>=prev-1e-6,('EM objective decreased',val,prev)
                if (val-prev)/n.sum()<tol:
                    converged=phase_i==len(phases)-1
                    break
            prev=val
            weighted=R*n[:,None]
            a=np.einsum('ik,ij->kj',weighted,mut)
            b=np.einsum('ik,ij->kj',weighted,wt)
            if phase=='soft':p=(a+.5)/(a+b+1)
            else:
                score1=a*np.log(errors[1])+b*np.log1p(-errors[1])
                score0=a*np.log(errors[0])+b*np.log1p(-errors[0])
                # Ties preserve current state; no arbitrary MUT/WT switch.
                g=np.where(score1>score0+1e-10,True,np.where(score0>score1+1e-10,False,g))
                if phase=='learned':
                    errors[0]=np.clip(((a*(~g)).sum(0)+.5)/(((a+b)*(~g)).sum(0)+1),.001,fp_cap)
                    errors[1]=np.clip(((a*g).sum(0)+.5)/(((a+b)*g).sum(0)+1),1-fn_cap,.999)
                p=np.where(g,errors[1],errors[0])
            w=np.stack([weighted[s==v].sum(0)+.5 for v in range(S)])
            w/=w.sum(1,keepdims=True)
        trace.append(dict(phase=phase,iterations=it+1,objective=val))
    R,norm=expectation(x,s,p,w)
    G=g.astype(np.int8) if mode!='soft' else (p>=.5).astype(np.int8)
    return dict(p=p,w=w,G=G,errors=errors,posterior=R[inv],loglik=float(n@norm),
                objective=objective(norm,n,p,w,mode,errors),converged=converged,trace=trace)

def predict(X,sample,model):
    x,s,inv,n=compress(X,sample)
    R,_=expectation(x,s,model['p'],model['w'])
    return np.einsum('ik,kj->ij',R,model['p'])[inv]

def scores(prob,truth,mask):
    y=truth[mask];p=np.clip(prob[mask],1e-12,1-1e-12)
    assert len(y)>0 and np.isin(y,[0,1]).all()
    return dict(n=int(len(y)),brier=float(np.mean((p-y)**2)),
                logloss=float(-np.mean(y*np.log(p)+(1-y)*np.log1p(-p))),
                accuracy=float(np.mean((p>=.5)==y)))

def baseline(train,sample):
    out=np.empty(train.shape,float)
    for s in np.unique(sample):
        ix=sample==s
        out[ix]=((train[ix]==1).sum(0)+.5)/((train[ix]>=0).sum(0)+1)
    return out
