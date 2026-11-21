"""Experimental gain-only phylogeny mixture; gene names never enter fitting.

One edge per event defines a complete event tree. Some ancestral states have
zero sampled prevalence, so a detectable clone can acquire multiple events.
All cells/events enter the likelihood. Original calls remain immutable.
"""
import numpy as np
from scipy.special import betaln, gammaln, logsumexp
from scipy.optimize import minimize


def evidence(X,A,T,mutant_mean,mode='counts',concentration=2.,fp=.01,fn=.30):
    X,A,T=np.asarray(X),np.asarray(A),np.asarray(T)
    assert X.shape==A.shape==T.shape and np.isin(X,[-1,0,1]).all()
    assert np.isfinite(A).all() and np.isfinite(T).all() and ((A>=0)&(T>=A)).all()
    assert (A==np.floor(A)).all() and (T==np.floor(T)).all()
    assert 0<fp<.5 and 0<fn<1 and concentration>0
    q=np.asarray(mutant_mean);assert q.shape==(X.shape[1],) and ((q>0)&(q<1)).all()
    l0=np.zeros(X.shape,float);l1=l0.copy()
    for call in [0,1]:
        mask=X==call
        l0[mask]=np.log(fp if call else 1-fp)
        l1[mask]=np.log(1-fn if call else fn)
    if mode=='counts':
        direct=T>0
        choose=gammaln(T+1)-gammaln(A+1)-gammaln(T-A+1)
        # Fixed observation assumptions, not fitted/empirically calibrated rates.
        for dest,mu,c in [(l0,np.full(X.shape[1],fp),100.),(l1,q,concentration)]:
            alpha=mu*c;beta=(1-mu)*c
            ll=choose+betaln(A+alpha,T-A+beta)-betaln(alpha,beta)
            dest[direct]=ll[direct]
    else:assert mode=='binary'
    assert np.isfinite(l0).all() and np.isfinite(l1).all()
    return l0,l1


def genotypes(parents):
    """Node 0 normal, node j+1 gains event j; arbitrary valid label ordering."""
    p=np.asarray(parents,dtype=int);n=len(p)
    assert p[0]==-1 and ((p[1:]>=0)&(p[1:]<n)).all()
    G=np.zeros((n,n-1),dtype=np.int8)
    for k in range(1,n):
        node=k;seen=set()
        while node!=0:
            assert node not in seen,'cycle';seen.add(node)
            G[k,node-1]=1;node=p[node]
    return G


def detector_evidence(X,A,T,mutant_mean,sample,Z,available,concentration=2.,fp=.01,fn=.30,strength=2.):
    """Partially pooled surrogate-detection probabilities, calibrated on covered cells.

    Absence of a positive surrogate is an observation of its detection process,
    NOT a negative genotype call. The original X matrix remains unchanged.
    Direct cells enter the genotype likelihood via counts once; their surrogate
    mask contributes only to calibration. No pseudo-read counts are constructed.
    """
    assert strength>0 and Z.shape==X.shape and available.shape==X.shape
    assert np.isin(Z,[0,1]).all() and not ((Z==1)&~available).any()
    l0,l1=evidence(X,A,T,mutant_mean,'counts',concentration,fp,fn)
    S=int(np.max(sample))+1;J=X.shape[1];stats=np.zeros((S,J,4));prevalence=np.full((S,J),np.nan)
    for s in range(S):
        for j in range(J):
            ix=(sample==s)&(T[:,j]>0)&available[:,j]
            if not ix.any():continue
            ll=np.column_stack([l0[ix,j],l1[ix,j]]);off=ll.max(1);L=np.exp(ll-off[:,None])
            f=mixture_weights(L,np.ones(ix.sum()),np.zeros(ix.sum(),int),1)
            w=f['weights'][0];q=(L[:,1]*w[1])/np.maximum(L@w,1e-250)
            stats[s,j]=[(q*Z[ix,j]).sum(),q.sum(),((1-q)*Z[ix,j]).sum(),(1-q).sum()]
            prevalence[s,j]=w[1]
    pooled=stats.sum(axis=0)
    pooled_tpr=(pooled[:,0]+.5)/(pooled[:,1]+1.)
    pooled_fpr=(pooled[:,2]+.5)/(pooled[:,3]+1.)
    records=[]
    for s in range(S):
        for j in range(J):
            if not available[sample==s,j].any():continue
            x=stats[s,j]
            tpr=(x[0]+strength*pooled_tpr[j])/(x[1]+strength)
            fpr=(x[2]+strength*pooled_fpr[j])/(x[3]+strength)
            assert 0<tpr<1 and 0<fpr<1
            ix=(sample==s)&(T[:,j]==0)&available[:,j]
            # This replaces, rather than double-counts, the categorical MUT call
            # for cells that carry a validated surrogate in the union mask.
            l0[ix,j]=np.where(Z[ix,j],np.log(fpr),np.log1p(-fpr))
            l1[ix,j]=np.where(Z[ix,j],np.log(tpr),np.log1p(-tpr))
            # Other accepted channels (DEG/DEI and explicit imputed WT) stay present.
            extra=ix&(Z[:,j]==0)&(X[:,j]>=0)
            l0[extra,j]+=np.where(X[extra,j]==1,np.log(fp),np.log1p(-fp))
            l1[extra,j]+=np.where(X[extra,j]==1,np.log1p(-fn),np.log(fn))
            records.append(dict(sample_index=s,event_index=j,TPR=tpr,FPR=fpr,
                                effective_mutant_validation_cells=x[1],effective_WT_validation_cells=x[3],
                                positive_effective_mutants=x[0],positive_effective_WT=x[2],
                                direct_only_mutant_fraction=prevalence[s,j],prior_strength=strength,
                                pooled_TPR=pooled_tpr[j],pooled_FPR=pooled_fpr[j],
                                total_cells=int((sample==s).sum()),surrogate_positive_cells=int(Z[sample==s,j].sum())))
    return l0,l1,records


def mixture_weights(L,n,sample,S,initial=None,tol=1e-8):
    """Concave likelihood over a simplex, with a directional KKT check."""
    K=L.shape[1];w=np.empty((S,K));ll=0.;gaps=[]
    for s in range(S):
        ix=sample==s;F=L[ix];v=n[ix];N=v.sum();freq=v/N
        start=np.full(K,1/K) if initial is None else np.maximum(initial[s],1e-10)
        start/=start.sum()
        def fun(a):
            den=np.maximum(F@a,1e-250)
            return -float(freq@np.log(den)),-(F.T@(freq/den))
        if K==1:a=np.ones(1)
        else:
            opt=minimize(fun,start,jac=True,method='SLSQP',bounds=[(0.,1.)]*K,
                         constraints={'type':'eq','fun':lambda a:a.sum()-1,'jac':lambda a:np.ones(K)},
                         options={'ftol':tol,'maxiter':150})
            a=np.clip(opt.x,0,1);a/=a.sum()
        val,neggrad=fun(a);gap=float(-neggrad.min()-1.)
        # Multiplicative EM repairs finite-tolerance SLSQP gaps without clipping data.
        for it in range(1000):
            if gap<2e-5:break
            a*= -neggrad;a/=a.sum();val,neggrad=fun(a);gap=float(-neggrad.min()-1.)
        # Gap is reported for every tree; uncertified fits cannot be final winners.
        w[s]=a;ll-=N*val;gaps.append(gap)
    return dict(weights=w,loglik=ll,max_gap=max(gaps))


class TreeSearch:
    def __init__(self,l0,l1,sample,prune=True):
        self.J=l0.shape[1];self.S=int(np.max(sample))+1;self.prune=prune
        # Exact compression of identical likelihood and specimen rows, no rounding.
        raw=np.column_stack([sample,l0,l1]);u,inv,n=np.unique(raw,axis=0,return_inverse=True,return_counts=True)
        self.sample=u[:,0].astype(int);self.l0=u[:,1:1+self.J];self.l1=u[:,1+self.J:]
        self.inverse=inv;self.n=n.astype(float);self.N=self.n.sum()
        self.offset=np.maximum(self.l0,self.l1).sum(1);self.constant=self.l0.sum(1)-self.offset
        self.delta=self.l1-self.l0;self.cache={};self.evaluations=0
    def evaluate(self,p):
        key=tuple(p)
        if key in self.cache:return self.cache[key]
        G=genotypes(p)
        raw=self.constant[:,None]+self.delta@G.T
        L=np.exp(np.maximum(raw,-575))
        fit=mixture_weights(L,self.n,self.sample,self.S)
        active=np.arange(len(G));best=fit
        # BIC parameter count uses the full specimen simplex, not a tuned positive-call cutoff.
        def penalized(f,k):return f['loglik']-.5*self.S*(k-1)*np.log(self.N)
        score=penalized(fit,len(active))
        if self.prune:
            # Delete a sampled state, never an event/tree node. WT stays available.
            while len(active)>1:
                proposals=[]
                for at in range(1,len(active)):
                    keep=np.delete(active,at);initial=np.delete(best['weights'],at,axis=1)
                    f=mixture_weights(L[:,keep],self.n,self.sample,self.S,initial)
                    proposals.append((penalized(f,len(keep)),keep,f))
                v,keep,f=max(proposals,key=lambda x:x[0])
                if v<=score+1e-6:break
                score,active,best=v,keep,f
        result=dict(parents=np.array(p),G=G,active=active,score=score,
                    loglik=best['loglik']+float(self.n@self.offset),weights=best['weights'],
                    max_gap=best['max_gap'])
        self.cache[key]=result;self.evaluations+=1
        return result
    def neighbors(self,p,rng):
        proposals=[];n=len(p)
        # Reattachment permits ancestry changes without recurrent mutation gains.
        for node in range(1,n):
            for parent in range(n):
                if parent in [node,p[node]]:continue
                q=p.copy();q[node]=parent
                try:genotypes(q)
                except AssertionError:continue
                proposals.append(q)
        # Swap event labels, including ancestor/descendant positions.
        for a in range(1,n):
            for b in range(a+1,n):
                perm=np.arange(n);perm[a],perm[b]=b,a
                q=np.full(n,-1,dtype=int)
                for child in range(1,n):q[perm[child]]=perm[p[child]]
                proposals.append(q)
        rng.shuffle(proposals)
        return proposals
    def run(self,seed,max_rounds=12):
        rng=np.random.default_rng(seed);order=rng.permutation(np.arange(1,self.J+1))
        p=np.zeros(self.J+1,dtype=int);p[0]=-1
        # Random recursive tree has no input from gene labels or expected genotypes.
        for k,node in enumerate(order):p[node]=rng.choice(np.r_[0,order[:k]])
        best=self.evaluate(p);trace=[];converged=False
        for r in range(max_rounds):
            before=best['score'];candidate=best
            for q in self.neighbors(best['parents'],rng):
                f=self.evaluate(q)
                if f['score']>candidate['score']+1e-6:candidate=f
            best=candidate;trace.append(dict(round=r,score=best['score'],evaluations=self.evaluations))
            print('round',r,'score',best['score'],'trees',self.evaluations,'K',len(best['active']),flush=True)
            if best['score']<=before+1e-6:converged=True;break
        best=dict(best);best.update(converged=converged,trace=trace,evaluations=self.evaluations,seed=seed)
        return best
    def posterior(self,fit):
        G=fit['G'][fit['active']]
        ll=self.constant[:,None]+self.delta@G.T
        joint=ll+np.log(np.maximum(fit['weights'][self.sample],1e-250))
        P=np.exp(joint-logsumexp(joint,axis=1)[:,None])
        return P[self.inverse],ll[self.inverse]
