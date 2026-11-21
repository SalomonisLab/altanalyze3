"""Exhaustive finite-genotype mixture: binary events, explicit missingness, no fixed clone list.

This is a custom comparator, NOT native scBayes. For genotype g and observed calls x,
L(x|g)=product_j P(x_j|g_j), omitting unknown entries. All 2**J configurations
are retained. Sample weights maximize sum_i log(sum_g w_g L_i,g), w>=0,sum(w)=1.
Identical observed rows are sufficient-statistic compressed without dropping cells.
The objective -sum(q*log(Lw))+sum(w) is convex; its optimum has sum(w)=1.
"""
import numpy as np
from scipy.optimize import minimize,linprog
from scipy.linalg import svd

def configurations(n_events):
 if not 1<=n_events<=16:raise ValueError('Exhaustive enumeration supports 1..16 events; no implicit feature reduction')
 return ((np.arange(2**n_events)[:,None]>>np.arange(n_events))&1).astype(np.int8)

def likelihood(patterns,G,fp,fn):
 if not 0<fp<.5 or not 0<fn<.5:raise ValueError('Fixed FP/FN must be between zero and .5')
 patterns=np.asarray(patterns);G=np.asarray(G)
 if not np.isin(patterns,[-1,0,1]).all() or not np.isin(G,[0,1]).all():raise ValueError('Invalid calls/genotypes')
 p=fp+(1-fn-fp)*G
 ll=(patterns==1)@np.log(p).T+(patterns==0)@np.log1p(-p).T
 reference=np.zeros(ll.shape)
 for j in range(G.shape[1]):reference+=np.where(patterns[:,j,None]==1,np.log(p[:,j]),np.where(patterns[:,j,None]==0,np.log1p(-p[:,j]),0))
 if not np.allclose(ll,reference,atol=1e-10,rtol=0):raise RuntimeError('Numerical backend failed independent likelihood gate; do not fit on this BLAS configuration')
 offset=ll.max(axis=1)
 return np.exp(ll-offset[:,None]),offset

def fit(X,fp=.01,fn=.3,seed=0,maxiter=20000):
 X=np.asarray(X)
 if X.ndim!=2 or len(X)==0 or not np.isin(X,[-1,0,1]).all():raise ValueError('Invalid observation matrix')
 G=configurations(X.shape[1]);patterns,inverse,counts=np.unique(X,axis=0,return_inverse=True,return_counts=True)
 L,offset=likelihood(patterns,G,fp,fn);q=counts/counts.sum()
 w0=np.full(len(G),1/len(G)) if seed==0 else np.random.default_rng(seed).dirichlet(np.ones(len(G)))
 def fun(w):
  v=L@w
  if (v<=0).any():return np.inf,np.full(len(w),-np.inf)
  return -q@np.log(v)+w.sum(),1-L.T@(q/v)
 result=minimize(fun,w0,jac=True,method='L-BFGS-B',bounds=[(0,None)]*len(G),options={'maxiter':maxiter,'ftol':1e-15,'gtol':1e-8,'maxls':50,'maxcor':20})
 w=np.maximum(result.x,0);w/=w.sum();val,grad=fun(w)
 # Dual directional gap on the simplex bounds suboptimality for this convex objective.
 gap=float(grad@w-grad.min())
 polish=0
 while np.isfinite(gap) and gap>1e-6 and polish<50000:
  # Frank-Wolfe vertex step with exact convex one-dimensional minimization.
  # Unlike a boundary-stalled L-BFGS iteration, this can revive any of 2**J states.
  j=int(np.argmin(grad));v=L@w;d=L[:,j]-v
  lo=0.;hi=1.
  if -q@(d/(v+d))<=0:step=1.
  else:
   for _ in range(45):
    step=(lo+hi)/2
    if -q@(d/(v+step*d))<0:lo=step
    else:hi=step
   step=(lo+hi)/2
  oldval=val;w*=1-step;w[j]+=step
  # Multiplicative EM updates retain all states and cannot lower the likelihood.
  for _ in range(5):
   w*=L.T@(q/(L@w));w/=w.sum()
  val,grad=fun(w);gap=float(grad@w-grad.min());polish+=1
  if val>oldval+1e-12:raise RuntimeError('Monotonicity check failed')
 if not np.isfinite(gap) or gap>1e-5:raise RuntimeError(f'Insufficient optimization: dual gap={gap}; {result.message}; polish={polish}')
 P=L*w[None,:];P/=P.sum(axis=1)[:,None]
 return dict(G=G,patterns=patterns,inverse=inverse,counts=counts,L=L,offset=offset,w=w,P=P,loglik=float(counts@(np.log(L@w)+offset)),gap=gap,iterations=result.nit,polish_iterations=polish,message=str(result.message),fp=fp,fn=fn)

def equivalent_bounds(fit,objective,tolerance=1e-9):
 """Linear bounds holding each observed-pattern likelihood fixed within tolerance.
 Not sampling confidence intervals. Endpoints retain all configurations and data.
 """
 L=fit['L'];v=L@fit['w'];c=np.asarray(objective,dtype=float)
 # An orthonormal basis removes only numerical row dependencies, not events.
 # Every original pattern constraint is independently verified after solving.
 C=np.vstack([np.ones((1,L.shape[1])),L])
 if '_equivalence_basis' not in fit:
  _,sing,Vh=svd(C,full_matrices=False)
  fit['_equivalence_basis']=Vh[sing>sing[0]*1e-12]
 Aeq=fit['_equivalence_basis'];beq=Aeq@fit['w']
 endpoints=[]
 for direction in [1,-1]:
  r=linprog(direction*c,A_eq=Aeq,b_eq=beq,bounds=(0,None),method='highs-ds',options={'presolve':False,'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
  if not r.success:raise RuntimeError('Likelihood-equivalence LP failed: '+r.message)
  residual=float(np.max(np.abs(L@r.x-v)))
  if residual>5e-9 or abs(r.x.sum()-1)>1e-8:raise RuntimeError('LP constraint validation failed')
  ll=float(fit['counts']@(np.log(L@r.x)+fit['offset']))
  if abs(ll-fit['loglik'])/fit['counts'].sum()>1e-6:raise RuntimeError('LP likelihood validation failed')
  endpoints.append(dict(value=float(c@r.x),max_pattern_difference=residual,loglik_difference=ll-fit['loglik'],weights=r.x))
 return endpoints
