"""Bounded-memory implementation of the notebook's five-stage cellHarmony transfer.

The notebook defines independent Louvain communities, reciprocal community matches,
Pearson-nearest reference cells within matched communities, query label centroids,
then Pearson-nearest query centroids. Unit-vector Euclidean nearest neighbors are
algebraically equivalent to Pearson argmax and avoid full cell × cell matrices.
igraph's multilevel Louvain backend is used rather than notebook vtraag-louvain;
community stochasticity/backend and limited shared-marker PCA are recorded explicitly.
"""
import numpy as np
from sklearn.neighbors import NearestNeighbors
from .transfer import _scale_both
from .normalize import kde_quantile_map,reference_spline


def unit_rows(X):
    X=np.asarray(X,dtype=np.float64);Z=X-X.mean(1,keepdims=True);norm=np.linalg.norm(Z,axis=1,keepdims=True)
    return np.divide(Z,norm,out=np.zeros_like(Z),where=norm>0)


def nearest_pearson(reference,query,batch_size=1024):
    ref=unit_rows(reference);q=unit_rows(query)
    nn=NearestNeighbors(n_neighbors=1,metric='euclidean',algorithm='brute',n_jobs=1).fit(ref)
    matches=np.empty(len(q),np.int64);rhos=np.empty(len(q))
    for start in range(0,len(q),batch_size):
        _,idx=nn.kneighbors(q[start:start+batch_size]);idx=idx[:,0];matches[start:start+len(idx)]=idx
        rhos[start:start+len(idx)]=np.einsum('ij,ij->i',q[start:start+len(idx)],ref[idx])
    return matches,rhos


def communities(X,seed=0,resolution=1):
    import scanpy as sc,anndata as ad,igraph as ig,random
    a=ad.AnnData(np.asarray(X,np.float32));npcs=min(30,X.shape[1]-1,len(X)-1)
    if npcs<2:raise ValueError('At least three shared markers are required')
    sc.tl.pca(a,n_comps=npcs,random_state=seed)
    sc.pp.neighbors(a,n_neighbors=min(15,len(X)-1),n_pcs=npcs,random_state=seed)
    graph=a.obsp['connectivities'].tocoo();keep=graph.row<graph.col
    g=ig.Graph(n=len(X),edges=list(zip(graph.row[keep],graph.col[keep])),directed=False)
    ig.set_random_number_generator(random.Random(seed))
    return np.asarray(g.community_multilevel(weights=graph.data[keep].tolist(),resolution=resolution).membership)


def assign_communities(C,F,labels,ref_groups,query_groups):
    """Stages 2–5, exposed separately to test against the notebook with frozen groups."""
    C=np.asarray(C,float);F=np.asarray(F,float);labels=np.asarray(labels).astype(str)
    z=lambda x:(x-x.mean(0))/np.where(x.std(0)>0,x.std(0),1)
    C,F=z(C),z(F);rlevels=np.unique(ref_groups);qlevels=np.unique(query_groups)
    rcent=np.vstack([C[ref_groups==g].mean(0) for g in rlevels]);qcent=np.vstack([F[query_groups==g].mean(0) for g in qlevels])
    correlations=unit_rows(qcent)@unit_rows(rcent).T
    mapping={i:{int(np.argmax(correlations[i]))} for i in range(len(qlevels))}
    for j in range(len(rlevels)):mapping[int(np.argmax(correlations[:,j]))].add(j)
    initial=np.full(len(F),'unassigned',object);cell_rho=np.zeros(len(F))
    for i,qg in enumerate(qlevels):
        qi=np.flatnonzero(query_groups==qg);ri=np.flatnonzero(np.isin(ref_groups,rlevels[list(mapping[i])]))
        matches,rhos=nearest_pearson(C[ri],F[qi]);initial[qi]=labels[ri[matches]];cell_rho[qi]=rhos
    levels=np.unique(initial);centroids=np.vstack([F[initial==l].mean(0) for l in levels])
    # Only O(batch × number of labels), rather than query × reference correlations.
    normalized_query,normalized_centroids=unit_rows(F),unit_rows(centroids)
    idx=np.empty(len(F),int);final_rho=np.empty(len(F))
    for start in range(0,len(F),1024):
        stop=min(start+1024,len(F));scores=normalized_query[start:stop]@normalized_centroids.T
        ii=np.argmax(scores,axis=1);idx[start:stop]=ii;final_rho[start:stop]=scores[np.arange(stop-start),ii]
    return levels[idx],dict(initial_labels=initial,reference_cell_rho=cell_rho,final_rho=final_rho,
                            reference_communities=len(rlevels),query_communities=len(qlevels))


def kde_cellharmony(cite,flow,labels,seed=0,lo=1,hi=99,n_ref=20000,return_details=False,**kwargs):
    C,F=_scale_both(cite,flow,lo,hi);rng=np.random.default_rng(seed);mapped=np.empty_like(C)
    for j in range(C.shape[1]):
        ref=F[:,j]
        if len(ref)>n_ref:ref=ref[rng.choice(len(ref),n_ref,replace=False)]
        mapped[:,j]=kde_quantile_map(C[:,j],None,spline=reference_spline(ref),ties='first')
    r=communities(mapped,seed);q=communities(F,seed)
    pred,details=assign_communities(mapped,F,labels,r,q)
    details['community_backend']='igraph multilevel Louvain; notebook uses vtraag Louvain; PCA min(30,n_markers-1)'
    return (pred,details) if return_details else pred
