"""Corrected graph-width angular-search diagnostic against exhaustive correlation.

Synthetic 30,000-fit/12,000-query/827-feature fixture only. This prototype was
rejected: reusing a graph without a search forest recovered ~21% of exact
neighbors. It is not evidence against all angular ANN implementations.
Run in a fresh process with NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4.
"""
from pathlib import Path
import json


def main():
    approval = json.loads(Path(__file__).with_name("discover_umap_transform_profile_20261007.json").read_text())["decision"]
    if approval["state"] != "approved" or "Yes, do that" not in approval["source"]:
        raise ValueError("Search-backend comparison requires the explicit human authorization.")
    import os,time,json,resource,sys
    import numpy as np
    from umap import UMAP
    from pynndescent import NNDescent
    from threadpoolctl import threadpool_limits
    threadpool_limits(4)
    rng=np.random.default_rng(20261007)
    f,n,q,g=827,30000,12000,51
    centers=rng.random((g,f),dtype=np.float32);centers[centers<.85]=0
    labels=np.arange(n+q)%g
    x=centers[labels].copy();x+=rng.random(x.shape,dtype=np.float32)*.2;x[rng.random(x.shape,dtype=np.float32)<.65]=0;np.log1p(x*20,out=x)
    m=UMAP(n_neighbors=15,min_dist=.75,metric='correlation',random_state=0)
    t=time.perf_counter();m.fit_transform(x[:n]);print('fit',time.perf_counter()-t,flush=True)
    original=m._knn_search_index
    class Exact:
     _angular_trees=False
     def __init__(self,train):
      a=np.array(train,dtype=np.float32,order='C',copy=True);a-=a.mean(axis=1,keepdims=True);a/=np.maximum(np.linalg.norm(a,axis=1,keepdims=True),1e-30);self.a=a
     def query(self,data,k,epsilon=.12):
      a=np.array(data,dtype=np.float32,order='C',copy=True);a-=a.mean(axis=1,keepdims=True);a/=np.maximum(np.linalg.norm(a,axis=1,keepdims=True),1e-30)
      idx=np.empty((len(a),k),dtype=np.int32);d=np.empty((len(a),k),dtype=np.float32)
      for i in range(0,len(a),256):
       sim=a[i:i+256]@self.a.T
       ix=np.argpartition(sim,-k,axis=1)[:,-k:];dd=1-np.take_along_axis(sim,ix,axis=1);order=np.argsort(dd,axis=1)
       idx[i:i+256]=np.take_along_axis(ix,order,axis=1);d[i:i+256]=np.take_along_axis(dd,order,axis=1)
      return idx,np.maximum(d,0)
    e=Exact(x[:n]);reports={};outputs={}
    for name,index in [('exact',e)]:
     m._knn_search_index=index;t=time.perf_counter();out=m.transform(x[n:]);elapsed=time.perf_counter()-t;outputs[name]=out
     reports[name]={'seconds':elapsed,'finite':bool(np.isfinite(out).all())};print(name,reports[name],flush=True)
    t=time.perf_counter();new=NNDescent(e.a,metric='cosine',n_neighbors=original.neighbor_graph[0].shape[1],init_graph=original.neighbor_graph[0],n_iters=0,random_state=0,n_jobs=4);new.prepare();print('cosine_index',time.perf_counter()-t,flush=True)
    class Angular:
     _angular_trees=original._angular_trees
     def query(self,data,k,epsilon=.12):
      a=np.array(data,dtype=np.float32,order='C',copy=True);a-=a.mean(axis=1,keepdims=True);a/=np.maximum(np.linalg.norm(a,axis=1,keepdims=True),1e-30)
      return new.query(a,k=k,epsilon=epsilon)
    m._knn_search_index=Angular();t=time.perf_counter();out=m.transform(x[n:]);reports['angular']={'seconds':time.perf_counter()-t,'finite':bool(np.isfinite(out).all())};print('angular',reports['angular'],flush=True)
    anchors=np.arange(0,q,50);target=e.query(x[n:][anchors],15)[0]
    for name,index in [('angular',Angular())]:
     ix=index.query(x[n:][anchors],15)[0];reports[name]['exact_neighbor_recall']=float(np.mean([len(set(a)&set(b))/15 for a,b in zip(ix,target)]))
    reports['peak_rss_gib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**3 if sys.platform=='darwin' else 1024**2)
    open('/tmp/benchmark_umap_search_corrected_angular.json','w').write(json.dumps(reports,indent=2));print(reports,flush=True)


if __name__ == "__main__":
    main()
