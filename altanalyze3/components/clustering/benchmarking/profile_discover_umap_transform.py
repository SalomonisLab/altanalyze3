"""Synthetic UMAP transform profiling only; no biological analysis or model changes.

Run in its own process with NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
OMP_NUM_THREADS=4. Writes /tmp/profile_discover_umap_transform_20261007.json.
The complete 827-feature input is generated from seed 20261007; 30,000 fit
rows evenly represent all 51 synthetic groups. This profiles existing UMAP,
not a replacement method or the visitor dataset.
"""
import json,time,resource,sys,os
import numpy as np
import numba
from umap import UMAP
import umap.umap_ as implementation
numba.set_num_threads(4)
# Synthetic computational fixture; it is not a substitute for the visitor input.
rng=np.random.default_rng(20261007)
features,fit_cells,query_cells,states=827,30000,12000,51
centers=rng.random((states,features),dtype=np.float32)
centers[centers<.85]=0
labels=np.arange(fit_cells+query_cells)%states
x=centers[labels].copy()
x+=rng.random(x.shape,dtype=np.float32)*.2
x[rng.random(x.shape,dtype=np.float32)<.65]=0
np.log1p(x*20,out=x)
model=UMAP(n_neighbors=15,min_dist=.75,metric='correlation',random_state=0)
start=time.perf_counter();model.fit_transform(x[:fit_cells]);fit_seconds=time.perf_counter()-start
print(json.dumps({'fit_seconds':fit_seconds,'pid':os.getpid()}),flush=True)
index=model._knn_search_index
phases=[]
query_original=index.query
for name in ['smooth_knn_dist','compute_membership_strengths','init_graph_transform','optimize_layout_euclidean']:
 original=getattr(implementation,name)
 def wrapper(*args,_name=name,_original=original,**kwargs):
  start=time.perf_counter();result=_original(*args,**kwargs);phases.append((_name,time.perf_counter()-start));return result
 setattr(implementation,name,wrapper)
def query(*args,**kwargs):
 start=time.perf_counter();result=query_original(*args,**kwargs);phases.append(('neighbor_search',time.perf_counter()-start));return result
index.query=query
reports=[]
for attempt in range(2):
 phases.clear();start=time.perf_counter();coordinates=model.transform(x[fit_cells:]);elapsed=time.perf_counter()-start
 assert coordinates.shape==(query_cells,2) and np.isfinite(coordinates).all()
 reports.append({'attempt':attempt+1,'total_seconds':elapsed,'phases':dict(phases),'finite_cells':len(coordinates)})
 print(json.dumps(reports[-1]),flush=True)
report={'scope':'Synthetic UMAP-only diagnostic, not the reported visitor dataset or full workflow','fit_cells':fit_cells,'query_cells':query_cells,'features':features,'states':states,'metric':'correlation','n_neighbors':15,'min_dist':.75,'random_state':0,'numba_threads':numba.get_num_threads(),'fit_seconds':fit_seconds,'transforms':reports,'peak_rss_gib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**3 if sys.platform=='darwin' else 1024**2)}
open('/tmp/profile_discover_umap_transform_20261007.json','w').write(json.dumps(report,indent=2)+'\n')
