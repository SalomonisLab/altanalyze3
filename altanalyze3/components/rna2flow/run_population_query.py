"""Reproducible example: actual marrow ML/CLP protein profiles, fixed-polarity gates and FlowSOM placement."""
import json,urllib.request,time
from pathlib import Path
import pandas as pd
import argparse
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--url',default='http://127.0.0.1:8085');p.add_argument('--out',required=True);p.add_argument('--methods',default='kde_knn,kde_k5_average,kde_cellharmony,minmax_knn,zscore_knn');a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
def post(path,q):
 r=urllib.request.Request(a.url+path,data=json.dumps(q).encode(),headers={'Content-Type':'application/json'})
 return json.load(urllib.request.urlopen(r))
pops=['ML-1a','ML-1b','CLP1-a','CLP1-b','CLP1-c'];markers=['Sca-1','CD11c','CD11b','CD27','CD117','CD25','CD4','CD8','RNA:Ly6a','RNA:Itgax','RNA:Itgam','RNA:Spi1','RNA:Irf4','RNA:Irf8']
profiles=post('/api/marker_profiles',dict(populations=pops,markers=markers));pd.DataFrame(profiles['rows']).to_csv(out/'measured_CITE_profiles.tsv',sep='\t',index=False);(out/'profiles.json').write_text(json.dumps(profiles,indent=2))
signs={'CD4':'<=','CD8':'<=','CD117':'>','CD25':'<=','CD11c':'<=','CD11b':'<=','CD27':'>','Sca-1':'>'};rows=[];somrows=[]
for panel in ['cite_marrow_ADT195','cite_marrow_ADT112']:
 for method in a.methods.split(','):
  label='transfer_'+panel+'_'+method
  for pop in ['ML-1a','CLP1-a','CLP1-b','CLP1-c']:
   t=time.time()
   try:
    q=post('/api/optimize',dict(space='flow',label_set=label,populations=[pop],signs=signs));best=q['best'];q['label_set']=label
    (out/(label+'_'+pop+'.json')).write_text(json.dumps(q,indent=2))
    rows.append(dict(panel=panel,method=method,population=pop,**best['test'],n_gated_all=best['all_events']['n_gated'],shuffled_f1=best['shuffled_test_labels']['f1'],conditions=json.dumps(best['conditions']),seconds=time.time()-t))
    # All-event FlowSOM composition of a population label, independent of gate thresholds.
    p={'space':'flow','a':label,'b':'FlowSOM_J8DW'};ct=post('/api/concordance',p);counts=ct['counts'][ct['rows'].index(pop)];total=sum(counts)
    for cluster,n in zip(ct['cols'],counts):somrows.append(dict(panel=panel,method=method,population=pop,FlowSOM_cluster=cluster,n_events=n,percent_of_population=100*n/max(total,1)))
   except Exception as e:
    error=e.read().decode() if hasattr(e,'read') else str(e)
    rows.append(dict(panel=panel,method=method,population=pop,error=error))
   print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv(out/'constrained_gate_validation.tsv',sep='\t',index=False);pd.DataFrame(somrows).to_csv(out/'population_FlowSOM_placement.tsv',sep='\t',index=False)
