"""Add optimized KDE variants without replacing existing benchmark predictions."""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ.setdefault(v,'1')
import argparse,json,time
from pathlib import Path
import numpy as np,pandas as pd
from .data_api import FlowBundle
from ...rna2flow.crosswalk import build_crosswalk
from ...rna2flow.transfer import kde_knn
from ...rna2flow.marker_concordance import concordance
from ...rna2flow.evaluate import score_transfer


def build(bundle,flow_npz):
 b=FlowBundle(bundle);root=Path(bundle);arr=root/'arrays';man=b.man
 z=np.load(flow_npz);F=z['X'];flow_names=list(z['channels'])
 som=np.array(b.levels('FlowSOM_J8DW'))[b.labels('FlowSOM_J8DW')];rows=[]
 for name in ['cite_grimes','cite_chinese']:
  features=b.features(name);adt=[f for f in features if not f.startswith('RNA:')]
  pairs,_,_=build_crosswalk(flow_names,adt,extra_aliases={'CD8':'CD8a'});pairs=pairs.drop_duplicates('cite').drop_duplicates('flow').reset_index(drop=True)
  C=b.feature_matrix(name)[:,[features.index(c) for c in pairs.cite]];Fx=F[:,[flow_names.index(f) for f in pairs.flow]]
  for annotation,prefix in b.space(name)['transfer_links'].items():
   if annotation=='celltype':continue
   levels=b.levels(annotation,name);labs=np.array(levels)[b.labels(annotation,name)];keep=labs!='unassigned'
   for ties in ['average','first']:
    method='kde_k5_'+ties;t=time.time();lo,hi=(5,95) if name=='cite_chinese' else (1,99)
    pred=kde_knn(C[keep],Fx,labs[keep],k=5,ties=ties,lo=lo,hi=hi)
    key=prefix+method;cat=pd.Categorical(pred);fn=key+'.i16';cat.codes.astype(np.int16).tofile(arr/fn)
    man['spaces']['flow']['labels'][key]=dict(file=fn,levels=list(map(str,cat.categories)),source='transferred',source_space=name,source_annotation=annotation,method=method,k=5,clip=[lo,hi],ties=ties,n_source_cells=int(keep.sum()))
    tab,summary=concordance(C[keep],labs[keep],Fx,pred,pairs);summary.update(score_transfer(pred,som))
    summary.update(dataset=name,annotation=annotation,method=method,seconds=time.time()-t,n_source_cells=int(keep.sum()),n_source_labels=len(set(labs[keep])),validation_kind='Same-marker consistency, not independent validation')
    rows.append(summary);print(summary,flush=True)
    tab.to_csv(root/'validation_20261002'/f'consistency_{name}_{annotation}_{method}.tsv',sep='\t',index=False)
 pd.DataFrame(rows).to_csv(root/'validation_20261002'/'optimized_kde_consistency.tsv',sep='\t',index=False)
 man['audit']['optimized_kde_variants']=rows
 tmp=root/'manifest.new.json';tmp.write_text(json.dumps(man,indent=1));os.replace(tmp,root/'manifest.json')


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',required=True);p.add_argument('--flow-npz',required=True);a=p.parse_args();build(a.bundle,a.flow_npz)
if __name__=='__main__':main()
