"""Adversarial transfer audit on antibodies excluded from label assignment.

MarkerFinder ranks on assignment markers are a consistency check, not independent
validation. This audit excludes a prespecified alternating half of shared antibodies
from each classifier, and evaluates their population rankings after transfer. It
retains the complete CITE label universe, including absent/too-small flow states.
No result establishes functional stem-cell identity or biological replication.
"""
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ.setdefault(v,'1')
import argparse,json,time
from pathlib import Path
import numpy as np,pandas as pd
from .transfer import METHODS,kde_knn
from .crosswalk import build_crosswalk
from .marker_concordance import concordance
from .evaluate import score_transfer
from ..visualization.flow_viewer.data_api import FlowBundle


def audit(bundle_root,flow_npz,out_dir,spaces=None,method_names=None,annotation_names=None):
    b=FlowBundle(bundle_root);root=Path(out_dir);root.mkdir(parents=True,exist_ok=True)
    # The audited RDS matrix is the exact input used by the original benchmarks.
    z=np.load(flow_npz);F=z['X'];flow_names=list(z['channels'])
    som=np.array(b.levels('FlowSOM_J8DW'))[b.labels('FlowSOM_J8DW')]
    methods=dict(METHODS)
    methods['kde_k5_average']=lambda c,f,l,**kw:kde_knn(c,f,l,k=5,ties='average',**kw)
    methods['kde_k5_first']=lambda c,f,l,**kw:kde_knn(c,f,l,k=5,ties='first',**kw)
    if method_names:methods={k:methods[k] for k in method_names}
    rows=[];population_rows=[]
    for space in spaces or ['cite_grimes','cite_chinese']:
        allfeatures=b.features(space);adt=[f for f in allfeatures if not f.startswith('RNA:')]
        pairs,_,_=build_crosswalk(flow_names,adt,extra_aliases={'CD8':'CD8a'})
        pairs=pairs.drop_duplicates('cite').drop_duplicates('flow').reset_index(drop=True)
        C=b.feature_matrix(space)[:,[allfeatures.index(c) for c in pairs.cite]]
        Fx=F[:,[flow_names.index(f) for f in pairs.flow]]
        annotations=annotation_names or (['Population'] if space.startswith('cite_marrow') else ['Mm-MarrowAtlas-L4','pruned'])
        for annotation in annotations:
            if annotation not in b.space(space)['labels']:continue
            levels=b.levels(annotation,space);codes=b.labels(annotation,space)
            labs=np.array(levels)[codes];keep=labs!='unassigned';c=np.asarray(C[keep]);l=labs[keep]
            for fold in [0,1]:
                held=np.arange(len(pairs))[np.arange(len(pairs))%2==fold]
                train=np.array([i for i in range(len(pairs)) if i not in held])
                for method,fn in methods.items():
                    start=time.time();row=dict(dataset=space,annotation=annotation,fold=fold,method=method,
                        n_assignment_markers=len(train),n_validation_markers=len(held),
                        assignment_markers=','.join(pairs.flow.iloc[train]),held_out_markers=','.join(pairs.flow.iloc[held]))
                    try:
                        kw={'lo':5.,'hi':95.} if method.startswith('kde_k5') and space=='cite_chinese' else {}
                        pred=fn(c[:,train],Fx[:,train],l,**kw)
                        row.update(score_transfer(pred,som))
                        table,summary=concordance(c[:,held],l,Fx[:,held],pred,pairs.iloc[held]);row.update(summary)
                        lookup=table.set_index('population').to_dict('index') if len(table) else {}
                        for population in sorted(set(l)):
                            nref=int((l==population).sum());nflow=int((pred==population).sum())
                            result=lookup.get(population,{});population_rows.append(dict(dataset=space,annotation=annotation,fold=fold,method=method,population=population,n_reference=nref,n_flow=nflow,status='compared' if result else ('absent_in_flow' if nflow==0 else 'insufficient_cells'),**result))
                        row['n_reference_labels']=len(set(l));row['n_labels_absent_in_flow']=len(set(l)-set(pred));row['seconds']=time.time()-start
                    except Exception as e:row['error']=str(e)
                    rows.append(row);pd.DataFrame(rows).to_csv(root/'held_out_marker_summary.tsv',sep='\t',index=False)
                    pd.DataFrame(population_rows).to_csv(root/'held_out_marker_populations.tsv',sep='\t',index=False)
                    print(json.dumps(row),flush=True)
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',required=True);p.add_argument('--flow-npz',required=True);p.add_argument('--out',required=True);p.add_argument('--space',action='append');p.add_argument('--annotation',action='append');p.add_argument('--methods',help='Comma-separated methods');a=p.parse_args();audit(a.bundle,a.flow_npz,a.out,a.space,a.methods.split(',') if a.methods else None,a.annotation)
if __name__=='__main__':main()
