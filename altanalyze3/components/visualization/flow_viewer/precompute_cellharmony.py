"""Add KDE + community/correlation cellHarmony predictions and correlation diagnostics."""
import os
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,json,time
from pathlib import Path
import numpy as np,pandas as pd
from .data_api import FlowBundle
from ...rna2flow.cellharmony import kde_cellharmony
from ...rna2flow.crosswalk import build_crosswalk
from ...rna2flow.marker_concordance import concordance
from ...rna2flow.evaluate import score_transfer

def build(bundle,flow_npz,spaces=None,annotations=None):
    root=Path(bundle);b=FlowBundle(root);man=b.man;arr=root/'arrays';z=np.load(flow_npz);names=list(z['channels']);F=z['X'];som=np.array(b.levels('FlowSOM_J8DW'))[b.labels('FlowSOM_J8DW')];rows=list(man.get('audit',{}).get('cellharmony_assignment',[]))
    for space in ['cite_marrow_ADT195','cite_marrow_ADT112','cite_grimes','cite_chinese']:
     if spaces and space not in spaces:continue
     features=b.features(space);adt=[f for f in features if not f.startswith('RNA:')];pairs,_,_=build_crosswalk(names,adt);pairs=pairs.drop_duplicates('cite').drop_duplicates('flow').reset_index(drop=True)
     C=b.feature_matrix(space)[:,[features.index(c) for c in pairs.cite]];Fx=F[:,[names.index(f) for f in pairs.flow]]
     anns=['Population'] if space.startswith('cite_marrow') else ['Mm-MarrowAtlas-L4','pruned', 'StJude' if space=='cite_grimes' else 'Author_celltype']
     for annotation in anns:
      if annotations and annotation not in annotations:continue
      labs=np.array(b.levels(annotation,space))[b.labels(annotation,space)];keep=labs!='unassigned';t=time.time();pred,details=kde_cellharmony(C[keep],Fx,labs[keep],return_details=True)
      prefix=b.space(space)['transfer_links'][annotation];key=prefix+'kde_cellharmony';cat=pd.Categorical(pred);file=key+'.i16';cat.codes.astype(np.int16).tofile(arr/file)
      man['spaces']['flow']['labels'][key]=dict(file=file,levels=list(map(str,cat.categories)),source='transferred',source_space=space,source_annotation=annotation,method='KDE + five-stage community/correlation cellHarmony',community_backend=details['community_backend'],assignment_crosswalk=pairs.to_dict('records'))
      table,summary=concordance(C[keep],labs[keep],Fx,pred,pairs);summary.update(score_transfer(pred,som));summary.update(dataset=space,annotation=annotation,method='kde_cellharmony',seconds=time.time()-t,n_reference_communities=details['reference_communities'],n_flow_communities=details['query_communities'],validation_kind='Same-marker consistency, not independent validation');rows=[r for r in rows if (r['dataset'],r['annotation'])!=(space,annotation)];rows.append(summary);print(summary,flush=True)
      table.to_csv(root/'validation_20261002'/('consistency_'+key+'.tsv'),sep='\t',index=False)
      np.savez(root/'validation_20261002'/(key+'_correlations.npz'),nearest_reference_rho=details['reference_cell_rho'],final_centroid_rho=details['final_rho'])
    man['audit']['cellharmony_assignment']=rows;pd.DataFrame(rows).to_csv(root/'validation_20261002/cellharmony_assignment_summary.tsv',sep='\t',index=False)
    tmp=root/'manifest.new.json';tmp.write_text(json.dumps(man,indent=1));os.replace(tmp,root/'manifest.json')


def main():
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',required=True);p.add_argument('--flow-npz',required=True);p.add_argument('--space',action='append');p.add_argument('--annotation',action='append');a=p.parse_args();build(a.bundle,a.flow_npz,a.space,a.annotation)
if __name__=='__main__':main()
