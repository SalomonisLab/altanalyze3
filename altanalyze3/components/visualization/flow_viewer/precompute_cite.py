"""Build a CITE-only viewer bundle, or add ADT cells to an existing flow bundle.

Input must be an ADT-feature AnnData (not an RNA matrix). Preserve its measured
normalization and cell-aligned embeddings; no inferred antibodies or expression
transformations are introduced here.
"""
import argparse
import json
from pathlib import Path
import re
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def export_cite(adata, bundle, *, name='cite', labels, scale, layer=None, overwrite=False):
    if not re.fullmatch(r'[A-Za-z0-9_-]+',name):raise ValueError('Space name must contain only letters, digits, hyphens and underscores')
    if not labels or not set(labels).issubset(adata.obs):raise ValueError('Select existing observation label columns')
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:raise ValueError('Cell and ADT identifiers must be unique')
    if adata.n_vars<2 or adata.n_obs<1:raise ValueError('Expected a nonempty cell × ADT matrix with at least two ADTs')
    if any(str(f).startswith('RNA:') for f in adata.var_names):raise ValueError('Supply an ADT-only AnnData; RNA features cannot be used for sort gates')
    root=Path(bundle);manifest=root/'manifest.json';root.mkdir(parents=True,exist_ok=True);arr=root/'arrays';arr.mkdir(exist_ok=True)
    man=json.loads(manifest.read_text()) if manifest.exists() else dict(spaces={},gatesets={})
    if 'spaces' not in man:
        from .data_api import FlowBundle
        man=FlowBundle(root).man
    if name in man['spaces'] and not overwrite:raise ValueError('Space already exists; choose a new name or explicitly overwrite')
    prefix='cite_'+name;feature_file=prefix+'_features.f32';temp=arr/(feature_file+'.tmp')
    matrix=adata.layers[layer] if layer else adata.X
    try:
        with temp.open('wb') as handle:
            for start in range(0,adata.n_obs,4096):
                values=matrix[start:start+4096];values=values.toarray() if sparse.issparse(values) else np.asarray(values)
                values=np.asarray(values,np.float32)
                if not np.isfinite(values).all():raise ValueError('ADT matrix contains nonfinite values')
                values.tofile(handle)
        temp.replace(arr/feature_file)
    finally:
        temp.unlink(missing_ok=True)
    spec=dict(n=adata.n_obs,features=list(map(str,adata.var_names)),features_file=feature_file,
              labels={},embeddings={},modality='CITE-seq ADT',feature_scale=scale,default_label=labels[0],
              provenance=dict(source='ADT AnnData',layer=layer or 'X',normalization=scale),cell_ids_file=prefix+'_cell_ids.txt')
    (arr/spec['cell_ids_file']).write_text('\n'.join(map(str,adata.obs_names))+'\n')
    for i,key in enumerate(labels):
        categorical=pd.Categorical(adata.obs[key].astype(str))
        if len(categorical.categories)>32767:raise ValueError('Too many annotation levels for the viewer')
        filename=prefix+'_labels_'+str(i)+'.i16';categorical.codes.astype(np.int16).tofile(arr/filename)
        spec['labels'][key]=dict(file=filename,levels=list(map(str,categorical.categories)))
    for i,(key,values) in enumerate(adata.obsm.items()):
        values=np.asarray(values)
        if values.shape!=(adata.n_obs,2) or not np.isfinite(values).all():continue
        filename=prefix+'_embedding_'+str(i)+'.f32';values.astype(np.float32).tofile(arr/filename)
        spec['embeddings'][key]=dict(file=filename,provenance='Existing cell-aligned AnnData coordinates; feature selection depends on source analysis')
    if spec['embeddings']:spec['default_embedding']=next(iter(spec['embeddings']))
    man['spaces'][name]=spec
    temporary=manifest.with_suffix('.json.tmp');temporary.write_text(json.dumps(man,indent=2));temporary.replace(manifest)
    return spec


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--adt-h5ad',required=True);parser.add_argument('--bundle',required=True);parser.add_argument('--name',default='cite');parser.add_argument('--labels',nargs='+',required=True);parser.add_argument('--scale',required=True,help='Measured ADT normalization/units, e.g. TotalVI denoised or DSB');parser.add_argument('--layer');parser.add_argument('--overwrite',action='store_true');args=parser.parse_args()
    data=ad.read_h5ad(args.adt_h5ad,backed='r')
    try:export_cite(data,args.bundle,name=args.name,labels=args.labels,scale=args.scale,layer=args.layer,overwrite=args.overwrite)
    finally:data.file.close()


if __name__=='__main__':main()
