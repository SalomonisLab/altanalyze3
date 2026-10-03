"""Add measured marrow RNA, a native surface-marker flow UMAP and labelled flow projections."""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
import sys,json
from pathlib import Path
import anndata as ad,numpy as np,pandas as pd
from .precompute_reference import aligned,canonical,dense,add_embedding
from ..scalable_viewer.precompute import compute_embedding
from ..approximate_umap import approximate_umap

def enhance(config,bundle):
    root=Path(bundle);arr=root/'arrays';m=json.loads((root/'manifest.json').read_text())
    # Add actual marrow RNA features to reference cells, joined by full library-qualified IDs.
    ref=m['spaces']['marrow_reference'];ids=pd.Index((arr/ref['cell_ids_file']).read_text().splitlines())
    a=ad.read_h5ad(config['reference']['expression'],backed='r')
    ix=canonical(a.obs_names).get_indexer(canonical(ids));present=ix>=0
    oldref=ref;ref=json.loads(json.dumps(oldref));ref['n']=int(present.sum());ref['description']='Bone marrow reference with measured RNA expression'
    ref['cell_ids_file']='marrow_reference_RNA_cell_ids.txt';ids=ids[present];ix=ix[present]
    (arr/ref['cell_ids_file']).write_text('\n'.join(ids)+'\n')
    for en,spec in ref['embeddings'].items():
     xy=np.fromfile(arr/spec['file'],np.float32).reshape(-1,2)[present];spec['file']='marrow_RNA_'+en+'.f32';xy.tofile(arr/spec['file'])
    for ln,spec in ref['labels'].items():
     codes=np.fromfile(arr/spec['file'],np.int16)[present];spec['file']='marrow_RNA_'+ln+'.i16';codes.tofile(arr/spec['file'])
    m['spaces']['marrow_reference_RNA']=ref;m['audit']['marrow_expression_alignment']={'n_coordinates':len(present),'n_expression_matched':int(present.sum()),'n_expression_missing':int((~present).sum())}
    for k in ['marrow_reference','marrow_reference_RNA']:m['spaces'][k]['transfer_links']={'Population':'transfer_grimes_Ferchen_L4_'}
    genes=[f[4:] for f in m['spaces']['cite_grimes']['features'] if f.startswith('RNA:')]
    cols=a.var_names.get_indexer(genes);out=np.empty((len(ids),len(cols)),np.float32)
    order=np.argsort(ix)
    for start in range(0,len(ids),4096):
     rows=order[start:start+4096];out[rows]=dense(a[ix[rows],:].X[:,cols])
    out.tofile(arr/'marrow_reference_RNA.f32');ref['features_file']='marrow_reference_RNA.f32';ref['features']=['RNA:'+g for g in genes];ref['feature_scale']='Source CP10k/log1p RNA; measured reference expression';ref['provenance']['expression']=config['reference']['expression'];a.file.close()
    print('Added measured marrow RNA genes',len(genes),flush=True)
    # The original Grimes RNA marker UMAP is joined to the matching scTriangulate subset.
    d=next(d for d in config['datasets'] if d['name']=='cite_grimes');df=pd.read_csv(d['marker_coordinates'],sep='\t',index_col=0)
    s=ad.read_h5ad(d['sct'],backed='r')
    xy=df.iloc[aligned(df.index,s.obs_names)][['UMAP_1','UMAP_2']].to_numpy(np.float32)
    add_embedding(m,arr,'cite_grimes_scTriangulate','RNA_MarkerFinder_UMAP',xy,dict(source=d['marker_coordinates'],method='Existing ICGS RNA MarkerFinder cell UMAP; annotation-informed feature selection'))
    s.file.close()
    # A native flow UMAP uses measured surface markers only, excluding scatter, viability, reporters.
    flow=m['spaces']['flow'];features=flow['features'];X=np.memmap(arr/flow['features_file'],np.float32,'r',shape=(flow['n'],len(features)))
    keep=[i for i,c in enumerate(features) if not c.startswith(('SSC','FSC')) and c not in ('DEAD','tdTomato-A')]
    existing=arr/'flow_Flow_measured_surface_UMAP.f32'
    if config.get('reuse_measured_flow_embedding',False) and existing.exists() and existing.stat().st_size==flow['n']*8:
        e=np.fromfile(existing,np.float32).reshape(-1,2)
    else:
        e,pcs,warnings=compute_embedding(np.array(X[:,keep]),min(15,len(keep)-1),30,.3,42)
    add_embedding(m,arr,'flow','Flow_measured_surface_UMAP',e,dict(method='Measured flow surface channels: z-score, PCA, UMAP',independent_of_labels=True,seed=42,channels=[features[i] for i in keep]))
    flow['default_embedding']='Flow_measured_surface_UMAP'
    # Show these same flow events within marrow geometry using transferred marrow annotations.
    refcoords=np.fromfile(arr/ref['embeddings']['Marrow_RNA_reference_UMAP']['file'],np.float32).reshape(-1,2)
    reflevels=ref['labels']['Population']['levels'];refcodes=np.fromfile(arr/ref['labels']['Population']['file'],np.int16)
    reference=ad.AnnData(np.zeros((len(refcodes),1)),obs=pd.DataFrame({'Population':np.array(reflevels)[refcodes]},index=ids));reference.obsm['X_umap']=refcoords
    for method in ['kde_knn','harmony_knn','xgboost']:
     label='transfer_grimes_Ferchen_L4_'+method;spec=flow['labels'][label];codes=np.fromfile(arr/spec['file'],np.int16);lab=np.array(spec['levels'])[codes]
     query=ad.AnnData(np.zeros((len(codes),1)),obs=pd.DataFrame({'Population':lab},index=[str(i) for i in range(len(lab))]))
     q=approximate_umap(query,reference,query_cluster_key='Population',random_state=42,jitter=.05)
     add_embedding(m,arr,'flow','Marrow_projection_Grimes_'+method,np.asarray(q.query_adata.obsm['X_umap']),dict(method='Flow events placed by CITE-seq → flow transferred marrow labels',source_label=label,independent_of_labels=False,limitation='Label-conditioned reference placement; not measured RNA in flow events.'))
    (root/'manifest.json').write_text(json.dumps(m,indent=1));print('Views complete',flush=True)


def main():
    import argparse
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',required=True);ap.add_argument('--bundle',required=True)
    a=ap.parse_args();enhance(json.load(open(a.config)),a.bundle)
if __name__=='__main__':main()
