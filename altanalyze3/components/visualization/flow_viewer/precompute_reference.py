"""Add cell-ID-aligned RNA, ADT and marrow-reference views to an existing bundle.

Configuration specifies actual source paths. RNA HVGs use no annotation labels; the
approximate marrow placement calls the existing approximate_umap implementation and
is explicitly labelled as annotation-conditioned, rather than independent evidence.
"""
from __future__ import annotations
import argparse, json, os
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from .precompute_spaces import add_space


def canonical(ids):
    # Match cellHarmony barcodes to the vendor convention without dropping library IDs.
    return pd.Index([str(s).rsplit('.',1)[1]+'_'+str(s).rsplit('.',1)[0]
                     if '.' in str(s) else str(s) for s in ids])


def aligned(source_ids, target_ids):
    s,t=canonical(source_ids),canonical(target_ids)
    if not s.is_unique or not t.is_unique:raise ValueError('Nonunique canonical cell IDs')
    ix=s.get_indexer(t)
    if (ix<0).any():raise ValueError(f'{int((ix<0).sum())} unmatched cells')
    return ix


def dense(x):return np.asarray(x.toarray() if sparse.issparse(x) else x,dtype=np.float32)


def add_embedding(man,arr,space,name,xy,metadata):
    if xy.shape!=(man['spaces'][space]['n'],2) or not np.isfinite(xy).all():
        raise ValueError('Embedding must contain finite coordinates for every cell')
    fn=f'{space}_{name}.f32';np.asarray(xy,np.float32).tofile(arr/fn)
    man['spaces'][space]['embeddings'][name]={**{k:v for k,v in metadata.items() if k != 'file'},'file':fn}


def coords(path,ids):
    df=pd.read_csv(path,sep='\t',index_col=0)
    columns=next((c for c in [['UMAP-X','UMAP-Y'],['UMAP_1','UMAP_2'],['UMAP1','UMAP2']] if set(c)<=set(df)),None)
    if columns is None:raise ValueError(f'No coordinate columns in {path}')
    return df.iloc[aligned(df.index,ids)][columns].to_numpy(np.float32)


def build(config,bundle):
    root=Path(bundle);arr=root/'arrays';man=json.loads((root/'manifest.json').read_text())
    refcfg=config['reference']
    rc=pd.read_csv(refcfg['coordinates'],sep='\t',index_col=0)
    n_reference_missing=int((~np.isfinite(rc[['UMAP1','UMAP2']].to_numpy()).all(1)).sum())
    rc=rc[np.isfinite(rc[['UMAP1','UMAP2']].to_numpy()).all(1)]
    rl=pd.read_csv(refcfg['labels'],sep='\t',index_col=0)
    rl=rl.iloc[aligned(rl.index,rc.index)]
    ref=ad.AnnData(sparse.csr_matrix((len(rc),1)),obs=pd.DataFrame({'Population':rl['Population'].astype(str).values},index=rc.index))
    ref.obsm['X_umap']=rc[['UMAP1','UMAP2']].to_numpy(np.float32)
    add_space(man,str(arr),'marrow_reference',np.empty((len(rc),0)),[],
              {'Marrow_RNA_reference_UMAP':ref.obsm['X_umap']},{'Population':ref.obs.Population},'marrow_reference')
    man['spaces']['marrow_reference'].update(description='Bone marrow RNA reference cells',modality='RNA reference',
        default_embedding='Marrow_RNA_reference_UMAP',default_label='Population',provenance=refcfg)
    (arr/'marrow_reference_cell_ids.txt').write_text('\n'.join(map(str,rc.index))+'\n')
    man['spaces']['marrow_reference']['cell_ids_file']='marrow_reference_cell_ids.txt'
    from ..approximate_umap import approximate_umap
    from ..scalable_viewer.precompute import compute_embedding
    import scanpy as sc
    audit={'reference':{'excluded_missing_coordinates':n_reference_missing}}
    for d in config['datasets']:
        name=d['name'];print('Building',name,flush=True)
        rna=ad.read_h5ad(d['rna'],backed='r');adt=ad.read_h5ad(d['adt'],backed='r');sct=ad.read_h5ad(d['sct'],backed='r')
        rids=canonical(rna.obs_names);aids=canonical(adt.obs_names)
        keep=np.flatnonzero(aids.isin(rids));ids=adt.obs_names[keep];ri=aligned(rna.obs_names,ids)
        # backed AnnData fancy indexes need sorted source row indices.
        sort=np.argsort(ri);inverse=np.argsort(sort)
        R=rna[ri[sort],:].to_memory()[inverse].copy()
        A=dense(adt[keep,:].X)
        labels={}
        for c in ['Mm-MarrowAtlas-L4','Author_celltype','Library','group']:
            if c in R.obs:labels[c]=R.obs[c].astype(str).values
        for c in ['celltype','library']:
            if c in adt.obs:labels[c]=adt.obs[c].iloc[keep].astype(str).values
        if 'Author_celltype' not in labels and 'celltype' in labels:
            labels['Author_celltype']=labels['celltype']
        si=canonical(sct.obs_names).get_indexer(canonical(ids));present=si>=0
        for c in ['pruned','cluster_name','hopach_cluster','ICGS3_RNA','Leiden_RNA_r2','ICGS3_ADT','Ferchen_L4','StJude']:
            if c in sct.obs:
                v=np.full(len(ids),'unassigned',object);v[present]=sct.obs[c].astype(str).values[si[present]];labels[c]=v
        labels['Tissue']=np.full(len(ids),'Thymus')
        # Source Condition is known to mislabel the Grimes Thymus library as marrow.
        if 'Condition' in adt.obs:labels['Source_Condition_unverified']=adt.obs.Condition.iloc[keep].astype(str).values
        add_space(man,str(arr),name,A,list(adt.var_names),{},labels,name)
        space=man['spaces'][name]
        space.update(description=d['description'],modality='CITE-seq (measured RNA and ADT)',
                     provenance={k:d[k] for k in ['rna','adt','sct']},feature_scale=d['feature_scale'],
                     transfer_links=d['transfer_links'],default_label=d['default_label'],
                     default_embedding='RNA_HVG_UMAP',cell_ids_file=f'{name}_cell_ids.txt')
        (arr/f'{name}_cell_ids.txt').write_text('\n'.join(map(str,ids))+'\n')
        # HVG selection on already-normalized RNA: no labels enter PCA/UMAP.
        sc.pp.highly_variable_genes(R,n_top_genes=min(2000,R.n_vars),flavor='seurat',inplace=True)
        hvg=np.flatnonzero(R.var.highly_variable.values)
        emb,pcs,warnings=compute_embedding(dense(R.X[:,hvg]),30,30,.3,42)
        add_embedding(man,arr,name,'RNA_HVG_UMAP',emb,dict(method='RNA HVG z-score, randomized PCA, UMAP',
            n_features=len(hvg),seed=42,source=d['rna'],independent_of_labels=True))
        (root/f'{name}_RNA_HVG_genes.txt').write_text('\n'.join(map(str,R.var_names[hvg]))+'\n')
        for e in d.get('coordinate_files',[]):
            add_embedding(man,arr,name,e['name'],coords(e['path'],ids),dict(source=e['path'],method=e['method']))
        # Existing scTriangulate UMAP covers a subset. NaNs preserve missingness instead
        # of inventing coordinates; serve it in a separate, matching subset space.
        sub_ids=sct.obs_names
        subix=aligned(ids,sub_ids)
        labs={k:np.asarray(v)[subix] for k,v in labels.items()}
        subname=name+'_scTriangulate'
        add_space(man,str(arr),subname,A[subix],list(adt.var_names),
                  {'scTriangulate_original_UMAP':np.asarray(sct.obsm['X_umap'],np.float32)},labs,subname)
        man['spaces'][subname].update(description=d['description']+'; scTriangulate subset',
              modality='CITE-seq',transfer_links=d['transfer_links'],feature_scale=d['feature_scale'],
              default_embedding='scTriangulate_original_UMAP',default_label='pruned')
        # RNA colour features are measured expression, distinguishable from ADT names.
        genes=[g for g in ['Kit','Gata2','Flt3','Il7r','Rag1','Rag2','Tcf7','Bcl11b','Lmo2','Hlf','Meis1','Cd34','Spi1','Gata1','Cebpa','Cd3d','Cd4','Cd8a','Ly6a','Dntt','Clec12a','Cd55'] if g in R.var_names]
        if not genes:raise ValueError('RNA gene symbols unavailable; configure a gene-name mapping')
        G=dense(R[:,genes].X)
        F=np.column_stack([A,G]);features=list(adt.var_names)+['RNA:'+g for g in genes]
        F.tofile(arr/space['features_file']);space['features']=features
        sub=man['spaces'][subname];F[subix].tofile(arr/sub['features_file']);sub['features']=features
        if set(labels['Mm-MarrowAtlas-L4'])-set(ref.obs.Population):raise ValueError('Marrow labels absent from reference')
        query=ad.AnnData(sparse.csr_matrix((len(ids),1)),obs=pd.DataFrame({'Population':labels['Mm-MarrowAtlas-L4']},index=ids))
        result=approximate_umap(query,ref,query_cluster_key='Population',reference_cluster_key='Population',random_state=42,jitter=.05)
        xy=np.asarray(result.query_adata.obsm['X_umap'],np.float32)
        meta=dict(method='approximate_umap: label-conditioned sampled reference coordinates + jitter',
                  independent_of_labels=False,seed=42,source=refcfg['coordinates'],
                  limitation='Placement follows assigned marrow labels; separation is not independent validation.')
        add_embedding(man,arr,name,'Marrow_approximate_UMAP',xy,meta)
        add_embedding(man,arr,subname,'Marrow_approximate_UMAP',xy[subix],meta)
        add_embedding(man,arr,subname,'RNA_HVG_UMAP',emb[subix],space['embeddings']['RNA_HVG_UMAP'])
        audit[name]=dict(adt_cells=adt.n_obs,rna_cells=rna.n_obs,joined_cells=len(ids),
                         unmatched_adt_cells=int(adt.n_obs-len(ids)),sct_cells=int(present.sum()),
                         sct_missing=int((~present).sum()),n_hvg=len(hvg),RNA_genes=genes,
                         source_condition_warning='Grimes source Condition says Bone Marrow for Thymus library' if 'Condition' in adt.obs else None)
        rna.file.close();adt.file.close();sct.file.close()
    man.setdefault('audit',{})['cell_alignment']=audit
    man['audit']['scientific_interpretation']='Labels transferred CITE-seq → flow. Gates approximate these labels; do not establish stem-cell function. Approximate reference placement uses labels. Event sampling is not biological replication.'
    tmp=root/'manifest.new.json';tmp.write_text(json.dumps(man,indent=1));os.replace(tmp,root/'manifest.json')
    return audit


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',required=True);ap.add_argument('--bundle',required=True)
    a=ap.parse_args();print(json.dumps(build(json.load(open(a.config)),a.bundle),indent=2))

if __name__=='__main__':main()
