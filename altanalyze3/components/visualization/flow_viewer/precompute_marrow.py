"""Actual bone-marrow CITE-seq panels and CITE-seq → flow mappings, excluding Thymus.

Panel-specific TotalVI scales remain separate. Original reference RNA coordinates are
joined by cell ID in a subset space; approximate placements are labelled as conditioned
on annotations. Each panel independently transfers marrow labels into the flow events.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ.setdefault(v,'1')
import argparse,json,time
from pathlib import Path
import anndata as ad,numpy as np,pandas as pd
from scipy import sparse
from .precompute_reference import aligned,canonical,dense,add_embedding
from .precompute_spaces import add_space
from .data_api import FlowBundle
from ..approximate_umap import approximate_umap
from ...rna2flow.transfer import kde_knn,minmax_knn,zscore_knn
from ...rna2flow.crosswalk import build_crosswalk
from ...rna2flow.marker_concordance import concordance
from ...rna2flow.evaluate import score_transfer

GENES=['Kit','Ly6a','Itgax','Itgam','Cd27','Il2ra','Cd4','Cd8a','Spi1','Irf1','Irf2','Irf3','Irf4','Irf7','Irf8','Irf9','Flt3','Il7r','Rag1','Rag2','Tcf7','Bcl11b','Gata2','Lmo2','Hlf','Meis1','Dntt','Clec12a','Cd55']


def gene_matrix(source,ids,genes):
    ix=aligned(source.obs_names,ids);cols=source.var_names.get_indexer(genes)
    if (cols<0).any():raise ValueError('RNA genes absent from source')
    order=np.argsort(ix);out=np.empty((len(ids),len(genes)),np.float32)
    for start in range(0,len(ids),2048):
        rows=order[start:start+2048];out[rows]=dense(source[ix[rows],:].X[:,cols])
    return out,ix


def build(config,bundle,flow_npz):
    root=Path(bundle);arr=root/'arrays';b=FlowBundle(root);man=b.man;RNA=ad.read_h5ad(config['reference']['expression'],backed='r')
    ref=man['spaces']['marrow_reference'];refids=pd.Index((arr/ref['cell_ids_file']).read_text().splitlines())
    refxy=b.embedding('Marrow_RNA_reference_UMAP','marrow_reference');reflabs=np.array(b.levels('Population','marrow_reference'))[b.labels('Population','marrow_reference')]
    refadata=ad.AnnData(sparse.csr_matrix((len(refids),1)),obs=pd.DataFrame({'Population':reflabs},index=refids));refadata.obsm['X_umap']=np.asarray(refxy)
    z=np.load(flow_npz);flow_names=list(z['channels']);F=z['X'];som=np.array(b.levels('FlowSOM_J8DW'))[b.labels('FlowSOM_J8DW')];summaries=[]
    # Add measured RNA, including PU.1/IRFs, to the existing reference/thymus cell views.
    for name in ['marrow_reference_RNA','cite_grimes','cite_chinese']:
        sp=man['spaces'][name];ids=pd.Index((arr/sp['cell_ids_file']).read_text().splitlines())
        source=RNA if name!='cite_chinese' else ad.read_h5ad(sp['provenance']['rna'],backed='r')
        G,ix=gene_matrix(source,ids,GENES)
        old=np.asarray(b.feature_matrix(name));adt=[f for f in sp['features'] if not f.startswith('RNA:')];cols=[sp['features'].index(f) for f in adt];X=np.column_stack([old[:,cols],G]);features=adt+['RNA:'+g for g in GENES]
        fn=name+'_PU1_IRF_features.f32';X.tofile(arr/fn);sp['features']=features;sp['features_file']=fn
        if name=='cite_chinese':source.file.close()
        subname=name+'_scTriangulate'
        if subname in man['spaces']:
            # The subset label arrays were built in sctriangulate.obs_names order.
            sct=ad.read_h5ad(sp['provenance']['sct'],backed='r');si=aligned(ids,sct.obs_names)
            sub=man['spaces'][subname];fn=subname+'_PU1_IRF_features.f32';X[si].tofile(arr/fn);sub['features']=features;sub['features_file']=fn;sct.file.close()
    for panel in config['marrow_panels']:
        name=panel['name'];A=ad.read_h5ad(panel['adt']);keep=A.obs.Library.astype(str)!='Thymus';A=A[keep].copy();ids=A.obs_names
        G,ri=gene_matrix(RNA,ids,GENES);labels={'Population':RNA.obs['Mm-MarrowAtlas-L4'].astype(str).values[ri],
            'Library':A.obs.Library.astype(str).values,'Tissue':np.full(A.n_obs,'Bone Marrow')}
        for col in ['groups','group']:
            if col in A.obs:labels['Sort_'+col]=A.obs[col].astype(str).values
        X=np.column_stack([dense(A.X),G]);features=list(A.var_names)+['RNA:'+g for g in GENES]
        add_space(man,str(arr),name,X,features,{},labels,name);sp=man['spaces'][name]
        sp.update(description=panel['description'],modality='Bone marrow CITE-seq: measured TotalVI ADT and RNA',default_label='Population',
            feature_scale='TotalVI batch-corrected non-log ADT; RNA source log-normalized',provenance={'adt':panel['adt'],'rna':config['reference']['expression']},
            transfer_links={'Population':'transfer_'+name+'_'},cell_ids_file=name+'_cell_ids.txt')
        (arr/sp['cell_ids_file']).write_text('\n'.join(ids)+'\n')
        present=canonical(ids).isin(canonical(refids));idx=np.flatnonzero(present);coordix=aligned(refids,ids[idx])
        if len(idx):
            subname=name+'_reference_coordinates';add_space(man,str(arr),subname,X[idx],features,{'Original_marrow_RNA_UMAP':refxy[coordix]}, {k:v[idx] for k,v in labels.items()},subname)
            man['spaces'][subname].update(description=panel['description']+'; original RNA-coordinate subset',modality=sp['modality'],default_label='Population',
                default_embedding='Original_marrow_RNA_UMAP',transfer_links=sp['transfer_links'],feature_scale=sp['feature_scale'])
        # A full-panel label-conditioned placement is available only for labels in the atlas.
        match=np.isin(labels['Population'],reflabs)
        if match.all():
            query=ad.AnnData(sparse.csr_matrix((A.n_obs,1)),obs=pd.DataFrame({'Population':labels['Population']},index=ids))
            q=approximate_umap(query,refadata,query_cluster_key='Population',random_state=42,jitter=.05)
            add_embedding(man,arr,name,'Marrow_approximate_UMAP',np.asarray(q.query_adata.obsm['X_umap']),dict(method='Label-conditioned approximate marrow placement',independent_of_labels=False))
            sp['default_embedding']='Marrow_approximate_UMAP'
        pairs,_,_=build_crosswalk(flow_names,A.var_names);pairs=pairs.drop_duplicates('cite').drop_duplicates('flow').reset_index(drop=True)
        C=dense(A[:,list(pairs.cite)].X);Fx=F[:,[flow_names.index(f) for f in pairs.flow]];L=labels['Population']
        methods={'kde_knn':lambda c,f,l:kde_knn(c,f,l,k=15),'kde_k5_average':lambda c,f,l:kde_knn(c,f,l,k=5,ties='average'),
                 'minmax_knn':minmax_knn,'zscore_knn':zscore_knn}
        for method,fn in methods.items():
            t=time.time();pred=fn(C,Fx,L);key='transfer_'+name+'_'+method;cat=pd.Categorical(pred);file=key+'.i16';cat.codes.astype(np.int16).tofile(arr/file)
            man['spaces']['flow']['labels'][key]=dict(file=file,levels=list(map(str,cat.categories)),source='transferred',source_space=name,source_annotation='Population',source_tissue='Bone Marrow',n_source_cells=A.n_obs)
            table,summary=concordance(C,L,Fx,pred,pairs);summary.update(score_transfer(pred,som));summary.update(source=name,method=method,n_source_cells=A.n_obs,n_source_labels=len(set(L)),n_shared_markers=len(pairs),n_reference_coordinates=int(present.sum()),seconds=time.time()-t,validation_kind='Same-marker consistency, not independent validation')
            summaries.append(summary);table.to_csv(root/'validation_20261002'/f'consistency_{key}.tsv',sep='\t',index=False);print(summary,flush=True)
            # Flow projection uses the transferred marrow labels, not measured flow RNA.
            if set(pred)<=set(reflabs) and method=='kde_k5_average':
                query=ad.AnnData(sparse.csr_matrix((len(pred),1)),obs=pd.DataFrame({'Population':pred},index=[str(i) for i in range(len(pred))]))
                q=approximate_umap(query,refadata,query_cluster_key='Population',random_state=42,jitter=.05)
                add_embedding(man,arr,'flow',name+'_projected_flow',np.asarray(q.query_adata.obsm['X_umap']),dict(method='Flow events placed by bone-marrow CITE-seq → flow labels',source_label=key,independent_of_labels=False))
        man.setdefault('audit',{}).setdefault('marrow_CITE_alignment',{})[name]=dict(n_cells=A.n_obs,n_reference_coordinates=int(present.sum()),n_without_original_coordinates=int((~present).sum()),excluded_libraries=['Thymus'],unknown_reference_states=sorted(set(L)-set(reflabs)))
    RNA.file.close();man['audit']['marrow_CITE_transfer']=summaries
    pd.DataFrame(summaries).to_csv(root/'validation_20261002'/'marrow_CITE_transfer_summary.tsv',sep='\t',index=False)
    tmp=root/'manifest.new.json';tmp.write_text(json.dumps(man,indent=1));os.replace(tmp,root/'manifest.json')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);p.add_argument('--bundle',required=True);p.add_argument('--flow-npz',required=True);a=p.parse_args();build(json.load(open(a.config)),a.bundle,a.flow_npz)
if __name__=='__main__':main()
