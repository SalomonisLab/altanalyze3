"""Shared cellHarmony analysis and sparse sampling for RNA-release evaluation."""
from pathlib import Path
import gc
import hashlib
import json
import os
import subprocess
import shutil
import time
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, io
from . import benchmark_nd20_167 as previous
from .. import ambient_subtract as ambient
from altanalyze3.components.cellHarmony import cellHarmony_lite as lite
from altanalyze3.components.cellHarmony import cellHarmony_differential as diff
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/'benchmarking/ND20_167_release_20261003'
REF = previous.REF
KEY = previous.KEY
RSCRIPT = '/Library/Frameworks/R.framework/Resources/bin/Rscript'
METHODS = ['contaminated','python_auto','SoupX_known','SoupX_auto']
def save_table(rows,path):
    pd.DataFrame(rows).to_csv(path,sep='\t',index=False)

def align(a,directory,cutoff=None):
    directory.mkdir(parents=True,exist_ok=True)
    timing=[]; original=ambient.process_anndata
    def timed(*args,**kwargs):
        start=time.perf_counter(); result=original(*args,**kwargs)
        timing.append(time.perf_counter()-start); return result
    ambient.process_anndata=timed
    try:
        _, result=lite.combine_and_align_h5([],str(REF),adata=a,output_dir=str(directory),
            min_genes=0,min_cells=None,min_counts=0,mit_percent=101,generate_umap=False,
            save_adata=False,unsupervised_cluster=False,alignment_mode='cosine',
            min_alignment_score=None,ambient_correct_cutoff=cutoff,
            ambient_memory_efficient=True,return_adata=True)
    finally:
        ambient.process_anndata=original
    result=result[a.obs_names].copy()
    assert 'counts' in result.layers and np.isfinite(result.X.data).all()
    result.raw=None
    for layer in list(result.layers):
        if layer!='counts': del result.layers[layer]
    return result,sum(timing)

def poisson_matrix(means,weights,rng):
    # Independent Poisson gene counts are equivalent to a Poisson total
    # followed by multinomial allocation; construct bounded dense blocks.
    support=np.flatnonzero(weights>0); blocks=[]
    for start in range(0,len(means),128):
        sample=rng.poisson(means[start:start+128,None]*weights[support][None,:])
        local=sparse.coo_matrix(sample)
        blocks.append(sparse.csr_matrix((local.data,(local.row,support[local.col])),
            shape=(sample.shape[0],len(weights)),dtype=np.float32))
    return sparse.vstack(blocks,format='csr')

def truth_metrics(corrected,base,contaminated,method,seconds):
    rows=[]
    for cap in ['HSC','MPP']:
        mask=base.obs.Library.to_numpy()==cap
        clean=base.layers['counts'][mask]; observed=contaminated.X[mask]; output=corrected.layers['counts'][mask]
        difference=(output-clean).tocsr(); added=previous.total(observed)-previous.total(clean)
        excess=float(difference.data[difference.data>0].sum(dtype=np.float64))
        deficit=float(-difference.data[difference.data<0].sum(dtype=np.float64))
        removed=previous.total(observed)-previous.total(output)
        rows.append(dict(method=method,capture=cap,baseline_counts=previous.total(clean),
            observed_counts=previous.total(observed),retained_counts=previous.total(output),
            added_counts=added,excess_counts=excess,baseline_deficit_counts=deficit,
            absolute_error_counts=excess+deficit,error_per_added_count=(excess+deficit)/added,
            recovery_percent=100*(1-(excess+deficit)/added),
            baseline_deficit_percent=100*deficit/previous.total(clean),
            loss_percent=100*removed/previous.total(observed),correction_seconds=seconds))
    return rows

def fingerprint(x):
    h=hashlib.sha256()
    for array in [x.data,x.indices,x.indptr]:h.update(array.tobytes())
    return h.hexdigest()

def evaluate(a,base,method,directory,shared,primary):
    result_file=directory/'metrics.json'
    if result_file.exists():return json.loads(result_file.read_text())
    keep=base.obs.analysis_included.to_numpy()
    b=a[keep].copy();b.obs['fixed_population']=base.obs.loc[keep,'fixed_population'].to_numpy()
    store=diff.run_de_for_comparisons(b[b.obs.fixed_population.isin(shared)].copy(),
        'fixed_population','Library','MPP','HSC',method='wilcoxon',alpha=.05,fc_thresh=1.2,
        min_cells_per_group=20,use_rawp=False)
    tables=[]
    for state,table in store['per_population_deg'].items():
        tables.append(table.reset_index().assign(population=state))
    stats=pd.concat(tables,ignore_index=True)
    stats.to_csv(directory/'all_gene_statistics.tsv.gz',sep='\t',index=False)
    continuous=[];gene_state=set()
    normalized=b.X.copy();normalized.data=np.expm1(normalized.data)
    for pop in shared:
        t=stats[stats.population==pop].copy()
        selected=(t.fdr<.05)&(t.log2fc.abs()>np.log2(1.2))
        eligible=(b.obs.fixed_population.to_numpy()==pop)
        case=eligible&(b.obs.Library.to_numpy()=='MPP');control=eligible&(b.obs.Library.to_numpy()=='HSC')
        cp_effect=np.log2((np.asarray(normalized[case].mean(0)).ravel()+1)/
                         (np.asarray(normalized[control].mean(0)).ravel()+1))
        mapping=pd.Series(cp_effect,index=b.var_names)
        t['CP10k_log2FC']=t.gene.map(mapping)
        normalized_deg=(t.fdr<.05)&(t.CP10k_log2FC.abs()>np.log2(1.2))
        continuous.append(dict(population=pop,DEGs=int(selected.sum()),CP10k_DEGs=int(normalized_deg.sum())))
        if pop in primary:gene_state.update((pop,gene)for gene in t.loc[selected,'gene'])
    save_table(continuous,directory/'DE_by_state.tsv')
    result=dict(method=method,primary_DEG_pairs=len(gene_state),
        primary_CP10k_DEG_pairs=sum(row['CP10k_DEGs'] for row in continuous if row['population'] in primary),
        population_reassignment_percent=100*float(np.mean(b.obs[KEY].astype(str).to_numpy()!=b.obs.fixed_population.to_numpy())))
    base_stats=OUT/'baseline/all_gene_statistics.tsv.gz'
    if method!='baseline' and base_stats.exists():
        t=pd.read_csv(base_stats,sep='\t');t=t[t.population.isin(primary)]
        baseline_set=set(zip(t.loc[(t.fdr<.05)&(t.log2fc.abs()>np.log2(1.2)),'population'],
                             t.loc[(t.fdr<.05)&(t.log2fc.abs()>np.log2(1.2)),'gene']))
        result.update(new_DEG_pairs=len(gene_state-baseline_set),lost_baseline_DEG_pairs=len(baseline_set-gene_state),
                      retained_baseline_DEG_pairs=len(gene_state&baseline_set))
    embeddings=[]
    for seed in (0,1,2):
        path=directory/f'embedding_seed{seed}.npz'
        if path.exists():
            e=np.load(path); xy,pcs=e['umap'],e['pca']
        else:
            xy,pcs,n_genes=previous.embed(b,seed)
            np.savez_compressed(path,umap=xy,pca=pcs,cell=b.obs_names.to_numpy(dtype=str))
        embeddings.extend(previous.distances(b,xy,pcs,method,seed,shared))
        if seed==0:
            pd.DataFrame(dict(cell=b.obs_names,capture=b.obs.Library,population=b.obs.fixed_population,
                              UMAP1=xy[:,0],UMAP2=xy[:,1])).to_csv(directory/'umap_coordinates.tsv',sep='\t',index=False)
    save_table(embeddings,directory/'embedding_distances.tsv')
    pred=pd.DataFrame(embeddings);pred=pred[pred.population.isin(primary)]
    result.update(UMAP_D_R=float(pred.UMAP_distance_within_scaled.mean()),
                  PCA_D_R=float(pred.PCA_distance_within_scaled.mean()),PCA_mixing=float(pred.PCA_neighbor_mixing.mean()))
    result_file.write_text(json.dumps(result,indent=2)+'\n')
    del b,store,normalized;gc.collect()
    return result
