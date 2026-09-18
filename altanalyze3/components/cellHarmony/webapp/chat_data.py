"""Sample-aware adapters and bounded caches for the shared Chat executors."""
from collections import OrderedDict
from functools import cached_property
from pathlib import Path
from threading import RLock
import json
import re
import numpy as np
import pandas as pd
from scipy import sparse
from .integration_data import UploadedIntegrationData

_DATA_LOCK=RLock()


class ProtocolData:
    def _init_cache(self):
        self._aggregates=OrderedDict();self._lock=RLock();self.aggregate_builds=0

    @cached_property
    def states(self):return list(dict.fromkeys(self.obs[self.cluster_key].astype(str)))
    @cached_property
    def state_code(self):return pd.Categorical(self.obs[self.cluster_key].astype(str),categories=self.states).codes
    @cached_property
    def state_n(self):return np.bincount(self.state_code,weights=self.cell_weights,minlength=len(self.states))
    @cached_property
    def _covariates(self):
        result={}
        identity={'donor','Donor','Library','sample','Sample','sample_id','meta_sample'}
        for col in self.obs:
            v=self.obs[col]
            if col not in identity and pd.api.types.is_numeric_dtype(v) and not isinstance(v.dtype,pd.CategoricalDtype):
                result[col]=('numeric',pd.to_numeric(v,errors='coerce').to_numpy(dtype=float),None)
            else:
                cat=pd.Categorical(v)
                result[col]=('categorical',cat.codes,[str(x) for x in cat.categories])
        return result
    def covariate_names(self):
        return {k:{'kind':v[0],'categories':v[2]} for k,v in self._covariates.items()}
    def covariate_values(self,name):return self._covariates[name]
    @cached_property
    def _gene_index(self):return {str(g).upper():i for i,g in enumerate(self.symbols)}
    def resolve_gene(self,g):return self._gene_index.get(str(g).upper())
    def ordered_covariate_values(self,name):
        kind,values,labels=self.covariate_values(name)
        if kind!='categorical':return kind,values,labels
        column=self.obs[name]
        if isinstance(column.dtype,pd.CategoricalDtype) and column.cat.ordered:
            return kind,values,labels
        roman={'i':1,'ii':2,'iii':3,'iv':4,'v':5}
        scores={}
        for i,label in enumerate(labels):
            text=str(label).lower().strip()
            if text.startswith(('not applicable','unknown','not recorded')):continue
            if text=='normal spirometry':scores[i]=0;continue
            match=re.fullmatch(r'(?:gold|stage|grade)?\s*(\d+(?:\.\d+)?|iv|iii|ii|i|v)',text)
            if not match:raise ValueError(f'{name} needs explicit ordered categories before a stage trend can be tested.')
            token=match.group(1);scores[i]=roman[token] if token in roman else float(token)
        order=sorted(scores,key=lambda i:scores[i]);mapping={old:i for i,old in enumerate(order)}
        return kind,np.array([mapping.get(int(v),-1) for v in values]),[labels[i] for i in order]

    def donor_axis(self):
        col=next((c for c in ('donor','Donor','Sample','sample_id','sample','Library','meta_sample') if c in self.obs),None)
        if col is None:return None,None
        labels=self.obs[col].astype('string')
        study=next((c for c in ('Study_internal','study') if c in self.obs),None)
        if study:labels=self.obs[study].astype('string')+'|'+labels
        cat=pd.Categorical(labels)
        return cat.codes,[str(x) for x in cat.categories]
    def donor_pseudobulk(self,rows,state='',min_cells=5):
        key=(tuple(rows),state,int(min_cells))
        with self._lock:
            if key in self._aggregates:
                self._aggregates.move_to_end(key);return self._aggregates[key]
            value=self._aggregate(rows,state,min_cells);self.aggregate_builds+=1
            self._aggregates[key]=value
            while len(self._aggregates)>12:self._aggregates.popitem(last=False)
            return value


class UploadedChatData(ProtocolData,UploadedIntegrationData):
    def __init__(self,app,meta):
        super().__init__(app,meta);self._init_cache();self.obs=self.rna['adata'].obs
        self.cluster_key=self.rna['cluster_key'];self.symbols=list(map(str,self.rna['var_names']))
        self.cell_weights=np.ones(len(self.obs));self.sv={'cluster_key':self.cluster_key}
        self.provenance='Real biological-sample count sums normalized to log2(1+CP10k); cells are not statistical replicates.'
    def comparison_covariate(self,contrast):
        run=self.runs[contrast];cfg=run['config'];values=self.obs[cfg['sample_field']].astype(str)
        case,control=run.get('case_label','Group 1'),run.get('control_label','Group 2')
        codes=np.where(values.isin(cfg['group1_samples']),0,np.where(values.isin(cfg['group2_samples']),1,-1))
        key='_chat_comparison_groups'
        self._covariates[key]=('categorical',codes,[case,control])
        return key,case,control

    def deg_table(self,comp_id,max_rows=200000,fdr_max=None,state=None):
        return super().deg_table(comp_id,max_rows=max_rows,fdr_max=fdr_max,state=state)
    def markers(self):
        f=self.web._chat_marker_table(self.meta,'rna')
        if f.empty:return []
        return f.rename(columns={'Gene':'gene'}).to_dict('records')
    def _aggregate(self,rows,state,min_cells):
        codes,labels=self.donor_axis()
        if codes is None:return None,[],[]
        a=self.rna['adata']
        if 'counts' not in a.layers:raise ValueError('Raw RNA counts are required for sample-level Chat analyses.')
        keep=codes>=0
        if state:keep &= self.obs[self.cluster_key].astype(str).to_numpy()==state
        counts=np.bincount(codes[keep],minlength=len(labels));usable=np.flatnonzero(counts>=min_cells)
        if not usable.size:return None,[],[]
        cells=np.flatnonzero(keep);mapping=np.full(len(labels),-1);mapping[usable]=np.arange(len(usable))
        cells=cells[mapping[codes[cells]]>=0]
        grouping=sparse.csr_matrix((np.ones(len(cells)),(mapping[codes[cells]],np.arange(len(cells)))),shape=(len(usable),len(cells)))
        # Sum the entire transcriptome before library-size normalization.
        summed=grouping @ a.layers['counts'][cells,:]
        totals=np.asarray(summed.sum(axis=1)).ravel();selected=summed[:,rows]
        selected=selected.toarray() if sparse.issparse(selected) else np.asarray(selected)
        if np.any(selected<0):raise ValueError('Raw RNA counts cannot be negative.')
        valid=totals>0
        matrix=np.log2(1+10000*selected[valid]/totals[valid,None]).T
        return matrix,[labels[i] for i in usable[valid]],counts[usable[valid]].tolist()


class BundleChatData(ProtocolData):
    def __init__(self,ds):
        self.ds=ds;self._init_cache();self.cluster_key='cell_state';self.sv=ds.sv
        root=Path(ds.sv.get('integrated_root',''));manifest=root/'manifest.json'
        info=json.loads(manifest.read_text()) if manifest.is_file() else {}
        self.has_sample_pseudobulks=bool(info.get('chat_pseudobulk'))
        if not self.has_sample_pseudobulks:
            from altanalyze3.components.visualization.scalable_viewer.bundle_meta import build_obs
            self.obs=build_obs(ds)[0];self.cluster_key=ds.cluster_key
            self.cell_weights=np.ones(len(self.obs));self.symbols=ds.symbols
            self.provenance='Existing precomputed results; sample-level RNA is not available.'
            return
        base=root/info['chat_pseudobulk'];self.obs=pd.read_csv(base/'observations.tsv',sep='\t')
        self.cell_weights=self.obs.n_cells.to_numpy(dtype=float)
        self.matrix=np.load(base/'expression.npy',mmap_mode='r');self.symbols=np.load(base/'genes.npy').tolist()
        self.provenance=info['chat_pseudobulk_scale']+'; one biological sample per cell state.'
    def __getattr__(self,k):return getattr(self.ds,k)
    def donor_axis(self):
        if not self.has_sample_pseudobulks:return None,None
        return super().donor_axis()
    def _aggregate(self,rows,state,min_cells):
        if not self.has_sample_pseudobulks:raise ValueError('Real sample pseudobulks are required; imputed meta-samples are not biological donors.')
        codes,labels=self.donor_axis();keep=(codes>=0)&(self.cell_weights>=min_cells)
        if state:keep &= self.obs.cell_state.astype(str).to_numpy()==state
        positions=np.flatnonzero(keep)
        if not len(positions):return None,[],[]
        if len(set(codes[positions]))!=len(positions):
            raise ValueError('Choose a cell state for sample pseudobulk expression analyses.')
        return np.asarray(self.matrix[positions][:,rows].T,dtype=float),[labels[codes[i]] for i in positions],self.cell_weights[positions].astype(int).tolist()


def protocol_adapter(ds):
    if isinstance(ds,ProtocolData):return ds
    with _DATA_LOCK:
        cached=getattr(ds,'_chat_protocol_adapter',None)
        if cached is None:
            cached=BundleChatData(ds);ds._chat_protocol_adapter=cached
        return cached


def data_for_job(app,meta):
    if hasattr(app.state.job_store,'dataset'):
        return protocol_adapter(app.state.job_store.dataset(meta['job_id']))
    # The expression cache already includes the h5ad generation. Include differential
    # history and file mtimes so recomputation or in-place artifact replacement expires it.
    paths=[(meta.get('artifacts') or {}).get('combined_h5ad','')]
    for run in (meta.get('differential_history') or {}).values():paths.extend((run.get('artifacts') or {}).values())
    stamp=[]
    for p in paths:
        if isinstance(p,str) and Path(p).is_file():
            st=Path(p).stat();stamp.append((p,st.st_mtime_ns,st.st_size))
    signature=(meta['job_id'],meta.get('updated_at'),json.dumps(meta.get('differential') or {},sort_keys=True),tuple(stamp))
    with _DATA_LOCK:
        if not hasattr(app.state,'chat_protocol_data'):app.state.chat_protocol_data=OrderedDict()
        cache=app.state.chat_protocol_data
        if signature not in cache:
            cache[signature]=UploadedChatData(app,meta)
            while len(cache)>4:cache.popitem(last=False)
        cache.move_to_end(signature)
        return cache[signature]
