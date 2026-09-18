"""Dataset adapters for Discover-compatible integrated views."""
from functools import lru_cache
from pathlib import Path
import json
import sqlite3
import numpy as np
from scipy import sparse
from altanalyze3.components.cellHarmony.grn_analysis import _sibling_comparison
from .grn_data import UploadedGrnData


class UploadedIntegrationData(UploadedGrnData):
    state_expression_scale = 'mean stored RNA expression in the cell state'

    def marker_rows(self, state):
        frame = self.web._chat_marker_table(self.meta, 'rna')
        if frame.empty:
            return []
        frame = frame.loc[frame['cluster'].astype(str) == state].rename(
            columns={'Gene': 'gene', 'Fold': 'log2fc', 'FDR p-value': 'fdr'})
        return json.loads(frame.to_json(orient='records'))

    def state_expression(self, state, features):
        cache = self.rna
        wanted = set(features)
        cols = [i for i, g in enumerate(cache['var_names']) if g in wanted]
        rows = np.flatnonzero(np.asarray(cache['populations'], dtype=str) == state)
        if not len(rows) or not cols:
            return {}
        means = np.asarray(cache['adata'].X[rows][:, cols].mean(axis=0)).ravel()
        return {str(cache['var_names'][i]): float(v) for i, v in zip(cols, means) if np.isfinite(v)}

    def comparison(self, contrast):
        return contrast or self.current_contrast

    def available(self, contrast, modality, state=''):
        selected=_sibling_comparison(self,self.comparison(contrast),modality)
        if not selected:return False
        # A completed run with zero calls is still an analysis, but another state is not.
        if state:
            snapshot=dict(self.meta,differential=self.runs[selected])
            tested=set(self.web._differential_result_populations(self.app,snapshot))
            tested.update(self.web._differential_heatmap_populations(self.app,snapshot))
            return state in tested
        return True

    def differentials(self, contrast, state, modality):
        selected=_sibling_comparison(self,self.comparison(contrast),modality)
        if not selected:return []
        return self.deg_table(selected,max_rows=1000000,state=state or None)['rows']

    def features(self, modality):
        if modality=='rna':return list(map(str,self.rna['var_names']))
        try:return self.modality(modality).features
        except (KeyError,FileNotFoundError):return []

    def edge_scores(self,state):
        store=self.modality('grn');col=self.states.index(state)
        return dict(zip(store.features,map(float,store.stats_mean[:,col])))

    def marker_features(self,state,modality='rna',limit=50):
        frame=self.web._chat_marker_table(self.meta,modality)
        if frame.empty:return []
        return frame.loc[frame['cluster'].astype(str)==state,'Gene'].astype(str).head(limit).tolist()

    def arm_expression(self,contrast,state,features):
        selected=_sibling_comparison(self,self.comparison(contrast),'rna')
        if not selected:return {},{}
        cfg=self.runs[selected]['config'];a=self.rna['adata'];obs=a.obs
        state_col=cfg['population_col'];group_col=cfg['sample_field']
        sample_col=next((c for c in ('sample','Sample','Library','library','Donor','donor') if c in obs),None)
        if not sample_col or 'counts' not in a.layers:
            raise ValueError('Integrated TF expression requires raw counts and biological sample identifiers for RNA pseudobulks.')
        idx=np.flatnonzero(obs[state_col].astype(str).to_numpy()==state)
        levels=[];names=list(map(str,a.var_names));lookup={g:i for i,g in enumerate(names)}
        for group in ('group1_samples','group2_samples'):
            arm=idx[obs.iloc[idx][group_col].astype(str).isin(cfg[group]).to_numpy()]
            means=[]
            for sample in obs.iloc[arm][sample_col].astype(str).unique():
                rows=arm[obs.iloc[arm][sample_col].astype(str).to_numpy()==sample]
                counts=np.asarray(a.layers['counts'][rows].sum(axis=0)).ravel().astype(float)
                if counts.sum()>0:means.append(np.log2(1+10000*counts/counts.sum()))
            mean=np.mean(means,axis=0) if means else None
            levels.append({g:float(mean[lookup[g]]) for g in features if g in lookup} if mean is not None else {})
        return tuple(levels)


class BundleIntegrationData:
    state_expression_scale = 'mean stored RNA expression in the cell state'

    def marker_rows(self, state):
        return [dict(gene=r['gene'], log2fc=r.get('fold'), fdr=r.get('p'))
                for r in self.ds.markers() if r['cluster'] == state]

    def state_expression(self, state, features):
        col = self.states.index(state)
        wanted = set(features)
        return {str(g): float(self.ds.stats_mean[i, col])
                for i, g in enumerate(self.ds.symbols)
                if g in wanted and np.isfinite(self.ds.stats_mean[i, col])}

    def __init__(self,ds,meta,assets):
        self.ds,self.meta=ds,meta
        self.current_contrast=(meta.get('differential') or {}).get('run_id','')
        self.states=ds.states
        root=assets.get('integrated_root') or ds.sv.get('integrated_root')
        self.root=Path(root) if root else None
        self.manifest=json.loads((self.root/'manifest.json').read_text()) if self.root and (self.root/'manifest.json').is_file() else {'comparisons':{}}

    def comparison(self,contrast):
        key=contrast or self.current_contrast
        if key in self.manifest['comparisons']:return key
        c=next((c for c in self.ds.deg_manifest().get('comparisons',[]) if c['id']==key),{})
        return c.get('comparison','')

    def available(self,contrast,modality,state=''):
        r=self.manifest['comparisons'].get(self.comparison(contrast),{}).get('modalities',{}).get(modality,{})
        if not r or str(r.get('status','')).startswith('not_applicable'):return False
        return not state or state in r.get('states',[])

    def differentials(self,contrast,state,modality):
        if not self.root:return []
        with sqlite3.connect(f"file:{self.root/self.manifest['differential_database']}?mode=ro",uri=True) as con:
            con.row_factory=sqlite3.Row
            query='SELECT * FROM differential WHERE comparison=? AND modality=?'
            params=[self.comparison(contrast),modality]
            if state:query+=' AND population=?';params.append(state)
            return [dict(row) for row in con.execute(query,params)]

    def features(self,modality):
        if modality=='rna':return self.ds.symbols
        try:return list(self.ds.modality(modality).features)
        except (KeyError,FileNotFoundError):return []

    def marker_features(self,state,modality='rna',limit=50):
        if modality!='rna':return []
        return self.ds.marker_genes_for([state],limit)

    def edge_scores(self,state):
        if self.root and self.manifest.get('edge_scores'):
            return _edge_index(str(self.root/self.manifest['edge_scores'])).get(state,{})
        store=self.ds.modality('grn');col=self.states.index(state)
        return dict(zip(store.features,map(float,store.stats_mean[:,col])))

    def arm_expression(self,contrast,state,features):
        r=self.manifest['comparisons'][self.comparison(contrast)]
        path=self.root/r['arm_means'];stamp=path.stat()
        return tuple({g:values[g] for g in features if g in values}
                     for values in _arm_means(str(path),stamp.st_mtime_ns,stamp.st_size,state))


@lru_cache(maxsize=8)
def _arm_means(path,mtime_ns,size,state):
    """Cache only one state's arm values; file generations invalidate replacements."""
    with np.load(path,allow_pickle=False) as a:
        states=list(a['states'])
        if state not in states:return {},{}
        i=states.index(state);genes=a['genes']
        return tuple({str(g):float(value) for g,value in zip(genes,a[arm+'_mean_log2cp10k'][i])
                      if np.isfinite(value)} for arm in ('case','control'))


@lru_cache(maxsize=2)
def _edge_index(path):
    return json.loads(Path(path).read_text())

def integration_data(app,meta):
    if hasattr(app.state.job_store,'dataset'):
        return BundleIntegrationData(app.state.job_store.dataset(meta['job_id']),meta,
                                     (getattr(app.state,'assets',{}) or {}).get(meta['job_id'],{}))
    return UploadedIntegrationData(app,meta)
