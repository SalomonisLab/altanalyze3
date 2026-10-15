"""Expandable methods and native per-state Interpret evidence from saved outputs.

No models are trained, no cells are excluded and no new cell annotations are made.
The GRN reader delegates gates and graph construction to the shared Discover code.
"""
from collections import Counter, defaultdict
from functools import cached_property
from types import SimpleNamespace
import json
from pathlib import Path
import re

import h5py
import numpy as np
from anndata.io import read_elem, sparse_dataset

from .. import analysis_record, discover_integration
from ..grn_analysis import _sibling_comparison
from .grn_data import completed_differentials, comparison_entry


def methods(meta):
    artifacts=meta.get('artifacts') or {}
    path=artifacts.get('combined_h5ad')
    # Resolve only the job's own directory, not an arbitrary bundle/reference location.
    directory=Path(path).parent.parent if path else Path(__file__).parent/'jobs'/meta['job_id']
    while directory.name != meta['job_id'] and directory != directory.parent:
        directory=directory.parent
    if directory.name != meta['job_id']:
        directory=Path(__file__).parent/'jobs'/meta['job_id']
    manifest,_=analysis_record.build_manifest(meta,directory)
    blocks=re.split(r'^## ',analysis_record.render_methods(manifest),flags=re.M)
    sections=[]
    for block in blocks[1:]:
        title,_,body=block.partition('\n')
        sections.append({'title':title.strip(),'text':body.strip()})
    sections.append({'title':'Biological summary and on-demand Discover interpretation',
        'text':('The summary reads completed differential tables, GO-Elite results and row-chunked observation/count metadata. '
                'Regulatory synthesis requires exact matching RNA and GRN-edge comparisons (state layer, sample field, unit and both arms). '
                'TF-activity differentials are included only when an exactly matching comparison was performed; otherwise their absence is reported. '
                'The native Discover network uses FDR ≤ 0.05, absolute RNA/edge fold change ≥ 1.2, '
                'and a factor-expression floor of 0.5 on sample-pseudobulk log2(1 + CP10k) in either arm. '
                'Up to 80 target features define each network; the five factors with the most retained targets are displayed as hubs. '
                'Hub-target overlaps with saved GO-Elite lists are distinguished from a new target-set enrichment. '
                'On-demand Interpret submits the full displayed hub-target set to native ToppGene GO Biological Process enrichment; unresolved mappings are reported. '
                'For matching RNA and lipid/metabolite comparisons, published WikiPathways diagrams are ranked by distinct reported RNA and molecular-feature nodes painted; '
                'the top three per state use native Discover interpretation, without an additional significance filter. '
                'Predicted modality nodes remain predictions, and diagram coverage is not a measurement of pathway flux. '
                'DrugMatch is run only when requested in Chat. It uses the saved RNA signature at FDR ≤ 0.05 and absolute log2 fold change ≥ 0, '
                'with the native 1,000-feature cap and measured/all-gene modes. Signed Spearman correlation, tied ranks and direction counts follow '
                'Discover cellref2-dart-like-signature-v1, one row per source treatment experiment; no compound-average signature or inferential FDR is generated. '
                'Library version and query coverage are returned with each screen. Drug Interpret retrieves exact-compound registry facts, '
                'ClinicalTrials.gov trial records using the annotated case condition, and PubMed citations. '
                'Literature checks use exact feature identities and documented state aliases, with explicitly broader disease-context searches. '
                'Publication hits require claim-level review before an exact disease/state/direction finding can be classified as verified. '
                'For Both mode, descriptive reference-gap screening reports every cluster and requires at least 50 cells, 70% unassigned within a candidate, '
                'and enrichment above the overall unassigned fraction; concentrated candidates additionally contain at least 50% of all unassigned cells. '
                'These are screening criteria, not significance tests or proof of cell-state novelty.')})
    models=[]
    for modality,record in manifest['model_versions'].items():
        models.append({'modality':modality,'methods_text':analysis_record.model_methods(modality,record),'execution_record':record,
            'training_and_validation':{key:record[key] for key in (
                'training_datasets','training_provenance','training_summary','validation',
                'validation_summary','evaluation','confidence','declared_method') if key in record},
            'confidence_note':'Training-fit metrics do not establish independent predictive accuracy. Missing training/holdout provenance is unresolved; current model metadata is not substituted for this job’s version.'})
    return {'sections':sections,'models':models,'provenance_gaps':manifest['provenance_gaps']}


def missing_states(meta,path):
    from .biology_summary import _column,_length,CHUNK
    if meta.get('analysis_mode') != 'both':return {'applicable':False,'rows':[]}
    if not path or not Path(path).is_file():return {'applicable':True,'status':'observations_unavailable','rows':[]}
    pairs=Counter();names=defaultdict(Counter)
    with h5py.File(path,'r') as f:
        obs=f['obs'];n=_length(obs[obs.attrs.get('_index','_index')])
        required=('supervised_state','unsupervised_cluster')
        if not all(k in obs for k in required):
            return {'applicable':True,'status':'required_layers_unavailable','required':list(required),'rows':[]}
        for start in range(0,n,CHUNK):
            supervised=_column(obs[required[0]],start,start+CHUNK).astype(str)
            clusters=_column(obs[required[1]],start,start+CHUNK).astype(str)
            states=_column(obs['unsupervised_state'],start,start+CHUNK).astype(str) if 'unsupervised_state' in obs else clusters
            pairs.update(zip(clusters,supervised))
            for cluster,state in zip(clusters,states):names[cluster][state]+=1
    nulls={'unaligned','unassigned','unknown','unk','nan','(missing)',''}
    total=sum(pairs.values());unaligned=sum(v for (_,s),v in pairs.items() if s.casefold() in nulls)
    base=unaligned/total if total else 0
    rows=[]
    for cluster in names:
        members={s:v for (c,s),v in pairs.items() if c==cluster}
        n=sum(members.values());u=sum(v for s,v in members.items() if s.casefold() in nulls)
        invalid=cluster.casefold() in nulls|{'not clustered'}
        fraction=u/n if n else 0;share=u/unaligned if unaligned else 0
        rows.append({'cluster':cluster,'names':dict(names[cluster]),'cells':n,'unassigned_cells':u,
            'unassigned_fraction':fraction,'share_of_all_unassigned':share,
            'enrichment_over_dataset':fraction/base if base else None,
            'supervised_assignments':members,
            'candidate':not invalid and n>=50 and fraction>=.7 and fraction>base,
            'concentrated':not invalid and n>=50 and fraction>=.7 and share>=.5})
    rows.sort(key=lambda r:(not r['concentrated'],not r['candidate'],-r['unassigned_cells'],r['cluster']))
    return {'applicable':True,'status':'completed','rows':rows,'dataset_unassigned_fraction':base,
        'criteria':{'min_cluster_cells':50,'min_within_cluster_unassigned':.7,'min_share_for_concentrated':.5,
                    'type':'descriptive screening criteria; no enrichment p-value or proven novelty'},
        'recommendation':'For supported concentrated candidate states, evaluate ICGS-integrate in AltAnalyze3 with the reference as primary and these ICGS3 clusters as secondary. ICGS-integrate tests marker distinctness and correlation, retains nonredundant source clusters with provenance, and records redundant clusters as support. It does not merge clusters or perform latent batch correction.',
        'interpretation':'A concentrated unassigned cluster is a candidate reference gap, not proof of a novel cell type. Review its GO-Elite name, markers, QC, sample composition and alternatives before integration. All cells remain in the combined H5AD.'}


class SavedInterpretData:
    """Native Discover adapter with row-chunked counts and saved GRN aggregates.

    No in-memory AnnData, full expression matrix, or model pickle is opened.
    """
    def __init__(self,meta,contrast):
        self.meta=meta;self.runs=completed_differentials(meta);self.current_contrast=contrast
        self.path=(meta.get('artifacts') or {}).get('combined_h5ad')
        self.entries=[comparison_entry(k,v) for k,v in reversed(list(self.runs.items()))]
        self.layer=self.runs[contrast]['config']['population_col'];self.tables={}
        self._edge_cache={};self._arms={}

    def deg_manifest(self):return {'comparisons':self.entries}
    def comparison(self,contrast):
        chosen=next((c for c in self.entries if c['id']==(contrast or self.current_contrast)),{})
        return chosen.get('comparison',contrast)
    def _run_id(self,contrast):return contrast if contrast in self.runs else self.current_contrast
    def features(self,modality):
        from ..modalities import modality_artifacts
        canonical='lipids' if modality=='lipid' else modality
        path=self.path if canonical=='rna' else modality_artifacts(self.meta).get(canonical,{}).get('h5ad')
        if not path:return []
        with h5py.File(path,'r') as f:return list(read_elem(f['var']).index.astype(str))

    def available(self,contrast,modality,state=''):
        return bool(_sibling_comparison(self,self._run_id(contrast),'lipids' if modality=='lipid' else modality))
    def differentials(self,contrast,state,modality):
        run=_sibling_comparison(self,self._run_id(contrast),'lipids' if modality=='lipid' else modality)
        if not run:return []
        if run not in self.tables:
            from .biology_summary import _rows,finite
            table=defaultdict(list)
            path=next((v for k,v in (self.runs[run].get('artifacts') or {}).items() if k.startswith('DEG_detailed_')),None)
            for row in _rows(path):
                row={k:(None if isinstance(v,float) and not np.isfinite(v) else v) for k,v in row.items()}
                for key in ('log2fc','fdr','pval'):row[key]=finite(row.get(key))
                table[str(row['population'])].append(row)
            self.tables[run]=table
        return self.tables[run].get(state,[]) if state else [r for rows in self.tables[run].values() for r in rows]

    @cached_property
    def states(self):
        from .biology_summary import _column,_length
        with h5py.File(self.path,'r') as f:
            obs=f['obs'];n=_length(obs[obs.attrs.get('_index','_index')])
            return list(dict.fromkeys(map(str,_column(obs[self.layer],0,n))))

    def edge_scores(self,state):
        from .biology_summary import _column,_length
        from ..modalities import modality_artifacts
        if state in self._edge_cache:return self._edge_cache[state]
        record=modality_artifacts(self.meta).get('grn',{})
        path=record.get('network_h5ad') or record.get('pseudobulk_h5ad')
        if not path:raise FileNotFoundError('Saved sample × state GRN aggregate not recorded')
        with h5py.File(path,'r') as f:
            obs=f['obs'];n=_length(obs[obs.attrs.get('_index','_index')])
            from .grn_data import aggregate_state_column
            layer=aggregate_state_column(self.meta,record,SimpleNamespace(obs=read_elem(obs)),self.layer)
            idx=np.flatnonzero(_column(obs[layer],0,n).astype(str)==state)
            names=read_elem(f['var']).index.astype(str);total=np.zeros(len(names))
            matrix=f['X'] if isinstance(f['X'],h5py.Dataset) else sparse_dataset(f['X'])
            for i in idx:
                block=matrix[i:i+1];total+=np.asarray(block.toarray() if hasattr(block,'toarray') else block).ravel()
        result=dict(zip(names,total/len(idx))) if len(idx) else {}
        # Only the current state's edge vector survives between calls.
        self._edge_cache={state:result};return result

    def arm_expression(self,contrast,state,features):
        """Same sample pseudobulk log2(1 + CP10k) as Discover; bounded CSR blocks."""
        from .biology_summary import _column,_length
        features=tuple(sorted(features));key=(state,features)
        if key in self._arms:return self._arms[key]
        run=_sibling_comparison(self,self._run_id(contrast),'rna');cfg=self.runs[run]['config']
        with h5py.File(self.path,'r') as f:
            obs=f['obs'];n=_length(obs[obs.attrs.get('_index','_index')]);names=read_elem(f['var']).index.astype(str)
            if 'counts' not in f.get('layers',{}):raise ValueError('Native GRN Interpret requires retained raw counts for sample pseudobulks')
            sample=next((s for s in ('sample','Sample','Library','library','Donor','donor') if s in obs),None)
            if sample is None:raise ValueError('Native GRN Interpret requires a sample observation field')
            labels=_column(obs[self.layer],0,n).astype(str);groups=_column(obs[cfg['sample_field']],0,n).astype(str)
            samples=_column(obs[sample],0,n).astype(str)
            node=f['layers/counts'];matrix=node if isinstance(node,h5py.Dataset) else sparse_dataset(node)
            selected=[i for i,g in enumerate(names) if g in features];sums={};totals={}
            for start in range(0,n,512):
                # Read only blocks holding a cell from this state and a comparison arm.
                stop=min(n,start+512);wanted=(labels[start:stop]==state)&np.isin(groups[start:stop],cfg['group1_samples']+cfg['group2_samples'])
                if not wanted.any():continue
                block=matrix[start:stop][wanted];btotal=np.asarray(block.sum(axis=1)).ravel()
                partial=block[:,selected];partial=partial.toarray() if hasattr(partial,'toarray') else partial
                for g,s,t,values in zip(groups[start:stop][wanted],samples[start:stop][wanted],btotal,partial):
                    ident=(g,s);sums.setdefault(ident,np.zeros(len(selected)));sums[ident]+=values;totals[ident]=totals.get(ident,0)+t
            out=[]
            for arm in ('group1_samples','group2_samples'):
                means=[np.log2(1+10000*v/totals[k]) for k,v in sums.items() if k[0] in cfg[arm] and totals[k]>0]
                mean=np.mean(means,axis=0) if means else None
                out.append(dict(zip(names[selected],map(float,mean))) if mean is not None else {})
        self._arms[key]=tuple(out);return tuple(out)


def regulatory_sections(meta,digest):
    """Evaluate matching RNA/edge contrasts; annotate TF evidence only if available."""
    groups=defaultdict(dict)
    for c in digest['comparisons']:groups[c['contrast']][c['modality']]=c
    result=[]
    for same in groups.values():
        if not all(m in same for m in ('rna','grn')):continue
        c=same['rna'];ds=SavedInterpretData(meta,c['id']);layer=c['config']['population_col']
        states=digest['observations']['layers'].get(layer,{}).get('counts',{})
        for state,n in states.items():
            item={'state':state,'layer':layer,'cells':n,'comparison':c['comparison'],'run_id':c['id'],
                'grn_run_id':same['grn']['id'],'tf_run_id':same.get('grn_tf',{}).get('id',''),
                'tf_activity_status':'matching differential available' if 'grn_tf' in same else 'No matching TF-activity differential was performed; TF-activity changes are not inferred.', 'hubs':[]}
            try:
                net=discover_integration.network(ds,contrast=c['id'],cell_state=state,limit=80,significance='fdr')
                by_source=defaultdict(list);nodes={v['id']:v for v in net.get('nodes',[])}
                for e in net.get('edges',[]):by_source[e['source']].append(e)
                for factor,edges in sorted(by_source.items(),key=lambda kv:(-len(kv[1]),kv[0]))[:5]:
                    targets=[e['target'] for e in edges];own=nodes.get(factor,{})
                    item['hubs'].append({'factor':factor,'n_targets':len(targets),'targets':targets,
                        'target_RNA_up':sum((nodes.get(t,{}).get('log2fc') or 0)>0 for t in targets),
                        'target_RNA_down':sum((nodes.get(t,{}).get('log2fc') or 0)<0 for t in targets),
                        'target_RNA_calls':{t:{'log2fc':nodes.get(t,{}).get('log2fc'),'fdr':nodes.get(t,{}).get('fdr')} for t in targets},
                        'edge_up':sum((e.get('log2fc') or 0)>0 for e in edges),
                        'edge_down':sum((e.get('log2fc') or 0)<0 for e in edges),
                        'RNA_log2fc':own.get('log2fc'),'RNA_fdr':own.get('fdr'),
                        'TF_log2fc':own.get('activity_log2fc'),'TF_fdr':own.get('activity_fdr'),
                        'saved_pathway_overlaps':[{'term':p['facts'].get('term_name'),'genes':sorted(set(targets)&set(re.split(r'[,;\s]+',str(p['facts'].get('overlap_genes',''))))),'fdr':p['facts'].get('fdr')}
                            for p in digest['evidence'] if p['kind']=='pathway' and p['run_id']==c['id'] and re.sub(r'__(up|down)$','',p['state'])==state and set(targets)&set(re.split(r'[,;\s]+',str(p['facts'].get('overlap_genes',''))))]})
                item.update(status='completed',interpretation=net.get('interpretation') or net.get('note'),
                    n_edges=len(net.get('edges',[])),method={'target_limit':80,'gene_fold':1.2,'edge_fold':1.2,'max_fdr':.05,'TF_expression_floor':.5,'implementation':'shared Discover network'})
            except (ValueError,KeyError,FileNotFoundError,OSError) as exc:
                item.update(status='unavailable',interpretation=str(exc))
            result.append(item)
    result.sort(key=lambda r:(not bool(r['hubs']),-max((h['n_targets'] for h in r['hubs']),default=0),r['state']))
    return {'rows':result,'states_examined':len(result),'selection':'Top hubs by number of native displayed targets; no causal importance or independent validation inferred.',
        'scope':'Requires matching RNA and GRN-edge comparisons: population layer, sample field, comparison unit and both arms. Matching TF-activity differentials are included when available; their absence does not prevent RNA/edge interpretation.'}


def pathway_sections(meta,digest):
    groups=defaultdict(dict)
    for c in digest['comparisons']:groups[c['contrast']][c['modality']]=c
    rows=[];coverage=[]
    for same in groups.values():
        if 'rna' not in same:continue
        c=same['rna'];ds=SavedInterpretData(meta,c['id']);layer=c['config']['population_col']
        for modality in ('lipids','metabolite'):
            if modality not in same:continue
            native_modality='lipid' if modality=='lipids' else modality
            for state in digest['observations']['layers'].get(layer,{}).get('counts',{}):
                record=dict(run_id=c['id'],state=state,modality=native_modality,layer=layer)
                try:
                    ranked=discover_integration.pathways(ds,contrast=c['id'],cell_state=state,modality=native_modality,significance='reported')
                    for pathway in ranked.get('pathways',[])[:3]:
                        diagram=discover_integration.pathway(ds,id=pathway['id'],contrast=c['id'],cell_state=state,modality=native_modality,significance='reported')
                        rows.append(dict(pathway,run_id=c['id'],modality=native_modality,comparison=c['comparison'],state=state,layer=layer,
                            interpretation=diagram.get('interpretation') or diagram.get('note') or diagram.get('error'),
                            nodes=[n for n in diagram.get('nodes',[]) if n.get('log2fc') is not None]))
                    record.update(status='completed',eligible_diagrams=len(ranked.get('pathways',[])),note=ranked.get('note'))
                except (ValueError,KeyError,FileNotFoundError,OSError) as exc:
                    record.update(status='unavailable',note=str(exc))
                coverage.append(record)
    return {'rows':rows,'coverage':coverage,'states_examined':len(coverage),'method':'Native Discover ranking and interpretation of published WikiPathways diagrams; top three per state and matching comparison. Saved calls are overlaid without an additional filter; no metabolic flux inferred.'}
