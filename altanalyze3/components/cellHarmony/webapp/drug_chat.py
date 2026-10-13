"""Uploaded RNA signatures screened by the unchanged native Discover DrugMatch engine."""
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
import json
import re
import threading
from pathlib import Path

from ..discover_resources import drug_matching as dm,drug_facts
from .grn_data import completed_differentials,comparison_entry
from .biology_summary import _rows,finite,source_identity

_LOCK=threading.Lock()


def resource_identity():
    """Invalidate response/score caches when mounted resources or adapters change."""
    paths=[Path(dm.__file__),Path(drug_facts.__file__),Path(__file__)]
    paths.append(Path(dm._store_dir())/'manifest.json')
    return json.dumps([(str(p),p.stat().st_mtime_ns,p.stat().st_size) for p in paths if p.is_file()])


def requested(question):
    return bool(re.search(r'\b(drug|drugs|drugmatch|ilincs|compound|compounds)\b',question,re.I) and re.search(r'revers|oppos|mimic|signature|drugmatch|ilincs|match',question,re.I))


def choices(meta):
    rows=[]
    for key,run in completed_differentials(meta).items():
        if (run.get('config') or {}).get('modality','rna') != 'rna':continue
        path=next((v for k,v in (run.get('artifacts') or {}).items() if k.startswith('DEG_detailed_')),None)
        states=list(dict.fromkeys(str(r['population']) for r in _rows(path)))
        rows.append(dict(comparison_entry(key,run),states=states))
    return rows


def read(meta,question):
    runs=choices(meta);matches=[]
    for run in runs:
        states=[s for s in run['states'] if re.search(r'(?<!\w)'+re.escape(s)+r'(?!\w)',question,re.I)]
        for state in states:matches.append((run,state))
    # A longest exact state match prevents AT2 matching AT2-prolif as a second state.
    if matches:
        length=max(len(s) for _,s in matches);matches=[m for m in matches if len(m[1])==length]
    named=[m for m in matches if m[0]['comparison'].casefold() in question.casefold() or m[0]['id'] in question]
    if named:matches=named
    if len(matches)!=1:
        return {'kind':'drug_match_choices','intent':'drug_match','answer':'Select the completed RNA comparison and cell state for this signature. No comparison or cell state is guessed.','choices':runs,
            'tail':'mimicking' if re.search(r'mimic',question,re.I) else 'opposing','gene_space':'measured' if re.search(r'ilincs|measured',question,re.I) else 'all'}
    run,state=matches[0]
    return {'run_id':run['id'],'cell_state':state,'tail':'mimicking' if re.search(r'mimic',question,re.I) else 'opposing',
            'gene_space':'measured' if re.search(r'ilincs|measured',question,re.I) else 'all'}


def prepare(meta,run_id,state,gene_space='all',max_fdr=.05,min_abs_log2fc=0,gene_cap=1000):
    run=completed_differentials(meta).get(run_id)
    if not run or (run.get('config') or {}).get('modality','rna')!='rna':raise ValueError('Choose a completed RNA differential comparison')
    if gene_space not in dm.GENE_SPACES:raise ValueError('Gene space must be all or measured')
    if not 0<=max_fdr<=1 or min_abs_log2fc<0 or not 1<=gene_cap<=10000:raise ValueError('Invalid signature selection parameters')
    path=next((v for k,v in (run.get('artifacts') or {}).items() if k.startswith('DEG_detailed_')),None)
    genes=[];effects=[]
    for row in _rows(path):
        if str(row['population'])!=state:continue
        fc,q=finite(row.get('log2fc')),finite(row.get('fdr'))
        if fc is not None and q is not None and q<=max_fdr and abs(fc)>=min_abs_log2fc:
            genes.append(str(row['gene']));effects.append(fc)
    if not genes:raise ValueError('This state has no saved RNA calls satisfying the declared signature criteria')
    if any(re.match(r'^ENS[A-Z]*G\d+',g) for g in genes):
        raise ValueError('Ensembl identifiers require a verified symbol mapping before DrugMatch. All query features are retained; no symbol subset is substituted.')
    return dm.prepare_query(genes,effects,job_id=meta['job_id'],contrast=run_id,cell_state=state,modality='rna',
        max_fdr=max_fdr,min_abs_log2fc=min_abs_log2fc,gene_cap=gene_cap,gene_space=gene_space,source_identity=source_identity(meta)),run


@lru_cache(maxsize=4)
def _screen(key,genes,effects,parameters,gene_space,store_version):
    prepared=dm.prepare_query(list(genes),list(effects),**json.loads(parameters))
    return dm.score_library(prepared,gene_space=gene_space)


def screen(meta,run_id,cell_state,tail='opposing',gene_space='all',offset=0,limit=25,compound='',max_fdr=.05,min_abs_log2fc=0,gene_cap=1000):
    prepared,run=prepare(meta,run_id,cell_state,gene_space,max_fdr,min_abs_log2fc,gene_cap)
    # Native partition memmaps and result arrays are shared safely; bounded score cache.
    with _LOCK:
        result=_screen(prepared.hash,tuple(prepared.genes),tuple(prepared.effects),json.dumps(prepared.params,sort_keys=True),gene_space,resource_identity())
        rows,n_tail,n_match=dm.page(result,tail,offset,limit,compound)
    payload={k:v for k,v in result.items() if not k.startswith('_')}
    payload.update(n_query_features_offered=prepared.n_offered,duplicate_query_features_resolved=prepared.duplicate_features,kind='drug_match',intent='drug_match',run_id=run_id,cell_state=cell_state,tail=tail,
        offset=offset,limit=limit,compound_search=compound,n_in_tail=n_tail,n_matched=n_match,rows=rows,
        comparison=comparison_entry(run_id,run)['comparison'],layer=run['config'].get('population_col'),
        species=meta.get('species'),query_parameters=prepared.params,
        provenance=dm.provenance(),
        answer=f"{n_tail:,} source signatures {tail} the saved RNA program in {cell_state}. One row per source experiment; cell line, dose, time and species are retained. Spearman correlation and direction counts are separate; reversal is not demonstrated treatment efficacy.")
    return payload


def interpret(meta,**params):
    payload=screen(meta,**params);run=completed_differentials(meta)[params['run_id']]
    disease=run.get('case_label') or ''
    # Prefer an exact disease vocabulary record; otherwise retain the annotated arm
    # verbatim, expose it, and never turn the reference name into a disease.
    from .biology_summary import knowledge,_normal
    for row in knowledge().get('diseases',[]):
        if _normal(disease) in {_normal(row['term']),_normal(row['name'])}:disease=row['name'];break
    named={}
    for r in payload['rows']:
        if drug_facts.nameable(r.get('compound','')):named.setdefault(r['compound'].casefold(),r)
    # As Discover does, extend the tail to find named compounds, explicitly mark
    # those outside the displayed page. No score aggregation occurs.
    if len(named)<8 and not params.get('compound'):
        more=screen(meta,**dict(params,offset=payload['offset']+len(payload['rows']),limit=500))
        for r in more['rows']:
            if drug_facts.nameable(r.get('compound','')):named.setdefault(r['compound'].casefold(),dict(r,on_page=False))
            if len(named)>=8:break
    rows=list(named.values())[:8]
    with ThreadPoolExecutor(max_workers=4) as pool:
        facts=list(pool.map(lambda r:drug_facts.facts(r['compound'],disease,str(r.get('pubchem_cid') or '')),rows))
    return {'kind':'drug_interpretation','comparison':payload['comparison'],'cell_state':payload['cell_state'],
        'disease_query':disease,'tail':payload['tail'],'gene_space':payload['gene_space'],
        'compounds':[{'compound':r['compound'],'rank':r['rank_in_tail'],'rho':r['rho'],
            'source_signature_id':r['source_signature_id'],'on_page':r.get('on_page',True),
            'context':{k:r.get(k) for k in ('cell_line','experimental_species','dose','time','source_url')},'facts':f,
            'summary':fact_summary(r,f,disease,payload['cell_state'],payload['tail'])} for r,f in zip(rows,facts)],
        'unnamed_on_page':[r['compound'] for r in payload['rows'] if not drug_facts.nameable(r.get('compound',''))],
        'interpretation':'Registry and clinical-trial facts describe these exact compounds and the explicitly displayed disease query. Trial inclusion, FDA approval in another indication and signature matching do not establish efficacy in this dataset.'}


def fact_summary(row,facts,disease,state,tail):
    """Readable rendering of retrieved facts; unavailable queries stay unavailable."""
    fd=facts.get('fda') or {};tr=facts.get('trials') or {};pm=facts.get('pubmed') or {}
    parts=[f"{row['compound']} {tail} the {state} signature (Spearman {row['rho']:+.3f}) in the source experiment: {row.get('cell_line') or 'cell line unrecorded'}, {row.get('dose') or 'dose unrecorded'}, {row.get('time') or 'time unrecorded'}. This is signature similarity, not efficacy in that cell state."]
    approval=fd.get('fda_approved')
    parts.append('An FDA approval record was retrieved.' if approval is True else 'No FDA approval record was found by the recorded identity query.' if approval is False else 'FDA approval status is unresolved.')
    if fd.get('boxed_warning'):parts.append('Recorded boxed warning: '+str(fd['boxed_warning']))
    if tr.get('errors'):parts.append('Trial retrieval is incomplete: '+ '; '.join(tr['errors']))
    elif 'n_in_condition' in tr:parts.append(f"{tr['n_in_condition']} registered trials match the explicit condition query '{disease}'; {tr.get('n_trials',0)} match the compound overall. Registration does not establish a positive trial result.")
    if pm.get('error'):parts.append(pm['error'])
    elif (pm.get('themes') or {}).get('disease'):parts.append(f"{pm['themes']['disease']['count']} PubMed records match compound and disease; these are literature matches, not verified benefit.")
    return ' '.join(parts)
