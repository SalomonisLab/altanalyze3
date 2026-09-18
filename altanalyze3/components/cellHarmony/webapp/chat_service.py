"""Fast, conservative local routing and shared execution for sample-level Chat."""
from collections import OrderedDict
from copy import deepcopy
import re
import time
from altanalyze3.components.cellHarmony import chat_protocols as protocols
from altanalyze3.components.cellHarmony.grn_analysis import _sibling_comparison
from .chat_data import data_for_job

PROTOCOLS={'severity_gradient','dose_response','composition_shift','coexpression_module','annotation_concordance','donor_heterogeneity','most_affected_state','pathway_program'}

def _norm(text):return re.sub(r'[^a-z0-9]+',' ',str(text).lower()).strip()
def _contains(text,name):return bool(re.search(r'(?<!\w)'+re.escape(str(name).lower())+r'(?!\w)',text))

def read_question(question,states,genes,covariates):
    q=question.lower();names=[s for s in sorted(states,key=len,reverse=True) if _contains(q,s)]
    index={str(g).lower():str(g) for g in genes}
    found=list(dict.fromkeys(index[w.lower()] for w in re.findall(r'[A-Za-z0-9_.-]+',question) if w.lower() in index))
    columns=[c for c in covariates if not c.endswith(('__n_obs','__homogeneous')) and not c.startswith('wmean_')]
    exact=[c for c in columns if _contains(q,c) or _contains(_norm(q),_norm(c))]
    cov=max(exact,key=len) if exact else ''
    aliases={'fev1':['fev1_percent_predicted','fev1_liters'], 'fev1/fvc':['fev1_fvc_ratio'], 'dlco':['dlco_percent_predicted'],
             'age':['age_years','age','Age'], 'bmi':['bmi'], 'pack years':['smoking_pack_years'],
             'gold':['gold_ordinal','copd_gold_stage'],'condition':['condition','copd_status'], 'copd':['copd_status']}
    if not cov:
        for alias,candidates in sorted(aliases.items(),key=lambda x:-len(x[0])):
            if _contains(q,alias):
                cov=next((c for c in candidates if c in columns),'')
                if cov:break
    intent=''
    if re.search(r'\b(concordance|annotations?|labell?ings?)\b',q) and re.search(r'agree|agreement|compare|concordance|consistent|confus|match|map',q):intent='annotation_concordance'
    elif re.search(r'co[- ]?express|co[- ]?vary|correlat.*with',q) and found:intent='coexpression_module'
    elif re.search(r'most affected|most changed|which (cell )?(types?|states?).*(respond|affected|differential)',q):intent='most_affected_state'
    elif re.search(r'heterogeneity|signature.*(donor|sample|subset)|(?:donor|sample).*(?:signature|outlier)',q):intent='donor_heterogeneity'
    elif re.search(r'composition|abundance|deplet|enrich.*cell (state|type)|cell.*proportion|cell.*frequenc',q):intent='composition_shift'
    elif re.search(r'dose.response|stepwise|monotonic|across.*(?:stages|levels|grades)|ordered stages',q):intent='dose_response'
    elif re.search(r'\b(track|tracks|gradient|severity|correlate|correlation|associated)\b|change with|vary with',q):intent='severity_gradient'
    elif re.search(r'\b(go|goelite|go-elite|ontology|processes)\b',q):intent='pathway_program'
    # Common single-modality requests also avoid the external router.
    elif re.search(r'\bmarkers?\b',q) and not re.search(r'network|pathway|grn',q):intent='markers'
    elif len(names)>1 and re.search(r'distinguish|difference|compare|separat',q):intent='compare'
    elif found and re.search(r'express|where|distribution',q):intent='expression'
    elif re.search(r'genes?.*(significant|differential|changed|upregulated|downregulated)|differential.*genes?',q):intent='differential'
    if not intent:return None
    if intent=='dose_response' and cov=='gold_ordinal' and 'copd_gold_stage' in columns:cov='copd_gold_stage'
    return {'intent':intent,'cell_state':names[0] if names else '', 'cell_state_2':names[1] if len(names)>1 else '',
            'genes':found,'covariate':cov,'modality':'rna','router':'local_protocol',
            'direction':'down' if re.search(r'\bdown\b',q) else 'up' if re.search(r'\bup\b',q) else 'both'}


def resolve_contrast(ds,meta,question):
    current=(meta.get('differential') or {}).get('run_id','')
    entries=ds.deg_manifest().get('comparisons',[])
    candidates=[c for c in entries if c.get('modality','rna')=='rna' and c.get('kind')=='per_cell_state']
    q=_norm(question)
    full=[c for c in candidates if _norm(c.get('comparison','')) and _contains(q,_norm(c['comparison']))]
    if len(full)==1:return full[0]['id']
    return _sibling_comparison(ds,current,'rna')


def execute(app,meta,question,reading):
    if reading.get('intent') not in PROTOCOLS:return None
    started=time.perf_counter();intent=reading['intent'];state=reading.get('cell_state','');cov=reading.get('covariate','');genes=reading.get('genes') or []
    try:ds=data_for_job(app,meta)
    except (ValueError,KeyError,FileNotFoundError) as exc:
        return {'question':question,'reading':reading,'intent':intent,'status':'not_covered','answer':str(exc)}
    needed=[]
    if intent in {'severity_gradient','dose_response','coexpression_module','donor_heterogeneity','pathway_program'} and not state:needed.append('cell state')
    if intent in {'severity_gradient','dose_response','composition_shift'} and not cov:needed.append('clinical variable or grouping field')
    if intent=='coexpression_module' and not genes:needed.append('seed gene')
    if state and state not in ds.states:needed.append('valid cell state')
    needs_contrast=intent in {'most_affected_state','donor_heterogeneity','pathway_program'}
    contrast=(reading.get('contrast') or resolve_contrast(ds,meta,question)) if needs_contrast else ''
    if needs_contrast and not contrast:needed.append('completed RNA comparison')
    if needed:
        return {'question':question,'reading':reading,'intent':intent,'status':'clarify','answer':'Please specify '+', '.join(needed)+'.',
                'choices':{'states':ds.states,'covariates':list(ds.covariate_names())}}
    key=(intent,state,cov,tuple(genes),contrast,reading.get('direction','both'),_norm(question) if intent=='composition_shift' else '')
    with ds._lock:
        if not hasattr(ds,'answers'):ds.answers=OrderedDict()
        hit=key in ds.answers
        if hit:result=deepcopy(ds.answers[key]);ds.answers.move_to_end(key)
        else:
            try:
                if intent=='severity_gradient':result=protocols._run_severity_gradient(ds,state,cov)
                elif intent=='coexpression_module':result=protocols._run_coexpression(ds,state,genes[0])
                elif intent=='composition_shift':result=protocols._run_composition(ds,cov,question)
                elif intent=='dose_response':result=protocols._run_dose_response(ds,state,cov,genes)
                elif intent=='annotation_concordance':
                    data=getattr(ds,'ds',ds)
                    result=protocols._run_concordance(data,state)
                elif intent=='most_affected_state':result=protocols._run_most_affected_state(ds,contrast)
                elif intent=='donor_heterogeneity':
                    if hasattr(ds,'runs'):ds.comparison_group_field=ds.comparison_covariate(contrast)[0]
                    result=protocols._run_donor_heterogeneity(ds,state,contrast)
                else:
                    if hasattr(ds,'runs'):
                        assets={'differential':{k:{'goelite_tsv':(v.get('artifacts') or {}).get('goelite_tsv','')} for k,v in ds.runs.items()}}
                    else:assets=(getattr(app.state,'assets',{}) or {}).get(meta['job_id'],{})
                    result=protocols._run_pathway_program(assets,ds,state,contrast,reading.get('direction','both'))
            except (ValueError,KeyError,FileNotFoundError) as exc:
                result={'status':'not_covered','answer':str(exc)}
            ds.answers[key]=deepcopy(result)
            while len(ds.answers)>64:ds.answers.popitem(last=False)
    return dict(result,question=question,reading=dict(reading,contrast=contrast),intent=intent,
                provenance=ds.provenance,performance={'cache_hit':hit,'execution_ms':round(1000*(time.perf_counter()-started),2)})
