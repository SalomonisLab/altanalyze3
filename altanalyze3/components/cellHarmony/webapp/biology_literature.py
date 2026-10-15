"""Auditable PubMed checks of job-derived hypotheses; search absence is not novelty."""
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor
import re
import urllib.parse

from ..discover_resources.gene_literature import _get,EUTILS,phrase_terms,pubmed_url


def query_terms(feature,state,disease,species):
    # Exact symbols are kept distinct (SPI1 is never SP1). No guessed synonyms.
    feature=re.sub(r'["\[\]()]',' ',feature).strip()
    states=state if isinstance(state,tuple) else (state,)
    states=tuple(re.sub(r'[_-][cC]\d+$','',s).strip() for s in states)
    pieces=[phrase_terms(feature),phrase_terms(disease),' OR '.join(phrase_terms(s) for s in states if s)]
    return ' AND '.join(f'({v})' for v in pieces if v)


@lru_cache(maxsize=128)
def search(feature,state,disease,species):
    query=query_terms(feature,state,disease,species)
    if not feature:return {'status':'not_searched','papers':[],'note':'No feature-specific literature query defined.'}
    # Explicit broader context searches prevent abbreviated atlas labels from
    # hiding related work. They retrieve candidates, never verify a claim.
    queries=[('feature_disease_state',query),('feature_disease',query_terms(feature,'',disease,species)),
             ('feature_state_background',query_terms(feature,state,'',species))]
    queries=list(dict((q,(scope,q)) for scope,q in queries).values())
    searches=[];ids=[];scopes={}
    for scope,q in queries:
        found=_get(EUTILS+'esearch.fcgi?'+urllib.parse.urlencode({'db':'pubmed','retmode':'json','retmax':3,'sort':'relevance','term':q}))['esearchresult']
        searches.append({'scope':scope,'query':q,'count':int(found.get('count',0)),'url':pubmed_url(q)})
        ids.extend(p for p in found.get('idlist',[]) if p not in ids)
        for pmid in found.get('idlist',[]):scopes.setdefault(pmid,[]).append(scope)
    papers=[]
    if ids:
        summaries=_get(EUTILS+'esummary.fcgi?'+urllib.parse.urlencode({'db':'pubmed','retmode':'json','id':','.join(ids)}))['result']
        for pmid in ids:
            p=summaries.get(pmid,{})
            papers.append({'pmid':pmid,'title':p.get('title',''),'year':p.get('pubdate',''),
                'journal':p.get('fulljournalname',''),'url':f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/','retrieval_scopes':scopes.get(pmid,[])})
    return {'status':'searched','query':query,'search_url':pubmed_url(query),'searches':searches,'count':searches[0]['count'],
        'papers':papers,'classification':'prior_literature_candidates' if papers else 'unverified_prediction',
        'note':'Retrieved publications are candidates for verifying the exact finding, not proof of its direction, species or mechanism. Zero hits do not establish novelty.'}


def review(digest,limit=8):
    comparisons={c['id']:c for c in digest['comparisons']};rows=[];seen=set()
    for p in digest.get('regulatory_review',{}).get('rows',[]):
        for hub in p['hubs'][:1]:
            key=(p['run_id'],p['state'],hub['factor'])
            if key not in seen:
                seen.add(key);rows.append({'run_id':p['run_id'],'state':p['state'],'feature':hub['factor'],'modality':'grn','observation':p['interpretation']})
    for e in digest['evidence']:
        if e['kind']!='differential_feature' or e['modality']=='cell_communication':continue
        feature=e['feature'].split('|')[0] if e['modality']=='grn' else e['feature']
        key=(e['run_id'],e['state'],feature)
        if key not in seen:
            seen.add(key);rows.append({'run_id':e['run_id'],'state':e['state'],'feature':feature,'modality':e['modality'],'evidence_id':e['id'],'observation':e['title'],
                'effect':e.get('facts',{}).get('log2fc')})
    # Modalities each get a place before subsequent results of one modality.
    selected=[]
    for m in digest['modalities']:
        selected.extend([r for r in rows if r['modality']==m][:1])
    # Give distinct states literature coverage before adding more genes from one
    # state. Full job evidence remains available regardless of prompt ranking.
    covered={r['state'] for r in selected}
    for r in rows:
        if r['state'] not in covered:selected.append(r);covered.add(r['state'])
    buckets={}
    for r in rows:
        if r not in selected:buckets.setdefault((r['state'],r['modality']),[]).append(r)
    while any(buckets.values()):
        for bucket in buckets.values():
            if bucket:selected.append(bucket.pop(0))
    def check(r):
        c=comparisons[r['run_id']];disease=c['comparison'].split(' vs ')[0]
        row=dict(r,comparison=c['comparison'],species=digest['species'],disease=disease)
        from .biology_summary import _normal
        search_disease=next((v['name'] for v in digest.get('covariate_biology',[])
                            if _normal(disease) in {_normal(v['term']),_normal(v['name'])}),disease)
        row['literature_disease_term']=search_disease
        background=next((v for v in digest.get('cell_state_biology',{}).get('matched',[]) if v['state']==r['state']),{})
        aliases=tuple(dict.fromkeys([r['state']]+background.get('card',{}).get('names',[])))
        row['literature_state_terms']=list(aliases)
        try:row['literature']=search(r['feature'],aliases,search_disease,str(digest['species']))
        except Exception as exc:row['literature']={'status':'unavailable','papers':[],'error':type(exc).__name__,'note':'Literature could not be checked; no novelty inference made.'}
        return row
    # Independent HTTP waits can overlap. Every request still passes through the
    # native shared NCBI rate limiter; map preserves the original result order.
    with ThreadPoolExecutor(max_workers=4) as pool:
        checked=list(pool.map(check,selected[:limit]))
    # Exact verified observations must have a reviewed claim/source mapping. The
    # existing CellCards are background, not proof of this disease differential.
    known=[];predictions=[]
    for row in checked:
        sources=verified_sources(row,digest.get('reviewed_observations',[]))
        if sources:known.append(dict(row,verification_sources=sources))
        else:predictions.append(row)
    return {'verified_known_biology':known, 'literature_candidates':[r for r in checked if r['literature'].get('papers')],
        'novel_predictions':predictions,'coverage':{'checked':len(checked),'eligible':len(rows),'limit':limit},
        'verification_note':'Exact disease × cell-state × feature × direction findings are not labelled verified from a title match. Novel predictions are testable hypotheses with novelty unconfirmed; literature queries and citation candidates are retained.'}


def verified_sources(row,records):
    """Claim-level provenance gate; never infer approval from a PubMed hit."""
    from .biology_summary import _normal,finite
    effect=finite(row.get('effect'))
    if effect is None or effect==0:return []
    identity={k:_normal(row.get(k,'')) for k in ('species','state','disease','feature','modality')}
    identity['direction']='up' if effect>0 else 'down'
    sources=[]
    for record in records:
        if record.get('status')!='reviewed':continue
        if not all(record.get(k) for k in ('reviewed_by','reviewed_at','claim','evidence_excerpt','citation_url')):continue
        url=urllib.parse.urlparse(record['citation_url'])
        if url.scheme!='https' or not url.netloc:continue
        if any(_normal(record.get(k,''))!=v for k,v in identity.items()):continue
        sources.append(record)
    return sources
