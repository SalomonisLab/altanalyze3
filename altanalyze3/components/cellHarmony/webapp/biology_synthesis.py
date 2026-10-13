"""Synthesize retained results; no fitting, imputation, tests or new annotations.

Cross-modality joins require the same saved contrast and cell-state identifier.
RNA-derived predictions count as additional hypotheses, not independent replication.
"""
from collections import defaultdict
import re

from .grn_data import completed_differentials, comparison_entry


def population(name):
    match=re.match(r'^(.*)__(up|down)$',str(name))
    return (match.group(1),match.group(2)) if match else (str(name),'mixed/unspecified')


def edge_genes(feature):
    parts=str(feature).split('|')
    return parts if len(parts)==2 and all(parts) else []


def enrich_regulatory_targets(digest,meta,rows,add):
    """Read exact target/TF RNA calls for retained edges, including weaker-ranked calls."""
    from .biology_summary import finite, clean
    comparisons={c['id']:c for c in digest['comparisons']}
    wanted=defaultdict(set)
    for e in digest['evidence']:
        if e['kind']=='pathway' and e['run_id'] in comparisons:
            genes=set(re.split(r'[,;\s]+',str(e['facts'].get('overlap_genes',''))))-{''}
            wanted[comparisons[e['run_id']]['contrast']].update((population(e['state'])[0],g) for g in genes)
        if e['modality']=='grn' and e['kind']=='differential_feature':
            genes=edge_genes(e['feature'])
            if genes:wanted[comparisons[e['run_id']]['contrast']].update((e['state'],g) for g in genes)
    existing={(e['run_id'],e['state'],e['feature']) for e in digest['evidence'] if e['kind']=='differential_feature'}
    for run_id,run in completed_differentials(meta).items():
        modality=(run.get('config') or {}).get('modality','rna')
        if modality not in ('rna','grn_tf'):continue
        c=comparisons[run_id];needed=wanted.get(c['contrast'],set())
        if not needed:continue
        path=next((v for k,v in (run.get('artifacts') or {}).items() if k.startswith('DEG_detailed_')),None)
        for row in rows(path):
            state,feature=clean(row.get('population')),clean(row.get('gene'))
            if (state,feature) not in needed or (run_id,state,feature) in existing:continue
            fc,fdr=finite(row.get('log2fc')),finite(row.get('fdr'))
            if fc is None or fdr is None or fdr>.05:continue
            add('differential_feature',f'{feature} in {state} ({c["comparison"]})',
                {'state':state,'feature':feature,'log2fc':fc,'fdr':fdr,'n_case':finite(row.get('n_case')),
                 'n_control':finite(row.get('n_control')),'selection':f'exact native {modality} support for a retained regulatory edge or pathway gene'},
                modality=modality,layer=c['config'].get('population_col',''),state=state,run_id=run_id,feature=feature,mode='volcano')
            existing.add((run_id,state,feature))


def change(item):
    value=item['facts'].get('delta_score') if item['kind']=='communication_contrast' else item['facts'].get('log2fc')
    return 'higher' if value is not None and value>0 else 'lower' if value is not None and value<0 else 'unchanged'


def synthesize(digest):
    comparisons={c['id']:c for c in digest['comparisons']}
    groups=defaultdict(list);patterns=defaultdict(dict);limits=[]
    for c in comparisons.values():
        cfg=c['config'];field=cfg.get('sample_field','')
        a,b=cfg.get('group1_samples',[]),cfg.get('group2_samples',[])
        if len(a)<2 or len(b)<2:
            text=f'{c["comparison"]}: fewer than two selected {field or "sample"} values in an arm. Cell-level FDR does not establish replication across biological samples; donor identity is not inferred from cell counts.'
            if text not in limits:limits.append(text)
    for e in digest['evidence']:
        if not e['run_id']:continue
        contrast=comparisons[e['run_id']]['contrast']
        if e['kind']=='differential_patterns':
            for state,counts in e['facts']['directions_by_state'].items():patterns[(contrast,state)][e['modality']]={'counts':counts,'evidence_id':e['id']}
        else:groups[(contrast,population(e['state'])[0])].append(e)
    programs=[]
    for (contrast,state),evidence in groups.items():
        by=defaultdict(list)
        for e in evidence:by[e['modality']].append(e)
        features={m:[e for e in records if e['kind']=='differential_feature'] for m,records in by.items()}
        paths=[e for e in evidence if e['kind']=='pathway']
        paths.sort(key=lambda e:(e['facts'].get('namespace')!='biological_process',e['facts'].get('fdr',1)))
        balanced=[]
        for direction in ('down','up','mixed/unspecified'):
            balanced.extend([e for e in paths if population(e['state'])[1]==direction][:2])
        paths=balanced[:4]
        rna={e['feature']:e for e in features.get('rna',[])};tf={e['feature']:e for e in features.get('grn_tf',[])}
        chains=[]
        for e in features.get('grn',[]):
            genes=edge_genes(e['feature'])
            if not genes:continue
            regulator,target=genes;target_call=rna.get(target);tf_call=tf.get(regulator)
            if not target_call:continue
            ids=[e['id'],target_call['id']]+([tf_call['id']] if tf_call else [])
            pathway=next((p for p in paths if target in re.split(r'[,;\s]+',str(p['facts'].get('overlap_genes','')))),None)
            if pathway:ids.append(pathway['id'])
            chains.append({'regulator':regulator,'target':target,'edge_change':change(e),'target_RNA_change':change(target_call),
                'TF_prediction_change':change(tf_call) if tf_call else 'no saved differential call',
                'same_direction':change(e)==change(target_call) and (not tf_call or change(tf_call)==change(e)),
                'pathway_overlap':pathway['facts']['term_name'] if pathway else '', 'edge_fdr':e['facts']['fdr'], 'evidence_ids':ids})
        chains.sort(key=lambda c:(c['TF_prediction_change']=='no saved differential call',c['edge_fdr'],-len(c['evidence_ids']),not c['same_direction']))
        # Each regulator gets a chance; many targets of one TF do not monopolize the narrative.
        distinct=[];seen=set()
        for chain in chains:
            if chain['regulator'] not in seen:distinct.append(chain);seen.add(chain['regulator'])
        chains=distinct[:3]
        supporting=[];parts=[]
        if paths:
            terms=[]
            for p in paths[:3]:
                _,direction=population(p['state'])
                terms.append(f'{p["facts"].get("term_name")} ({direction}-regulated gene list)')
                supporting.append(p['id'])
            parts.append('GO-Elite enrichment: '+ '; '.join(terms)+'.')
        if chains:
            for c in chains[:2]:
                parts.append(f'{c["regulator"]}→{c["target"]}: predicted edge score {c["edge_change"]}, target RNA {c["target_RNA_change"]}; predicted TF activity {c["TF_prediction_change"]}. This is directional association, not demonstrated activation or repression.')
                supporting.extend(c['evidence_ids'])
        for mod,label in [('rna','RNA'),('adt','predicted surface proteins'),('lipids','predicted lipids'),('grn_tf','predicted TF activities')]:
            records=features.get(mod,[])[:3]
            if records:
                parts.append(label+': '+', '.join(f'{e["feature"]} {change(e)}' for e in records)+'.')
                supporting.extend(e['id'] for e in records)
        communication=[e for e in evidence if e['kind']=='communication_contrast'][:2]
        for e in communication:
            parts.append(f'Predicted communication {e["feature"]}: score {change(e)} in the case arm. '+
                ('A saved statistical test is recorded.' if e['facts'].get('test_applied') else 'Descriptive contrast only; no statistical test was applied.'))
            supporting.append(e['id'])
        insights=[]
        sterol=next((p for p in paths if population(p['state'])[1]=='down' and
                     'cholesterol biosynthetic' in str(p['facts'].get('term_name','')).lower()),None)
        ester=next((e for e in features.get('lipids',[]) if e['feature'].startswith('CE(') and change(e)=='higher'),None)
        if sterol and ester:
            insights.append({'text':f'{state}: the down-regulated RNA list is enriched for cholesterol biosynthesis, while predicted {ester["feature"]} is higher. This discordance nominates sterol synthesis versus ester remodeling for direct measurement; it does not establish cholesterol storage, metabolic flux or model accuracy.',
                             'evidence_ids':[sterol['id'],ester['id']]})
        if len(communication)>1 and change(communication[0])!=change(communication[1]):
            insights.append({'text':f'{state}: predicted intercellular inputs change in opposite directions ({communication[0]["feature"]} {change(communication[0])}; {communication[1]["feature"]} {change(communication[1])}). Test a change in the balance of incoming signals rather than a uniform signaling increase. These expression-based scores do not establish physical contact or a spatial niche.',
                             'evidence_ids':[e['id'] for e in communication]})
        if not parts:continue
        c=next(comparisons[e['run_id']] for e in evidence)
        composition=[]
        layer=c['config'].get('population_col','');field=c['config'].get('sample_field','')
        source=digest['observations']['layers'].get(layer,{})
        for row in source.get('composition',{}).get(field,[]):
            if row['state']==state and row['value'] in c['config'].get('group1_samples',[])+c['config'].get('group2_samples',[]):composition.append(row)
        modalities=sorted({e['modality'] for e in evidence if e['kind'] in ('differential_feature','communication_contrast')})
        if state not in digest['observations']['layers'].get(layer,{}).get('counts',{}):
            # Coregulated/global pathway groups are not invented cell states.
            continue
        programs.append({'state':state,'layer':layer,'contrast':contrast,'comparison':c['comparison'],'modalities':modalities,
            'text':f'{state} — '+ ' '.join(parts),'regulatory_chains':chains,'composition':composition,
            'direction_counts':patterns[(contrast,state)],'evidence_ids':list(dict.fromkeys(supporting)),
            'ranking':{'covered_modalities':len(modalities),'pathways':len(paths),'RNA_supported_regulators':len(chains)},
            'pathway_ids':[e['id'] for e in paths]})
        programs[-1]['integrative_insights']=insights
    represented={(p['contrast'],p['state']) for p in programs}
    for c in comparisons.values():
        layer=c['config'].get('population_col','');field=c['config'].get('sample_field','')
        source=digest['observations']['layers'].get(layer,{})
        selected=c['config'].get('group1_samples',[])+c['config'].get('group2_samples',[])
        for state in source.get('counts',{}):
            if (c['contrast'],state) in represented:continue
            programs.append({'state':state,'layer':layer,'contrast':c['contrast'],'comparison':c['comparison'],
                'modalities':[],'text':f'{state}: no retained differential-feature or pathway evidence supports a program for this contrast. This does not establish absence of an effect.',
                'regulatory_chains':[],'composition':[row for row in source.get('composition',{}).get(field,[]) if row['state']==state and row['value'] in selected],
                'direction_counts':patterns[(c['contrast'],state)],'evidence_ids':[],'pathway_ids':[],'integrative_insights':[],
                'ranking':{'covered_modalities':0,'pathways':0,'RNA_supported_regulators':0}})
            represented.add((c['contrast'],state))
    programs.sort(key=lambda p:(-p['ranking']['covered_modalities'],-p['ranking']['RNA_supported_regulators'],-p['ranking']['pathways'],p['state']))
    return {'programs':programs,'design_limits':limits,'method':'Descriptive evidence assembly. Rank by modality coverage, RNA-supported regulator chains and saved pathway support. No new statistical significance; RNA-derived modalities are dependent evidence.'}


def narrative_from_synthesis(digest):
    evidence={e['id']:e for e in digest['evidence']};programs=digest['synthesis']['programs'];hypotheses=[]
    overview=[{'text':p['text'],'evidence_ids':p['evidence_ids'],'origin':'computed cross-modality synthesis'} for p in programs[:4]]
    overview.extend(dict(i,origin='computed integrative hypothesis') for p in programs[:4] for i in p['integrative_insights'])
    for modality in digest['modalities']:
        for p in [p for p in programs if modality in p['modalities']][:2]:
            refs=[e for e in (evidence[i] for i in p['evidence_ids']) if e['modality']==modality and e['kind'] in ('differential_feature','communication_contrast')]
            if not refs:continue
            paths=[evidence[i] for i in p['pathway_ids']]
            context=paths[0]['facts'].get('term_name') if paths else 'the recorded state phenotype'
            targets=', '.join(f'{e["feature"]} ({change(e)})' for e in refs[:2])
            if modality=='cell_communication':
                title=f'Test altered intercellular input to {p["state"]}: {refs[0]["feature"]}'
                validation='Repeat across biological samples, measure ligand and receptor proteins, and test receptor blockade with a downstream response assay. Expression compatibility alone does not demonstrate a spatial niche.'
            elif modality in ('grn','grn_tf'):
                title=f'Test the regulatory contribution to {context} in {p["state"]}'
                validation='Perturb the candidate TF and measure the linked target transcripts and relevant cellular response; verify the state assignment and replicate the contrast across biological samples.'
            elif modality=='lipids':
                title=f'Test whether lipid remodeling accompanies {context} in {p["state"]}'
                validation='Directly measure the named lipid species in independently collected, matched cell populations. RNA-derived lipid predictions do not measure lipid abundance or metabolic flux.'
            elif modality=='adt':
                title=f'Test the surface phenotype accompanying {context} in {p["state"]}'
                validation='Measure the named proteins by flow cytometry or matched protein profiling in independent samples; check state specificity and RNA/protein discordance.'
            else:
                title=f'Test {context} in {p["state"]}'
                validation='Replicate the state-specific transcript changes across biological samples and use a functional assay appropriate to the named process.'
            ids=[e['id'] for e in refs]+p['pathway_ids'][:2]
            for chain in p['regulatory_chains'][:1]:ids+=chain['evidence_ids']
            link_context=f'Consider these alongside the saved {context} enrichment in this state; they are derived from the same RNA experiment, not independent validation.' if paths else 'No matching saved pathway enrichment is available for this state.'
            if p['regulatory_chains']:
                chain=p['regulatory_chains'][0]
                link_context+=f' The matched {chain["regulator"]}→{chain["target"]} prediction has {chain["target_RNA_change"]} target RNA in the same contrast.'
            hypotheses.append({'modality':modality,'title':title,'rationale':f'{p["comparison"]}: {targets}. '+link_context,
                'validation':validation,'evidence_ids':list(dict.fromkeys(ids)),'origin':'computed evidence-linked hypothesis; requires experimental validation'})
    from .biology_questions import overview as question_overview
    return {'overview':question_overview(digest) or overview,'hypotheses':hypotheses}
