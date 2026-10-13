"""Explicit biological questions evaluated against native, completed results.

This is an interpretation contract, not a replacement inference/statistics
pipeline. Prediction agreement never counts as independent validation.
"""
from collections import defaultdict
import hashlib
import re

from .biology_synthesis import population, change

# These are question categories, not gene sets or asserted disease mechanisms.
# Membership is determined by the saved enriched GO term, not by the LLM.
DOMAINS = (
    ('host_defense', 'Host defense, inflammatory recruitment and antigen presentation',
     r'immune|inflamm|chemotaxis|antigen|phagocyt|defen[cs]e|cytokine',
     'Measure the implicated immune response or phagocytic function in matched cell populations.'),
    ('surfactant', 'Surfactant production and turnover', r'surfactant|lamellar|phospholipid',
     'Measure surfactant composition and secretion or turnover in the implicated cell population.'),
    ('lipid_metabolism', 'Lipid synthesis, handling and remodeling', r'lipid|sterol|cholesterol|fatty acid',
     'Measure the implicated lipid species and biosynthetic or uptake flux; discriminate synthesis from remodeling.'),
    ('matrix', 'Extracellular matrix production and remodeling', r'collagen|extracellular matrix|matrix organization',
     'Measure matrix deposition or remodeling and test the candidate regulator or receptor by perturbation.'),
    ('barrier', 'Vascular and epithelial barrier maintenance', r'junction|permeability|blood vessel|angiogen|endotheli',
     'Measure permeability or junction organization in matched cells and perturb the nominated regulator.'),
    ('airway_clearance', 'Mucociliary clearance and airway secretion', r'cili|mucus|mucin|secretory',
     'Measure ciliary transport or the implicated secretory function rather than inferring function from transcripts.'),
    ('development_repair', 'Differentiation, development and injury repair',
     r'differentiation|development|morphogenesis|proliferation|cell cycle|regeneration|wound',
     'Test differentiation or repair in a time course or perturbation; compare matched developmental stages where recorded.'),
    ('cellular_state', 'Biosynthetic demand and cellular stress',
     r'translation|ribosom|mitochond|oxidative|unfolded|apoptosis|hypoxi',
     'Measure the implicated biosynthetic or stress response and distinguish viable state adaptation from technical stress.'),
)


def evaluate(digest):
    evidence = {e['id']: e for e in digest['evidence']}
    comparisons = {c['id']: c for c in digest['comparisons']}
    candidates, shared, domains = [], defaultdict(list), defaultdict(list)
    matched = {r['state']: r for r in digest['cell_state_biology']['matched']}
    for p in digest['synthesis']['programs']:
        if not p['evidence_ids']:
            continue
        records = [evidence[i] for i in p['evidence_ids']]
        paths = [evidence[i] for i in p['pathway_ids']]
        features = [e for e in evidence.values() if e['kind']=='differential_feature'
                    and e['state']==p['state'] and e['layer']==p['layer']
                    and comparisons[e['run_id']]['contrast']==p['contrast']]
        rna = [e for e in features if e['modality'] == 'rna']
        predictions = [e for e in features if e['modality'] != 'rna']
        communication = [e for e in records if e['kind'] == 'communication_contrast']
        c = next(comparisons[e['run_id']] for e in records if e['run_id'])
        cfg = c['config']
        for path in paths:
            term = str(path['facts'].get('term_name', ''))
            domain = next((d for d in DOMAINS if re.search(d[2], term, re.I)),
                ('other_process', 'Other annotated cellular process', '',
                 'Use a functional assay specific to the recorded enriched process.'))
            _, direction = population(path['state'])
            genes = set(re.split(r'[,;\s]+', str(path['facts'].get('overlap_genes', '')))) - {''}
            targets = [e for e in rna if e['feature'] in genes]
            chains = [ch for ch in p['regulatory_chains'] if ch['target'] in genes]
            ids = [path['id']] + [e['id'] for e in targets]
            ids += [i for ch in chains for i in ch['evidence_ids']]
            # Same-state predictions are follow-ups; no feature-specific biochemical
            # link to this process is asserted solely from coexistence.
            followups = predictions + communication
            ids += [e['id'] for e in followups]
            ids = list(dict.fromkeys(ids))
            key = f"{p['contrast']}|{p['state']}|{term}|{direction}"
            candidate = {
                'id': 'Q' + hashlib.sha256(key.encode()).hexdigest()[:12],
                'domain': domain[0], 'domain_name': domain[1], 'state': p['state'],
                'layer': p['layer'], 'contrast': p['contrast'], 'comparison': p['comparison'],
                'question': f"Does {p['state']} show a covariate-associated difference in {term}, and which predicted molecular phenotypes warrant testing alongside it?",
                'observation': f"The {direction}-regulated gene list is enriched for {term}. This describes list enrichment, not measured pathway activity.",
                'pathway_direction': direction, 'pathway_evidence_id': path['id'],
                'RNA_anchor_ids': [e['id'] for e in targets],
                'regulatory_chains': chains,
                'predicted_followup_ids': [e['id'] for e in followups],
                'evidence_ids': ids,
                'evidence_level': 'RNA-list-anchored hypothesis; predicted modalities require direct validation',
                'background_source': matched.get(p['state'], {}).get('card', {}).get('url', ''),
                'established_state_function': matched.get(p['state'], {}).get('card', {}).get('facts', {}).get('Function', {}),
                'background_role': 'Established state function provides context, not evidence of a covariate effect or a predicted phenotype.',
                'alternatives': [
                    'Within-state biological change versus altered mixtures or annotation within this state.',
                    'RNA-derived prediction dependence or model error versus an actual molecular phenotype.',
                    'Covariate biology versus library, donor or technical differences not separated by this contrast.'],
                'disagreements': [ch for ch in chains if not ch['same_direction']],
                'integrative_followups': p['integrative_insights'],
                'validation': domain[3] + ' Replicate across independent biological samples; directly measure predicted modalities.',
                'comparison_unit': cfg.get('comparison_type', ''),
                'replication_established': False,
                'prioritization': {'pathway_FDR': path['facts'].get('fdr', 1),
                    'RNA_supported_regulator_chains': len(chains), 'named_RNA_support': len(targets)}
            }
            candidates.append(candidate)
            shared[(p['contrast'], p['layer'], term, direction)].append(candidate)
            domains[(p['contrast'],p['layer'],domain[0])].append(candidate)
    # A hierarchy, not a significance score or a vote among dependent modalities.
    candidates.sort(key=lambda q: (-bool(q['prioritization']['RNA_supported_regulator_chains']),
        -bool(q['prioritization']['named_RNA_support']), q['domain']=='other_process',
        q['prioritization']['pathway_FDR'], q['id']))
    for rank, candidate in enumerate(candidates, 1):
        candidate['priority'] = rank
    cross_state = [{'contrast': contrast, 'layer': layer, 'term': term,
                    'direction': direction, 'states': sorted({q['state'] for q in rows}),
                    'question_ids': [q['id'] for q in rows],
                    'interpretation': 'Shared enriched process in the same contrast; does not establish cell-to-cell coupling.'}
        for (contrast, layer, term, direction), rows in shared.items() if len({q['state'] for q in rows}) > 1]
    contrasts = defaultdict(list)
    for c in digest['comparisons']:
        contrasts[c['contrast']].append(c)
    plan = []
    for contrast, runs in contrasts.items():
        cfg = runs[0]['config']
        completed = sorted({c['modality'] for c in runs})
        plan.append({'contrast': contrast, 'layer': cfg.get('population_col'),
            'group_values_from': cfg.get('sample_field'),
            'case_values': cfg.get('group1_samples', []), 'control_values': cfg.get('group2_samples', []),
            'comparison_unit': cfg.get('comparison_type'), 'completed_modalities': completed,
            'missing_modality_comparisons': sorted(set(digest['modalities']) - set(completed)),
            'pathway_analysis': 'available' if any(e['kind']=='pathway' and e['run_id'] in {c['id'] for c in runs}
                                                    for e in digest['evidence']) else 'not available',
            'replication': 'Not established by selected group values or cell counts; verify biological sample identities.',
            'adjusted_confounder_analysis': 'Not inferred from an unadjusted saved comparison.'})
    framework = {
        'version': 1,
        'role': 'Code defines questions, eligible contrasts, evidence joins and priority. Code authors the propositions. The model selects experimental priorities from defined alternatives; it cannot create analyses or evidence.',
        'catalogue': [{'id': d[0], 'name': d[1]} for d in DOMAINS],
        'analysis_contract': [
            'Compare the same state layer within explicitly selected covariate arms; retain the original analysis unit and parameters.',
            'Anchor functional questions in saved RNA differential results and GO-Elite up/down enrichment; marker identity enrichment is separate.',
            'Link regulators only through named target RNA in the same state/contrast; link to a pathway only through its recorded overlap genes.',
            'Treat surface proteins, lipids, metabolites and regulatory predictions as testable molecular phenotypes; their common RNA origin is not replication.',
            'Communication contrasts without statistical tests are descriptive. Score compatibility does not establish signaling, colocalization or a niche.',
            'Evaluate concordance and discordance, composition and competing state/technical explanations; specify a distinguishing experiment.',
            'Use covariate biology only for recorded/completed comparisons; disease, age and sex labels do not automatically establish causation or adjustment.',
            'Prioritize RNA-supported regulator-target-pathway chains, named RNA support, the explicit physiological question catalogue, then saved pathway FDR; this is an editorial hierarchy, not a new significance score or a vote among predicted modalities.'],
        'analysis_plan': plan, 'questions': candidates, 'cross_state_processes': cross_state,
        'cross_state_domain_questions': [{'contrast':contrast,'layer':layer,'domain':domain,
            'question':'Do related transcript programs show different directions across these cell states rather than a uniform tissue-wide response?',
            'state_processes':[{'state':q['state'],'observation':q['observation'],'direction':q['pathway_direction'],
                               'question_id':q['id'],'evidence_id':q['pathway_evidence_id']} for q in rows],
            'opposing_directions':len({q['pathway_direction'] for q in rows})>1,
            'interpretation':'Shared question category, not an identical pathway/activity score or demonstrated cross-cell coupling.'}
            for (contrast,layer,domain),rows in domains.items() if domain!='other_process' and
            len({q['state'] for q in rows})>1],
        'state_coverage': [{'state':p['state'],'layer':p['layer'],'contrast':p['contrast'],
            'RNA_functional_question_available':any(q['state']==p['state'] and q['contrast']==p['contrast'] for q in candidates),
            'available_modalities':p['modalities'],
            'interpretation':'No RNA-anchored functional question is not evidence of no biological change.'}
            for p in digest['synthesis']['programs']],
        'design_limits': digest['synthesis']['design_limits'],
        'coverage': 'Questions evaluate retained summary evidence, not an exhaustive pathway screen. Full native analyses remain downloadable.',
        'missing_requirements': ([] if plan else ['No completed covariate comparison: evaluate state identity/markers only; no disease/development effect can be inferred.'])
    }
    tasks=[]
    for modality in digest['modalities']:
        eligible=[q for q in candidates if any(evidence[i]['modality']==modality and
            evidence[i]['kind'] in ('differential_feature','communication_contrast') for i in q['evidence_ids'])]
        if modality=='lipids':eligible.sort(key=lambda q:(q['domain']!='lipid_metabolism',q['priority']))
        if not eligible:continue
        q=eligible[0]
        records=[evidence[i] for i in q['evidence_ids'] if evidence[i]['modality']==modality and
                 evidence[i]['kind'] in ('differential_feature','communication_contrast')]
        refs=[q['pathway_evidence_id']]+q['RNA_anchor_ids'][:2]
        observations=[q['observation']]
        if modality in ('grn','grn_tf') and q['regulatory_chains']:
            ch=q['regulatory_chains'][0];refs+=ch['evidence_ids']
            observations.append(f'{ch["regulator"]}→{ch["target"]}: predicted edge {ch["edge_change"]}, target RNA {ch["target_RNA_change"]}, predicted TF activity {ch["TF_prediction_change"]}. Directional association, not demonstrated activation/repression.')
        refs+=[e['id'] for e in records[:2]]
        observations += [f'{e["feature"]}: {change(e)} in the case arm.' for e in records[:2]]
        if modality=='lipids':
            for insight in q['integrative_followups']:
                if any(evidence[i]['modality']=='lipids' for i in insight['evidence_ids']):
                    observations.append(insight['text']);refs+=insight['evidence_ids']
        if modality=='cell_communication':
            observations.append('These saved contrasts have no statistical test unless explicitly recorded; investigate incoming-signal balance, not a proven niche.')
        validation=q['validation']
        if modality=='adt':validation='Measure the named surface proteins in matched cell populations; distinguish a true protein phenotype from prediction error and within-state mixtures.'
        if modality in ('grn','grn_tf'):validation='Perturb the named regulator, measure its target transcripts and the relevant cell function in independent biological samples.'
        if modality=='cell_communication':validation='Measure the named ligand and receptor proteins, establish proximity and test receptor blockade with a downstream-response assay in independent samples.'
        if modality=='rna':
            logic='Test whether the recorded transcript differences reflect a distinct functional state versus altered mixtures within this annotation. Gene-list enrichment does not measure the function itself.'
        elif modality=='adt':
            logic='Test whether the predicted surface phenotype is present despite the RNA-program direction, versus a prediction or annotation-mixture explanation. A surface-marker prediction cannot establish increased antigen presentation or inflammatory recruitment.'
        elif modality=='lipids':
            logic='Test the supplied RNA/lipid concordance or discordance by direct molecular measurement. Distinguish synthesis, uptake/remodeling and prediction error; do not assert lipid storage or flux from predictions.'
        elif modality in ('grn','grn_tf'):
            logic='Test the named regulator-target association by perturbation. Same-direction TF/edge/RNA predictions nominate a coupled program; they do not identify activation versus repression, a feedback loop, or a causal regulator.'
        elif modality=='cell_communication':
            logic='Test whether opposite incoming expression-compatibility score changes correspond to different functional inputs. A score change does not demonstrate ligand expression change, signaling, physical proximity or a spatial niche.'
        else:
            logic='Test the named predicted molecular phenotype by direct measurement and distinguish a biological association from RNA-derived model error.'
        tasks.append({'task_id':f'H{len(tasks)+1}','question_id':q['id'],'modality':modality,
            'title':f'{q["state"]}: {q["domain_name"]} — {modality}',
            'state':q['state'],'comparison':q['comparison'],'observations':observations,
            'alternatives':q['alternatives'],'validation':validation,
            'logical_proposition':logic,
            'experiment_options':[validation,
                'Replicate this exact state/contrast in independently collected biological samples before assigning a covariate effect.',
                'Check annotation purity and within-state mixtures using the recorded RNA markers, then directly measure the named molecular phenotype.'],
            'evidence_ids':list(dict.fromkeys(refs))})
    framework['hypothesis_tasks']=tasks
    return framework


def overview(digest):
    """Lead with relationships and distinguishing questions, not feature lists."""
    f=digest.get('question_framework',{});evidence={e['id']:e for e in digest['evidence']};rows=[]
    labels={c['contrast']:c['comparison'] for c in digest['comparisons']}
    for q in f.get('cross_state_domain_questions',[]):
        if not q['opposing_directions']:continue
        parts=[f'{p["state"]}: {evidence[p["evidence_id"]]["facts"]["term_name"]} ({p["direction"]}-regulated RNA list)'
               for p in q['state_processes']]
        rows.append({'text':labels.get(q['contrast'],'Recorded comparison')+': '+ '; '.join(parts)+
            '. Test a compartment-specific response versus differences in within-state mixtures. Related process categories do not establish a common pathway activity or cell-to-cell coupling.',
            'evidence_ids':list(dict.fromkeys(p['evidence_id'] for p in q['state_processes'])),
            'origin':'code-defined cross-state question'})
    seen=set()
    for q in f.get('questions',[]):
        for insight in q['integrative_followups']:
            key=tuple(insight['evidence_ids'])
            if key not in seen:
                rows.append(dict(insight,origin='code-defined integrative question'));seen.add(key)
    for q in f.get('questions',[]):
        if not q['regulatory_chains']:continue
        c=q['regulatory_chains'][0];key=tuple(c['evidence_ids'])
        if key in seen:continue
        seen.add(key)
        rows.append({'text':f'{q["state"]} ({q["comparison"]}): {c["regulator"]}→{c["target"]} has predicted edge {c["edge_change"]}, target RNA {c["target_RNA_change"]} and predicted TF activity {c["TF_prediction_change"]}; the target occurs in the saved {evidence[q["pathway_evidence_id"]]["facts"]["term_name"]} enrichment. Perturb the regulator and measure the target and relevant function to distinguish regulatory involvement from co-associated predictions. Directional agreement does not identify activation/repression or prove causation.',
            'evidence_ids':list(dict.fromkeys(c['evidence_ids']+[q['pathway_evidence_id']])),
            'origin':'code-defined regulator-target-process question'})
        if len(rows)>=6:break
    return rows[:6]
