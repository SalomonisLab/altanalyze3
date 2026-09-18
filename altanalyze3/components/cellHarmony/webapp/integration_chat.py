"""Local routing for integrated views; statistics always come from stored analyses."""
import re
from altanalyze3.components.cellHarmony import discover_integration as views
from altanalyze3.components.cellHarmony.grn_analysis import read_regulatory_question
from .integration_data import integration_data


def read_question(question, states, genes):
    q=question.lower()
    pathway=bool(re.search(r'\b(pathway|pathways)\b',q)) and not re.search(r'\b(go|goelite|go-elite|ontology)\b',q)
    network=bool(re.search(r'\b(network|networks|grn|regulatory|regulators)\b',q)) and not re.search(r'\b(communication|ligand|receptor)\b',q)
    if not pathway and not network:return None
    if network and re.search(r'\b(differential|changed|significant|upregulated|downregulated)\b',q) and not re.search(r'\b(network|networks|impacts)\b',q):return None
    if network and re.search(r'\b(activity|activities)\b',q) and not re.search(r'\b(network|networks|edges)\b',q):return None
    base=read_regulatory_question('Show regulatory network '+question,states,genes) or {}
    return dict(base,intent='integrated_pathway' if pathway else 'integrated_network',
                modality='metabolite' if 'metabol' in q else 'lipid' if pathway else 'grn',
                source='marker' if 'marker' in q else 'differential',
                limit=base.get('limit',50) if re.search(r'(?:top|show)\s+\d+',q) else 50)


def answer(app,meta,question,reading):
    ds=integration_data(app,meta)
    state=reading.get('cell_state') or ''
    contrast=reading.get('contrast') or ds.current_contrast
    # Resolve a named comparison within this dataset; never substitute an unrelated run.
    entries=ds.ds.deg_manifest().get('comparisons',[]) if hasattr(ds,'ds') else ds.deg_manifest().get('comparisons',[])
    norm=lambda s:re.sub(r'[^a-z0-9]+',' ',s.lower()).strip()
    named=next((c['id'] for c in entries if c.get('comparison') and norm(c['comparison']) in norm(question)),None)
    contrast=named or contrast
    if reading['modality']=='lipid' and not ds.available(contrast,'lipid',state) and ds.available(contrast,'lipids',state):
        reading=dict(reading,modality='lipids')
    spec=dict(kind=reading['intent'],cell_state=state,contrast=contrast,modality=reading['modality'],
              source=reading['source'],features=reading.get('genes') or [],limit=int(reading.get('limit') or 50))
    if not state:
        return {'question':question,'reading':reading,'intent':reading['intent'],'answer':'Choose a cell state for the integrated view.',
                'choices':{'states':ds.states},'status':'clarify'}
    try:
        if reading['intent']=='integrated_network':
            result=views.network(ds,cell_state=state,contrast=contrast,features=spec['features'],source=spec['source'],limit=spec['limit'])
        else:
            result=views.pathways(ds,cell_state=state,contrast=contrast,modality=spec['modality'],source=spec['source'],features=spec['features'])
    except (ValueError,KeyError,FileNotFoundError) as exc:
        result={'available':False,'note':str(exc)}
    return dict(result,question=question,reading=dict(reading,contrast=contrast),intent=reading['intent'],plot=spec,
                answer=result.get('note') or ('Select a pathway to see the matching gene and metabolite differentials.'
                if reading['intent']=='integrated_pathway' else 'Regulatory network from matching expression and edge differentials.'))
