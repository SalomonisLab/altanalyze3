"""Section synthesis uses native records and cannot borrow identities or secrets."""
import json
import pytest
from altanalyze3.components.cellHarmony.webapp import biology_interpretation as section


def digest():
    hub={'factor':'STAT3','targets':['COL1A1'],'target_RNA_up':1,'target_RNA_down':0,'saved_pathway_overlaps':[{'term':'matrix organization','genes':['COL1A1']}],'TF_log2fc':1}
    return {'species':'human','reference':'test','observations':{'layers':{'state':{'counts':{'AM':20,'AT2':20}}}},
        'evidence':[{'kind':'differential_feature','feature':'STAT3','state':'AM','modality':'rna'},
                    {'kind':'differential_feature','feature':'SFTPC','state':'AT2','modality':'rna'}],
        'comparisons':[],'synthesis':{'design_limits':['one library per arm']},
        'regulatory_review':{'rows':[dict(state='AM',layer='state',comparison='Case vs Control',run_id='r',grn_run_id='g',tf_run_id='t',hubs=[hub])]}}


def response(monkeypatch, text, refs=None, experiment=None):
    monkeypatch.setenv('SCALABLE_SUMMARY_URL','http://localhost/test')
    monkeypatch.setenv('SCALABLE_SUMMARY_MODEL','fixture')
    reply={'choices':[{'message':{'content':json.dumps({'sections':{'grn':[{'implication':text,'alternative':'A change in captured cell identity could account for the association.','experiment':experiment or 'Assess the named targets in independent libraries to distinguish identity from a functional response.','source_ids':refs or ['G1']}]}})}}]}
    class Reply:
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def read(self,*args):return json.dumps(reply).encode()
    monkeypatch.setattr(section.urllib.request,'urlopen',lambda *a,**k:Reply())


def test_interpretation_can_explain_consequences_as_hypotheses(monkeypatch):
    text='Could the predicted matrix-remodeling program in AM affect tissue repair rather than reflect a change in cell identity?'
    response(monkeypatch,text)
    r=section.interpret_sections(digest())
    assert text in r['sections']['grn'][0]['text'] and r['sources']['G1']['hubs'][0]['factor']=='STAT3'


@pytest.mark.parametrize('text,refs',[
    ('Could the predicted STAT3 network in AT2 is consistent with matrix remodeling; measure COL1A1 in independent libraries?',['G1']),
    ('Could the predicted SFTPC program in AM could reflect matrix remodeling; independent libraries should distinguish the alternatives.',['G1']),
    ('Could the predicted STAT3 factor represses COL1A1 in AM and could therefore modify matrix remodeling in this library comparison?',['G1']),
    ('Could the predicted STAT3 program in AM could affect matrix remodeling, requiring independent validation against composition effects?',['G999'])])
def test_rejects_borrowed_state_feature_signed_regulation_and_unbound_sources(monkeypatch,text,refs):
    response(monkeypatch,text,refs)
    r=section.interpret_sections(digest())
    assert r['status']=='rejected' and r['rejected_sections'] and not r['sections']['grn']


def test_no_model_does_not_fabricate_narrative(monkeypatch):
    monkeypatch.delenv('SCALABLE_SUMMARY_URL',raising=False)
    assert section.interpret_sections(digest())['status']=='unavailable'


def test_model_validation_suggestion_is_audited_but_not_promoted(monkeypatch):
    text='Could the predicted matrix-remodeling program in AM affect tissue repair rather than reflect a change in cell identity?'
    invalid='Validate RNA-derived lipid predictions with independent RNA-seq.'
    response(monkeypatch,text,experiment=invalid)
    row=section.interpret_sections(digest())['sections']['grn'][0]
    assert invalid not in row['text'] and row['model_suggestions']['experiment']==invalid
    assert 'matrix deposition' in row['functional_checks'][0]['functional_readout']
    assert 'code-defined' in row['text']


def test_rejects_unperformed_significance_test(monkeypatch):
    response(monkeypatch,'Could the predicted STAT3 network demonstrate a significant matrix-remodeling response in AM?')
    result=section.interpret_sections(digest())
    assert result['status']=='rejected'


def test_rejects_model_replacing_the_canonical_state_with_a_guess(monkeypatch):
    response(monkeypatch,'Could the predicted matrix-remodeling program in airway macrophages affect tissue repair?')
    result=section.interpret_sections(digest())
    assert result['status']=='rejected'
    assert 'canonical cell-state identity' in result['rejected_sections'][0]['reason']


def test_section_prompts_are_separate_and_keep_native_directions(monkeypatch):
    d=digest()
    d['comparisons']=[{'id':'r','comparison':'Case vs Control'}]
    d['evidence'].append({'kind':'pathway','state':'AM__down','layer':'state','run_id':'r','id':'E1','feature':'',
        'facts':{'term_name':'matrix organization','overlap_genes':'COL1A1','fdr':0.01},'title':'matrix organization'})
    monkeypatch.setenv('SCALABLE_SUMMARY_URL','http://localhost/test')
    monkeypatch.setenv('SCALABLE_SUMMARY_MODEL','fixture')
    sent=[]
    class Reply:
        def __init__(self,body):self.body=body
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def read(self,*args):return json.dumps(self.body).encode()
    def call(request,**kwargs):
        payload=json.loads(request.data)
        data=json.loads(payload['messages'][1]['content']);sent.append(data)
        section_name=data['section']
        row={'implication':'Could the predicted matrix-remodeling program in AM reflect altered tissue repair?',
             'alternative':'Different cell identities may account for this association.',
             'experiment':'Measure matrix deposition with identity controls in independent donors.'}
        return Reply({'choices':[{'message':{'content':json.dumps({'sections':{section_name:[row]}})}}]})
    monkeypatch.setattr(section.urllib.request,'urlopen',call)
    result=section.interpret_sections(d)
    assert result['status']=='completed'
    assert [p['section'] for p in sent]==['grn','pathways']
    assert 'predicted_regulatory_programs' in sent[0]['findings'][0]
    assert sent[1]['findings'][0]['gene_list_direction']=='down'
    assert all('hubs' not in p['findings'][0] for p in sent)


def test_prompt_removes_geometry_without_losing_distinct_molecular_nodes():
    d=digest();node={'id':'glyph1','label':'Cholesterol','type':'Metabolite','x':10,'y':20,'link_feature':'CE(22:4)','log2fc':1}
    d['integrated_pathway_review']={'rows':[{'state':'AM','layer':'state','comparison':'Case vs Control','run_id':'r','modality':'lipid','id':'WP1','name':'Lipid pathway','nodes':[node,dict(node,id='glyph2',x=30),dict(node,id='glyph3',link_feature='CE(18:0)')]}]}
    payload,sources=section.context(d)
    assert len(payload['sections']['multimodal_pathways'][0]['nodes'])==2
    assert len(sources['M1']['nodes'])==3
    assert all('x' not in n for n in payload['sections']['multimodal_pathways'][0]['nodes'])


def test_state_name_prefix_is_not_a_borrowed_cell_state(monkeypatch):
    d=digest();d['observations']['layers']['state']['counts']['AM-prolif']=20
    d['regulatory_review']['rows'][0]['state']='AM-prolif'
    text='Could the predicted matrix-remodeling program in AM-prolif affect tissue repair rather than reflect a change in cell identity?'
    response(monkeypatch,text)
    r=section.interpret_sections(d)
    assert text in r['sections']['grn'][0]['text']
