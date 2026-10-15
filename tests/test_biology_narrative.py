"""Integrative interpretations retain actual model prose and bind saved sources."""
import copy
import json
import pytest
from altanalyze3.components.cellHarmony.webapp import biology_narrative as narrative
from altanalyze3.components.cellHarmony.webapp import biology_literature,summary_runtime


def digest():
    return {'species':'human','reference':'test','modalities':['rna','grn'],
      'comparisons':[{'id':'r','contrast':'case','comparison':'Case vs Control'}],
      'observations':{'layers':{'supervised_state':{'counts':{'AM':20,'AT2':20}}}},
      'synthesis':{'design_limits':['one library per arm']},
      'evidence':[{'id':'E1','run_id':'r','layer':'supervised_state','state':'AM','kind':'differential_feature','modality':'rna','feature':'STAT3','facts':{'log2fc':1}},
                  {'id':'E2','run_id':'r','layer':'supervised_state','state':'AT2','kind':'differential_feature','modality':'rna','feature':'SFTPC','facts':{'log2fc':-1}}],
      'regulatory_review':{'rows':[{'run_id':'r','layer':'supervised_state','state':'AM','tf_activity_status':'No matched TF activity differential',
           'hubs':[{'factor':'STAT3','targets':['COL1A1'],'saved_pathway_overlaps':[],'target_RNA_up':1,'target_RNA_down':0}]}]}}


def stub(monkeypatch,rows):
    monkeypatch.setenv('SCALABLE_SUMMARY_URL','http://localhost/fixture')
    monkeypatch.setenv('SCALABLE_SUMMARY_MODEL','fixture')
    monkeypatch.setattr(biology_literature,'review',lambda d:{'literature_candidates':[{'state':'AM','feature':'STAT3','comparison':'Case vs Control','literature':{'papers':[{'pmid':'123','title':'STAT3 and repair','url':'https://pubmed.ncbi.nlm.nih.gov/123/'}]}}]})
    monkeypatch.setattr(narrative,'abstracts',lambda pmids:{'123':'Related repair background, not this dataset.'})
    captured=[]
    def completion(request,*args):
        captured.append(json.loads(request.data))
        return {'model':'actual-fixture-model','choices':[{'finish_reason':'stop','message':{'content':json.dumps({'sections':rows})}}]}
    monkeypatch.setattr(summary_runtime,'completion',completion)
    return captured


def row(**kw):
    return dict({'kind':'hypothesis','source_ids':['S1'],'citation_ids':['L1'],
        'text':'Could the STAT3-associated program in AM modify matrix remodeling and tissue repair?',
        'alternative':'Differences in captured cellular states could explain this association.',
        'experiment':'Measure matrix deposition after perturbation in independent donor cultures.'},**kw)


def test_model_prose_and_real_citation_are_retained(monkeypatch):
    item=row();calls=stub(monkeypatch,[item]);r=narrative.synthesize(digest())
    assert r['status']=='completed' and r['sections'][0]['text']==item['text']
    assert r['citations']['L1']['url']=='https://pubmed.ncbi.nlm.nih.gov/123/'
    assert r['model']=='actual-fixture-model' and len(calls)==1
    assert 'No matched TF' in calls[0]['messages'][1]['content']
    assert 'abstract_excerpt' in calls[0]['messages'][1]['content']


@pytest.mark.parametrize('item',[
    row(source_ids=['S99']),row(citation_ids=['L99']),
    row(text='Could the STAT3 program in AT2 modify tissue repair during the comparison?'),
    row(text='Could the SFTPC program in AM modify tissue repair during the comparison?'),
    row(text='Could the STAT3 program in AM explain 80% of cells in this comparison?'),
    row(experiment='None')])
def test_invalid_prose_is_explicitly_rejected_without_static_replacement(monkeypatch,item):
    stub(monkeypatch,[item]);r=narrative.synthesize(digest())
    assert r['status']=='rejected' and not r['sections'] and r['rejected_sections']


def test_focus_preserves_regulatory_context_without_inventing_tf_measurements():
    d=digest();e,s=narrative.packet(d,'rna')
    assert len(s)==2 and s['S1']['hubs'][0]['factor']=='STAT3'
    assert s['S1']['TF_activity_evidence']=='No matched TF activity differential'
    assert 'TF_log2fc' in s['S1']['hubs'][0] and s['S1']['hubs'][0]['TF_log2fc'] is None


def test_unconfigured_model_is_explicit(monkeypatch):
    stub(monkeypatch,[row()]);monkeypatch.delenv('SCALABLE_SUMMARY_URL')
    r=narrative.synthesize(digest())
    assert r['status']=='unavailable' and not r['sections'] and r['citations']


def test_mixed_layers_are_not_integrated(monkeypatch):
    stub(monkeypatch,[row(source_ids=['S1','S2'])]);d=digest()
    d['observations']['layers']['unsupervised_cluster']={'counts':{'AT2':20}}
    d['evidence'][1]['layer']='unsupervised_cluster'
    r=narrative.synthesize(d)
    assert r['status']=='rejected' and 'Mixed cell-state' in r['rejected_sections'][0]['reason']


def test_literature_ids_cannot_be_used_for_unrelated_paragraphs(monkeypatch):
    stub(monkeypatch,[row(text='Could the COL1A1 program in AM modify matrix remodeling during tissue repair?')])
    r=narrative.synthesize(digest())
    assert r['status']=='rejected' and 'Literature citation' in r['rejected_sections'][0]['reason']


def test_shared_rna_claim_requires_records_in_each_state(monkeypatch):
    stub(monkeypatch,[row(source_ids=['S1','S2'],citation_ids=[],
         text='Could a shared increase of STAT3 and SFTPC in both AM and AT2 modify tissue repair?')])
    r=narrative.synthesize(digest())
    assert r['status']=='rejected' and 'Shared RNA assertion' in r['rejected_sections'][0]['reason']


def test_detail_dispatches_each_requested_focus_without_preloading_other_sections(monkeypatch):
    from altanalyze3.components.cellHarmony.webapp.biology_summary import build_summary_detail
    calls=[]
    monkeypatch.setattr(narrative,'synthesize',lambda d,focus:calls.append(focus) or {'sections':[],'literature_review':{'checked':True}})
    surface={'digest':digest()}
    assert build_summary_detail(surface,'synthesis')['kind']=='summary_synthesis'
    assert calls==['systems']
    build_summary_detail(surface,'synthesis','rna');assert calls==['systems','rna']
    assert build_summary_detail(surface,'literature')['review']=={'checked':True}
    assert calls==['systems','rna','literature']


def test_gene_suffix_is_not_a_borrowed_b_cell_state(monkeypatch):
    stub(monkeypatch,[row(text='Could HLA-B-associated changes in AM alter antigen presentation during tissue repair?',citation_ids=[])])
    d=digest();d['observations']['layers']['supervised_state']['counts']['B']=10
    d['evidence'][0]['feature']='HLA-B'
    r=narrative.synthesize(d)
    assert r['status']=='completed'


def test_literature_focus_requests_source_bound_interpretation(monkeypatch):
    calls=stub(monkeypatch,[row(kind='literature',text='AM: STAT3-related repair literature provides a functional context for this saved association.')])
    result=narrative.synthesize(digest(),'literature')
    assert result['status']=='completed'
    context=json.loads(calls[0]['messages'][1]['content'])
    assert context['literature_questions'][0]['source_ids']==['S1']
    assert set(context['sources'])=={'S1'}
    assert 'AT2' not in json.dumps(context) and 'SFTPC' not in json.dumps(context)
    assert context['literature_questions'][0]['required_text_prefix']=='AM: '
    assert calls[0]['response_format']['schema']['properties']['sections']['minItems']==2


def test_literature_symbol_case_is_preserved_without_losing_its_source_binding(monkeypatch):
    stub(monkeypatch,[row(kind='literature',text='AM: Stat3 was studied in a distinct experimental repair model, providing a testable context for this association.')])
    result=narrative.synthesize(digest(),'literature')
    assert result['status']=='completed' and 'Stat3' in result['sections'][0]['text']
    assert result['sources']['S1']['features'][0]=='STAT3'


def test_literature_cannot_present_a_dataset_association_as_functional_validation(monkeypatch):
    stub(monkeypatch,[row(kind='literature',text='AM: STAT3 is lower in Case, which limits inflammation and tissue injury in these cells.')])
    result=narrative.synthesize(digest(),'literature')
    assert result['status']=='rejected' and 'established functional consequence' in result['rejected_sections'][0]['reason']
