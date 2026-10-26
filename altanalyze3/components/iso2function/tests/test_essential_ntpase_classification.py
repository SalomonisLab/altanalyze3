"""A broad nucleotide-binding fold must not confer a helicase-specific label."""
from altanalyze3.components.iso2function.plotting.protein_features import essential_rule,load_essential_whitelist

def test_p_loop_superfamily_keeps_generic_native_annotation_class():
    row=dict(provider='protein_features',name='SSF52540',accession='IPR027417',description='P-loop_NTPase')
    match=essential_rule(row,load_essential_whitelist())
    assert match['cls']=='NTPASE' and match['label']=='P-loop NTPase'

def test_helicase_specific_and_gtpase_specific_domains_remain_distinct():
    rules=load_essential_whitelist()
    helicase=dict(provider='protein_features',name='PF00270',accession='IPR011545',description='DEAD/DEAH box helicase')
    gtpase=dict(provider='protein_features',name='PF00071',accession='IPR001806',description='Small_GTPase')
    assert essential_rule(helicase,rules)['cls']=='HELICASE'
    assert essential_rule(gtpase,rules)['cls']=='GTPASE'

def test_gpi_attachment_is_visible_without_relabelling_other_lipidations():
    rules=load_essential_whitelist()
    gpi=dict(provider='uniprot_features',name='',accession='',description='LIPID-GPI-anchor amidated glycin')
    other=dict(provider='uniprot_features',name='',accession='',description='LIPID-N-myristoyl glycin')
    assert essential_rule(gpi,rules)['cls']=='GPI'
    assert essential_rule(other,rules) is None
