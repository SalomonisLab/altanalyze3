import pytest
from altanalyze3.components.annotation.splice_survey_report import check_expectations,html_table


def test_literature_match_fails_on_changed_native_coding_signature():
    expected=dict(aa_1=512,aa_2=491,unmapped_aa_1=[[23,43]])
    check_expectations(dict(expected),expected,'LYN')
    with pytest.raises(ValueError,match='LYN aa_2'):check_expectations(dict(expected,aa_2=490),expected,'LYN')
    assert '&lt;native&gt;' in html_table([{'isoform':'<native>'}],['isoform'])


def test_missing_gene_mapping_is_resolved_from_source_instead_of_zero_filled(tmp_path):
    from altanalyze3.components.annotation.splice_survey_report import resolve_screen_symbols
    screen=[dict(gene_symbol='FOXP3',genes='')]
    with pytest.raises(ValueError,match='gene mapping'):resolve_screen_symbols(screen,None)
    p=tmp_path/'symbols.tsv';p.write_text('ENSG_FIRST\tOTHER\nLRG_62\tFOXP3\nENSG00000049768\tFOXP3\n')
    resolved=resolve_screen_symbols(screen,p)
    assert resolved[0]['genes']=='ENSG00000049768'
    assert 'absent from usage source' in resolved[0]['gene_mapping_source']


def test_shared_state_support_is_strict_and_does_not_conflate_independent_maxima():
    from altanalyze3.components.annotation.splice_survey_report import shared_read_support
    assert not shared_read_support({'A':101,'B':0},{'A':100,'B':102},100)['same_state_read_supported']
    r=shared_read_support({'A':101,'B':200},{'A':102,'B':300},100)
    assert r['n_states_both_gt_threshold']==2 and r['shared_high_read_state']=='B'
    assert r['reads_1_in_shared_state']==200
