from altanalyze3.components.annotation.splice_survey_correspondence import compare_peptides,ranges


def test_internal_splice_gap_is_not_terminal_truncation():
    full='MABCDEFGHIKLMNPQRSTVWY';short=full[:8]+full[11:]
    r=compare_peptides(full,short)
    assert r['unmapped_aa_1']==[[9,11]] and r['unmapped_aa_2']==[] and r['substitutions']==0
    r=compare_peptides(full,full[5:])
    assert r['unmapped_aa_1']==[[1,5]] and r['identical_aligned_residues']==len(full)-5


def test_equal_protein_lengths_can_contain_mutually_exclusive_sequence():
    r=compare_peptides('MAAAKLMNPQRSTVWY','MAAAKLMNACRSTVWY')
    assert r['aa_1']==r['aa_2'] and r['substitutions']==2
    assert ranges([0,1,4])==[[1,2],[5,5]]
