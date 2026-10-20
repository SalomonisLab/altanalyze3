from altanalyze3.components.iso2function.protein_isoform_audit import project_feature, align_mapping
from altanalyze3.components.iso2function.protein_isoform_audit import validate_native_inputs
from altanalyze3.components.iso2function.protein_isoform_audit import matched_activation
import copy
import pytest


def native_fixture():
    rows = [dict(gene_symbol='G', isoform_id='I', protein_sequence='MA',
                 cds_sequence='ATGGCT', cdna_sequence='GGATGGCTTAA')]
    decision = dict(pending_question=None, actual_user_answer='User supplied native database',
                    required_roster=[dict(gene_symbol='G', isoform_id='I')])
    return rows, {('G', 'I')}, decision


def test_native_gate_accepts_exact_roster_and_translation():
    validate_native_inputs(*native_fixture())


def test_blank_clone_ids_cannot_transfer_activation():
    rows = [dict(gene_symbol='G', isoform_id='I', tfiso_exact_clones='')]
    measurements = [dict(clone_id='', gene_symbol='OTHER', M1H_mean='10')]
    assert matched_activation(rows, measurements) == []


def test_exact_activation_match_preserves_source_gene():
    rows = [dict(gene_symbol='G', isoform_id='I', tfiso_exact_clones='G-1')]
    measurements = [dict(clone_id='G-1', gene_symbol='G', M1H_mean='2')]
    result = matched_activation(rows, measurements)
    assert len(result) == 1 and result[0]['isoform_id'] == 'I'
    measurements[0]['gene_symbol'] = 'OTHER'
    with pytest.raises(ValueError, match='gene identity'):
        matched_activation(rows, measurements)


def test_native_gate_rejects_missing_and_duplicate_identities():
    rows, rendered, decision = native_fixture()
    with pytest.raises(ValueError, match='roster'):
        validate_native_inputs([], rendered, decision)
    with pytest.raises(ValueError, match='Duplicate'):
        validate_native_inputs(rows + rows, rendered, decision)


def test_native_gate_rejects_identity_swap_even_when_count_unchanged():
    rows, rendered, decision = native_fixture()
    rows[0]['isoform_id'] = 'OTHER'
    with pytest.raises(ValueError, match='Expected roster'):
        validate_native_inputs(rows, {('G', 'OTHER')}, decision)


def test_native_gate_rejects_pending_decision_and_sequence_substitution():
    rows, rendered, decision = native_fixture()
    pending = copy.deepcopy(decision); pending['pending_question'] = 'Which source?'
    with pytest.raises(ValueError, match='Unresolved source'):
        validate_native_inputs(rows, rendered, pending)
    rows[0]['protein_sequence'] = 'MM'
    with pytest.raises(ValueError, match='translation differs'):
        validate_native_inputs(rows, rendered, decision)


def test_intact_feature_keeps_one_based_inclusive_coordinates():
    reference = "MALWMRLLPLLALLALWGPDPAAA"
    query = reference[4:]
    f = project_feature(reference, query, 6, 16, align_mapping(reference, query))
    assert f["status"] == "intact_sequence"
    assert (f["query_start"], f["query_end"]) == (2, 12)
    assert f["exact_fraction"] == 1


def test_partial_feature_is_never_intact():
    reference = "MALWMRLLPLLALLALWGPDPAAA"
    query = reference[:10]
    f = project_feature(reference, query, 6, 16, align_mapping(reference, query))
    assert f["status"] == "altered_or_partial"
    assert f["exact_fraction"] < 1


def test_inserted_residues_do_not_count_as_intact_domain():
    reference = "MALWMRLLPLLALLALWGPDPAAA"
    query = reference[:10] + "EEEE" + reference[10:]
    f = project_feature(reference, query, 6, 16, align_mapping(reference, query))
    assert f["status"] == "altered_or_partial"


def test_repeated_feature_cannot_be_placed_unambiguously():
    f = project_feature("MABCDE", "MABCDEABCDE", 2, 6, {i:i for i in range(6)})
    assert f["status"] == "ambiguous_repeated_sequence"
    assert f["query_start"] == ""
