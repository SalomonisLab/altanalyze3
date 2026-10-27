"""Integrated RNA/protein views must not mislabel an audited literal protein."""
from copy import deepcopy

import pytest

from altanalyze3.components.visualization.isv_web.isv_web_remote import _architecture_structure_length


def test_reported_length_discrepancy_is_preserved_and_literal_sequence_is_displayed():
    model = {'protein_length': 23, 'count': 158}
    track = {'isoform_id': 'ENST00000412051', 'reported_aa': 23, 'protein_aa': 115}
    before = deepcopy(model)
    audit = _architecture_structure_length(model, track)
    assert audit['displayed_literal_sequence_aa'] == 115
    assert audit['source_RNA_protein_length'] == audit['source_architecture_reported_aa'] == 23
    assert model == before


def test_unreconciled_source_length_still_fails():
    with pytest.raises(ValueError, match='differs from audited architecture'):
        _architecture_structure_length({'protein_length': 99},
                                       {'isoform_id': 'X', 'reported_aa': 23, 'protein_aa': 115})


def test_matching_literal_length_and_missing_source_length_are_supported():
    track = {'isoform_id': 'X', 'reported_aa': None, 'protein_aa': 198}
    for source in (198, None):
        audit = _architecture_structure_length({'protein_length': source}, track)
        assert audit['displayed_literal_sequence_aa'] == 198
        assert audit['source_RNA_protein_length'] == source
