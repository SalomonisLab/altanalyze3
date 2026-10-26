"""Prevent a missing transcript from borrowing another RNA's sequence."""

import pytest

from altanalyze3.components.annotation.isoform_predict import (
    _lookup_sequence,
    load_sequence_fasta,
)


def test_missing_pb_does_not_use_generic_pb_sequence():
    sequences = {"PB": "FOREIGN", "PB.1.1": "LITERAL"}
    assert _lookup_sequence(sequences, "PB.22997.203") == ""
    assert _lookup_sequence(sequences, "PB.1.1") == "LITERAL"


def test_physical_library_ids_remain_distinct(tmp_path):
    path = tmp_path / "proteins.fa"
    path.write_text(
        ">101.first__source010\nMAAA\n"
        ">101.second__source010\nMCCC\n"
    )
    sequences = load_sequence_fasta(str(path))
    assert _lookup_sequence(sequences, "101.first__source010") == "MAAA"
    assert _lookup_sequence(sequences, "101.second__source010") == "MCCC"
    assert _lookup_sequence(sequences, "101.missing__source010") == ""


def test_gene_qualified_pb_retains_complete_pb_identity(tmp_path):
    path = tmp_path / "proteins.fa"
    path.write_text(
        ">ENSG00000115935_PB.24879.86\nMAAA\n"
        ">ENSG00000140983_PB.146926.188\nMCCC\n"
    )
    sequences = load_sequence_fasta(str(path))
    assert _lookup_sequence(sequences, "PB.24879.86") == "MAAA"
    assert _lookup_sequence(sequences, "PB.146926.188") == "MCCC"
    assert _lookup_sequence(sequences, "PB.22997.203") == ""


def test_documented_ensembl_version_alias_remains_supported(tmp_path):
    path = tmp_path / "proteins.fa"
    path.write_text(">ENSG00000135363_ENST00000257818.3\nMAAA\n")
    sequences = load_sequence_fasta(str(path))
    assert _lookup_sequence(sequences, "ENST00000257818.3") == "MAAA"
    assert _lookup_sequence(sequences, "ENST00000257818") == "MAAA"


def test_conflicting_source_identity_stops_before_sequence_substitution(tmp_path):
    path = tmp_path / "proteins.fa"
    path.write_text(
        ">ENSG00000115935_PB.1.1\nMAAA\n"
        ">ENSG00000140983_PB.1.1\nMCCC\n"
    )
    with pytest.raises(ValueError, match="Ambiguous sequence identifier 'PB.1.1'"):
        load_sequence_fasta(str(path))
