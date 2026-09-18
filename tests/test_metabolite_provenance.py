"""Provenance must not turn ambiguous MS annotations into compound identities."""
import csv
import json
import zipfile
from xml.sax.saxutils import escape

import pytest

from altanalyze3.components.rna2metabolite.provenance import audit, workbook_rows


def workbook(path):
    ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    sheets = {
        "Table21": [
            ["Metabolite", "mz", "Rt.min", "Tags", "Synonyms.or.Isomers", "Formula"],
            ["Known", "101", "2", "MS2", "", "CH4"],
            ["Unknown 1", "202", "0", "MS2", "Candidate only", "C2H4"],
        ],
        "Table22": [
            ["Metabolite", "mz", "Rt.min", "Tags"],
            ["Known", "303", "4", "MS2"],
        ],
        "Table24": [["feature", "score"], ["HILIC:Unknown 1", "1"]],
    }
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("xl/workbook.xml", f'<workbook xmlns="{ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>' + ''.join(
            f'<sheet name="{name}" r:id="r{i}"/>' for i, name in enumerate(sheets)) + '</sheets></workbook>')
        z.writestr("xl/_rels/workbook.xml.rels", '<Relationships>' + ''.join(
            f'<Relationship Id="r{i}" Target="worksheets/s{i}.xml"/>' for i in range(3)) + '</Relationships>')
        for i, rows in enumerate(sheets.values()):
            xml = f'<worksheet xmlns="{ns}"><sheetData>'
            for n, values in enumerate(rows, 1):
                xml += f'<row r="{n}">' + ''.join(
                    f'<c r="{chr(65+c)}{n}" t="inlineStr"><is><t>{escape(value)}</t></is></c>'
                    for c, value in enumerate(values) if value != '') + '</row>'
            z.writestr(f"xl/worksheets/s{i}.xml", xml + '</sheetData></worksheet>')


def test_assay_specific_join_preserves_unknowns_and_source_coordinates(tmp_path):
    source = tmp_path / "source.xlsx"
    workbook(source)
    selection = tmp_path / "selection.tsv"
    selection.write_text("\tprotocol\tdetection\tmedian_signal\nKnown\tRP\t1\t2\nUnknown 1\tHILIC\t0.9\t1\n")
    pdc = tmp_path / "pdc.csv"
    pdc.write_text(',sample\nRP:Known,1\nHILIC:Known,2\nHILIC:Unknown 1,3\n')
    output = tmp_path / "audit"
    report = audit(source, selection, pdc, output)
    with (output / "aml_metabolite_provenance.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert rows[0]["mz"] == "303"  # the selected RP row, not the first named row
    assert rows[1]["model_feature_id"] == "Unknown 1"
    assert rows[1]["identification_status"] == "unidentified_in_source"
    assert rows[1]["source_synonyms_or_isomers"] == "Candidate only"
    assert rows[1]["source_cell"] == "Table21!A3"
    assert rows[1]["retention_time_min"] == "0"
    assert json.loads(rows[1]["raw_source_annotations_json"])["Formula"] == "C2H4"
    assert report["confirmed_unknown_renamings"] == 0
    assert report["source_pdc_feature_ids_identical"]
    assert report["selected_measurements_in_paper_ica"] == 1
    # A changed deposit must fail rather than silently assigning a different ID.
    pdc.write_text(',sample\nRP:Known,1\n')
    with pytest.raises(ValueError, match="missing from deposited"):
        audit(source, selection, pdc, output)


def test_sparse_inline_cells_keep_original_columns(tmp_path):
    path = tmp_path / "source.xlsx"
    workbook(path)
    number, values = list(workbook_rows(path, "Table21"))[1]
    assert number == 2
    assert "E" not in values
    assert values["F"] == "CH4"
