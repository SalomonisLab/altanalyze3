from altanalyze3.components.rna2metabolite.annotations import (
    catalog, display_name, for_dataset, label_tsv, pdf_payload,
)


def aml_meta():
    return {"modalities": {"available": [{"id": "metabolite", "label": "Metabolite (AML)"}]},
            "differential": {"config": {"modality": "metabolite"}}}


def test_annotations_are_study_and_modality_scoped():
    assert not for_dataset({}, "metabolite")
    assert not for_dataset(aml_meta(), "rna")
    unrelated = {"modalities": {"available": [{"id": "metabolite", "label": "Other study"}]}}
    assert not for_dataset(unrelated)
    unrelated["modalities"]["available"][0]["annotation_source"] = "PDC000561"
    assert len(for_dataset(unrelated)) == 2115


def test_known_source_values_and_candidates_are_not_invented_names():
    annotations = for_dataset(aml_meta())
    item = annotations["Unknown 0235"]
    assert item["label"] == "Unknown 0235 · m/z 234.09706 · RP+ · RT 1.221 min"
    assert item["source_cell"] == "Table22!A464"
    assert "Sarcosine" not in annotations["Unknown0 1621"]["label"]
    assert annotations["Unknown0 1621"]["status"] == "Unidentified feature"
    assert display_name("Lactate", annotations) == "Lactate"
    assert display_name("Unknown from another study", annotations) == "Unknown from another study"


def test_presentation_copy_preserves_lookup_ids_and_numerical_results():
    original = {"gene": "Unknown 0235", "rows": [{"gene": "Unknown 0235", "values": [None, 0, 1.25]}]}
    labeled = pdf_payload(original, aml_meta())
    assert "m/z" in labeled["gene"]
    assert labeled["rows"][0]["values"] == original["rows"][0]["values"]
    assert original["gene"] == original["rows"][0]["gene"] == "Unknown 0235"
    assert pdf_payload(original, aml_meta(), "rna") == original


def test_marker_export_labels_keep_cluster_prefix_and_exact_measurements():
    text = "\tcell1\tcell2\nHSC-1:Unknown 0235\t-1.234\t0\nHSC-1:ATP\t2\t3\n"
    labeled = label_tsv(text, catalog())
    assert "HSC-1:Unknown 0235 · m/z 234.09706" in labeled
    assert [line.partition("\t")[2] for line in text.splitlines()] == [line.partition("\t")[2] for line in labeled.splitlines()]


def test_literal_metabolite_names_remain_selectable_in_gene_set_plots():
    from altanalyze3.components.cellHarmony.webapp.app import _split_expression_features
    import numpy as np
    cache = {"var_names": np.array(["Unknown 0235", "Unknown 0309", "1,3-Diphenylurea", "TP53", "MPO"])}
    assert _split_expression_features("Unknown 0235", cache) == ["Unknown 0235"]
    assert _split_expression_features("1,3-Diphenylurea\nUnknown 0235", cache) == ["1,3-Diphenylurea", "Unknown 0235"]
    assert _split_expression_features("Unknown 0235, Unknown 0309", cache) == ["Unknown 0235", "Unknown 0309"]
    assert _split_expression_features("TP53 MPO", cache) == ["TP53", "MPO"]
