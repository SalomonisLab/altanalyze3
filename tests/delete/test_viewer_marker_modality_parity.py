"""scALABLE-viewer parity: per-modality MarkerFinder registration.

flask/pipeline.py gives an uploaded job one `marker_analysis_by_modality` entry per
imputed modality. A bundle hard-coded `{"rna": ...}`, so `cross_pathways.marker_tables`
scored one modality of five on the COPD viewer. `_marker_analysis_by_modality` now
registers what the bundle ships and stays RNA-only when it ships nothing else.
"""
import pandas as pd
import pytest

from altanalyze3.components.visualization.scalable_viewer import bundle_meta
from altanalyze3.components.cellHarmony.webapp.modality_markers import marker_tables


RNA = {"enabled": True, "status": "completed", "markers_tsv": "/rna/markers.tsv",
       "populations": ["AT1", "AT2"]}


def write_markers(path, cluster, genes):
    frame = pd.DataFrame({"Gene": genes, "cluster": [cluster] * len(genes),
                          "rho": [0.9] * len(genes)})
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


def test_no_modality_assets_keeps_rna_only():
    for assets in (None, {}, {"adt": {}}):
        assert bundle_meta._marker_analysis_by_modality(RNA, assets) == {"rna": RNA}


def test_declared_but_absent_file_is_not_registered(tmp_path):
    assets = {"adt": {"markers_tsv": str(tmp_path / "missing.tsv")}}
    assert bundle_meta._marker_analysis_by_modality(RNA, assets) == {"rna": RNA}


def test_present_modality_tables_are_registered(tmp_path):
    adt = write_markers(tmp_path / "adt_markers.tsv", "AT1", ["CD31", "CD45"])
    lipid = write_markers(tmp_path / "lipid_markers.tsv", "AT1", ["PC(34:1)"])
    absent = str(tmp_path / "gone.tsv")
    result = bundle_meta._marker_analysis_by_modality(RNA, {
        "ADT": {"markers_tsv": adt, "redundant_markers_tsv": absent},
        "metabolites": {"markers_tsv": lipid, "populations": ["AT1"]},
        "grn_tf": {"markers_tsv": absent},
    })
    # "ADT" and "metabolites" reach their canonical ids through _normalize_modality_id.
    assert sorted(result) == ["adt", "metabolite", "rna"]
    assert result["rna"] == RNA
    assert result["adt"]["markers_tsv"] == adt
    assert "redundant_markers_tsv" not in result["adt"]        # absent file is dropped
    assert result["adt"]["populations"] == RNA["populations"]  # inherited when unstated
    assert result["metabolite"]["populations"] == ["AT1"]


def test_cross_modality_marker_tables_see_every_registered_modality(tmp_path):
    """The measured defect: one modality of five reaching cross_pathways."""
    rna_tsv = write_markers(tmp_path / "rna_markers.tsv", "AT1", ["SFTPC", "AGER"])
    adt_tsv = write_markers(tmp_path / "adt_markers.tsv", "AT1", ["CD31"])
    rna = dict(RNA, markers_tsv=rna_tsv)
    available = [{"id": m, "label": m.upper()} for m in ("rna", "adt", "lipid")]

    before = {"modalities": {"available": available},
              "marker_analysis": rna, "marker_analysis_by_modality": {"rna": rna}}
    after = {"modalities": {"available": available}, "marker_analysis": rna,
             "marker_analysis_by_modality": bundle_meta._marker_analysis_by_modality(
                 rna, {"adt": {"markers_tsv": adt_tsv}})}

    assert sorted(marker_tables(before)[0]) == ["rna"]
    assert sorted(marker_tables(after)[0]) == ["adt", "rna"]
    # A modality the bundle does not ship is still reported, never silently dropped.
    coverage = {c["modality"]: c["status"] for c in marker_tables(after)[1]}
    assert coverage["lipid"] == "no retained MarkerFinder correlation scores"


def test_build_meta_accepts_the_new_argument():
    import inspect
    params = inspect.signature(bundle_meta.build_meta).parameters
    assert params["marker_assets_by_modality"].default is None


# --- producer side: prepare_assets discovery -------------------------------------

from altanalyze3.components.visualization.scalable_viewer.prepare_assets import (  # noqa: E402
    find_modality_marker_tables)


def build_run_outputs(tmp_path, modalities):
    """Mimic flask/pipeline.py `_emit_modality_marker_heatmap` output layout."""
    outputs = tmp_path / "outputs"
    rna_dir = outputs / "marker_heatmap"
    rna_dir.mkdir(parents=True)
    rna = write_markers(rna_dir / "cell_state_marker_heatmap_markers.tsv", "AT1", ["SFTPC"])
    for modality in modalities:
        d = outputs / f"marker_heatmap_{modality}"
        d.mkdir()
        write_markers(d / f"cell_state_{modality}_marker_heatmap_markers.tsv", "AT1", ["F1"])
        write_markers(d / f"cell_state_{modality}_marker_heatmap_redundant_markers.tsv",
                      "AT1", ["F1", "F2"])
    return rna


def test_discovery_finds_every_imputed_modality(tmp_path):
    rna = build_run_outputs(tmp_path, ["adt", "metabolite", "lipid", "grn_tf"])
    found = find_modality_marker_tables(rna)
    assert sorted(found) == ["adt", "grn_tf", "lipid", "metabolite"]
    for entry in found.values():
        assert entry["markers_tsv"].endswith("_marker_heatmap_markers.tsv")
        assert entry["redundant_markers_tsv"].endswith("_redundant_markers.tsv")
        # The unique table must never be mistaken for the redundant one.
        assert entry["markers_tsv"] != entry["redundant_markers_tsv"]


def test_discovery_is_empty_without_imputed_modalities(tmp_path):
    assert find_modality_marker_tables(build_run_outputs(tmp_path, [])) == {}


def test_discovery_skips_a_directory_with_no_unique_table(tmp_path):
    rna = build_run_outputs(tmp_path, ["adt"])
    empty = tmp_path / "outputs" / "marker_heatmap_lipid"
    empty.mkdir()
    (empty / "cell_state_lipid_marker_heatmap_redundant_markers.tsv").write_text("x")
    assert sorted(find_modality_marker_tables(rna)) == ["adt"]


def test_producer_and_consumer_agree(tmp_path):
    """prepare_assets discovery feeds bundle_meta registration end to end."""
    rna = build_run_outputs(tmp_path, ["adt", "lipid"])
    registered = bundle_meta._marker_analysis_by_modality(
        dict(RNA, markers_tsv=str(rna)), find_modality_marker_tables(rna))
    assert sorted(registered) == ["adt", "lipid", "rna"]
