"""`_absolutise_asset_paths` must resolve paths written relative to the PROJECT.

stage_discover writes `scalable_viewer/fastComm_v7/...` and a relative `bundle_dir`.
Tried under the release folder alone those became scalable_viewer/scalable_viewer/...,
so fastComm answered "fastComm scores are unavailable" and marker networks came back
empty while every file existed (COPD, 2026-09-28).
"""
from pathlib import Path

from altanalyze3.components.visualization.scalable_viewer.scalable_app import _absolutise_asset_paths


def layout(tmp_path):
    project = tmp_path / "COPD-atlas"
    release = project / "scalable_viewer"
    (release / "assets_integrated").mkdir(parents=True)
    (release / "bundles_integrated" / "COPD-metacells").mkdir(parents=True)
    (release / "fastComm_v7").mkdir()
    (release / "fastComm_v7" / "scores.tsv").write_text("x")
    (release / "in_release.tsv").write_text("x")
    return project, release


def test_project_relative_paths_resolve_one_level_above_the_release(tmp_path):
    project, release = layout(tmp_path)
    data = {"bundle_dir": "bundles_integrated/COPD-metacells",
            "fastcomm_analysis": {"scores_tsv": "scalable_viewer/fastComm_v7/scores.tsv"},
            "networks": [{"id": "AF", "tsv": "scalable_viewer/fastComm_v7/scores.tsv"}]}
    out = _absolutise_asset_paths(data, release_root=str(release))
    assert out["fastcomm_analysis"]["scores_tsv"] == str(release / "fastComm_v7" / "scores.tsv")
    assert Path(out["networks"][0]["tsv"]).is_file()


def test_release_relative_and_absolute_paths_are_untouched(tmp_path):
    project, release = layout(tmp_path)
    absolute = str(release / "in_release.tsv")
    data = {"bundle_dir": "bundles_integrated/COPD-metacells",
            "a": "in_release.tsv", "b": absolute}
    out = _absolutise_asset_paths(data, release_root=str(release))
    assert out["a"] == str(release / "in_release.tsv")
    assert out["b"] == absolute


def test_unresolvable_value_is_returned_unchanged_for_the_caller_to_report(tmp_path):
    project, release = layout(tmp_path)
    data = {"bundle_dir": "bundles_integrated/COPD-metacells", "x": "nowhere/at/all.tsv"}
    assert _absolutise_asset_paths(data, release_root=str(release))["x"] == "nowhere/at/all.tsv"


def test_relative_bundle_dir_still_yields_the_project_root(tmp_path):
    """The bundle's third parent is the project; a relative bundle_dir must not drop it."""
    project, release = layout(tmp_path)
    (project / "analysis").mkdir()
    (project / "analysis" / "t.tsv").write_text("x")
    data = {"bundle_dir": "bundles_integrated/COPD-metacells", "x": "analysis/t.tsv"}
    assert _absolutise_asset_paths(data, release_root=str(release))["x"] == str(project / "analysis" / "t.tsv")
