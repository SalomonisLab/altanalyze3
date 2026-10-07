import copy
import json
import pickle
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.model_registry.registry import describe_model, write_run_provenance
from altanalyze3.components.model_registry.proposals import validate_proposal


def test_identity_tracks_bytes_and_code_but_not_location(tmp_path):
    bundle = tmp_path / "model.pkl"
    bundle.write_bytes(b"aaaa")
    code = tmp_path / "api.py"
    code.write_text("version one")
    first = describe_model("test", {"bundle": bundle}, code_paths=[code])
    renamed = tmp_path / "renamed.pkl"
    renamed.write_bytes(bundle.read_bytes())
    assert describe_model("test", {"bundle": renamed}, code_paths=[code])["model_version_id"] == first["model_version_id"]
    bundle.write_bytes(b"bbbb")  # Same length, same path; no stat-based hash cache.
    changed = describe_model("test", {"bundle": bundle}, code_paths=[code])
    assert changed["model_version_id"] != first["model_version_id"]
    code.write_text("version two")
    code_changed = describe_model("test", {"bundle": bundle}, code_paths=[code])
    assert changed["model_version_id"] == code_changed["model_version_id"]
    assert changed["analysis_model_version_id"] != code_changed["analysis_model_version_id"]
    assert changed["inference_version_id"] != code_changed["inference_version_id"]


def test_missing_artifacts_cannot_receive_a_version(tmp_path):
    with pytest.raises(FileNotFoundError):
        describe_model("test", {"bundle": tmp_path / "missing"}, code_paths=[])


def test_catalog_lookup_is_by_hash_not_name(tmp_path):
    artifact = tmp_path / "model"
    artifact.write_bytes(b"model")
    identity = describe_model("test", {"bundle": artifact}, code_paths=[])
    catalog = tmp_path / "catalog.json"
    catalog.write_text(json.dumps({"models": [{"model_version_id": identity["model_version_id"],
                                              "status": "observed-default", "name": "known"}]}))
    assert describe_model("test", {"bundle": artifact}, code_paths=[], catalog_path=catalog)["registry_status"] == "observed-default"
    artifact.write_bytes(b"replacement")
    assert describe_model("test", {"bundle": artifact}, code_paths=[], catalog_path=catalog)["registry_status"] == "unregistered"


class SyntheticRegressor:
    def predict(self, matrix, **kwargs):
        return np.asarray(matrix)[:, :1] * 2


@pytest.mark.parametrize("component", ["rna2lipid", "rna2adt", "rna2grn", "rna2metabolite", "rna2lipid_aml"])
def test_loader_and_prediction_keep_provenance_and_values(tmp_path, component):
    from importlib import import_module
    module = component.replace("rna2lipid_aml", "rna2lipid.aml")
    api = import_module("altanalyze3.components." + module + ".api")
    bundle_path = tmp_path / "synthetic.pkl"
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler(with_mean=False, with_std=False).fit(pd.DataFrame(np.ones((2, 2)), columns=["GENE1", "GENE2"]))
    record = {"model": SyntheticRegressor(), "scaler_x": scaler, "scaler_y": None,
              "X_columns": ["GENE1", "GENE2"], "Y_columns": ["target"],
              "mu": np.zeros(2), "sd": np.ones(2), "sel_idx": [[0]],
              "coef": [[2.0]], "intercept": [0.0]}
    bundle_path.write_bytes(pickle.dumps(record))
    bundle = api.load_bundle(bundle_path)
    original_version = bundle.model_provenance["model_version_id"]
    frame = pd.DataFrame([[1., 2.], [3., 4.]], columns=record["X_columns"], index=["s1", "s2"])
    if component in ("rna2metabolite", "rna2lipid_aml", "rna2grn"):
        result = bundle.predict_from_dataframe(frame, normalized=True)
    else:
        result = bundle.predict_from_dataframe(frame)
    np.testing.assert_allclose(result.predictions.to_numpy().ravel(), [2., 6.])
    assert result.summary["model_version_id"] == original_version
    assert result.summary["component"] == component
    # A loaded object retains the identity of its loaded model, even if the path changes later.
    bundle_path.write_bytes(b"replacement")
    assert bundle.model_info()["model_version_id"] == original_version
    matrix = ad.AnnData(X=result.predictions.to_numpy(), uns={"prediction_summary": result.summary})
    output = tmp_path / "result.h5ad"
    matrix.write_h5ad(output)
    assert ad.read_h5ad(output).uns["prediction_summary"]["model_version_id"] == original_version


def test_run_manifest_serializes_h5ad_values(tmp_path):
    path = write_run_provenance(tmp_path, {"lipids": {"matched_genes": np.int64(3)}}, application="scALABLE")
    assert json.loads(path.read_text())["models"]["lipids"]["matched_genes"] == 3


def example():
    return json.loads((Path(__file__).resolve().parents[3] / "model_repository/proposal.example.json").read_text())


def test_proposal_cannot_self_promote_or_omit_rosters():
    record = example()
    validate_proposal(record)
    for field, value in [("status", "default"), ("approval", {"approved": True}),
                         ("input_roster", {}), ("model_version_id", "fake")]:
        altered = copy.deepcopy(record)
        altered[field] = value
        with pytest.raises(ValueError):
            validate_proposal(altered)


def test_fastcomm_writes_resource_versions(tmp_path):
    from altanalyze3.components.fastComm.api import FastCommParams, run_fastcomm
    lr = tmp_path / "lr.tsv"
    lr.write_text("ligand\treceptor\nLIG\tREC\n")
    expression = tmp_path / "expression.tsv"
    pd.DataFrame([[2, 1], [1, 2]], columns=["LIG", "REC"], index=["c1", "c2"]).to_csv(expression, sep="\t")
    metadata = tmp_path / "metadata.tsv"
    pd.DataFrame({"cell_state": ["A", "B"]}, index=["c1", "c2"]).to_csv(metadata, sep="\t")
    output = tmp_path / "scores.tsv"
    result = run_fastcomm(FastCommParams(expression=expression, metadata=metadata, lr_table=lr,
                         response_matrix=None, output=output, state_pair_output=tmp_path / "pairs.tsv"))
    companion = json.loads(Path(str(output) + ".model_provenance.json").read_text())
    assert companion["models"]["model_version_id"] == result.summary["model_version_id"]
    assert set(result.summary["artifacts"]) == {"ligand_receptor"}
    exported = pd.read_csv(output, sep="\t")
    assert list(exported.columns) == list(result.scores.columns)
    np.testing.assert_allclose(exported["fastcomm_score"], result.scores["fastcomm_score"])
    assert Path(str(tmp_path / "pairs.tsv") + ".model_provenance.json").exists()


def test_static_defaults_never_execute_code():
    from altanalyze3.components.model_registry.static_defaults import api_selection
    source = 'BUNDLE_DIR = Path(__file__).parent\nREFERENCE_BUNDLES = {"lung": BUNDLE_DIR / "lung.pkl"}\nDEFAULT_REFERENCE = "lung"\nDEFAULT_BUNDLE_PATH = REFERENCE_BUNDLES[DEFAULT_REFERENCE]\nraise RuntimeError("must not execute")'
    assert api_selection(source, "components/rna2grn/api.py") == "components/rna2grn/lung.pkl"


def test_result_provenance_survives_exports_without_inventing_old_ids(tmp_path):
    from altanalyze3.components.model_registry.registry import read_result_provenance, write_provenance
    source = tmp_path / "predictions.csv"
    source.write_text("gene,value\nG1,1\n")
    assert read_result_provenance(source) == {}
    record = {"model_version_id": "test:sha256:" + "a" * 64}
    write_provenance(source, record)
    assert read_result_provenance(source) == record
    matrix = ad.AnnData(X=np.ones((2, 2)), uns={"prediction_summary": record})
    h5ad = tmp_path / "predictions.h5ad"
    matrix.write_h5ad(h5ad)
    assert read_result_provenance(h5ad) == record


def test_viewer_modality_retains_recorded_identity(tmp_path):
    from altanalyze3.components.visualization.scalable_viewer import bundle, precompute
    record = {"model_version_id": "test:sha256:" + "a" * 64}
    matrix = ad.AnnData(X=np.array([[1., 2.], [2., 1.]], dtype=np.float32),
                       obs=pd.DataFrame(index=["c1", "c2"]),
                       var=pd.DataFrame(index=["F1", "F2"]),
                       uns={"prediction_summary": record})
    source = tmp_path / "result.h5ad"
    matrix.write_h5ad(source)
    output = tmp_path / "bundle"
    output.mkdir()
    info = precompute.ingest_modality("test", str(source), paths=bundle.BundlePaths(str(output), "test"),
            barcodes=["c1", "c2"], states=["A", "B"], state_code=np.array([0, 1]), state_n=np.array([1, 1]))
    assert info["model_provenance"] == record


def test_static_default_can_follow_a_release_module():
    from altanalyze3.components.model_registry.static_defaults import api_selection
    source = 'from .release import RELEASE_BUNDLE_PATH\nDEFAULT_BUNDLE_PATH = RELEASE_BUNDLE_PATH'
    release = 'HERE = Path(__file__).resolve().parent\nRELEASE_BUNDLE_PATH = HERE / "release.pkl"'
    assert api_selection(source, "components/rna2lipid/api.py", lambda path: release) == "components/rna2lipid/release.pkl"


def test_catalog_refresh_retains_prior_artifacts_and_code():
    from altanalyze3.components.model_registry.snapshot import merge_catalog
    def record(mid, iid):
        return {"model_version_id": mid, "inference_version_id": iid,
                "analysis_model_version_id": mid + iid, "inference_code_sha256": {"api": iid},
                "status": "observed-default"}
    original = record("original", "v1")
    unchanged_artifact = record("retained", "v1")
    current = {"models": [record("new", "v1"), record("retained", "v2")], "defaults": []}
    result = merge_catalog({"models": [original, unchanged_artifact]}, current)
    by_id = {r["model_version_id"]: r for r in result["models"]}
    assert set(by_id) == {"original", "retained", "new"}
    assert by_id["original"]["status"] == "historical-default"
    assert [v["inference_version_id"] for v in by_id["retained"]["inference_versions"]] == ["v1", "v2"]
    assert original["status"] == "observed-default"  # Do not mutate the archived record.
    assert merge_catalog(result, current) == result


def test_archived_catalog_is_stable(tmp_path):
    from altanalyze3.components.model_registry.snapshot import archive_catalog
    record = {"models": [], "defaults": []}
    name = archive_catalog(record, tmp_path)
    assert archive_catalog(record, tmp_path) == name
    assert json.loads((tmp_path / "snapshots" / name).read_text()) == record


def test_downstream_method_has_its_own_version(tmp_path):
    from altanalyze3.components.model_registry.registry import describe_application
    code = tmp_path / "scale.py"
    code.write_text("version1")
    before = describe_application("scALABLE", code_paths=[code])
    code.write_text("version2")
    after = describe_application("scALABLE", code_paths=[code])
    assert before["application_method_version_id"] != after["application_method_version_id"]


def test_refresh_preserves_download_and_release_evidence():
    from altanalyze3.components.model_registry.snapshot import merge_catalog
    model = {"model_version_id": "test", "inference_version_id": "v1",
             "analysis_model_version_id": "av1", "inference_code_sha256": {},
             "artifacts": {"bundle": {"sha256": "hash", "url": "https://example.org/model"}},
             "release_evidence": {"manifest": "release.json"}}
    current = copy.deepcopy(model)
    current.pop("release_evidence")
    current["artifacts"]["bundle"].pop("url")
    result = merge_catalog({"registry_release": "release1", "models": [model]},
                           {"models": [current], "defaults": []})
    assert result["registry_release"] == "release1"
    assert result["models"][0]["artifacts"]["bundle"]["url"] == "https://example.org/model"
    assert result["models"][0]["release_evidence"] == model["release_evidence"]


def test_readable_versions_keep_artifact_ids_and_default_links(tmp_path):
    from altanalyze3.components.model_registry.layout import apply_labels, write_layout
    mid = "rna2lipid:sha256:original"
    catalog = {"models": [{"model_version_id": mid, "artifacts": {}, "status": "observed-default"}],
               "defaults": [{"model_version_id": mid}]}
    apply_labels(catalog, {mid: {"modality": "rna2lipid", "version": "v1.1", "variant": "lung"}})
    assert catalog["models"][0]["model_version_id"] == mid
    assert catalog["defaults"][0]["model_display_id"] == "rna2lipid/lung/v1.1"
    write_layout(catalog, tmp_path)
    assert (tmp_path / "models/rna2lipid/v1.1/lung/model.json").exists()
    assert (tmp_path / "models/rna2lipid/v1.1/README.md").exists()


def test_two_artifacts_cannot_reuse_a_readable_version():
    from altanalyze3.components.model_registry.layout import apply_labels
    catalog = {"models": [{"model_version_id": "a"}, {"model_version_id": "b"}], "defaults": []}
    label = {"modality": "rna2lipid", "version": "v1.0", "variant": "lung"}
    with pytest.raises(ValueError, match="multiple artifacts"):
        apply_labels(catalog, {"a": label, "b": label})


def test_result_summary_exposes_human_version(tmp_path):
    artifact = tmp_path / "model"
    artifact.write_bytes(b"model")
    identity = describe_model("test", {"bundle": artifact}, code_paths=[])
    catalog = tmp_path / "catalog.json"
    catalog.write_text(json.dumps({"models": [{"model_version_id": identity["model_version_id"],
        "status": "observed-default", "name": "test", "model_modality": "rna2lipid",
        "model_version": "v1.1", "model_variant": "lung", "model_display_id": "rna2lipid/lung/v1.1"}]}))
    result = describe_model("test", {"bundle": artifact}, code_paths=[], catalog_path=catalog)
    assert result["model_version"] == "v1.1"
    assert result["model_display_id"] == "rna2lipid/lung/v1.1"


def test_assigned_model_version_cannot_be_renumbered():
    from altanalyze3.components.model_registry.layout import apply_labels
    catalog = {"models": [{"model_version_id": "known", "model_version": "v1.0"}], "defaults": []}
    with pytest.raises(ValueError, match="Cannot relabel"):
        apply_labels(catalog, {"known": {"modality": "rna2lipid", "version": "v1.1", "variant": "lung"}})
