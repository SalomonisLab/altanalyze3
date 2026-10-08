"""Analysis-record export uses native persisted evidence, never scientific data reads."""
import ast
import hashlib
import io
import json
import zipfile
from pathlib import Path
from importlib import import_module

import pytest
from fastapi.testclient import TestClient

from altanalyze3.components.cellHarmony.analysis_record import (
    QC_DEFAULTS, build_archive, build_manifest, capture_execution, record_effective, render_methods)
from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
from altanalyze3.components.cellHarmony.webapp.app import QCSettings, create_app
from altanalyze3.components.fastComm.api import FastCommParams, _resolve_default_resource_paths
from altanalyze3.components.fastComm import upstream_resources as resources
from altanalyze3.components.fastComm.verify_packaged_resources import verify_resources


def fixture_job(tmp_path, mode="supervised"):
    store = JobStore(tmp_path / "jobs")
    job = store.create_job("human", "recorded-reference", None,
                           [{"filename": "input.h5ad", "sample_name": "a"}])["job_id"]
    # Deliberately invalid H5AD: exporting methods must never read expression.
    (store.uploads_dir(job) / "input.h5ad").write_bytes(b"do not read expression")
    qc = dict(QC_DEFAULTS, analysis_mode=mode, min_genes=700, ambient_correction="yes")
    store.update_job(job, status="completed", qc=qc, analysis_mode=mode,
        model_versions={"grn": {"model_version_id": "original-bundle-hash", "runtime_versions": {"numpy": "original-version"}}},
        fastcomm_analysis={"enabled": False, "status": "failed", "message": "No built fastComm ligand-receptor bundle"})
    store.append_log(job, "[scale] skipping normalize_total and log1p")
    store.append_log(job, "...skipping QC: input is already scaled and log-transformed")
    store.append_log(job, "[params] alignment min_genes=700 ambient_correction=yes")
    return store, job


def read_archive(store, job):
    with build_archive(store.get_job(job), store.logs_dir(job).parent) as spool:
        return spool.read()


def test_report_defaults_match_interface_schema():
    assert QC_DEFAULTS == QCSettings().model_dump()


@pytest.mark.parametrize("mode", ["supervised", "unsupervised", "both"])
def test_zip_is_deterministic_complete_and_notebook_unexecuted(tmp_path, mode, monkeypatch):
    store, job = fixture_job(tmp_path, mode)
    import anndata
    monkeypatch.setattr(anndata, "read_h5ad", lambda *a, **k: pytest.fail("export read expression"))
    ref = tmp_path / "reference.tsv"
    ref.write_text("gene\tstate\nA\t1\n")
    snapshot = capture_execution(store, job, reference={
        "label": "Verified selected study", "states_tsv": str(ref),
        "study_citation": "Selected study citation", "study_url": "https://example.org/study"})
    assert snapshot["reference_files"]["states_tsv"]["sha256"] == hashlib.sha256(ref.read_bytes()).hexdigest()
    first = read_archive(store, job)
    assert first == read_archive(store, job)
    with zipfile.ZipFile(io.BytesIO(first)) as archive:
        assert set(archive.namelist()) == {"methods.md", "analysis.ipynb", "run_manifest.json",
                                         "model_provenance.json", "references.json", "logs/pipeline.log"}
        assert archive.read("logs/pipeline.log") == (store.logs_dir(job) / "pipeline.log").read_bytes()
        methods = archive.read("methods.md").decode()
        assert "| min_genes | 700 | 500 | yes |" in methods
        assert "skipping QC: input is already scaled" in methods
        assert "No built fastComm ligand-receptor bundle" in methods
        assert "original-bundle-hash" in methods and "Selected study citation" in methods
        assert "ICGS.run_icgs3" in methods if mode != "supervised" else "ICGS.run_icgs3" not in methods
        notebook = json.loads(archive.read("analysis.ipynb"))
        assert notebook["nbformat"] == 4 and notebook["nbformat_minor"] == 5
        assert notebook["metadata"]["kernelspec"]["name"] == "python3"
        assert len({cell["id"] for cell in notebook["cells"]}) == len(notebook["cells"])
        assert all(cell["cell_type"] in {"markdown", "code"} and isinstance(cell["source"], list)
                   for cell in notebook["cells"])
        for cell in notebook["cells"]:
            if cell["cell_type"] == "code":
                assert cell["execution_count"] is None and cell["outputs"] == []
                ast.parse("".join(cell["source"]))


def test_historical_report_does_not_invent_original_versions_or_defaults(tmp_path, monkeypatch):
    store, job = fixture_job(tmp_path)
    monkeypatch.setattr(import_module("altanalyze3.components.cellHarmony.analysis_record").importlib.metadata,
                        "version", lambda *a: pytest.fail("download used current versions"))
    with zipfile.ZipFile(io.BytesIO(read_archive(store, job))) as archive:
        manifest = json.loads(archive.read("run_manifest.json"))
        assert manifest["interface_defaults_at_execution"] == {}
        assert len(manifest["provenance_gaps"]) == 2
        assert "original-version" in archive.read("methods.md").decode()


def test_all_differential_runs_preserve_layer_groups_and_actual_dispatch(tmp_path):
    store, job = fixture_job(tmp_path, "both")
    capture_execution(store, job)
    runs = {"sup": {"status": "completed", "config": {"population_col": "supervised_state", "modality": "rna",
                    "comparison_type": "cells", "group1_samples": ["a"], "group2_samples": ["b"]}},
            "unsup": {"status": "completed", "config": {"population_col": "unsupervised_cluster", "modality": "lipids",
                      "comparison_type": "pseudobulk", "max_cells_per_state_sample": 500}}}
    capture_execution(store, job, key="differential:unsup", effective={"configured_method": "wilcoxon"})
    record_effective(store, job, "differential:unsup", actual_population_test="limma_like_moderated_t")
    store.update_job(job, differential_history=runs,
                     differential={"run_id": "failed", "status": "failed", "config": {"modality": "adt"}})
    manifest, _ = build_manifest(store.get_job(job), store.logs_dir(job).parent)
    assert set(manifest["differential_runs"]) == {"sup", "unsup", "failed"}
    assert manifest["analysis_records"]["main"]["differential:unsup"]["effective"] == {
        "configured_method": "wilcoxon", "actual_population_test": "limma_like_moderated_t"}
    methods = render_methods(manifest)
    assert "supervised_state" in methods and "unsupervised_cluster" in methods
    assert "limma_like_moderated_t" in methods and '"status": "failed"' in methods


def test_unstarted_differential_is_not_reported_as_executed(tmp_path):
    store, job = fixture_job(tmp_path)
    store.update_job(job, differential={"status": "not_run", "config": {}})
    manifest, _ = build_manifest(store.get_job(job), store.logs_dir(job).parent)
    assert manifest["differential_runs"] == {}
    assert "No differential comparison was recorded" in render_methods(manifest)


def test_branch_records_and_logs_retained_without_external_paths(tmp_path):
    store, job = fixture_job(tmp_path, "both")
    branch_root = store.outputs_dir(job) / "branches" / "unsupervised" / job
    (branch_root / "logs").mkdir(parents=True)
    (branch_root / "job.json").write_text(json.dumps({"analysis_records": {"icgs3": {"effective": {"umap_fit_cells": 30000}}}}))
    (branch_root / "logs" / "pipeline.log").write_text("[params] icgs3 input_normalized=True\n")
    outside = tmp_path / "private.log"
    outside.write_text("must not be included")
    (store.logs_dir(job) / "external.log").symlink_to(outside)
    with zipfile.ZipFile(io.BytesIO(read_archive(store, job))) as archive:
        assert not any("external" in name for name in archive.namelist())
        assert any("branches/unsupervised" in name for name in archive.namelist())
        manifest = json.loads(archive.read("run_manifest.json"))
        assert "input_normalized=True" in archive.read("methods.md").decode()
        assert len(manifest["analysis_records"]) == 2


def test_endpoint_returns_zip_and_is_independent_of_selected_layer(tmp_path):
    store, job = fixture_job(tmp_path, "both")
    store.update_job(job, cluster_key="unsupervised_cluster", fastcomm_analysis={"tag": "main"}, cell_state_layers={
        "default": "unsupervised_cluster", "layers": [
            {"key": "unsupervised_cluster", "fastcomm_analysis": {"tag": "main"}},
            {"key": "supervised_state", "fastcomm_analysis": {"tag": "request-specific"}}]})
    app = create_app({"JOB_STORAGE": str(store.root), "ISOLATE_JOBS": False})
    with TestClient(app) as client:
        one = client.get(f"/api/jobs/{job}/log?state_layer=supervised_state")
        two = client.get(f"/api/jobs/{job}/log?state_layer=unsupervised_cluster")
        assert one.status_code == 200 and one.content == two.content
        assert one.headers["content-type"] == "application/zip"
        assert one.headers["content-disposition"].endswith('_logs_methods.zip"')
        with zipfile.ZipFile(io.BytesIO(one.content)) as archive:
            assert json.loads(archive.read("run_manifest.json"))["fastcomm_analysis"] == {"tag": "main"}
        assert client.get("/api/jobs/does-not-exist/log").status_code == 404
        (store.logs_dir(job) / "pipeline.log").unlink()
        assert client.get(f"/api/jobs/{job}/log").status_code == 404


@pytest.mark.parametrize("species", ["human", "mouse"])
def test_complete_packaged_lr_matches_existing_native_resource_and_fallback(tmp_path, monkeypatch, species):
    verified = verify_resources()
    assert verified[species]["interactions"] == (7340 if species == "human" else 8091)
    generated = resources.bundle_paths_for_species(species)
    packaged = resources.bundle_paths_for_species(species, root=resources.PACKAGED_BUNDLE_ROOT)
    # Native local resource, when present, remains the preferred resource and is byte-identical.
    if generated["lr_table"].is_file():
        assert resources.inference_bundle_paths_for_species(species) == generated
        assert generated["lr_table"].read_bytes() == packaged["lr_table"].read_bytes()
    monkeypatch.setattr(resources, "GENERATED_BUNDLE_ROOT", tmp_path / "absent-in-clean-checkout")
    assert resources.inference_bundle_paths_for_species(species) == packaged
    lr, response, resolved = _resolve_default_resource_paths(FastCommParams(species=species, response_matrix=None))
    assert lr == packaged["lr_table"] and response is None and resolved == species
