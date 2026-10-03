"""Wrong 10x artifacts receive actionable errors before allocating expression."""
import io

import h5py
import numpy as np
import pytest
from fastapi.testclient import TestClient

from altanalyze3.components.cellHarmony.input_validation import validate_10x_h5
from altanalyze3.components.cellHarmony.merge_inputs import concat_matching_10x
from altanalyze3.components.cellHarmony.webapp.app import create_app


def _write_h5(path, kind):
    with h5py.File(path, "w") as handle:
        if kind == "molecule":
            handle.create_dataset("count", data=np.array([1, 2], dtype=np.uint32))
            handle.create_dataset("barcode_idx", data=np.array([0, 1], dtype=np.uint64))
            handle.create_dataset("feature_idx", data=np.array([0, 1], dtype=np.uint32))
        else:
            group = handle.create_group("matrix" if kind == "modern" else "mm10")
            for key, data in {"data": [1, 2], "indices": [0, 1], "indptr": [0, 1, 2],
                              "shape": [2, 2]}.items():
                group.create_dataset(key, data=data)
            group.create_dataset("barcodes", data=[b"cell1", b"cell2"])
            if kind == "modern":
                features = group.create_group("features")
                for key, data in {"name": [b"G1", b"G2"], "id": [b"G1", b"G2"],
                                  "feature_type": [b"Gene Expression"] * 2}.items():
                    features.create_dataset(key, data=data)
            else:
                group.create_dataset("genes", data=[b"G1", b"G2"])
                group.create_dataset("gene_names", data=[b"G1", b"G2"])


@pytest.mark.parametrize("filename", ["molecule_info.h5", "renamed_upload.h5"])
def test_upload_rejects_molecule_info_before_job_creation(tmp_path, filename):
    path = tmp_path / filename
    _write_h5(path, "molecule")
    root = tmp_path / "jobs"
    app = create_app({"JOB_STORAGE": str(root)})
    with TestClient(app) as client:
        response = client.post("/api/jobs", data={"species": "mouse", "reference": "marrow",
                                                  "sample_names": "sample1"},
                               files={"files": (filename, path.read_bytes(), "application/octet-stream")})
    assert response.status_code == 400, response.text
    message = response.json()["detail"]
    assert filename in message and "molecule_info.h5" in message
    assert "filtered_feature_bc_matrix.h5" in message and "H5AD" in message
    assert not list(root.glob("*/job.json"))


@pytest.mark.parametrize("kind", ["modern", "legacy"])
def test_valid_matrix_upload_preserves_bytes(tmp_path, kind):
    path = tmp_path / "matrix.h5"
    _write_h5(path, kind)
    app = create_app({"JOB_STORAGE": str(tmp_path / "jobs")})
    with TestClient(app) as client:
        response = client.post("/api/jobs", data={"species": "mouse", "reference": "marrow",
                                                  "sample_names": "sample1"},
                               files={"files": (path.name, path.read_bytes(), "application/octet-stream")})
    assert response.status_code == 200, response.text
    store = app.state.job_store
    job = response.json()["job_id"]
    saved = store.uploads_dir(job) / store.get_job(job)["files"][0]["filename"]
    assert saved.read_bytes() == path.read_bytes()


def test_layout_check_reads_no_arrays_and_preserves_upload_position(tmp_path, monkeypatch):
    path = tmp_path / "molecule.h5"
    _write_h5(path, "molecule")
    upload = io.BytesIO(path.read_bytes())
    upload.seek(7)
    def forbidden(*args, **kwargs):
        raise AssertionError("validation must not read HDF5 arrays")
    monkeypatch.setattr(h5py.Dataset, "__getitem__", forbidden)
    with pytest.raises(ValueError, match="per-molecule"):
        validate_10x_h5(upload, "renamed.h5")
    assert upload.tell() == 7
    with pytest.raises(ValueError, match="filtered_feature_bc_matrix"):
        concat_matching_10x([str(path)], forbidden)


@pytest.mark.parametrize("stream", [False, True])
def test_existing_jobs_and_cli_get_same_error_before_scanpy(tmp_path, monkeypatch, stream):
    from altanalyze3.components.cellHarmony import cellHarmony_lite as chl

    path = tmp_path / "renamed.h5"
    _write_h5(path, "molecule")
    def forbidden(*args, **kwargs):
        raise AssertionError("Scanpy must not interpret per-molecule files as matrices")
    monkeypatch.setattr(chl.sc, "read_10x_h5", forbidden)
    reference = tmp_path / "reference.tsv"
    reference.write_text("gene\tstate1\tstate2\nG1\t1\t2\nG2\t2\t1\n")
    with pytest.raises(ValueError, match="molecule_info.h5.*filtered_feature_bc_matrix.h5"):
        chl.combine_and_align_h5(h5_files=[str(path)], cellharmony_ref=str(reference),
                                 output_dir=str(tmp_path / "output"), stream_10x_inputs=stream,
                                 concat_on_disk=stream)
