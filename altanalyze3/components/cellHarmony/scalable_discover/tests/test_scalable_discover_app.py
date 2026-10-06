"""Fast checks of the scALABLE-discover app; the end-to-end runs live in validation/."""
import io
from pathlib import Path
import tempfile

import anndata as ad
import numpy as np
import pytest
from fastapi.testclient import TestClient

from altanalyze3.components.cellHarmony.scalable_discover import app as discover
from altanalyze3.components.cellHarmony.scalable_discover.pipeline import discover_registry


@pytest.fixture()
def client(tmp_path):
    application = discover.create_discover_app({"JOB_STORAGE": str(tmp_path), "ISOLATE_JOBS": False})
    return TestClient(application)


def test_page_is_rewritten_from_the_shared_template(client):
    page = client.get("/").text
    assert "<h1>scALABLE-discover</h1>" in page
    assert 'data-tab="differential"' not in page
    assert '<option value="relative">UMAP broad</option>' not in page
    assert "discover-static/discover.js" in page
    assert '"references": [{"id": "icgs3"' in page


def test_rewrite_fails_loudly_when_an_anchor_moves():
    with pytest.raises(RuntimeError, match="expected 1"):
        discover.rewrite_index("<html><body></body></html>", "", 0)


def test_species_menu_is_human_and_mouse_only(client):
    assert client.get("/api/meta/species").json() == discover_registry()
    assert [s["id"] for s in discover_registry()["species"]] == ["human", "mouse"]


@pytest.mark.parametrize("path", ["/api/jobs/x/differential/status", "/api/jobs/x/differential/interactive/summary",
                                  "/api/meta/reference-preview"])
def test_differential_and_reference_routes_are_gone(client, path):
    assert client.get(path).status_code == 404


def _upload(client, species="human", reference="icgs3"):
    # Upload validation inspects the H5AD header before accepting a job.
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "s1.h5ad"
        ad.AnnData(np.array([[1, 2], [3, 4]], dtype=np.float32)).write_h5ad(path)
        payload = path.read_bytes()
    return client.post("/api/jobs", data={"species": species, "reference": reference, "sample_names": ["s1"]},
                       files={"files": ("s1.h5ad", io.BytesIO(payload), "application/octet-stream")})


def test_upload_refuses_other_species_and_references(client):
    assert _upload(client, species="zebrafish").status_code == 400
    assert _upload(client, reference="hs_lung_cellref2_reference").status_code == 400


def test_qc_refuses_imputation_and_stores_no_alignment_cutoff(client):
    job = _upload(client).json()["job_id"]
    assert client.post(f"/api/jobs/{job}/qc", json={"impute_modalities": ["adt"]}).status_code == 400
    saved = client.post(f"/api/jobs/{job}/qc", json={"min_genes": 300, "ambient_correction": "yes"}).json()["qc"]
    assert saved["min_genes"] == 300 and saved["impute_modalities"] == [] and "align_cutoff" not in saved


def test_runner_is_the_discover_runner(client):
    runner = client.app.state.job_runner
    assert type(runner).__name__ == "DiscoverJobRunner"
    assert runner.WORKER_MODULE.endswith("scalable_discover.worker")
