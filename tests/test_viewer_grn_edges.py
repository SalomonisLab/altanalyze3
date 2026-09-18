"""scALABLE-viewer: GRN edges served from the bundle's own per-cell-state store.

`webapp/app.py` `_grn_edges_adata` opens
`meta['modality_artifacts']['grn']['network_h5ad']`, the edge-level object
cellHarmony-differential writes per contrast. A bundle has none, so the GRN edges panel
answered HTTP 404 "GRN edge output is unavailable for this job." while the edges sat in
the bundle's `grn` store.
"""
from types import SimpleNamespace

import numpy as np
import pytest

from altanalyze3.components.visualization.scalable_viewer import scalable_app as SA


EDGES = ["SPI1|SLC2A9", "SPI1|ZFYVE26", "GATA1|HBB"]
STATES = ["AF", "AM", "AT1", "AT2"]
# features x states, the layout `<prefix>_grn_stats_mean.npy` uses.
MEAN = np.array([[3.5, 0.1, 0.0, 2.0],
                 [2.3, 0.2, 1.0, 0.5],
                 [0.0, 4.4, 0.3, 0.1]], dtype=np.float32)


def fake_ds(mean=MEAN, edges=EDGES, states=STATES):
    store = SimpleNamespace(stats_mean=mean, features=list(edges))
    return SimpleNamespace(states=list(states), modality=lambda name: store)


def test_adata_is_cell_states_by_edges():
    a = SA._bundle_grn_edge_adata(fake_ds(), "cell_type")
    assert a.shape == (len(STATES), len(EDGES))          # transposed from the store
    assert list(a.var_names) == EDGES
    assert list(a.obs["cell_type"]) == STATES


def test_each_state_row_holds_that_states_stored_scores():
    """The defect this guards: a transpose or reorder would attribute another state."""
    a = SA._bundle_grn_edge_adata(fake_ds(), "cell_type")
    for col, state in enumerate(STATES):
        row = np.asarray(a[a.obs["cell_type"] == state].X).ravel()
        assert np.allclose(row, MEAN[:, col]), f"{state} row does not match stored column"


def test_cluster_key_names_the_obs_column():
    a = SA._bundle_grn_edge_adata(fake_ds(), "Population")
    assert "Population" in a.obs.columns and "cell_type" not in a.obs.columns
    # A blank key still yields a usable default rather than an unnamed column.
    assert "cell_type" in SA._bundle_grn_edge_adata(fake_ds(), "").obs.columns


def test_shape_disagreement_raises_rather_than_mislabelling():
    bad = fake_ds(mean=np.zeros((2, len(STATES)), dtype=np.float32))   # 2 rows, 3 names
    with pytest.raises(ValueError, match="GRN store is"):
        SA._bundle_grn_edge_adata(bad, "cell_type")


def test_uploaded_job_still_uses_the_original_reader():
    """A cellHarmony job has no `scalable_viewer` block and must not be diverted."""
    calls = []
    original = SA.W._grn_edges_adata
    try:
        SA.W._grn_edges_adata = lambda meta: calls.append(meta) or "original"
        SA._install_grn_edge_adata(SimpleNamespace(state=SimpleNamespace(catalog=None)),
                                   SimpleNamespace())
        assert SA.W._grn_edges_adata({"job": "upload"}) == "original"
        assert calls == [{"job": "upload"}]
    finally:
        SA.W._grn_edges_adata = original


# --- GRN edge detail: sidecar edge set vs ragged edge set -------------------------

import numpy as np
from importlib import import_module

# `from ...webapp import app` binds the FastAPI instance the package __init__ re-exports,
# not the module. Import the module explicitly, as bundle_meta.py does.
WEBAPP = import_module("altanalyze3.components.cellHarmony.webapp.app")


def test_ragged_reader_returns_none_for_an_edge_it_lacks(tmp_path):
    """The guard the detail panel now depends on: absent edge -> None, never a fabrication."""
    meta = {"grn_ragged_dir": str(tmp_path)}
    # No store directory entry at all.
    assert WEBAPP._grn_ragged_values(meta, "AT1", "KLF6|NTN4", ["a", "b"]) is None
    # Blank arguments stay None rather than reading an arbitrary file.
    assert WEBAPP._grn_ragged_values(meta, "", "KLF6|NTN4", ["a"]) is None
    assert WEBAPP._grn_ragged_values(meta, "AT1", "", ["a"]) is None
    assert WEBAPP._grn_ragged_values({"grn_ragged_dir": ""}, "AT1", "E", ["a"]) is None


def test_ragged_reader_aligns_values_to_the_requested_barcodes(tmp_path):
    """A barcode the ragged store lacks must stay NaN, not borrow a neighbour's value."""
    import anndata as ad
    import pandas as pd

    x = np.array([[1.5], [2.5]], dtype=np.float32)
    a = ad.AnnData(x, obs=pd.DataFrame(index=["bc1", "bc2"]),
                   var=pd.DataFrame(index=["KLF6|NTN4"]))
    a.write_h5ad(tmp_path / "AT1.h5ad")
    out = WEBAPP._grn_ragged_values({"grn_ragged_dir": str(tmp_path)}, "AT1",
                                    "KLF6|NTN4", ["bc2", "missing", "bc1"])
    assert out is not None
    assert out[0] == 2.5 and out[2] == 1.5      # aligned by barcode, not by position
    assert np.isnan(out[1])                      # unknown barcode is not filled in


def test_population_name_is_sanitised_into_the_store_filename(tmp_path):
    """Cell states carry spaces and slashes; the reader must find the same file."""
    import anndata as ad
    import pandas as pd

    a = ad.AnnData(np.array([[3.0]], dtype=np.float32),
                   obs=pd.DataFrame(index=["bc1"]), var=pd.DataFrame(index=["E|T"]))
    a.write_h5ad(tmp_path / "AT2-AT1_int..h5ad")
    out = WEBAPP._grn_ragged_values({"grn_ragged_dir": str(tmp_path)},
                                    "AT2-AT1 int.", "E|T", ["bc1"])
    assert out is not None and out[0] == 3.0
