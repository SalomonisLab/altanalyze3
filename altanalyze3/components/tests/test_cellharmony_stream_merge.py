"""Streaming 10x merge preserves the normal loader and concatenation contract."""
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp

from altanalyze3.components.cellHarmony.merge_inputs import concat_matching_10x
from altanalyze3.components.cellHarmony import cellHarmony_lite as chl


def _write_10x(path, matrix, names, feature_types=None):
    matrix = sp.csr_matrix(matrix, dtype=np.int32)
    with h5py.File(path, "w") as handle:
        group = handle.create_group("matrix")
        for key, values in dict(data=matrix.data, indices=matrix.indices, indptr=matrix.indptr).items():
            group.create_dataset(key, data=values)
        group.create_dataset("shape", data=[matrix.shape[1], matrix.shape[0]])
        group.create_dataset("barcodes", data=[f"c{i}".encode() for i in range(matrix.shape[0])])
        features = group.create_group("features")
        for key, values in dict(name=[g.encode() for g in names], id=[g.encode() for g in names],
                                feature_type=feature_types or [b"Gene Expression"] * len(names),
                                genome=[b"GRCh38"] * len(names)).items():
            features.create_dataset(key, data=values)


def _load(path, sample_name_override=None):
    sample = sc.read_10x_h5(path)
    name = sample_name_override or Path(path).stem
    sample.var_names_make_unique()
    sample.obs_names = [f"{cell}.{name}" for cell in sample.obs_names]
    sample.obs["sample"] = name
    sample.obs["Library"] = name
    return sample, name


def test_streamed_merge_preserves_duplicate_gene_names_and_metadata(tmp_path):
    names = ["g0", "g0", "g2", "g3"]
    files = []
    for i in range(3):
        path = tmp_path / f"sample{i}.h5"
        _write_10x(path, [[0, 1, 3, i], [4, 0, 2, 0]], names)
        files.append((str(path), f"override{i}"))
    expected = ad.concat([_load(path, name)[0] for path, name in files], label="sample", join="outer", fill_value=0)
    actual = concat_matching_10x(files, _load)
    np.testing.assert_array_equal(actual.X.toarray(), expected.X.toarray())
    pd.testing.assert_frame_equal(actual.obs, expected.obs)
    pd.testing.assert_frame_equal(actual.var, expected.var)


@pytest.mark.parametrize("reason", ["different", "reordered", "multimodal"])
def test_preflight_falls_back_before_loading_expression(tmp_path, reason):
    names = ["g0", "g1", "g2", "g3"]
    left, right = tmp_path / "left.h5", tmp_path / "right.h5"
    _write_10x(left, np.ones((3, 4)), names)
    second_names = names.copy()
    types = None
    if reason == "different":
        second_names[2] = "different"
    elif reason == "reordered":
        second_names.reverse()
    else:
        types = [b"Gene Expression"] * 3 + [b"Antibody Capture"]
    _write_10x(right, np.ones((3, 4)), second_names, types)
    def unexpected_load(*args, **kwargs):
        raise AssertionError("preflight must precede loading expression")
    assert concat_matching_10x([str(left), str(right)], unexpected_load) is None


def test_stream_pipeline_matches_normal_import_after_qc(tmp_path):
    names = [f"g{i}" for i in range(12)]
    files = []
    rng = np.random.default_rng(3)
    for i in range(3):
        path = tmp_path / f"sample{i}.h5"
        matrix = rng.integers(0, 50, size=(30, 12), dtype=np.int32)
        matrix[0] = 0
        _write_10x(path, matrix, names)
        files.append(str(path))
    reference = tmp_path / "reference.tsv"
    pd.DataFrame(rng.random((12, 3)), index=names, columns=["a", "b", "c"]).to_csv(reference, sep="\t")
    outputs = []
    for stream in (False, True):
        outputs.append(chl.combine_and_align_h5(
            h5_files=files, cellharmony_ref=str(reference), output_dir=str(tmp_path / str(stream)),
            ambient_correct_cutoff="auto", ambient_memory_efficient=True,
            concat_on_disk=stream, concat_batch_size=1 if stream else None,
            stream_10x_inputs=stream, min_genes=2, min_counts=10, min_cells=0,
            mit_percent=100, generate_umap=False, export_h5ad=False, save_adata=False,
            return_adata=True,
        ))
    (expected_assignments, expected), (actual_assignments, actual) = outputs
    pd.testing.assert_frame_equal(actual_assignments, expected_assignments, atol=1e-7, rtol=1e-7)
    assert actual.n_obs == 87
    pd.testing.assert_index_equal(actual.obs_names, expected.obs_names)
    np.testing.assert_array_equal(actual.X.toarray(), expected.X.toarray())
    np.testing.assert_array_equal(actual.layers["soupx_raw"].toarray(), expected.layers["soupx_raw"].toarray())
