"""Prove h5ad_group_obs maps the named labels, sends the rest to Other, and rejects a bad mapping."""

import os
import tempfile

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")
sparse = pytest.importorskip("scipy.sparse")

from altanalyze3.components.aggregate.h5ad_group_obs import group_obs, read_mapping

MAPPING = """out_key\tsource_key\tgroup_name\tsource_label
Ref_group\tRef\tCLP-1\tCLP1-a
Ref_group\tRef\tCLP-1\tCLP1-b
Ref_group\tRef\tEILP\tEILP
"""


def _source(tmp, labels):
    path = os.path.join(tmp, "src.h5ad")
    anndata.AnnData(
        X=sparse.csr_matrix(np.zeros((len(labels), 3), dtype=np.float32)),
        obs=pd.DataFrame({"Ref": labels}, index=[f"BC{i:03d}" for i in range(len(labels))]),
        var=pd.DataFrame(index=["a", "b", "c"]),
    ).write(path)
    return path


def _mapping(tmp, text=MAPPING):
    path = os.path.join(tmp, "map.tsv")
    open(path, "w").write(text)
    return path


def test_named_labels_group_and_the_rest_fall_to_other():
    with tempfile.TemporaryDirectory() as tmp:
        labels = ["CLP1-a"] * 5 + ["CLP1-b"] * 3 + ["EILP"] * 4 + ["ML-1a"] * 8
        src = _source(tmp, labels)
        out = os.path.join(tmp, "side.h5ad")
        counts = group_obs(src, _mapping(tmp), out, log=lambda *_: None)

        side = anndata.read_h5ad(out)
        assert side.n_obs == 20
        got = side.obs["Ref_group"].astype(str)
        assert list(got[:5]) == ["CLP-1"] * 5
        assert list(got[5:8]) == ["CLP-1"] * 3
        assert list(got[8:12]) == ["EILP"] * 4
        assert list(got[12:]) == ["Other"] * 8
        assert got.value_counts().to_dict() == {"CLP-1": 8, "EILP": 4, "Other": 8}
        # the category order follows the mapping file, with Other last
        assert list(side.obs["Ref_group"].cat.categories) == ["CLP-1", "EILP", "Other"]
        # the source column survives, so one sidecar serves every covariate
        assert "Ref" in side.obs.columns
        assert int(counts["cells"].sum()) == 20
        assert os.path.exists(os.path.join(tmp, "side_group_counts.tsv"))


def test_a_label_the_column_does_not_hold_fails_by_default():
    with tempfile.TemporaryDirectory() as tmp:
        src = _source(tmp, ["EILP"] * 4 + ["ML-1a"] * 4)   # no CLP1-a, no CLP1-b
        with pytest.raises(ValueError, match="holds no cell labelled"):
            group_obs(src, _mapping(tmp), os.path.join(tmp, "o.h5ad"), log=lambda *_: None)


def test_allow_missing_labels_keeps_going_and_drops_the_empty_group():
    with tempfile.TemporaryDirectory() as tmp:
        src = _source(tmp, ["EILP"] * 4 + ["ML-1a"] * 4)
        out = os.path.join(tmp, "o.h5ad")
        group_obs(src, _mapping(tmp), out, allow_missing_labels=True, log=lambda *_: None)
        side = anndata.read_h5ad(out)
        assert list(side.obs["Ref_group"].cat.categories) == ["EILP", "Other"]
        assert side.obs["Ref_group"].value_counts().to_dict() == {"EILP": 4, "Other": 4}


def test_one_label_in_two_groups_is_rejected():
    with tempfile.TemporaryDirectory() as tmp:
        bad = MAPPING + "Ref_group\tRef\tEILP-2\tEILP\n"
        with pytest.raises(ValueError, match="two groups"):
            read_mapping(_mapping(tmp, bad))


def test_one_out_key_reading_two_source_columns_is_rejected():
    with tempfile.TemporaryDirectory() as tmp:
        bad = MAPPING + "Ref_group\tOtherCol\tX\tY\n"
        with pytest.raises(ValueError, match="more than one source_key"):
            read_mapping(_mapping(tmp, bad))


def test_missing_obs_column_is_rejected():
    with tempfile.TemporaryDirectory() as tmp:
        src = _source(tmp, ["EILP"] * 4)
        bad = "out_key\tsource_key\tgroup_name\tsource_label\nG\tNope\tA\tEILP\n"
        with pytest.raises(KeyError):
            group_obs(src, _mapping(tmp, bad), os.path.join(tmp, "o.h5ad"), log=lambda *_: None)
