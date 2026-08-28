"""Prove build_multimodal_input attaches labels to the right barcodes and keeps only complete cells."""

import os
import tempfile

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")
sparse = pytest.importorskip("scipy.sparse")

from altanalyze3.components.sctriangulate.build_multimodal_input import build, read_annotation_table


def _fixture(tmp):
    rng = np.random.default_rng(0)
    barcodes = [f"BC{i:03d}-1.Lib" for i in range(60)]
    rna = anndata.AnnData(
        X=sparse.csr_matrix(rng.random((60, 30)).astype(np.float32)),
        obs=pd.DataFrame({"Ref_L4": ["StateA"] * 30 + ["StateB"] * 30}, index=barcodes),
        var=pd.DataFrame(index=[f"Gene{j}" for j in range(30)]),
    )
    rna_path = os.path.join(tmp, "rna.h5ad")
    rna.write(rna_path)

    # ADT covers the first 50 barcodes only
    adt = anndata.AnnData(
        X=(rng.random((50, 8)) * 200).astype(np.float32),
        obs=pd.DataFrame(index=barcodes[:50]),
        var=pd.DataFrame(index=[f"CD{j}" for j in range(8)]),
    )
    adt_path = os.path.join(tmp, "adt.h5ad")
    adt.write(adt_path)

    # a cluster table covering the first 55 barcodes, label column not second
    tsv = os.path.join(tmp, "clusters.tsv")
    pd.DataFrame({
        "barcode": barcodes[:55],
        "junk": ["x"] * 55,
        "ICGS3_cluster": [f"C{i % 3}" for i in range(55)],
    }).to_csv(tsv, sep="\t", index=False)
    return rna_path, adt_path, tsv, barcodes


def test_labels_land_on_the_right_barcodes_and_complete_cells_survive():
    with tempfile.TemporaryDirectory() as tmp:
        rna_path, adt_path, tsv, barcodes = _fixture(tmp)
        out = os.path.join(tmp, "combined.h5ad")
        report = build(
            rna_path, out, adt_h5ad=adt_path,
            file_annotations=[("ICGS3_RNA", tsv, "ICGS3_cluster")],
            obs_annotations=[("Ferchen", "Ref_L4")],
            require_complete=True, log=lambda *_: None,
        )
        combined = anndata.read_h5ad(out)
        # ADT covers 50, the TSV covers 55 -> the intersection is 50
        assert combined.n_obs == 50 == report["cells_kept"]
        assert combined.n_vars == 30 + 8
        assert sum(str(v).startswith("AB_") for v in combined.var_names) == 8
        assert list(combined.obs_names) == barcodes[:50]
        expected = [f"C{i % 3}" for i in range(50)]
        assert list(combined.obs["ICGS3_RNA"].astype(str)) == expected
        assert list(combined.obs["Ferchen"].astype(str)) == ["StateA"] * 30 + ["StateB"] * 20
        assert not combined.obs[["ICGS3_RNA", "Ferchen"]].isna().any().any()
        assert report["query"] == "ICGS3_RNA,Ferchen"
        assert os.path.exists(os.path.join(tmp, "combined_annotations.tsv"))
        assert os.path.exists(os.path.join(tmp, "combined_build_report.json"))


def test_adt_log1p_is_applied_before_concatenation():
    with tempfile.TemporaryDirectory() as tmp:
        rna_path, adt_path, tsv, _ = _fixture(tmp)
        raw = anndata.read_h5ad(adt_path)
        for mode, transform in [("log1p", np.log1p), ("none", lambda a: a)]:
            out = os.path.join(tmp, f"c_{mode}.h5ad")
            build(rna_path, out, adt_h5ad=adt_path, adt_normalization=mode,
                  obs_annotations=[("Ferchen", "Ref_L4")], require_complete=True,
                  log=lambda *_: None)
            combined = anndata.read_h5ad(out)
            ab = combined[:, [v for v in combined.var_names if str(v).startswith("AB_")]]
            got = ab.X.toarray() if sparse.issparse(ab.X) else np.asarray(ab.X)
            assert np.allclose(got, transform(np.asarray(raw.X)), atol=1e-5), mode


def test_read_annotation_table_defaults_to_the_second_column():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "two.tsv")
        pd.DataFrame({"CellBarcode": ["a", "b"], "Leiden": ["1", "2"]}).to_csv(
            path, sep="\t", index=False)
        s = read_annotation_table(path)
        assert s.to_dict() == {"a": "1", "b": "2"}


def test_duplicate_barcodes_are_rejected():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "dup.tsv")
        pd.DataFrame({"bc": ["a", "a"], "c": ["1", "2"]}).to_csv(path, sep="\t", index=False)
        with pytest.raises(ValueError, match="repeats"):
            read_annotation_table(path)


def test_missing_obs_column_is_rejected():
    with tempfile.TemporaryDirectory() as tmp:
        rna_path, _, _, _ = _fixture(tmp)
        with pytest.raises(KeyError):
            build(rna_path, os.path.join(tmp, "o.h5ad"),
                  obs_annotations=[("X", "not_a_column")], log=lambda *_: None)
