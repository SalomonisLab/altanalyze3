import anndata as ad
import numpy as np
import pandas as pd
import pytest

fitz = pytest.importorskip("fitz")

from altanalyze3.components.visualization import violin


def _toy(path):
    rng = np.random.default_rng(1)
    obs = pd.DataFrame({"state": np.repeat(["A", "B"], 6),
                        "age_group": pd.Categorical(np.tile(["young", "aged"], 6), categories=["young", "aged"])},
                       index=[f"c{i}" for i in range(12)])
    ad.AnnData(X=rng.normal(size=(12, 1)).astype(np.float32), obs=obs,
               var=pd.DataFrame(index=["g1"])).write_h5ad(path)


def _text(pdf):
    return "".join(page.get_text() for page in fitz.open(pdf))


def test_default_ylabel_is_expression(tmp_path):
    _toy(tmp_path / "toy.h5ad")
    out = violin.write_lineage_violin_pdf(tmp_path / "toy.h5ad", output_pdf=tmp_path / "d.pdf", genes=["g1"],
                                          cluster_key="state", mode="population-obs", population="A",
                                          compare_obs_field="age_group")
    assert "Expression" in _text(out)


def test_ylabel_replaces_only_the_axis_label(tmp_path):
    _toy(tmp_path / "toy.h5ad")
    a = _text(violin.write_lineage_violin_pdf(tmp_path / "toy.h5ad", output_pdf=tmp_path / "d.pdf", genes=["g1"],
                                              cluster_key="state", mode="population-obs", population="A",
                                              compare_obs_field="age_group"))
    b = _text(violin.write_lineage_violin_pdf(tmp_path / "toy.h5ad", output_pdf=tmp_path / "u.pdf", genes=["g1"],
                                              cluster_key="state", mode="population-obs", population="A",
                                              compare_obs_field="age_group", ylabel="Isoform usage"))
    assert "Expression" not in b and "Isoform usage" in b
    assert a.replace("Expression", "") == b.replace("Isoform usage", "")
