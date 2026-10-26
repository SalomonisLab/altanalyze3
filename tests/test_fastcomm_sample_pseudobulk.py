"""fastComm per-sample scores -> sample x receiver-state pseudobulk (sample-pseudobulk)."""
import json
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.fastComm.sample_pseudobulk import build_sample_pseudobulk, write_covariates


def _long(tmp_path, rows):
    p = tmp_path / "split_scores_long.tsv"
    pd.DataFrame(rows, columns=["split", "sender_state", "receiver_state", "ligand", "receptor", "fastcomm_score",
                                "pathway"]).to_csv(p, sep="\t", index=False)
    return p


def test_shape_values_and_duplicate_collapse(tmp_path):
    p = _long(tmp_path, [
        ["D1", "MSC", "HSC-1", "CXCL12", "CXCR4", 0.8, "CXCL"],
        ["D1", "MSC", "HSC-1", "CXCL12", "CXCR4", 0.6, "CXCL"],
        ["D1", "Mono", "HSC-1", "TNF", "TNFRSF1A", 0.3, "TNF"],
        ["D2", "MSC", "HSC-1", "CXCL12", "CXCR4", 0.5, "CXCL"],
        ["D2", "MSC", "MPP-1", "KITLG", "KIT", 0.4, "KIT"],
        ["D3", "MSC", "MPP-1", "KITLG", "KIT", 0.0, "KIT"],
    ])
    out = tmp_path / "pb.h5ad"
    a, s = build_sample_pseudobulk(p, out, log=lambda *x: None)
    b = ad.read_h5ad(out)
    assert list(b.obs_names) == ["D1|HSC-1", "D2|HSC-1", "D2|MPP-1", "D3|MPP-1"]
    assert list(b.var_names) == ["MSC|CXCL12->CXCR4", "MSC|KITLG->KIT", "Mono|TNF->TNFRSF1A"]
    X = b.X.toarray()
    assert np.allclose(X, [[0.8, 0, 0.3], [0.5, 0, 0], [0, 0.4, 0], [0, 0, 0]])
    assert b.uns["pseudobulk_method"] == "pseudobulk" and b.uns["value"] == "fastcomm_score"
    assert list(b.obs.Sample) == ["D1", "D2", "D2", "D3"] and list(b.obs.cell_type) == ["HSC-1", "HSC-1", "MPP-1", "MPP-1"]
    assert s["records"] == 6 and s["collapsed_duplicates"] == 1 and s["samples"] == 3 and s["nnz"] == 4


def test_rejects_pipe_in_keys_and_bad_scores(tmp_path):
    p = _long(tmp_path, [["D|1", "MSC", "HSC-1", "CXCL12", "CXCR4", 0.8, "CXCL"]])
    with pytest.raises(ValueError, match="contain"):
        build_sample_pseudobulk(p, tmp_path / "x.h5ad", log=lambda *x: None)
    p = _long(tmp_path, [["D1", "MSC", "HSC-1", "CXCL12", "CXCR4", "high", "CXCL"]])
    with pytest.raises(ValueError, match="non-numeric"):
        build_sample_pseudobulk(p, tmp_path / "x.h5ad", log=lambda *x: None)


def test_covariates_require_every_sample(tmp_path):
    cov = tmp_path / "cov.tsv"
    pd.DataFrame({"Sample": ["D1", "D2", "D9"], "Condition": ["aged", "young", "young"]}).to_csv(cov, sep="\t", index=False)
    out = write_covariates(["D1", "D2", "D1"], cov, tmp_path / "c.txt", log=lambda *x: None)
    assert list(out.Sample) == ["D1", "D2"]
    with pytest.raises(ValueError, match="without a covariate"):
        write_covariates(["D1", "D3"], cov, tmp_path / "c2.txt", log=lambda *x: None)


def test_cli_subcommand(tmp_path):
    p = _long(tmp_path, [["D1", "MSC", "HSC-1", "CXCL12", "CXCR4", 0.8, "CXCL"],
                         ["D2", "MSC", "HSC-1", "CXCL12", "CXCR4", 0.5, "CXCL"]])
    cov = tmp_path / "cov.tsv"
    pd.DataFrame({"Sample": ["D1", "D2"], "Condition": ["aged", "young"]}).to_csv(cov, sep="\t", index=False)
    out = tmp_path / "pb.h5ad"
    r = subprocess.run([sys.executable, "-m", "altanalyze3.components.fastComm.cli", "sample-pseudobulk",
                        "--split-scores", str(p), "--output", str(out), "--covariates", str(cov)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert out.exists() and (tmp_path / "pb_covariates.txt").exists()
    assert json.load(open(tmp_path / "pb.summary.json"))["rows"] == 2
