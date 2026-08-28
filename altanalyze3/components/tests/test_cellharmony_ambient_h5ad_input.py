#!/usr/bin/env python3
"""--ambient_correct_cutoff must apply to h5ad input, not only to multi-h5 input.

The ambient correction block used to sit inside the multi-h5 branch of
combine_and_align_h5. A run started from a single --h5ad accepted the option,
never applied it, printed no warning, and returned uncorrected counts.
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from altanalyze3.components.cellHarmony import cellHarmony_lite as chl


def _tiny_dataset(tmp_path, n_cells=400, n_genes=150):
    rng = np.random.default_rng(0)
    depth = rng.integers(3000, 12000, size=n_cells)
    profile = rng.dirichlet(np.ones(n_genes) * 0.4)
    X = np.vstack([rng.multinomial(int(d), profile) for d in depth]).astype(np.float32)
    genes = [f"Gene{i}" for i in range(n_genes)]
    obs = pd.DataFrame({"Library": ["L1"] * (n_cells // 2) + ["L2"] * (n_cells - n_cells // 2)},
                       index=[f"C{i}" for i in range(n_cells)])
    adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs,
                       var=pd.DataFrame(index=pd.Index(genes)))
    h5ad_path = os.path.join(tmp_path, "query.h5ad")
    adata.write_h5ad(h5ad_path)

    ref = pd.DataFrame(
        rng.random((40, 3)) * 5,
        index=pd.Index(genes[:40], name="uid"),
        columns=["StateA", "StateB", "StateC"],
    )
    ref_path = os.path.join(tmp_path, "ref.txt")
    ref.to_csv(ref_path, sep="\t")
    return h5ad_path, ref_path


def _run(tmp_path, rho, capsys):
    h5ad_path, ref_path = _tiny_dataset(str(tmp_path))
    outdir = os.path.join(str(tmp_path), f"out_{rho}")
    chl.combine_and_align_h5(
        h5_files=[], cellharmony_ref=ref_path, h5ad_file=h5ad_path,
        output_dir=outdir, min_genes=0, min_cells=0, min_counts=0, mit_percent=100,
        generate_umap=False, save_adata=False, export_h5ad=False, export_cptt=False,
        ambient_correct_cutoff=rho, return_adata=False,
    )
    return capsys.readouterr().out


def test_ambient_correction_runs_for_h5ad_input(tmp_path, capsys):
    out = _run(tmp_path, "0.2", capsys)
    assert "Running ambient RNA correction (rho=0.2)" in out, out[-2000:]
    assert "Ambient RNA correction complete" in out


def test_no_ambient_correction_when_not_requested(tmp_path, capsys):
    out = _run(tmp_path, None, capsys)
    assert "Running ambient RNA correction" not in out
