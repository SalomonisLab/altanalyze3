from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from altanalyze3.components.fastCNV.main import FastCNVParams, run_fastcnv


def test_fastcnv_detects_broad_gain(tmp_path: Path):
    rng = np.random.default_rng(7)
    n_cells = 120
    n_genes = 180
    genes = [f"G{i:03d}" for i in range(n_genes)]
    counts = rng.poisson(2.0, size=(n_cells, n_genes)).astype(np.float32)
    clone_cells = np.arange(90, 120)
    counts[clone_cells, 20:95] += rng.poisson(6.0, size=(len(clone_cells), 75)).astype(np.float32)

    obs = pd.DataFrame(
        {
            "cell_state": ["epithelial"] * n_cells,
            "sample": ["sample1"] * n_cells,
        },
        index=[f"cell{i:03d}" for i in range(n_cells)],
    )
    adata = ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=genes))
    adata.layers["counts"] = counts
    h5ad_path = tmp_path / "synthetic.h5ad"
    adata.write_h5ad(h5ad_path)

    coords = pd.DataFrame(
        {
            "gene": genes,
            "chr": ["chr1"] * 120 + ["chr2"] * 60,
            "start": np.arange(n_genes) * 1000 + 1,
            "end": np.arange(n_genes) * 1000 + 900,
        }
    )
    coords_path = tmp_path / "coords.tsv"
    coords.to_csv(coords_path, sep="\t", index=False)

    outputs = run_fastcnv(
        FastCNVParams(
            h5ad=h5ad_path,
            gene_coordinates=coords_path,
            output_prefix=tmp_path / "fastcnv",
            state_key="cell_state",
            sample_key="sample",
            window_genes=21,
            stride_genes=5,
            min_interval_genes=30,
            min_run_windows=2,
            high_threshold=2.2,
            low_threshold=1.4,
            min_mean_score=1.5,
            cnv_burden_threshold=1.4,
        )
    )

    cells = pd.read_csv(outputs["cells"], sep="\t")
    intervals = pd.read_csv(outputs["intervals"], sep="\t")
    clone_intervals = pd.read_csv(outputs["clone_intervals"], sep="\t")

    clone_status = cells.set_index("CellBarcode").loc[[f"cell{i:03d}" for i in clone_cells], "cnv_status"]
    assert (clone_status == "CNV").mean() >= 0.75
    assert "state_clone_id" in cells.columns
    assert "global_clone_id" in cells.columns
    assert set(cells.loc[cells["cnv_status"] == "CNV", "global_clone_id"]) - {"WT"}
    assert not intervals.empty
    assert not clone_intervals.empty
    assert "gain" in set(intervals["call"])
    assert "chr1" in set(intervals["chr"])
    assert outputs["clone_pdf"].exists()
