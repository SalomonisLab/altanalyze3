import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from altanalyze3.components.clustering.ICGS import (
    ICGS3Config,
    apply_expression_batch_adjustment,
    compute_sparse_graph,
    _batch_keys,
    _expression_batch_keys,
    legacy_rnaseq_gene_filter,
    read_inputs,
    run_icgs3,
)


def _write_sparse_h5ad(path):
    rng = np.random.default_rng(7)
    n_cells = 90
    n_genes = 45
    groups = np.repeat(["A", "B", "C"], n_cells // 3)
    X = rng.poisson(0.12, size=(n_cells, n_genes)).astype(float)
    X[:30, :5] += rng.poisson(3.0, size=(30, 5))
    X[30:60, 5:10] += rng.poisson(3.0, size=(30, 5))
    X[60:, 10:15] += rng.poisson(3.0, size=(30, 5))
    adata = AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame({"truth": groups}, index=[f"cell{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(n_genes)]),
    )
    adata.write_h5ad(path)


def test_ICGS_sparse_pipeline(tmp_path):
    input_path = tmp_path / "synthetic.h5ad"
    outdir = tmp_path / "icgs3"
    _write_sparse_h5ad(input_path)

    result = run_icgs3(
        ICGS3Config(
            input_paths=[str(input_path)],
            output_dir=str(outdir),
            input_normalized=False,
            min_genes=0,
            min_cells=0,
            min_counts=0,
            mito_percent=None,
            target_cells=None,
            pagerank_cells=60,
            louvain_downsample_cutoff=1000,
            downsample_key="truth",
            n_top_features=0,
            n_pcs=10,
            n_neighbors=8,
            marker_top_n=5,
            marker_min_per_cluster=1,
            rank=3,
            species="Hs",
            generate_umap=False,
            write_h5ad=False,
            heatmap_cells_per_cluster=10,
        )
    )

    assert result.adata.n_obs <= 90
    assert sp.issparse(result.adata.X)
    assert "ICGS3_cluster" in result.adata.obs
    assert not result.clusters.empty
    assert (outdir / "icgs3_clusters.tsv").exists()
    assert (outdir / "icgs3_cell_barcode_clusters.tsv").exists()
    cluster_df = pd.read_csv(outdir / "icgs3_cell_barcode_clusters.tsv", sep="\t")
    assert "ICGS3_SVM_score" in cluster_df.columns
    assert "ICGS3_original_NMF_cluster" in cluster_df.columns
    assert "ICGS3_cell_state_prediction" in cluster_df.columns
    for _, group in cluster_df.groupby("ICGS3_cluster", sort=False):
        scores = group["ICGS3_SVM_score"].to_numpy(dtype=float)
        assert np.all(scores[:-1] >= scores[1:])
    assert (outdir / "MarkerFinder" / "icgs3_markers.tsv").exists()
    assert (outdir / "sNMF" / "icgs3_nmf_variable_features.tsv").exists()
    assert (outdir / "sNMF" / "icgs3_nmf_rank_selection.tsv").exists()
    assert (outdir / "sNMF" / "icgs3_svm_reclassification_scores.final.tsv").exists()
    assert (outdir / "sNMF" / "icgs3_pagerank_downsampling.tsv").exists()
    assert (outdir / "MarkerFinder" / "icgs3_marker_heatmap.pdf").exists()
    assert (outdir / "MarkerFinder" / "icgs3_marker_heatmap.svg").exists()
    assert (outdir / "MarkerFinder" / "icgs3_marker_heatmap_fold_matrix.npz").exists()
    assert not (outdir / "MarkerFinder" / "icgs3_marker_heatmap_exp_matrix.tsv").exists()
    logs = list((outdir / "logs").glob("icgs3_*.log"))
    assert logs
    log_text = logs[0].read_text()
    assert "cli equivalent:" in log_text
    assert "NMF variable features:" in log_text
    assert not result.markers.empty


def test_legacy_rnaseq_gene_filter_exclusions():
    genes = ["ACTB", "RPL10", "RPS3", "MT-CO1", "A.B", "Gm123", "XIST", "TSIX", "RSPH1", "HLA-A", "DDX3Y", "CD34"]
    adata = AnnData(
        X=sp.csr_matrix(np.ones((3, len(genes)))),
        var=pd.DataFrame({"gene_symbols": genes}, index=genes),
    )
    keep = legacy_rnaseq_gene_filter(adata, protein_coding_genes=set(genes))
    kept = set(np.array(genes)[keep])
    assert kept == {"ACTB", "CD34"}


def test_read_single_h5ad_preserves_existing_sample_metadata(tmp_path):
    input_path = tmp_path / "metadata.h5ad"
    adata = AnnData(
        X=sp.csr_matrix(np.ones((4, 3))),
        obs=pd.DataFrame(
            {
                "Library": ["L1", "L1", "L2", "L2"],
                "sample": ["S1", "S1", "S2", "S2"],
            },
            index=["cell1", "cell2", "cell3", "cell4"],
        ),
        var=pd.DataFrame(index=["gene1", "gene2", "gene3"]),
    )
    adata.write_h5ad(input_path)
    loaded = read_inputs([str(input_path)])
    assert loaded.obs["Library"].astype(str).tolist() == ["L1", "L1", "L2", "L2"]
    assert loaded.obs["sample"].astype(str).tolist() == ["S1", "S1", "S2", "S2"]
    assert loaded.obs_names.tolist() == ["cell1", "cell2", "cell3", "cell4"]


def test_harmony_batch_correction_replaces_graph_pcs():
    try:
        import harmonypy  # noqa: F401
    except Exception:
        return
    rng = np.random.default_rng(9)
    n_cells = 60
    n_genes = 20
    batch = np.repeat(["ADT112", "ADT195"], n_cells // 2)
    X = rng.normal(0, 1, size=(n_cells, n_genes))
    X[batch == "ADT195", :6] += 2.0
    adata = AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame({"ADT_panel": batch}, index=[f"cell{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"feature{i}" for i in range(n_genes)]),
    )
    graph = compute_sparse_graph(
        adata,
        ICGS3Config(
            input_paths=["dummy.h5ad"],
            output_dir="/tmp/icgs3_test",
            batch_correction="harmony",
            batch_key="ADT_panel",
            batch_correction_use="graph",
            n_pcs=8,
            n_neighbors=8,
        ),
    )
    assert "X_pca_harmony" in graph.obsm
    assert "X_pca_uncorrected" in graph.obsm
    assert graph.uns["icgs3_graph"]["representation"] == "X_pca_harmony"
    assert "connectivities" in graph.obsp


def test_expression_batch_adjustment_is_nonnegative_and_batch_centered(tmp_path):
    X = np.array(
        [
            [10.0, 1.0],
            [12.0, 2.0],
            [2.0, 9.0],
            [3.0, 11.0],
        ]
    )
    adata = AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame({"panel": ["A", "A", "B", "B"]}, index=[f"cell{i}" for i in range(4)]),
        var=pd.DataFrame(index=["CD1", "CD2"]),
    )
    outdir = tmp_path / "icgs3"
    (outdir / "sNMF").mkdir(parents=True)
    adjusted = apply_expression_batch_adjustment(
        adata,
        ICGS3Config(
            input_paths=["dummy.h5ad"],
            output_dir=str(outdir),
            batch_adjust_expression=True,
            expression_batch_key="panel",
        ),
        str(outdir),
    )
    dense = adjusted.X.toarray()
    assert np.all(dense >= 0)
    assert np.allclose(dense[:2].mean(axis=0), dense[2:].mean(axis=0))
    assert not np.allclose(adata.X.toarray(), dense)
    assert (outdir / "sNMF" / "icgs3_expression_batch_adjustment_groups.tsv").exists()


def test_batch_key_resolvers_combine_generic_technology_and_tissue_keys():
    config = ICGS3Config(
        input_paths=["dummy.h5ad"],
        output_dir="/tmp/icgs3_test",
        batch_key="Library",
        technology_batch_key="ADT_panel",
        tissue_batch_key="TissueGroup",
    )
    assert _batch_keys(config) == ["Library", "ADT_panel", "TissueGroup"]
    assert _expression_batch_keys(config) == ["Library", "ADT_panel", "TissueGroup"]

    config.expression_batch_key = "ADT_panel"
    assert _expression_batch_keys(config) == ["ADT_panel", "TissueGroup"]


def test_nmf_rowsum_assignment_normalization_balances_component_scale():
    from altanalyze3.components.clustering.ICGS import _add_udon_path

    _add_udon_path()
    from nmf import binarize_nmf

    # Shape mirrors nimfa basis output: samples x components. Raw argmax is
    # dominated by the first component scale; rowsum normalization recovers the
    # component-specific winners.
    w = np.array(
        [
            [100.0, 1.0, 1.0],
            [90.0, 20.0, 1.0],
            [80.0, 1.0, 30.0],
        ]
    )
    raw = binarize_nmf(w)["cluster"].tolist()
    rowsum = binarize_nmf(w, normalization="rowsum")["cluster"].tolist()
    assert raw == [0, 0, 0]
    assert rowsum == [0, 1, 2]
