import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from altanalyze3.components.clustering.ICGS import (
    ADT_OCCUPANCY_KEY,
    ADT_OCCUPANCY_STATUS_KEY,
    ICGS3Config,
    apply_adt_occupancy_filter,
    apply_expression_batch_adjustment,
    apply_modality_defaults,
    compute_adt_occupancy,
    estimate_udon_nmf_rank,
    resolve_normalization,
    resolve_small_feature_rank,
    _trim_occupancy_within_groups,
    build_arg_parser,
    classify_adata_with_scores_chunked,
    classify_with_scores,
    cli_equivalent,
    compute_sparse_graph,
    _batch_keys,
    _expression_batch_keys,
    legacy_rnaseq_gene_filter,
    pagerank_downsample_adata,
    prepare_expression,
    read_inputs,
    run_icgs3,
    write_retention_audit,
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


def test_prepare_expression_rounds_normalized_sparse_values_to_four_decimals():
    X = sp.csr_matrix(np.array([[0.123456, 1.987654], [0.0, 2.222229]], dtype=float))
    adata = AnnData(X=X.copy(), obs=pd.DataFrame(index=["c1", "c2"]), var=pd.DataFrame(index=["g1", "g2"]))
    out = prepare_expression(
        adata,
        ICGS3Config(
            input_paths=["dummy.h5ad"],
            output_dir="/tmp/icgs3_test",
            input_normalized=True,
            normalized_decimals=4,
        ),
    )
    assert sp.issparse(out.X)
    np.testing.assert_allclose(out.X.toarray(), np.array([[0.1235, 1.9877], [0.0, 2.2222]], dtype=np.float32))
    np.testing.assert_allclose(out.layers["counts"].toarray(), X.toarray())
    assert out.uns["icgs3_normalized_decimals"] == 4


def test_prepare_expression_does_not_round_by_default():
    X = sp.csr_matrix(np.array([[0.123456, 1.987654]], dtype=float))
    adata = AnnData(X=X.copy(), obs=pd.DataFrame(index=["c1"]), var=pd.DataFrame(index=["g1", "g2"]))
    out = prepare_expression(
        adata,
        ICGS3Config(
            input_paths=["dummy.h5ad"],
            output_dir="/tmp/icgs3_test",
            input_normalized=True,
        ),
    )
    np.testing.assert_allclose(out.X.toarray(), X.toarray())
    assert "icgs3_normalized_decimals" not in out.uns


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


def test_automatic_pc_selection_records_elbow_metadata():
    rng = np.random.default_rng(11)
    X = rng.normal(0, 1, size=(50, 20))
    X[:25, :4] += 3
    adata = AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(50)]),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(20)]),
    )
    graph = compute_sparse_graph(
        adata,
        ICGS3Config(
            input_paths=["dummy.h5ad"],
            output_dir="/tmp/icgs3_test",
            n_pcs=0,
            max_auto_pcs=12,
            n_neighbors=8,
        ),
    )
    assert graph.uns["icgs3_graph"]["pc_selection"] == "elbow"
    assert 1 <= graph.uns["icgs3_graph"]["n_pcs"] <= 12
    assert graph.uns["icgs3_graph"]["pca_computed"] == 12


def test_downsampling_graph_can_use_sparse_expression_without_pca():
    rng = np.random.default_rng(12)
    X = rng.poisson(0.2, size=(40, 15)).astype(float)
    adata = AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(40)]),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(15)]),
    )
    graph = compute_sparse_graph(
        adata,
        ICGS3Config(input_paths=["dummy.h5ad"], output_dir="/tmp/icgs3_test", n_neighbors=6),
        use_pca=False,
    )
    assert graph.uns["icgs3_graph"]["representation"] == "X"
    assert graph.uns["icgs3_graph"]["pc_selection"] == "none"
    assert graph.uns["icgs3_graph"]["n_pcs"] == 0
    assert "X_pca" not in graph.obsm
    assert "connectivities" in graph.obsp


def test_pagerank_downsample_uses_icgs2_annoy_louvain_fields():
    rng = np.random.default_rng(13)
    X = rng.normal(0, 1, size=(45, 18)).astype(float)
    X[:15, :4] += 3
    X[15:30, 4:8] += 3
    X[30:, 8:12] += 3
    adata = AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(45)]),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(18)]),
    )
    sampled, scores = pagerank_downsample_adata(
        adata,
        ICGS3Config(
            input_paths=["dummy.h5ad"],
            output_dir="/tmp/icgs3_test",
            pagerank_cells=12,
            louvain_downsample_cutoff=20,
            downsample_var_genes=10,
        ),
    )
    assert sampled.n_obs <= 12
    assert "precommunity" in scores.columns
    assert "icgs2_community_size" in scores.columns
    assert "icgs2_chunk_start" in scores.columns
    assert "pagerank" in scores.columns


def test_rna_pagerank_downsample_filters_gene_space_before_dispersion_selection():
    rng = np.random.default_rng(31)
    genes = ["ACTB", "CD34", "RPL10", "RPS3", "MT-CO1", "DDX3Y", "HLA-A", "A.B"]
    X = rng.poisson(1.0, size=(30, len(genes))).astype(float)
    X[:, genes.index("RPL10")] += np.arange(30) * 10.0
    adata = AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(30)]),
        var=pd.DataFrame({"gene_symbols": genes}, index=genes),
    )
    sampled, _ = pagerank_downsample_adata(
        adata,
        ICGS3Config(
            input_paths=["dummy.h5ad"],
            output_dir="/tmp/icgs3_test",
            modality="rna",
            species="Hs",
            pagerank_cells=12,
            louvain_downsample_cutoff=1000,
            downsample_var_genes=2,
        ),
    )
    assert sampled.var_names.tolist() == ["ACTB", "CD34"]
    assert "RPL10" not in sampled.var_names
    assert "MT-CO1" not in sampled.var_names
    assert "DDX3Y" not in sampled.var_names


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


def test_marker_heatmap_annotations_render(tmp_path):
    from altanalyze3.components.visualization.marker_heatmap_h5ad import _plot_heatmap

    heatmap = pd.DataFrame(
        np.array(
            [
                [2.0, 1.5, -0.5, -1.0],
                [1.8, 1.2, -0.8, -1.2],
                [-1.0, -0.8, 1.4, 2.0],
                [-1.2, -0.4, 1.2, 1.8],
            ]
        ),
        index=["CD3", "CD4", "CD19", "MS4A1"],
        columns=["cell1", "cell2", "cell3", "cell4"],
    )
    covariates = pd.DataFrame(
        {
            "TissueGroup": ["BoneMarrow", "Thymus", "BoneMarrow", "Thymus"],
            "Library": ["L1", "L1", "L2", "L2"],
        },
        index=heatmap.columns,
    )
    out = tmp_path / "annotated_heatmap.pdf"
    _plot_heatmap(
        heatmap,
        str(out),
        [("C1", 2), ("C2", 2)],
        ["C1", "C2"],
        ["C1", "C1", "C2", "C2"],
        ["C1", "C1", "C2", "C2"],
        covariate_df=covariates,
        go_terms={"C1": [("T cell activation", 0.001)], "C2": [("B cell receptor signaling", 0.002)]},
        go_terms_max=2,
    )
    assert out.exists()
    assert out.with_suffix(".svg").exists()


def test_marker_heatmap_covariate_category_limit_and_palettes():
    from altanalyze3.components.visualization.marker_heatmap_h5ad import (
        _build_covariate_track,
        _filter_covariate_columns_for_heatmap,
    )

    cells = [f"cell{i}" for i in range(13)]
    obs = pd.DataFrame(
        {
            "small": ["A", "B"] * 6 + ["A"],
            "large": [f"L{i}" for i in range(13)],
        },
        index=cells,
    )
    retained, skipped = _filter_covariate_columns_for_heatmap(obs, cells, ["small", "large"])
    assert retained == ["small"]
    assert skipped == [("large", 13)]

    _, legends = _build_covariate_track(obs.loc[cells[:4], ["small", "large"]])
    assert legends["small"]["colors"]["A"] != legends["large"]["colors"]["L0"]


def test_auto_nmf_rank_uses_udon_rank_without_default_cap():
    from altanalyze3.components.clustering.ICGS import _resolve_auto_nmf_rank

    config = ICGS3Config(input_paths=["dummy.h5ad"], output_dir="/tmp/icgs3_test")
    candidate, rank, cap = _resolve_auto_nmf_rank(38, config)
    assert candidate == 38
    assert rank == 38
    assert cap is None

    config.max_auto_nmf_k = 30
    candidate, rank, cap = _resolve_auto_nmf_rank(38, config)
    assert candidate == 38
    assert rank == 30
    assert cap == 30


def test_svm_reclassify_all_cells_cli_default_and_false_alias():
    parser = build_arg_parser()
    args = parser.parse_args(["--input", "dummy.h5ad", "--output-dir", "/tmp/icgs3_test"])
    assert args.svm_reclassify_all_cells is True

    args = parser.parse_args(
        ["--input", "dummy.h5ad", "--output-dir", "/tmp/icgs3_test", "--no-svm-reclassify-all-cells"]
    )
    assert args.svm_reclassify_all_cells is False

    config = ICGS3Config(input_paths=["dummy.h5ad"], output_dir="/tmp/icgs3_test", svm_reclassify_all_cells=False)
    assert "--no-svm-reclassify-all-cells" in cli_equivalent(config)


def test_default_mito_percent_is_30():
    parser = build_arg_parser()
    args = parser.parse_args(["--input", "dummy.h5ad", "--output-dir", "/tmp/icgs3_test"])
    assert args.mito_percent == 30.0
    assert ICGS3Config(input_paths=["dummy.h5ad"], output_dir="/tmp/icgs3_test").mito_percent == 30.0
    assert args.pre_pagerank_cells == 0
    # The marker gate carries no parser-level default; apply_modality_defaults resolves it.
    # test_marker_gate_defaults_resolve_by_modality asserts the resolved values.
    assert args.marker_rho is None
    assert args.marker_min_per_cluster is None


def test_marker_gate_defaults_resolve_by_modality():
    """RNA keeps rho 0.3 / 2 markers; ADT gets the 0.25 / 1 gate; users still win.

    The superseded implementation detected an unset value by comparing against the literal
    defaults of the day, so a later default change silently disabled the ADT gate.
    """
    rna = apply_modality_defaults(
        ICGS3Config(input_paths=["dummy.h5ad"], output_dir="/tmp/icgs3_test", modality="rna")
    )
    assert rna.marker_rho == 0.3
    assert rna.marker_min_per_cluster == 2

    adt = apply_modality_defaults(
        ICGS3Config(input_paths=["dummy.h5ad"], output_dir="/tmp/icgs3_test", modality="adt")
    )
    assert adt.marker_rho == 0.25
    assert adt.marker_min_per_cluster == 1
    assert adt.min_genes == 0
    assert adt.min_counts == 0
    assert adt.mito_percent is None
    assert adt.n_top_features == 0

    explicit = apply_modality_defaults(
        ICGS3Config(
            input_paths=["dummy.h5ad"],
            output_dir="/tmp/icgs3_test",
            modality="adt",
            marker_rho=0.4,
            marker_min_per_cluster=5,
        )
    )
    assert explicit.marker_rho == 0.4
    assert explicit.marker_min_per_cluster == 5

    parser = build_arg_parser()
    args = parser.parse_args(
        ["--input", "dummy.h5ad", "--output-dir", "/tmp/icgs3_test", "--modality", "adt", "--marker-rho", "0.25"]
    )
    resolved = apply_modality_defaults(
        ICGS3Config(
            input_paths=list(args.input),
            output_dir=args.output_dir,
            modality=args.modality,
            marker_rho=args.marker_rho,
            marker_min_per_cluster=args.marker_min_per_cluster,
        )
    )
    assert resolved.marker_rho == 0.25
    assert resolved.marker_min_per_cluster == 1


def test_chunked_sparse_svm_matches_dense_classifier():
    genes = ["g1", "g2", "g3"]
    train = pd.DataFrame(
        {
            "P1": [3.0, 0.1, 0.0],
            "P2": [0.0, 0.2, 3.0],
        },
        index=genes,
    )
    X = np.array(
        [
            [4.0, 0.0, 0.0],
            [3.5, 0.1, 0.0],
            [0.0, 0.0, 4.0],
            [0.1, 0.0, 3.5],
        ]
    )
    cells = [f"cell{i}" for i in range(X.shape[0])]
    adata = AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(index=genes),
    )
    dense = classify_with_scores(
        train=train,
        expression=pd.DataFrame(X.T, index=genes, columns=cells),
        cluster_key="cluster",
        min_decision_score=0.0,
    )
    chunked = classify_adata_with_scores_chunked(
        train=train,
        adata=adata,
        genes=genes,
        cluster_key="cluster",
        min_decision_score=0.0,
        chunk_size=2,
    )
    pd.testing.assert_series_equal(dense["cluster"], chunked["cluster"])
    np.testing.assert_allclose(dense["svm_score"].values, chunked["svm_score"].values)


def test_retention_audit_flags_rare_population_loss(tmp_path):
    obs = pd.DataFrame(
        {"HLCA": ["rare"] * 30 + ["common"] * 70},
        index=[f"cell{i}" for i in range(100)],
    )
    var = pd.DataFrame(index=["gene1", "gene2"])
    full = AnnData(X=sp.csr_matrix(np.ones((100, 2))), obs=obs, var=var)
    candidates = full.copy()
    sampled = full[[f"cell{i}" for i in range(10)] + [f"cell{i}" for i in range(30, 70)]].copy()
    outdir = tmp_path / "icgs3"
    (outdir / "sNMF").mkdir(parents=True)
    audit = write_retention_audit(
        full=full,
        candidates=candidates,
        sampled=sampled,
        final_clusters=None,
        config=ICGS3Config(input_paths=["dummy.h5ad"], output_dir=str(outdir), retention_audit_obs="HLCA"),
        outdir=str(outdir),
    )
    rare = audit[audit["annotation"] == "rare"].iloc[0]
    assert rare["post_qc_cells"] == 30
    assert rare["pagerank_louvain_sampled_cells"] == 10
    assert bool(rare["sampled_below_20"])
    assert (outdir / "sNMF" / "icgs3_rare_population_retention_audit.tsv").exists()


# ---------------------------------------------------------------------------
# ADT occupancy filter (AltAnalyze2 removeCellsWithHighADTOccupancy)
# ---------------------------------------------------------------------------


def _legacy_occupancy(values: np.ndarray, sd_multiplier: float = 1.0) -> np.ndarray:
    """Direct transcription of AltAnalyze2 clustering.py:10973-10986.

    values is features x cells, exactly the orientation the legacy function read from the
    ADT value file. statistics.stdev there is the SAMPLE standard deviation, so ddof=1.
    """
    marks = []
    for adt_vector in values:
        stdev = float(np.std(adt_vector, ddof=1))
        avg = float(np.mean(adt_vector))
        threshold = avg + sd_multiplier * stdev
        marks.append([1 if v > threshold else 0 for v in adt_vector])
    return np.asarray(marks).transpose().sum(axis=1)


def _adt_adata(values_cells_by_features: np.ndarray, sparse: bool = True):
    n_obs, n_vars = values_cells_by_features.shape
    X = sp.csr_matrix(values_cells_by_features) if sparse else values_cells_by_features.astype(np.float32)
    return AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=[f"ADT{j}" for j in range(n_vars)]),
    )


def _occupancy_config(**kwargs):
    base = dict(input_paths=["dummy.h5ad"], output_dir="/tmp/icgs3_test", modality="adt")
    base.update(kwargs)
    return ICGS3Config(**base)


def test_adt_occupancy_matches_legacy_binarization():
    """compute_adt_occupancy must equal the legacy mean+1SD mark count, cell for cell."""
    rng = np.random.default_rng(0)
    values = rng.gamma(shape=2.0, scale=5.0, size=(400, 24)).astype(np.float64)
    adata = _adt_adata(values, sparse=True)
    out = compute_adt_occupancy(adata, _occupancy_config(), stage="test")
    expected = _legacy_occupancy(values.T, sd_multiplier=1.0)
    np.testing.assert_array_equal(out.obs[ADT_OCCUPANCY_KEY].to_numpy(), expected)
    assert int(out.uns["icgs3_adt_occupancy"]["cells"]) == 400
    assert int(out.uns["icgs3_adt_occupancy"]["features"]) == 24


def test_adt_occupancy_sparse_and_dense_paths_agree():
    """The sparse-nonzero path and the dense-block path must return identical scores."""
    rng = np.random.default_rng(1)
    values = rng.gamma(shape=1.5, scale=3.0, size=(250, 16))
    values[values < 1.0] = 0.0  # force real sparsity; the two paths must still agree
    sparse_out = compute_adt_occupancy(_adt_adata(values, sparse=True), _occupancy_config(), stage="sparse")
    dense_out = compute_adt_occupancy(_adt_adata(values, sparse=False), _occupancy_config(), stage="dense")
    np.testing.assert_array_equal(
        sparse_out.obs[ADT_OCCUPANCY_KEY].to_numpy(),
        dense_out.obs[ADT_OCCUPANCY_KEY].to_numpy(),
    )
    np.testing.assert_array_equal(
        sparse_out.obs[ADT_OCCUPANCY_KEY].to_numpy(),
        _legacy_occupancy(values.T),
    )


def test_occupancy_trim_removes_both_tails_inside_each_group():
    """20 percent per tail, computed inside the group, not across the whole matrix."""
    # Group A occupancies 0..49, group B occupancies 100..149. A global percentile would
    # delete all of A's low cells and all of B's high cells; a per-group trim must not.
    barcodes = [f"c{i}" for i in range(100)]
    occupancy = pd.Series(list(range(50)) + list(range(100, 150)), index=barcodes, dtype=float)
    groups = pd.Series(["A"] * 50 + ["B"] * 50, index=barcodes)
    keep_mask, removed, summary = _trim_occupancy_within_groups(
        occupancy, groups, percent=20.0, min_group=10
    )
    assert int(keep_mask.sum()) == 60, "50 cells per group, 10 removed per tail, 30 kept each"
    assert removed.shape[0] == 40
    for group in ("A", "B"):
        sub = removed[removed["broad_group"] == group]
        assert int((sub[ADT_OCCUPANCY_STATUS_KEY] == "removed_understained").sum()) == 10
        assert int((sub[ADT_OCCUPANCY_STATUS_KEY] == "removed_overstained").sum()) == 10
    # The retained cells of each group are its middle 60 percent, by that group's own range.
    kept = set(np.asarray(barcodes)[keep_mask])
    assert kept.issuperset({f"c{i}" for i in range(10, 40)})
    assert kept.issuperset({f"c{i}" for i in range(60, 90)})
    assert kept.isdisjoint({f"c{i}" for i in range(10)})
    assert kept.isdisjoint({f"c{i}" for i in range(90, 100)})
    assert summary.shape[0] == 2
    assert list(summary["n_kept"]) == [30, 30]


def test_occupancy_trim_keeps_small_group_that_legacy_emptied():
    """AltAnalyze2 scores[int(l*p):-1*int(l*p)] is scores[0:0] when int(l*p) is 0.

    A 4-cell group at 20 percent gives int(4*0.2)=0, so the legacy slice returned the EMPTY
    list and deleted the whole group. ICGS3 must keep it and label it.
    """
    barcodes = [f"s{i}" for i in range(4)]
    occupancy = pd.Series([1.0, 2.0, 3.0, 40.0], index=barcodes)
    groups = pd.Series(["small"] * 4, index=barcodes)
    legacy_slice = list(range(4))[int(4 * 0.2): -1 * int(4 * 0.2)]
    assert legacy_slice == [], "the legacy expression really does empty this group"
    keep_mask, removed, summary = _trim_occupancy_within_groups(
        occupancy, groups, percent=20.0, min_group=2
    )
    assert int(keep_mask.sum()) == 4
    assert removed.empty
    assert bool(summary["kept_whole_group"].iloc[0]) is True


def test_occupancy_trim_is_deterministic_under_ties():
    """Equal occupancies must resolve the same way on every run."""
    barcodes = [f"t{i}" for i in range(20)]
    occupancy = pd.Series([5.0] * 20, index=barcodes)
    groups = pd.Series(["g"] * 20, index=barcodes)
    first = _trim_occupancy_within_groups(occupancy, groups, percent=20.0, min_group=5)[0]
    second = _trim_occupancy_within_groups(occupancy, groups, percent=20.0, min_group=5)[0]
    np.testing.assert_array_equal(first, second)
    assert int(first.sum()) == 12


def test_occupancy_percent_out_of_range_is_rejected():
    occupancy = pd.Series([1.0, 2.0], index=["a", "b"])
    groups = pd.Series(["g", "g"], index=["a", "b"])
    with pytest.raises(ValueError, match="adt-occupancy-percent"):
        _trim_occupancy_within_groups(occupancy, groups, percent=50.0, min_group=1)


def test_normalization_resolves_by_modality_and_honours_legacy_flag():
    adt = apply_modality_defaults(_occupancy_config())
    rna = ICGS3Config(input_paths=["d.h5ad"], output_dir="/tmp/icgs3_test", modality="rna")
    assert resolve_normalization(adt) == "log1p"
    assert resolve_normalization(rna) == "cp10k-log1p"
    assert resolve_normalization(_occupancy_config(input_normalized=True)) == "none"
    assert resolve_normalization(_occupancy_config(normalization="cp10k-log1p")) == "cp10k-log1p"
    with pytest.raises(ValueError, match="conflicts"):
        resolve_normalization(_occupancy_config(input_normalized=True, normalization="log1p"))


def test_log1p_normalization_skips_library_size_rescaling():
    values = np.array([[1.0, 9.0], [4.0, 4.0]], dtype=np.float64)
    out = prepare_expression(_adt_adata(values), _occupancy_config(normalization="log1p"))
    np.testing.assert_allclose(out.X.toarray(), np.log1p(values).astype(np.float32), rtol=1e-6)
    assert out.uns["icgs3_normalization_mode"] == "log1p"
    # CP10K would have forced both rows to the same total; log1p alone must not.
    cp10k = prepare_expression(_adt_adata(values), _occupancy_config(normalization="cp10k-log1p"))
    assert not np.allclose(cp10k.X.toarray(), out.X.toarray())


def test_log1p_normalization_rejects_negative_values():
    values = np.array([[-1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
    with pytest.raises(ValueError, match="log1p"):
        prepare_expression(_adt_adata(values, sparse=False), _occupancy_config(normalization="log1p"))


def test_small_feature_rank_branch_resolves_by_modality():
    """rank_rel_threshold existed as a dead CLI parameter; it must now reach the estimator."""
    assert resolve_small_feature_rank(_occupancy_config()) is True
    assert resolve_small_feature_rank(
        ICGS3Config(input_paths=["d.h5ad"], output_dir="/tmp/icgs3_test", modality="rna")
    ) is False
    assert resolve_small_feature_rank(_occupancy_config(small_feature_rank="false")) is False
    assert resolve_small_feature_rank(
        ICGS3Config(input_paths=["d.h5ad"], output_dir="/tmp/icgs3_test",
                    modality="rna", small_feature_rank="true")
    ) is True


def test_estimate_udon_nmf_rank_forwards_small_feature_and_rel_threshold(monkeypatch):
    """rank_rel_threshold was declared, echoed by cli_equivalent, and read by nothing.

    determine_nmf_ranks uses it ONLY inside its small_feature branch, and ICGS.py:2073 passed
    small_feature=False for every modality, so no run ever reached that branch and the
    parameter was dead. Capturing the call proves both values now arrive, and that ADT turns
    the branch on while RNA leaves it off.
    """
    from altanalyze3.components.clustering.ICGS import _add_udon_path

    _add_udon_path()
    import nmf as udon_nmf

    captured = {}

    def fake_determine_nmf_ranks(df, small_feature=False, rel_threshold=0.1, max_rank=30):
        captured["small_feature"] = small_feature
        captured["rel_threshold"] = rel_threshold
        return 7

    monkeypatch.setattr(udon_nmf, "determine_nmf_ranks", fake_determine_nmf_ranks)
    expr = pd.DataFrame(np.ones((3, 4)), index=list("abc"), columns=list("wxyz"))

    rank = estimate_udon_nmf_rank(expr, _occupancy_config(rank_rel_threshold=0.42))
    assert rank == 7
    assert captured["small_feature"] is True, "ADT must reach the small_feature branch"
    assert captured["rel_threshold"] == 0.42, "--rank-rel-threshold must reach the estimator"

    estimate_udon_nmf_rank(
        expr,
        ICGS3Config(input_paths=["d.h5ad"], output_dir="/tmp/icgs3_test", modality="rna",
                    rank_rel_threshold=0.42),
    )
    assert captured["small_feature"] is False, "RNA must stay on the Tracy-Widom boundary"


def test_adt_occupancy_filter_end_to_end_removes_cells_from_every_output(tmp_path):
    """A removed cell must reach no NMF, MarkerFinder or SVM step and appear in no output.

    The matrix carries three real ADT blocks plus two planted technical populations: cells
    stained high for the whole panel and cells stained low for the whole panel. Both must be
    gone from the cluster table and from the returned AnnData, not merely flagged.
    """
    rng = np.random.default_rng(11)
    n_per_block, n_features = 120, 20
    blocks = []
    for b in range(3):
        values = rng.gamma(2.0, 2.0, size=(n_per_block, n_features))
        values[:, b * 5:(b + 1) * 5] += rng.gamma(6.0, 4.0, size=(n_per_block, 5))
        blocks.append(values)
    values = np.vstack(blocks)
    names = [f"cell{i}" for i in range(values.shape[0])]
    # 30 over-stained and 30 under-stained cells, spread across the three blocks.
    over = np.arange(0, 90, 3)[:30]
    under = np.arange(1, 90, 3)[:30]
    values[over] *= 12.0
    values[under] *= 0.02
    adata = AnnData(
        X=sp.csr_matrix(values.astype(np.float32)),
        obs=pd.DataFrame(index=names),
        var=pd.DataFrame(index=[f"ADT{j}" for j in range(n_features)]),
    )
    path = tmp_path / "adt.h5ad"
    adata.write_h5ad(path)

    result = run_icgs3(
        ICGS3Config(
            input_paths=[str(path)],
            output_dir=str(tmp_path / "out"),
            modality="adt",
            normalization="log1p",
            adt_occupancy_filter=True,
            adt_occupancy_percent=20.0,
            adt_occupancy_scale="input",
            adt_occupancy_scope="all",
            adt_occupancy_min_group=5,
            generate_umap=False,
            write_h5ad=False,
        )
    )

    removed_path = tmp_path / "out" / "sNMF" / "icgs3_adt_occupancy_removed_cells.tsv"
    summary_path = tmp_path / "out" / "sNMF" / "icgs3_adt_occupancy_group_summary.tsv"
    assert removed_path.exists() and summary_path.exists()
    removed = pd.read_csv(removed_path, sep="\t")
    assert removed.shape[0] > 0, "the filter removed nothing on a matrix with planted outliers"

    removed_barcodes = set(removed["barcode"].astype(str))
    surviving = set(result.adata.obs_names.astype(str))
    assert removed_barcodes.isdisjoint(surviving), "a removed cell survived into the AnnData"
    assert removed_barcodes.isdisjoint(set(result.clusters.index.astype(str))), \
        "a removed cell received a cluster"
    assert len(surviving) + len(removed_barcodes) == values.shape[0]

    # Both technical tails must be caught, not just one.
    over_names = {names[i] for i in over}
    under_names = {names[i] for i in under}
    assert len(removed_barcodes & over_names) > 0, "no over-stained cell was removed"
    assert len(removed_barcodes & under_names) > 0, "no under-stained cell was removed"
    assert set(removed[ADT_OCCUPANCY_STATUS_KEY]) <= {"removed_overstained", "removed_understained"}
