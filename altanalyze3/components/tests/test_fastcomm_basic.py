from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from altanalyze3.components.fastComm.api import FastCommParams, run_fastcomm
from altanalyze3.components.fastComm.benchmark import FastCommBenchmarkParams, run_benchmark
from altanalyze3.components.fastComm.reporting import ExemplarReportParams, build_exemplar_table, write_exemplar_report
from altanalyze3.components.fastComm.training import ResponseTrainingParams, build_response_matrix
from altanalyze3.components.fastComm.upstream_resources import (
    build_cellchat_lr_table,
    build_nichenet_lr_table,
    build_nichenet_response_matrix,
    merge_ligand_receptor_tables,
)


def test_fastcomm_prioritizes_response_supported_edge(tmp_path: Path):
    cells = [f"cell{i:03d}" for i in range(12)]
    genes = ["TGFB1", "TGFBR1", "TGFBR2", "SERPINE1", "SMAD7", "TNF", "TNFRSF1A"]
    expression = pd.DataFrame(0.0, index=cells, columns=genes)
    metadata = pd.DataFrame(
        {"cell_state": ["sender"] * 4 + ["receiver_base"] * 4 + ["receiver_activated"] * 4},
        index=cells,
    )

    expression.loc[cells[:4], "TGFB1"] = 8.0
    expression.loc[cells[4:], ["TGFBR1", "TGFBR2"]] = 5.0
    expression.loc[cells[8:], ["SERPINE1", "SMAD7"]] = 6.0
    expression.loc[cells[:4], "TNF"] = 6.0
    expression.loc[cells[4:], "TNFRSF1A"] = 4.0

    expression_path = tmp_path / "expression.tsv"
    metadata_path = tmp_path / "metadata.tsv"
    response_path = tmp_path / "response.tsv"
    output_path = tmp_path / "fastcomm.tsv"
    expression.to_csv(expression_path, sep="\t")
    metadata.to_csv(metadata_path, sep="\t")

    response = pd.DataFrame(
        {
            "SERPINE1": [1.0, 0.0],
            "SMAD7": [1.0, 0.0],
            "TNF": [0.0, 1.0],
        },
        index=["TGFB1", "TNF"],
    )
    response.to_csv(response_path, sep="\t")

    result = run_fastcomm(
        FastCommParams(
            expression=expression_path,
            metadata=metadata_path,
            response_matrix=response_path,
            output=output_path,
            state_key="cell_state",
            baseline_state="receiver_base",
            include_self_edges=False,
        )
    )

    assert output_path.exists()
    assert not result.scores.empty
    top = result.scores.iloc[0]
    assert top["sender_state"] == "sender"
    assert top["receiver_state"] == "receiver_activated"
    assert top["ligand"] == "TGFB1"
    assert top["receiver_response_score"] > 0.9
    assert "SERPINE1" in top["response_support_genes"]
    assert 0.0 <= top["fastcomm_percentile"] <= 1.0
    assert np.isfinite(result.scores["fastcomm_score"]).all()


def test_fastcomm_h5ad_loads_only_required_genes(tmp_path: Path):
    cells = [f"cell{i:03d}" for i in range(8)]
    genes = ["TGFB1", "TGFBR1", "TGFBR2", "SERPINE1", "SMAD7"] + [f"NOISE{i:03d}" for i in range(50)]
    X = np.zeros((len(cells), len(genes)), dtype=np.float32)
    gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}
    X[:4, gene_to_idx["TGFB1"]] = 5.0
    X[4:, gene_to_idx["TGFBR1"]] = 4.0
    X[4:, gene_to_idx["TGFBR2"]] = 4.0
    X[4:, gene_to_idx["SERPINE1"]] = 3.0
    X[4:, gene_to_idx["SMAD7"]] = 3.0

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame({"cell_state": ["sender"] * 4 + ["receiver"] * 4}, index=cells),
        var=pd.DataFrame(index=genes),
    )
    h5ad_path = tmp_path / "mini.h5ad"
    adata.write_h5ad(h5ad_path)

    lr_path = tmp_path / "lr.tsv"
    pd.DataFrame(
        {
            "ligand": ["TGFB1"],
            "receptor": ["TGFBR1+TGFBR2"],
            "pathway": ["TGFB1"],
            "evidence_weight": [1.0],
        }
    ).to_csv(lr_path, sep="\t", index=False)
    response_path = tmp_path / "response.tsv"
    pd.DataFrame({"SERPINE1": [1.0], "SMAD7": [1.0]}, index=["TGFB1"]).to_csv(response_path, sep="\t")

    result = run_fastcomm(
        FastCommParams(
            h5ad=h5ad_path,
            lr_table=lr_path,
            response_matrix=response_path,
            output=tmp_path / "scores.tsv",
            include_self_edges=False,
        )
    )

    assert result.summary["n_genes"] == 5
    assert result.summary["n_loaded_genes"] == 5
    assert result.scores.iloc[0]["ligand"] == "TGFB1"


def test_fastcomm_scores_mouse_title_case_symbols(tmp_path: Path):
    cells = [f"cell{i:03d}" for i in range(10)]
    genes = ["Tgfb1", "Tgfbr1", "Tgfbr2", "Serpine1", "Smad7", "NoiseGene"]
    X = np.zeros((len(cells), len(genes)), dtype=np.float32)
    gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}
    X[:5, gene_to_idx["Tgfb1"]] = 6.0
    X[5:, gene_to_idx["Tgfbr1"]] = 4.0
    X[5:, gene_to_idx["Tgfbr2"]] = 4.0
    X[5:, gene_to_idx["Serpine1"]] = 5.0
    X[5:, gene_to_idx["Smad7"]] = 4.0

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame({"cell_state": ["sender"] * 5 + ["receiver"] * 5}, index=cells),
        var=pd.DataFrame(index=genes),
    )
    h5ad_path = tmp_path / "mouse.h5ad"
    adata.write_h5ad(h5ad_path)

    lr_path = tmp_path / "lr.tsv"
    pd.DataFrame(
        {
            "ligand": ["TGFB1"],
            "receptor": ["TGFBR1+TGFBR2"],
            "pathway": ["TGFB1"],
            "evidence_weight": [1.0],
        }
    ).to_csv(lr_path, sep="\t", index=False)
    response_path = tmp_path / "response.tsv"
    pd.DataFrame({"SERPINE1": [1.0], "SMAD7": [1.0]}, index=["TGFB1"]).to_csv(response_path, sep="\t")

    result = run_fastcomm(
        FastCommParams(
            h5ad=h5ad_path,
            lr_table=lr_path,
            response_matrix=response_path,
            output=tmp_path / "scores.tsv",
            include_self_edges=False,
        )
    )

    assert result.summary["n_matched_required_genes"] == 5
    assert result.summary["n_missing_required_genes"] == 0
    assert set(result.state_expression.columns) == {"TGFB1", "TGFBR1", "TGFBR2", "SERPINE1", "SMAD7"}
    top = result.scores.iloc[0]
    assert top["ligand"] == "TGFB1"
    assert top["receptor"] == "TGFBR1+TGFBR2"
    assert top["receiver_response_score"] > 0.9


def test_fastcomm_benchmark_writes_split_stability(tmp_path: Path):
    cells = [f"cell{i:03d}" for i in range(16)]
    genes = ["TGFB1", "TGFBR1", "TGFBR2", "SERPINE1", "SMAD7"]
    X = np.zeros((len(cells), len(genes)), dtype=np.float32)
    gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}
    X[:8, gene_to_idx["TGFB1"]] = 5.0
    X[8:, gene_to_idx["TGFBR1"]] = 4.0
    X[8:, gene_to_idx["TGFBR2"]] = 4.0
    X[8:, gene_to_idx["SERPINE1"]] = 3.0
    X[8:, gene_to_idx["SMAD7"]] = 3.0

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(
            {
                "cell_state": ["sender"] * 8 + ["receiver"] * 8,
                "Donor": ["D1"] * 4 + ["D2"] * 4 + ["D1"] * 4 + ["D2"] * 4,
            },
            index=cells,
        ),
        var=pd.DataFrame(index=genes),
    )
    h5ad_path = tmp_path / "benchmark.h5ad"
    adata.write_h5ad(h5ad_path)

    lr_path = tmp_path / "lr.tsv"
    pd.DataFrame(
        {
            "ligand": ["TGFB1"],
            "receptor": ["TGFBR1+TGFBR2"],
            "pathway": ["TGFB1"],
            "evidence_weight": [1.0],
        }
    ).to_csv(lr_path, sep="\t", index=False)
    response_path = tmp_path / "response.tsv"
    pd.DataFrame({"SERPINE1": [1.0], "SMAD7": [1.0]}, index=["TGFB1"]).to_csv(response_path, sep="\t")

    output_dir = tmp_path / "bench"
    summary = run_benchmark(
        FastCommBenchmarkParams(
            h5ad=h5ad_path,
            output_dir=output_dir,
            lr_table=lr_path,
            response_matrix=response_path,
            min_cells=2,
            top_n_stability=10,
        )
    )

    assert summary["n_splits"] == 2
    assert (output_dir / "full_scores.tsv").exists()
    stability = pd.read_csv(output_dir / "split_stability.tsv", sep="\t")
    assert set(stability["status"]) == {"ok"}


def test_fastcomm_builds_response_matrix_from_long_table(tmp_path: Path):
    long = pd.DataFrame(
        {
            "signature": ["TGFB1", "TGFB1", "TGFB1", "TNF"],
            "gene": ["SERPINE1", "SMAD7", "NOISE", "NFKBIA"],
            "score": [3.0, 2.0, 0.1, -4.0],
        }
    )
    input_path = tmp_path / "long.tsv"
    output_path = tmp_path / "response.tsv"
    manifest_path = tmp_path / "manifest.json"
    long.to_csv(input_path, sep="\t", index=False)

    manifest = build_response_matrix(
        ResponseTrainingParams(
            input=input_path,
            output=output_path,
            manifest=manifest_path,
            top_genes=2,
            min_abs_score=0.5,
        )
    )

    matrix = pd.read_csv(output_path, sep="\t", index_col=0)
    assert manifest["output_shape"] == [2, 3]
    assert "NOISE" not in matrix.columns
    assert np.isclose(np.sqrt((matrix.loc["TGFB1"] ** 2).sum()), 1.0)
    assert manifest_path.exists()


def test_fastcomm_exemplar_report(tmp_path: Path):
    scores = pd.DataFrame(
        {
            "sender_state": ["sender", "sender"],
            "receiver_state": ["receiver", "receiver"],
            "ligand": ["TGFB1", "TNF"],
            "receptor": ["TGFBR1+TGFBR2", "TNFRSF1A"],
            "pathway": ["TGF_beta", "TNF"],
            "interaction_class": ["secreted", "secreted"],
            "ligand_expr": [4.0, 3.0],
            "receptor_expr": [5.0, 2.0],
            "ligand_detection": [1.0, 0.8],
            "receptor_detection": [0.9, 0.7],
            "complex_completeness": [0.9, 0.7],
            "lr_expression_score_scaled": [0.8, 0.4],
            "receiver_response_score": [0.7, 0.2],
            "state_promotion_score": [0.7, 0.2],
            "response_key": ["TGFB1", "TNF"],
            "response_support_genes": ["SERPINE1:1;SMAD7:0.5", "NFKBIA:0.2"],
            "fastcomm_score": [0.75, 0.31],
            "fastcomm_percentile": [0.99, 0.5],
            "lr_within_lr_percentile": [0.95, 0.4],
            "score_rank": [1, 2],
        }
    )
    scores_path = tmp_path / "scores.tsv"
    output_tsv = tmp_path / "exemplars.tsv"
    output_md = tmp_path / "exemplars.md"
    scores.to_csv(scores_path, sep="\t", index=False)

    params = ExemplarReportParams(
        scores=scores_path,
        output_tsv=output_tsv,
        output_md=output_md,
        top_n=5,
    )
    table = build_exemplar_table(params)
    manifest = write_exemplar_report(params)

    assert table.iloc[0]["confidence"] == "high"
    assert table.iloc[0]["prototype_probability"] == 0.75
    assert manifest["n_exemplars"] == 1
    assert output_tsv.exists()
    assert "TGFB1" in output_md.read_text(encoding="utf-8")


def test_fastcomm_builds_lr_table_from_cellchat_and_nichenet():
    cellchat_db = {
        "interaction": pd.DataFrame(
            {
                "interaction_name": ["TGFB1_TGFBR1_TGFBR2", "CXCL12_CXCR4"],
                "pathway_name": ["TGFb", "CXCL"],
                "ligand": ["TGFB1", "CXCL12"],
                "receptor": ["TGFbR1_R2", "CXCR4"],
                "annotation": ["Secreted Signaling", "Secreted Signaling"],
                "evidence": ["KEGG: hsa04350", "PMID:123"],
                "ligand.symbol": ["TGFB1", "CXCL12"],
                "receptor.symbol": ["TGFBR1, TGFBR2", "CXCR4"],
                "agonist": ["", ""],
                "antagonist": ["", ""],
                "co_A_receptor": ["", ""],
                "co_I_receptor": ["", ""],
                "version": ["CellChatDB v2", "CellChatDB v2"],
            }
        ),
        "complex": pd.DataFrame(
            {"subunit_1": ["TGFBR1"], "subunit_2": ["TGFBR2"]},
            index=["TGFbR1_R2"],
        ),
    }
    nichenet_lr = pd.DataFrame(
        {
            "from": ["CXCL12", "FLT3LG"],
            "to": ["CXCR4", "FLT3"],
            "database": ["omnipath", "omnipath"],
            "source": ["omnipath", "omnipath"],
        }
    )

    cellchat_lr = build_cellchat_lr_table(cellchat_db)
    nichenet_only = build_nichenet_lr_table(nichenet_lr)
    merged = merge_ligand_receptor_tables(cellchat_lr, nichenet_only)
    observed_pairs = set(merged[["ligand", "receptor"]].apply(tuple, axis=1))

    assert ("TGFB1", "TGFBR1+TGFBR2") in observed_pairs
    assert ("CXCL12", "CXCR4") in observed_pairs
    assert ("FLT3LG", "FLT3") in observed_pairs
    assert merged.loc[merged["ligand"].eq("FLT3LG"), "source"].iloc[0] == "NicheNet"


def test_fastcomm_builds_nichenet_response_matrix_with_pathway_rows():
    lr_table = pd.DataFrame(
        {
            "ligand": ["TGFB1", "TGFB2", "CXCL12"],
            "receptor": ["TGFBR1+TGFBR2", "TGFBR1+TGFBR2", "CXCR4"],
            "pathway": ["TGFb", "TGFb", "CXCL"],
            "evidence_weight": [1.0, 1.0, 1.0],
            "interaction_class": ["secreted", "secreted", "secreted"],
        }
    )
    ligand_target = pd.DataFrame(
        {
            "TGFB1": [0.9, 0.8, 0.05],
            "TGFB2": [0.7, 0.6, 0.04],
            "CXCL12": [0.2, 0.1, 0.95],
        },
        index=["SERPINE1", "SMAD7", "CXCL10"],
    )

    response = build_nichenet_response_matrix(
        lr_table,
        ligand_target,
        top_targets_per_signature=2,
        min_target_weight=0.1,
    )

    assert "TGFB1" in response.index
    assert "CXCL12" in response.index
    assert "TGFb" in response.index
    assert "CXCL" in response.index
    assert "CXCL10" in response.columns
    assert np.isclose(np.sqrt((response.loc["TGFB1"] ** 2).sum()), 1.0)
