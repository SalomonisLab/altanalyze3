"""Tests for the post-pruning annotate default (steps 4-6)."""

import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from altanalyze3.components.sctriangulate import annotate as A


def _toy_adata(n_per=40, n_genes=60, n_adt=6, seed=0):
    """Three well-separated clusters, plus antibody features so the ADT branch runs."""
    rng = np.random.default_rng(seed)
    blocks, labels = [], []
    for i in range(3):
        x = rng.poisson(0.4, size=(n_per, n_genes)).astype(np.float32)
        x[:, i * 12:(i + 1) * 12] += rng.poisson(9.0, size=(n_per, 12))
        blocks.append(x)
        labels += [f"src{i}@C{i}"] * n_per
    rna = np.vstack(blocks)
    adt = np.zeros((rna.shape[0], n_adt), dtype=np.float32)
    for i in range(3):
        adt[i * n_per:(i + 1) * n_per, i * 2:(i + 1) * 2] += rng.poisson(7.0, size=(n_per, 2))
    X = np.hstack([rna, adt])
    adata = ad.AnnData(X=sp.csr_matrix(np.log1p(X)))
    adata.var_names = ([f"Gene{i}" for i in range(n_genes)]
                       + [f"{A.DEFAULT_ADT_PREFIX}Ab{i}_TotalA" for i in range(n_adt)])
    adata.obs_names = [f"cell{i}" for i in range(rna.shape[0])]
    adata.obs["pruned"] = pd.Categorical(labels)
    adata.obs["curated"] = pd.Categorical(
        ["Alpha"] * n_per + ["unassigned"] * n_per + ["Gamma"] * n_per)
    return adata


def test_species_maps_to_icgs_codes():
    assert A._icgs_species("mouse") == "Mm"
    assert A._icgs_species("Mm") == "Mm"
    assert A._icgs_species("human") == "Hs"
    with pytest.raises(ValueError):
        A._icgs_species("zebrafish")


def test_uninformative_labels_cover_the_no_call_vocabulary():
    for label in ("unassigned", "Unknown", "NA", "none", ""):
        assert label.strip().lower() in A.UNINFORMATIVE_LEAD_LABELS
    # QC calls DO name a cluster and must not be treated as missing
    for label in ("doublet", "dying", "erythrocyte"):
        assert label not in A.UNINFORMATIVE_LEAD_LABELS


def test_fallback_label_strips_the_antibody_decoration():
    markers = pd.DataFrame({
        "Gene": [f"{A.DEFAULT_ADT_PREFIX}CD4_TotalA", "Foxp3"],
        "Fold": [9.0, 1.0],
        "cluster": ["c1", "c1"],
    })
    assert A._fallback_label("c1", markers, A.DEFAULT_ADT_PREFIX) == "CD4-high"
    assert A._fallback_label("absent", markers, A.DEFAULT_ADT_PREFIX) == "Unresolved"


def test_enrichment_gates_reject_a_weak_term(monkeypatch, tmp_path):
    """A term at FDR 0.03 with an overlap of 2 must not name a cluster."""
    import altanalyze3.components.clustering.ICGS as ICGS

    labels = pd.DataFrame({
        "cluster": ["c1", "c2"],
        "term_name": ["Adult Marrow granulocyte", "Adult Retina rods"],
        "fdr": [1e-40, 3e-2],
        "overlap": [30, 2],
    })
    monkeypatch.setattr(ICGS, "_default_biomarker_file", lambda species: __file__)
    monkeypatch.setattr(ICGS, "biomarker_enrichment", lambda *a, **k: labels)

    markers = pd.DataFrame({"Gene": ["A", "B"], "cluster": ["c1", "c2"],
                            "Fold": [3.0, 3.0], "_is_adt": [False, False]})
    kept = A._predictions_from_enrichment(
        markers, ["A", "B"], ["c1", "c2"], "mouse", str(tmp_path), None,
        1e-5, 5, log=lambda *_: None,
    )
    assert "c1" in kept and "c2" not in kept

    loose = A._predictions_from_enrichment(
        markers, ["A", "B"], ["c1", "c2"], "mouse", str(tmp_path), None,
        0.05, 1, log=lambda *_: None,
    )
    assert "c2" in loose, "a loosened gate must still let the weak term through"


def test_annotate_end_to_end_without_a_biomarker_file(tmp_path):
    """No BioMarkers table: every name falls back to the marker, and the invariants still hold."""
    adata = _toy_adata()
    out = tmp_path / "run"
    summary = A.annotate_pruned_clusters(
        adata, str(out), cluster_key="pruned", species="mouse",
        top_n=6, cells_per_cluster=10, hopach_mincluster=2,
        biomarker_file=str(tmp_path / "does_not_exist.txt"),
        log=lambda *_: None,
    )

    assert set(summary["name_source"].values()) == {"top-marker"}
    assert len(summary["names"]) == 3
    assert adata.obs["cluster_name"].notna().all()
    assert adata.obs["cluster_name"].nunique() == 3
    assert adata.obs["hopach_cluster"].notna().all()
    assert list(adata.uns["lineage_order"]) == list(
        pd.read_csv(summary["cluster_table"], sep="\t")["named_label"])
    # every name carries its source cluster, so a reader can trace it back
    for cluster, name in summary["names"].items():
        assert name.endswith(f"({cluster})")
    # sizes survive the renaming
    assert (sorted(adata.obs["cluster_name"].value_counts().tolist())
            == sorted(adata.obs["pruned"].value_counts().tolist()))
    for path in (summary["cluster_table"], summary["ordered_centroids"]):
        assert os.path.exists(path)
    assert summary["marker_outputs_adt"] is not None, "ADT companion must run when AB_ features exist"


def test_annotate_lead_skips_an_uninformative_label(tmp_path):
    adata = _toy_adata()
    summary = A.annotate_pruned_clusters(
        adata, str(tmp_path / "lead"), cluster_key="pruned", species="mouse",
        top_n=6, cells_per_cluster=10, hopach_mincluster=2,
        lead_annotation="curated",
        biomarker_file=str(tmp_path / "does_not_exist.txt"),
        log=lambda *_: None,
    )
    sources = summary["name_source"]
    assert sources["src0@C0"] == "lead"
    assert sources["src2@C2"] == "lead"
    # the 'unassigned' cluster must NOT be named 'unassigned'
    assert sources["src1@C1"] == "top-marker"
    assert not summary["names"]["src1@C1"].lower().startswith("unassigned")


def test_annotate_rejects_a_missing_cluster_key(tmp_path):
    adata = _toy_adata()
    with pytest.raises(KeyError):
        A.annotate_pruned_clusters(adata, str(tmp_path / "x"), cluster_key="absent",
                                   log=lambda *_: None)


def test_annotate_rejects_a_missing_lead_column(tmp_path):
    adata = _toy_adata()
    with pytest.raises(KeyError):
        A.annotate_pruned_clusters(adata, str(tmp_path / "y"), cluster_key="pruned",
                                   lead_annotation="not_a_column", log=lambda *_: None)


def test_cli_defaults_annotate_on_and_flags_resolve():
    from altanalyze3.components.sctriangulate import cli

    parser = cli.build_parser()
    base = ["--h5ad", "x.h5ad", "--query", "a,b", "--outdir", "o"]
    args = parser.parse_args(base)
    assert args.annotate is True
    assert args.annotate_top_n == 60
    assert args.annotate_cells_per_cluster == 100
    assert args.annotate_hopach_distance == "cor"
    assert args.annotate_max_fdr == 1e-5
    assert args.annotate_min_overlap == 5
    assert args.annotate_lead is None
    assert parser.parse_args(base + ["--no-annotate"]).annotate is False
