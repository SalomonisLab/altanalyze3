"""Prove the --annotate-lead-mode enrichment naming picks the enriched label, not the biggest.

The failure the mode exists to stop: a reference state that covers most of the object wins every
cluster under dominant mode, because a cluster's most frequent label is that state even when the
cluster holds no more of it than chance predicts.
"""

import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.sctriangulate.annotate import (
    _benjamini_hochberg,
    _predictions_from_lead_enrichment,
    lead_enrichment_table,
)


def _fixture():
    """900 cells. Ferchen 'ILC1_3+NKP' covers 70% of everything, so it is every cluster's mode.

    Cluster C1 is genuinely 100% St Jude 'DN'. Cluster C2 is genuinely 100% Ferchen 'CLP1-c'.
    Cluster C3 holds the background mixture and is enriched for nothing.
    """
    n = 900
    cluster = np.array(["C1"] * 150 + ["C2"] * 150 + ["C3"] * 600)
    ferchen = np.array(
        ["ILC1_3+NKP"] * 150            # C1: the background state, at its background share
        + ["CLP1-c"] * 150              # C2: a small state, wholly inside C2
        + ["ILC1_3+NKP"] * 480 + ["EILP"] * 120  # C3: background mixture
    )
    stjude = np.array(
        ["DN"] * 150                    # C1: a small state, wholly inside C1
        + ["DP-blast"] * 150
        + ["DP-blast"] * 450 + ["ISP"] * 150
    )
    obs = pd.DataFrame({"Ferchen": ferchen, "StJude": stjude},
                       index=[f"BC{i:04d}" for i in range(n)])
    return pd.Series(cluster, index=obs.index), obs


def test_bh_matches_a_hand_computed_example():
    p = np.array([0.01, 0.02, 0.03, 0.04, 0.05])
    adj = _benjamini_hochberg(p)
    assert np.allclose(adj, [0.05, 0.05, 0.05, 0.05, 0.05])
    assert np.all(np.diff(adj) >= -1e-12)
    assert _benjamini_hochberg(np.array([])).size == 0
    assert np.all(_benjamini_hochberg(np.array([0.9, 0.9])) <= 1.0)


def test_enrichment_beats_the_background_state():
    labels, obs = _fixture()
    clusters = ["C1", "C2", "C3"]
    table = lead_enrichment_table(labels, obs, ["Ferchen", "StJude"], clusters=clusters)
    assert set(table["reference"]) == {"Ferchen", "StJude"}
    assert (table["overlap"] <= table["n_cluster"]).all()
    assert (table["overlap"] <= table["n_reference_label"]).all()

    picks = _predictions_from_lead_enrichment(
        table, clusters, max_fdr=0.05, min_overlap=10, log=lambda *_: None)

    # Dominant mode would name C1 by its most frequent Ferchen label, the background state.
    assert obs.loc[labels == "C1", "Ferchen"].value_counts().index[0] == "ILC1_3+NKP"
    # Enrichment mode names it from the St Jude state that is exclusive to it.
    assert picks["C1"]["label"] == "DN"
    assert picks["C1"]["reference"] == "StJude"
    assert picks["C2"]["label"] == "CLP1-c"
    assert picks["C2"]["reference"] == "Ferchen"
    assert picks["C1"]["fold_enrichment"] > 1.0
    assert picks["C2"]["fold_enrichment"] > 1.0


def test_uninformative_reference_labels_are_excluded():
    labels = pd.Series(["C1"] * 50 + ["C2"] * 50, index=[f"B{i}" for i in range(100)])
    obs = pd.DataFrame({"Ref": ["unassigned"] * 50 + ["ProperState"] * 50}, index=labels.index)
    table = lead_enrichment_table(labels, obs, ["Ref"], clusters=["C1", "C2"])
    assert "unassigned" not in set(table["reference_label"])
    assert set(table["cluster"]) == {"C2"}
    assert int(table["n_annotated"].iloc[0]) == 50


def test_min_overlap_rejects_a_tiny_but_significant_label():
    labels, obs = _fixture()
    table = lead_enrichment_table(labels, obs, ["Ferchen", "StJude"], clusters=["C1", "C2", "C3"])
    strict = _predictions_from_lead_enrichment(
        table, ["C1"], max_fdr=0.05, min_overlap=10_000, log=lambda *_: None)
    assert strict == {}


def test_fold_enrichment_is_observed_over_expected():
    labels, obs = _fixture()
    table = lead_enrichment_table(labels, obs, ["StJude"], clusters=["C1"])
    row = table[table["reference_label"] == "DN"].iloc[0]
    expected = row["n_cluster"] * row["n_reference_label"] / row["n_annotated"]
    assert np.isclose(row["expected_overlap"], expected)
    assert np.isclose(row["fold_enrichment"], row["overlap"] / expected)


# --- MarkerFinder scaling guard on a CITE-seq object -----------------------------------------

def _cite_object(n_cells=200, n_rna=400, n_adt=20, depth_normalized=True, seed=0):
    import anndata
    from scipy.sparse import csr_matrix
    import scanpy as sc
    rng = np.random.default_rng(seed)
    counts = rng.poisson(1.2, size=(n_cells, n_rna)).astype(np.float64)
    counts[:, 0] += rng.integers(1, 50, n_cells)          # per-cell depth spread
    rna = anndata.AnnData(X=csr_matrix(counts))
    if depth_normalized:
        sc.pp.normalize_total(rna, target_sum=1e4)
        sc.pp.log1p(rna)
    adt = np.log1p(rng.gamma(2.0, 20.0, size=(n_cells, n_adt)))
    X = np.hstack([rna.X.toarray(), adt])
    var = pd.DataFrame(index=[f"Gene{i}" for i in range(n_rna)]
                             + [f"AB_CD{i}" for i in range(n_adt)])
    return anndata.AnnData(X=csr_matrix(X), var=var,
                           obs=pd.DataFrame(index=[f"BC{i}" for i in range(n_cells)]))


def test_guard_is_switched_off_only_when_the_rna_block_alone_passes():
    from altanalyze3.components.sctriangulate.annotate import _resolve_marker_validate_scaling
    from altanalyze3.components.cellHarmony.markerFinder import detect_input_scaling

    good = _cite_object(depth_normalized=True)
    rna_only = good[:, [v for v in good.var_names if not str(v).startswith("AB_")]]
    assert detect_input_scaling(rna_only.X)["status"] == "ok"
    assert detect_input_scaling(good.X)["status"] != "ok"       # the AB_ block breaks it
    assert _resolve_marker_validate_scaling(good, "AB_", None, False, lambda *_: None) is False


def test_guard_stays_on_when_the_rna_block_itself_is_not_depth_normalized():
    from altanalyze3.components.sctriangulate.annotate import _resolve_marker_validate_scaling
    bad = _cite_object(depth_normalized=False)
    assert _resolve_marker_validate_scaling(bad, "AB_", None, False, lambda *_: None) is True


def test_guard_stays_on_for_an_rna_only_object():
    import anndata
    from scipy.sparse import csr_matrix
    from altanalyze3.components.sctriangulate.annotate import _resolve_marker_validate_scaling
    rng = np.random.default_rng(0)
    a = anndata.AnnData(X=csr_matrix(rng.random((30, 20)).astype(np.float32)),
                        var=pd.DataFrame(index=[f"G{i}" for i in range(20)]))
    assert _resolve_marker_validate_scaling(a, "AB_", None, False, lambda *_: None) is True
