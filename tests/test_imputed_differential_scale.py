"""Known synthetic abundances verify folds; never rerun a biological analysis."""
import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests

from altanalyze3.components.cellHarmony.cellHarmony_differential import (
    _compute_pseudobulk_log2fc, _prepare_scanpy_rank_input, _rank_genes_scanpy,
    _moderated_t_test, run_de_for_comparisons)
from altanalyze3.components.cellHarmony.flask.pipeline import (
    _attach_imputed_expression_metadata, _read_differential_h5ad,
    _build_pseudobulk_differential_adata)
from altanalyze3.components.cellHarmony.flask.disk_predictions import PredictionWriter
from altanalyze3.components.cellHarmony.imputed_scale import (
    prediction_encoding, encode_predictions, inverse_predictions)
from altanalyze3.components.cellHarmony.imputed_pseudobulk import aggregate_imputed_predictions

ENCODINGS = [("linear", None, 0), ("log1p", np.e, 1), ("log1p", 10, 1),
             ("log2p1", 2, 1), ("log2", 2, 0), ("log2", 2, 10)]


def data(encoding, *, n_features=3):
    scale, base, offset = encoding
    rng = np.random.default_rng(417)
    linear = rng.uniform(.01, .09, (24, n_features))
    # Every corresponding case value is exactly 2x its control counterpart.
    linear[:12] = linear[12:] * 2
    obs = pd.DataFrame({"condition": pd.Categorical(["case"] * 12 + ["control"] * 12),
                        "state": "state", "sample": [f"donor{i // 2}" for i in range(24)]},
                       index=[f"c{i}" for i in range(24)])
    a = ad.AnnData(encode_predictions(linear, encoding), obs=obs,
                   var=pd.DataFrame(index=[f"feature{i}" for i in range(n_features)]),
                   uns={"modality": "adt" if scale != "linear" else "grn_tf", "expression_scale": scale})
    if base is not None:
        a.uns.update(log_base=base, log_pseudocount=offset)
    if scale == "log1p":
        a.uns["log1p"] = {"base": base}
    return a, linear


@pytest.mark.parametrize("encoding", ENCODINGS)
def test_exact_folds_use_source_inverse_without_adding_one(encoding):
    a, linear = data(encoding)
    # Deliberately stale auxiliaries must not override declared prediction X.
    a.raw = ad.AnnData(np.full(a.shape, 800.), obs=a.obs.copy(), var=a.var.copy())
    a.layers["counts"] = np.full(a.shape, 900.)
    before = a.X.copy()
    result = _compute_pseudobulk_log2fc(a, "condition", "case", "control")
    np.testing.assert_allclose(result.log2fc, 1., atol=1e-12)
    np.testing.assert_allclose(result.case_mean_linear, linear[:12].mean(0), atol=1e-12)
    np.testing.assert_allclose(result.control_mean_linear, linear[12:].mean(0), atol=1e-12)
    assert result.attrs["fold_convention"]["fold_pseudocount"] == 0
    np.testing.assert_array_equal(a.X, before)


@pytest.mark.parametrize("encoding", ENCODINGS)
def test_raw_tests_unchanged_and_reported_folds_are_exact(encoding):
    a, _ = data(encoding)
    reference = a.copy()
    sc.tl.rank_genes_groups(reference, groupby="condition", groups=["case"],
                           reference="control", method="wilcoxon", use_raw=False)
    names, _, fc, p = _rank_genes_scanpy(a, "condition", "case", "control", "wilcoxon")
    expected_p = pd.Series(reference.uns["rank_genes_groups"]["pvals"]["case"],
                           index=reference.uns["rank_genes_groups"]["names"]["case"])
    np.testing.assert_array_equal(p.sort_index(), expected_p.sort_index())
    np.testing.assert_allclose(fc, 1., atol=1e-12)
    np.testing.assert_allclose(a.uns["rank_genes_groups"]["logfoldchanges"]["case"], 1., atol=1e-12)
    assert "scanpy_logfoldchanges_approximate" in a.uns["rank_genes_groups"]
    # Independently reconstruct the existing moderated test from source X.
    case, control = a.X[:12], a.X[12:]
    s2 = (case.var(0, ddof=1) + control.var(0, ddof=1)) / 2
    shrunk = .8 * s2 + .2 * np.median(s2)
    t = (case.mean(0) - control.mean(0)) / np.maximum(np.sqrt(shrunk / 6), 1e-12)
    expected = 2 * stats.t.sf(np.abs(t), df=22)
    result, _ = _moderated_t_test(a, "condition", "case", "control", "state")
    np.testing.assert_allclose(result.pval, expected, atol=1e-12)
    np.testing.assert_allclose(result.log2fc, 1., atol=1e-12)


@pytest.mark.parametrize("encoding", ENCODINGS)
def test_mean_aggregation_then_fold_matches_known_abundances(encoding):
    a, linear = data(encoding)
    pb, audit = aggregate_imputed_predictions(a, population_col="state", sample_col="sample",
                                             covariate_col="condition", min_cells=2)
    expected = linear.reshape(12, 2, 3).mean(1)
    # Existing H5AD prediction storage is float32; inversion with a source
    # offset of 10 amplifies its quantization at abundances below 0.1.
    np.testing.assert_allclose(inverse_predictions(pb.X, prediction_encoding(pb)), expected, atol=1e-6)
    np.testing.assert_allclose(_compute_pseudobulk_log2fc(pb, "condition", "case", "control").log2fc,
                               1., atol=5e-5)
    assert audit["model_recomputed"] is False
    assert audit["aggregated_cells"] == a.n_obs
    assert list(pb.var_names) == list(a.var_names)


@pytest.mark.parametrize("encoding", ENCODINGS)
def test_ram_export_and_stream_export_have_identical_linear_layers(tmp_path, encoding):
    a, _ = data(encoding)
    a.X = a.X.astype(np.float32)
    info = _attach_imputed_expression_metadata(a, a.X, expression_scale=encoding[0],
                                               log_base=encoding[1], log_pseudocount=encoding[2])
    assert info["target_encoding_complete"]
    path = tmp_path / "predictions.h5ad"
    writer = PredictionWriter(path, a.obs, a.var, dict(a.uns), {},
                              expression_scale=encoding[0], log_base=encoding[1])
    try:
        writer.append(a.X[:7]); writer.append(a.X[7:]); writer.close()
    except BaseException:
        writer.abort()
        raise
    saved = ad.read_h5ad(path)
    np.testing.assert_array_equal(saved.X, a.X)
    np.testing.assert_array_equal(saved.layers["counts"], a.layers["counts"])


@pytest.mark.parametrize("encoding", ENCODINGS)
def test_disk_feature_blocks_match_ram_and_never_renormalize(tmp_path, encoding):
    a, _ = data(encoding, n_features=520)  # Exercise multiple feature blocks.
    if encoding[0] == "linear":
        a.X *= 10000  # High linear TF scores must not become CP10k log1p.
    assert _prepare_scanpy_rank_input(a) is a
    before = a.X.copy()
    path = tmp_path / "source.h5ad"
    a.write_h5ad(path)
    disk = _read_differential_h5ad(path, disk_backed=True)
    try:
        ram = _rank_genes_scanpy(a, "condition", "case", "control", "wilcoxon")
        backed = _rank_genes_scanpy(disk, "condition", "case", "control", "wilcoxon")
        for x, y in zip(ram[1:], backed[1:]):
            np.testing.assert_allclose(x.sort_index(), y.sort_index(), atol=1e-12)
        pb_ram = _moderated_t_test(a, "condition", "case", "control", "state")[0]
        pb_disk = _moderated_t_test(disk, "condition", "case", "control", "state")[0]
        pd.testing.assert_frame_equal(pb_ram, pb_disk, atol=1e-12)
    finally:
        disk._analysis_h5_handle.close()
    np.testing.assert_array_equal(a.X, before)


@pytest.mark.parametrize("encoding", [ENCODINGS[0], ENCODINGS[1], ENCODINGS[3]])
def test_zero_means_are_honest_and_all_feature_identities_retained(encoding):
    a, linear = data(encoding)
    linear[:12, 0] = 0; linear[12:, 1] = 0; linear[:, 2] = 0
    a.X = encode_predictions(linear, encoding)
    result = _compute_pseudobulk_log2fc(a, "condition", "case", "control")
    assert np.isneginf(result.log2fc.iloc[0])
    assert np.isposinf(result.log2fc.iloc[1])
    assert np.isnan(result.log2fc.iloc[2])
    assert list(result.fold_status) == ["case_mean_zero", "control_mean_zero", "both_means_zero"]
    assert list(result.index) == list(a.var_names)


@pytest.mark.parametrize("modality", ["adt", "metabolite", "lipid", "grn", "grn_tf"])
@pytest.mark.parametrize("encoding", [ENCODINGS[0], ENCODINGS[4]])
def test_all_imputed_modalities_have_full_panel_bh_despite_low_or_negative_X(modality, encoding):
    a, _ = data(encoding); a.uns["modality"] = modality
    if encoding[0] == "linear":
        a.X *= .1
    assert (a.X < .1).all()  # Every feature would fail the old eligibility rule.
    names, fdr, _, p = _rank_genes_scanpy(a, "condition", "case", "control", "wilcoxon")
    assert set(names) == set(a.var_names)
    np.testing.assert_allclose(fdr, multipletests(p, method="fdr_bh")[1])
    result, tested = _moderated_t_test(a, "condition", "case", "control", "state")
    assert tested == a.n_vars
    assert list(result.gene) == list(a.var_names)
    np.testing.assert_allclose(result.fdr, multipletests(result.pval, method="fdr_bh")[1])


def test_pooled_algorithm_remains_moderated_for_cell_comparisons():
    a, _ = data(ENCODINGS[1])
    expected, _ = _moderated_t_test(a, "condition", "case", "control", "pooled_overall", store_full_log=False)
    result = run_de_for_comparisons(a, "state", "condition", "case", "control", method="wilcoxon",
                                    alpha=1.0, fc_thresh=1.5, min_cells_per_group=4, use_rawp=True)
    np.testing.assert_allclose(result["pooled_overall"].sort_index().fdr,
                               expected.set_index("gene").sort_index().fdr)


def test_undeclared_aml_inverse_blocks_both_tests_and_aggregation(tmp_path):
    a, _ = data(ENCODINGS[0]); a.uns = {"modality": "metabolite", "expression_scale": "log2"}
    _attach_imputed_expression_metadata(a, a.X, expression_scale="log2", log_base=2)
    # Export can preserve source X for display, without a fabricated inverse.
    assert "counts" not in a.layers
    writer = PredictionWriter(tmp_path / "ambiguous.h5ad", a.obs, a.var, dict(a.uns), {}, expression_scale="log2")
    writer.append(a.X); writer.close()
    saved = ad.read_h5ad(tmp_path / "ambiguous.h5ad")
    np.testing.assert_array_equal(saved.X, a.X.astype(np.float32))
    assert "counts" not in saved.layers
    for fn in [lambda: _rank_genes_scanpy(a, "condition", "case", "control", "wilcoxon"),
               lambda: _moderated_t_test(a, "condition", "case", "control", "state"),
               lambda: aggregate_imputed_predictions(a, population_col="state", sample_col="sample",
                                                      covariate_col="condition", min_cells=2)]:
        with pytest.raises(ValueError, match="preprocessing protocol"):
            fn()


def test_signed_linear_export_is_preserved_but_abundance_fold_is_rejected(tmp_path):
    a, _ = data(ENCODINGS[0]); a.X[0, 0] = -2
    _attach_imputed_expression_metadata(a, a.X, expression_scale="linear")
    assert a.layers["counts"][0, 0] == -2
    writer = PredictionWriter(tmp_path / "signed.h5ad", a.obs, a.var, dict(a.uns), {})
    writer.append(a.X); writer.close()
    assert ad.read_h5ad(tmp_path / "signed.h5ad").layers["counts"][0, 0] == -2
    with pytest.raises(ValueError, match="Signed activity scores"):
        _compute_pseudobulk_log2fc(a, "condition", "case", "control")


@pytest.mark.parametrize("metadata", [
    {"expression_scale": "native_relative_log2", "log1p": {"base": None}},
    {"expression_scale": "log2p1", "log_base": 10},
    {"expression_scale": "log1p", "log_base": 0},
    {"expression_scale": "log1p", "log_pseudocount": 0},
    {"expression_scale": "log2", "log_pseudocount": -1},
])
def test_conflicting_encodings_stop_instead_of_guessing(metadata):
    a, _ = data(ENCODINGS[0]); a.uns = metadata
    with pytest.raises(ValueError):
        prediction_encoding(a)


def test_auxiliary_grn_pseudobulks_declare_linear_scale():
    a, _ = data(ENCODINGS[0]); a.obs["Library"] = a.obs["sample"]
    keys = [f"donor{i}|state" for i in range(12)]
    predictions = pd.DataFrame(np.arange(36).reshape(12, 3) + .01, index=keys, columns=a.var_names)
    pb = _build_pseudobulk_differential_adata(a, predictions, "state", feature_type="GRN-TF", modality_id="grn_tf")
    assert pb.uns["expression_scale"] == "linear"
    np.testing.assert_array_equal(pb.layers["counts"], pb.X)
    assert _prepare_scanpy_rank_input(pb) is pb


@pytest.mark.parametrize("pseudobulk", [False, True])
def test_full_comparison_outputs_use_exact_fold_and_annotate_means(pseudobulk):
    a, _ = data(ENCODINGS[1])
    if pseudobulk:
        a, _ = aggregate_imputed_predictions(a, population_col="state", sample_col="sample",
                                             covariate_col="condition", min_cells=2)
    result = run_de_for_comparisons(a, "state", "condition", "case", "control", method="wilcoxon",
                                    alpha=1.0, fc_thresh=1.5, min_cells_per_group=4, use_rawp=True)
    detail = result["detailed_deg"]
    assert len(detail) == 3
    np.testing.assert_allclose(detail.log2fc, 1., atol=3e-6)
    np.testing.assert_allclose(result["pooled_overall"].log2fc, 1., atol=3e-6)
    np.testing.assert_allclose(detail.case_mean_linear / detail.control_mean_linear, 2., atol=3e-6)
    assert result["prediction_scale"]["fold_pseudocount"] == 0
