#!/usr/bin/env python3
"""Regression tests for the ICGS3 -> MarkerFinder scaling interface.

Two defects motivate these tests.

1. markerFinder.detect_input_scaling classifies a matrix as raw counts only when
   frac_integer == 1.0 exactly. AltAnalyze3 ambient RNA correction subtracts a
   fractional ambient profile and therefore writes NON-integer counts. Those
   counts were classified 'unknown', and enforce_input_scaling then refused to
   scale them even when the caller passed scale_data=True.

2. ICGS.run_canonical_heatmap branches on icgs3_normalization_mode and had no
   branch for 'cp10k-log1p', which resolve_normalization returns for every RNA
   run. The run therefore reached the final else, passed a log-scale X unscaled,
   and died at the last step.

The tests below must keep every existing refusal in force. A matrix that is
log-transformed without depth normalization, or mean-centred, still carries no
recoverable counts and must still be refused.
"""
import numpy as np
import pytest
import scipy.sparse as sp

from altanalyze3.components.cellHarmony import markerFinder as mf


RNG = np.random.default_rng(0)


def _counts(n_cells=400, n_genes=300, seed=0):
    rng = np.random.default_rng(seed)
    depth = rng.integers(2000, 20000, size=n_cells)
    profile = rng.dirichlet(np.ones(n_genes) * 0.3)
    X = np.vstack([rng.multinomial(int(d), profile) for d in depth]).astype(np.float64)
    return X


def _ambient_corrected(X, rho=0.2):
    """Reproduce the ambient_subtract model: subtract rho * n_j * b_g, floor at 0."""
    totals = X.sum(axis=1, keepdims=True)
    ambient = X.sum(axis=0) / X.sum()
    corrected = np.maximum(X - rho * totals * ambient[None, :], 0.0)
    return corrected


# --------------------------------------------------------------------------
# Behaviour that must NOT change
# --------------------------------------------------------------------------
def test_integer_counts_still_classified_as_counts():
    report = mf.detect_input_scaling(sp.csr_matrix(_counts()))
    assert report["status"] == "counts"
    assert report["frac_integer"] == 1.0


def test_log2_cp10k_still_ok():
    X = _counts()
    cp10k = X / X.sum(axis=1, keepdims=True) * 1e4
    report = mf.detect_input_scaling(sp.csr_matrix(np.log2(cp10k + 1.0)))
    assert report["status"] == "ok"


def test_linear_cp10k_still_linear_normalized():
    X = _counts()
    cp10k = X / X.sum(axis=1, keepdims=True) * 1e4
    report = mf.detect_input_scaling(sp.csr_matrix(cp10k))
    assert report["status"] == "linear_normalized"


def test_log_without_depth_normalization_still_refused():
    """The case the guard exists for: a log matrix carries no recoverable counts.

    A log of counts reports best_transform == 'identity', because taking a log compresses
    the depth variation and leaves the untransformed row sums the most constant of the
    candidates. Any repair rule keyed on best_transform would therefore accept this matrix
    by mistake. The rule must be the magnitude ceiling instead.
    """
    logged = sp.csr_matrix(np.log1p(_counts()))
    report = mf.detect_input_scaling(logged)
    assert report["max"] <= mf._LOG_INVERSE_CEILING
    assert report["best_transform"] == "identity"
    assert report["status"] == "unknown"
    with pytest.raises(mf.MarkerFinderInputError):
        mf.enforce_input_scaling(logged, scale_data=True, verbose=False)


def test_mean_centred_still_refused():
    X = _counts()
    centred = sp.csr_matrix(X - X.mean(axis=0, keepdims=True))
    report = mf.detect_input_scaling(centred)
    assert report["status"] == "prescaled"
    with pytest.raises(mf.MarkerFinderInputError):
        mf.enforce_input_scaling(centred, scale_data=True, verbose=False)


def test_refusal_without_scale_data_still_stands():
    """Un-normalized input must still be refused when the caller does not ask to scale."""
    corrected = sp.csr_matrix(_ambient_corrected(_counts()))
    with pytest.raises(mf.MarkerFinderInputError):
        mf.enforce_input_scaling(corrected, scale_data=False, verbose=False)


# --------------------------------------------------------------------------
# Defect 1: ambient-corrected, non-integer counts
# --------------------------------------------------------------------------
def test_ambient_corrected_counts_are_recognised_as_linear_counts():
    corrected = sp.csr_matrix(_ambient_corrected(_counts()))
    report = mf.detect_input_scaling(corrected)
    assert report["frac_integer"] < 1.0, "test fixture must be non-integer"
    assert report["max"] > mf._LOG_INVERSE_CEILING, "no log could produce values this large"
    assert report["status"] == "counts_like", report


def test_ambient_corrected_counts_scale_under_scale_data():
    corrected = sp.csr_matrix(_ambient_corrected(_counts()))
    scaled, report = mf.enforce_input_scaling(corrected, scale_data=True, verbose=False)
    assert mf.detect_input_scaling(scaled)["status"] == "ok"


def test_scaling_matches_explicit_pre_scaling():
    """The repaired path must equal calling scale_expression_matrix directly."""
    corrected = sp.csr_matrix(_ambient_corrected(_counts()))
    auto, _ = mf.enforce_input_scaling(corrected, scale_data=True, verbose=False)
    manual, _ = mf.scale_expression_matrix(corrected)
    a = auto.toarray() if sp.issparse(auto) else np.asarray(auto)
    m = manual.toarray() if sp.issparse(manual) else np.asarray(manual)
    assert np.allclose(a, m, atol=1e-10)


def test_marker_finder_runs_on_ambient_corrected_counts():
    """End to end: MarkerFinder must score ambient-corrected counts under scale_data."""
    counts = _counts(n_cells=200, n_genes=120)
    counts[:100, :20] *= 6          # give group A a real marker block
    corrected = _ambient_corrected(counts)
    genes = [f"G{i}" for i in range(corrected.shape[1])]
    groups = ["A"] * 100 + ["B"] * 100
    r_df, p_df = mf.marker_finder(sp.csr_matrix(corrected), groups,
                                  gene_names=genes, scale_data=True)
    # MarkerFinder drops zero-variance features; ambient subtraction can empty a few genes.
    assert r_df.shape[1] == 2
    assert 100 <= r_df.shape[0] <= 120
    assert np.isfinite(r_df.to_numpy()).any()
    top_a = r_df["A"].sort_values(ascending=False).index[:20]
    assert sum(g in {f"G{i}" for i in range(20)} for g in top_a) >= 15


def test_marker_finder_refuses_ambient_corrected_counts_without_scale_data():
    """The guard must still stop an un-normalized matrix when the caller does not opt in."""
    corrected = _ambient_corrected(_counts(n_cells=120, n_genes=80))
    with pytest.raises(mf.MarkerFinderInputError):
        mf.marker_finder(sp.csr_matrix(corrected), ["A"] * 60 + ["B"] * 60,
                         scale_data=False, validate_scaling=True)


# --------------------------------------------------------------------------
# Defect 2: every normalization mode ICGS3 can emit must be handled
# --------------------------------------------------------------------------
def test_every_normalization_mode_is_handled_by_the_heatmap_branch():
    import inspect
    from altanalyze3.components.clustering import ICGS

    modes = [m for m in ICGS.NORMALIZATION_MODES if m != "auto"]
    source = inspect.getsource(ICGS.run_canonical_heatmap)
    for mode in modes:
        assert f'"{mode}"' in source, (
            f"run_canonical_heatmap has no branch for normalization mode {mode!r}; "
            f"resolve_normalization can return it for any RNA run"
        )


# --------------------------------------------------------------------------
# Root cause: ICGS3 must recognise ambient-corrected counts AS counts, so that
# it preserves layers['counts'] and MarkerFinder receives the log2 depth-
# normalized values it asks for.
# --------------------------------------------------------------------------
def test_icgs3_calls_raw_integer_counts_counts():
    from altanalyze3.components.clustering.ICGS import infer_expression_scale
    r = infer_expression_scale(sp.csr_matrix(_counts()))
    assert r["verdict"] == "counts"


def test_icgs3_calls_ambient_corrected_counts_counts():
    """The defect: these were called TotalVI-like denoised expression."""
    from altanalyze3.components.clustering.ICGS import infer_expression_scale, DEPTH_CV_MIN
    r = infer_expression_scale(sp.csr_matrix(_ambient_corrected(_counts())))
    assert r["frac_integer"] < 0.99, "fixture must be non-integer"
    assert r["max"] > 50
    assert r["cv_cell_total"] >= DEPTH_CV_MIN
    assert r["verdict"] == "counts", r["reason"]


def test_icgs3_still_protects_depth_corrected_model_output():
    """Depth-corrected linear values must NOT be renormalized. Protection preserved."""
    from altanalyze3.components.clustering.ICGS import infer_expression_scale, DEPTH_CV_MIN
    X = _counts()
    denoised = X / X.sum(axis=1, keepdims=True) * 1e4   # depth removed, max well above 50
    r = infer_expression_scale(sp.csr_matrix(denoised))
    assert r["max"] > 50
    assert r["cv_cell_total"] < DEPTH_CV_MIN
    assert r["verdict"] == "linear_non_integer", r["reason"]


def test_icgs3_log_and_centered_verdicts_unchanged():
    from altanalyze3.components.clustering.ICGS import infer_expression_scale
    X = _counts()
    cp10k = X / X.sum(axis=1, keepdims=True) * 1e4
    assert infer_expression_scale(sp.csr_matrix(np.log2(cp10k + 1.0)))["verdict"] == "log"
    assert infer_expression_scale(sp.csr_matrix(X - X.mean(axis=0, keepdims=True)))["verdict"] == "centered"
