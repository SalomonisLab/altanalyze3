"""Compare optimized ICGS computations with the original scientific calculations."""
import warnings

import numpy as np
import pandas as pd
import pytest
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from altanalyze3.components.clustering import ICGS as icgs


def original_counts(values, threshold):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return (np.nan_to_num(np.corrcoef(values)) > threshold).sum(axis=0)


def expression(dtype=np.float64, order="C"):
    rng = np.random.default_rng(491)
    values = rng.normal(size=(173, 191))
    values[1:15] = values[0] + rng.normal(0, .1, (14, values.shape[1]))
    values[20:30] = -values[0] + rng.normal(0, .1, (10, values.shape[1]))
    values[45] = 0
    values[46] = 7
    if dtype == np.int64:
        values = np.rint(values * 5)
    return np.array(values, dtype=dtype, order=order)


@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.int64])
@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize("threshold", [.4, -.4, 0., .99, -1., 1.])
def test_all_neighbor_counts_match_original_across_tiles(dtype, order, threshold):
    values = expression(dtype, order)
    before = values.copy()
    counts, _ = icgs._pearson_neighbor_counts(values, threshold, block_genes=37)
    np.testing.assert_array_equal(counts, original_counts(values, threshold))
    np.testing.assert_array_equal(values, before)


def test_near_threshold_pairs_are_recalculated_without_changing_selection():
    values = expression()
    x = values[0] - values[0].mean()
    x /= np.linalg.norm(x)
    y = values[50] - values[50].mean()
    y -= np.dot(x, y) * x
    y /= np.linalg.norm(y)
    values[60] = .400002 * x + np.sqrt(1 - .400002**2) * y
    values[61] = .399998 * x + np.sqrt(1 - .399998**2) * y
    counts, refined = icgs._pearson_neighbor_counts(values, .4, block_genes=37)
    assert refined > 0
    np.testing.assert_array_equal(counts, original_counts(values, .4))


def test_exact_float64_boundary_uses_the_original_whole_matrix(monkeypatch):
    values = expression()
    threshold = np.corrcoef(values)[0, 55]
    messages = []
    monkeypatch.setattr(icgs, "_log", messages.append)
    counts, _ = icgs._pearson_neighbor_counts(values, threshold, block_genes=37)
    np.testing.assert_array_equal(counts, original_counts(values, threshold))
    assert any("float64 threshold boundary" in message for message in messages)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
@pytest.mark.parametrize("threshold", [.4, -.4])
def test_nonfinite_rows_retain_original_zero_correlation_semantics(bad, threshold):
    values = expression()
    values[2, 3] = bad
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        counts, _ = icgs._pearson_neighbor_counts(values, threshold)
    np.testing.assert_array_equal(counts, original_counts(values, threshold))


@pytest.mark.parametrize("scale,offset", [(1e-140, 0), (1e140, 0), (1, 1e14)])
def test_extreme_scales_match_original(scale, offset):
    values = expression() * scale + offset
    counts, _ = icgs._pearson_neighbor_counts(values, .4, block_genes=37)
    np.testing.assert_array_equal(counts, original_counts(values, .4))


@pytest.mark.parametrize("scale", [1e-155, 1e-200])
def test_subnormal_variances_retain_original_correlations(scale):
    values = expression()
    values[0] *= scale
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        counts, _ = icgs._pearson_neighbor_counts(values, .4)
    np.testing.assert_array_equal(counts, original_counts(values, .4))


def test_feature_names_order_and_full_input_are_preserved():
    values = expression(np.float32, "F")
    frame = pd.DataFrame(values, index=[f"gene_{i}" for i in range(len(values))])
    before = frame.copy(deep=True)
    names = frame.index[::-1]
    expected = names[original_counts(frame.loc[names].to_numpy(), .4) > 5]
    actual = icgs.fast_intercorrelation(frame, names)
    assert actual.equals(expected)
    pd.testing.assert_frame_equal(frame, before)


def test_vectorized_biomarkers_preserve_probabilities_fdr_labels_and_background(tmp_path):
    catalog = tmp_path / "catalog.tsv"
    pd.DataFrame({
        "Ensembl": ["ENSG000001", "ENSG000002", "ENSG000003", "ENSG000001", "ENSG000004", "ENSG000003"],
        "Gene": ["GATA1", "SPI1", "LST1", "GATA1", "UNKNOWN", "LST1"],
        "Term": ["State A"] * 3 + ["State B"] * 3,
    }).to_csv(catalog, sep="\t", index=False)
    background = ["ENSG000001.1", "ENSG000002", "GATA1", "SPI1", "LST1", "UNMAPPED"]
    markers = pd.DataFrame({"marker": ["GATA1", "ENSG000001.1", "UNMAPPED", "SPI1", "LST1"],
                            "top_cluster": ["c1"] * 3 + ["c2"] * 2})
    before = markers.copy(deep=True)
    cfg = icgs.ICGS3Config(input_paths=[], output_dir=str(tmp_path), biomarker_file=str(catalog))
    labels = icgs.biomarker_enrichment(markers, background, cfg, str(tmp_path))
    observed = pd.read_csv(tmp_path / "GO-Elite/icgs3_biomarker_enrichment.tsv", sep="\t", float_precision="round_trip")
    rows = []
    for cluster, size, overlap_a, overlap_b, genes_a, genes_b in [
        ("c1", 3, 2, 2, "ENSG000001.1,GATA1", "ENSG000001.1,GATA1"),
        ("c2", 2, 2, 1, "LST1,SPI1", "LST1"),
    ]:
        for term, term_size, overlap, genes in [("State A", 5, overlap_a, genes_a), ("State B", 3, overlap_b, genes_b)]:
            rows.append({"cluster": cluster, "term_name": term,
                         "p_value": hypergeom(M=6, n=term_size, N=size).sf(overlap-1),
                         "overlap": overlap, "query_size": size, "term_size": term_size, "overlap_genes": genes})
    expected = pd.DataFrame(rows)
    expected["fdr"] = multipletests(expected.p_value, method="fdr_bh")[1]
    expected = expected.sort_values(["cluster", "fdr", "p_value", "term_name"])
    pd.testing.assert_frame_equal(observed, expected.reset_index(drop=True), check_exact=True)
    expected_labels = icgs.clean_biomarker_prediction_labels(expected.groupby("cluster").head(1)[["cluster", "term_name", "fdr", "overlap"]])
    pd.testing.assert_frame_equal(labels, expected_labels, check_exact=True)
    assert pd.read_csv(tmp_path / "GO-Elite/icgs3_biomarker_background.tsv", sep="\t").feature.tolist() == background
    pd.testing.assert_frame_equal(markers, before)
