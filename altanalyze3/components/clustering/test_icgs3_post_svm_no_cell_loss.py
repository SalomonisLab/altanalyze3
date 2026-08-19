"""The post-SVM marker pass selects genes; it must never delete cells (ICGS2 fidelity).

ICGS_NMF.py:1134-1141 runs MarkerFinder after the SVM and feeds the result to
ExpandSampleClusters.filterRows, which filters ROWS (genes) of the expression matrix. The
cluster set `name` is never revised after the SVM. A regression in ICGS3 subsetted the
per-barcode table instead, deleting 4,575 of 18,388 cells on a mouse thymus run.
"""

import inspect
import re

from altanalyze3.components.clustering import ICGS


def _post_svm_block() -> str:
    src = inspect.getsource(ICGS.run_nmf_marker_svm)
    start = src.index("post-SVM MarkerFinder gene pool")
    return src[start:]


def test_post_svm_pass_does_not_subset_the_cell_table():
    block = _post_svm_block()
    offending = re.findall(
        r"final_clusters\s*=\s*final_clusters\[[^\]]*isin\(\s*robust2\s*\)[^\]]*\]", block
    )
    assert not offending, (
        "the post-SVM marker pass subsets final_clusters, which is indexed by barcode, so it "
        f"deletes cells: {offending}"
    )


def test_pre_svm_cluster_fitness_gate_is_still_applied():
    """The gate ICGS2 DOES apply runs before the SVM, where it costs no cell."""
    src = inspect.getsource(ICGS.run_nmf_marker_svm)
    pre = src[: src.index("SVM reclassification target")]
    assert "nmf_robust = nmf_clusters[nmf_clusters[\"cluster\"].isin(robust)]" in pre
    assert "generate_train_data(expr_markers_sampled, nmf_robust)" in pre


def test_starved_clusters_are_reported_not_removed():
    block = _post_svm_block()
    assert "starved" in block, "clusters without markers must still be reported"
    assert "ICGS2 keeps them" in block


def test_svm_min_decision_score_remains_the_only_cell_level_gate():
    """One cell-level filter after classification, matching ExpandSampleClusters.py:307."""
    src = inspect.getsource(ICGS._classification_frame_from_scores)
    assert "keep = win_score > float(min_decision_score)" in src
