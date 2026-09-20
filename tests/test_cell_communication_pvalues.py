"""Cell communication reports a p-value only when it ran a test.

With fewer than 2 samples in a group the code used to write a normalised rank of the
effect size into `pval`, then pass it to Benjamini-Hochberg. The `fdr` column that came
out was a monotone transform of the ranking and carried no error rate, while looking
like one a reader could filter at 0.05.

Run:  python -m pytest tests/test_cell_communication_pvalues.py -q
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from altanalyze3.components.cellHarmony.flask.pipeline import _bh_fdr  # noqa: E402


def test_bh_keeps_a_missing_pvalue_missing():
    out = _bh_fdr([0.01, float("nan"), 0.04])
    assert math.isnan(out[1])
    assert all(np.isfinite(out[i]) for i in (0, 2))


def test_bh_ranks_only_the_tested_features():
    """A missing value must not take a BH rank and shift everything else."""
    with_gap = _bh_fdr([0.01, float("nan"), 0.04])
    without = _bh_fdr([0.01, 0.04])
    assert with_gap[0] == without[0]
    assert with_gap[2] == without[1]


def test_bh_on_all_missing_returns_all_missing():
    out = _bh_fdr([float("nan"), float("nan")])
    assert len(out) == 2 and all(math.isnan(v) for v in out)


def test_bh_still_adjusts_a_normal_vector():
    out = _bh_fdr([0.01, 0.02, 0.03])
    assert out == sorted(out)
    assert all(0.0 <= v <= 1.0 for v in out)


def test_the_rank_no_longer_lives_in_pval():
    source = (REPO / "altanalyze3" / "components" / "cellHarmony" / "flask"
              / "pipeline.py").read_text()
    assert 'detailed["pval"] = (\n            detailed["abs_delta_score"].rank' not in source
    assert 'detailed["effect_rank"] = (' in source
    assert 'detailed["test_applied"] = bool(tested)' in source
    assert 'no_test_reason' in source


def test_no_placeholder_pvalue_of_one():
    source = (REPO / "altanalyze3" / "components" / "cellHarmony" / "flask"
              / "pipeline.py").read_text()
    block = source[source.index("def _run_cell_communication_differential"):]
    block = block[:block.index("def ", 100)]
    assert "pval = 1.0" not in block
