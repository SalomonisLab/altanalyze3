"""The population selector offers cell states, never the pooled whole-sample test.

`Pooled overall` answers a different question from a per-cell-state test. The pipeline
used to substitute it whenever the per-state table came back empty, which put it in the
population selector as though it were a cell state. The Regulatory network then had no
population to draw and reported a missing modality instead of a missing result.

Run:  python -m pytest tests/test_cell_state_populations_only.py -q
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from altanalyze3.components.cellHarmony.webapp.app import (  # noqa: E402
    _POOLED_OVERALL_LABEL,
    _cell_state_populations,
)


def test_pooled_label_never_reaches_a_population_list():
    ordered = ["AT2-prolif", _POOLED_OVERALL_LABEL, "CAP2", "AM"]
    assert _cell_state_populations(ordered) == ["AT2-prolif", "CAP2", "AM"]


def test_a_pooled_only_run_offers_no_population():
    # What job f17d2aed's GRN run produced: the pooled row and nothing else.
    assert _cell_state_populations([_POOLED_OVERALL_LABEL]) == []


def test_empty_and_blank_values_are_dropped():
    assert _cell_state_populations(["", "AT1", ""]) == ["AT1"]
    assert _cell_state_populations([]) == []


def test_order_is_preserved():
    ordered = ["Ciliated-TB", "AT2", _POOLED_OVERALL_LABEL, "PArEC"]
    assert _cell_state_populations(ordered) == ["Ciliated-TB", "AT2", "PArEC"]


def test_the_pipeline_no_longer_substitutes_a_pooled_result():
    source = (REPO / "altanalyze3" / "components" / "cellHarmony" / "flask"
              / "pipeline.py").read_text()
    assert "detailed_deg = pooled_detail" not in source
    assert 'pooled_detail["population"] = _POOLED_OVERALL_LABEL' not in source
    assert "No cell state produced a differential result" in source
