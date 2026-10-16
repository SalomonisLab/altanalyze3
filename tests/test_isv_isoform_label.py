"""Every isoform name carries its protein length.

Nathan, 2026-10-15: "include the AA length the isoform names in parentheses for standard ISV
views". The label is built in one place so the read-panel group labels and the structure rows
agree, and `--no-protein-length` writes bare ids.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from altanalyze3.components.visualization.isv_web.isv_web_remote import (  # noqa: E402
    isoform_label,
)


def test_length_in_parentheses():
    assert isoform_label({"isoform_id": "35292108.AML-7", "protein_length": 568}) == \
        "35292108.AML-7 (568aa)"


def test_nmd_sits_inside_the_brackets():
    assert isoform_label({"isoform_id": "X", "protein_length": 144, "nmd_status": "NMD"}) == \
        "X (144aa, NMD)"
    assert isoform_label({"isoform_id": "X", "protein_length": 144,
                          "nmd_status": "Potential-NMD"}) == "X (144aa, NMD)"


def test_final_isoform_id_wins():
    assert isoform_label({"isoform_id": "read", "final_isoform_id": "final",
                          "protein_length": 10}) == "final (10aa)"


def test_no_length_known():
    assert isoform_label({"isoform_id": "X"}) == "X"
    assert isoform_label({"isoform_id": "X", "protein_length": ""}) == "X"
    assert isoform_label({"isoform_id": "X", "nmd_status": "NMD"}) == "X (NMD)"


def test_a_float_or_string_length_still_prints_as_an_integer():
    assert isoform_label({"isoform_id": "X", "protein_length": "568.0"}) == "X (568aa)"
    assert isoform_label({"isoform_id": "X", "protein_length": 568.0}) == "X (568aa)"


def test_switching_it_off_gives_the_bare_id():
    assert isoform_label({"isoform_id": "X", "protein_length": 568}, with_length=False) == "X"


# ---------------------------------------------------------------------------------------------
# '+' in --isoform_filter means a read must contain BOTH patterns. Nathan, 2026-10-15: "generate
# the ISV for reads that contain ENSG00000173621:E2.1-E3.1_66857513 and I1 combined".
from altanalyze3.components.visualization.isoform_structure_view import (  # noqa: E402
    _isoform_matches_filters,
    _normalize_isoform_filter_spec,
)

LRFN4_35292108 = ("U1.1_66856538|E1.1_66856657|E1.1_66856776|E1.1|I1.1|E2.1|E3.1_66857513|"
                  "E3.1_66857544|E3.3_66858374|E3.3|E3.4|E4.1|E4.1_66859679|E4.1_66859746|"
                  "E4.2_66860473").split("|")


def test_plus_ands_the_patterns():
    f = _normalize_isoform_filter_spec(["E2.1|E3.1_66857513 + I1.1"])
    assert f == [[["E2.1", "E3.1_66857513"], ["I1.1"]]]
    assert _isoform_matches_filters(LRFN4_35292108, f)
    assert not _isoform_matches_filters([t for t in LRFN4_35292108 if t != "I1.1"], f)
    assert not _isoform_matches_filters(
        [t for t in LRFN4_35292108 if t != "E3.1_66857513"], f)


def test_separate_arguments_still_or():
    """Without '+', each argument is its own clause, so either one alone keeps the read."""
    f = _normalize_isoform_filter_spec(["E2.1|E3.1_66857513", "I1.1"])
    assert f == [[["E2.1", "E3.1_66857513"]], [["I1.1"]]]
    assert _isoform_matches_filters([t for t in LRFN4_35292108 if t != "I1.1"], f)
    assert _isoform_matches_filters(
        [t for t in LRFN4_35292108 if t != "E3.1_66857513"], f)


def test_only_the_read_panel_label_uses_parentheses():
    """Nathan, 2026-10-15: "the AA length only goes in parentheses for the view that doesn't
    already have the AA". The structure rows already carried <id>|<aa> and keep it."""
    m = {"isoform_id": "52751007.WF86", "protein_length": 635}
    assert isoform_label(m) == "52751007.WF86 (635aa)"
    assert isoform_label(m, style="bar") == "52751007.WF86|635"


def test_nmd_marking_differs_by_style():
    m = {"isoform_id": "X", "protein_length": 39, "nmd_status": "Potential-NMD"}
    assert isoform_label(m) == "X (39aa, NMD)"
    assert isoform_label(m, style="bar") == "X|39 (NMD)"


def test_the_id_is_shortened_but_the_length_survives():
    """Truncating the whole label hid the very thing it was added for."""
    lab = isoform_label({"isoform_id": "52751007.WF86-MDS-0609-6020_CD34",
                         "protein_length": 635}, id_max=26)
    assert lab.endswith(" (635aa)")
    assert lab.startswith("52751007.WF86-MDS-06")
    assert "…" in lab


def test_a_short_id_is_left_alone():
    assert isoform_label({"isoform_id": "X", "protein_length": 82}, id_max=26) == "X (82aa)"
