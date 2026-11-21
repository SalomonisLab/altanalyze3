"""Standard variant-evidence cutoffs: values of record and the rule functions (MDS-AML-KINNEX-5 cases)."""
import argparse
import sys

from altanalyze3.components.bam import variant_evidence_defaults as d


def test_values_of_record():
    assert d.NUCLEAR_PASSENGER_MAX_BACKGROUND == 0.15
    assert d.MT_PASSENGER_MAX_BACKGROUND == 0.15
    assert d.PASSENGER_MAX_Q == 0.05 and d.PASSENGER_MIN_MUT_CELLS == 10
    assert d.CLASSIFIER_MIN_CELLS_PER_CLASS == 15 and d.CLASSIFIER_MIN_ACCURACY == 0.85


def test_germline_rule_on_known_cases():
    assert d.germline_call(1210, 1211, 34, 34) == "germline_like"        # TET2 c.86C>G, AM75
    assert d.germline_call(5343, 5350, 4, 4) == "undetermined"           # GATA1 c.212A>G, AM75: 4 lymphoid reads
    assert d.germline_call(580, 609, 1, 1) == "undetermined"             # BCOR c.2257_2258insG
    assert d.germline_call(1690, 1821, 0, 2) == "somatic_pattern"        # U2AF1-S34F: patient VAF 0.93
    assert d.germline_call(99, 99, 30, 30) == "somatic_pattern"          # fewer than 100 patient reads
    assert d.germline_call(1000, 1000, 10, 40) == "somatic_pattern"      # lymphoid cells hold reference reads


def test_self_marker():
    assert d.is_self_marker("NPM1_chr5:171410542", "NPM1_chr5:171410539_C>CTCTG")
    assert not d.is_self_marker("NPM1_chr5:171410542", "NPM1_chr5:171410729_CT>C")
    assert not d.is_self_marker("DDX10::SKA3_chr11:108715038", "PTMA_chr2:231713269_G>GA")
    assert not d.is_self_marker("NPM1_chr5:171410542", "chrM:9570_A>AC")


def test_loh_gate():
    assert d.loh_is_evidence(819, 57, 324, 872)          # 5801 chr21, RUNX1 p.W279*
    assert not d.loh_is_evidence(7, 22, 22, 106)         # AM75 chr4, TET2 c.3646C>T


def test_module_cli_defaults_follow_the_cutoff_file(monkeypatch):
    from altanalyze3.components.bam import mt_variants
    monkeypatch.setattr(sys, "argv", ["mt_variants", "--output-dir", "x", "--nominate-fdr"])
    assert mt_variants.parse_arguments().max_background == d.MT_PASSENGER_MAX_BACKGROUND
    monkeypatch.setattr(sys, "argv", ["mt_variants", "--output-dir", "x"])
    assert mt_variants.parse_arguments().max_background == 0.10


def test_inherited_by_clone_absent_sample():
    others = [(6, 11195), (0, 366), (0, 496), (2, 111)]                      # AM75 on therapy: SRSF2, EZH2, RUNX1, TET2 c.3646
    assert d.inherited_by_clone_absent_sample(1611, 1611, others)            # GATA1 c.212A>G
    assert d.inherited_by_clone_absent_sample(333, 333, others)              # TET2 c.86C>G
    assert not d.inherited_by_clone_absent_sample(1289, 1291, [(4373, 9127)])    # clone present: no call
    assert not d.inherited_by_clone_absent_sample(100, 100, [(0, 10)])           # no other variant with enough reads
    assert not d.inherited_by_clone_absent_sample(40, 100, others)               # variant itself not near-fixed
