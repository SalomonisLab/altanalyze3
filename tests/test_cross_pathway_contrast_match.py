"""A comparison named in a question must route to that contrast.

The stored label is `cancer_vs_no_cancer`; a reader types "in cancer vs no cancer".
The verbatim check answered "Select a saved comparison", so the contrast view, which
carries the imputed modalities, never ran (COPD viewer, 2026-09-28).
"""
import pytest

from altanalyze3.components.cellHarmony.webapp.cross_pathways import comparison_matches, _label_tokens

LABELS = ['COPD_vs_non-COPD', 'DLCO_high_vs_DLCO_low', 'GOLD_IV_vs_GOLD_I_II', 'M_vs_F',
          'cancer_vs_no_cancer', 'current_vs_never', 'former_vs_never']


def hits(question):
    return [label for label in LABELS if comparison_matches(question, label)]


def test_reader_phrasing_matches_the_underscored_label():
    assert hits("best cross-modality AT2 representation in cancer vs no cancer?") == ['cancer_vs_no_cancer']


def test_versus_and_hyphens_and_case_are_equivalent():
    assert hits("AT1 pathways in COPD versus non-COPD") == ['COPD_vs_non-COPD']
    assert hits("AT1 in gold iv vs gold i ii") == ['GOLD_IV_vs_GOLD_I_II']


def test_a_shorter_label_does_not_shadow_a_longer_one():
    assert hits("AT1 pathways former vs never") == ['former_vs_never']
    assert hits("AT1 pathways current vs never") == ['current_vs_never']


def test_verbatim_label_still_matches():
    assert hits("AT2 in cancer_vs_no_cancer") == ['cancer_vs_no_cancer']


def test_no_label_means_no_match():
    assert hits("What pathways have the best cross-modality AT2 representation?") == []
    assert not comparison_matches("anything", "")


def test_tokens_split_on_underscore_hyphen_and_map_versus():
    assert _label_tokens("COPD_vs_non-COPD") == ['copd', 'vs', 'non', 'copd']
    assert _label_tokens("cancer versus no cancer") == ['cancer', 'vs', 'no', 'cancer']
