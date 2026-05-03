from __future__ import annotations

import math
from typing import Any, Mapping


def _f(row: Mapping[str, Any], key: str, default: float = 0.0) -> float:
    value = row.get(key, default)
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _clamp01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))


def _ratio(numer: float, denom: float, default: float = 0.0) -> float:
    if abs(denom) < 1e-12:
        return default
    return numer / denom


def _similarity_delta(delta: float, scale: float) -> float:
    return _clamp01(1.0 - abs(delta) / max(scale, 1e-6))


def _max2(a: float, b: float) -> float:
    return a if a >= b else b


def _signal_support(row: Mapping[str, Any]) -> float:
    alt_signal = _max2(
        _f(row, "alternative_predicted_signal_candidate"),
        _clamp01(_f(row, "alternative_predicted_signal_score") / 3.0),
    )
    ref_signal = _max2(
        _f(row, "reference_has_signal_peptide"),
        _max2(_clamp01(_f(row, "reference_signal_feature_count")), _f(row, "reference_predicted_signal_candidate")),
    )
    return _clamp01(0.65 * alt_signal + 0.35 * ref_signal)


def _tm_similarity(row: Mapping[str, Any]) -> float:
    ref_count = _f(row, "reference_predicted_tm_count")
    ref_span = _f(row, "reference_predicted_tm_total_span")
    count_sim = _similarity_delta(_f(row, "predicted_tm_count_delta"), max(ref_count, 1.0))
    span_sim = _similarity_delta(_f(row, "predicted_tm_total_span_delta"), max(ref_span, 20.0))
    hydropathy = _clamp01(_f(row, "alternative_predicted_tm_max_hydropathy_19") / 2.5)
    return _clamp01(0.45 * count_sim + 0.35 * span_sim + 0.20 * hydropathy)


def _tm_disruption(row: Mapping[str, Any]) -> float:
    return _clamp01(
        0.45 * _f(row, "overlap_tm_frac")
        + 0.20 * _f(row, "is_complex_change")
        + 0.20 * _f(row, "is_truncation")
        + 0.15 * _clamp01(_f(row, "alt_changed_frac"))
    )


def score_transmembrane_tasks(row: Mapping[str, Any]) -> dict[str, float]:
    signal = _signal_support(row)
    nterm_change_penalty = _clamp01(_f(row, "is_n_term_change") * _f(row, "alt_changed_frac"))
    truncation_penalty = _clamp01(_f(row, "is_truncation") * (1.0 - _f(row, "protein_length_ratio", 1.0)))
    tm_similarity = _tm_similarity(row)
    tm_disruption = _tm_disruption(row)
    extracellular_prior = _clamp01(
        0.55 * _clamp01(_ratio(_f(row, "reference_extracellular_topology_aa"), 80.0))
        + 0.20 * _f(row, "gene_hpa_has_extracellular_annotation")
        + 0.25 * _clamp01(_f(row, "reference_tm_feature_count"))
    )
    extracellular_access = _clamp01(
        0.45 * _f(row, "overlap_extracellular_frac")
        + 0.20 * extracellular_prior
        + 0.15 * (1.0 - _f(row, "overlap_tm_frac"))
        + 0.10 * (1.0 - _f(row, "overlap_cytoplasmic_frac"))
        + 0.10 * _clamp01(_ratio(_f(row, "alt_changed_len"), 30.0))
    )
    glyco_support = _clamp01(
        0.60 * _clamp01(_ratio(_f(row, "alternative_glyco_motif_count"), 3.0))
        + 0.40 * _clamp01(_ratio(_f(row, "reference_glycosylation_count"), 3.0))
    )
    disulfide_support = _clamp01(_ratio(_f(row, "reference_disulfide_count"), 4.0))
    cysteine_support = _clamp01(_ratio(_f(row, "alternative_cysteine_count"), 8.0))
    secretory_pathway_compatibility = _clamp01(
        0.45 * signal
        + 0.20 * (1.0 - nterm_change_penalty)
        + 0.15 * (1.0 - _f(row, "is_complex_change"))
        + 0.20 * (1.0 - truncation_penalty)
    )
    membrane_insertion_topology = _clamp01(
        0.35 * tm_similarity
        + 0.20 * secretory_pathway_compatibility
        + 0.20 * (1.0 - tm_disruption)
        + 0.15 * (1.0 - _f(row, "overlap_cytoplasmic_frac"))
        + 0.10 * _clamp01(_ratio(_f(row, "reference_cytoplasmic_topology_aa"), 80.0))
    )
    folding_stability_qc_escape = _clamp01(
        0.30 * membrane_insertion_topology
        + 0.20 * secretory_pathway_compatibility
        + 0.15 * glyco_support
        + 0.15 * disulfide_support
        + 0.10 * cysteine_support
        + 0.10 * (1.0 - truncation_penalty)
    )
    cell_surface_localization = _clamp01(
        0.25 * secretory_pathway_compatibility
        + 0.25 * membrane_insertion_topology
        + 0.25 * folding_stability_qc_escape
        + 0.25 * extracellular_prior
    )
    antibody_targetability = _clamp01(
        0.45 * cell_surface_localization
        + 0.35 * extracellular_access
        + 0.10 * _clamp01(_ratio(_f(row, "alt_changed_len"), 25.0))
        + 0.10 * (1.0 - _f(row, "overlap_cytoplasmic_frac"))
    )
    stable_functional = _clamp01(
        0.18 * secretory_pathway_compatibility
        + 0.22 * membrane_insertion_topology
        + 0.15 * extracellular_access
        + 0.20 * folding_stability_qc_escape
        + 0.25 * cell_surface_localization
    )
    non_stable_non_functional = _clamp01(
        1.0
        - (
            0.20 * secretory_pathway_compatibility
            + 0.25 * membrane_insertion_topology
            + 0.25 * folding_stability_qc_escape
            + 0.30 * cell_surface_localization
        )
    )
    return {
        "tm_secretory_pathway_compatibility_score": secretory_pathway_compatibility,
        "tm_membrane_insertion_topology_score": membrane_insertion_topology,
        "tm_extracellular_altered_segment_accessibility_score": extracellular_access,
        "tm_folding_stability_qc_escape_score": folding_stability_qc_escape,
        "tm_cell_surface_localization_score": cell_surface_localization,
        "tm_antibody_targetability_score": antibody_targetability,
        "tm_stable_functional_score": stable_functional,
        "tm_non_stable_non_functional_score": non_stable_non_functional,
    }


def _nls_support(row: Mapping[str, Any]) -> float:
    ref_nls = _f(row, "ref_full_motif_nls_basic_count")
    alt_nls = _f(row, "alt_full_motif_nls_basic_count")
    retained = _clamp01(_ratio(alt_nls, max(ref_nls, 1.0), 1.0 if alt_nls > 0 else 0.0))
    return _clamp01(
        0.55 * retained
        + 0.25 * (1.0 - _f(row, "overlap_nls_region_frac"))
        + 0.20 * (1.0 - _f(row, "overlap_nes_region_frac"))
    )


def score_tf_tasks(row: Mapping[str, Any]) -> dict[str, float]:
    truncation_penalty = _clamp01(_f(row, "is_truncation") * (1.0 - _f(row, "protein_length_ratio", 1.0)))
    nuclear_localization = _clamp01(
        0.45 * _nls_support(row)
        + 0.20 * (1.0 - _f(row, "is_n_term_change") * _f(row, "alt_changed_frac"))
        + 0.15 * (1.0 - _f(row, "overlap_nes_region_frac"))
        + 0.20 * (1.0 - truncation_penalty)
    )
    dna_ref_weight = _clamp01(
        _ratio(
            _f(row, "reference_dna_binding_feature_count")
            + _f(row, "reference_dna_region_feature_count")
            + _f(row, "reference_zinc_finger_feature_count"),
            6.0,
        )
    )
    dbd_integrity_specificity = _clamp01(
        0.45 * (1.0 - _f(row, "overlap_dna_interface_frac"))
        + 0.20 * dna_ref_weight
        + 0.15 * (1.0 - truncation_penalty)
        + 0.20 * (1.0 - _f(row, "is_complex_change"))
    )
    cofactor_ppi_rewiring_risk = _clamp01(
        0.35 * _f(row, "overlap_ppi_interface_frac")
        + 0.20 * _f(row, "overlap_dimerization_frac")
        + 0.15 * _f(row, "overlap_activation_region_frac")
        + 0.10 * _f(row, "overlap_repression_region_frac")
        + 0.20 * truncation_penalty
    )
    activation_repression_competence = _clamp01(
        0.30 * nuclear_localization
        + 0.30 * dbd_integrity_specificity
        + 0.20 * (1.0 - _f(row, "overlap_activation_region_frac"))
        + 0.20 * (1.0 - _f(row, "overlap_repression_region_frac"))
    )
    dominant_negative_neomorphic_risk = _clamp01(
        0.30 * nuclear_localization
        + 0.20 * (1.0 - _f(row, "overlap_dimerization_frac"))
        + 0.20 * (1.0 - _f(row, "overlap_dna_interface_frac"))
        + 0.30 * (1.0 - activation_repression_competence)
    )
    overall_regulatory_functionality = _clamp01(
        0.25 * nuclear_localization
        + 0.25 * dbd_integrity_specificity
        + 0.15 * (1.0 - cofactor_ppi_rewiring_risk)
        + 0.20 * activation_repression_competence
        + 0.15 * (1.0 - dominant_negative_neomorphic_risk)
    )
    non_stable_non_functional = _clamp01(
        1.0
        - (
            0.22 * nuclear_localization
            + 0.22 * dbd_integrity_specificity
            + 0.16 * (1.0 - cofactor_ppi_rewiring_risk)
            + 0.20 * activation_repression_competence
            + 0.20 * (1.0 - dominant_negative_neomorphic_risk)
        )
    )
    return {
        "tf_nuclear_localization_score": nuclear_localization,
        "tf_dbd_integrity_specificity_score": dbd_integrity_specificity,
        "tf_cofactor_ppi_rewiring_risk_score": cofactor_ppi_rewiring_risk,
        "tf_activation_repression_competence_score": activation_repression_competence,
        "tf_dominant_negative_neomorphic_risk_score": dominant_negative_neomorphic_risk,
        "tf_overall_regulatory_functionality_score": overall_regulatory_functionality,
        "tf_non_stable_non_functional_score": non_stable_non_functional,
    }


def score_phase_d_objectives(row: Mapping[str, Any]) -> dict[str, float | None]:
    family = str(row.get("family_class", ""))
    scores: dict[str, float | None] = {}
    if family == "membrane":
        scores.update(score_transmembrane_tasks(row))
    else:
        scores.update(
            {
                "tm_secretory_pathway_compatibility_score": None,
                "tm_membrane_insertion_topology_score": None,
                "tm_extracellular_altered_segment_accessibility_score": None,
                "tm_folding_stability_qc_escape_score": None,
                "tm_cell_surface_localization_score": None,
                "tm_antibody_targetability_score": None,
                "tm_stable_functional_score": None,
                "tm_non_stable_non_functional_score": None,
            }
        )
    if family == "transcription_factor":
        scores.update(score_tf_tasks(row))
    else:
        scores.update(
            {
                "tf_nuclear_localization_score": None,
                "tf_dbd_integrity_specificity_score": None,
                "tf_cofactor_ppi_rewiring_risk_score": None,
                "tf_activation_repression_competence_score": None,
                "tf_dominant_negative_neomorphic_risk_score": None,
                "tf_overall_regulatory_functionality_score": None,
                "tf_non_stable_non_functional_score": None,
            }
        )
    if family == "membrane":
        stable = float(scores["tm_stable_functional_score"])
        unstable = float(scores["tm_non_stable_non_functional_score"])
    elif family == "transcription_factor":
        stable = float(scores["tf_overall_regulatory_functionality_score"])
        unstable = float(scores["tf_non_stable_non_functional_score"])
    else:
        stable = None
        unstable = None
    if stable is None or unstable is None:
        scores["aggregate_stable_functional_score"] = None
        scores["aggregate_non_stable_non_functional_score"] = None
        scores["aggregate_functional_annotation"] = None
    else:
        scores["aggregate_stable_functional_score"] = stable
        scores["aggregate_non_stable_non_functional_score"] = unstable
        if stable >= 0.70 and unstable <= 0.35:
            label = "stable_functional"
        elif unstable >= 0.70 and stable <= 0.35:
            label = "non_stable_non_functional"
        else:
            label = "indeterminate"
        scores["aggregate_functional_annotation"] = label
    return scores
