#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PHASE_A = Path(__file__).resolve().parents[2] / "phase_a" / "data" / "interim"
EXTERNAL = Path(__file__).resolve().parents[1] / "data" / "external"
IN_PATH = PHASE_A / "isoform_pair_candidates.with_splits.tsv"
OUT_PATH = EXTERNAL / "proxy_supervised_tasks.tsv"


def _ratio(num: pd.Series, den: pd.Series) -> pd.Series:
    den_safe = den.replace(0, np.nan)
    return (num / den_safe).fillna(0.0)


def _append_task(rows: list[pd.DataFrame], base: pd.DataFrame, task: str, **extra: object) -> None:
    frame = base.copy()
    frame["task"] = task
    for key, value in extra.items():
        frame[key] = value
    rows.append(frame)


def main() -> int:
    EXTERNAL.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(IN_PATH, sep="\t")
    id_cols = ["species", "gene_name", "reference_transcript_id", "alternative_transcript_id"]
    rows: list[pd.DataFrame] = []

    ppi = df[df["reference_ppi_binding_feature_count"] > 0][id_cols + ["reference_ppi_binding_feature_count", "alternative_ppi_binding_feature_count", "split"]].copy()
    if not ppi.empty:
        _append_task(
            rows,
            ppi[id_cols + ["split"]],
            "ppi_impact",
            ppi_preserved=(ppi["alternative_ppi_binding_feature_count"] > 0).astype(int).to_numpy(),
            ppi_delta_score=_ratio(ppi["alternative_ppi_binding_feature_count"], ppi["reference_ppi_binding_feature_count"]).to_numpy(),
            evidence_strength=np.log1p(ppi["reference_ppi_binding_feature_count"]).to_numpy(),
        )

    tf_mask = df["family_class"] == "transcription_factor"
    tf_ref_support = (
        df["reference_dna_binding_feature_count"]
        + df["reference_dna_region_feature_count"]
        + df["reference_zinc_finger_feature_count"]
    )
    tf_alt_support = (
        df["alternative_dna_binding_feature_count"]
        + df["alternative_dna_region_feature_count"]
        + df["alternative_zinc_finger_feature_count"]
    )
    tf_pairs = df[tf_mask & (tf_ref_support > 0)].copy()
    if not tf_pairs.empty:
        ref_support = (
            tf_pairs["reference_dna_binding_feature_count"]
            + tf_pairs["reference_dna_region_feature_count"]
            + tf_pairs["reference_zinc_finger_feature_count"]
        )
        alt_support = (
            tf_pairs["alternative_dna_binding_feature_count"]
            + tf_pairs["alternative_dna_region_feature_count"]
            + tf_pairs["alternative_zinc_finger_feature_count"]
        )
        _append_task(
            rows,
            tf_pairs[id_cols + ["split"]],
            "pdi_impact",
            pdi_preserved=(alt_support > 0).astype(int).to_numpy(),
            pdi_delta_score=_ratio(alt_support, ref_support).to_numpy(),
            evidence_strength=np.log1p(ref_support).to_numpy(),
        )
        _append_task(
            rows,
            tf_pairs[id_cols + ["split"]],
            "tf_dna_binding",
            pdi_preserved=(alt_support > 0).astype(int).to_numpy(),
            pdi_delta_score=_ratio(alt_support, ref_support).to_numpy(),
            evidence_strength=np.log1p(ref_support).to_numpy(),
        )
        ref_activity = tf_pairs["reference_domain_count"] + tf_pairs["reference_motif_count"] + ref_support
        alt_activity = tf_pairs["alternative_domain_count"] + tf_pairs["alternative_motif_count"] + alt_support
        _append_task(
            rows,
            tf_pairs[id_cols + ["split"]],
            "tf_transcriptional_activity",
            activity_preserved=(alt_activity > 0).astype(int).to_numpy(),
            transcriptional_activity_score=_ratio(alt_activity, ref_activity).to_numpy(),
            evidence_strength=np.log1p(ref_activity).to_numpy(),
        )
        _append_task(
            rows,
            tf_pairs[id_cols + ["split"]],
            "tf_perturbation_response",
            bootstrap_significant=(alt_activity > 0).astype(int).to_numpy(),
            perturbation_effect_size=_ratio(alt_activity, ref_activity).to_numpy(),
            evidence_strength=np.log1p(ref_activity).to_numpy(),
        )

    loc_ref = (
        (df["reference_is_membrane"] > 0)
        | (df["reference_has_signal_peptide"] > 0)
        | (df["reference_extracellular_topology_aa"] > 0)
    )
    loc_pairs = df[loc_ref].copy()
    if not loc_pairs.empty:
        alt_loc = (
            (loc_pairs["alternative_is_membrane"] > 0)
            | (loc_pairs["alternative_has_signal_peptide"] > 0)
            | (loc_pairs["alternative_extracellular_topology_aa"] > 0)
        )
        loc_class = np.where(
            alt_loc & (loc_pairs["alternative_extracellular_topology_aa"] > 0),
            "surface_retained",
            np.where(alt_loc, "membrane_like", "surface_lost"),
        )
        _append_task(
            rows,
            loc_pairs[id_cols + ["split"]],
            "localization_shift",
            localization_preserved=alt_loc.astype(int).to_numpy(),
            localization_class=loc_class,
            evidence_strength=np.log1p(
                loc_pairs["reference_extracellular_topology_aa"] + loc_pairs["reference_tm_feature_count"]
            ).to_numpy(),
        )

    kinase = df[(df["family_class"] == "kinase") & ((df["reference_phospho_feature_count"] > 0) | (df["reference_domain_count"] > 0))].copy()
    if not kinase.empty:
        kinase_ref = kinase["reference_phospho_feature_count"] + kinase["reference_domain_count"]
        kinase_alt = kinase["alternative_phospho_feature_count"] + kinase["alternative_domain_count"]
        _append_task(
            rows,
            kinase[id_cols + ["split"]],
            "kinase_signaling",
            kinase_activity_preserved=((kinase["alternative_phospho_feature_count"] > 0) & (kinase["alternative_domain_count"] > 0)).astype(int).to_numpy(),
            kinase_activity_score=_ratio(kinase_alt, kinase_ref).to_numpy(),
            evidence_strength=np.log1p(kinase_ref).to_numpy(),
        )

    signal_pairs = df[(df["reference_has_signal_peptide"] > 0) | (df["reference_signal_feature_count"] > 0)].copy()
    if not signal_pairs.empty:
        _append_task(
            rows,
            signal_pairs[id_cols + ["split"]],
            "signal_retention",
            signal_retained=((signal_pairs["alternative_has_signal_peptide"] > 0) | (signal_pairs["alternative_signal_feature_count"] > 0)).astype(int).to_numpy(),
            evidence_strength=np.log1p(signal_pairs["reference_signal_feature_count"] + signal_pairs["reference_has_signal_peptide"]).to_numpy(),
        )

    tm_pairs = df[df["reference_tm_feature_count"] > 0].copy()
    if not tm_pairs.empty:
        _append_task(
            rows,
            tm_pairs[id_cols + ["split"]],
            "tm_insertion",
            tm_inserted=(tm_pairs["alternative_tm_feature_count"] > 0).astype(int).to_numpy(),
            tm_insertion_score=_ratio(tm_pairs["alternative_tm_feature_count"], tm_pairs["reference_tm_feature_count"]).to_numpy(),
            evidence_strength=np.log1p(tm_pairs["reference_tm_feature_count"]).to_numpy(),
        )
        tm_folded = (
            (tm_pairs["alternative_tm_feature_count"] > 0)
            & (
                (tm_pairs["alternative_extracellular_topology_aa"] > 0)
                | (tm_pairs["alternative_glycosylation_count"] > 0)
                | (tm_pairs["alternative_disulfide_count"] > 0)
            )
        )
        fold_score = (
            _ratio(tm_pairs["alternative_extracellular_topology_aa"], tm_pairs["reference_extracellular_topology_aa"].replace(0, 1))
            + _ratio(tm_pairs["alternative_glycosylation_count"], tm_pairs["reference_glycosylation_count"].replace(0, 1))
            + _ratio(tm_pairs["alternative_disulfide_count"], tm_pairs["reference_disulfide_count"].replace(0, 1))
        ) / 3.0
        _append_task(
            rows,
            tm_pairs[id_cols + ["split"]],
            "tm_fold",
            tm_folded=tm_folded.astype(int).to_numpy(),
            tm_fold_score=fold_score.to_numpy(),
            evidence_strength=np.log1p(
                tm_pairs["reference_tm_feature_count"]
                + tm_pairs["reference_extracellular_topology_aa"]
                + tm_pairs["reference_glycosylation_count"]
                + tm_pairs["reference_disulfide_count"]
            ).to_numpy(),
        )

    out = pd.concat(rows, axis=0, ignore_index=True)
    out.to_csv(OUT_PATH, sep="\t", index=False)
    print(f"[ok] Wrote {OUT_PATH} rows={len(out)}")
    for task, n in out["task"].value_counts().sort_index().items():
        print(f"[ok] task={task} rows={n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
