#!/usr/bin/env python3
"""Consensus-ensemble binary analysis on k-fold OOF predictions.

Rule: an isoform is called PM only if >= K ensemble members predict PM
(prob >= 0.5). Otherwise it is flagged as NON-PM (excluded as a likely
true negative). Scans K from 1..N and reports per-group outcomes.

Reads:
    classical-kfold predictions (tm_isoform_kfold_predictions.tsv)
    residual-mlp kfold predictions (tm_isoform_kfold_predictions.tsv)

Writes:
    phase_d/checkpoints/ensemble_consensus_binary.tsv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

_DAEDALUS = Path(__file__).resolve().parents[2]
CLASSICAL_DIR = _DAEDALUS / "training" / "checkpoints" / "kfold_classical"
RESIDUAL_DIR = _DAEDALUS / "training" / "checkpoints" / "kfold_residual_mlp"

DEFAULT_MODELS = [
    # Top-performing per the single-threshold 0.5 binary analysis
    "prob_phase_d_logistic",
    "prob_phase_c_adaboost",
    "prob_residual_mlp_ae_strongrecon",
    "prob_phase_c_elasticnet",
    "prob_residual_mlp_ae_narrow",
    "prob_residual_mlp_ae_wide",
]


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--threshold", type=float, default=0.5,
                   help="Per-model threshold for calling PM")
    p.add_argument("--models", default=",".join(DEFAULT_MODELS),
                   help="Comma-separated prob_<model> columns to pool")
    p.add_argument("--out", type=Path,
                   default=Path("phase_d/checkpoints/ensemble_consensus_binary.tsv"))
    args = p.parse_args()

    classical = pd.read_csv(CLASSICAL_DIR / "tm_isoform_kfold_predictions.tsv", sep="\t")
    residual = pd.read_csv(RESIDUAL_DIR / "tm_isoform_kfold_predictions.tsv", sep="\t")

    join_cols = ["species", "gene_name", "primary_accession",
                 "reference_isoform_id", "alternative_isoform_id",
                 "benchmark_group", "label_binary", "alt_isoform_locations"]
    r_probs = [c for c in residual.columns if c.startswith("prob_")]
    merged = classical.merge(residual[join_cols + r_probs], on=join_cols, how="outer")

    wanted = [m.strip() for m in args.models.split(",") if m.strip()]
    for col in wanted:
        if col not in merged.columns:
            raise SystemExit(f"Missing column {col}. Available: {[c for c in merged.columns if c.startswith('prob_')]}")

    # Valid rows = predictions exist for every selected model
    prob_mat = merged[wanted].to_numpy(dtype=np.float64)
    valid = ~np.isnan(prob_mat).any(axis=1)
    prob_mat = prob_mat[valid]
    sub = merged.loc[valid].reset_index(drop=True)
    votes_pm = (prob_mat >= args.threshold).sum(axis=1)
    n_models = prob_mat.shape[1]

    neg_mask = sub["benchmark_group"] == "tm_retained_non_surface_true_negative"
    pos_mask = sub["benchmark_group"] == "uniprot_tm_positive"
    ctrl_mask = sub["benchmark_group"] == "tm_lost_no_tm_control"

    rows = []
    for k in range(1, n_models + 1):
        pred_pm = (votes_pm >= k).astype(int)   # 1 = called PM, 0 = excluded
        fp = int(((pred_pm == 1) & neg_mask).sum())
        tn = int(((pred_pm == 0) & neg_mask).sum())
        tp = int(((pred_pm == 1) & pos_mask).sum())
        fn = int(((pred_pm == 0) & pos_mask).sum())
        ctrl_fp = int(((pred_pm == 1) & ctrl_mask).sum())
        ctrl_tn = int(((pred_pm == 0) & ctrl_mask).sum())
        rows.append({
            "consensus_k": int(k),
            "rule": f">={k} of {n_models} models call PM",
            "tm_retained_negatives_n": int(neg_mask.sum()),
            "wrongly_called_PM": fp,
            "correctly_called_nonPM": tn,
            "pct_wrongly_called_PM": round(100 * fp / max(fp + tn, 1), 2),
            "positives_n": int(pos_mask.sum()),
            "pos_correctly_called_PM": tp,
            "pos_wrongly_excluded": fn,
            "pos_sensitivity_pct": round(100 * tp / max(tp + fn, 1), 2),
            "tm_lost_controls_n": int(ctrl_mask.sum()),
            "tm_lost_wrongly_called_PM": ctrl_fp,
            "tm_lost_correctly_called_nonPM": ctrl_tn,
            "tm_lost_pct_wrongly_called_PM": round(100 * ctrl_fp / max(ctrl_fp + ctrl_tn, 1), 2),
        })

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"Ensemble members ({n_models}): {wanted}")
    print(f"Wrote: {args.out}\n")
    for r in rows:
        print(f"  k>={r['consensus_k']}:  wrongly-PM (neg) {r['wrongly_called_PM']}/{r['tm_retained_negatives_n']} "
              f"({r['pct_wrongly_called_PM']}%)  |  sens {r['pos_sensitivity_pct']}%  "
              f"| ctrl wrong {r['tm_lost_pct_wrongly_called_PM']}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
