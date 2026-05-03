#!/usr/bin/env python3
"""Summarize out-of-fold (OOF) k-fold predictions: at binary threshold 0.5,
how many true-negative TM isoforms are wrongly predicted as Cell Membrane?

Reads:
    {kfold_dir}/tm_isoform_kfold_predictions.tsv  (one row per benchmark pair,
                                                   columns `prob_<model>` per model)

Emits:
    {kfold_dir}/tm_isoform_kfold_binary_negative_detection.tsv
    (one row per model; counts and rates on the tm_retained true-negative set)

Usage:
    python summarize_kfold_negative_detection.py \
        --kfold-dir phase_d/checkpoints/kfold_classical \
        --threshold 0.5
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

NEG_GROUP_TM_RETAINED = "tm_retained_non_surface_true_negative"
NEG_GROUP_TM_LOST = "tm_lost_no_tm_control"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--kfold-dir", type=Path, required=True)
    p.add_argument("--threshold", type=float, default=0.5)
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args()

    in_tsv = args.kfold_dir / "tm_isoform_kfold_predictions.tsv"
    if not in_tsv.exists():
        sys.exit(f"missing: {in_tsv}")
    df = pd.read_csv(in_tsv, sep="\t")
    model_cols = [c for c in df.columns if c.startswith("prob_")]
    if not model_cols:
        sys.exit("no prob_<model> columns found")

    # tm_retained true-negative set
    neg_mask = df["benchmark_group"] == NEG_GROUP_TM_RETAINED
    pos_mask = df["benchmark_group"] == "uniprot_tm_positive"
    ctrl_mask = df["benchmark_group"] == NEG_GROUP_TM_LOST
    rows = []
    for col in model_cols:
        model = col.replace("prob_", "", 1)
        prob = df[col].to_numpy(dtype=np.float64)
        valid = ~np.isnan(prob)
        # Per-scope counts
        def _counts(mask):
            m = mask & valid
            probs = prob[m]
            labels = df.loc[m, "label_binary"].to_numpy(dtype=np.int32)
            pred = (probs >= args.threshold).astype(int)
            tp = int(((pred == 1) & (labels == 1)).sum())
            fn = int(((pred == 0) & (labels == 1)).sum())
            tn = int(((pred == 0) & (labels == 0)).sum())
            fp = int(((pred == 1) & (labels == 0)).sum())
            return tp, fn, tn, fp, int(m.sum())

        tp_n, fn_n, tn_n, fp_n, n_neg = _counts(neg_mask)
        tp_p, fn_p, tn_p, fp_p, n_pos = _counts(pos_mask)
        tp_c, fn_c, tn_c, fp_c, n_ctrl = _counts(ctrl_mask)

        total_neg = fp_n + tn_n  # all are label=0
        miss_neg_rate = (fp_n / total_neg) if total_neg else float("nan")

        rows.append({
            "model": model,
            "threshold": float(args.threshold),
            "tm_retained_n_negatives": total_neg,
            "tm_retained_wrongly_called_PM": fp_n,
            "tm_retained_correctly_called_nonPM": tn_n,
            "tm_retained_pct_wrongly_called_PM": round(100 * miss_neg_rate, 2),
            "tm_lost_n": (fp_c + tn_c),
            "tm_lost_wrongly_called_PM": fp_c,
            "tm_lost_pct_wrongly_called_PM": round(
                100 * fp_c / max(fp_c + tn_c, 1), 2
            ),
            "pos_n": (tp_p + fn_p),
            "pos_correctly_called_PM": tp_p,
            "pos_pct_correctly_called_PM": round(
                100 * tp_p / max(tp_p + fn_p, 1), 2
            ),
        })

    out = args.out or (args.kfold_dir / "tm_isoform_kfold_binary_negative_detection.tsv")
    pd.DataFrame(rows).sort_values(
        "tm_retained_pct_wrongly_called_PM"
    ).to_csv(out, sep="\t", index=False)
    print(f"Wrote: {out}")
    for r in rows:
        print(f"  {r['model']:<40s}  wrongly-PM={r['tm_retained_wrongly_called_PM']}/{r['tm_retained_n_negatives']}  ({r['tm_retained_pct_wrongly_called_PM']}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
