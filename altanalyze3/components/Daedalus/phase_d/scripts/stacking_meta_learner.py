#!/usr/bin/env python3
"""Logistic-regression stacking meta-learner on OOF probability vectors.

Joins classical and residual-MLP OOF predictions, uses all prob_<model>
columns as meta-features, trains a gene-grouped logistic regression via
inner 5-fold CV to produce meta-OOF predictions, then sweeps thresholds
to report the NPV / sensitivity tradeoff — compared against the K=6
unanimous consensus baseline.

Reads:
    phase_d/checkpoints/kfold_classical/tm_isoform_kfold_predictions.tsv
    phase_d/checkpoints/kfold_residual_mlp/tm_isoform_kfold_predictions.tsv

Writes:
    phase_d/checkpoints/stacking/meta_oof_predictions.tsv
    phase_d/checkpoints/stacking/meta_threshold_sweep.tsv
    phase_d/checkpoints/stacking/meta_coefficients.tsv
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

CLASSICAL_DIR = Path("phase_d/checkpoints/uniprot_kfold_classical_tier1")
RESIDUAL_DIR = Path("phase_d/checkpoints/uniprot_kfold_residual_tier1")
OUT_DIR = Path("phase_d/checkpoints/uniprot_stacking_tier1")

JOIN_COLS = [
    "species", "gene_name", "primary_accession",
    "reference_isoform_id", "alternative_isoform_id",
    "benchmark_group", "label_binary", "alt_isoform_locations",
]

NEG_GROUP = "tm_retained_non_surface_true_negative"
POS_GROUP = "uniprot_tm_positive"
CTRL_GROUP = "tm_lost_no_tm_control"


def _gene_group(row: pd.Series) -> str:
    key = f"{row['species']}|{row['gene_name']}"
    return hashlib.md5(key.encode()).hexdigest()


def _threshold_metrics(sub: pd.DataFrame, meta_prob: np.ndarray, threshold: float) -> dict:
    pred = (meta_prob >= threshold).astype(int)
    neg = sub["benchmark_group"] == NEG_GROUP
    pos = sub["benchmark_group"] == POS_GROUP
    ctrl = sub["benchmark_group"] == CTRL_GROUP

    fp = int(((pred == 1) & neg).sum())
    tn = int(((pred == 0) & neg).sum())
    tp = int(((pred == 1) & pos).sum())
    fn = int(((pred == 0) & pos).sum())
    ctrl_fp = int(((pred == 1) & ctrl).sum())
    ctrl_tn = int(((pred == 0) & ctrl).sum())

    sensitivity = tp / max(tp + fn, 1)
    specificity = tn / max(tn + fp, 1)
    npv = tn / max(tn + fn, 1)
    ppv = tp / max(tp + fp, 1)

    return {
        "threshold": round(threshold, 4),
        "neg_n": int(neg.sum()),
        "wrongly_called_PM": fp,
        "correctly_called_nonPM": tn,
        "pct_wrongly_called_PM": round(100 * fp / max(fp + tn, 1), 2),
        "pos_n": int(pos.sum()),
        "pos_correctly_called_PM": tp,
        "pos_wrongly_excluded": fn,
        "sensitivity_pct": round(100 * sensitivity, 2),
        "specificity_pct": round(100 * specificity, 2),
        "npv": round(npv, 4),
        "ppv": round(ppv, 4),
        "ctrl_n": int(ctrl.sum()),
        "ctrl_wrongly_called_PM": ctrl_fp,
        "ctrl_correctly_called_nonPM": ctrl_tn,
        "ctrl_pct_wrongly_called_PM": round(100 * ctrl_fp / max(ctrl_fp + ctrl_tn, 1), 2),
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--classical-dir", type=Path, default=CLASSICAL_DIR)
    p.add_argument("--residual-dir", type=Path, default=RESIDUAL_DIR)
    p.add_argument("--out-dir", type=Path, default=OUT_DIR)
    p.add_argument("--n-folds", type=int, default=5,
                   help="Inner gene-grouped CV folds for meta-OOF")
    p.add_argument("--c", type=float, default=1.0,
                   help="LogisticRegression C (regularisation inverse strength)")
    p.add_argument("--sensitivity-floor", type=float, default=0.85,
                   help="Minimum sensitivity when selecting operating threshold")
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    classical = pd.read_csv(args.classical_dir / "tm_isoform_kfold_predictions.tsv", sep="\t")
    residual = pd.read_csv(args.residual_dir / "tm_isoform_kfold_predictions.tsv", sep="\t")

    r_probs = [c for c in residual.columns if c.startswith("prob_")]
    merged = classical.merge(residual[JOIN_COLS + r_probs], on=JOIN_COLS, how="inner")

    feat_cols = [c for c in merged.columns if c.startswith("prob_")]
    print(f"Meta-features ({len(feat_cols)}): {feat_cols}")

    X_full = merged[feat_cols].to_numpy(dtype=np.float64)
    y_full = merged["label_binary"].to_numpy(dtype=np.int32)

    # Drop rows missing ANY meta-feature
    valid = ~np.isnan(X_full).any(axis=1)
    X = X_full[valid]
    y = y_full[valid]
    sub = merged.loc[valid].reset_index(drop=True)

    groups = sub.apply(_gene_group, axis=1).to_numpy()
    print(f"Valid rows: {valid.sum()} / {len(merged)}  |  "
          f"unique genes: {len(np.unique(groups))}")

    # Gene-grouped inner CV to produce meta-OOF probabilities
    meta_prob = np.full(len(sub), np.nan)
    fold_coefs = []

    gkf = GroupKFold(n_splits=args.n_folds)
    for fold_idx, (train_idx, val_idx) in enumerate(gkf.split(X, y, groups)):
        X_tr, X_va = X[train_idx], X[val_idx]
        y_tr = y[train_idx]

        scaler = StandardScaler()
        X_tr_s = scaler.fit_transform(X_tr)
        X_va_s = scaler.transform(X_va)

        clf = LogisticRegression(
            C=args.c,
            class_weight="balanced",
            max_iter=1000,
            solver="lbfgs",
        )
        clf.fit(X_tr_s, y_tr)
        meta_prob[val_idx] = clf.predict_proba(X_va_s)[:, 1]

        fold_val = sub.iloc[val_idx]
        neg_mask = fold_val["benchmark_group"] == NEG_GROUP
        pos_mask = fold_val["benchmark_group"] == POS_GROUP
        pv = meta_prob[val_idx]
        pred05 = (pv >= 0.5).astype(int)
        fp = int(((pred05 == 1) & neg_mask.values).sum())
        tn = int(((pred05 == 0) & neg_mask.values).sum())
        tp = int(((pred05 == 1) & pos_mask.values).sum())
        fn = int(((pred05 == 0) & pos_mask.values).sum())
        print(f"  fold {fold_idx}: neg_wrong={fp}/{fp+tn}  "
              f"sens={round(100*tp/max(tp+fn,1),1)}%  "
              f"n_train={len(train_idx)}  n_val={len(val_idx)}")
        fold_coefs.append(clf.coef_[0])

    # Save meta-OOF predictions
    out_oof = sub.copy()
    out_oof["meta_prob"] = meta_prob
    out_oof.to_csv(args.out_dir / "meta_oof_predictions.tsv", sep="\t", index=False)

    # Threshold sweep
    thresholds = np.round(np.arange(0.05, 0.96, 0.01), 3)
    sweep_rows = [_threshold_metrics(sub, meta_prob, t) for t in thresholds]
    sweep_df = pd.DataFrame(sweep_rows)
    sweep_df.to_csv(args.out_dir / "meta_threshold_sweep.tsv", sep="\t", index=False)

    # Mean coefficients across folds (interpretability)
    mean_coef = np.mean(fold_coefs, axis=0)
    coef_df = pd.DataFrame({"feature": feat_cols, "mean_coef": mean_coef})
    coef_df = coef_df.sort_values("mean_coef", ascending=False)
    coef_df.to_csv(args.out_dir / "meta_coefficients.tsv", sep="\t", index=False)

    # Pick operating threshold: maximize specificity (min wrongly-PM) where sensitivity >= floor
    floor = args.sensitivity_floor
    eligible = sweep_df[sweep_df["sensitivity_pct"] >= floor * 100]
    if eligible.empty:
        print(f"\nNo threshold achieves sensitivity >= {floor*100:.0f}% — relaxing to best specificity overall.")
        eligible = sweep_df

    best_row = eligible.loc[eligible["pct_wrongly_called_PM"].idxmin()]

    # Print summary
    print("\n" + "=" * 72)
    print(f"Meta-learner OOF threshold sweep  (sensitivity floor >= {floor*100:.0f}%)")
    print("=" * 72)
    print(f"{'Threshold':>9}  {'WrongPM%':>8}  {'WrongPM/N':>10}  "
          f"{'Sens%':>6}  {'NPV':>5}  {'ctrlWrong%':>10}")
    for _, r in sweep_df.iterrows():
        marker = " <<< selected" if abs(r["threshold"] - best_row["threshold"]) < 1e-6 else ""
        print(f"  {r['threshold']:7.2f}  {r['pct_wrongly_called_PM']:8.1f}  "
              f"{int(r['wrongly_called_PM']):4d}/{int(r['neg_n']):<5d}  "
              f"{r['sensitivity_pct']:6.1f}  {r['npv']:5.3f}  "
              f"{r['ctrl_pct_wrongly_called_PM']:10.1f}{marker}")

    print("\n--- Selected operating point ---")
    for k, v in best_row.items():
        print(f"  {k}: {v}")

    # Compare to K=6 consensus (hardcoded reference from prior analysis)
    print("\n--- Comparison vs K=6 unanimous consensus (from prior run) ---")
    print("  consensus K=6:  wrongly-PM ~14%  |  sensitivity ~76.7%")
    print(f"  meta-learner @{best_row['threshold']}:  wrongly-PM {best_row['pct_wrongly_called_PM']}%  "
          f"|  sensitivity {best_row['sensitivity_pct']}%  |  NPV {best_row['npv']}")

    summary = {
        "n_meta_features": len(feat_cols),
        "n_valid_rows": int(valid.sum()),
        "operating_threshold": float(best_row["threshold"]),
        "sensitivity_floor": floor,
        "sensitivity_pct": float(best_row["sensitivity_pct"]),
        "specificity_pct": float(best_row["specificity_pct"]),
        "npv": float(best_row["npv"]),
        "ppv": float(best_row["ppv"]),
        "neg_wrongly_called_PM": int(best_row["wrongly_called_PM"]),
        "neg_n": int(best_row["neg_n"]),
        "pct_wrongly_called_PM": float(best_row["pct_wrongly_called_PM"]),
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nWrote: {args.out_dir}/meta_oof_predictions.tsv")
    print(f"Wrote: {args.out_dir}/meta_threshold_sweep.tsv")
    print(f"Wrote: {args.out_dir}/meta_coefficients.tsv")
    print(f"Wrote: {args.out_dir}/summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
