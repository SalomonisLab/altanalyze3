#!/usr/bin/env python3
"""Gene-grouped K-fold cross-validation over the full TM-isoform benchmark.

Motivation: the single-split run_tm_isoform_model_comparison.py evaluates
on a fixed 10% gene-hash test split (~69 tm_retained negatives). This
script instead pools train+val+test, performs GroupKFold on (species,
gene_name), and reports metrics aggregated over every fold so that
EVERY negative isoform (and every positive pair) is held out exactly
once. This maximizes evaluation coverage of the 427 true-negative
catalog.

Outputs:
    checkpoints/tm_isoform_kfold_fold_metrics.tsv
    checkpoints/tm_isoform_kfold_summary.tsv
    checkpoints/tm_isoform_kfold_predictions.tsv   (per-sample OOF predictions)
    checkpoints/tm_isoform_kfold_summary.json

Usage:
    python -m phase_d.scripts.run_tm_isoform_kfold_benchmark \
        --models phase_c_residual_mlp_autoencoder,residual_mlp_ae_wide,phase_d_xgboost \
        --k 5 \
        --tune-threshold \
        --literature-tsv phase_a/data/interim/tm_negatives/gene_literature_evidence.tsv
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.ensemble import AdaBoostClassifier, HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.svm import SVC

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

PHASE_D = Path(__file__).resolve().parents[1]
PROCESSED = PHASE_D / "data" / "processed"
CHECKPOINTS = PHASE_D / "checkpoints"
IN_PATH = PROCESSED / "tm_isoform_benchmark_instances.parquet"
SCHEMA_PATH = PROCESSED / "tm_isoform_benchmark_schema.json"

POS_GROUP = "uniprot_tm_positive"
NEG_GROUP_TM_RETAINED = "tm_retained_non_surface_true_negative"
NEG_GROUP_TM_LOST = "tm_lost_no_tm_control"

# Reuse the single-split implementation for all the model training and
# threshold-tuning logic; we just slice the data differently.
from Daedalus.phase_d.scripts.run_tm_isoform_model_comparison import (  # noqa: E402
    MODEL_ORDER,
    SplitArrays,
    _classical_prob,
    _safe_metric,
    _train_deep_binary,
    _tune_threshold_for_negatives,
)


def _prepare_arrays_from_frames(
    train_df: pd.DataFrame,
    val_df: pd.DataFrame,
    test_df: pd.DataFrame,
    feature_columns: list[str],
    categorical_columns: list[str],
) -> tuple[dict[str, SplitArrays], dict[str, object]]:
    scaler = StandardScaler()
    train_numeric = scaler.fit_transform(train_df[feature_columns].to_numpy(dtype=np.float32)).astype(np.float32)
    val_numeric = scaler.transform(val_df[feature_columns].to_numpy(dtype=np.float32)).astype(np.float32)
    test_numeric = scaler.transform(test_df[feature_columns].to_numpy(dtype=np.float32)).astype(np.float32)

    encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    train_cat = encoder.fit_transform(train_df[categorical_columns].astype(str))
    val_cat = encoder.transform(val_df[categorical_columns].astype(str))
    test_cat = encoder.transform(test_df[categorical_columns].astype(str))

    def _map(series: pd.Series) -> dict[str, int]:
        return {v: i for i, v in enumerate(sorted(series.astype(str).unique().tolist()))}

    family_map = _map(train_df["family_class"])
    refsrc_map = _map(train_df["reference_source"])
    change_map = _map(train_df["change_class"])

    def encode(series: pd.Series, mapping: dict[str, int]) -> np.ndarray:
        return series.astype(str).map(lambda x: mapping.get(x, 0)).to_numpy(dtype=np.int64)

    def bundle(frame: pd.DataFrame, x_num: np.ndarray, x_cat: np.ndarray) -> SplitArrays:
        return SplitArrays(
            x_numeric=x_num,
            x_dense=np.concatenate([x_num, x_cat.astype(np.float32)], axis=1).astype(np.float32),
            y=frame["label_binary"].fillna(0).to_numpy(dtype=np.int64),
            family_id=encode(frame["family_class"], family_map),
            refsrc_id=encode(frame["reference_source"], refsrc_map),
            change_id=encode(frame["change_class"], change_map),
            task_id=np.zeros(len(frame), dtype=np.int64),
        )

    arrays = {
        "train": bundle(train_df, train_numeric, train_cat),
        "val": bundle(val_df, val_numeric, val_cat),
        "test": bundle(test_df, test_numeric, test_cat),
    }
    meta = {
        "family_map_size": len(family_map),
        "refsrc_map_size": len(refsrc_map),
        "change_map_size": len(change_map),
        "numeric_dim": int(train_numeric.shape[1]),
        "dense_dim": int(arrays["train"].x_dense.shape[1]),
    }
    return arrays, meta


def _train_and_predict(
    model_name: str,
    arrays: dict[str, SplitArrays],
    meta: dict[str, object],
    device: Any,
    max_epochs: int,
    patience: int,
    batch_size: int,
    seed: int,
    sample_weight: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    train, val, test = arrays["train"], arrays["val"], arrays["test"]
    x_train, y_train = train.x_dense, train.y
    x_val, x_test = val.x_dense, test.x_dense

    pos = float(y_train.sum())
    neg = float(len(y_train) - pos)

    if model_name == "phase_d_logistic":
        return _classical_prob(
            LogisticRegression(solver="saga", C=1.0, max_iter=1000,
                               class_weight="balanced", random_state=seed),
            x_train, y_train, x_val, val.y, x_test, sample_weight=sample_weight,
        )
    if model_name == "phase_c_elasticnet":
        return _classical_prob(
            SGDClassifier(loss="log_loss", penalty="elasticnet", alpha=1e-4,
                          l1_ratio=0.15, max_iter=3000, tol=1e-3,
                          class_weight="balanced", random_state=seed),
            x_train, y_train, x_val, val.y, x_test, sample_weight=sample_weight,
        )
    if model_name == "phase_d_hist_gradient_boosting":
        hgb_sw = np.ones(len(y_train), dtype=np.float32)
        if neg > 0:
            hgb_sw[y_train == 0] = pos / neg
        if sample_weight is not None:
            hgb_sw = hgb_sw * sample_weight
        return _classical_prob(
            HistGradientBoostingClassifier(max_depth=6, learning_rate=0.05,
                                           max_iter=400, early_stopping=True,
                                           validation_fraction=0.15,
                                           random_state=seed),
            x_train, y_train, x_val, val.y, x_test, sample_weight=hgb_sw,
        )
    if model_name == "phase_d_xgboost":
        import xgboost as xgb
        return _classical_prob(
            xgb.XGBClassifier(n_estimators=400, max_depth=5, learning_rate=0.05,
                              subsample=0.8, colsample_bytree=0.8, reg_lambda=1.0,
                              min_child_weight=1.0,
                              scale_pos_weight=neg / max(pos, 1.0),
                              objective="binary:logistic", eval_metric="aucpr",
                              tree_method="hist", random_state=seed, n_jobs=1),
            x_train, y_train, x_val, val.y, x_test, sample_weight=sample_weight,
        )
    if model_name == "phase_c_random_forest":
        return _classical_prob(
            RandomForestClassifier(n_estimators=400, class_weight="balanced_subsample",
                                   min_samples_leaf=2, random_state=seed, n_jobs=-1),
            x_train, y_train, x_val, val.y, x_test, sample_weight=sample_weight,
        )
    if model_name == "phase_c_svm_rbf":
        return _classical_prob(
            SVC(kernel="rbf", probability=True, class_weight="balanced",
                gamma="scale", C=1.0, random_state=seed),
            x_train, y_train, x_val, val.y, x_test, sample_weight=sample_weight,
        )
    if model_name == "phase_c_adaboost":
        ada_sw = np.ones(len(y_train), dtype=np.float32)
        if neg > 0:
            ada_sw[y_train == 0] = pos / neg
        if sample_weight is not None:
            ada_sw = ada_sw * sample_weight
        return _classical_prob(
            AdaBoostClassifier(n_estimators=200, learning_rate=0.5, random_state=seed),
            x_train, y_train, x_val, val.y, x_test, sample_weight=ada_sw,
        )
    if model_name == "phase_c_knn":
        return _classical_prob(
            KNeighborsClassifier(n_neighbors=25, weights="distance"),
            x_train, y_train, x_val, val.y, x_test, sample_weight=None,
        )
    # Deep models
    if model_name == "phase_c_deepimmuno_cnn":
        from Daedalus.phase_c.models import DeepImmunoStyleCNN
        m = DeepImmunoStyleCNN(
            n_numeric=meta["numeric_dim"], n_families=meta["family_map_size"],
            n_reference_sources=meta["refsrc_map_size"],
            n_change_classes=meta["change_map_size"],
        ).to(device)
        return _train_deep_binary(m, train, val, test, device, max_epochs, patience,
                                  batch_size, sample_weight=sample_weight)
    if model_name == "phase_c_deepimmuno_cnn_autoencoder":
        from Daedalus.phase_c.models import DeepImmunoStyleCNNAutoencoder
        m = DeepImmunoStyleCNNAutoencoder(
            n_numeric=meta["numeric_dim"], n_families=meta["family_map_size"],
            n_reference_sources=meta["refsrc_map_size"],
            n_change_classes=meta["change_map_size"],
        ).to(device)
        return _train_deep_binary(m, train, val, test, device, max_epochs, patience,
                                  batch_size, recon_weight=0.10, sample_weight=sample_weight)
    if model_name == "phase_c_residual_mlp_autoencoder":
        from Daedalus.phase_c.models import ResidualMLPAutoencoder
        m = ResidualMLPAutoencoder(
            n_numeric=meta["numeric_dim"], n_families=meta["family_map_size"],
            n_reference_sources=meta["refsrc_map_size"],
            n_change_classes=meta["change_map_size"],
        ).to(device)
        return _train_deep_binary(m, train, val, test, device, max_epochs, patience,
                                  batch_size, recon_weight=0.10, sample_weight=sample_weight)
    if model_name == "phase_d_family_multitask":
        from Daedalus.phase_d.models.family_multitask_net import FamilyMultitaskNet
        m = FamilyMultitaskNet(
            n_numeric=meta["numeric_dim"], n_tasks=1,
            n_families=meta["family_map_size"],
            n_reference_sources=meta["refsrc_map_size"],
            n_change_classes=meta["change_map_size"],
        ).to(device)
        return _train_deep_binary(m, train, val, test, device, max_epochs, patience,
                                  batch_size, recon_weight=0.10, use_task_ids=True,
                                  sample_weight=sample_weight)
    # variants
    variant_configs = {
        "residual_mlp_ae_wide":        (dict(hidden_dim=256, latent_dim=64, dropout=0.15), 0.10),
        "residual_mlp_ae_deep":        (dict(hidden_dim=192, latent_dim=48, dropout=0.15, n_blocks=3), 0.10),
        "residual_mlp_ae_highdrop":    (dict(hidden_dim=192, latent_dim=48, dropout=0.30), 0.10),
        "residual_mlp_ae_strongrecon": (dict(hidden_dim=192, latent_dim=48, dropout=0.15), 0.30),
        "residual_mlp_ae_narrow":      (dict(hidden_dim=128, latent_dim=32, dropout=0.15), 0.10),
    }
    if model_name in variant_configs:
        cfg, recon_w = variant_configs[model_name]
        from Daedalus.phase_c.models.residual_mlp_variants import ResidualMLPAutoencoderVariant
        m = ResidualMLPAutoencoderVariant(
            n_numeric=meta["numeric_dim"],
            n_families=meta["family_map_size"],
            n_reference_sources=meta["refsrc_map_size"],
            n_change_classes=meta["change_map_size"],
            **cfg,
        ).to(device)
        return _train_deep_binary(m, train, val, test, device, max_epochs, patience,
                                  batch_size, recon_weight=recon_w, sample_weight=sample_weight)
    raise ValueError(f"Unknown model: {model_name}")


def _metrics(y_true: np.ndarray, y_prob: np.ndarray, threshold: float) -> dict[str, float]:
    y_true = np.nan_to_num(y_true, nan=0).astype(int)
    pred = (y_prob >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, pred, labels=[0, 1]).ravel()
    return {
        "threshold": float(threshold),
        "ap": _safe_metric(average_precision_score, y_true, y_prob),
        "auroc": _safe_metric(roc_auc_score, y_true, y_prob),
        "balanced_accuracy": _safe_metric(balanced_accuracy_score, y_true, pred),
        "specificity": float(tn / max(tn + fp, 1)),
        "sensitivity": float(tp / max(tp + fn, 1)),
        "npv": float(tn / max(tn + fn, 1)),
        "ppv": float(tp / max(tp + fp, 1)),
        "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp),
        "n": int(len(y_true)),
        "n_pos": int(y_true.sum()),
        "n_neg": int(len(y_true) - y_true.sum()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Gene-grouped K-fold CV over the full TM-isoform benchmark.")
    parser.add_argument("--models", default="phase_c_residual_mlp_autoencoder,residual_mlp_ae_wide,residual_mlp_ae_deep,residual_mlp_ae_highdrop,residual_mlp_ae_strongrecon,residual_mlp_ae_narrow,phase_d_xgboost")
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--max-epochs", type=int, default=24)
    parser.add_argument("--patience", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--tune-threshold", action="store_true")
    parser.add_argument("--sensitivity-floor", type=float, default=0.85)
    parser.add_argument("--literature-tsv", type=Path, default=None)
    parser.add_argument(
        "--hard-neg-weight",
        type=float,
        default=4.0,
        help="Sample weight multiplier applied to tm_retained_non_surface_true_negative "
             "training instances relative to all other samples (weight=1.0). "
             "Upweights hard negatives independently of class-balance weighting.",
    )
    parser.add_argument("--out-dir", type=Path, default=CHECKPOINTS)
    args = parser.parse_args()

    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    warnings.filterwarnings("ignore", message=".*'n_jobs' has no effect.*")

    requested = [m.strip() for m in args.models.split(",") if m.strip()]
    unknown = [m for m in requested if m not in MODEL_ORDER]
    if unknown:
        raise ValueError(f"Unknown model ids: {unknown}")

    deep_ids = {
        "phase_c_deepimmuno_cnn", "phase_c_deepimmuno_cnn_autoencoder",
        "phase_c_residual_mlp_autoencoder", "phase_d_family_multitask",
        "residual_mlp_ae_wide", "residual_mlp_ae_deep", "residual_mlp_ae_highdrop",
        "residual_mlp_ae_strongrecon", "residual_mlp_ae_narrow",
    }
    if any(m in deep_ids for m in requested):
        from Daedalus.phase_d.models.family_multitask_net import choose_device
        device = choose_device(args.device)
    else:
        class _Cpu:
            type = "cpu"
        device = _Cpu()

    with SCHEMA_PATH.open() as h:
        schema = json.load(h)
    df = pd.read_parquet(IN_PATH).reset_index(drop=True)
    feature_columns = schema["feature_columns"]
    categorical_columns = schema["categorical_columns"]

    # Literature weighting
    lit_weight_by_gene: dict[str, float] = {}
    if args.literature_tsv is not None:
        lit = pd.read_csv(args.literature_tsv, sep="\t")
        weight_by_evidence = {"HIGH": 1.0, "MEDIUM": 0.7, "LOW": 0.4,
                              "AMBIGUOUS": 0.3, "NONE": 0.2, "NO": 0.0}
        for _, row in lit.iterrows():
            conf = str(row.get("confidence", "NONE"))
            lit_weight_by_gene[str(row["gene_name"])] = weight_by_evidence.get(conf, 0.5)

    # Gene-grouped K-fold: group = species|gene_name
    groups = (df["species"].astype(str) + "|" + df["gene_name"].astype(str)).to_numpy()
    gkf = GroupKFold(n_splits=args.k)

    per_fold_rows: list[dict[str, object]] = []
    oof_rows: list[dict[str, object]] = []
    y_all = df["label_binary"].fillna(0).to_numpy(dtype=np.int64)

    # Accumulator for OOF predictions per model
    oof_by_model: dict[str, np.ndarray] = {m: np.full(len(df), np.nan, dtype=np.float64) for m in requested}

    for fold_idx, (train_val_idx, test_idx) in enumerate(gkf.split(df, y_all, groups)):
        test_df = df.iloc[test_idx].reset_index(drop=True)
        tv_df = df.iloc[train_val_idx]
        tv_groups = groups[train_val_idx]
        # Inner gene-grouped split for val (~10% of train+val)
        inner = GroupKFold(n_splits=9)
        _, (inner_train_idx, inner_val_idx) = next(enumerate(inner.split(tv_df, tv_df["label_binary"].fillna(0), tv_groups)))
        train_df = tv_df.iloc[inner_train_idx].reset_index(drop=True)
        val_df = tv_df.iloc[inner_val_idx].reset_index(drop=True)

        arrays, meta = _prepare_arrays_from_frames(
            train_df, val_df, test_df, feature_columns, categorical_columns
        )

        # Per-sample weights for this fold's training set.
        # Hard negatives (tm_retained) are upweighted by --hard-neg-weight
        # relative to all other samples, then optionally modulated by
        # literature evidence confidence for those same genes.
        w = np.ones(len(train_df), dtype=np.float32)
        if "benchmark_group" in train_df.columns:
            hard_neg_mask = (
                train_df["benchmark_group"] == NEG_GROUP_TM_RETAINED
            ).to_numpy()
            w[hard_neg_mask] = float(args.hard_neg_weight)
        if lit_weight_by_gene:
            for i, row in train_df.iterrows():
                if int(row["label_binary"]) == 1:
                    continue
                g = str(row.get("gene_name", ""))
                if g in lit_weight_by_gene:
                    w[i] *= lit_weight_by_gene[g]
        train_sw = w if (w != 1.0).any() else None

        for model_name in requested:
            try:
                val_prob, test_prob = _train_and_predict(
                    model_name, arrays, meta, device,
                    args.max_epochs, args.patience, args.batch_size, args.seed,
                    sample_weight=train_sw,
                )
            except Exception as exc:
                print(f"[fold={fold_idx} model={model_name}] FAILED: {exc}", flush=True)
                continue

            # Store OOF test predictions
            oof_by_model[model_name][test_idx] = test_prob

            # Threshold = 0.5
            row = _metrics(test_df["label_binary"].fillna(0).to_numpy(), test_prob, 0.5)
            row["fold"] = int(fold_idx)
            row["model"] = model_name
            row["threshold_mode"] = "fixed_0.5"
            row["scope"] = "overall"
            per_fold_rows.append(row)
            # Per-scope
            for scope in [NEG_GROUP_TM_RETAINED, NEG_GROUP_TM_LOST]:
                mask = test_df["benchmark_group"].isin([POS_GROUP, scope]).to_numpy()
                if mask.sum() == 0:
                    continue
                m = _metrics(
                    test_df.loc[mask, "label_binary"].fillna(0).to_numpy(),
                    test_prob[mask], 0.5,
                )
                m["fold"] = int(fold_idx)
                m["model"] = model_name
                m["threshold_mode"] = "fixed_0.5"
                m["scope"] = scope
                per_fold_rows.append(m)

            # Tuned threshold
            if args.tune_threshold:
                thr = _tune_threshold_for_negatives(
                    y_val=arrays["val"].y, p_val=val_prob,
                    sensitivity_floor=float(args.sensitivity_floor),
                )
                row_t = _metrics(test_df["label_binary"].fillna(0).to_numpy(), test_prob, thr)
                row_t["fold"] = int(fold_idx)
                row_t["model"] = model_name
                row_t["threshold_mode"] = f"tuned_npv_sens>={args.sensitivity_floor:.2f}"
                row_t["scope"] = "overall"
                per_fold_rows.append(row_t)
                for scope in [NEG_GROUP_TM_RETAINED, NEG_GROUP_TM_LOST]:
                    mask = test_df["benchmark_group"].isin([POS_GROUP, scope]).to_numpy()
                    if mask.sum() == 0:
                        continue
                    m = _metrics(
                        test_df.loc[mask, "label_binary"].to_numpy(),
                        test_prob[mask], thr,
                    )
                    m["fold"] = int(fold_idx)
                    m["model"] = model_name
                    m["threshold_mode"] = f"tuned_npv_sens>={args.sensitivity_floor:.2f}"
                    m["scope"] = scope
                    per_fold_rows.append(m)
            print(f"[fold={fold_idx} model={model_name}] overall AP={row['ap']:.4f} "
                  f"AUROC={row['auroc']:.4f} Spec={row['specificity']:.3f} "
                  f"NPV={row['npv']:.3f}", flush=True)

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    fold_df = pd.DataFrame(per_fold_rows)
    fold_df.to_csv(out_dir / "tm_isoform_kfold_fold_metrics.tsv", sep="\t", index=False)

    # Aggregate: mean +/- std across folds per (model, scope, threshold_mode)
    grp_cols = ["model", "scope", "threshold_mode"]
    metric_cols = ["ap", "auroc", "balanced_accuracy", "specificity", "sensitivity", "npv", "ppv"]
    agg_rows = []
    for (m, s, tm), g in fold_df.groupby(grp_cols):
        row = {"model": m, "scope": s, "threshold_mode": tm, "n_folds": int(len(g))}
        for c in metric_cols:
            row[f"{c}_mean"] = float(g[c].mean())
            row[f"{c}_std"] = float(g[c].std())
        row["tp_sum"] = int(g["tp"].sum())
        row["tn_sum"] = int(g["tn"].sum())
        row["fp_sum"] = int(g["fp"].sum())
        row["fn_sum"] = int(g["fn"].sum())
        row["n_sum"] = int(g["n"].sum())
        row["n_pos_sum"] = int(g["n_pos"].sum())
        row["n_neg_sum"] = int(g["n_neg"].sum())
        agg_rows.append(row)
    agg_df = pd.DataFrame(agg_rows).sort_values(
        ["scope", "threshold_mode", "npv_mean", "specificity_mean"],
        ascending=[True, True, False, False],
    )
    agg_df.to_csv(out_dir / "tm_isoform_kfold_summary.tsv", sep="\t", index=False)

    # OOF predictions per sample for later analysis
    oof_frame = df[["species", "gene_name", "reference_transcript_id",
                    "alternative_transcript_id", "benchmark_group", "label_binary",
                    "negative_isoform_ids"]].copy()
    for m, probs in oof_by_model.items():
        oof_frame[f"prob_{m}"] = probs
    oof_frame.to_csv(out_dir / "tm_isoform_kfold_predictions.tsv", sep="\t", index=False)

    summary = {
        "k": int(args.k),
        "models": requested,
        "tune_threshold": bool(args.tune_threshold),
        "sensitivity_floor": float(args.sensitivity_floor),
        "literature_tsv": str(args.literature_tsv) if args.literature_tsv else None,
        "n_total_rows": int(len(df)),
        "n_positives": int((df["label_binary"].fillna(0) == 1).sum()),
        "n_negatives": int((df["label_binary"].fillna(0) == 0).sum()),
        "n_tm_retained_negatives": int((df["benchmark_group"] == NEG_GROUP_TM_RETAINED).sum()),
        "best_per_scope": {},
    }
    for scope in ["overall", NEG_GROUP_TM_RETAINED, NEG_GROUP_TM_LOST]:
        sub = agg_df[agg_df["scope"] == scope]
        if sub.empty:
            continue
        # Prefer tuned threshold if available
        tuned_sub = sub[sub["threshold_mode"].str.startswith("tuned_")]
        chosen = tuned_sub if not tuned_sub.empty else sub
        chosen = chosen.sort_values(["npv_mean", "specificity_mean"], ascending=[False, False])
        top = chosen.iloc[0]
        summary["best_per_scope"][scope] = {
            "model": str(top["model"]),
            "threshold_mode": str(top["threshold_mode"]),
            "npv_mean": float(top["npv_mean"]),
            "specificity_mean": float(top["specificity_mean"]),
            "sensitivity_mean": float(top["sensitivity_mean"]),
            "auroc_mean": float(top["auroc_mean"]),
            "ap_mean": float(top["ap_mean"]),
        }

    (out_dir / "tm_isoform_kfold_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[ok] Wrote {out_dir / 'tm_isoform_kfold_fold_metrics.tsv'}")
    print(f"[ok] Wrote {out_dir / 'tm_isoform_kfold_summary.tsv'}")
    print(f"[ok] Wrote {out_dir / 'tm_isoform_kfold_predictions.tsv'}")
    print(f"[ok] Wrote {out_dir / 'tm_isoform_kfold_summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
