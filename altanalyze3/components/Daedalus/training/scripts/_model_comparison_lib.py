#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import random
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.ensemble import AdaBoostClassifier, HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import average_precision_score, balanced_accuracy_score, confusion_matrix, roc_auc_score
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
OUT_TSV = CHECKPOINTS / "tm_isoform_model_comparison_metrics.tsv"
OUT_MD = CHECKPOINTS / "tm_isoform_model_comparison_report.md"
OUT_JSON = CHECKPOINTS / "tm_isoform_model_comparison_summary.json"

POS_GROUP = "uniprot_tm_positive"
NEG_GROUP_TM_RETAINED = "tm_retained_non_surface_true_negative"
NEG_GROUP_TM_LOST = "tm_lost_no_tm_control"
MODEL_ORDER = [
    "phase_d_logistic",
    "phase_c_elasticnet",
    "phase_d_hist_gradient_boosting",
    "phase_d_xgboost",
    "phase_c_random_forest",
    "phase_c_svm_rbf",
    "phase_c_adaboost",
    "phase_c_knn",
    "phase_c_deepimmuno_cnn",
    "phase_c_deepimmuno_cnn_autoencoder",
    "phase_c_residual_mlp_autoencoder",
    "residual_mlp_ae_wide",
    "residual_mlp_ae_deep",
    "residual_mlp_ae_highdrop",
    "residual_mlp_ae_strongrecon",
    "residual_mlp_ae_narrow",
    "phase_d_family_multitask",
]


@dataclass
class SplitArrays:
    x_numeric: np.ndarray
    x_dense: np.ndarray
    y: np.ndarray
    family_id: np.ndarray
    refsrc_id: np.ndarray
    change_id: np.ndarray
    task_id: np.ndarray


def _set_seeds(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except Exception:
        pass


def _safe_metric(fn, y_true: np.ndarray, y_score: np.ndarray) -> float:
    try:
        return float(fn(y_true, y_score))
    except Exception:
        return float("nan")


def _tune_threshold_for_negatives(
    y_val: np.ndarray,
    p_val: np.ndarray,
    sensitivity_floor: float = 0.85,
) -> float:
    """Select the probability cutoff that maximises NPV (TN/(TN+FN))
    subject to sensitivity >= sensitivity_floor on validation data.

    Positive class = 1 (PM-localized). Negative class = 0 (true-negative TM
    isoform). "Call negative when prob < threshold."
    """
    if len(y_val) == 0:
        return 0.5
    best_thr = 0.5
    best_npv = -1.0
    # Sweep candidate thresholds from the actual predicted probs
    candidates = np.unique(np.concatenate([
        np.linspace(0.05, 0.95, 19),
        np.quantile(p_val, np.linspace(0.05, 0.95, 19)),
    ]))
    for thr in candidates:
        pred = (p_val >= thr).astype(int)
        tp = int(((pred == 1) & (y_val == 1)).sum())
        fn = int(((pred == 0) & (y_val == 1)).sum())
        tn = int(((pred == 0) & (y_val == 0)).sum())
        fp = int(((pred == 1) & (y_val == 0)).sum())
        sens = tp / max(tp + fn, 1)
        if sens < sensitivity_floor:
            continue
        npv = tn / max(tn + fn, 1)
        if npv > best_npv or (npv == best_npv and thr < best_thr):
            best_npv = npv
            best_thr = float(thr)
    return best_thr


def _scope_metrics(
    frame: pd.DataFrame,
    probs: np.ndarray,
    scope: str,
    threshold: float = 0.5,
) -> dict[str, object]:
    if scope == "overall":
        sub = frame.copy()
    else:
        sub = frame[
            frame["benchmark_group"].isin([POS_GROUP, scope])
        ].copy()
    y_true = sub["label_binary"].to_numpy(dtype=np.int32)
    y_prob = probs[sub.index.to_numpy()]
    pred = (y_prob >= threshold).astype(int)
    labels = [0, 1]
    tn, fp, fn, tp = confusion_matrix(y_true, pred, labels=labels).ravel()
    specificity = tn / max(tn + fp, 1)
    sensitivity = tp / max(tp + fn, 1)
    npv = tn / max(tn + fn, 1)
    ppv = tp / max(tp + fp, 1)
    return {
        "scope": scope,
        "n": int(len(sub)),
        "positives": int(y_true.sum()),
        "negatives": int(len(sub) - y_true.sum()),
        "threshold": float(threshold),
        "average_precision": _safe_metric(average_precision_score, y_true, y_prob),
        "auroc": _safe_metric(roc_auc_score, y_true, y_prob),
        "balanced_accuracy": _safe_metric(balanced_accuracy_score, y_true, pred),
        "specificity": float(specificity),
        "sensitivity": float(sensitivity),
        "npv": float(npv),
        "ppv": float(ppv),
        "tp": int(tp),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
    }


def _encode_map(train_values: pd.Series) -> dict[str, int]:
    return {value: idx for idx, value in enumerate(sorted(train_values.astype(str).unique().tolist()))}


def _encode_with_map(series: pd.Series, mapping: dict[str, int]) -> np.ndarray:
    return series.astype(str).map(lambda x: mapping.get(x, 0)).to_numpy(dtype=np.int64)


def _prepare_arrays(df: pd.DataFrame, feature_columns: list[str], categorical_columns: list[str]) -> tuple[dict[str, SplitArrays], dict[str, object]]:
    train_df = df[df["split"] == "train"].copy()
    val_df = df[df["split"] == "val"].copy()
    test_df = df[df["split"] == "test"].copy()

    scaler = StandardScaler()
    train_numeric = scaler.fit_transform(train_df[feature_columns].to_numpy(dtype=np.float32)).astype(np.float32)
    val_numeric = scaler.transform(val_df[feature_columns].to_numpy(dtype=np.float32)).astype(np.float32)
    test_numeric = scaler.transform(test_df[feature_columns].to_numpy(dtype=np.float32)).astype(np.float32)

    encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    train_cat = encoder.fit_transform(train_df[categorical_columns].astype(str))
    val_cat = encoder.transform(val_df[categorical_columns].astype(str))
    test_cat = encoder.transform(test_df[categorical_columns].astype(str))

    family_map = _encode_map(train_df["family_class"])
    refsrc_map = _encode_map(train_df["reference_source"])
    change_map = _encode_map(train_df["change_class"])
    task_map = {"cell_surface_localization": 0}

    def bundle(frame: pd.DataFrame, x_num: np.ndarray, x_cat: np.ndarray) -> SplitArrays:
        return SplitArrays(
            x_numeric=x_num,
            x_dense=np.concatenate([x_num, x_cat.astype(np.float32)], axis=1).astype(np.float32),
            y=frame["label_binary"].to_numpy(dtype=np.int64),
            family_id=_encode_with_map(frame["family_class"], family_map),
            refsrc_id=_encode_with_map(frame["reference_source"], refsrc_map),
            change_id=_encode_with_map(frame["change_class"], change_map),
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
        "task_map_size": len(task_map),
        "numeric_dim": int(train_numeric.shape[1]),
        "dense_dim": int(arrays["train"].x_dense.shape[1]),
    }
    return arrays, meta


def _classical_prob(
    model,
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_val: np.ndarray | None,
    y_val: np.ndarray | None,
    x_test: np.ndarray,
    sample_weight: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit model and return (val_prob, test_prob) for threshold tuning + test eval."""
    fit_kwargs: dict[str, np.ndarray] = {}
    if sample_weight is not None:
        try:
            model.fit(x_train, y_train, sample_weight=sample_weight)
        except TypeError:
            model.fit(x_train, y_train)
    else:
        model.fit(x_train, y_train)

    def _prob(x: np.ndarray) -> np.ndarray:
        if hasattr(model, "predict_proba"):
            return model.predict_proba(x)[:, 1]
        score = model.decision_function(x)
        return 1.0 / (1.0 + np.exp(-score))

    val_prob = _prob(x_val) if x_val is not None else np.array([])
    test_prob = _prob(x_test)
    return val_prob, test_prob


def _train_deep_binary(
    model: Any,
    train: SplitArrays,
    val: SplitArrays,
    test: SplitArrays,
    device: Any,
    max_epochs: int,
    patience: int,
    batch_size: int,
    recon_weight: float = 0.0,
    use_task_ids: bool = False,
    sample_weight: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Train a deep model and return (val_prob, test_prob)."""
    import torch
    from torch.utils.data import DataLoader, TensorDataset

    pos = float(train.y.sum())
    neg = float(len(train.y) - pos)
    pos_weight = torch.tensor([neg / max(pos, 1.0)], device=device, dtype=torch.float32)
    # Use per-sample reduction="none" when sample_weight is provided so we
    # can multiply by the weight before mean.
    per_sample = sample_weight is not None
    criterion = torch.nn.BCEWithLogitsLoss(
        pos_weight=pos_weight,
        reduction="none" if per_sample else "mean",
    )
    recon_criterion = torch.nn.MSELoss(reduction="none" if per_sample else "mean")
    optimizer = torch.optim.AdamW(model.parameters(), lr=7e-4, weight_decay=2e-4)

    sw = (
        torch.tensor(sample_weight, dtype=torch.float32)
        if sample_weight is not None
        else torch.ones(len(train.y), dtype=torch.float32)
    )

    if use_task_ids:
        dataset = TensorDataset(
            torch.tensor(train.x_numeric, dtype=torch.float32),
            torch.tensor(train.task_id, dtype=torch.long),
            torch.tensor(train.family_id, dtype=torch.long),
            torch.tensor(train.refsrc_id, dtype=torch.long),
            torch.tensor(train.change_id, dtype=torch.long),
            torch.tensor(train.y, dtype=torch.float32),
            sw,
        )
    else:
        dataset = TensorDataset(
            torch.tensor(train.x_numeric, dtype=torch.float32),
            torch.tensor(train.family_id, dtype=torch.long),
            torch.tensor(train.refsrc_id, dtype=torch.long),
            torch.tensor(train.change_id, dtype=torch.long),
            torch.tensor(train.y, dtype=torch.float32),
            sw,
        )
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    def predict(bundle: SplitArrays) -> np.ndarray:
        model.eval()
        with torch.no_grad():
            if use_task_ids:
                logits, _, _ = model(
                    torch.tensor(bundle.x_numeric, dtype=torch.float32, device=device),
                    torch.tensor(bundle.task_id, dtype=torch.long, device=device),
                    torch.tensor(bundle.family_id, dtype=torch.long, device=device),
                    torch.tensor(bundle.refsrc_id, dtype=torch.long, device=device),
                    torch.tensor(bundle.change_id, dtype=torch.long, device=device),
                )
            else:
                out = model(
                    torch.tensor(bundle.x_numeric, dtype=torch.float32, device=device),
                    torch.tensor(bundle.family_id, dtype=torch.long, device=device),
                    torch.tensor(bundle.refsrc_id, dtype=torch.long, device=device),
                    torch.tensor(bundle.change_id, dtype=torch.long, device=device),
                )
                logits = out[0]
            return torch.sigmoid(logits).detach().cpu().numpy()

    best_state = None
    best_ap = -math.inf
    stale = 0
    for _ in range(max_epochs):
        model.train()
        for batch in loader:
            optimizer.zero_grad(set_to_none=True)
            if use_task_ids:
                numeric, task_id, family_id, refsrc_id, change_id, labels, weights = batch
                numeric = numeric.to(device)
                task_id = task_id.to(device)
                family_id = family_id.to(device)
                refsrc_id = refsrc_id.to(device)
                change_id = change_id.to(device)
                labels = labels.to(device)
                weights = weights.to(device)
                logits, recon, _ = model(numeric, task_id, family_id, refsrc_id, change_id)
            else:
                numeric, family_id, refsrc_id, change_id, labels, weights = batch
                numeric = numeric.to(device)
                family_id = family_id.to(device)
                refsrc_id = refsrc_id.to(device)
                change_id = change_id.to(device)
                labels = labels.to(device)
                weights = weights.to(device)
                out = model(numeric, family_id, refsrc_id, change_id)
                logits = out[0]
                recon = out[1] if len(out) > 2 else None
            loss = criterion(logits, labels)
            if per_sample:
                loss = (loss * weights).sum() / weights.sum().clamp_min(1e-6)
            if recon is not None and recon_weight > 0.0:
                rloss = recon_criterion(recon, numeric)
                if per_sample:
                    rloss = rloss.mean(dim=1)
                    rloss = (rloss * weights).sum() / weights.sum().clamp_min(1e-6)
                loss = loss + recon_weight * rloss
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=2.0)
            optimizer.step()

        val_prob = predict(val)
        val_ap = _safe_metric(average_precision_score, val.y, val_prob)
        if val_ap > best_ap:
            best_ap = val_ap
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
            if stale >= patience:
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    return predict(val), predict(test)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run model comparison on the dedicated TM isoform benchmark.")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--max-epochs", type=int, default=24)
    parser.add_argument("--patience", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--models",
        default="all",
        help="Comma-separated model ids to run, or 'all'.",
    )
    parser.add_argument(
        "--tune-threshold",
        action="store_true",
        help="In addition to the fixed-0.5 threshold metrics, emit metrics at "
        "a per-model threshold chosen on val to maximise NPV subject to "
        "sensitivity >= --sensitivity-floor.",
    )
    parser.add_argument(
        "--sensitivity-floor",
        type=float,
        default=0.85,
        help="Sensitivity floor for NPV-optimised threshold selection.",
    )
    parser.add_argument(
        "--literature-tsv",
        type=Path,
        default=None,
        help="Path to gene_literature_evidence.tsv. When provided, negative "
        "training instances are weighted by literature confidence.",
    )
    parser.add_argument("--out-tsv", type=Path, default=OUT_TSV)
    parser.add_argument("--out-md", type=Path, default=OUT_MD)
    parser.add_argument("--out-json", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    _set_seeds(args.seed)
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    warnings.filterwarnings("ignore", message=".*'penalty' was deprecated.*")
    warnings.filterwarnings("ignore", message=".*'n_jobs' has no effect.*")
    if str(args.models).strip().lower() == "all":
        requested_models = list(MODEL_ORDER)
    else:
        requested_models = [token.strip() for token in str(args.models).split(",") if token.strip()]
    unknown_models = sorted(set(requested_models) - set(MODEL_ORDER))
    if unknown_models:
        raise ValueError(f"Unknown model ids: {', '.join(unknown_models)}")
    deep_model_ids = {
        "phase_c_deepimmuno_cnn",
        "phase_c_deepimmuno_cnn_autoencoder",
        "phase_c_residual_mlp_autoencoder",
        "phase_d_family_multitask",
        "residual_mlp_ae_wide",
        "residual_mlp_ae_deep",
        "residual_mlp_ae_highdrop",
        "residual_mlp_ae_strongrecon",
        "residual_mlp_ae_narrow",
    }
    if any(model_name in deep_model_ids for model_name in requested_models):
        from Daedalus.models.family_multitask_net import choose_device

        device = choose_device(args.device)
    else:
        class _CpuDevice:
            type = "cpu"

        device = _CpuDevice()
    with SCHEMA_PATH.open() as handle:
        schema = json.load(handle)
    df = pd.read_parquet(IN_PATH).reset_index(drop=True)
    feature_columns = schema["feature_columns"]
    categorical_columns = schema["categorical_columns"]
    arrays, meta = _prepare_arrays(df, feature_columns, categorical_columns)

    train = arrays["train"]
    val = arrays["val"]
    test = arrays["test"]
    train_frame = df[df["split"] == "train"].copy().reset_index(drop=True)
    val_frame = df[df["split"] == "val"].copy().reset_index(drop=True)
    test_frame = df[df["split"] == "test"].copy().reset_index(drop=True)

    # --- literature-weighted sample weights ----------------------------
    train_sample_weight: np.ndarray | None = None
    if args.literature_tsv is not None:
        lit = pd.read_csv(args.literature_tsv, sep="\t")
        lit_map = {row["gene_name"]: row for _, row in lit.iterrows()}
        weight_by_evidence = {"HIGH": 1.0, "MEDIUM": 0.7, "LOW": 0.4,
                              "AMBIGUOUS": 0.3, "NONE": 0.2, "NO": 0.0}
        weights = np.ones(len(train_frame), dtype=np.float32)
        n_neg = 0
        n_down = 0
        for idx, row in train_frame.iterrows():
            if int(row["label_binary"]) == 1:
                continue
            n_neg += 1
            gene = str(row.get("gene_name", ""))
            ev = lit_map.get(gene)
            if ev is None:
                continue
            conf = str(ev.get("confidence", "NONE"))
            w = weight_by_evidence.get(conf, 0.5)
            if w < 1.0:
                n_down += 1
            weights[idx] = w
        train_sample_weight = weights
        print(f"[lit-weight] negatives={n_neg} downweighted={n_down}")

    models: list[tuple[str, np.ndarray, np.ndarray]] = []
    x_train = train.x_dense
    y_train = train.y
    x_val = val.x_dense
    y_val = val.y
    x_test = test.x_dense

    if "phase_d_logistic" in requested_models:
        vp, tp = _classical_prob(
            LogisticRegression(
                solver="saga",
                C=1.0,
                max_iter=1000,
                class_weight="balanced",
                random_state=args.seed,
            ),
            x_train, y_train, x_val, y_val, x_test,
            sample_weight=train_sample_weight,
        )
        models.append(("phase_d_logistic", vp, tp))
    if "phase_c_elasticnet" in requested_models:
        vp, tp = _classical_prob(
            SGDClassifier(
                loss="log_loss",
                penalty="elasticnet",
                alpha=1e-4,
                l1_ratio=0.15,
                max_iter=3000,
                tol=1e-3,
                class_weight="balanced",
                random_state=args.seed,
            ),
            x_train, y_train, x_val, y_val, x_test,
            sample_weight=train_sample_weight,
        )
        models.append(("phase_c_elasticnet", vp, tp))
    if "phase_d_hist_gradient_boosting" in requested_models:
        # HGB has no class_weight; simulate via sample_weight
        hgb_sw = np.ones(len(y_train), dtype=np.float32)
        pos = float(y_train.sum())
        neg = float(len(y_train) - pos)
        if neg > 0:
            hgb_sw[y_train == 0] = pos / max(neg, 1.0)
        if train_sample_weight is not None:
            hgb_sw = hgb_sw * train_sample_weight
        vp, tp = _classical_prob(
            HistGradientBoostingClassifier(
                max_depth=6,
                learning_rate=0.05,
                max_iter=400,
                early_stopping=True,
                validation_fraction=0.15,
                random_state=args.seed,
            ),
            x_train, y_train, x_val, y_val, x_test,
            sample_weight=hgb_sw,
        )
        models.append(("phase_d_hist_gradient_boosting", vp, tp))
    if "phase_d_xgboost" in requested_models:
        import xgboost as xgb

        pos = float(y_train.sum())
        neg = float(len(y_train) - pos)
        scale_pos_weight = neg / max(pos, 1.0)
        vp, tp = _classical_prob(
            xgb.XGBClassifier(
                n_estimators=400,
                max_depth=5,
                learning_rate=0.05,
                subsample=0.8,
                colsample_bytree=0.8,
                reg_lambda=1.0,
                min_child_weight=1.0,
                scale_pos_weight=scale_pos_weight,
                objective="binary:logistic",
                eval_metric="aucpr",
                tree_method="hist",
                random_state=args.seed,
                n_jobs=1,
            ),
            x_train, y_train, x_val, y_val, x_test,
            sample_weight=train_sample_weight,
        )
        models.append(("phase_d_xgboost", vp, tp))
    if "phase_c_random_forest" in requested_models:
        vp, tp = _classical_prob(
            RandomForestClassifier(
                n_estimators=400,
                class_weight="balanced_subsample",
                min_samples_leaf=2,
                random_state=args.seed,
                n_jobs=-1,
            ),
            x_train, y_train, x_val, y_val, x_test,
            sample_weight=train_sample_weight,
        )
        models.append(("phase_c_random_forest", vp, tp))
    if "phase_c_svm_rbf" in requested_models:
        vp, tp = _classical_prob(
            SVC(
                kernel="rbf",
                probability=True,
                class_weight="balanced",
                gamma="scale",
                C=1.0,
                random_state=args.seed,
            ),
            x_train, y_train, x_val, y_val, x_test,
            sample_weight=train_sample_weight,
        )
        models.append(("phase_c_svm_rbf", vp, tp))
    if "phase_c_adaboost" in requested_models:
        # AdaBoost has no class_weight parameter; use sample_weight rebalance.
        ada_sw = np.ones(len(y_train), dtype=np.float32)
        pos = float(y_train.sum())
        neg = float(len(y_train) - pos)
        if neg > 0:
            ada_sw[y_train == 0] = pos / max(neg, 1.0)
        if train_sample_weight is not None:
            ada_sw = ada_sw * train_sample_weight
        vp, tp = _classical_prob(
            AdaBoostClassifier(
                n_estimators=200,
                learning_rate=0.5,
                random_state=args.seed,
            ),
            x_train, y_train, x_val, y_val, x_test,
            sample_weight=ada_sw,
        )
        models.append(("phase_c_adaboost", vp, tp))
    if "phase_c_knn" in requested_models:
        vp, tp = _classical_prob(
            KNeighborsClassifier(
                n_neighbors=25,
                weights="distance",
            ),
            x_train, y_train, x_val, y_val, x_test,
            sample_weight=None,
        )
        models.append(("phase_c_knn", vp, tp))

    if "phase_c_deepimmuno_cnn" in requested_models:
        from Daedalus.models import DeepImmunoStyleCNN

        cnn = DeepImmunoStyleCNN(
            n_numeric=meta["numeric_dim"],
            n_families=meta["family_map_size"],
            n_reference_sources=meta["refsrc_map_size"],
            n_change_classes=meta["change_map_size"],
        ).to(device)
        vp, tp = _train_deep_binary(
            cnn, train, val, test, device, args.max_epochs, args.patience,
            args.batch_size, sample_weight=train_sample_weight,
        )
        models.append(("phase_c_deepimmuno_cnn", vp, tp))
    if "phase_c_deepimmuno_cnn_autoencoder" in requested_models:
        from Daedalus.models import DeepImmunoStyleCNNAutoencoder

        cnn_ae = DeepImmunoStyleCNNAutoencoder(
            n_numeric=meta["numeric_dim"],
            n_families=meta["family_map_size"],
            n_reference_sources=meta["refsrc_map_size"],
            n_change_classes=meta["change_map_size"],
        ).to(device)
        vp, tp = _train_deep_binary(
            cnn_ae, train, val, test, device, args.max_epochs, args.patience,
            args.batch_size, recon_weight=0.10, sample_weight=train_sample_weight,
        )
        models.append(("phase_c_deepimmuno_cnn_autoencoder", vp, tp))
    if "phase_c_residual_mlp_autoencoder" in requested_models:
        from Daedalus.models import ResidualMLPAutoencoder

        res_ae = ResidualMLPAutoencoder(
            n_numeric=meta["numeric_dim"],
            n_families=meta["family_map_size"],
            n_reference_sources=meta["refsrc_map_size"],
            n_change_classes=meta["change_map_size"],
        ).to(device)
        vp, tp = _train_deep_binary(
            res_ae, train, val, test, device, args.max_epochs, args.patience,
            args.batch_size, recon_weight=0.10, sample_weight=train_sample_weight,
        )
        models.append(("phase_c_residual_mlp_autoencoder", vp, tp))

    # ResidualMLPAutoencoder variants -- explore hidden/latent/dropout/recon
    _variant_configs = [
        ("residual_mlp_ae_wide",       dict(hidden_dim=256, latent_dim=64, dropout=0.15), 0.10),
        ("residual_mlp_ae_deep",       dict(hidden_dim=192, latent_dim=48, dropout=0.15, n_blocks=3), 0.10),
        ("residual_mlp_ae_highdrop",   dict(hidden_dim=192, latent_dim=48, dropout=0.30), 0.10),
        ("residual_mlp_ae_strongrecon",dict(hidden_dim=192, latent_dim=48, dropout=0.15), 0.30),
        ("residual_mlp_ae_narrow",     dict(hidden_dim=128, latent_dim=32, dropout=0.15), 0.10),
    ]
    for v_name, v_cfg, v_recon in _variant_configs:
        if v_name not in requested_models:
            continue
        from Daedalus.models.residual_mlp_variants import ResidualMLPAutoencoderVariant
        model_v = ResidualMLPAutoencoderVariant(
            n_numeric=meta["numeric_dim"],
            n_families=meta["family_map_size"],
            n_reference_sources=meta["refsrc_map_size"],
            n_change_classes=meta["change_map_size"],
            **v_cfg,
        ).to(device)
        vp, tp = _train_deep_binary(
            model_v, train, val, test, device, args.max_epochs, args.patience,
            args.batch_size, recon_weight=v_recon, sample_weight=train_sample_weight,
        )
        models.append((v_name, vp, tp))

    if "phase_d_family_multitask" in requested_models:
        from Daedalus.models.family_multitask_net import FamilyMultitaskNet

        family_model = FamilyMultitaskNet(
            n_numeric=meta["numeric_dim"],
            n_tasks=1,
            n_families=meta["family_map_size"],
            n_reference_sources=meta["refsrc_map_size"],
            n_change_classes=meta["change_map_size"],
        ).to(device)
        vp, tp = _train_deep_binary(
            family_model, train, val, test, device, args.max_epochs, args.patience,
            args.batch_size, recon_weight=0.10, use_task_ids=True,
            sample_weight=train_sample_weight,
        )
        models.append(("phase_d_family_multitask", vp, tp))

    scopes = ["overall", NEG_GROUP_TM_RETAINED, NEG_GROUP_TM_LOST]
    rows: list[dict[str, object]] = []
    val_y = val.y
    for model_name, val_prob, test_prob in models:
        # Fixed-0.5 threshold metrics
        for scope in scopes:
            metric_row = _scope_metrics(test_frame, test_prob, scope, threshold=0.5)
            metric_row["model"] = model_name
            metric_row["threshold_mode"] = "fixed_0.5"
            rows.append(metric_row)
        # NPV-tuned threshold metrics
        if args.tune_threshold and len(val_y) > 0:
            tuned_thr = _tune_threshold_for_negatives(
                y_val=val_y,
                p_val=val_prob,
                sensitivity_floor=float(args.sensitivity_floor),
            )
            for scope in scopes:
                metric_row = _scope_metrics(test_frame, test_prob, scope, threshold=tuned_thr)
                metric_row["model"] = model_name
                metric_row["threshold_mode"] = f"tuned_npv_sens>={args.sensitivity_floor:.2f}"
                rows.append(metric_row)

    out = pd.DataFrame(rows)
    sort_cols = ["scope", "threshold_mode", "npv", "specificity", "balanced_accuracy"]
    sort_asc = [True, True, False, False, False]
    out = out.sort_values(
        [c for c in sort_cols if c in out.columns],
        ascending=[sort_asc[i] for i, c in enumerate(sort_cols) if c in out.columns],
    )
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)

    summary = {
        "benchmark": schema["benchmark_name"],
        "device": device.type,
        "train_rows": int(len(train.y)),
        "val_rows": int(len(val.y)),
        "test_rows": int(len(test.y)),
        "test_group_counts": test_frame["benchmark_group"].value_counts().to_dict(),
        "best_by_scope": {},
    }
    for scope in scopes:
        sub = out[out["scope"] == scope].copy()
        if "threshold_mode" in sub.columns:
            sub = sub.sort_values(
                ["threshold_mode", "npv", "specificity"], ascending=[True, False, False]
            )
        else:
            sub = sub.sort_values(
                ["balanced_accuracy", "specificity", "auroc"], ascending=False
            )
        if sub.empty:
            continue
        top = sub.iloc[0]
        summary["best_by_scope"][scope] = {
            "model": str(top["model"]),
            "threshold_mode": str(top.get("threshold_mode", "fixed_0.5")),
            "threshold": float(top.get("threshold", 0.5)),
            "balanced_accuracy": float(top["balanced_accuracy"]),
            "specificity": float(top["specificity"]),
            "sensitivity": float(top["sensitivity"]),
            "npv": float(top.get("npv", float("nan"))),
            "ppv": float(top.get("ppv", float("nan"))),
            "auroc": float(top["auroc"]),
            "average_precision": float(top["average_precision"]),
        }
    args.out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# TM Isoform Benchmark Model Comparison",
        "",
        f"Benchmark: `{schema['benchmark_name']}`",
        f"Device: `{device.type}`",
        f"Train rows: `{len(train.y)}`",
        f"Validation rows: `{len(val.y)}`",
        f"Test rows: `{len(test.y)}`",
        "",
        "Test composition:",
        f"- `{POS_GROUP}`: `{int((test_frame['benchmark_group'] == POS_GROUP).sum())}`",
        f"- `{NEG_GROUP_TM_RETAINED}`: `{int((test_frame['benchmark_group'] == NEG_GROUP_TM_RETAINED).sum())}`",
        f"- `{NEG_GROUP_TM_LOST}`: `{int((test_frame['benchmark_group'] == NEG_GROUP_TM_LOST).sum())}`",
        "",
    ]
    for scope in scopes:
        lines.append(f"## {scope}")
        lines.append("")
        sub = out[out["scope"] == scope].copy()
        if "threshold_mode" in sub.columns:
            sub = sub.sort_values(
                ["threshold_mode", "npv", "specificity"], ascending=[True, False, False]
            )
        else:
            sub = sub.sort_values(
                ["balanced_accuracy", "specificity", "auroc"], ascending=False
            )
        for row in sub.itertuples(index=False):
            tmode = getattr(row, "threshold_mode", "fixed_0.5")
            thr = getattr(row, "threshold", 0.5)
            npv_val = getattr(row, "npv", float("nan"))
            lines.append(
                f"- `{row.model}` [{tmode} thr={thr:.3f}] BalAcc={row.balanced_accuracy:.4f} "
                f"Spec={row.specificity:.4f} Sens={row.sensitivity:.4f} "
                f"NPV={npv_val:.4f} AUROC={row.auroc:.4f} AP={row.average_precision:.4f}"
            )
        lines.append("")
    args.out_md.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")

    print(f"[ok] Wrote {args.out_tsv}")
    print(f"[ok] Wrote {args.out_md}")
    print(f"[ok] Wrote {args.out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
