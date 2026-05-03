#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from sklearn.decomposition import PCA
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.svm import SVC
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Daedalus.phase_c.models import (
    DeepImmunoStyleCNN,
    DeepImmunoStyleCNNAutoencoder,
    ResidualMLPAutoencoder,
    choose_device,
)


PHASE_C = Path(__file__).resolve().parents[1]
PROCESSED = PHASE_C / "data" / "processed"
CHECKPOINTS = PHASE_C / "checkpoints"
IN_PATH = PROCESSED / "biochem_multitask_instances.parquet"
SCHEMA_PATH = PROCESSED / "biochem_multitask_schema.json"
OUT_TSV = CHECKPOINTS / "deepimmuno_style_benchmark_metrics.tsv"
OUT_MD = CHECKPOINTS / "deepimmuno_style_benchmark_report.md"
OUT_JSON = CHECKPOINTS / "deepimmuno_style_benchmark_summary.json"


@dataclass
class SplitArrays:
    x_numeric: np.ndarray
    family_id: np.ndarray
    refsrc_id: np.ndarray
    change_id: np.ndarray
    y: np.ndarray


def _safe_metric(fn, y_true, y_score) -> float:
    try:
        return float(fn(y_true, y_score))
    except Exception:
        return float("nan")


def _metrics(y_true: np.ndarray, probs: np.ndarray) -> dict[str, float]:
    pred = (probs >= 0.5).astype(int)
    return {
        "n": int(len(y_true)),
        "positives": int(y_true.sum()),
        "negatives": int(len(y_true) - y_true.sum()),
        "average_precision": _safe_metric(average_precision_score, y_true, probs),
        "auroc": _safe_metric(roc_auc_score, y_true, probs),
        "balanced_accuracy": _safe_metric(balanced_accuracy_score, y_true, pred),
    }


def _group_map(numeric_columns: list[str]) -> dict[str, list[str]]:
    segment = {
        "common_prefix_len",
        "common_suffix_len",
        "ref_changed_len",
        "alt_changed_len",
        "ref_changed_frac",
        "alt_changed_frac",
        "is_n_term_change",
        "is_c_term_change",
        "is_internal_change",
        "is_complex_change",
        "is_truncation",
        "is_extension",
        "ref_change_start",
        "ref_change_end",
    }
    groups: dict[str, list[str]] = {
        "base_reference_delta": [],
        "segment_change": [],
        "full_sequence_physchem": [],
        "terminal_sequence_physchem": [],
        "motifs": [],
        "typed_overlap": [],
    }
    for col in numeric_columns:
        if col in segment:
            groups["segment_change"].append(col)
        elif col.startswith(("ref_full_", "alt_full_")):
            groups["full_sequence_physchem"].append(col)
        elif col.startswith(("ref_seq_nterm_", "ref_seq_cterm_", "alt_seq_nterm_", "alt_seq_cterm_")):
            groups["terminal_sequence_physchem"].append(col)
        elif "_motif_" in col:
            groups["motifs"].append(col)
        elif col.startswith("overlap_"):
            groups["typed_overlap"].append(col)
        else:
            groups["base_reference_delta"].append(col)
    return groups


def _encode_map(train_values: pd.Series, all_values: list[str]) -> dict[str, int]:
    categories = sorted(set(train_values.astype(str).tolist()) | set(all_values))
    return {v: i for i, v in enumerate(categories)}


def _fit_block_pca(train: np.ndarray, val: np.ndarray, test: np.ndarray, max_components: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, float]:
    n_features = train.shape[1]
    if n_features == 0:
        return train, val, test, 0, 0.0
    if n_features <= 4:
        return train, val, test, n_features, 1.0
    scaler = StandardScaler()
    train_s = scaler.fit_transform(train)
    val_s = scaler.transform(val)
    test_s = scaler.transform(test)
    pca_full = PCA().fit(train_s)
    cum = np.cumsum(pca_full.explained_variance_ratio_)
    n_comp = int(np.searchsorted(cum, 0.95) + 1)
    n_comp = max(2, min(max_components, n_comp, n_features))
    pca = PCA(n_components=n_comp)
    return pca.fit_transform(train_s), pca.transform(val_s), pca.transform(test_s), n_comp, float(cum[n_comp - 1])


def _build_task_splits(task_df: pd.DataFrame, numeric_columns: list[str]) -> tuple[dict[str, SplitArrays], dict[str, object]]:
    groups = _group_map(numeric_columns)
    train_df = task_df[task_df["split"] == "train"].copy()
    val_df = task_df[task_df["split"] == "val"].copy()
    test_df = task_df[task_df["split"] == "test"].copy()

    family_map = _encode_map(train_df["family_class"], task_df["family_class"].astype(str).unique().tolist())
    refsrc_map = _encode_map(train_df["reference_source"], task_df["reference_source"].astype(str).unique().tolist())
    change_map = _encode_map(train_df["change_class"], task_df["change_class"].astype(str).unique().tolist())

    pca_meta: dict[str, dict[str, object]] = {}
    reduced_parts: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    for name, cols in groups.items():
        if name in {"full_sequence_physchem", "terminal_sequence_physchem", "typed_overlap"}:
            tr, va, te, n_comp, explained = _fit_block_pca(
                train_df[cols].to_numpy(dtype=np.float32),
                val_df[cols].to_numpy(dtype=np.float32),
                test_df[cols].to_numpy(dtype=np.float32),
                max_components=24 if name != "typed_overlap" else 16,
            )
            reduced_parts[name] = (tr, va, te)
            pca_meta[name] = {
                "original_dim": len(cols),
                "reduced_dim": int(n_comp),
                "explained_variance": explained,
            }
        else:
            reduced_parts[name] = (
                train_df[cols].to_numpy(dtype=np.float32),
                val_df[cols].to_numpy(dtype=np.float32),
                test_df[cols].to_numpy(dtype=np.float32),
            )
            pca_meta[name] = {
                "original_dim": len(cols),
                "reduced_dim": len(cols),
                "explained_variance": 1.0,
            }

    def stack(idx: int) -> np.ndarray:
        mats = [reduced_parts[name][idx] for name in groups]
        return np.concatenate(mats, axis=1).astype(np.float32)

    x_train = stack(0)
    x_val = stack(1)
    x_test = stack(2)
    scaler = StandardScaler()
    x_train = scaler.fit_transform(x_train).astype(np.float32)
    x_val = scaler.transform(x_val).astype(np.float32)
    x_test = scaler.transform(x_test).astype(np.float32)

    def encode(frame: pd.DataFrame, mapping: dict[str, int], column: str) -> np.ndarray:
        return frame[column].astype(str).map(lambda x: mapping.get(x, 0)).to_numpy(dtype=np.int64)

    bundles = {
        "train": SplitArrays(
            x_numeric=x_train,
            family_id=encode(train_df, family_map, "family_class"),
            refsrc_id=encode(train_df, refsrc_map, "reference_source"),
            change_id=encode(train_df, change_map, "change_class"),
            y=train_df["label_preservation_seed"].to_numpy(dtype=np.int64),
        ),
        "val": SplitArrays(
            x_numeric=x_val,
            family_id=encode(val_df, family_map, "family_class"),
            refsrc_id=encode(val_df, refsrc_map, "reference_source"),
            change_id=encode(val_df, change_map, "change_class"),
            y=val_df["label_preservation_seed"].to_numpy(dtype=np.int64),
        ),
        "test": SplitArrays(
            x_numeric=x_test,
            family_id=encode(test_df, family_map, "family_class"),
            refsrc_id=encode(test_df, refsrc_map, "reference_source"),
            change_id=encode(test_df, change_map, "change_class"),
            y=test_df["label_preservation_seed"].to_numpy(dtype=np.int64),
        ),
    }
    meta = {
        "pca": pca_meta,
        "family_map_size": len(family_map),
        "refsrc_map_size": len(refsrc_map),
        "change_map_size": len(change_map),
        "reduced_numeric_dim": int(x_train.shape[1]),
        "train_n": int(len(train_df)),
        "val_n": int(len(val_df)),
        "test_n": int(len(test_df)),
    }
    return bundles, meta


def _classical_matrix(bundle: SplitArrays, n_family: int, n_refsrc: int, n_change: int) -> np.ndarray:
    family = np.eye(max(n_family, 1), dtype=np.float32)[bundle.family_id]
    refsrc = np.eye(max(n_refsrc, 1), dtype=np.float32)[bundle.refsrc_id]
    change = np.eye(max(n_change, 1), dtype=np.float32)[bundle.change_id]
    return np.concatenate([bundle.x_numeric, family, refsrc, change], axis=1)


def _limit_train(bundle: SplitArrays, max_n: int) -> SplitArrays:
    if max_n <= 0 or len(bundle.y) <= max_n:
        return bundle
    rng = np.random.default_rng(0)
    idx = rng.choice(len(bundle.y), size=max_n, replace=False)
    idx.sort()
    return SplitArrays(
        x_numeric=bundle.x_numeric[idx],
        family_id=bundle.family_id[idx],
        refsrc_id=bundle.refsrc_id[idx],
        change_id=bundle.change_id[idx],
        y=bundle.y[idx],
    )


def _run_classical_models(task: str, splits: dict[str, SplitArrays]) -> list[dict[str, object]]:
    train = splits["train"]
    val = splits["val"]
    test = splits["test"]

    n_family = int(max(train.family_id.max(), val.family_id.max(), test.family_id.max()) + 1)
    n_refsrc = int(max(train.refsrc_id.max(), val.refsrc_id.max(), test.refsrc_id.max()) + 1)
    n_change = int(max(train.change_id.max(), val.change_id.max(), test.change_id.max()) + 1)
    x_train = _classical_matrix(train, n_family, n_refsrc, n_change)
    x_val = _classical_matrix(val, n_family, n_refsrc, n_change)
    x_test = _classical_matrix(test, n_family, n_refsrc, n_change)
    y_train, y_val, y_test = train.y, val.y, test.y

    rows: list[dict[str, object]] = []

    print(f"[task={task}] fitting elasticnet", flush=True)
    elastic = SGDClassifier(
        loss="log_loss",
        penalty="elasticnet",
        alpha=1e-4,
        l1_ratio=0.15,
        max_iter=2000,
        class_weight="balanced",
        random_state=0,
    )
    elastic.fit(x_train, y_train)
    probs = elastic.predict_proba(x_test)[:, 1]
    rows.append({"model": "elasticnet", "task": task, "feature_mode": "pca_reduced", **_metrics(y_test, probs)})

    print(f"[task={task}] fitting knn", flush=True)
    knn_train = _limit_train(train, 10000)
    knn = make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=25, weights="distance"))
    knn.fit(knn_train.x_numeric, knn_train.y)
    probs = knn.predict_proba(test.x_numeric)[:, 1]
    rows.append({"model": "knn", "task": task, "feature_mode": "pca_reduced", **_metrics(y_test, probs)})

    print(f"[task={task}] fitting svm_rbf", flush=True)
    svm_train = _limit_train(train, 6000)
    svm = make_pipeline(
        StandardScaler(),
        SVC(C=2.0, kernel="rbf", gamma="scale", probability=True, class_weight="balanced", random_state=0),
    )
    svm.fit(svm_train.x_numeric, svm_train.y)
    probs = svm.predict_proba(test.x_numeric)[:, 1]
    rows.append({"model": "svm_rbf", "task": task, "feature_mode": "pca_reduced", **_metrics(y_test, probs)})

    print(f"[task={task}] fitting random_forest", flush=True)
    rf_train = _limit_train(train, 15000)
    x_rf = _classical_matrix(rf_train, n_family, n_refsrc, n_change)
    rf = RandomForestClassifier(
        n_estimators=250,
        max_depth=None,
        min_samples_leaf=2,
        n_jobs=8,
        class_weight="balanced_subsample",
        random_state=0,
    )
    rf.fit(x_rf, rf_train.y)
    probs = rf.predict_proba(x_test)[:, 1]
    rows.append({"model": "random_forest", "task": task, "feature_mode": "pca_reduced", **_metrics(y_test, probs)})

    print(f"[task={task}] fitting adaboost", flush=True)
    ada_train = _limit_train(train, 12000)
    x_ada = _classical_matrix(ada_train, n_family, n_refsrc, n_change)
    ada = AdaBoostClassifier(n_estimators=200, learning_rate=0.05, random_state=0)
    ada.fit(x_ada, ada_train.y)
    probs = ada.predict_proba(x_test)[:, 1]
    rows.append({"model": "adaboost", "task": task, "feature_mode": "pca_reduced", **_metrics(y_test, probs)})

    print(f"[task={task}] fitting hist_gradient_boosting", flush=True)
    hgb_train = _limit_train(train, 18000)
    x_hgb = _classical_matrix(hgb_train, n_family, n_refsrc, n_change)
    hgb = HistGradientBoostingClassifier(
        learning_rate=0.05,
        max_depth=6,
        max_iter=300,
        early_stopping=True,
        validation_fraction=0.1,
        random_state=0,
    )
    hgb.fit(x_hgb, hgb_train.y)
    probs = hgb.predict_proba(x_test)[:, 1]
    rows.append({"model": "hist_gradient_boosting", "task": task, "feature_mode": "pca_reduced", **_metrics(y_test, probs)})

    return rows


def _predict_torch(model: nn.Module, bundle: SplitArrays, device: torch.device) -> np.ndarray:
    model.eval()
    with torch.no_grad():
        numeric = torch.tensor(bundle.x_numeric, dtype=torch.float32, device=device)
        family = torch.tensor(bundle.family_id, dtype=torch.long, device=device)
        refsrc = torch.tensor(bundle.refsrc_id, dtype=torch.long, device=device)
        change = torch.tensor(bundle.change_id, dtype=torch.long, device=device)
        outputs = model(numeric, family, refsrc, change)
        logits = outputs[0]
        probs = torch.sigmoid(logits).detach().cpu().numpy()
    return probs


def _fit_torch_model(
    model_name: str,
    model: nn.Module,
    train: SplitArrays,
    val: SplitArrays,
    device: torch.device,
    max_epochs: int = 12,
    patience: int = 3,
    batch_size: int = 512,
) -> nn.Module:
    model = model.to(device)
    train_ds = TensorDataset(
        torch.tensor(train.x_numeric, dtype=torch.float32),
        torch.tensor(train.family_id, dtype=torch.long),
        torch.tensor(train.refsrc_id, dtype=torch.long),
        torch.tensor(train.change_id, dtype=torch.long),
        torch.tensor(train.y, dtype=torch.float32),
    )
    loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
    pos = float(train.y.sum())
    neg = float(len(train.y) - pos)
    pos_weight = torch.tensor([neg / max(pos, 1.0)], dtype=torch.float32, device=device)
    bce = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    mse = nn.MSELoss()
    optimizer = torch.optim.AdamW(model.parameters(), lr=8e-4, weight_decay=2e-4)
    best_state = None
    best_ap = -1.0
    stale = 0
    has_recon = "Autoencoder" in model_name

    for _ in range(max_epochs):
        model.train()
        for numeric, family, refsrc, change, labels in loader:
            numeric = numeric.to(device)
            family = family.to(device)
            refsrc = refsrc.to(device)
            change = change.to(device)
            labels = labels.to(device)
            optimizer.zero_grad(set_to_none=True)
            outputs = model(numeric, family, refsrc, change)
            logits = outputs[0]
            loss = bce(logits, labels)
            if has_recon:
                recon = outputs[1]
                loss = loss + 0.10 * mse(recon, numeric)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=2.0)
            optimizer.step()
        val_probs = _predict_torch(model, val, device)
        val_ap = average_precision_score(val.y, val_probs)
        if val_ap > best_ap:
            best_ap = float(val_ap)
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
            if stale >= patience:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model.cpu()


def _run_torch_models(
    task: str,
    splits: dict[str, SplitArrays],
    meta: dict[str, object],
    selected_models: set[str] | None = None,
    max_epochs: int = 12,
    patience: int = 3,
) -> list[dict[str, object]]:
    device = choose_device("auto")
    train, val, test = splits["train"], splits["val"], splits["test"]
    rows: list[dict[str, object]] = []

    models = {
        "residual_mlp_autoencoder": ResidualMLPAutoencoder(
            n_numeric=train.x_numeric.shape[1],
            n_families=int(meta["family_map_size"]),
            n_reference_sources=int(meta["refsrc_map_size"]),
            n_change_classes=int(meta["change_map_size"]),
        ),
        "deepimmuno_cnn": DeepImmunoStyleCNN(
            n_numeric=train.x_numeric.shape[1],
            n_families=int(meta["family_map_size"]),
            n_reference_sources=int(meta["refsrc_map_size"]),
            n_change_classes=int(meta["change_map_size"]),
        ),
        "deepimmuno_cnn_autoencoder": DeepImmunoStyleCNNAutoencoder(
            n_numeric=train.x_numeric.shape[1],
            n_families=int(meta["family_map_size"]),
            n_reference_sources=int(meta["refsrc_map_size"]),
            n_change_classes=int(meta["change_map_size"]),
        ),
    }
    for name, model in models.items():
        if selected_models is not None and name not in selected_models:
            continue
        print(f"[task={task}] fitting {name}", flush=True)
        fitted = _fit_torch_model(name, model, train, val, device, max_epochs=max_epochs, patience=patience)
        probs = _predict_torch(fitted, test, torch.device("cpu"))
        rows.append({"model": name, "task": task, "feature_mode": "pca_reduced", **_metrics(test.y, probs)})
    return rows


def _write_outputs(rows: list[dict[str, object]], pca_summary: dict[str, object], tested_models: list[str]) -> None:
    out = pd.DataFrame(rows).sort_values(["task", "average_precision"], ascending=[True, False]).reset_index(drop=True)
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)
    summary = {
        "tested_models": tested_models,
        "pca_summary": pca_summary,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    lines = [
        "# Phase C DeepImmuno-Style Benchmark",
        "",
        "Implemented to mirror the DeepImmuno logic for this problem:",
        "- systematic ML/DL model comparison",
        "- PCA reduction on the high-dimensional physicochemical blocks",
        "- explicit comparison of deep learning versus classical methods",
        "- autoencoder + classifier joint training as one deep-learning option",
        "",
        "Tested models:",
    ]
    for model in tested_models:
        lines.append(f"- `{model}`")
    lines.extend(["", "## Test metrics by task", ""])
    for task in sorted(out["task"].astype(str).unique().tolist()):
        lines.append(f"### {task}")
        lines.append("")
        sub = out[out["task"] == task].copy().sort_values("average_precision", ascending=False)
        for row in sub.itertuples(index=False):
            lines.append(
                f"- `{row.model}` AP={row.average_precision:.6f} "
                f"AUROC={row.auroc:.6f} BalAcc={row.balanced_accuracy:.6f}"
            )
        meta = pca_summary[task]
        lines.append("")
        lines.append("PCA reduction:")
        for block_name, block_meta in meta["pca"].items():
            lines.append(
                f"- `{block_name}` {block_meta['original_dim']} -> {block_meta['reduced_dim']} "
                f"(explained={block_meta['explained_variance']:.4f})"
            )
        lines.append("")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tasks", nargs="*", default=None)
    parser.add_argument("--models", nargs="*", default=None)
    parser.add_argument("--max-epochs", type=int, default=12)
    parser.add_argument("--patience", type=int, default=3)
    args = parser.parse_args()

    df = pd.read_parquet(IN_PATH)
    schema = json.loads(SCHEMA_PATH.read_text())
    numeric_columns = schema["numeric_columns"]
    tasks = sorted(df["task"].astype(str).unique().tolist())
    if args.tasks:
        keep = set(args.tasks)
        tasks = [task for task in tasks if task in keep]

    rows: list[dict[str, object]] = []
    pca_summary: dict[str, object] = {}
    tested_models = [
        "elasticnet",
        "knn",
        "svm_rbf",
        "random_forest",
        "adaboost",
        "hist_gradient_boosting",
        "residual_mlp_autoencoder",
        "deepimmuno_cnn",
        "deepimmuno_cnn_autoencoder",
    ]
    selected_models = set(args.models) if args.models else None
    classical_models = {"elasticnet", "knn", "svm_rbf", "random_forest", "adaboost", "hist_gradient_boosting"}
    torch_models = {"residual_mlp_autoencoder", "deepimmuno_cnn", "deepimmuno_cnn_autoencoder"}

    for task in tasks:
        print(f"[benchmark] start task={task}", flush=True)
        task_df = df[df["task"] == task].copy()
        splits, meta = _build_task_splits(task_df, numeric_columns)
        pca_summary[task] = meta
        if selected_models is None or selected_models & classical_models:
            rows.extend(_run_classical_models(task, splits))
        if selected_models is None or selected_models & torch_models:
            rows.extend(
                _run_torch_models(
                    task,
                    splits,
                    meta,
                    selected_models=selected_models,
                    max_epochs=args.max_epochs,
                    patience=args.patience,
                )
            )
        output_models = tested_models if selected_models is None else [m for m in tested_models if m in selected_models]
        _write_outputs(rows, pca_summary, output_models)
    print(f"[ok] Wrote {OUT_TSV}")
    print(f"[ok] Wrote {OUT_MD}")
    print(f"[ok] Wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
