#!/usr/bin/env python3
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    roc_auc_score,
)
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


PROCESSED = Path(__file__).resolve().parents[1] / "data" / "processed"
IN_PATH = PROCESSED / "baseline_feature_matrix.tsv"
METRICS_OUT = PROCESSED / "torch_seed_baseline_metrics.tsv"
REPORT_OUT = PROCESSED / "torch_seed_baseline_report.md"

TASKS = [
    "task_global_seed",
    "task_membrane_seed",
    "task_kinase_seed",
    "task_tf_seed",
]

FEATURE_COLUMNS = [
    "feature_protein_ratio",
    "feature_protein_ratio_distance",
    "feature_protein_delta_norm",
    "feature_transcript_delta_norm",
    "feature_alt_has_protein",
    "feature_ref_alt_same_protein_length",
    "feature_gene_has_appris_principal_1",
    "feature_gene_clinvar_splice",
    "feature_gene_hpa_extracellular",
    "feature_is_membrane_pair",
    "feature_is_surface_pair",
    "feature_is_kinase_pair",
    "feature_is_tf_pair",
    "feature_reference_source_mane",
    "feature_reference_source_appris",
    "feature_reference_source_uniprot_longest",
    "feature_reference_source_longest_protein",
]


def _device() -> torch.device:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


class TabularMLP(nn.Module):
    def __init__(self, input_dim: int) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(),
            nn.Dropout(0.15),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Dropout(0.10),
            nn.Linear(32, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(-1)


@dataclass
class SplitData:
    x: np.ndarray
    y: np.ndarray


def _safe_metric(metric_fn, y_true, y_score):
    try:
        return float(metric_fn(y_true, y_score))
    except Exception:
        return float("nan")


def _metrics_dict(model: str, task: str, split: str, y_true, y_prob, y_pred) -> dict[str, str]:
    return {
        "model": model,
        "task": task,
        "split": split,
        "n": str(len(y_true)),
        "positives": str(int(np.sum(y_true))),
        "negatives": str(int(len(y_true) - np.sum(y_true))),
        "auroc": f"{_safe_metric(roc_auc_score, y_true, y_prob):.6f}",
        "average_precision": f"{_safe_metric(average_precision_score, y_true, y_prob):.6f}",
        "brier": f"{_safe_metric(brier_score_loss, y_true, y_prob):.6f}",
        "accuracy": f"{accuracy_score(y_true, y_pred):.6f}",
        "balanced_accuracy": f"{balanced_accuracy_score(y_true, y_pred):.6f}",
    }


def _prepare_splits(task_df: pd.DataFrame) -> tuple[SplitData, SplitData, SplitData]:
    train_df = task_df[task_df["split"] == "train"].copy()
    val_df = task_df[task_df["split"] == "val"].copy()
    test_df = task_df[task_df["split"] == "test"].copy()
    train_x = train_df[FEATURE_COLUMNS].to_numpy(dtype=np.float32)
    val_x = val_df[FEATURE_COLUMNS].to_numpy(dtype=np.float32)
    test_x = test_df[FEATURE_COLUMNS].to_numpy(dtype=np.float32)
    mean = train_x.mean(axis=0, keepdims=True)
    std = train_x.std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    train_x = (train_x - mean) / std
    val_x = (val_x - mean) / std
    test_x = (test_x - mean) / std
    return (
        SplitData(train_x, train_df["label_preservation_seed"].to_numpy(dtype=np.float32)),
        SplitData(val_x, val_df["label_preservation_seed"].to_numpy(dtype=np.float32)),
        SplitData(test_x, test_df["label_preservation_seed"].to_numpy(dtype=np.float32)),
    )


def _predict(model: nn.Module, x: np.ndarray, device: torch.device) -> np.ndarray:
    model.eval()
    with torch.no_grad():
        logits = model(torch.tensor(x, dtype=torch.float32, device=device))
        probs = torch.sigmoid(logits).detach().cpu().numpy()
    return probs


def _train_model(train: SplitData, val: SplitData, device: torch.device) -> TabularMLP:
    model = TabularMLP(train.x.shape[1]).to(device)
    pos = float(train.y.sum())
    neg = float(len(train.y) - pos)
    pos_weight = torch.tensor([neg / max(pos, 1.0)], device=device, dtype=torch.float32)
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)

    loader = DataLoader(
        TensorDataset(
            torch.tensor(train.x, dtype=torch.float32),
            torch.tensor(train.y, dtype=torch.float32),
        ),
        batch_size=512,
        shuffle=True,
    )
    best_state = None
    best_ap = -1.0
    patience = 10
    stale = 0
    for _epoch in range(80):
        model.train()
        for xb, yb in loader:
            xb = xb.to(device)
            yb = yb.to(device)
            optimizer.zero_grad(set_to_none=True)
            loss = criterion(model(xb), yb)
            loss.backward()
            optimizer.step()

        val_prob = _predict(model, val.x, device)
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
    return model


def main() -> int:
    df = pd.read_csv(IN_PATH, sep="\t")
    device = _device()
    metric_rows: list[dict[str, str]] = []
    report_lines = [
        "# Torch Seed Baseline Report",
        "",
        f"Device: `{device.type}`",
        "",
        "Tabular MLP trained on structural-only features.",
        "",
    ]

    for task in TASKS:
        task_df = df[(df[task] == 1) & (df["label_preservation_seed"].isin([0, 1]))].copy()
        if task_df.empty:
            continue
        train, val, test = _prepare_splits(task_df)
        if len(train.y) == 0 or len(val.y) == 0 or len(test.y) == 0:
            continue
        model = _train_model(train, val, device)
        val_prob = _predict(model, val.x, device)
        test_prob = _predict(model, test.x, device)
        val_pred = (val_prob >= 0.5).astype(int)
        test_pred = (test_prob >= 0.5).astype(int)
        val_metrics = _metrics_dict("torch_mlp_structural", task, "val", val.y, val_prob, val_pred)
        test_metrics = _metrics_dict("torch_mlp_structural", task, "test", test.y, test_prob, test_pred)
        metric_rows.extend([val_metrics, test_metrics])
        report_lines.append(
            f"- `{task}` val AUROC={val_metrics['auroc']} AP={val_metrics['average_precision']} | "
            f"test AUROC={test_metrics['auroc']} AP={test_metrics['average_precision']} BalAcc={test_metrics['balanced_accuracy']}"
        )

    with METRICS_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "model",
                "task",
                "split",
                "n",
                "positives",
                "negatives",
                "auroc",
                "average_precision",
                "brier",
                "accuracy",
                "balanced_accuracy",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(metric_rows)

    with REPORT_OUT.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(report_lines) + "\n")

    print(f"[ok] Wrote {METRICS_OUT} rows={len(metric_rows)}")
    print(f"[ok] Wrote {REPORT_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
