#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Daedalus.phase_a.models import StructuralMLP, choose_device
PROCESSED = Path(__file__).resolve().parents[1] / "data" / "processed"
CHECKPOINTS = Path(__file__).resolve().parents[1] / "checkpoints"
IN_PATH = PROCESSED / "baseline_feature_matrix.tsv"

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


@dataclass
class SplitData:
    x: np.ndarray
    y: np.ndarray


def _safe_metric(metric_fn, y_true, y_score) -> float:
    try:
        return float(metric_fn(y_true, y_score))
    except Exception:
        return float("nan")


def _prepare_splits(task: str) -> tuple[SplitData, SplitData, SplitData, np.ndarray, np.ndarray]:
    df = pd.read_csv(IN_PATH, sep="\t")
    df = df[(df[task] == 1) & (df["label_preservation_seed"].isin([0, 1]))].copy()
    train_df = df[df["split"] == "train"]
    val_df = df[df["split"] == "val"]
    test_df = df[df["split"] == "test"]
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
        mean,
        std,
    )


def _predict(model: nn.Module, x: np.ndarray, device: torch.device) -> np.ndarray:
    model.eval()
    with torch.no_grad():
        logits = model(torch.tensor(x, dtype=torch.float32, device=device))
        probs = torch.sigmoid(logits).detach().cpu().numpy()
    return probs


def _fit(
    train: SplitData,
    val: SplitData,
    device: torch.device,
    epochs: int,
    batch_size: int,
    lr: float,
    weight_decay: float,
    patience: int,
) -> tuple[StructuralMLP, dict]:
    model = StructuralMLP(train.x.shape[1]).to(device)
    pos = float(train.y.sum())
    neg = float(len(train.y) - pos)
    pos_weight = torch.tensor([neg / max(pos, 1.0)], device=device, dtype=torch.float32)
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    loader = DataLoader(
        TensorDataset(
            torch.tensor(train.x, dtype=torch.float32),
            torch.tensor(train.y, dtype=torch.float32),
        ),
        batch_size=batch_size,
        shuffle=True,
    )
    best_state = None
    best = {"epoch": -1, "val_ap": -1.0}
    stale = 0
    for epoch in range(epochs):
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
        if val_ap > best["val_ap"]:
            best = {"epoch": epoch, "val_ap": val_ap}
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
            if stale >= patience:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, best


def main() -> int:
    parser = argparse.ArgumentParser(description="Train a structural MLP on one Daedalus seed task.")
    parser.add_argument("--task", required=True)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--patience", type=int, default=10)
    parser.add_argument("--output-prefix", default=None)
    args = parser.parse_args()

    device = choose_device(args.device)
    train, val, test, mean, std = _prepare_splits(args.task)
    model, best = _fit(
        train=train,
        val=val,
        device=device,
        epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        weight_decay=args.weight_decay,
        patience=args.patience,
    )

    val_prob = _predict(model, val.x, device)
    test_prob = _predict(model, test.x, device)
    val_pred = (val_prob >= 0.5).astype(int)
    test_pred = (test_prob >= 0.5).astype(int)

    metrics = {
        "task": args.task,
        "device": device.type,
        "best_epoch": best["epoch"],
        "val_ap": _safe_metric(average_precision_score, val.y, val_prob),
        "val_auroc": _safe_metric(roc_auc_score, val.y, val_prob),
        "val_balanced_accuracy": _safe_metric(balanced_accuracy_score, val.y, val_pred),
        "test_ap": _safe_metric(average_precision_score, test.y, test_prob),
        "test_auroc": _safe_metric(roc_auc_score, test.y, test_prob),
        "test_balanced_accuracy": _safe_metric(balanced_accuracy_score, test.y, test_pred),
        "n_train": int(len(train.y)),
        "n_val": int(len(val.y)),
        "n_test": int(len(test.y)),
    }

    prefix = args.output_prefix or args.task
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)
    ckpt_path = CHECKPOINTS / f"{prefix}.pt"
    metrics_path = CHECKPOINTS / f"{prefix}.metrics.json"
    torch.save(
        {
            "task": args.task,
            "feature_columns": FEATURE_COLUMNS,
            "state_dict": model.state_dict(),
            "mean": mean.astype(np.float32),
            "std": std.astype(np.float32),
            "metrics": metrics,
        },
        ckpt_path,
    )
    with metrics_path.open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)

    print(json.dumps({"checkpoint": str(ckpt_path), "metrics": metrics}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
