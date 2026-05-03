#!/usr/bin/env python3
from __future__ import annotations

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

from Daedalus.phase_b.models import ReferenceDeltaNet, choose_device


PROCESSED = Path(__file__).resolve().parents[1] / "data" / "processed"
CHECKPOINTS = Path(__file__).resolve().parents[1] / "checkpoints"
IN_PATH = PROCESSED / "multitask_instances.parquet"
SCHEMA_PATH = PROCESSED / "multitask_schema.json"
REPORT_PATH = CHECKPOINTS / "reference_delta_multitask_report.md"
METRICS_PATH = CHECKPOINTS / "reference_delta_multitask_metrics.json"
CKPT_PATH = CHECKPOINTS / "reference_delta_multitask.pt"


@dataclass
class SplitBundle:
    numeric: np.ndarray
    labels: np.ndarray
    task_ids: np.ndarray
    species_ids: np.ndarray
    family_ids: np.ndarray
    refsrc_ids: np.ndarray
    tasks: np.ndarray


def _safe_metric(metric_fn, y_true, y_score) -> float:
    try:
        return float(metric_fn(y_true, y_score))
    except Exception:
        return float("nan")


def _encode_categories(series: pd.Series) -> tuple[np.ndarray, dict[str, int]]:
    categories = sorted(series.astype(str).unique().tolist())
    mapping = {value: idx for idx, value in enumerate(categories)}
    encoded = series.astype(str).map(mapping).to_numpy(dtype=np.int64)
    return encoded, mapping


def _prepare() -> tuple[SplitBundle, SplitBundle, SplitBundle, list[str], dict]:
    df = pd.read_parquet(IN_PATH)
    with SCHEMA_PATH.open() as handle:
        schema = json.load(handle)
    numeric_columns = schema["numeric_columns"]
    train_df = df[df["split"] == "train"].copy()
    val_df = df[df["split"] == "val"].copy()
    test_df = df[df["split"] == "test"].copy()

    task_train, task_map = _encode_categories(train_df["task"])
    species_train, species_map = _encode_categories(train_df["species"])
    family_train, family_map = _encode_categories(train_df["family_class"])
    refsrc_train, refsrc_map = _encode_categories(train_df["reference_source"])

    def encode_with_map(series: pd.Series, mapping: dict[str, int]) -> np.ndarray:
        return series.astype(str).map(lambda x: mapping.get(x, 0)).to_numpy(dtype=np.int64)

    mean = train_df[numeric_columns].to_numpy(dtype=np.float32).mean(axis=0, keepdims=True)
    std = train_df[numeric_columns].to_numpy(dtype=np.float32).std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0

    def bundle(frame: pd.DataFrame, task_ids, species_ids, family_ids, refsrc_ids) -> SplitBundle:
        numeric = frame[numeric_columns].to_numpy(dtype=np.float32)
        numeric = (numeric - mean) / std
        labels = frame["label_preservation_seed"].to_numpy(dtype=np.float32)
        return SplitBundle(
            numeric=numeric,
            labels=labels,
            task_ids=task_ids,
            species_ids=species_ids,
            family_ids=family_ids,
            refsrc_ids=refsrc_ids,
            tasks=frame["task"].astype(str).to_numpy(),
        )

    train = bundle(train_df, task_train, species_train, family_train, refsrc_train)
    val = bundle(
        val_df,
        encode_with_map(val_df["task"], task_map),
        encode_with_map(val_df["species"], species_map),
        encode_with_map(val_df["family_class"], family_map),
        encode_with_map(val_df["reference_source"], refsrc_map),
    )
    test = bundle(
        test_df,
        encode_with_map(test_df["task"], task_map),
        encode_with_map(test_df["species"], species_map),
        encode_with_map(test_df["family_class"], family_map),
        encode_with_map(test_df["reference_source"], refsrc_map),
    )
    metadata = {
        "task_map": task_map,
        "species_map": species_map,
        "family_map": family_map,
        "reference_source_map": refsrc_map,
        "numeric_columns": numeric_columns,
        "mean": mean.astype(np.float32),
        "std": std.astype(np.float32),
    }
    return train, val, test, numeric_columns, metadata


def _predict(model: nn.Module, bundle: SplitBundle, device: torch.device) -> np.ndarray:
    model.eval()
    with torch.no_grad():
        logits = model(
            torch.tensor(bundle.numeric, dtype=torch.float32, device=device),
            torch.tensor(bundle.task_ids, dtype=torch.long, device=device),
            torch.tensor(bundle.species_ids, dtype=torch.long, device=device),
            torch.tensor(bundle.family_ids, dtype=torch.long, device=device),
            torch.tensor(bundle.refsrc_ids, dtype=torch.long, device=device),
        )
        probs = torch.sigmoid(logits).detach().cpu().numpy()
    return probs


def _task_metrics(tasks: np.ndarray, labels: np.ndarray, probs: np.ndarray) -> dict[str, dict[str, float]]:
    result: dict[str, dict[str, float]] = {}
    preds = (probs >= 0.5).astype(int)
    for task in sorted(set(tasks.tolist())):
        mask = tasks == task
        y = labels[mask]
        p = probs[mask]
        yhat = preds[mask]
        result[task] = {
            "n": int(mask.sum()),
            "positives": int(y.sum()),
            "negatives": int(mask.sum() - y.sum()),
            "average_precision": _safe_metric(average_precision_score, y, p),
            "auroc": _safe_metric(roc_auc_score, y, p),
            "balanced_accuracy": _safe_metric(balanced_accuracy_score, y, yhat),
        }
    return result


def _macro_ap(task_metrics: dict[str, dict[str, float]]) -> float:
    values = [m["average_precision"] for m in task_metrics.values() if not np.isnan(m["average_precision"])]
    return float(np.mean(values)) if values else float("nan")


def main() -> int:
    device = choose_device("auto")
    train, val, test, numeric_columns, metadata = _prepare()

    model = ReferenceDeltaNet(
        n_numeric=len(numeric_columns),
        n_tasks=len(metadata["task_map"]),
        n_species=len(metadata["species_map"]),
        n_families=len(metadata["family_map"]),
        n_reference_sources=len(metadata["reference_source_map"]),
    ).to(device)

    pos = float(train.labels.sum())
    neg = float(len(train.labels) - pos)
    pos_weight = torch.tensor([neg / max(pos, 1.0)], device=device, dtype=torch.float32)
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    dataset = TensorDataset(
        torch.tensor(train.numeric, dtype=torch.float32),
        torch.tensor(train.task_ids, dtype=torch.long),
        torch.tensor(train.species_ids, dtype=torch.long),
        torch.tensor(train.family_ids, dtype=torch.long),
        torch.tensor(train.refsrc_ids, dtype=torch.long),
        torch.tensor(train.labels, dtype=torch.float32),
    )
    loader = DataLoader(dataset, batch_size=1024, shuffle=True)

    best_state = None
    best = {"epoch": -1, "val_macro_ap": -1.0, "val_metrics": {}}
    stale = 0
    patience = 8
    for epoch in range(40):
        model.train()
        for numeric, task_ids, species_ids, family_ids, refsrc_ids, labels in loader:
            numeric = numeric.to(device)
            task_ids = task_ids.to(device)
            species_ids = species_ids.to(device)
            family_ids = family_ids.to(device)
            refsrc_ids = refsrc_ids.to(device)
            labels = labels.to(device)
            optimizer.zero_grad(set_to_none=True)
            logits = model(numeric, task_ids, species_ids, family_ids, refsrc_ids)
            loss = criterion(logits, labels)
            loss.backward()
            optimizer.step()

        val_probs = _predict(model, val, device)
        val_metrics = _task_metrics(val.tasks, val.labels, val_probs)
        val_macro = _macro_ap(val_metrics)
        if val_macro > best["val_macro_ap"]:
            best = {"epoch": epoch, "val_macro_ap": val_macro, "val_metrics": val_metrics}
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
            if stale >= patience:
                break

    if best_state is not None:
        model.load_state_dict(best_state)

    val_probs = _predict(model, val, device)
    test_probs = _predict(model, test, device)
    final = {
        "device": device.type,
        "best_epoch": best["epoch"],
        "val_macro_ap": _macro_ap(_task_metrics(val.tasks, val.labels, val_probs)),
        "test_macro_ap": _macro_ap(_task_metrics(test.tasks, test.labels, test_probs)),
        "val": _task_metrics(val.tasks, val.labels, val_probs),
        "test": _task_metrics(test.tasks, test.labels, test_probs),
    }

    CHECKPOINTS.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "state_dict": model.state_dict(),
            "metadata": metadata,
            "metrics": final,
        },
        CKPT_PATH,
    )
    with METRICS_PATH.open("w", encoding="utf-8") as handle:
        json.dump(final, handle, indent=2)
    lines = [
        "# Reference Delta Multitask Report",
        "",
        f"Device: `{final['device']}`",
        f"Best epoch: `{final['best_epoch']}`",
        f"Validation macro AP: `{final['val_macro_ap']:.6f}`",
        f"Test macro AP: `{final['test_macro_ap']:.6f}`",
        "",
        "## Test metrics",
        "",
    ]
    for task, metrics in final["test"].items():
        lines.append(
            f"- `{task}` n={metrics['n']} AP={metrics['average_precision']:.6f} "
            f"AUROC={metrics['auroc']:.6f} BalAcc={metrics['balanced_accuracy']:.6f}"
        )
    with REPORT_PATH.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")

    print(json.dumps({"checkpoint": str(CKPT_PATH), "metrics_path": str(METRICS_PATH), "report_path": str(REPORT_PATH), "summary": final}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
