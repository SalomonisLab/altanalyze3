"""
Evaluation metrics for Proteus.

Includes AUROC, AUPRC, Brier score, balanced accuracy,
and within-gene ranking accuracy.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd


def compute_task_metrics(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    task_name: str = "task",
) -> Dict[str, Any]:
    """
    Compute standard binary classification metrics for a single task.

    Handles degenerate cases (all-same label, empty inputs) gracefully.

    Parameters
    ----------
    y_true : np.ndarray
        Binary labels (0/1). Labels of -1 are excluded.
    y_prob : np.ndarray
        Predicted probabilities in [0, 1].
    task_name : str
        Task name (for logging only).

    Returns
    -------
    dict
        Keys: auroc, auprc, brier, balanced_accuracy, n, positives, negatives.
        Values are float or None if metric cannot be computed.
    """
    from scipy.stats import rankdata

    # Filter out unlabeled (-1) samples
    valid = y_true != -1
    y_true = y_true[valid]
    y_prob = y_prob[valid]

    n = len(y_true)
    n_pos = int(y_true.sum())
    n_neg = n - n_pos

    metrics: Dict[str, Any] = {
        "task": task_name,
        "n": n,
        "positives": n_pos,
        "negatives": n_neg,
        "auroc": None,
        "auprc": None,
        "brier": None,
        "balanced_accuracy": None,
    }

    if n == 0:
        return metrics

    # Brier score (works with any label distribution)
    metrics["brier"] = float(np.mean((y_prob - y_true) ** 2))

    if n_pos == 0 or n_neg == 0:
        # Degenerate: can't compute AUROC/AUPRC
        return metrics

    try:
        from sklearn.metrics import (
            roc_auc_score,
            average_precision_score,
            balanced_accuracy_score,
        )

        metrics["auroc"] = float(roc_auc_score(y_true, y_prob))
        metrics["auprc"] = float(average_precision_score(y_true, y_prob))

        y_pred = (y_prob >= 0.5).astype(int)
        metrics["balanced_accuracy"] = float(balanced_accuracy_score(y_true, y_pred))

    except Exception as e:
        metrics["_error"] = str(e)

    return metrics


def within_gene_ranking_accuracy(
    pair_df: pd.DataFrame,
    probs: np.ndarray,
    gene_col: str = "gene_id",
    label_col: str = "weak_label",
    source_col: str = "reference_source",
) -> float:
    """
    Compute within-gene ranking accuracy.

    For each gene that has multiple isoform pairs, rank pairs by predicted
    global_preservation probability (descending). Measure what fraction of
    genes have the reference-preferred pair ranked highest.

    "Reference-preferred" means:
    - weak_label == "reference_preferred" OR "1" OR 1, OR
    - reference_source in ("MANE_SELECT", "APPRIS_PRINCIPAL_1", "MANE_CLINICAL")

    Parameters
    ----------
    pair_df : pd.DataFrame
        DataFrame with gene_id, weak_label, and optionally reference_source columns.
    probs : np.ndarray
        Predicted global_preservation probabilities, aligned with pair_df rows.
    gene_col : str
        Column name for gene identifier.
    label_col : str
        Column name for weak label.
    source_col : str
        Column name for reference source annotation.

    Returns
    -------
    float
        Fraction of multi-isoform genes where the reference-preferred pair
        is ranked highest. Returns 0.0 if no eligible genes.
    """
    if len(pair_df) == 0:
        return 0.0

    df = pair_df.copy()
    df["_prob"] = probs

    # Determine which rows are reference-preferred
    preferred_sources = {"MANE_SELECT", "APPRIS_PRINCIPAL_1", "MANE_CLINICAL", "MANE"}
    preferred_labels = {"reference_preferred", "1", "preserved"}

    def is_preferred(row) -> bool:
        label = str(row.get(label_col, "")).strip().lower()
        source = str(row.get(source_col, "")).strip()
        if label in preferred_labels or label == "1":
            return True
        if source in preferred_sources:
            return True
        try:
            return int(float(label)) == 1
        except (ValueError, TypeError):
            return False

    df["_is_preferred"] = df.apply(is_preferred, axis=1)

    genes_correct = 0
    genes_eligible = 0

    for gene_id, group in df.groupby(gene_col):
        if len(group) < 2:
            continue  # need at least 2 isoforms to rank

        preferred_rows = group[group["_is_preferred"]]
        if len(preferred_rows) == 0:
            continue  # no ground truth for this gene

        genes_eligible += 1

        # Rank by predicted probability (highest = rank 1)
        best_pred_idx = group["_prob"].idxmax()
        if best_pred_idx in preferred_rows.index:
            genes_correct += 1

    if genes_eligible == 0:
        return 0.0

    return float(genes_correct) / float(genes_eligible)


def compute_all_metrics(
    outputs: Dict[str, Any],
    batch: Dict[str, Any],
) -> Dict[str, Dict[str, Any]]:
    """
    Compute metrics for all active tasks in the batch.

    Parameters
    ----------
    outputs : dict
        Model outputs (logits) from TaskHeads.
    batch : dict
        Batch dict containing labels and task masks.

    Returns
    -------
    dict
        Maps task_name -> metrics dict from compute_task_metrics().
    """
    import torch

    from ..heads.task_heads import BINARY_TASKS

    label_preservation = batch.get("label_preservation")
    task_masks = batch.get("task_masks", {})

    if label_preservation is None:
        return {}

    if hasattr(label_preservation, "numpy"):
        y_true_base = label_preservation.cpu().numpy()
    else:
        y_true_base = np.array(label_preservation)

    all_metrics: Dict[str, Dict[str, Any]] = {}

    for task in BINARY_TASKS:
        if task not in outputs:
            continue

        logits = outputs[task]
        if hasattr(logits, "detach"):
            probs = torch.sigmoid(logits).detach().cpu().numpy().squeeze(-1)
        else:
            probs = np.array(logits).squeeze(-1)

        # Apply task mask if available
        task_mask = task_masks.get(task)
        if task_mask is not None:
            if hasattr(task_mask, "numpy"):
                mask_np = task_mask.cpu().numpy().astype(bool)
            else:
                mask_np = np.array(task_mask).astype(bool)
            y_true_masked = y_true_base.copy()
            y_true_masked[~mask_np] = -1  # mark as unlabeled
        else:
            y_true_masked = y_true_base

        all_metrics[task] = compute_task_metrics(y_true_masked, probs, task_name=task)

    # Deviation class accuracy (if available)
    dev_labels = batch.get("deviation_class")
    if dev_labels is not None and "deviation_class" in outputs:
        import torch
        dev_logits = outputs["deviation_class"]
        if hasattr(dev_logits, "detach"):
            dev_probs = torch.softmax(dev_logits, dim=-1).detach().cpu().numpy()
            dev_preds = dev_probs.argmax(axis=-1)
        else:
            dev_probs = np.array(dev_logits)
            dev_preds = dev_probs.argmax(axis=-1)

        if hasattr(dev_labels, "numpy"):
            dev_true = dev_labels.cpu().numpy()
        else:
            dev_true = np.array(dev_labels)

        valid = dev_true != -1
        if valid.sum() > 0:
            from sklearn.metrics import balanced_accuracy_score
            try:
                bal_acc = balanced_accuracy_score(dev_true[valid], dev_preds[valid])
            except Exception:
                bal_acc = None
            all_metrics["deviation_class"] = {
                "task": "deviation_class",
                "n": int(valid.sum()),
                "balanced_accuracy": bal_acc,
            }

    return all_metrics
