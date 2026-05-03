"""
Loss functions for Proteus supervised training and SSL pretraining.
"""

from __future__ import annotations

from typing import Dict, Optional

import torch
import torch.nn as nn
from torch import Tensor

BINARY_TASKS = [
    "global_preservation",
    "topology_preserved",
    "surface_retained",
    "kinase_competent",
    "tf_competent",
    "disorder_preserved",
]


class ProteusLoss(nn.Module):
    """
    Multi-task loss for Proteus supervised training.

    Computes masked BCE loss for each binary task and cross-entropy for
    the deviation_class multiclass task. Masking ensures that task-specific
    losses (e.g., kinase_competent) are only computed on relevant samples.

    Labels of -1 indicate unlabeled samples and are excluded from all losses.

    Parameters
    ----------
    task_weights : dict
        Maps task name -> float loss weight.
    pos_weights : dict | None
        Maps task name -> positive class weight (for class imbalance correction).
        Passed as pos_weight to BCEWithLogitsLoss.
    """

    def __init__(
        self,
        task_weights: Optional[Dict[str, float]] = None,
        pos_weights: Optional[Dict[str, float]] = None,
    ) -> None:
        super().__init__()

        default_weights = {
            "global_preservation": 1.0,
            "topology_preserved": 0.8,
            "surface_retained": 0.6,
            "kinase_competent": 0.8,
            "tf_competent": 0.8,
            "disorder_preserved": 0.6,
            "deviation_class": 0.5,
        }
        self.task_weights = task_weights or default_weights

        # Per-task BCE losses with optional pos_weight
        self.binary_losses: nn.ModuleDict = nn.ModuleDict()
        for task in BINARY_TASKS:
            pw = None
            if pos_weights and task in pos_weights:
                pw = torch.tensor([pos_weights[task]])
            self.binary_losses[task] = nn.BCEWithLogitsLoss(
                pos_weight=pw, reduction="none"
            )

        # Cross-entropy for deviation class (ignore_index=-1 for unlabeled)
        self.deviation_loss = nn.CrossEntropyLoss(ignore_index=-1, reduction="mean")

    def forward(
        self,
        outputs: Dict[str, Tensor],
        batch: Dict[str, Tensor],
    ) -> Dict[str, Tensor]:
        """
        Compute per-task losses and weighted total.

        Parameters
        ----------
        outputs : dict
            Model outputs from TaskHeads.forward() — contains raw logits.
        batch : dict
            Batch dict from DataLoader. Expected keys:
            - label_preservation [B] (long, -1 for unlabeled)
            - deviation_class    [B] (long, -1 for unlabeled)
            - task_masks         dict[str, Tensor[B]] — 1.0 if task applies

        Returns
        -------
        dict[str, Tensor]
            Per-task scalar losses + "total" weighted sum.
        """
        label_preservation = batch["label_preservation"].float()  # [B]
        deviation_labels = batch.get("deviation_class", torch.full_like(batch["label_preservation"], -1))
        task_masks = batch.get("task_masks", {})

        losses: Dict[str, Tensor] = {}
        total_loss = torch.tensor(0.0, device=label_preservation.device)

        # Binary tasks
        for task in BINARY_TASKS:
            if task not in outputs:
                continue

            logits = outputs[task].squeeze(-1)  # [B]
            labels = label_preservation  # [B], 0/1/-1

            # Task mask: 1.0 if task applies to this sample
            task_mask = task_masks.get(task, torch.ones_like(labels))  # [B]

            # Combine: mask out unlabeled (label=-1) AND task-irrelevant samples
            valid_mask = (labels != -1).float() * task_mask  # [B]

            if valid_mask.sum() < 1:
                losses[task] = torch.tensor(0.0, device=logits.device)
                continue

            # BCE on valid samples
            # Replace -1 labels with 0 before BCE (they'll be masked out anyway)
            clamped_labels = labels.clamp(min=0)
            per_sample_loss = self.binary_losses[task](logits, clamped_labels)  # [B]
            masked_loss = (per_sample_loss * valid_mask).sum() / valid_mask.sum()
            losses[task] = masked_loss
            total_loss = total_loss + self.task_weights.get(task, 1.0) * masked_loss

        # Deviation class (multiclass)
        if "deviation_class" in outputs:
            dev_logits = outputs["deviation_class"]  # [B, 6]
            dev_labels_int = deviation_labels.long()  # [B]
            dev_loss = self.deviation_loss(dev_logits, dev_labels_int)
            if not torch.isnan(dev_loss):
                losses["deviation_class"] = dev_loss
                total_loss = total_loss + self.task_weights.get("deviation_class", 0.5) * dev_loss
            else:
                losses["deviation_class"] = torch.tensor(0.0, device=total_loss.device)

        losses["total"] = total_loss
        return losses


class PretrainLoss(nn.Module):
    """
    SSL pretraining loss for ProteusPretrainModel.

    Binary cross-entropy for domain_disrupted, nmd_triggered, tm_lost.
    MSE for disorder_delta regression.

    Parameters
    ----------
    loss_weights : dict | None
        Maps task name -> float weight.
    """

    def __init__(self, loss_weights: Optional[Dict[str, float]] = None) -> None:
        super().__init__()

        self.loss_weights = loss_weights or {
            "domain_disrupted": 1.0,
            "nmd_triggered": 1.0,
            "tm_lost": 1.0,
            "disorder_delta": 1.0,
        }

        self.bce = nn.BCEWithLogitsLoss()
        self.mse = nn.MSELoss()

    def forward(
        self,
        outputs: Dict[str, Tensor],
        ssl_labels: Dict[str, Tensor],
    ) -> Dict[str, Tensor]:
        """
        Compute SSL pretraining losses.

        Parameters
        ----------
        outputs : dict
            ProteusPretrainModel outputs (logits).
        ssl_labels : dict
            Keys: domain_disrupted [B], nmd_triggered [B], tm_lost [B],
            disorder_delta [B].

        Returns
        -------
        dict[str, Tensor]
            Per-task losses + "total".
        """
        losses: Dict[str, Tensor] = {}
        total = torch.tensor(0.0)

        for task in ("domain_disrupted", "nmd_triggered", "tm_lost"):
            if task in outputs and task in ssl_labels:
                logits = outputs[task].squeeze(-1)
                labels = ssl_labels[task].float()
                loss = self.bce(logits, labels)
                losses[task] = loss
                total = total + self.loss_weights.get(task, 1.0) * loss

        if "disorder_delta" in outputs and "disorder_delta" in ssl_labels:
            pred = outputs["disorder_delta"].squeeze(-1)
            target = ssl_labels["disorder_delta"].float()
            loss = self.mse(pred, target)
            losses["disorder_delta"] = loss
            total = total + self.loss_weights.get("disorder_delta", 1.0) * loss

        losses["total"] = total
        return losses
