"""
PretrainObjective: Loss computation for SSL pretraining.

Wraps BCE (binary tasks) and MSE (disorder_delta regression) into
a single interface with configurable task weights.
"""

from __future__ import annotations

from typing import Dict, Optional

import torch
import torch.nn as nn
from torch import Tensor


class PretrainObjective(nn.Module):
    """
    SSL pretraining objective for ProteusPretrainModel.

    Computes:
    - Binary cross-entropy for: domain_disrupted, nmd_triggered, tm_lost
    - Mean squared error for: disorder_delta (regression)

    Returns per-task losses and a weighted total.

    Parameters
    ----------
    loss_weights : dict | None
        Maps task name -> float weight. Defaults to 1.0 for all tasks.
    """

    BINARY_TASKS = ("domain_disrupted", "nmd_triggered", "tm_lost")
    REGRESSION_TASKS = ("disorder_delta",)

    def __init__(self, loss_weights: Optional[Dict[str, float]] = None) -> None:
        super().__init__()

        self.loss_weights = loss_weights or {
            "domain_disrupted": 1.0,
            "nmd_triggered": 1.0,
            "tm_lost": 1.0,
            "disorder_delta": 1.0,
        }

        # Per-task BCE (reduction=mean by default, but we allow filtering NaN)
        self.bce = nn.BCEWithLogitsLoss(reduction="none")
        self.mse = nn.MSELoss(reduction="none")

    def compute(
        self,
        outputs: Dict[str, Tensor],
        labels: Dict[str, Tensor],
    ) -> Dict[str, Tensor]:
        """
        Compute SSL losses.

        Parameters
        ----------
        outputs : dict
            ProteusPretrainModel outputs (raw logits / regression predictions).
        labels : dict
            SSL labels: domain_disrupted [B], nmd_triggered [B], tm_lost [B],
            disorder_delta [B].

        Returns
        -------
        dict[str, Tensor]
            Per-task scalar losses and "total" weighted sum.
        """
        losses: Dict[str, Tensor] = {}
        total = torch.tensor(0.0, requires_grad=False)

        for task in self.BINARY_TASKS:
            if task not in outputs or task not in labels:
                continue
            logits = outputs[task].squeeze(-1).float()
            y = labels[task].float()

            # Filter any NaN labels
            valid = ~torch.isnan(y)
            if valid.sum() < 1:
                losses[task] = torch.tensor(0.0, device=logits.device)
                continue

            per_sample = self.bce(logits[valid], y[valid])
            loss = per_sample.mean()
            losses[task] = loss
            total = total.to(loss.device) + self.loss_weights.get(task, 1.0) * loss

        for task in self.REGRESSION_TASKS:
            if task not in outputs or task not in labels:
                continue
            pred = outputs[task].squeeze(-1).float()
            target = labels[task].float()

            valid = ~torch.isnan(target)
            if valid.sum() < 1:
                losses[task] = torch.tensor(0.0, device=pred.device)
                continue

            per_sample = self.mse(pred[valid], target[valid])
            loss = per_sample.mean()
            losses[task] = loss
            total = total.to(loss.device) + self.loss_weights.get(task, 1.0) * loss

        losses["total"] = total
        return losses

    def forward(
        self,
        outputs: Dict[str, Tensor],
        labels: Dict[str, Tensor],
    ) -> Dict[str, Tensor]:
        """Alias for compute() to support nn.Module interface."""
        return self.compute(outputs, labels)
