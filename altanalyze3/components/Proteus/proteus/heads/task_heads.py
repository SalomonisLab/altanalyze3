"""
TaskHeads: Multi-task prediction heads for Proteus.

Ten prediction heads operating on the fused 256-dim representation:
nine binary heads and one 8-class deviation classification head.

New in v1.1:
  - localization_preserved   : subcellular localization maintained
  - tm_insertion_competent   : novel TM isoform can properly insert/fold into membrane
  - signaling_competent      : transmembrane receptor can transmit signal
  - ppi_interface_preserved  : ≥1 PPI binding interface retained
  - deviation_class expanded : 8 classes (added localization_change, ppi_disrupted)
"""

from __future__ import annotations

from typing import Dict, Optional

import torch
import torch.nn as nn
from torch import Tensor

# Deviation class taxonomy (8 classes)
DEVIATION_CLASSES: Dict[int, str] = {
    0: "preserved",
    1: "truncation_partial",
    2: "truncation_nmd",
    3: "topology_loss",
    4: "domain_loss",
    5: "novel_fusion",
    6: "localization_change",
    7: "ppi_disrupted",
}

# All binary task names (used by model.predict() and loss computation)
BINARY_TASKS = [
    # Original tasks
    "global_preservation",
    "topology_preserved",
    "surface_retained",
    "kinase_competent",
    "tf_competent",
    "disorder_preserved",
    # New in v1.1
    "localization_preserved",
    "tm_insertion_competent",
    "signaling_competent",
    "ppi_interface_preserved",
]

# Task → mask column in dataset (None = always applicable)
TASK_MASK_COLUMNS: Dict[str, Optional[str]] = {
    "global_preservation":      None,           # always computed
    "topology_preserved":       "is_membrane_pair",
    "surface_retained":         "is_surface_pair",
    "kinase_competent":         "is_kinase_pair",
    "tf_competent":             "is_tf_pair",
    "disorder_preserved":       None,            # always computed
    "localization_preserved":   "is_localization_pair",
    "tm_insertion_competent":   "is_membrane_pair",  # membrane proteins only
    "signaling_competent":      "is_signaling_pair",
    "ppi_interface_preserved":  "is_ppi_pair",
}


class _BinaryHead(nn.Module):
    """Small MLP binary head: Linear → GELU → LayerNorm → Linear → scalar logit."""

    def __init__(self, fused_dim: int, hidden_dim: int = 64) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(fused_dim, hidden_dim),
            nn.GELU(),
            nn.LayerNorm(hidden_dim),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, x: Tensor) -> Tensor:
        return self.net(x)  # [B, 1]


class TaskHeads(nn.Module):
    """
    Multi-task prediction heads operating on the Proteus fused representation.

    All binary heads output raw logits [B, 1] — sigmoid is applied during
    loss computation or in .predict(), not in forward().
    The deviation_class head outputs [B, n_classes] logits.

    Parameters
    ----------
    fused_dim : int
        Dimensionality of the fused representation from CrossModalFusion.
    n_deviation_classes : int
        Number of deviation classes (default 8).
    head_hidden_dim : int
        Hidden dim for each binary head MLP (default 64).
    """

    def __init__(
        self,
        fused_dim: int = 256,
        n_deviation_classes: int = 8,
        head_hidden_dim: int = 64,
    ) -> None:
        super().__init__()

        self.fused_dim = fused_dim
        self.n_deviation_classes = n_deviation_classes

        # ---- Original binary heads ----
        self.global_preservation     = _BinaryHead(fused_dim, head_hidden_dim)
        self.topology_preserved      = _BinaryHead(fused_dim, head_hidden_dim)
        self.surface_retained        = _BinaryHead(fused_dim, head_hidden_dim)
        self.kinase_competent        = _BinaryHead(fused_dim, head_hidden_dim)
        self.tf_competent            = _BinaryHead(fused_dim, head_hidden_dim)
        self.disorder_preserved      = _BinaryHead(fused_dim, head_hidden_dim)

        # ---- New v1.1 binary heads ----
        self.localization_preserved  = _BinaryHead(fused_dim, head_hidden_dim)
        self.tm_insertion_competent  = _BinaryHead(fused_dim, head_hidden_dim)
        self.signaling_competent     = _BinaryHead(fused_dim, head_hidden_dim)
        self.ppi_interface_preserved = _BinaryHead(fused_dim, head_hidden_dim)

        # ---- Deviation classification head (8-class) ----
        self.deviation_class = nn.Sequential(
            nn.Linear(fused_dim, 256),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(256, 128),
            nn.GELU(),
            nn.Linear(128, n_deviation_classes),
        )

    def forward(
        self,
        fused_repr: Tensor,
        task_masks: Optional[Dict[str, Tensor]] = None,
    ) -> Dict[str, Tensor]:
        """
        Compute raw logits for all tasks.

        Parameters
        ----------
        fused_repr : Tensor [B, fused_dim]
            Fused representation from CrossModalFusion.
        task_masks : dict[str, Tensor] | None
            Per-task boolean/float masks [B]. Stored in output dict as
            "{task}_mask" keys for use in loss computation.

        Returns
        -------
        dict[str, Tensor]
            Task name → raw logits. Binary tasks: [B, 1].
            deviation_class: [B, n_classes].
            "{task}_mask" entries if task_masks provided.
        """
        outputs: Dict[str, Tensor] = {}

        # Binary heads
        for task in BINARY_TASKS:
            head = getattr(self, task)
            outputs[task] = head(fused_repr)

        # Multi-class deviation head
        outputs["deviation_class"] = self.deviation_class(fused_repr)

        # Store masks
        if task_masks is not None:
            for task, mask in task_masks.items():
                outputs[f"{task}_mask"] = mask

        return outputs

    def predict(self, fused_repr: Tensor) -> Dict[str, Tensor]:
        """
        Calibrated predictions (sigmoid/softmax applied) for inference.

        Parameters
        ----------
        fused_repr : Tensor [B, fused_dim]

        Returns
        -------
        dict[str, Tensor]
            Binary tasks: probabilities [B, 1] ∈ [0, 1].
            deviation_class: class probabilities [B, n_classes] summing to 1.
        """
        with torch.no_grad():
            raw = self.forward(fused_repr)

        preds: Dict[str, Tensor] = {}
        for task in BINARY_TASKS:
            preds[task] = torch.sigmoid(raw[task])

        preds["deviation_class"] = torch.softmax(raw["deviation_class"], dim=-1)
        return preds

    def extra_repr(self) -> str:
        return (
            f"fused_dim={self.fused_dim}, "
            f"n_deviation_classes={self.n_deviation_classes}, "
            f"n_binary_heads={len(BINARY_TASKS)}"
        )
