"""
ProteusPretrainModel: Self-supervised pretraining model.

Uses the same delta + fusion architecture as ProteusModel but with
SSL task heads suited for synthetic exon-skip pretraining.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import torch
import torch.nn as nn
import yaml
from torch import Tensor

from .delta.reference_delta import ReferenceDelta
from .delta.cross_modal import CrossModalFusion


def _load_yaml(path: Path) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


class ProteusPretrainModel(nn.Module):
    """
    Proteus self-supervised pretraining model.

    Shares the ReferenceDelta + CrossModalFusion architecture with ProteusModel,
    but replaces the supervised task heads with SSL pretraining heads suited
    for synthetic exon-skip objectives:

    - domain_disrupted: binary — did the skip disrupt a domain/active site?
    - nmd_triggered: binary — does the skip create a premature stop codon triggering NMD?
    - tm_lost: binary — was a TM helix lost?
    - disorder_delta: regression — signed change in IDR content

    After pretraining, the delta + fusion weights can be transferred to
    ProteusModel via get_transferable_state_dict().

    Parameters
    ----------
    cfg : dict | None
        Model configuration. If None, loads from configs/model.yaml.
    """

    def __init__(self, cfg: Optional[Dict[str, Any]] = None) -> None:
        super().__init__()

        if cfg is None:
            cfg_path = Path(__file__).parent.parent / "configs" / "model.yaml"
            if cfg_path.exists():
                cfg = _load_yaml(cfg_path)
            else:
                cfg = {}

        self.cfg = cfg

        delta_cfg = cfg.get("delta", {})
        fusion_cfg = cfg.get("fusion", {})

        rna_dim = delta_cfg.get("rna_dim", 512)
        protein_dim = delta_cfg.get("protein_dim", 1280)
        disorder_proj_dim = delta_cfg.get("disorder_dim", 64)
        structural_proj_dim = delta_cfg.get("structural_dim", 64)
        hidden_dim = delta_cfg.get("delta_hidden_dim", 256)
        out_dim = delta_cfg.get("delta_out_dim", 256)
        n_heads = delta_cfg.get("n_heads", 8)
        dropout = delta_cfg.get("dropout", 0.1)
        common_dim = fusion_cfg.get("common_dim", 256)

        # Input projections (same as ProteusModel)
        self.disorder_proj = nn.Linear(50, disorder_proj_dim)
        self.structural_proj = nn.Linear(96, structural_proj_dim)

        # Shared delta modules
        self.rna_delta = ReferenceDelta(
            input_dim=rna_dim,
            hidden_dim=hidden_dim,
            out_dim=out_dim,
            n_heads=n_heads,
            dropout=dropout,
        )
        self.protein_delta = ReferenceDelta(
            input_dim=protein_dim,
            hidden_dim=hidden_dim,
            out_dim=out_dim,
            n_heads=n_heads,
            dropout=dropout,
        )
        self.disorder_delta = ReferenceDelta(
            input_dim=disorder_proj_dim,
            hidden_dim=hidden_dim // 2,
            out_dim=out_dim,
            n_heads=max(1, n_heads // 2),
            dropout=dropout,
        )
        self.structural_delta = ReferenceDelta(
            input_dim=structural_proj_dim,
            hidden_dim=hidden_dim // 2,
            out_dim=out_dim,
            n_heads=max(1, n_heads // 2),
            dropout=dropout,
        )

        # Shared fusion
        self.fusion = CrossModalFusion(
            common_dim=common_dim,
            mlp_hidden=fusion_cfg.get("mlp_hidden", 512),
            mlp_layers=fusion_cfg.get("mlp_layers", 3),
            dropout=fusion_cfg.get("dropout", 0.15),
            use_cross_modal_attention=fusion_cfg.get("use_cross_modal_attention", True),
            modality_dropout_prob=fusion_cfg.get("modality_dropout_prob", 0.1),
        )

        # SSL pretraining heads
        self.domain_disrupted_head = nn.Linear(common_dim, 1)    # binary
        self.nmd_triggered_head = nn.Linear(common_dim, 1)       # binary
        self.tm_lost_head = nn.Linear(common_dim, 1)             # binary
        self.disorder_delta_head = nn.Linear(common_dim, 1)      # regression

    def forward(self, batch: Dict[str, Any]) -> Dict[str, Tensor]:
        """
        Forward pass for SSL pretraining.

        Parameters
        ----------
        batch : dict
            Same embedding fields as ProteusModel.forward().

        Returns
        -------
        dict[str, Tensor]
            SSL task outputs:
            - domain_disrupted [B, 1] raw logits
            - nmd_triggered    [B, 1] raw logits
            - tm_lost          [B, 1] raw logits
            - disorder_delta   [B, 1] regression predictions
        """
        ref_rna = batch["ref_rna_emb"]
        alt_rna = batch["alt_rna_emb"]
        ref_prot = batch["ref_protein_emb"]
        alt_prot = batch["alt_protein_emb"]
        ref_dis_raw = batch["ref_disorder_raw"]
        alt_dis_raw = batch["alt_disorder_raw"]
        ref_str_raw = batch["ref_structural_raw"]
        alt_str_raw = batch["alt_structural_raw"]

        has_alt_protein = batch.get(
            "has_alt_protein",
            torch.ones(ref_rna.size(0), dtype=torch.bool, device=ref_rna.device),
        )

        # Zero out absent alt proteins
        alt_prot = alt_prot.clone()
        absent_mask = ~has_alt_protein.bool()
        if absent_mask.any():
            alt_prot[absent_mask] = 0.0

        # Project raw features
        ref_dis = self.disorder_proj(ref_dis_raw)
        alt_dis = self.disorder_proj(alt_dis_raw)
        ref_str = self.structural_proj(ref_str_raw)
        alt_str = self.structural_proj(alt_str_raw)

        # Compute deltas
        rna_delta = self.rna_delta(ref_rna, alt_rna)
        protein_delta = self.protein_delta(ref_prot, alt_prot)
        dis_delta = self.disorder_delta(ref_dis, alt_dis)
        str_delta = self.structural_delta(ref_str, alt_str)

        # Fuse
        fused = self.fusion(rna_delta, protein_delta, dis_delta, str_delta)

        # SSL heads
        return {
            "domain_disrupted": self.domain_disrupted_head(fused),   # [B, 1]
            "nmd_triggered": self.nmd_triggered_head(fused),          # [B, 1]
            "tm_lost": self.tm_lost_head(fused),                      # [B, 1]
            "disorder_delta": self.disorder_delta_head(fused),        # [B, 1]
        }

    def get_transferable_state_dict(self) -> Dict[str, Tensor]:
        """
        Returns only the delta + fusion weights suitable for transfer
        to a ProteusModel (excludes SSL task heads).

        Returns
        -------
        dict[str, Tensor]
            State dict containing delta and fusion parameters.
        """
        exclude_prefixes = (
            "domain_disrupted_head",
            "nmd_triggered_head",
            "tm_lost_head",
            "disorder_delta_head",
        )
        return {
            k: v
            for k, v in self.state_dict().items()
            if not any(k.startswith(pfx) for pfx in exclude_prefixes)
        }

    def save_checkpoint(
        self,
        path: Path,
        epoch: int,
        metrics: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Save pretrain checkpoint with transferable state dict."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        torch.save(
            {
                "model_state_dict": self.state_dict(),
                "transferable_state_dict": self.get_transferable_state_dict(),
                "cfg": self.cfg,
                "epoch": epoch,
                "metrics": metrics or {},
                "model_class": "ProteusPretrainModel",
            },
            path,
        )
        print(f"[ProteusPretrainModel] Checkpoint saved: {path} (epoch {epoch})")

    def count_parameters(self) -> int:
        """Return total number of trainable parameters."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
