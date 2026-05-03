"""
ProteusModel: Reference-conditioned isoform function prediction model.

The main trainable model. Uses pre-cached embeddings from frozen foundation
models (Orthrus, ESM2, disorder encoders) as input.
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
from .heads.task_heads import TaskHeads


def _load_yaml(path: Path) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


class ProteusModel(nn.Module):
    """
    Proteus: Reference-conditioned isoform function prediction model.

    Predicts whether an alternative RNA isoform preserves the core function
    of its reference isoform across seven biological dimensions.

    Architecture:
    - Four ReferenceDelta modules operating on frozen foundation model embeddings:
      * RNA delta: from Orthrus embeddings (512-dim)
      * Protein delta: from ESM2 embeddings (1280-dim)
      * Disorder delta: from STARLING/metapredict features (64-dim projected)
      * Structural delta: from UniProt/topology annotations (64-dim projected)
    - CrossModalFusion: fuses all four deltas into a 256-dim representation
    - TaskHeads: seven multi-task prediction heads

    The encoders (Orthrus, ESM2) are NOT part of this model — they are used
    in extraction scripts and their outputs are cached. Only the delta
    architecture is trainable.

    Parameters
    ----------
    cfg : dict | None
        Model configuration dict. If None, loads from configs/model.yaml.
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
        heads_cfg = cfg.get("heads", {})

        rna_dim = delta_cfg.get("rna_dim", 512)
        protein_dim = delta_cfg.get("protein_dim", 1280)
        disorder_proj_dim = delta_cfg.get("disorder_dim", 64)
        structural_proj_dim = delta_cfg.get("structural_dim", 64)
        hidden_dim = delta_cfg.get("delta_hidden_dim", 256)
        out_dim = delta_cfg.get("delta_out_dim", 256)
        n_heads = delta_cfg.get("n_heads", 8)
        dropout = delta_cfg.get("dropout", 0.1)

        common_dim = fusion_cfg.get("common_dim", 256)

        # Input projections for disorder and structural (raw features → projected)
        self.disorder_proj = nn.Linear(50, disorder_proj_dim)
        self.structural_proj = nn.Linear(96, structural_proj_dim)

        # Four ReferenceDelta modules
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

        # CrossModalFusion
        self.fusion = CrossModalFusion(
            common_dim=common_dim,
            mlp_hidden=fusion_cfg.get("mlp_hidden", 512),
            mlp_layers=fusion_cfg.get("mlp_layers", 3),
            dropout=fusion_cfg.get("dropout", 0.15),
            use_cross_modal_attention=fusion_cfg.get("use_cross_modal_attention", True),
            modality_dropout_prob=fusion_cfg.get("modality_dropout_prob", 0.1),
        )

        # Task heads
        n_dev_classes = heads_cfg.get("deviation_class", {}).get("n_classes", 8)
        head_hidden = heads_cfg.get("head_hidden_dim", 64)
        self.heads = TaskHeads(
            fused_dim=common_dim,
            n_deviation_classes=n_dev_classes,
            head_hidden_dim=head_hidden,
        )

    def forward(self, batch: Dict[str, Any]) -> Dict[str, Tensor]:
        """
        Forward pass.

        Parameters
        ----------
        batch : dict
            Expected keys:
            - ref_rna_emb      [B, 512]  — Orthrus RNA embedding of reference
            - alt_rna_emb      [B, 512]  — Orthrus RNA embedding of alternative
            - ref_protein_emb  [B, 1280] — ESM2 embedding of reference protein
            - alt_protein_emb  [B, 1280] — ESM2 embedding of alternative protein
            - ref_disorder_raw [B, 50]   — disorder features of reference protein
            - alt_disorder_raw [B, 50]   — disorder features of alternative protein
            - ref_structural_raw [B, 96] — structural features of reference protein
            - alt_structural_raw [B, 96] — structural features of alternative protein
            - has_alt_protein  [B]       — bool, False if alt protein is absent (NMD)
            - task_masks       dict      — per-task boolean masks [B]

        Returns
        -------
        dict[str, Tensor]
            Outputs from TaskHeads: logits per task + mask entries.
        """
        # Extract and move to device
        ref_rna = batch["ref_rna_emb"]
        alt_rna = batch["alt_rna_emb"]
        ref_prot = batch["ref_protein_emb"]
        alt_prot = batch["alt_protein_emb"]
        ref_dis_raw = batch["ref_disorder_raw"]
        alt_dis_raw = batch["alt_disorder_raw"]
        ref_str_raw = batch["ref_structural_raw"]
        alt_str_raw = batch["alt_structural_raw"]
        has_alt_protein = batch.get("has_alt_protein", torch.ones(ref_rna.size(0), dtype=torch.bool))
        task_masks = batch.get("task_masks", None)

        # Zero out protein embeddings for NMD/absent alternative proteins
        # This forces the protein delta to represent "complete loss"
        if has_alt_protein is not None:
            alt_prot = alt_prot.clone()
            absent_mask = ~has_alt_protein.bool()
            if absent_mask.any():
                alt_prot[absent_mask] = 0.0

        # Project disorder and structural raw features to common dim
        ref_dis = self.disorder_proj(ref_dis_raw)    # [B, 64]
        alt_dis = self.disorder_proj(alt_dis_raw)    # [B, 64]
        ref_str = self.structural_proj(ref_str_raw)  # [B, 64]
        alt_str = self.structural_proj(alt_str_raw)  # [B, 64]

        # Compute four asymmetric deltas
        rna_delta = self.rna_delta(ref_rna, alt_rna)         # [B, 256]
        protein_delta = self.protein_delta(ref_prot, alt_prot)  # [B, 256]
        dis_delta = self.disorder_delta(ref_dis, alt_dis)    # [B, 256]
        str_delta = self.structural_delta(ref_str, alt_str)  # [B, 256]

        # Fuse across modalities
        fused = self.fusion(rna_delta, protein_delta, dis_delta, str_delta)  # [B, 256]

        # Multi-task prediction
        outputs = self.heads(fused, task_masks=task_masks)

        return outputs

    def predict(self, batch: Dict[str, Any]) -> Dict[str, Tensor]:
        """
        Run forward pass and return calibrated probabilities (inference mode).

        Parameters
        ----------
        batch : dict
            Same format as forward().

        Returns
        -------
        dict[str, Tensor]
            Calibrated probabilities per task.
        """
        self.eval()
        with torch.no_grad():
            raw_outputs = self.forward(batch)
            # Extract logits and apply sigmoid/softmax
            import torch.nn.functional as F
            from .heads.task_heads import BINARY_TASKS
            preds: Dict[str, Tensor] = {}
            for task in BINARY_TASKS:
                if task in raw_outputs:
                    preds[task] = torch.sigmoid(raw_outputs[task])
            if "deviation_class" in raw_outputs:
                preds["deviation_class"] = F.softmax(raw_outputs["deviation_class"], dim=-1)
        return preds

    def count_parameters(self) -> int:
        """Return total number of trainable parameters."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    def save_checkpoint(
        self,
        path: Path,
        epoch: int,
        metrics: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Save model checkpoint.

        Parameters
        ----------
        path : Path
            Save path (.pt file).
        epoch : int
            Current epoch number.
        metrics : dict | None
            Validation metrics to store alongside weights.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        checkpoint = {
            "model_state_dict": self.state_dict(),
            "cfg": self.cfg,
            "epoch": epoch,
            "metrics": metrics or {},
            "model_class": "ProteusModel",
        }
        torch.save(checkpoint, path)
        print(f"[ProteusModel] Checkpoint saved: {path} (epoch {epoch})")

    @classmethod
    def load_checkpoint(cls, path: Path) -> "ProteusModel":
        """
        Load a ProteusModel from a checkpoint file.

        Parameters
        ----------
        path : Path
            Checkpoint .pt file path.

        Returns
        -------
        ProteusModel
            Model with loaded weights.
        """
        path = Path(path)
        checkpoint = torch.load(path, map_location="cpu", weights_only=False)
        cfg = checkpoint.get("cfg", {})
        model = cls(cfg=cfg)
        model.load_state_dict(checkpoint["model_state_dict"])
        epoch = checkpoint.get("epoch", 0)
        metrics = checkpoint.get("metrics", {})
        print(f"[ProteusModel] Loaded checkpoint: {path} (epoch {epoch}, metrics: {metrics})")
        return model

    @classmethod
    def from_pretrain(cls, pretrain_checkpoint_path: Path, cfg: Optional[Dict] = None) -> "ProteusModel":
        """
        Create ProteusModel and initialize delta + fusion weights from a
        ProteusPretrainModel checkpoint.

        Parameters
        ----------
        pretrain_checkpoint_path : Path
            Path to pretrain checkpoint (saved by ProteusPretrainModel).
        cfg : dict | None
            Model config override.

        Returns
        -------
        ProteusModel
            Model with pretrained delta/fusion weights.
        """
        from .pretrain_model import ProteusPretrainModel

        path = Path(pretrain_checkpoint_path)
        pretrain_checkpoint = torch.load(path, map_location="cpu", weights_only=False)

        # Create ProteusModel
        model = cls(cfg=cfg)

        # Get transferable weights from pretrain checkpoint
        if "transferable_state_dict" in pretrain_checkpoint:
            transferable = pretrain_checkpoint["transferable_state_dict"]
        else:
            # Try to extract from full state dict
            pretrain_state = pretrain_checkpoint.get("model_state_dict", pretrain_checkpoint)
            transferable = {
                k: v for k, v in pretrain_state.items()
                if not k.startswith("domain_disrupted_head")
                and not k.startswith("nmd_triggered_head")
                and not k.startswith("tm_lost_head")
                and not k.startswith("disorder_delta_head")
            }

        # Load with strict=False to ignore missing task head weights
        missing, unexpected = model.load_state_dict(transferable, strict=False)
        print(f"[ProteusModel.from_pretrain] Loaded from {path}")
        if missing:
            print(f"  Missing keys (will be randomly initialized): {len(missing)}")
        if unexpected:
            print(f"  Unexpected keys (ignored): {len(unexpected)}")

        return model
