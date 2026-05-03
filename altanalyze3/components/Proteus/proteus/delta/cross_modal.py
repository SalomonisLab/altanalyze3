"""
CrossModalFusion: Fuses RNA, protein, disorder, and structural delta
representations into a unified embedding via cross-modal attention and MLP fusion.
"""

from __future__ import annotations

import torch
import torch.nn as nn
from torch import Tensor


class CrossModalFusion(nn.Module):
    """
    Fuses four modality-specific delta representations (RNA, protein,
    disorder, structural) into a single unified representation.

    Cross-modal attention between RNA and protein deltas allows the model
    to reason about joint consequences: "this exon skip removed the exact
    region encoding the kinase activation loop".

    Architecture:
    1. Optional cross-modal attention: RNA <-> Protein
    2. Concatenate all four modality vectors: [B, 4 * common_dim]
    3. MLP with LayerNorm + GELU + Dropout per layer
    4. Output: [B, common_dim]

    During training, each modality is independently dropped out with
    probability modality_dropout_prob to improve robustness to missing data.

    Parameters
    ----------
    common_dim : int
        Dimensionality of each input modality delta (all must match).
    mlp_hidden : int
        Hidden dimension of the fusion MLP.
    mlp_layers : int
        Number of MLP layers (excluding final projection).
    dropout : float
        Dropout probability in MLP layers.
    use_cross_modal_attention : bool
        If True, applies cross-attention between RNA and protein deltas.
    modality_dropout_prob : float
        Probability of zeroing out each modality independently during training.
    """

    def __init__(
        self,
        common_dim: int = 256,
        mlp_hidden: int = 512,
        mlp_layers: int = 3,
        dropout: float = 0.15,
        use_cross_modal_attention: bool = True,
        modality_dropout_prob: float = 0.1,
    ) -> None:
        super().__init__()

        self.common_dim = common_dim
        self.use_cross_modal = use_cross_modal_attention
        self.modality_dropout_prob = modality_dropout_prob

        if use_cross_modal_attention:
            # RNA attends to protein
            self.rna_cross_attn = nn.MultiheadAttention(
                common_dim, num_heads=4, batch_first=True
            )
            # Protein attends to RNA
            self.protein_cross_attn = nn.MultiheadAttention(
                common_dim, num_heads=4, batch_first=True
            )
            self.norm_rna = nn.LayerNorm(common_dim)
            self.norm_protein = nn.LayerNorm(common_dim)

        # MLP fusion: 4*common_dim → mlp_hidden → ... → common_dim
        mlp_layers_list: list[nn.Module] = []

        in_dim = 4 * common_dim
        for i in range(mlp_layers - 1):
            mlp_layers_list.extend([
                nn.Linear(in_dim, mlp_hidden),
                nn.LayerNorm(mlp_hidden),
                nn.GELU(),
                nn.Dropout(dropout),
            ])
            in_dim = mlp_hidden

        # Final projection to common_dim
        mlp_layers_list.extend([
            nn.Linear(in_dim, common_dim),
            nn.LayerNorm(common_dim),
        ])

        self.mlp = nn.Sequential(*mlp_layers_list)

    def _modality_dropout(self, x: Tensor) -> Tensor:
        """
        Zero out tensor with probability modality_dropout_prob during training.

        Parameters
        ----------
        x : Tensor [B, D]

        Returns
        -------
        Tensor [B, D]
        """
        if not self.training or self.modality_dropout_prob <= 0.0:
            return x
        # Bernoulli per sample in batch
        mask = torch.bernoulli(
            torch.full((x.size(0), 1), 1.0 - self.modality_dropout_prob, device=x.device)
        )  # [B, 1]
        return x * mask

    def forward(
        self,
        rna_delta: Tensor,
        protein_delta: Tensor,
        disorder_delta: Tensor,
        structural_delta: Tensor,
    ) -> Tensor:
        """
        Fuse four modality delta vectors.

        Parameters
        ----------
        rna_delta : Tensor [B, common_dim]
            RNA ReferenceDelta output.
        protein_delta : Tensor [B, common_dim]
            Protein ReferenceDelta output.
        disorder_delta : Tensor [B, common_dim]
            Disorder ReferenceDelta output.
        structural_delta : Tensor [B, common_dim]
            Structural ReferenceDelta output.

        Returns
        -------
        Tensor [B, common_dim]
            Fused representation.
        """
        # 1. Modality dropout (training only)
        rna_delta = self._modality_dropout(rna_delta)
        protein_delta = self._modality_dropout(protein_delta)
        disorder_delta = self._modality_dropout(disorder_delta)
        structural_delta = self._modality_dropout(structural_delta)

        # 2. Cross-modal attention between RNA and protein
        if self.use_cross_modal:
            # RNA attends to protein
            rna_seq = rna_delta.unsqueeze(1)           # [B, 1, D]
            protein_seq = protein_delta.unsqueeze(1)   # [B, 1, D]

            rna_ctx, _ = self.rna_cross_attn(
                rna_seq,       # query = RNA
                protein_seq,   # key = protein
                protein_seq,   # value = protein
            )  # [B, 1, D]
            rna_attended = self.norm_rna(
                rna_delta + rna_ctx.squeeze(1)
            )  # [B, D] — residual connection

            protein_ctx, _ = self.protein_cross_attn(
                protein_seq,   # query = protein
                rna_seq,       # key = RNA
                rna_seq,       # value = RNA
            )  # [B, 1, D]
            protein_attended = self.norm_protein(
                protein_delta + protein_ctx.squeeze(1)
            )  # [B, D]
        else:
            rna_attended = rna_delta
            protein_attended = protein_delta

        # 3. Concatenate all four modalities
        fused = torch.cat(
            [rna_attended, protein_attended, disorder_delta, structural_delta],
            dim=-1,
        )  # [B, 4 * common_dim]

        # 4. MLP fusion
        out = self.mlp(fused)  # [B, common_dim]

        return out

    def extra_repr(self) -> str:
        return (
            f"common_dim={self.common_dim}, "
            f"use_cross_modal={self.use_cross_modal}, "
            f"modality_dropout={self.modality_dropout_prob}"
        )
