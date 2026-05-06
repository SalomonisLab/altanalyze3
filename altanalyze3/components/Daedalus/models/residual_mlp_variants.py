"""
Variants of ResidualMLPAutoencoder for TM-negative detection.

All variants share the same IO contract as the baseline so they plug into
_train_deep_binary in run_tm_isoform_model_comparison.py. Configurable
knobs:
    * hidden_dim      width of residual blocks
    * latent_dim      bottleneck width (autoencoder)
    * dropout         dropout inside blocks / head
    * n_blocks        number of residual blocks (1-4)
    * decoder_hidden  decoder hidden width
    * head_hidden     classification head hidden width
"""
from __future__ import annotations

import torch
from torch import nn

from .deepimmuno_style import _CatContext


class _ResBlock(nn.Module):
    def __init__(self, hidden_dim: int, dropout: float) -> None:
        super().__init__()
        self.body = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.post = nn.Sequential(nn.LayerNorm(hidden_dim), nn.GELU())

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.post(x + self.body(x))


class ResidualMLPAutoencoderVariant(nn.Module):
    """Generalised residual-MLP autoencoder with configurable depth/width."""

    def __init__(
        self,
        n_numeric: int,
        n_families: int,
        n_reference_sources: int,
        n_change_classes: int,
        hidden_dim: int = 192,
        latent_dim: int = 48,
        dropout: float = 0.15,
        n_blocks: int = 2,
        decoder_hidden: int | None = None,
        head_hidden: int = 96,
    ) -> None:
        super().__init__()
        self.context = _CatContext(n_families, n_reference_sources, n_change_classes)
        ctx_dim = 8 + 6 + 6  # matches _CatContext concat width
        self.input_proj = nn.Sequential(
            nn.Linear(n_numeric, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
        )
        self.blocks = nn.ModuleList(
            [_ResBlock(hidden_dim, dropout) for _ in range(max(1, n_blocks))]
        )
        self.lat_proj = nn.Sequential(
            nn.Linear(hidden_dim + ctx_dim, latent_dim),
            nn.LayerNorm(latent_dim),
            nn.GELU(),
        )
        dec_hidden = decoder_hidden or hidden_dim
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, dec_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(dec_hidden, n_numeric),
        )
        self.head = nn.Sequential(
            nn.Linear(latent_dim + ctx_dim, head_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(head_hidden, 1),
        )

    def forward(
        self,
        numeric: torch.Tensor,
        family_id: torch.Tensor,
        refsrc_id: torch.Tensor,
        change_id: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        x = self.input_proj(numeric)
        for block in self.blocks:
            x = block(x)
        ctx = self.context(family_id, refsrc_id, change_id)
        latent = self.lat_proj(torch.cat([x, ctx], dim=1))
        recon = self.decoder(latent)
        logits = self.head(torch.cat([latent, ctx], dim=1)).squeeze(1)
        return logits, recon, latent
