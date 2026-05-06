from __future__ import annotations

import torch
from torch import nn


def choose_device(requested: str = "auto") -> torch.device:
    requested = (requested or "auto").lower()
    if requested != "auto":
        return torch.device(requested)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


class _CatContext(nn.Module):
    def __init__(self, n_families: int, n_reference_sources: int, n_change_classes: int) -> None:
        super().__init__()
        self.family = nn.Embedding(n_families, 8)
        self.refsrc = nn.Embedding(n_reference_sources, 6)
        self.change = nn.Embedding(n_change_classes, 6)

    def forward(self, family_id: torch.Tensor, refsrc_id: torch.Tensor, change_id: torch.Tensor) -> torch.Tensor:
        return torch.cat(
            [self.family(family_id), self.refsrc(refsrc_id), self.change(change_id)],
            dim=1,
        )


class ResidualMLPAutoencoder(nn.Module):
    def __init__(
        self,
        n_numeric: int,
        n_families: int,
        n_reference_sources: int,
        n_change_classes: int,
        hidden_dim: int = 192,
        latent_dim: int = 48,
        dropout: float = 0.15,
    ) -> None:
        super().__init__()
        self.context = _CatContext(n_families, n_reference_sources, n_change_classes)
        self.input_proj = nn.Sequential(
            nn.Linear(n_numeric, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
        )
        self.block1 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.block2 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.post = nn.Sequential(nn.LayerNorm(hidden_dim), nn.GELU())
        self.lat_proj = nn.Sequential(
            nn.Linear(hidden_dim + 20, latent_dim),
            nn.LayerNorm(latent_dim),
            nn.GELU(),
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_numeric),
        )
        self.head = nn.Sequential(
            nn.Linear(latent_dim + 20, 96),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(96, 1),
        )

    def forward(
        self,
        numeric: torch.Tensor,
        family_id: torch.Tensor,
        refsrc_id: torch.Tensor,
        change_id: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        x = self.input_proj(numeric)
        x = self.post(x + self.block1(x))
        x = self.post(x + self.block2(x))
        ctx = self.context(family_id, refsrc_id, change_id)
        latent = self.lat_proj(torch.cat([x, ctx], dim=1))
        recon = self.decoder(latent)
        logits = self.head(torch.cat([latent, ctx], dim=1)).squeeze(1)
        return logits, recon, latent


class DeepImmunoStyleCNN(nn.Module):
    def __init__(
        self,
        n_numeric: int,
        n_families: int,
        n_reference_sources: int,
        n_change_classes: int,
        conv_channels: int = 48,
        hidden_dim: int = 96,
        dropout: float = 0.15,
    ) -> None:
        super().__init__()
        self.context = _CatContext(n_families, n_reference_sources, n_change_classes)
        self.conv = nn.Sequential(
            nn.Conv1d(1, conv_channels, kernel_size=5, padding=2),
            nn.GELU(),
            nn.BatchNorm1d(conv_channels),
            nn.Conv1d(conv_channels, conv_channels, kernel_size=3, padding=1),
            nn.GELU(),
            nn.BatchNorm1d(conv_channels),
            nn.MaxPool1d(kernel_size=2),
            nn.Conv1d(conv_channels, conv_channels * 2, kernel_size=3, padding=1),
            nn.GELU(),
            nn.BatchNorm1d(conv_channels * 2),
            nn.AdaptiveAvgPool1d(8),
        )
        self.proj = nn.Sequential(
            nn.Linear((conv_channels * 2) * 8 + 20, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
        )
        self.head = nn.Linear(hidden_dim, 1)

    def forward(
        self,
        numeric: torch.Tensor,
        family_id: torch.Tensor,
        refsrc_id: torch.Tensor,
        change_id: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        x = numeric.unsqueeze(1)
        feat = self.conv(x).flatten(1)
        ctx = self.context(family_id, refsrc_id, change_id)
        hidden = self.proj(torch.cat([feat, ctx], dim=1))
        logits = self.head(hidden).squeeze(1)
        return logits, hidden


class DeepImmunoStyleCNNAutoencoder(nn.Module):
    def __init__(
        self,
        n_numeric: int,
        n_families: int,
        n_reference_sources: int,
        n_change_classes: int,
        conv_channels: int = 48,
        hidden_dim: int = 96,
        latent_dim: int = 48,
        dropout: float = 0.15,
    ) -> None:
        super().__init__()
        self.context = _CatContext(n_families, n_reference_sources, n_change_classes)
        self.conv = nn.Sequential(
            nn.Conv1d(1, conv_channels, kernel_size=5, padding=2),
            nn.GELU(),
            nn.BatchNorm1d(conv_channels),
            nn.Conv1d(conv_channels, conv_channels, kernel_size=3, padding=1),
            nn.GELU(),
            nn.BatchNorm1d(conv_channels),
            nn.MaxPool1d(kernel_size=2),
            nn.Conv1d(conv_channels, conv_channels * 2, kernel_size=3, padding=1),
            nn.GELU(),
            nn.BatchNorm1d(conv_channels * 2),
            nn.AdaptiveAvgPool1d(8),
        )
        self.lat_proj = nn.Sequential(
            nn.Linear((conv_channels * 2) * 8 + 20, latent_dim),
            nn.LayerNorm(latent_dim),
            nn.GELU(),
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_numeric),
        )
        self.head = nn.Sequential(
            nn.Linear(latent_dim + 20, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

    def forward(
        self,
        numeric: torch.Tensor,
        family_id: torch.Tensor,
        refsrc_id: torch.Tensor,
        change_id: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        x = numeric.unsqueeze(1)
        feat = self.conv(x).flatten(1)
        ctx = self.context(family_id, refsrc_id, change_id)
        latent = self.lat_proj(torch.cat([feat, ctx], dim=1))
        recon = self.decoder(latent)
        logits = self.head(torch.cat([latent, ctx], dim=1)).squeeze(1)
        return logits, recon, latent
