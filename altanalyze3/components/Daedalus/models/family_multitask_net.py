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


class FamilyMultitaskNet(nn.Module):
    def __init__(
        self,
        n_numeric: int,
        n_tasks: int,
        n_families: int,
        n_reference_sources: int,
        n_change_classes: int,
        shared_hidden: int = 192,
        expert_hidden: int = 96,
        latent_dim: int = 64,
        dropout: float = 0.16,
    ) -> None:
        super().__init__()
        self.task_embed = nn.Embedding(n_tasks, 12)
        self.family_embed = nn.Embedding(n_families, 8)
        self.refsrc_embed = nn.Embedding(n_reference_sources, 8)
        self.change_embed = nn.Embedding(n_change_classes, 6)

        self.shared = nn.Sequential(
            nn.Linear(n_numeric, shared_hidden),
            nn.LayerNorm(shared_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(shared_hidden, shared_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
        )
        self.membrane_expert = nn.Sequential(
            nn.Linear(shared_hidden, expert_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(expert_hidden, expert_hidden),
            nn.GELU(),
        )
        self.tf_expert = nn.Sequential(
            nn.Linear(shared_hidden, expert_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(expert_hidden, expert_hidden),
            nn.GELU(),
        )
        self.other_expert = nn.Sequential(
            nn.Linear(shared_hidden, expert_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(expert_hidden, expert_hidden),
            nn.GELU(),
        )
        self.lat_proj = nn.Sequential(
            nn.Linear(shared_hidden + expert_hidden + 8 + 8 + 6, latent_dim),
            nn.LayerNorm(latent_dim),
            nn.GELU(),
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, shared_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(shared_hidden, n_numeric),
        )
        self.task_heads = nn.ModuleList(
            [
                nn.Sequential(
                    nn.Linear(latent_dim + expert_hidden + 12 + 8 + 6, 96),
                    nn.GELU(),
                    nn.Dropout(dropout),
                    nn.Linear(96, 1),
                )
                for _ in range(n_tasks)
            ]
        )

    def _route(self, shared: torch.Tensor, family_id: torch.Tensor) -> torch.Tensor:
        batch = shared.shape[0]
        out = shared.new_zeros((batch, self.membrane_expert[0].out_features))
        membrane_mask = family_id == 0
        tf_mask = family_id == 1
        other_mask = ~(membrane_mask | tf_mask)
        if torch.any(membrane_mask):
            out[membrane_mask] = self.membrane_expert(shared[membrane_mask])
        if torch.any(tf_mask):
            out[tf_mask] = self.tf_expert(shared[tf_mask])
        if torch.any(other_mask):
            out[other_mask] = self.other_expert(shared[other_mask])
        return out

    def forward(
        self,
        numeric: torch.Tensor,
        task_id: torch.Tensor,
        family_id: torch.Tensor,
        refsrc_id: torch.Tensor,
        change_id: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        shared = self.shared(numeric)
        expert = self._route(shared, family_id)
        fam = self.family_embed(family_id)
        refsrc = self.refsrc_embed(refsrc_id)
        change = self.change_embed(change_id)
        latent = self.lat_proj(torch.cat([shared, expert, fam, refsrc, change], dim=1))
        recon = self.decoder(latent)
        task_context = torch.cat([latent, expert, self.task_embed(task_id), fam, change], dim=1)
        task_logits = torch.cat([head(task_context) for head in self.task_heads], dim=1)
        logits = task_logits.gather(1, task_id.unsqueeze(1)).squeeze(1)
        return logits, recon, latent
