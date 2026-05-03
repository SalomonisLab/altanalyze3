"""
ReferenceDelta: Asymmetric reference-conditioned delta module.

This is the core novel architectural contribution of Proteus.
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor


class ReferenceDelta(nn.Module):
    """
    Asymmetric reference-conditioned delta module.

    Given (ref_emb, alt_emb) both of shape [B, input_dim], produces
    an asymmetric delta representation [B, out_dim] encoding what the
    alternative isoform lost or gained relative to the reference.

    Unlike a simple vector difference (alt - ref), this module captures:
    1. What the reference had that the alternative LOST (loss gate)
    2. What the alternative has that the reference LACKED (gain gate)
    3. How the alternative 'sees' the reference (cross-attention, asymmetric)

    Asymmetry is critical: (ref=kinase, alt=truncated) should produce a
    very different delta than (ref=truncated, alt=kinase).

    Parameters
    ----------
    input_dim : int
        Dimensionality of input embeddings (e.g., 512 for RNA, 1280 for protein).
    hidden_dim : int
        Internal hidden dimension for gating and projection.
    out_dim : int
        Output delta dimensionality.
    n_heads : int
        Number of attention heads for cross-attention.
    dropout : float
        Dropout probability.
    """

    def __init__(
        self,
        input_dim: int,
        hidden_dim: int = 256,
        out_dim: int = 256,
        n_heads: int = 8,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()

        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.out_dim = out_dim

        # Reference and alternative projection heads (separate — asymmetric by design)
        self.ref_proj = nn.Linear(input_dim, hidden_dim)
        self.alt_proj = nn.Linear(input_dim, hidden_dim)

        # Cross-attention: alt attends to ref (asymmetric direction)
        self.cross_attn = nn.MultiheadAttention(
            hidden_dim,
            n_heads,
            batch_first=True,
            dropout=dropout,
        )

        # Loss gate: what ref had that alt lost
        self.loss_gate_fc = nn.Linear(2 * hidden_dim, hidden_dim)
        nn.init.xavier_uniform_(self.loss_gate_fc.weight, gain=0.1)
        nn.init.zeros_(self.loss_gate_fc.bias)

        # Gain gate: what alt has that ref lacked
        self.gain_gate_fc = nn.Linear(2 * hidden_dim, hidden_dim)
        nn.init.xavier_uniform_(self.gain_gate_fc.weight, gain=0.1)
        nn.init.zeros_(self.gain_gate_fc.bias)

        # Final projection: 4H → out_dim
        self.delta_proj = nn.Linear(4 * hidden_dim, out_dim)
        self.layer_norm = nn.LayerNorm(out_dim)
        self.dropout = nn.Dropout(dropout)

        # Learned shortcut scalar — starts at 0 (shortcut initially off)
        self.shortcut_alpha = nn.Parameter(torch.zeros(1))

        # Shortcut projection for residual connection
        if input_dim != out_dim:
            self.shortcut_proj = nn.Linear(input_dim, out_dim)
        else:
            self.shortcut_proj = nn.Identity()

    def forward(self, ref_emb: Tensor, alt_emb: Tensor) -> Tensor:
        """
        Compute asymmetric delta representation.

        Parameters
        ----------
        ref_emb : Tensor
            Reference isoform embeddings, shape [B, input_dim].
        alt_emb : Tensor
            Alternative isoform embeddings, shape [B, input_dim].

        Returns
        -------
        Tensor
            Asymmetric delta representation, shape [B, out_dim].
        """
        # 1. Project to hidden space
        ref_h = self.ref_proj(ref_emb)  # [B, H]
        alt_h = self.alt_proj(alt_emb)  # [B, H]

        # 2. Cross-attention: alt (query) attends to ref (key/value)
        #    This captures what the alt "needs from" or "differs from" the ref.
        alt_h_seq = alt_h.unsqueeze(1)   # [B, 1, H]
        ref_h_seq = ref_h.unsqueeze(1)   # [B, 1, H]
        attended, _ = self.cross_attn(
            alt_h_seq,  # query = alt
            ref_h_seq,  # key = ref
            ref_h_seq,  # value = ref
        )  # [B, 1, H]
        attended = self.dropout(attended.squeeze(1))  # [B, H]

        # 3. Loss gate: what ref had that alt lost
        #    ref_h - alt_h is positive where ref was "stronger"
        loss_input = torch.cat([ref_h, alt_h], dim=-1)  # [B, 2H]
        loss_gate = torch.sigmoid(self.loss_gate_fc(loss_input))  # [B, H]
        diff_rl = ref_h - alt_h  # [B, H] — positive where ref > alt
        # Normalize direction but preserve magnitude
        diff_rl_norm = (
            F.normalize(diff_rl, dim=-1)
            * diff_rl.norm(dim=-1, keepdim=True).clamp(min=1e-8)
        )
        loss_signal = loss_gate * diff_rl_norm  # [B, H]

        # 4. Gain gate: what alt has that ref lacked
        #    alt_h - ref_h is positive where alt "added something"
        gain_input = torch.cat([alt_h, ref_h], dim=-1)  # [B, 2H]
        gain_gate = torch.sigmoid(self.gain_gate_fc(gain_input))  # [B, H]
        diff_ar = alt_h - ref_h  # [B, H] — positive where alt > ref
        diff_ar_norm = (
            F.normalize(diff_ar, dim=-1)
            * diff_ar.norm(dim=-1, keepdim=True).clamp(min=1e-8)
        )
        gain_signal = gain_gate * diff_ar_norm  # [B, H]

        # 5. Concatenate all four signals
        # [loss_signal, gain_signal, cross-attended, raw_diff] → [B, 4H]
        delta_input = torch.cat(
            [loss_signal, gain_signal, attended, ref_h - alt_h],
            dim=-1,
        )  # [B, 4H]

        # 6. Project and normalize
        delta = self.layer_norm(self.delta_proj(delta_input))  # [B, out_dim]

        # 7. Learned shortcut: alpha starts at 0, so shortcut is initially off.
        #    As training proceeds, the model can enable it if helpful.
        shortcut = self.shortcut_proj(ref_emb)  # [B, out_dim]
        delta = delta + torch.tanh(self.shortcut_alpha) * shortcut

        return delta

    def extra_repr(self) -> str:
        return (
            f"input_dim={self.input_dim}, hidden_dim={self.hidden_dim}, "
            f"out_dim={self.out_dim}"
        )


if __name__ == "__main__":
    # Smoke test for the ReferenceDelta module
    import sys

    print("Testing ReferenceDelta...")

    B = 8
    input_dim = 512
    hidden_dim = 256
    out_dim = 256

    module = ReferenceDelta(
        input_dim=input_dim,
        hidden_dim=hidden_dim,
        out_dim=out_dim,
        n_heads=8,
        dropout=0.0,  # disable dropout for determinism
    )
    module.eval()

    ref = torch.randn(B, input_dim, requires_grad=True)
    alt = torch.randn(B, input_dim, requires_grad=True)

    # Test 1: Output shape
    out = module(ref, alt)
    assert out.shape == (B, out_dim), f"Expected [{B}, {out_dim}], got {out.shape}"
    print(f"  [PASS] Output shape: {out.shape}")

    # Test 2: Asymmetry — swapping ref and alt should produce different output
    out_swapped = module(alt, ref)
    diff = (out - out_swapped).abs().sum().item()
    assert diff > 1e-3, (
        f"Asymmetry test FAILED: ReferenceDelta(ref, alt) == ReferenceDelta(alt, ref) "
        f"(diff={diff:.6f}). Module is not truly asymmetric."
    )
    print(f"  [PASS] Asymmetry: diff between (ref,alt) and (alt,ref) = {diff:.4f}")

    # Test 3: Gradient flow to both inputs
    out.sum().backward()
    assert ref.grad is not None and ref.grad.abs().sum() > 0, "No gradient to ref_emb"
    assert alt.grad is not None and alt.grad.abs().sum() > 0, "No gradient to alt_emb"
    print(f"  [PASS] Gradient to ref: {ref.grad.abs().mean():.4e}")
    print(f"  [PASS] Gradient to alt: {alt.grad.abs().mean():.4e}")

    # Test 4: NaN-free outputs
    assert not torch.isnan(out).any(), "NaN detected in output!"
    print("  [PASS] No NaNs in output")

    # Test 5: Different input dims
    for dim in [512, 1280, 64]:
        m = ReferenceDelta(input_dim=dim, hidden_dim=128, out_dim=128, n_heads=4)
        r = torch.randn(4, dim)
        a = torch.randn(4, dim)
        o = m(r, a)
        assert o.shape == (4, 128), f"Dim {dim}: expected [4, 128], got {o.shape}"
        print(f"  [PASS] input_dim={dim}: output shape {o.shape}")

    print("\nAll ReferenceDelta tests passed!")
    sys.exit(0)
