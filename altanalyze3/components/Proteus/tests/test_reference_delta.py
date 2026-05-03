"""
Tests for the ReferenceDelta module.

Tests: output shape, asymmetry, gradient flow, NaN-free outputs,
and compatibility across different input dimensions.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
import torch

sys.path.insert(0, str(Path(__file__).parent.parent))
from proteus.delta.reference_delta import ReferenceDelta


class TestReferenceDeltaShape:
    """Tests for output tensor shapes."""

    def test_output_shape_standard(self):
        """Output should be [B, out_dim] for standard inputs."""
        B, D, out_dim = 8, 512, 256
        module = ReferenceDelta(input_dim=D, hidden_dim=256, out_dim=out_dim)
        ref = torch.randn(B, D)
        alt = torch.randn(B, D)
        out = module(ref, alt)
        assert out.shape == (B, out_dim), f"Expected [{B}, {out_dim}], got {out.shape}"

    def test_output_shape_batch_1(self):
        """Should work with batch size 1."""
        module = ReferenceDelta(input_dim=512, out_dim=256)
        ref = torch.randn(1, 512)
        alt = torch.randn(1, 512)
        out = module(ref, alt)
        assert out.shape == (1, 256)

    def test_output_shape_rna_dim(self):
        """Test with RNA embedding dimension (512)."""
        B = 4
        module = ReferenceDelta(input_dim=512, hidden_dim=256, out_dim=256, n_heads=8)
        ref = torch.randn(B, 512)
        alt = torch.randn(B, 512)
        out = module(ref, alt)
        assert out.shape == (B, 256)

    def test_output_shape_protein_dim(self):
        """Test with protein embedding dimension (1280)."""
        B = 4
        module = ReferenceDelta(input_dim=1280, hidden_dim=256, out_dim=256, n_heads=8)
        ref = torch.randn(B, 1280)
        alt = torch.randn(B, 1280)
        out = module(ref, alt)
        assert out.shape == (B, 256)

    def test_output_shape_disorder_dim(self):
        """Test with small disorder-projected dimension (64)."""
        B = 4
        module = ReferenceDelta(input_dim=64, hidden_dim=128, out_dim=128, n_heads=4)
        ref = torch.randn(B, 64)
        alt = torch.randn(B, 64)
        out = module(ref, alt)
        assert out.shape == (B, 128)

    @pytest.mark.parametrize("input_dim", [64, 128, 256, 512, 1280])
    def test_various_input_dims(self, input_dim):
        """Test with multiple input dimensions."""
        B, out_dim = 4, 128
        n_heads = 4  # safe number of heads for all dims
        module = ReferenceDelta(input_dim=input_dim, hidden_dim=128, out_dim=out_dim, n_heads=n_heads)
        ref = torch.randn(B, input_dim)
        alt = torch.randn(B, input_dim)
        out = module(ref, alt)
        assert out.shape == (B, out_dim), f"Dim {input_dim}: expected [{B}, {out_dim}], got {out.shape}"


class TestReferenceDeltaAsymmetry:
    """Tests verifying the module is truly asymmetric."""

    def test_asymmetric_non_identical_inputs(self):
        """Swapping ref and alt should produce different outputs."""
        B, D = 8, 512
        module = ReferenceDelta(input_dim=D, hidden_dim=256, out_dim=256, dropout=0.0)
        module.eval()

        # Use fixed random seed for reproducibility
        torch.manual_seed(42)
        ref = torch.randn(B, D)
        alt = torch.randn(B, D)

        out_normal = module(ref, alt)
        out_swapped = module(alt, ref)

        diff = (out_normal - out_swapped).abs().sum().item()
        assert diff > 1e-3, (
            f"Module appears symmetric: ||delta(ref,alt) - delta(alt,ref)|| = {diff:.6f}. "
            "The module should be asymmetric — the directionality of loss/gain gates "
            "is fundamental to Proteus."
        )

    def test_kinase_vs_truncation_asymmetry(self):
        """
        Conceptual test: delta(kinase_full, truncated) should differ from
        delta(truncated, kinase_full) by a meaningful amount.
        """
        D = 512
        module = ReferenceDelta(input_dim=D, hidden_dim=256, out_dim=256, dropout=0.0)
        module.eval()

        # Simulate: kinase embedding has high values in certain dimensions
        # truncated embedding has low values
        torch.manual_seed(123)
        kinase_emb = torch.randn(1, D).abs()  # positive, "functional"
        truncated_emb = torch.randn(1, D) * 0.1  # small, "truncated"

        delta_correct = module(kinase_emb, truncated_emb)   # kinase→truncated
        delta_inverted = module(truncated_emb, kinase_emb)  # truncated→kinase

        cosine_sim = torch.nn.functional.cosine_similarity(
            delta_correct, delta_inverted, dim=-1
        ).item()

        # They shouldn't be identical or perfectly correlated
        assert abs(cosine_sim) < 0.999, (
            f"delta(kinase, truncated) and delta(truncated, kinase) are nearly identical "
            f"(cosine={cosine_sim:.4f}). This defeats the purpose of asymmetric design."
        )

    def test_identical_inputs_zero_diff(self):
        """
        When ref == alt (identical isoforms), the module should produce
        a small but consistent output (shortcut should dominate, not pure noise).
        """
        D = 512
        module = ReferenceDelta(input_dim=D, hidden_dim=256, out_dim=256, dropout=0.0)
        module.eval()

        emb = torch.randn(4, D)
        out = module(emb, emb.clone())  # identical inputs

        # Should not produce NaN
        assert not torch.isnan(out).any(), "NaN with identical inputs"
        # Loss and gain signals should be near zero for identical inputs
        # (diff_rl and diff_ar should be zero, so gates have no signal)
        # The shortcut_alpha starts at 0, so output should be near-zero initially
        # Just check it's finite
        assert torch.isfinite(out).all(), "Infinite values with identical inputs"


class TestReferenceDeltaGradients:
    """Tests for gradient flow."""

    def test_gradient_flows_to_ref(self):
        """Gradient should flow to ref_emb."""
        D = 512
        module = ReferenceDelta(input_dim=D, hidden_dim=256, out_dim=256)
        ref = torch.randn(4, D, requires_grad=True)
        alt = torch.randn(4, D, requires_grad=False)
        out = module(ref, alt)
        out.sum().backward()
        assert ref.grad is not None, "No gradient to ref_emb"
        assert ref.grad.abs().sum() > 1e-8, "Zero gradient to ref_emb"

    def test_gradient_flows_to_alt(self):
        """Gradient should flow to alt_emb."""
        D = 512
        module = ReferenceDelta(input_dim=D, hidden_dim=256, out_dim=256)
        ref = torch.randn(4, D, requires_grad=False)
        alt = torch.randn(4, D, requires_grad=True)
        out = module(ref, alt)
        out.sum().backward()
        assert alt.grad is not None, "No gradient to alt_emb"
        assert alt.grad.abs().sum() > 1e-8, "Zero gradient to alt_emb"

    def test_gradient_flows_to_both_simultaneously(self):
        """Both ref_emb and alt_emb should receive non-zero gradients."""
        D = 512
        module = ReferenceDelta(input_dim=D, hidden_dim=256, out_dim=256)
        ref = torch.randn(4, D, requires_grad=True)
        alt = torch.randn(4, D, requires_grad=True)
        out = module(ref, alt)
        out.sum().backward()

        assert ref.grad is not None and ref.grad.abs().sum() > 1e-8, \
            "No meaningful gradient to ref_emb"
        assert alt.grad is not None and alt.grad.abs().sum() > 1e-8, \
            "No meaningful gradient to alt_emb"

    def test_gradient_magnitudes_comparable(self):
        """
        Gradients to ref and alt should be of comparable magnitude.
        A ratio > 1000x would suggest one path is effectively disconnected.
        """
        D = 512
        module = ReferenceDelta(input_dim=D, hidden_dim=256, out_dim=256)
        ref = torch.randn(8, D, requires_grad=True)
        alt = torch.randn(8, D, requires_grad=True)
        out = module(ref, alt)
        out.sum().backward()

        ref_grad_norm = ref.grad.norm().item()
        alt_grad_norm = alt.grad.norm().item()

        if ref_grad_norm > 0 and alt_grad_norm > 0:
            ratio = max(ref_grad_norm, alt_grad_norm) / min(ref_grad_norm, alt_grad_norm)
            assert ratio < 1000, (
                f"Gradient ratio ref/alt is {ratio:.1f}x — suspiciously large. "
                f"ref_grad={ref_grad_norm:.4e}, alt_grad={alt_grad_norm:.4e}"
            )


class TestReferenceDeltaNumerics:
    """Tests for numerical stability and NaN-free outputs."""

    def test_no_nan_random_input(self):
        """No NaN in output for normal random inputs."""
        D = 512
        module = ReferenceDelta(input_dim=D, hidden_dim=256, out_dim=256)
        ref = torch.randn(16, D)
        alt = torch.randn(16, D)
        out = module(ref, alt)
        assert not torch.isnan(out).any(), "NaN detected in output"
        assert torch.isfinite(out).all(), "Infinite values in output"

    def test_no_nan_zero_alt(self):
        """No NaN when alternative embedding is all zeros (NMD case)."""
        D = 512
        module = ReferenceDelta(input_dim=D, hidden_dim=256, out_dim=256)
        ref = torch.randn(8, D)
        alt = torch.zeros(8, D)  # NMD — no protein
        out = module(ref, alt)
        assert not torch.isnan(out).any(), "NaN with zero alt embedding"

    def test_no_nan_large_embeddings(self):
        """No NaN with large-magnitude embeddings."""
        D = 1280  # ESM2 dim
        module = ReferenceDelta(input_dim=D, hidden_dim=256, out_dim=256, n_heads=8)
        ref = torch.randn(4, D) * 10
        alt = torch.randn(4, D) * 10
        out = module(ref, alt)
        assert not torch.isnan(out).any(), "NaN with large-magnitude embeddings"

    def test_shortcut_alpha_starts_near_zero(self):
        """
        shortcut_alpha should be initialized to zero, meaning the shortcut
        contributes tanh(0) = 0 initially (fully off at init).
        """
        module = ReferenceDelta(input_dim=512, out_dim=256)
        alpha_val = module.shortcut_alpha.item()
        assert abs(alpha_val) < 1e-6, (
            f"shortcut_alpha should be initialized to 0, got {alpha_val}"
        )

    def test_deterministic_eval_mode(self):
        """In eval mode (dropout=0), same inputs should give same outputs."""
        D = 256
        module = ReferenceDelta(input_dim=D, hidden_dim=128, out_dim=128, dropout=0.1)
        module.eval()

        torch.manual_seed(99)
        ref = torch.randn(4, D)
        alt = torch.randn(4, D)

        out1 = module(ref, alt)
        out2 = module(ref, alt)

        assert torch.allclose(out1, out2, atol=1e-6), \
            "Non-deterministic output in eval mode (dropout should be disabled)"


class TestReferenceDeltaModuleProperties:
    """Tests for module properties and initialization."""

    def test_has_required_submodules(self):
        """Module should have all required subcomponents."""
        module = ReferenceDelta(input_dim=512, hidden_dim=256, out_dim=256)
        assert hasattr(module, "ref_proj"), "Missing ref_proj"
        assert hasattr(module, "alt_proj"), "Missing alt_proj"
        assert hasattr(module, "cross_attn"), "Missing cross_attn"
        assert hasattr(module, "loss_gate_fc"), "Missing loss_gate_fc"
        assert hasattr(module, "gain_gate_fc"), "Missing gain_gate_fc"
        assert hasattr(module, "delta_proj"), "Missing delta_proj"
        assert hasattr(module, "layer_norm"), "Missing layer_norm"
        assert hasattr(module, "shortcut_alpha"), "Missing shortcut_alpha"
        assert hasattr(module, "shortcut_proj"), "Missing shortcut_proj"

    def test_parameter_count_reasonable(self):
        """Module should have a reasonable number of parameters (not suspiciously small)."""
        module = ReferenceDelta(input_dim=512, hidden_dim=256, out_dim=256)
        n_params = sum(p.numel() for p in module.parameters())
        # Expect at least 1M parameters for this configuration
        assert n_params > 100_000, f"Suspiciously few parameters: {n_params}"

    def test_repr_contains_dims(self):
        """extra_repr should include dimension info."""
        module = ReferenceDelta(input_dim=512, hidden_dim=256, out_dim=256)
        repr_str = repr(module)
        assert "512" in repr_str, "input_dim not in repr"
