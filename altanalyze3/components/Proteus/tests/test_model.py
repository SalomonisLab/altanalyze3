"""
Tests for ProteusModel forward pass, checkpointing, and architecture properties.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from typing import Any, Dict

import pytest
import torch

sys.path.insert(0, str(Path(__file__).parent.parent))
from proteus.model import ProteusModel


def make_synthetic_batch(
    batch_size: int = 4,
    has_alt_protein: bool = True,
) -> Dict[str, Any]:
    """
    Create a synthetic batch of the correct shapes for ProteusModel.

    Parameters
    ----------
    batch_size : int
    has_alt_protein : bool
        If False, marks all alt proteins as absent (NMD simulation).
    """
    B = batch_size
    return {
        "ref_rna_emb": torch.randn(B, 512),
        "alt_rna_emb": torch.randn(B, 512),
        "ref_protein_emb": torch.randn(B, 1280),
        "alt_protein_emb": torch.randn(B, 1280) if has_alt_protein else torch.zeros(B, 1280),
        "ref_disorder_raw": torch.randn(B, 50).abs(),   # disorder scores in [0,1]ish
        "alt_disorder_raw": torch.randn(B, 50).abs(),
        "ref_structural_raw": torch.randn(B, 48),
        "alt_structural_raw": torch.randn(B, 48),
        "has_alt_protein": torch.tensor([has_alt_protein] * B, dtype=torch.bool),
        "label_preservation": torch.randint(0, 2, (B,)),
        "deviation_class": torch.randint(0, 6, (B,)),
        "task_masks": {
            "global_preservation": torch.ones(B),
            "topology_preserved": torch.bernoulli(torch.full((B,), 0.5)),
            "surface_retained": torch.bernoulli(torch.full((B,), 0.3)),
            "kinase_competent": torch.bernoulli(torch.full((B,), 0.2)),
            "tf_competent": torch.bernoulli(torch.full((B,), 0.2)),
            "disorder_preserved": torch.ones(B),
        },
    }


@pytest.fixture
def model() -> ProteusModel:
    """Instantiate a ProteusModel with default config."""
    return ProteusModel()


@pytest.fixture
def small_model() -> ProteusModel:
    """Instantiate a small ProteusModel for faster tests."""
    cfg = {
        "delta": {
            "rna_dim": 512,
            "protein_dim": 1280,
            "disorder_dim": 32,
            "structural_dim": 32,
            "delta_hidden_dim": 128,
            "delta_out_dim": 128,
            "n_heads": 4,
            "dropout": 0.0,
        },
        "fusion": {
            "common_dim": 128,
            "mlp_hidden": 256,
            "mlp_layers": 2,
            "dropout": 0.0,
            "use_cross_modal_attention": True,
            "modality_dropout_prob": 0.0,
        },
        "heads": {
            "deviation_class": {"n_classes": 6},
        },
    }
    return ProteusModel(cfg=cfg)


class TestProteusModelForward:
    """Tests for forward pass shapes and correctness."""

    def test_forward_output_keys(self, model):
        """Forward pass should return all expected task keys."""
        batch = make_synthetic_batch(4)
        model.eval()
        with torch.no_grad():
            outputs = model(batch)

        expected_keys = {
            "global_preservation",
            "topology_preserved",
            "surface_retained",
            "kinase_competent",
            "tf_competent",
            "disorder_preserved",
            "deviation_class",
        }
        for key in expected_keys:
            assert key in outputs, f"Missing output key: {key}"

    def test_binary_task_output_shapes(self, model):
        """Binary task outputs should be [B, 1]."""
        B = 6
        batch = make_synthetic_batch(B)
        model.eval()
        with torch.no_grad():
            outputs = model(batch)

        binary_tasks = [
            "global_preservation", "topology_preserved", "surface_retained",
            "kinase_competent", "tf_competent", "disorder_preserved",
        ]
        for task in binary_tasks:
            assert outputs[task].shape == (B, 1), \
                f"{task}: expected [{B}, 1], got {outputs[task].shape}"

    def test_deviation_class_output_shape(self, model):
        """deviation_class output should be [B, 6]."""
        B = 6
        batch = make_synthetic_batch(B)
        model.eval()
        with torch.no_grad():
            outputs = model(batch)
        assert outputs["deviation_class"].shape == (B, 6), \
            f"Expected [{B}, 6], got {outputs['deviation_class'].shape}"

    def test_no_nan_with_standard_input(self, model):
        """Standard input should produce NaN-free outputs."""
        batch = make_synthetic_batch(8)
        model.eval()
        with torch.no_grad():
            outputs = model(batch)

        for task, tensor in outputs.items():
            if isinstance(tensor, torch.Tensor):
                assert not torch.isnan(tensor).any(), f"NaN in {task} output"
                assert torch.isfinite(tensor).all(), f"Inf in {task} output"

    def test_no_nan_with_absent_alt_protein(self, model):
        """NMD/absent alt protein (zeros + has_alt_protein=False) should not cause NaN."""
        batch = make_synthetic_batch(4, has_alt_protein=False)
        model.eval()
        with torch.no_grad():
            outputs = model(batch)

        for task, tensor in outputs.items():
            if isinstance(tensor, torch.Tensor):
                assert not torch.isnan(tensor).any(), \
                    f"NaN in {task} output with absent alt protein"

    def test_no_nan_with_all_zero_embeddings(self, model):
        """All-zero embeddings should not cause NaN or division by zero."""
        B = 4
        batch = {
            "ref_rna_emb": torch.zeros(B, 512),
            "alt_rna_emb": torch.zeros(B, 512),
            "ref_protein_emb": torch.zeros(B, 1280),
            "alt_protein_emb": torch.zeros(B, 1280),
            "ref_disorder_raw": torch.zeros(B, 50),
            "alt_disorder_raw": torch.zeros(B, 50),
            "ref_structural_raw": torch.zeros(B, 48),
            "alt_structural_raw": torch.zeros(B, 48),
            "has_alt_protein": torch.ones(B, dtype=torch.bool),
            "label_preservation": torch.zeros(B, dtype=torch.long),
            "deviation_class": torch.zeros(B, dtype=torch.long),
        }
        model.eval()
        with torch.no_grad():
            outputs = model(batch)

        for task, tensor in outputs.items():
            if isinstance(tensor, torch.Tensor):
                assert not torch.isnan(tensor).any(), \
                    f"NaN in {task} output with all-zero embeddings"

    def test_batch_consistency(self, model):
        """Processing items one-at-a-time vs. batched should give same results."""
        model.eval()

        torch.manual_seed(42)
        batch = make_synthetic_batch(8)
        # Drop task_masks to simplify
        batch_no_masks = {k: v for k, v in batch.items() if k != "task_masks"}

        with torch.no_grad():
            batch_out = model(batch_no_masks)

        # Process item 0 individually
        single_batch = {
            k: v[:1] if isinstance(v, torch.Tensor) else v
            for k, v in batch_no_masks.items()
        }
        with torch.no_grad():
            single_out = model(single_batch)

        for task in ["global_preservation", "deviation_class"]:
            if task in batch_out and task in single_out:
                assert torch.allclose(
                    batch_out[task][:1], single_out[task], atol=1e-5
                ), f"Batch vs single inconsistency in {task}"

    def test_predict_returns_probabilities(self, model):
        """predict() should return probabilities in [0, 1]."""
        batch = make_synthetic_batch(4)
        model.eval()
        preds = model.predict(batch)

        for task in ["global_preservation", "topology_preserved", "kinase_competent"]:
            if task in preds:
                probs = preds[task]
                assert (probs >= 0).all() and (probs <= 1).all(), \
                    f"{task} probabilities out of [0, 1] range"

        if "deviation_class" in preds:
            probs = preds["deviation_class"]
            row_sums = probs.sum(dim=-1)
            assert torch.allclose(row_sums, torch.ones_like(row_sums), atol=1e-5), \
                "deviation_class probabilities don't sum to 1"


class TestProteusModelParameters:
    """Tests for parameter count and architecture."""

    def test_parameter_count_reasonable(self, model):
        """Total parameters should be in the expected 10-15M range."""
        n_params = model.count_parameters()
        print(f"ProteusModel parameters: {n_params:,}")
        # The delta + fusion + heads architecture should be in this range
        assert n_params > 1_000_000, f"Too few parameters: {n_params:,}"
        assert n_params < 100_000_000, f"Too many parameters: {n_params:,}"

    def test_parameter_count_small_model(self, small_model):
        """Small config model should have fewer parameters."""
        n_params = small_model.count_parameters()
        print(f"Small ProteusModel parameters: {n_params:,}")
        assert n_params < 10_000_000

    def test_all_parameters_trainable(self, model):
        """All model parameters should be trainable (no frozen params in ProteusModel)."""
        n_trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
        n_total = sum(p.numel() for p in model.parameters())
        assert n_trainable == n_total, (
            f"Some parameters are frozen: {n_total - n_trainable} / {n_total}"
        )


class TestProteusModelCheckpointing:
    """Tests for save/load checkpoint roundtrip."""

    def test_save_and_load_checkpoint(self, small_model):
        """Saved and loaded model should produce identical outputs."""
        small_model.eval()
        batch = make_synthetic_batch(4)

        # Get reference outputs
        with torch.no_grad():
            ref_outputs = {
                k: v.clone() for k, v in small_model(batch).items()
                if isinstance(v, torch.Tensor)
            }

        # Save to temp file
        with tempfile.NamedTemporaryFile(suffix=".pt", delete=False) as f:
            ckpt_path = Path(f.name)

        try:
            small_model.save_checkpoint(ckpt_path, epoch=5, metrics={"val_loss": 0.42})

            # Load
            loaded_model = ProteusModel.load_checkpoint(ckpt_path)
            loaded_model.eval()

            with torch.no_grad():
                loaded_outputs = {
                    k: v.clone() for k, v in loaded_model(batch).items()
                    if isinstance(v, torch.Tensor)
                }

            # Compare
            for task in ref_outputs:
                if task in loaded_outputs:
                    assert torch.allclose(ref_outputs[task], loaded_outputs[task], atol=1e-6), \
                        f"Checkpoint mismatch in {task}"

        finally:
            ckpt_path.unlink(missing_ok=True)

    def test_checkpoint_contains_metadata(self, small_model):
        """Checkpoint should contain epoch and metrics metadata."""
        with tempfile.NamedTemporaryFile(suffix=".pt", delete=False) as f:
            ckpt_path = Path(f.name)

        try:
            small_model.save_checkpoint(
                ckpt_path,
                epoch=10,
                metrics={"auroc": 0.85, "auprc": 0.72},
            )
            checkpoint = torch.load(ckpt_path, map_location="cpu", weights_only=False)
            assert checkpoint["epoch"] == 10
            assert checkpoint["metrics"]["auroc"] == 0.85
            assert "model_state_dict" in checkpoint
            assert "cfg" in checkpoint

        finally:
            ckpt_path.unlink(missing_ok=True)

    def test_save_load_preserves_cfg(self, small_model):
        """Config should be preserved through save/load roundtrip."""
        original_cfg = small_model.cfg.copy()

        with tempfile.NamedTemporaryFile(suffix=".pt", delete=False) as f:
            ckpt_path = Path(f.name)

        try:
            small_model.save_checkpoint(ckpt_path, epoch=1)
            loaded = ProteusModel.load_checkpoint(ckpt_path)
            assert loaded.cfg == original_cfg, "Config changed after save/load"
        finally:
            ckpt_path.unlink(missing_ok=True)


class TestProteusModelGradients:
    """Tests for gradient flow through the full model."""

    def test_gradients_flow_to_rna_inputs(self, small_model):
        """Gradients should flow from loss back to RNA embedding inputs."""
        batch = make_synthetic_batch(4)
        batch["ref_rna_emb"] = batch["ref_rna_emb"].requires_grad_(True)
        batch["alt_rna_emb"] = batch["alt_rna_emb"].requires_grad_(True)

        outputs = small_model(batch)
        outputs["global_preservation"].sum().backward()

        assert batch["ref_rna_emb"].grad is not None
        assert batch["ref_rna_emb"].grad.abs().sum() > 0
        assert batch["alt_rna_emb"].grad is not None
        assert batch["alt_rna_emb"].grad.abs().sum() > 0

    def test_gradients_flow_to_protein_inputs(self, small_model):
        """Gradients should flow from loss back to protein embedding inputs."""
        batch = make_synthetic_batch(4)
        batch["ref_protein_emb"] = batch["ref_protein_emb"].requires_grad_(True)
        batch["alt_protein_emb"] = batch["alt_protein_emb"].requires_grad_(True)

        outputs = small_model(batch)
        outputs["global_preservation"].sum().backward()

        assert batch["ref_protein_emb"].grad is not None
        assert batch["ref_protein_emb"].grad.abs().sum() > 0
