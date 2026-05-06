"""Daedalus model classes and shared training-objective utilities."""

from .deepimmuno_style import (
    DeepImmunoStyleCNN,
    DeepImmunoStyleCNNAutoencoder,
    ResidualMLPAutoencoder,
    choose_device,
)
from .residual_mlp_variants import ResidualMLPAutoencoderVariant
from .family_multitask_net import FamilyMultitaskNet

__all__ = [
    "DeepImmunoStyleCNN",
    "DeepImmunoStyleCNNAutoencoder",
    "ResidualMLPAutoencoder",
    "ResidualMLPAutoencoderVariant",
    "FamilyMultitaskNet",
    "choose_device",
]
