"""
Proteus: Reference-conditioned isoform function prediction model.

Predicts whether an alternative RNA isoform preserves the core function
of its reference isoform across seven biological dimensions.
"""

from .model import ProteusModel
from .pretrain_model import ProteusPretrainModel

__version__ = "0.1.0"
__all__ = ["ProteusModel", "ProteusPretrainModel"]
