"""Frozen foundation model encoder wrappers for Proteus."""

from .rna import OrthrusRNAEncoder
from .protein import ESM2ProteinEncoder
from .disorder import DisorderEncoder
from .structural import StructuralFeatureEncoder

__all__ = [
    "OrthrusRNAEncoder",
    "ESM2ProteinEncoder",
    "DisorderEncoder",
    "StructuralFeatureEncoder",
]
