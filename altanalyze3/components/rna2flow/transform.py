"""Flow-cytometry value transforms.

The Logicle transform is the step that makes a flow panel comparable to a CITE-seq panel: a
raw FCS channel spans ~1e6 with compensation-induced negatives, and FlowJo clusters the
Logicle-transformed values, not the raw ones.

pyInfinityFlow is the reference implementation and a declared dependency of rna2flow (see
requirements.txt). `logicle_anndata` delegates to it so the published numbers are reproduced
exactly. `logicle` is a vectorized path over a plain array for the case where the per-channel
T/W/M/A are already known; `tests/test_logicle_parity.py` holds it to the reference.
"""
from __future__ import annotations

import numpy as np

# Standard external dependency. pyInfinityFlow owns the Logicle definition and the FCS keyword
# handling that supplies T, W, M and A; rna2flow does not re-derive either.
from pyInfinityFlow.InfinityFlow_Utilities import (  # noqa: F401
    anndata_to_df,
    apply_logicle_to_anndata,
    move_features_to_silent,
    read_fcs_into_anndata,
)
from pyInfinityFlow.Transformations import apply_logicle, apply_inverse_logicle  # noqa: F401

__all__ = [
    "read_fcs_anndata",
    "logicle_anndata",
    "logicle",
    "inverse_logicle",
    "anndata_to_df",
    "move_features_to_silent",
]


def read_fcs_anndata(fcs_path: str, obs_prefix: str = "", batch_key: str = ""):
    """FCS -> AnnData carrying the per-channel Logicle parameters in .var.

    Thin pass-through to pyInfinityFlow so .var keeps T/W/M/A; `io.read_fcs` is the faster
    reader for the matrix alone and does not populate those.
    """
    return read_fcs_into_anndata(fcs_path, obs_prefix=obs_prefix, batch_key=batch_key)


def logicle_anndata(adata, in_place: bool = True):
    """Apply the published Logicle transform, using each channel's own T/W/M/A from .var."""
    return apply_logicle_to_anndata(adata, in_place=in_place)


def logicle(x, T: float = 3000000, W: float = 0, M: float = 3, A: float = 1) -> np.ndarray:
    """Logicle over an array with explicit parameters. Same defaults as pyInfinityFlow."""
    return np.asarray(apply_logicle(np.asarray(x, dtype=np.float64), T=T, W=W, M=M, A=A))


def inverse_logicle(x, T: float = 3000000, W: float = 0, M: float = 3, A: float = 1) -> np.ndarray:
    return np.asarray(apply_inverse_logicle(np.asarray(x, dtype=np.float64), T=T, W=W, M=M, A=A))
