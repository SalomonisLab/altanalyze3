"""
Device management utilities for Proteus.
"""

from __future__ import annotations

from typing import Any, Dict

import torch


def get_device(preference: str = "auto") -> torch.device:
    """
    Select compute device based on availability.

    Parameters
    ----------
    preference : str
        "auto" selects: CUDA > MPS > CPU.
        "cuda" forces CUDA (raises if unavailable).
        "mps" forces MPS.
        "cpu" forces CPU.

    Returns
    -------
    torch.device
    """
    preference = preference.lower().strip()

    if preference == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")

    if preference == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA requested but not available.")
        return torch.device("cuda")

    if preference == "mps":
        if not (hasattr(torch.backends, "mps") and torch.backends.mps.is_available()):
            raise RuntimeError("MPS requested but not available.")
        return torch.device("mps")

    return torch.device("cpu")


def get_device_info() -> Dict[str, Any]:
    """
    Return information about the available compute device.

    Returns
    -------
    dict
        Keys: device_type, device_name, memory_gb (GPU only),
        compute_capability (CUDA only), n_gpus, cuda_version.
    """
    info: Dict[str, Any] = {
        "device_type": "cpu",
        "device_name": "CPU",
        "memory_gb": None,
        "compute_capability": None,
        "n_gpus": 0,
        "cuda_version": None,
    }

    if torch.cuda.is_available():
        info["device_type"] = "cuda"
        info["n_gpus"] = torch.cuda.device_count()
        info["cuda_version"] = torch.version.cuda

        try:
            props = torch.cuda.get_device_properties(0)
            info["device_name"] = props.name
            info["memory_gb"] = round(props.total_memory / (1024 ** 3), 2)
            info["compute_capability"] = (
                f"{props.major}.{props.minor}"
            )
        except Exception:
            pass

    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        info["device_type"] = "mps"
        info["device_name"] = "Apple Silicon MPS"

    return info
