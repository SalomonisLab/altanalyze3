from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_JOB_STORAGE = BASE_DIR / "jobs"
DEFAULT_REFERENCE_REGISTRY = BASE_DIR / "reference_config.json"


def load_config(overrides: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Produce a configuration dictionary for the Flask application.

    Parameters
    ----------
    overrides:
        Optional dictionary with keys that should replace the defaults. This is
        mainly intended for unit tests.
    """

    cfg: Dict[str, Any] = {
        "MAX_CONTENT_LENGTH": 1024 * 1024 * 1024,  # 1 GB uploads
        "JOB_STORAGE": str(DEFAULT_JOB_STORAGE),
        "REFERENCE_REGISTRY": str(DEFAULT_REFERENCE_REGISTRY),
        "ALLOWED_EXTENSIONS": {"h5", "h5ad"},
        "MAX_FILES_PER_JOB": 7,
        "SECRET_KEY": "cellharmony-lite-demo",
    }
    if overrides:
        cfg.update(overrides)
    return cfg
