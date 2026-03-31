from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Optional


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_JOB_STORAGE = BASE_DIR / "jobs"
DEFAULT_REFERENCE_REGISTRY = BASE_DIR.parent / "flask" / "reference_config.json"


def load_config(overrides: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    cfg: Dict[str, Any] = {
        "APP_TITLE": "cellHarmony-lite Portal",
        "MAX_CONTENT_LENGTH": 1024 * 1024 * 1024,
        "ROOT_PATH": os.getenv("CELLHARMONY_ROOT_PATH", "").strip(),
        "JOB_STORAGE": os.getenv("CELLHARMONY_JOB_STORAGE", str(DEFAULT_JOB_STORAGE)),
        "REFERENCE_REGISTRY": os.getenv("CELLHARMONY_REFERENCE_REGISTRY", str(DEFAULT_REFERENCE_REGISTRY)),
        "ALLOWED_EXTENSIONS": {"h5", "h5ad"},
        "MAX_FILES_PER_JOB": int(os.getenv("CELLHARMONY_MAX_FILES", "7")),
        "JOB_WORKERS": int(os.getenv("CELLHARMONY_JOB_WORKERS", "2")),
    }
    if overrides:
        cfg.update(overrides)
    return cfg
