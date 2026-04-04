from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_JOB_STORAGE = BASE_DIR / "jobs"


def load_config(overrides: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    cfg: Dict[str, Any] = {
        "APP_TITLE": "cellHarmony Explorer Sandbox",
        "JOB_STORAGE": str(DEFAULT_JOB_STORAGE),
        "TEMPLATE_DIR": str(BASE_DIR / "templates"),
        "STATIC_DIR": str(BASE_DIR / "static"),
        "INDEX_TEMPLATE": "index.html",
    }
    if overrides:
        cfg.update(overrides)
    return cfg
