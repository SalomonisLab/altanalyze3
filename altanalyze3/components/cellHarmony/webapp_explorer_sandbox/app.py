from __future__ import annotations

from altanalyze3.components.cellHarmony.webapp.app import create_app as create_base_app

from .config import load_config


def create_app(test_config: dict | None = None):
    cfg = load_config(test_config)
    return create_base_app(cfg)


app = create_app()
