from __future__ import annotations

from pathlib import Path

from flask import Flask

from .config import load_config
from .job_manager import JobStore
from .routes import register_routes
from .tasks import JobRunner


def create_app(test_config: dict | None = None) -> Flask:
    """
    Flask application factory for the cellHarmony-lite web experience.

    Parameters
    ----------
    test_config:
        Optional dictionary to override configuration values (useful for unit
        tests).
    """

    app = Flask(
        __name__,
        template_folder=str(Path(__file__).parent / "templates"),
        static_folder=str(Path(__file__).parent / "static"),
    )
    cfg = load_config(test_config)
    app.config.update(cfg)

    job_root = Path(app.config["JOB_STORAGE"])
    job_root.mkdir(parents=True, exist_ok=True)
    registry_path = Path(app.config["REFERENCE_REGISTRY"])
    app.job_store = JobStore(job_root)
    app.job_runner = JobRunner(
        app.job_store, registry_path, isolate_jobs=True,
        max_workers=cfg.get("JOB_WORKERS", 2),
        worker_memory_limit_gib=cfg.get("WORKER_MEMORY_LIMIT_GIB", 15),
        total_memory_limit_gib=cfg.get("TOTAL_MEMORY_LIMIT_GIB", 27),
    )

    register_routes(app)
    return app
