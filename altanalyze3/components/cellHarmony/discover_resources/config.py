"""Portable resource locations; no external checkout import at deployment."""
import os
from pathlib import Path

def cellref2_dir():
    return Path(os.environ.get("SCALABLE_DISCOVER_RESOURCES", Path(__file__).parents[1]/"resources/discover"))

def state_dir():
    path=Path(os.environ.get("SCALABLE_INTERPRET_CACHE", Path(__file__).parents[1]/"webapp/jobs/.interpret-cache"))
    path.mkdir(parents=True,exist_ok=True)
    return path
