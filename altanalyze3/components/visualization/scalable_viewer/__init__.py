"""scalable_viewer - a multi-dataset, precomputed single-cell browser.

Precompute a bundle once (precompute.py), then serve many bundles from one FastAPI
process (server.py). The server never opens an h5ad.
"""
from __future__ import annotations

__all__ = ["bundle", "data_api", "server", "precompute", "run"]
__version__ = "1.0"
