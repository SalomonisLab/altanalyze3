"""CLI launcher for the scalable_viewer.

  PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 \
  /opt/homebrew/opt/python@3.11/bin/python3.11 \
    -m altanalyze3.components.visualization.scalable_viewer.run \
    --root /path/to/bundles --assets /path/to/assets --state-dir /path/to/runtime \
    [--port 8062]

--root is walked once for */*_metadata.json files carrying a "scalable_viewer" block.
--catalog reads an explicit JSON instead. Both may be given.
--assets points at the directory prepare_assets.py wrote (<assets>/<dataset id>/...).
--state-dir holds per-dataset logs; nothing is written into a bundle.

The server is scALABLE's own FastAPI application (cellHarmony/webapp/app.py) with a
bundle-backed job store. See scalable_app.py.
"""
from __future__ import annotations

import argparse
import os
import sys


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Launch the scalable_viewer server.")
    ap.add_argument("--root", default=None, help="directory tree holding precomputed bundles")
    ap.add_argument("--catalog", default=None,
                    help='JSON: {"datasets":[{"bundle_dir":..,"prefix":..}]}')
    ap.add_argument("--assets", default=None, help="asset root written by prepare_assets.py")
    ap.add_argument("--state-dir", default=None, help="runtime/log directory (default <root>/../viewer_runtime)")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8062)
    a = ap.parse_args(argv)

    if not a.root and not a.catalog:
        ap.error("give --root or --catalog")

    from . import data_api as da
    from .scalable_app import create_scalable_app
    import uvicorn

    catalog = da.build_catalog(root=a.root, catalog_file=a.catalog)
    print(f"[scalable_viewer] {len(catalog.entries)} dataset(s) in the catalog", flush=True)
    for e in catalog.entries:
        print(f"[scalable_viewer]   {e['id']}: {e['n_cells']} cells x {e['n_genes']} genes, "
              f"{e['n_states']} states, bundle={e['bundle_dir']}", flush=True)
    for err in catalog.load_errors:
        print(f"[scalable_viewer]   ERROR {err['prefix']}: {err['error']}", flush=True)
    if not catalog.entries:
        print("[scalable_viewer] no servable bundle found. Run precompute.py first.", flush=True)
        return 2

    state_dir = a.state_dir or os.path.join(os.path.dirname(os.path.abspath(a.root or ".")),
                                            "viewer_runtime")
    os.makedirs(state_dir, exist_ok=True)
    print(f"[scalable_viewer] state dir  {state_dir}", flush=True)
    print(f"[scalable_viewer] assets     {a.assets}", flush=True)

    app = create_scalable_app(catalog, state_dir=state_dir, assets_root=a.assets)
    print(f"[scalable_viewer] serving scALABLE on http://{a.host}:{a.port}", flush=True)
    uvicorn.run(app, host=a.host, port=a.port, log_level="info")
    return 0


if __name__ == "__main__":
    sys.exit(main())
