#!/usr/bin/env python3
"""Launch the flow_viewer server. Entry point, in the shape of scalable_viewer.run."""
import argparse


def main():
    ap = argparse.ArgumentParser(description="Launch the rna2flow interactive viewer.")
    ap.add_argument("--bundle", required=True, help="directory written by precompute_interactive.py")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8085)
    a = ap.parse_args()
    from .app import create_app
    app = create_app(a.bundle)
    print("rna2flow viewer: http://%s:%d  (bundle %s)" % (a.host, a.port, a.bundle))
    app.run(host=a.host, port=a.port, debug=False, threaded=True)


if __name__ == "__main__":
    main()
