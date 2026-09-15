#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "$0")/../../.." && pwd)
python3 -m build --wheel --outdir "$repo/altanalyze3/deployment/containers/wheels" "$repo"
docker build -f "$repo/altanalyze3/deployment/containers/Dockerfile" -t altanalyze3-snaf:0.1.3 "$repo"
