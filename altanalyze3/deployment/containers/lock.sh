#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
# uv is build tooling only. CPU PyTorch avoids bundling CUDA libraries.
for arch in x86_64 aarch64; do
  uv pip compile requirements.in --python-version 3.11 \
    --python-platform "${arch}-manylinux_2_28" \
    --extra-index-url https://download.pytorch.org/whl/cpu \
    --index-strategy unsafe-best-match --generate-hashes --emit-index-url --no-emit-package altanalyze3 \
    --output-file "requirements-linux-${arch}.txt"
done
