#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  cat <<'EOF'
Usage:
  Daedalus/container/run_apptainer.sh /path/to/daedalus.sif [command...]

Examples:
  Daedalus/container/run_apptainer.sh daedalus.sif
  Daedalus/container/run_apptainer.sh daedalus.sif bash
  Daedalus/container/run_apptainer.sh daedalus.sif \
    python Daedalus/phase_b/scripts/train_reference_delta_multitask.py

Environment variables:
  DAEDALUS_BIND_SRC   Source directory to bind into the container.
                      Default: current working directory
  DAEDALUS_BIND_DST   Destination mount point inside the container.
                      Default: /workspace
  DAEDALUS_NV         GPU passthrough flag. Default: 1 (enabled)
EOF
  exit 1
fi

IMAGE_PATH="$1"
shift || true

BIND_SRC="${DAEDALUS_BIND_SRC:-$(pwd)}"
BIND_DST="${DAEDALUS_BIND_DST:-/workspace}"
USE_NV="${DAEDALUS_NV:-1}"

APPTAINER_ARGS=()
if [[ "${USE_NV}" == "1" ]]; then
  APPTAINER_ARGS+=(--nv)
fi
APPTAINER_ARGS+=(--bind "${BIND_SRC}:${BIND_DST}")

if [[ $# -eq 0 ]]; then
  exec apptainer exec "${APPTAINER_ARGS[@]}" "${IMAGE_PATH}" bash
fi

exec apptainer exec "${APPTAINER_ARGS[@]}" "${IMAGE_PATH}" "$@"
