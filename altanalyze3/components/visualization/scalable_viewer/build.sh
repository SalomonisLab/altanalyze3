#!/bin/bash
# Build the scalable_viewer image. Run from anywhere.
#
#   ./build.sh
#
# The context is altanalyze3/components, two levels up, because the viewer imports
# the cellHarmony web app. IMAGE and TAG override the name.
set -euo pipefail

cd "$(dirname "$0")" || exit 1

IMAGE="${IMAGE:-scalable-viewer}"
TAG="${TAG:-latest}"

docker build \
    -f Dockerfile \
    -t "$IMAGE:$TAG" \
    ../..

# The import closure is the thing that breaks silently when a component is left
# out of the context, and it breaks at request time rather than at start up. Ask
# the image to import the app factory now.
docker run --rm --entrypoint python "$IMAGE:$TAG" -c \
    'from altanalyze3.components.visualization.scalable_viewer.scalable_app import create_scalable_app'

size="$(docker image inspect -f '{{.Size}}' "$IMAGE:$TAG" | awk '{printf "%.0f MB", $1/1048576}')"
echo "$IMAGE:$TAG imports the viewer app, $size"
echo "Run it  ./run.sh --bundles /path/to/bundles --assets /path/to/assets"
