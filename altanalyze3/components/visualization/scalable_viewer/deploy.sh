#!/bin/bash
# Build the image and push it to ECR. Run from anywhere.
#
#   REGISTRY=<account>.dkr.ecr.us-east-1.amazonaws.com ./deploy.sh
#   ./deploy.sh --registry <account>.dkr.ecr.us-east-1.amazonaws.com --tag 2026-09-19
#   ./deploy.sh --skip-build    # push what is already built locally
#
# REGISTRY is required and names the ECR registry to push to; REGION, IMAGE and TAG
# override the defaults below. Nothing on the target host is restarted: pull the tag
# there and run ./run.sh, which replaces the container in place.
set -euo pipefail

cd "$(dirname "$0")" || exit 1

# No default: an account's own registry is not something a public repository should
# carry. Set it in the environment, or pass --registry.
REGISTRY="${REGISTRY:-}"
REGION="${REGION:-us-east-1}"
IMAGE="${IMAGE:-scalable-viewer}"
TAG="${TAG:-latest}"
SKIP_BUILD=0

while [ $# -gt 0 ]; do
    case "$1" in
        --tag) TAG="${2:-}"; shift 2 ;;
        --tag=*) TAG="${1#*=}"; shift ;;
        --image) IMAGE="${2:-}"; shift 2 ;;
        --image=*) IMAGE="${1#*=}"; shift ;;
        --registry) REGISTRY="${2:-}"; shift 2 ;;
        --registry=*) REGISTRY="${1#*=}"; shift ;;
        --region) REGION="${2:-}"; shift 2 ;;
        --region=*) REGION="${1#*=}"; shift ;;
        --skip-build) SKIP_BUILD=1; shift ;;
        -h | --help) sed -n '2,10p' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 1 ;;
    esac
done

if [ -z "$REGISTRY" ]; then
    echo "Set REGISTRY (or pass --registry) to the ECR registry to push to." >&2
    exit 1
fi

if [ "$SKIP_BUILD" -eq 0 ]; then
    IMAGE="$IMAGE" TAG="$TAG" ./build.sh
elif ! docker image inspect "$IMAGE:$TAG" >/dev/null 2>&1; then
    echo "No local $IMAGE:$TAG to push; drop --skip-build." >&2
    exit 1
fi

# Already there on a second deploy; the error is not a failure.
aws ecr create-repository --region "$REGION" --repository-name "$IMAGE" >/dev/null 2>&1 \
    || echo "Repository $IMAGE exists or could not be created; continuing"

aws ecr get-login-password --region "$REGION" \
    | docker login --username AWS --password-stdin "$REGISTRY"

docker tag "$IMAGE:$TAG" "$REGISTRY/$IMAGE:$TAG"
docker push "$REGISTRY/$IMAGE:$TAG"

echo
echo "Pushed   $REGISTRY/$IMAGE:$TAG"
echo "Deploy   docker pull $REGISTRY/$IMAGE:$TAG"
echo "         IMAGE=$REGISTRY/$IMAGE:$TAG ./run.sh --bundles /path/to/bundles --assets /path/to/assets"
