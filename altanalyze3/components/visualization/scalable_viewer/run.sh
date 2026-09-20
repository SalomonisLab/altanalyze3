#!/bin/bash
# Serve the scalable_viewer on port 8003. Run from anywhere.
#
#   ./run.sh --bundles /path/to/bundles --assets /path/to/assets
#   ./run.sh --bundles /path/to/bundles                  # asset-driven views empty
#   ./run.sh --catalog /path/to/catalog.json --assets /path/to/assets
#
# The image holds the code only. A bundle tree (--bundles, walked for
# */*_metadata.json) or an explicit catalog JSON (--catalog) must be given, and goes
# in read-only. --state holds the per-dataset logs and is the only directory the
# container writes to; it defaults to ./viewer_runtime beside this script.
#
# The paths a bundle records go in read-only at those same paths: the directories a
# catalog names, and each bundle's own source h5ad. The viewer serves without that h5ad;
# a host that has it also gets the Gene Detail fallback that opens it.
#
# PORT, BIND, IMAGE, NAME, NETWORK and STATE override the defaults below. NETWORK is a
# second docker network to join, for reaching another container by service name. SCALABLE_ROOT_PATH,
# SCALABLE_ASSISTANT_URL, LUNGMAP_SITE_DB, LUNGMAP_SOURCE_TABLES, LUNGMAP_STUDY_IDS and
# LUNGMAP_SITE_BASE pass through when they are set in the environment; a SITE_DB or
# SOURCE_TABLES path is mounted read-only at the same path inside the container.
#
# SCALABLE_ROOT_PATH is the path prefix the viewer answers on, for a proxy that passes
# the prefix through: SCALABLE_ROOT_PATH=/scalable-viewer behind
#   ProxyPass /scalable-viewer/ http://127.0.0.1:8005/scalable-viewer/
set -euo pipefail

cd "$(dirname "$0")" || exit 1

PORT="${PORT:-8003}"
BIND="${BIND:-0.0.0.0}"
IMAGE="${IMAGE:-scalable-viewer}"
NAME="${NAME:-scalable-viewer}"
NETWORK="${NETWORK:-}"
BUNDLES="${BUNDLES:-}"
ASSETS="${ASSETS:-}"
CATALOG="${CATALOG:-}"
STATE="${STATE:-$PWD/viewer_runtime}"

while [ $# -gt 0 ]; do
    case "$1" in
        --bundles) BUNDLES="${2:-}"; shift 2 ;;
        --bundles=*) BUNDLES="${1#*=}"; shift ;;
        --assets) ASSETS="${2:-}"; shift 2 ;;
        --assets=*) ASSETS="${1#*=}"; shift ;;
        --catalog) CATALOG="${2:-}"; shift 2 ;;
        --catalog=*) CATALOG="${1#*=}"; shift ;;
        --state) STATE="${2:-}"; shift 2 ;;
        --state=*) STATE="${1#*=}"; shift ;;
        --port) PORT="${2:-}"; shift 2 ;;
        --port=*) PORT="${1#*=}"; shift ;;
        -h | --help) sed -n '2,21p' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 1 ;;
    esac
done

if [ -z "$BUNDLES" ] && [ -z "$CATALOG" ]; then
    echo "Give --bundles <dir> or --catalog <file>; the server has nothing to serve without one." >&2
    exit 1
fi

mounts=()
envs=()

if [ -n "$BUNDLES" ]; then
    if [ ! -d "$BUNDLES" ]; then
        echo "No bundle tree at $BUNDLES" >&2
        exit 1
    fi
    BUNDLES="$(cd "$BUNDLES" && pwd)"
    mounts+=(-v "$BUNDLES:/data/bundles:ro")
fi

if [ -n "$ASSETS" ]; then
    if [ ! -d "$ASSETS" ]; then
        echo "No asset root at $ASSETS" >&2
        exit 1
    fi
    ASSETS="$(cd "$ASSETS" && pwd)"
    mounts+=(-v "$ASSETS:/data/assets:ro")
fi

if [ -n "$CATALOG" ]; then
    if [ ! -f "$CATALOG" ]; then
        echo "No catalog file at $CATALOG" >&2
        exit 1
    fi
    CATALOG="$(cd "$(dirname "$CATALOG")" && pwd)/$(basename "$CATALOG")"
    # A relative bundle_dir in the catalog resolves against the file, so the file
    # goes in at its own path and the tree it names has to be reachable there too.
    mounts+=(-v "$CATALOG:$CATALOG:ro")
    envs+=(-e "VIEWER_CATALOG=$CATALOG")
    # The release directory around that file, read-only. A release is more than its
    # bundles: the integrated pseudobulk tier sits beside them, and the Regulatory
    # network and integrated pathway views read it from there. Mounting only the
    # bundles left those views answering "differentials are missing" for data that was
    # on the host all along.
    RELEASE_DIR="$(dirname "$CATALOG")"
    mounts+=(-v "$RELEASE_DIR:$RELEASE_DIR:ro")
fi

# Two paths recorded inside a bundle resolve on the host, so they go in at those same
# paths read-only.
#
# A catalog's `bundle_dir` is the dataset itself: --catalog mounts the file at its own
# path, and relative entries resolve against it, so the directories it names have to be
# reachable there.
#
# `source_h5ad` is the h5ad precompute read, which the bundle does not need: the
# expression cache is bundle-owned since 2e98942. app.py still opens it for one thing,
# the Gene Detail fallback for a gene the comparison does not carry, so mount it where
# the host has it and say what is lost where it does not.
while IFS="	" read -r kind path; do
    [ -n "$path" ] || continue
    if [ "$kind" = "refused" ]; then
        echo "Refused to mount $path: a bundle names it, but it is outside the bundle root (or is not an .h5ad file). Nothing that depends on it will load." >&2
        continue
    fi
    if [ ! -e "$path" ]; then
        case "$kind" in
            "source h5ad")
                echo "Note: $path, the h5ad a bundle was built from, is not on this host. The bundle serves without it; Gene Detail loses its fallback for a gene outside the comparison." >&2
                ;;
            *)
                echo "A bundle directory the catalog names, $path, is not on this host; that dataset will not load." >&2
                ;;
        esac
        continue
    fi
    case " ${mounts[*]} " in
        *" $path:$path:ro "*) continue ;;
    esac
    mounts+=(-v "$path:$path:ro")
done <<< "$(python3 - "${BUNDLES:-}" "${CATALOG:-}" <<'PYEOF'
import json, os, sys

root, catalog = sys.argv[1], sys.argv[2]


def read(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


bundle_dirs = []
if catalog:
    cat = read(catalog) or {}
    for entry in cat.get("datasets", []):
        bundle_dir = entry.get("bundle_dir", "")
        if bundle_dir and not os.path.isabs(bundle_dir):
            bundle_dir = os.path.join(os.path.dirname(os.path.abspath(catalog)), bundle_dir)
        if bundle_dir:
            bundle_dirs.append(os.path.abspath(bundle_dir))

emitted = set()

# A bundle is a portable artifact: it is built elsewhere and copied here, so the paths
# inside its metadata are a stranger's input, not this host's configuration. Mounting
# one unchecked let a *_metadata.json choose what the container sees - name the docker
# socket, or a home directory, and an unauthenticated internet-facing service gets it.
# Nothing is mounted now unless it sits inside a root the operator named on the command
# line, and a source h5ad must additionally be a real .h5ad file.
roots = [os.path.realpath(p) for p in ([root] if root else [])]
if catalog:
    roots.append(os.path.realpath(os.path.dirname(os.path.abspath(catalog))))


def contained(path):
    real = os.path.realpath(path)
    return any(real == r or real.startswith(r + os.sep) for r in roots)


def emit(kind, path):
    if not path or path in emitted:
        return
    if not contained(path):
        print("refused\t" + path, file=sys.stderr)
        return
    if kind == "source h5ad" and not (os.path.isfile(path) and path.endswith(".h5ad")):
        print("refused\t" + path, file=sys.stderr)
        return
    emitted.add(path)
    print(kind + "\t" + path)


# The bundle tree goes in whole at /data/bundles, so only a catalog's directories
# need a mount of their own.
for bundle_dir in bundle_dirs:
    if not root or os.path.commonpath([os.path.abspath(root), bundle_dir]) != os.path.abspath(root):
        emit("bundle directory", bundle_dir)

for tree in filter(None, [root] + bundle_dirs):
    for dirpath, _dirnames, filenames in os.walk(tree):
        for name in sorted(filenames):
            if not name.endswith("_metadata.json"):
                continue
            meta = read(os.path.join(dirpath, name))
            if not isinstance(meta, dict) or "scalable_viewer" not in meta:
                continue
            emit("source h5ad", str(meta.get("source_h5ad") or "").strip())
PYEOF
)"

mkdir -p "$STATE"
STATE="$(cd "$STATE" && pwd)"
mounts+=(-v "$STATE:/srv/scalable_viewer/state")

for var in SCALABLE_ROOT_PATH SCALABLE_ASSISTANT_URL LUNGMAP_STUDY_IDS LUNGMAP_SITE_BASE FORWARDED_ALLOW_IPS; do
    if [ -n "${!var:-}" ]; then
        envs+=(-e "$var=${!var}")
    fi
done

# These two name paths on the host. Pass the value only when the path exists, and
# mount it where the container expects to find it.
for var in LUNGMAP_SITE_DB LUNGMAP_SOURCE_TABLES; do
    path="${!var:-}"
    [ -n "$path" ] || continue
    if [ ! -e "$path" ]; then
        echo "$var is set to $path, which does not exist; skipping it." >&2
        continue
    fi
    mounts+=(-v "$path:$path:ro")
    envs+=(-e "$var=$path")
done

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
    IMAGE="$IMAGE" ./build.sh
fi

docker rm -f "$NAME" >/dev/null 2>&1 || true

# An extra network, for reaching another container by its service name - the chat
# assistant runs in the LungMAP site container, which publishes on the host's loopback
# only and so cannot be reached from here through the host. Connecting after `run`
# rather than with --network keeps the default bridge, and with it the published port.
if [ -n "$NETWORK" ] && ! docker network inspect "$NETWORK" >/dev/null 2>&1; then
    echo "No docker network named $NETWORK" >&2
    exit 1
fi

# The viewer answers the public internet and only ever reads its data. Running it as
# root bought nothing and made every other weakness worse: a file write anywhere in the
# container, and root-owned directories left on the host through the one writable mount.
# The state directory belongs to whoever runs this script, and everything else goes in
# read-only, so the container has no reason to be anyone else.
docker run -d \
    --name "$NAME" \
    --restart unless-stopped \
    --user "$(id -u):$(id -g)" \
    --cap-drop ALL \
    --security-opt no-new-privileges \
    -p "$BIND:$PORT:8003" \
    ${envs[@]+"${envs[@]}"} \
    ${mounts[@]+"${mounts[@]}"} \
    "$IMAGE" >/dev/null

if [ -n "$NETWORK" ]; then
    docker network connect "$NETWORK" "$NAME"
fi

echo "Serving  http://127.0.0.1:$PORT/  from ${BUNDLES:-$CATALOG}"
echo "Health   http://127.0.0.1:$PORT/fast/healthz"
echo "Logs     docker logs -f $NAME"
echo "Stop     docker rm -f $NAME"
