#!/bin/sh
# Turn the container's environment into run.py's flags. Arguments given to
# `docker run` replace them entirely, so
#   docker run ... scalable-viewer --root /data/bundles --port 9000
# works, and so does `docker run ... scalable-viewer python -c ...` for a shell.
set -eu

if [ "$#" -gt 0 ]; then
    case "$1" in
        -*) ;;                       # flags for run.py, fall through
        *) exec "$@" ;;              # a command of its own
    esac
else
    set -- --host "$VIEWER_HOST" --port "$VIEWER_PORT" --state-dir "$VIEWER_STATE_DIR"

    if [ -n "${VIEWER_CATALOG:-}" ]; then
        set -- "$@" --catalog "$VIEWER_CATALOG"
    fi
    # The image creates /data/bundles, so existence proves nothing; content does.
    if [ -n "${VIEWER_ROOT:-}" ] && [ -d "$VIEWER_ROOT" ] && [ -n "$(ls -A "$VIEWER_ROOT" 2>/dev/null)" ]; then
        set -- "$@" --root "$VIEWER_ROOT"
    fi
    if [ -n "${VIEWER_ASSETS:-}" ] && [ -d "$VIEWER_ASSETS" ]; then
        set -- "$@" --assets "$VIEWER_ASSETS"
    fi

    # run.py exits 2 on an empty catalog, which reads as a crash loop. Say what is
    # actually missing instead.
    case " $* " in
        *" --root "* | *" --catalog "*) ;;
        *)
            echo "No bundles: mount a bundle tree at $VIEWER_ROOT, or set VIEWER_CATALOG." >&2
            exit 2
            ;;
    esac
fi

mkdir -p "$VIEWER_STATE_DIR"

exec python -m altanalyze3.components.visualization.scalable_viewer.run "$@"
