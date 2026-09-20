# scalable_viewer: one FastAPI process serving precomputed bundles.
#
# The build context is altanalyze3/components, not this directory, because
# scalable_app.py imports altanalyze3.components.cellHarmony.webapp.app, which in
# turn pulls in cellHarmony/flask and several rna2* components. The image holds the
# code only; bundles and assets are mounted at run time (see run.sh).
#
#   docker build -f visualization/scalable_viewer/Dockerfile -t scalable-viewer ../..
#   ./build.sh    does exactly that
FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib \
    NUMBA_CACHE_DIR=/tmp/numba \
    CELLHARMONY_JOB_STORAGE=/tmp/cellharmony-jobs \
    PYTHONPATH=/app \
    VIEWER_ROOT=/data/bundles \
    VIEWER_ASSETS=/data/assets \
    VIEWER_CATALOG= \
    VIEWER_STATE_DIR=/srv/scalable_viewer/state \
    VIEWER_HOST=0.0.0.0 \
    VIEWER_PORT=8003 \
    FORWARDED_ALLOW_IPS=*

# FORWARDED_ALLOW_IPS: the templates build absolute URLs with Starlette's `url_for`, so
# the app has to know the scheme the browser used. uvicorn reads X-Forwarded-Proto only
# from a peer it trusts, and its default, 127.0.0.1, is never the peer here: a proxy on
# the host reaches the container as the docker gateway, so every generated URL came out
# http:// and an https page blocked its own stylesheet and scripts as mixed content.
# The container is meant to sit behind that proxy - publish it on loopback, as run.sh
# does - and the proxy sets the header itself.

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    pkg-config \
    libcairo2 \
    libcairo2-dev \
    libffi-dev \
    libgomp1 \
    libhdf5-dev \
    libjpeg62-turbo-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# The viewer serves the web app's own routes, so it needs the web app's own
# dependency set. README.md names this file as the one to install; keeping a second
# copy here would only drift from it.
COPY cellHarmony/webapp/requirements.docker.txt /tmp/requirements.docker.txt
RUN python -m pip install --upgrade pip setuptools wheel \
    && python -m pip install -r /tmp/requirements.docker.txt

COPY . /app/altanalyze3/components

# altanalyze3 itself is a package; components is copied under it, so the top level
# needs an __init__ for `python -m altanalyze3.components...` to resolve.
# Importing the web app builds a second, unused FastAPI at module scope, and that one
# creates its job directory inside site-packages; CELLHARMONY_JOB_STORAGE moves it
# somewhere writable. The viewer's own app uses the state directory regardless.
# numba caches its compiled functions beside the library source unless told
# otherwise, and site-packages belongs to root: scanpy's import died on that as soon as
# the container stopped running as root. 1777 on the directories the process writes to: run.sh starts the container as the
# invoking user rather than root, and these are created here, by root, at build time.
RUN mkdir -p /app/altanalyze3 /data/bundles /data/assets /srv/scalable_viewer/state /tmp/matplotlib /tmp/numba /tmp/cellharmony-jobs \
    && chmod 1777 /srv/scalable_viewer/state /tmp/matplotlib /tmp/numba /tmp/cellharmony-jobs \
    && touch /app/altanalyze3/__init__.py \
    && chmod +x /app/altanalyze3/components/visualization/scalable_viewer/docker-entrypoint.sh

EXPOSE 8003

# server.py stays mounted at /fast, so /fast/healthz reports the catalog without
# rendering a page. Behind a path prefix the app answers only under it, so the check
# carries SCALABLE_ROOT_PATH too.
HEALTHCHECK --interval=30s --timeout=10s --start-period=180s --retries=3 \
  CMD python -c "import os, urllib.request; urllib.request.urlopen('http://127.0.0.1:' + os.environ['VIEWER_PORT'] + os.environ.get('SCALABLE_ROOT_PATH', '') + '/fast/healthz', timeout=5)"

ENTRYPOINT ["/app/altanalyze3/components/visualization/scalable_viewer/docker-entrypoint.sh"]
