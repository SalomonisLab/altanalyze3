# scALABLE-discover: its own image and container, separate from scALABLE-web.
# Build context: altanalyze3/components (see docker-compose.discover.yml). The Python
# requirements are scALABLE-web's list, cellHarmony/webapp/requirements.docker.txt, which
# carries python-louvain for ICGS3; altanalyze3/components/.dockerignore keeps job folders
# and per-module virtual environments out of the context.
FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib \
    PYTHONPATH=/app \
    SCALABLE_DISCOVER_ROOT_PATH= \
    SCALABLE_DISCOVER_JOB_STORAGE=/srv/scalable-discover/jobs \
    SCALABLE_DISCOVER_EXPORTS=minimal \
    CELLHARMONY_JOB_WORKERS=2

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

COPY cellHarmony/webapp/requirements.docker.txt /tmp/requirements.docker.txt
RUN python -m pip install --upgrade pip setuptools wheel \
    && python -m pip install -r /tmp/requirements.docker.txt

COPY . /app/altanalyze3/components

RUN python /app/altanalyze3/components/fastComm/verify_packaged_resources.py

RUN mkdir -p /app/altanalyze3 /srv/scalable-discover/jobs /tmp/matplotlib \
    && touch /app/altanalyze3/__init__.py

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/docs', timeout=5)"

CMD ["uvicorn", "altanalyze3.components.cellHarmony.scalable_discover.app:app", "--host", "0.0.0.0", "--port", "8000", "--proxy-headers", "--forwarded-allow-ips", "*"]
