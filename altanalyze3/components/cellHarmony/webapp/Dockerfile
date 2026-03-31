FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib \
    PYTHONPATH=/app \
    CELLHARMONY_ROOT_PATH= \
    CELLHARMONY_JOB_STORAGE=/srv/cellharmony/jobs \
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

RUN mkdir -p /app/altanalyze3 /srv/cellharmony/jobs /tmp/matplotlib \
    && touch /app/altanalyze3/__init__.py

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/docs', timeout=5)"

CMD ["uvicorn", "altanalyze3.components.cellHarmony.webapp.app:app", "--host", "0.0.0.0", "--port", "8000"]
