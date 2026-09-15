# Build from repository root after running deployment/containers/build.sh.
ARG PYTHON_IMAGE=python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534
FROM ${PYTHON_IMAGE}
ARG TARGETARCH
LABEL org.opencontainers.image.title="AltAnalyze3 SNAF" \
      org.opencontainers.image.version="0.1.3" \
      org.opencontainers.image.source="https://github.com/SalomonisLab/altanalyze3"
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 SNAF_OFFLINE=1 \
    MPLCONFIGDIR=/tmp/matplotlib XDG_CACHE_HOME=/tmp/cache
COPY altanalyze3/deployment/containers/requirements-linux-*.txt /tmp/locks/
RUN case "${TARGETARCH}" in amd64) platform=x86_64;; arm64) platform=aarch64;; *) exit 2;; esac \
    && pip install --require-hashes -r "/tmp/locks/requirements-linux-${platform}.txt"
COPY altanalyze3/deployment/containers/wheels/*.whl /tmp/wheels/
RUN pip install --no-deps /tmp/wheels/*.whl && pip check \
    && python -m altanalyze3.components.neoantigen --help \
    && rm -rf /tmp/wheels /tmp/locks
WORKDIR /work
CMD ["altanalyze3-neo", "--help"]
