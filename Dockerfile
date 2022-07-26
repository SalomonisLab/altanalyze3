#####################################################################################
# Docker image with AltAnalyze (currently mostly for running tests)                 #
# with custom UBUNTU_VERSION, PYTHON_VERSION, ALTANALYZE_VERSION, ALTANALYZE_URL    #
#####################################################################################
# Build Cmd:        docker build --no-cache --rm -t altanalyze:latest .             #
#####################################################################################

# can be provided through --build-arg PARAM=value
ARG UBUNTU_VERSION="20.04"
ARG PYTHON_VERSION="3.8.10"
ARG ALTANALYZE_VERSION="master"
ARG ALTANALYZE_URL="https://github.com/SalomonisLab/altanalyze3.git"

FROM ubuntu:${UBUNTU_VERSION}
LABEL maintainer="misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

WORKDIR /tmp

ARG PYTHON_VERSION
ARG ALTANALYZE_VERSION
ARG ALTANALYZE_URL

ENV PYTHON_URL "https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz"

ENV VERSION_HTSLIB 1.15.1
ENV URL_HTSLIB "https://github.com/samtools/htslib/releases/download/${VERSION_HTSLIB}/htslib-${VERSION_HTSLIB}.tar.bz2"

ENV VERSION_SAMTOOLS 1.15.1
ENV URL_SAMTOOLS "https://github.com/samtools/samtools/releases/download/${VERSION_SAMTOOLS}/samtools-${VERSION_SAMTOOLS}.tar.bz2"


RUN echo "Installing dependencies" && \
    apt-get update && \
    apt-get install -y gcc build-essential git wget curl zlib1g-dev libffi-dev libssl-dev libsqlite3-dev libbz2-dev liblzma-dev libncurses5-dev && \
    echo "Installing Python" && \
    wget ${PYTHON_URL} && \
    tar xzf Python-${PYTHON_VERSION}.tgz && \
    cd Python-${PYTHON_VERSION} && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    pip3 install -U pip && \
    echo "Installing htslib" && \
    wget -q -O - $URL_HTSLIB | tar -jxv && \
    cd htslib-${VERSION_HTSLIB} && \
    ./configure --prefix=/usr/local/ && \
    make -j 4 && \
    make install && \
    cd .. && \
    echo "Installing samtools" && \
    wget -q -O - $URL_SAMTOOLS | tar -jxv && \
    cd samtools-${VERSION_SAMTOOLS} && \
    ./configure --prefix=/usr/local/ && \
    make -j 4 && \
    make install && \
    cd .. && \
    echo "Installing AltAnalyze" && \
    git clone ${ALTANALYZE_URL} && \
    cd altanalyze3 && \
    git checkout ${ALTANALYZE_VERSION} && \
    pip3 install . && \
    cd .. && \
    echo "Cleaning up" && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true