#!/bin/bash
UBUNTU_VERSION=${1:-"20.04"}
PYTHON_VERSION=${2:-"3.8.10"}
ALTANALYZE_VERSION=${3:-`git rev-parse --abbrev-ref HEAD`}

WORKING_DIR=$( cd ../"$( dirname "${BASH_SOURCE[0]}" )" && pwd )
TMP_DIR="${WORKING_DIR}/tmp"

echo "Running unit tests with Python ${PYTHON_VERSION} in dockerized Ubuntu $UBUNTU_VERSION"
echo "Using AltAnalyze (${ALTANALYZE_VERSION})"
echo "Working directory $WORKING_DIR"
echo "Temporary directory $TMP_DIR"

mkdir $TMP_DIR
docker rmi altanalyze:latest --force
docker build --no-cache --build-arg UBUNTU_VERSION=$UBUNTU_VERSION \
                        --build-arg PYTHON_VERSION=$PYTHON_VERSION \
                        --build-arg ALTANALYZE_VERSION=$CWL_AIRFLOW_VERSION \
                        --rm -t altanalyze:latest $WORKING_DIR

docker run --rm -it -v ${WORKING_DIR}:${WORKING_DIR} \
                    --workdir ${WORKING_DIR} \
                    --env TMPDIR=${TMP_DIR} \
                    altanalyze:latest \
                    ${WORKING_DIR}/test/unit_tests.sh

echo "Cleaning temporary directory $TMP_DIR"
rm -rf $TMP_DIR