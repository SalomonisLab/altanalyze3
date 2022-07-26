#!/bin/bash
UBUNTU_VERSION=${1:-"20.04"}
PYTHON_VERSION=${2:-"3.8.10"}
ALTANALYZE_VERSION=${3:-`git rev-parse --abbrev-ref HEAD`}
ALTANALYZE_URL=${4:-`git config --get remote.origin.url`}

WORKING_DIR=$( cd ../"$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo "Running unit tests with Python ${PYTHON_VERSION} in dockerized Ubuntu $UBUNTU_VERSION"
echo "Using AltAnalyze downloaded from ${ALTANALYZE_URL}, version ${ALTANALYZE_VERSION}"
echo "Working directory $WORKING_DIR"

docker rmi altanalyze:latest --force
docker build --no-cache --build-arg UBUNTU_VERSION=$UBUNTU_VERSION \
                        --build-arg PYTHON_VERSION=$PYTHON_VERSION \
                        --build-arg ALTANALYZE_VERSION=$ALTANALYZE_VERSION \
                        --build-arg ALTANALYZE_URL=$ALTANALYZE_URL \
                        --rm -t altanalyze:latest $WORKING_DIR

docker run --rm -it -v ${WORKING_DIR}:${WORKING_DIR} \
                    --workdir ${WORKING_DIR}/tests \
                    altanalyze:latest \
                    ${WORKING_DIR}/tests/unit_tests.sh

echo "Cleaning temporary directory $TMP_DIR"
rm -rf $TMP_DIR