#!/usr/bin/env bash
UBUNTU_VERSION=${1:-"20.04"}
CROMWELL_VERSION=${2:-"84"}

WORKING_DIR=$( cd ../"$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo "Running wdl tests with Cromwell ${CROMWELL_VERSION} in dockerized Ubuntu $UBUNTU_VERSION"
echo "Working directory $WORKING_DIR"

docker rmi cromwell:latest --force
docker build --no-cache --build-arg UBUNTU_VERSION=$UBUNTU_VERSION \
                        --build-arg CROMWELL_VERSION=$CROMWELL_VERSION \
                        --rm -t cromwell:latest -f ${WORKING_DIR}/tests/wdl_tests/cromwell-Dockerfile .

docker run --rm -it -v /var/run/docker.sock:/var/run/docker.sock \
                    -v ${WORKING_DIR}:${WORKING_DIR} \
                    --workdir ${WORKING_DIR}/tests \
                    cromwell:latest \
                    ${WORKING_DIR}/tests/wdl_tests.sh \
                    ${WORKING_DIR}/tests/data/wdl_tmp/tmp