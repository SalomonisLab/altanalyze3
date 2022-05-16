#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
pip3 install -r $DIR/test_requirements.txt

$DIR/../prepare_data_for_tests.sh

pytest --cov=altanalyze --cov-append --forked $DIR/test_helpers.py
pytest --cov=altanalyze --cov-append --forked $DIR/test_data.py
pytest --cov=altanalyze --cov-append --forked $DIR/test_io.py
