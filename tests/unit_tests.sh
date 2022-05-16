#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
pip3 install -r $DIR/unit_tests/requirements.txt

$DIR/data_for_tests.sh

pytest --cov=altanalyze3 --cov-append --cov-report xml:cov.xml --forked $DIR/unit_tests/test_helpers.py
pytest --cov=altanalyze3 --cov-append --cov-report xml:cov.xml --forked $DIR/unit_tests/test_data.py
pytest --cov=altanalyze3 --cov-append --cov-report xml:cov.xml --forked $DIR/unit_tests/test_io.py
pytest --cov=altanalyze3 --cov-append --cov-report xml:cov.xml --forked $DIR/unit_tests/test_parser.py
pytest --cov=altanalyze3 --cov-append --cov-report xml:cov.xml --forked $DIR/unit_tests/test_intron_count.py
pytest --cov=altanalyze3 --cov-append --cov-report xml:cov.xml --forked $DIR/unit_tests/test_junction_count.py
