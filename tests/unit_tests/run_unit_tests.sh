#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
pip3 install -r $DIR/test_requirements.txt

pytest --cov=altanalyze --cov-append --forked $DIR/test_example.py