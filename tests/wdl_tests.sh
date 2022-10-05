#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

TMP_DIR="$1"
mkdir -p $TMP_DIR
echo -e '{\n  "final_workflow_outputs_dir": "TMP_DIR/outputs",\n  "final_workflow_log_dir": "TMP_DIR/wf_logs",\n  "final_call_logs_dir": "TMP_DIR/call_logs"\n}' > $TMP_DIR/options.json
sed -i "s|TMP_DIR|${TMP_DIR}|g" $TMP_DIR/options.json

$DIR/data_for_tests.sh

cd $DIR/wdl_tests

java -jar /usr/local/bin/cromwell.jar run $DIR/../wdls/altanalyze-intcount.wdl --inputs altanalyze-intcount-1.json --options $TMP_DIR/options.json
java -jar /usr/local/bin/cromwell.jar run $DIR/../wdls/altanalyze-juncount.wdl --inputs altanalyze-juncount-1.json --options $TMP_DIR/options.json
java -jar /usr/local/bin/cromwell.jar run $DIR/../wdls/altanalyze-juncount.wdl --inputs altanalyze-juncount-2.json --options $TMP_DIR/options.json
rm -rf cromwell-executions cromwell-workflow-logs

