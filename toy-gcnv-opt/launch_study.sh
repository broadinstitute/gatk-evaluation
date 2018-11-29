#!/usr/bin/env bash

# stop script if there is an error
set -e

ADVISOR_SERVER="127.0.0.1:8001"
CROMWELL_SERVER="https://cromwell-v34.dsde-methods.broadinstitute.org"
STUDY_NAME="toy-gCNV-opt"
ALGORITHM="BayesianOptimization"

docker run \
    --network host \
    -v $PWD/../pipeline-optimizer/launch_study.py:/launch_study.py \
    -v $PWD/workflow.wdl:/workflow.wdl \
    -v $PWD/template.json:/template.json \
    -v $PWD/scan.json:/scan.json \
    -it gatk-evaluation/pipeline-optimizer:test \
    python2 /launch_study.py \
        --advisor_server ${ADVISOR_SERVER} \
        --cromwell_server ${CROMWELL_SERVER} \
        --study_name ${STUDY_NAME} \
        --algorithm ${ALGORITHM} \
        --workflow_wdl /workflow.wdl \
        --template_json /template.json \
        --scan_json /scan.json \
        --womtool_path /usr/local/bin/womtool/womtool-36.jar
