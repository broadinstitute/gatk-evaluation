#!/usr/bin/env bash

# stop script if there is an error
set -e

ADVISOR_SERVER="http://127.0.0.1:8001"
CROMWELL_SERVER="https://cromwell-v36.dsde-methods.broadinstitute.org"
STUDY_NAME="toy-gCNV-opt"
ALGORITHM="BayesianOptimization"

docker run \
    --network host \
    -v $PWD/../pipeline-optimizer/launch_study.py:/launch_study.py \
    -v $PWD/workflow.wdl:/workflow.wdl \
    -v $PWD/template.json:/template.json \
    -v $PWD/scan.json:/scan.json \
    -it us.gcr.io/broad-dsde-methods/gatk-evaluation:pipeline-optimizer \
    python2 /launch_study.py \
        --advisor_server ${ADVISOR_SERVER} \
        --cromwell_server ${CROMWELL_SERVER} \
        --study_name ${STUDY_NAME} \
        --algorithm ${ALGORITHM} \
        --workflow_wdl /workflow.wdl \
        --template_json /template.json \
        --scan_json /scan.json \
        --womtool_path /usr/local/bin/womtool/womtool-36.jar
