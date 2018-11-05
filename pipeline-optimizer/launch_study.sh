#!/usr/bin/env bash

# stop script if there is an error
set -e

# parser arguments
while getopts "d:a:c:w:t:h:" option; do
	case "$option" in
		d) DOCKER="$OPTARG" ;;
		a) ADVISOR_SERVER="$OPTARG" ;;
		c) CROMWELL_SERVER="$OPTARG" ;;
		w) WORKFLOW_WDL="$OPTARG" ;;
		t) TEMPLATE_JSON="$OPTARG" ;;
		h) HYPERPARAMETERS_JSON="$OPTARG" ;;
	esac
done

if [ ! $DOCKER ] | [ ! $ADVISOR_SERVER ] | [ ! $CROMWELL_SERVER ] | [ ! $WORKFLOW_WDL ] | [ ! $TEMPLATE_JSON ] | [ ! $HYPERPARAMETERS_JSON ]; then
	echo "All arguments are required."
	exit 1
fi

echo "DOCKER: ${DOCKER}"
echo "ADVISOR_SERVER: ${ADVISOR_SERVER}"
echo "CROMWELL_SERVER: ${CROMWELL_SERVER}"
echo "WORKFLOW_WDL: ${WORKFLOW_WDL}"
echo "TEMPLATE_JSON: ${TEMPLATE_JSON}"
echo "HYPERPARAMETERS_JSON: ${HYPERPARAMETERS_JSON}"

docker run \
    --network host \
    -v $PWD/launch_study.py:/launch_study.py \
    -v $PWD/${WORKFLOW_WDL}:/workflow.wdl \
    -v $PWD/${TEMPLATE_JSON}:/template.json \
    -v $PWD/${HYPERPARAMETERS_JSON}:/hyperparameters.json \
    -it ${DOCKER} \
    python2 /launch_study.py \
        --advisor_server ${ADVISOR_SERVER} \
        --cromwell_server ${CROMWELL_SERVER} \
        --workflow_wdl /workflow.wdl \
        --template_json /template.json \
        --hyperparameters_json /hyperparameters.json \
        --womtool_path /usr/local/bin/womtool/womtool-35.jar
