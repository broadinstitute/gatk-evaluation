#!/usr/bin/env bash

# stop script if there is an error
set -e

# parser arguments
while getopts "d:a:c:n:l:w:t:s:" option; do
	case "$option" in
		d) DOCKER="$OPTARG" ;;
		a) ADVISOR_SERVER="$OPTARG" ;;
		c) CROMWELL_SERVER="$OPTARG" ;;
        n) STUDY_NAME="$OPTARG" ;;
        l) ALGORITHM="$OPTARG" ;;
		w) WORKFLOW_WDL="$OPTARG" ;;
		t) TEMPLATE_JSON="$OPTARG" ;;
		s) SCAN_JSON="$OPTARG" ;;
	esac
done

if [ ! $DOCKER ] | [ ! $ADVISOR_SERVER ] | [ ! $CROMWELL_SERVER ] | [ ! $STUDY_NAME ] | [ ! $ALGORITHM ] | [ ! $WORKFLOW_WDL ] | [ ! $TEMPLATE_JSON ] | [ ! $SCAN_JSON ]; then
	echo "All arguments are required."
	exit 1
fi

echo "DOCKER: ${DOCKER}"
echo "ADVISOR_SERVER: ${ADVISOR_SERVER}"
echo "CROMWELL_SERVER: ${CROMWELL_SERVER}"
echo "STUDY_NAME: ${STUDY_NAME}"
echo "ALGORITHM: ${ALGORITHM}"
echo "WORKFLOW_WDL: ${WORKFLOW_WDL}"
echo "TEMPLATE_JSON: ${TEMPLATE_JSON}"
echo "SCAN_JSON: ${SCAN_JSON}"

docker run \
    --network host \
    -v $PWD/launch_study.py:/launch_study.py \
    -v $PWD/${WORKFLOW_WDL}:/workflow.wdl \
    -v $PWD/${TEMPLATE_JSON}:/template.json \
    -v $PWD/${SCAN_JSON}:/scan.json \
    -it ${DOCKER} \
    python2 /launch_study.py \
        --advisor_server ${ADVISOR_SERVER} \
        --cromwell_server ${CROMWELL_SERVER} \
        --study_name ${STUDY_NAME} \
        --algorithm ${ALGORITHM} \
        --workflow_wdl /workflow.wdl \
        --template_json /template.json \
        --scan_json /scan.json \
        --womtool_path /usr/local/bin/womtool/womtool-36.jar
