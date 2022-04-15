#!/usr/bin/env bash

set -e

REPO=us.gcr.io/broad-dsde-methods/gatk-evaluation
PROJECT=cnv_somatic_1p5
EVALUATION_NAME=crsp_and_wgs
VERSION=1.2   #1.1
FINAL_GCR_TAG=${REPO}/${PROJECT}:${EVALUATION_NAME}_${VERSION}
echo "Building ${FINAL_GCR_TAG}..."
#gcloud docker -- build -t ${FINAL_GCR_TAG} -f Dockerfile .
docker build -t ${FINAL_GCR_TAG} .
echo "Pushing to GCR..."
#gcloud docker -- push ${FINAL_GCR_TAG}
docker push ${FINAL_GCR_TAG}
