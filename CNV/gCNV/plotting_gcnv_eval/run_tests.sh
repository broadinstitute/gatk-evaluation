#!/bin/bash

RED='\033[0;31m'
LIGHT_PURPLE='\033[1;35m'
NC='\033[0m' # No Color

PYTHON=$(which python)
TEMP_DIR=~/temp_dir_asmirnov/

#Replace this with this directory path
CODE_PATH=/Users/asmirnov/Desktop/evaluations/gCNV_theano/gatk-evaluation/CNV/gCNV/plotting_gcnv_eval
PLT_GCNV_EVAL_PATH=${CODE_PATH}/plt_gcnv_eval/

export PYTHONPATH="${PYTHONPATH}:${PLT_GCNV_EVAL_PATH}"

mkdir -p $TEMP_DIR

$PYTHON ${CODE_PATH}/tests/interval_collection__test.py ${TEMP_DIR} && 
$PYTHON ${CODE_PATH}/tests/interval__test.py ${TEMP_DIR}

if [[ "$?" -ne "0" ]]; then
    printf "${RED}Testing FAILED${NC} \n"
else
    printf "${LIGHT_PURPLE}Testing SUCCEDED${NC} \n"
fi

