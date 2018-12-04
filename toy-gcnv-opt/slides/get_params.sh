#!/bin/bash

echo -e "workflow_id\tgcnv_init_ard_rel_unexplained_variance\tgcnv_sample_psi_scale\tgcnv_interval_psi_scale\tgcnv_log_mean_bias_standard_deviation" > $2
while read w;
do
    metadata=$(curl -s -X GET "https://cromwell-v33.dsde-methods.broadinstitute.org/api/workflows/v1/$w/metadata")
    gcnv_init_ard_rel_unexplained_variance=$(echo $metadata | jq -r '.inputs["EvaluateCNVGermlineCohortWorkflow.gcnv_init_ard_rel_unexplained_variance"]')
    gcnv_sample_psi_scale=$(echo $metadata | jq -r '.inputs["EvaluateCNVGermlineCohortWorkflow.gcnv_sample_psi_scale"]')
    gcnv_interval_psi_scale=$(echo $metadata | jq -r '.inputs["EvaluateCNVGermlineCohortWorkflow.gcnv_interval_psi_scale"]')
    gcnv_log_mean_bias_standard_deviation=$(echo $metadata | jq -r '.inputs["EvaluateCNVGermlineCohortWorkflow.gcnv_log_mean_bias_standard_deviation"]')
    echo -e "$w\t$gcnv_init_ard_rel_unexplained_variance\t$gcnv_sample_psi_scale\t$gcnv_interval_psi_scale\t$gcnv_log_mean_bias_standard_deviation"
done < $1 >> $2
