#!/usr/bin/env python2

import argparse
import json
import os
import tempfile
import requests
import time
import advisor_client.client as AdvisorClient
from cromwell_tools.cromwell_api import CromwellAPI
from cromwell_tools.cromwell_auth import CromwellAuth

parser = argparse.ArgumentParser()
parser.add_argument('--advisor_server')
parser.add_argument('--cromwell_server')
parser.add_argument('--study_name')
parser.add_argument('--algorithm')
parser.add_argument('--workflow_wdl')
parser.add_argument('--template_json')
parser.add_argument('--scan_json')
parser.add_argument('--womtool_path')
args = parser.parse_args()

def calculate_metric(template_values, scan_values):
    merged_values = {k: v for (k, v) in (template_values.items() + scan_values.items())}
    
    merged_json_file, merged_json_path = tempfile.mkstemp()
    with open(merged_json_path, 'w') as f:
        json.dump(merged_values, f)

    cromwell_auth = CromwellAuth(url=args.cromwell_server, header={'Authorization': 'bearer fake_token'}, auth=None)
    with open(args.workflow_wdl, 'r') as w, open(merged_json_path, 'r') as j:
        submit = CromwellAPI.submit(cromwell_auth, w, j)
        
    workflow_id = submit.json()['id']
    print 'Submitted workflow:', workflow_id
    
    time.sleep(5)
    print 'Waiting for workflow to complete...'
    CromwellAPI.wait([workflow_id], cromwell_auth, timeout_minutes=1440, poll_interval_seconds=120)
    
    metadata = requests.post(url=cromwell_auth.url + CromwellAPI._metadata_endpoint.format(uuid=workflow_id),
                             auth=cromwell_auth.auth,
                             headers=cromwell_auth.header)

    workflow_name = metadata.json()['workflowName']
    objective_value = metadata.json()['outputs']['{}.objective_value'.format(workflow_name)]
    return objective_value 


def main():
    print 'Validating workflow...'
    CromwellAPI.validate_workflow(args.workflow_wdl, args.womtool_path)

    client = AdvisorClient.AdvisorClient(endpoint=args.advisor_server)

    with open(args.scan_json) as f:
        study_configuration = json.load(f)
    max_num_trials = study_configuration['maxTrials']

    with open(args.template_json) as f:
        template_values = json.load(f)

    study = client.get_or_create_study(args.study_name, study_configuration, algorithm=args.algorithm)
    print study

    for i in range(max_num_trials):
        trial = client.get_suggestions(study.name)[0]
        scan_values = json.loads(trial.parameter_values)
        metric = calculate_metric(template_values, scan_values)
        trial = client.complete_trial_with_one_metric(trial, metric)
        print trial

    best_trial = client.get_best_trial(study.name)
    print 'Best trial: {}'.format(best_trial)

if __name__ == '__main__':
    main()
