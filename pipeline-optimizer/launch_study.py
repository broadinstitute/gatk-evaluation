#!/usr/bin/env python2

import argparse
import json
import os
import tempfile
import requests
import time
import sys
import logging
import advisor_client.client as AdvisorClient
from cromwell_tools.cromwell_api import CromwellAPI
from cromwell_tools.cromwell_auth import CromwellAuth
from cromwell_tools.cromwell_api import WorkflowFailedException
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter

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

bad_value = None

# https://stackoverflow.com/questions/49121365/implementing-retry-for-requests-in-python
def retry_session(retries, session=None, backoff_factor=0.3, status_forcelist=(500, 502, 503, 504)):
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session

# https://stackoverflow.com/questions/28330317/print-timestamp-for-logging-in-python
def setup_custom_logger(name):
    formatter = logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')
    handler = logging.FileHandler('log.txt', mode='w')
    handler.setFormatter(formatter)
    screen_handler = logging.StreamHandler(stream=sys.stdout)
    screen_handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    logger.addHandler(screen_handler)
    return logger

logger = setup_custom_logger('logger')

def calculate_metric(template_values, scan_values):
    merged_values = {k: v for (k, v) in (template_values.items() + scan_values.items())}
    
    merged_json_file, merged_json_path = tempfile.mkstemp()
    with open(merged_json_path, 'w') as f:
        json.dump(merged_values, f)

    cromwell_auth = CromwellAuth(url=args.cromwell_server, header={'Authorization': 'bearer fake_token'}, auth=None)
    with open(args.workflow_wdl, 'r') as w, open(merged_json_path, 'r') as j:
        submit = CromwellAPI.submit(cromwell_auth, w, j)
        
    workflow_id = submit.json()['id']
    logger.info('Submitted workflow: ' + workflow_id)
    
    time.sleep(5)
    logger.info('Waiting for workflow to complete...')
    
    # Query workflow status indefinitely until success or failure returned.
    # If success returned, attempt to retrieve objective_value from metadata and return.
    # If failure returned or if exception raised during metadata retreival, return bad_value.
    try:
        while True:
            try:
                CromwellAPI.wait([workflow_id], cromwell_auth, timeout_minutes=600, poll_interval_seconds=20, verbose=False)
                response = CromwellAPI.status(workflow_id, cromwell_auth)
                status = response.json()['status']
                if status == 'Succeeded':
                    logger.info('Workflow succeeded...')
                    break
            except WorkflowFailedException:
                logger.info('Workflow failed, returning bad value...')
                return bad_value
            except Exception as e:
                logger.info(e)
                logger.info('Cromwell exception, retrying wait and status check...')
        logger.info('Getting metadata...')
        session = retry_session(retries=10)
        metadata = session.post(url=cromwell_auth.url + CromwellAPI._metadata_endpoint.format(uuid=workflow_id),
                                auth=cromwell_auth.auth,
                                headers=cromwell_auth.header)
        workflow_name = metadata.json()['workflowName']
        objective_value = metadata.json()['outputs']['{}.objective_value'.format(workflow_name)]
        return objective_value
    except Exception as e:
        logger.info(e)
        logger.info('Cromwell exception during metadata retrieval, returning bad value...')
        return bad_value

def main():
    logger.info('Validating workflow...')
    CromwellAPI.validate_workflow(args.workflow_wdl, args.womtool_path)

    client = AdvisorClient.AdvisorClient(endpoint=args.advisor_server)

    with open(args.scan_json) as f:
        study_configuration = json.load(f)
    max_num_trials = study_configuration['maxTrials']

    with open(args.template_json) as f:
        template_values = json.load(f)

    study = client.get_or_create_study(args.study_name, study_configuration, algorithm=args.algorithm)
    logger.info(study)

    for i in range(max_num_trials):
        try:
            trial = client.get_suggestions(study.name)[0]
            logger.info(trial)
            scan_values = json.loads(trial.parameter_values)
            metric = calculate_metric(template_values, scan_values)
            if metric is bad_value:
                logger.info('Trial returned bad value, skipping...')
                continue
            trial = client.complete_trial_with_one_metric(trial, metric)
            logger.info('Objective value: ' + str(metric))
            logger.info('Trial completed.')
        except Exception as e:
            logger.info(e)
            logger.info('Problem with trial, skipping...')
            continue

    best_trial = client.get_best_trial(study.name)
    logger.info('Best trial: {}'.format(best_trial))

if __name__ == '__main__':
    main()
