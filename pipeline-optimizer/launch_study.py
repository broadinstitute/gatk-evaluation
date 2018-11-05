#!/usr/bin/env python2

import argparse
import advisor_client.client as AdvisorClient
from cromwell_tools.cromwell_api import CromwellAPI
from cromwell_tools.cromwell_auth import CromwellAuth

parser = argparse.ArgumentParser()
parser.add_argument("--advisor_server")
parser.add_argument("--cromwell_server")
parser.add_argument("--workflow_wdl")
parser.add_argument("--template_json")
parser.add_argument("--hyperparameters_json")
parser.add_argument("--womtool_path")
args = parser.parse_args()


def main():
    print "Validating workflow..."
    CromwellAPI.validate_workflow(args.workflow_wdl, args.womtool_path)
    
    client = AdvisorClient.AdvisorClient(endpoint=args.advisor_server)
    
    cromwell_auth = CromwellAuth(url=args.cromwell_server, header={"Authorization": "bearer fake_token"}, auth=None)
    
    print CromwellAPI.health(cromwell_auth)
    
    print CromwellAPI.submit(cromwell_auth, args.workflow_wdl, args.template_json)


if __name__ == "__main__":
    main()
