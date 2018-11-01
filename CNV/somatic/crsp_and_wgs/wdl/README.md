
### What is deploy.sh.template?

You can make a copy of `deploy.sh.template` and modify the parameters to create the dependency zip file needed to run the evaluation in cromwell.

### How do I run this WDL, json, etc?

That is out of scope for this document, since that will depend on backend.  Note that the example json file should work using PAPI and having access to the `broad-dsde-methods` buckets in GCP.

Also, use `deploy.sh.template`, since that can get you to the point of being able to use `cromshell` (https://github.com/broadinstitute/cromshell)

For example, submitting to a cromwell server: 
```bash
cd /path/to/my/workdir/in/deploy.sh/
cromshell submit multi_cnv_validation.wdl multi_cnv_validation.json [path_to_options_json_file] eval_cnv_wgs_validation.zip
```

If you need to submit to a different server:
```bash
export CROMWELL_URL="https://my_server:cromwell_port"
```

### Is there json for creating the PoNs?

Yes, see `pon_input_*.json` files for PoN creation for the ICE and WGS portions of this evaluation.