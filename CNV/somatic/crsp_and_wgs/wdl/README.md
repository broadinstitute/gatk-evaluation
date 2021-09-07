# Validation of the somatic CNV pipeline
The pipeline consists of GATK ModelSegments and other data prep tools. 

## 2021 documentation
### Preparing to run this pipeline
Helper WDLs from GATK repo (in the `gatk_wdls` folder) were pulled in from [this folder](https://github.com/broadinstitute/gatk/tree/master/scripts/cnv_wdl/somatic), after [this commit](https://github.com/broadinstitute/gatk/commit/b4b63baf6f1af8dabc106a98530f86a9dea3c9a6).  

JSON of arguments cannot be committed to this repo because it is public; they are committed to the [hydro.gen repo](https://github.com/broadinstitute/hydro.gen) \[add details\].

## 2018 documentation
Probably broken after 2021 changes
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