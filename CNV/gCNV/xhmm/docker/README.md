## Docker instructions

### Make sure to have following resources in Dockerfile directory, which we ADD to the docker image:
- GenomeAnalysisTK-3.8-0.jar that can be downloaded here: https://software.broadinstitute.org/gatk/download/archive
- PLINK/Seq(plinkseq-0.10) binary executable file that can be downloaded using http://psychgen.u.hpc.mssm.edu/plinkseq_downloads/plinkseq-x86_64-latest.zip
- XHMM sourcecode that can be found here: https://bitbucket.org/statgen/xhmm/get/master.zip

Then unzip the downloaded resources and make sure that the unzipped directories' names are the same as the ones specified
in Dockerfile. Then run: ```docker build -t xhmm .``` to build docker image
