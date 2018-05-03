# Validation

These scripts are based on the ReCapSeg validation, but run using the GATK4 ModelSegments workflow (GATK CNV 1.5).  
A test for WGS concordance is also performed.

This implements:
- WGS Concordance analysis
- CRSP Clinical Sensitivity
- CRSP Reproducibility (Segmentation)
- CRSP Reproducibility (Caller)
- CRSP Purity Series


## Which samples are TCGA?

```
SM-74P44 	TCGA-HT-A5RC
SM-74P3P 	TCGA-FG-A4MU
SM-74P4H 	TCGA-DU-A5TY
SM-74P5I 	TCGA-55-8615
SM-74P47 	LUAD-TCGA-05-4390
```

These can be found in `gs://broad-dsde-methods/testdata/crsp-private-bams/`

## What is the list of validations?

Ideally:

Exome
-[x] Reproducibility RMSE plot
-[x] Reproducibility RMSE for 3*ploidy (Passing: 0.15)
-[x] Reproducibility call concordance (Passing: 0.97)

-[x] Purity Series (Precision and Sensitivity) Amplificiations
-[x] Purity Series (Precision and Sensitivity) Deletions

-[x] Analytical Sensitivity with TCGA samples (Sensitivity only) Amplificiations
-[x] Analytical Sensitivity with TCGA samples (Sensitivity only) Deletions
-[x] Analytical Sensitivity with TCGA samples (Sensitivity only) Small Amplificiations
-[x] Analytical Sensitivity with TCGA samples (Sensitivity only) Small Deletions

WGS
-[x] Base-pair concordance with PCAWG consensus segments



## How was the centromere file generated?

Step 1 was to go to UCSC genome browser and download the appropriate tracks for the appropriate genome build (hg19 being used here).
Table Browser
Gap: Centromeres and the like


```bash
# Generates the centromere file.
# This script is not as reusable as you might think
INPUT=Gap_hg19.tsv
GATK_JAR=~/IdeaProjects/gatk/build/libs/gatk.jar
FINAL_SEG_FILE=final_centromere_hg19.tsv

set -e

# Not doing XY
grep -Pv "\tchrUn_|chrX|chrY|chr[0-9]+_gl" $INPUT > tmp1_$INPUT

# chr1 --> 1
sed -r "s/\tchr([0-9]+)\t/\t\1\t/g" tmp1_$INPUT | sed -r "s/\#bin/bin/g" > tmp1a_$INPUT

# Remove the segments that start with zero.  
grep -Pv "\t[0-9XY]{1,2}\t0\t[0-9]+\t" tmp1a_$INPUT > tmp2_$INPUT

# alter the header (This step is now optional)
sed -r "s/\tchrom\t/\tCONTIG\t/g" tmp2_$INPUT | sed -r "s/\tstart\t|\tStart_Position\t|\tchromStart\t/\tSTART\t/g" | \
 sed -r "s/\tend\t|\tEnd_Position\t|\tchromEnd\t/\tEND\t/g" | egrep "END|centromere" > $FINAL_SEG_FILE
```

## Example json file 

See `multi_cnv_validation.wdl` in the wdl directory for a working example.