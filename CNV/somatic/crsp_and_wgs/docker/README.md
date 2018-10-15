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



## How was the centromere file generated (hg19)?

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

## How was the centromere file generated (hg38)?

Step 1 was to go to UCSC genome browser and download the appropriate tracks for the appropriate genome build (hg38).
Table Browser
Group: Mapping and Sequencing | Track: Centromere | Output format of BED


```bash
# Generates the centromere file.
# This script is not as reusable as you might think
INPUT=centromeres_hg38.bed
FINAL_SEG_FILE=final_centromere_hg38.seg
HG38_DICT=/home/lichtens/broad_oncotator_configs/ref/hg38/Homo_sapiens_assembly38.dict
set -e

# Not doing XY
grep -Pv "\tchrUn_|chrX|chrY|chr[0-9]+_gl" $INPUT > tmp1_$INPUT

# Remove the segments that start with zero.  
grep -Pv "\t[0-9XY]{1,2}\t0\t[0-9]+\t" tmp1_$INPUT > tmp2_$INPUT

# sort bed 
sort -k 1.4,1n -k 2,2n -k 3,3n tmp2_$INPUT > tmp2a_$INPUT

# Add header
cat $HG38_DICT >tmp3_$INPUT
echo -e "CONTIG\tSTART\tEND\ttype" >>tmp3_$INPUT
cat tmp2a_$INPUT >>tmp3_$INPUT

# alter the names:  E.g. GJ211954.1  --> "centromere"
sed -r "s/\tGJ[0-9]+\.[0-9]+/\tcentromere/g" tmp3_$INPUT  > $FINAL_SEG_FILE


```


## Example json file 

See `multi_cnv_validation.json` in the wdl directory for a working example.


## Converting GATK seg file (b37) into interval list and lifting over to hg38

```bash

INPUT_FILE=CNV.hg19.bypos.v1.CR1_event_added.mod.seg
OUTPUT_FILE=CNV.hg38liftover.bypos.v1.CR1_event_added.mod.seg

# ftp://gsapubftp-anonymous@ftp.broadinstitute.org/Liftover_Chain_Files/
B37TOHG19CHAIN=b37tohg19.chain

CHAIN_FILE=hg19ToHg38.over.chain
GATK_JAR=gatk.jar
TARGET_REF1=ucsc.hg19.dict
TARGET_REF2=Homo_sapiens_assembly38.dict

# WARNING -- this script will drop some of the files in the original input file.


# First we save the CNV file as an interval_list
egrep "^\@" ${INPUT_FILE} > tmp_old_dict
cut -f1-4 ${INPUT_FILE} | egrep -v "^\@|CONTIG" >tmp1_$INPUT_FILE
awk '{ print $2 "\t" $3 "\t" $4 "\t+\t" $1}' tmp1_$INPUT_FILE > tmp1a_$INPUT_FILE
egrep CONTIG ${INPUT_FILE} | awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' > new_header
cat tmp_old_dict > ${INPUT_FILE}.interval_list
cat tmp1a_$INPUT_FILE >> ${INPUT_FILE}.interval_list


# Liftover the interval list (twice -- once b37->hg19, then hg19->hg38)
java -jar ${GATK_JAR} LiftOverIntervalList \
-I ${INPUT_FILE}.interval_list \
-O tmp_${INPUT_FILE}.hg19.interval_list \
-SD ${TARGET_REF1} \
-CHAIN ${B37TOHG19CHAIN}

java -jar ${GATK_JAR} LiftOverIntervalList \
-I tmp_${INPUT_FILE}.hg19.interval_list \
-O tmp_${OUTPUT_FILE}.hg38.interval_list \
-SD ${TARGET_REF2} \
-CHAIN ${CHAIN_FILE}

egrep "^@" tmp_${OUTPUT_FILE}.hg38.interval_list > ${OUTPUT_FILE}
cat new_header >> ${OUTPUT_FILE}
egrep -v "^@" tmp_${OUTPUT_FILE}.hg38.interval_list  >> ${OUTPUT_FILE}


echo "Entries in the input file:"
egrep -v "^@" ${INPUT_FILE} | wc -l

echo "Entries in the output file:"
egrep -v "^@" ${OUTPUT_FILE} | wc -l
```