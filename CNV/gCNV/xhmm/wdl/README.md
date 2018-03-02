## WDL notes

### Optional steps
Two steps are optional in XHMM pipeline: 
- GetExtremeGCContentTargets which compiles a list of high/low GC content targets and will only run if gc_low_high_filter_params_final
parameter is specified and is an array that contains exactly two values. 
Note that WDL does not validate that array values are between 0.0 and 1.0
- GetLowComplexityTargets which compiles a list of low complexity targets and will be run if seq_db is specified which
is a resource that PLINK/SEQ uses and it can be downloaded here:
http://psychgen.u.hpc.mssm.edu/plinkseq_resources/hg19/seqdb.hg19.gz

### Parameters
The example json is included in the directory.
Here is a list of filtering parameters along with the default values that XHMM tutorial uses:
- "XHMM.interval_complexity_filter_threshold": 0.25
- "XHMM.gc_low_high_filter_params": [0.1, 0.9]
- "XHMM.min_target_size_filter": 10
- "XHMM.max_target_size_filter": 10000
- "XHMM.min_mean_target_rd": 10
- "XHMM.max_mean_target_rd": 500,
- "XHMM.min_mean_sample_rd": 25,
- "XHMM.max_mean_sample_rd": 200,
- "XHMM.max_sd_sample_rd": 150
