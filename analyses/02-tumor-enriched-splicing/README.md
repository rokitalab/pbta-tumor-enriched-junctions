# Tumor-enriched splice junction assessment

Module authors: Ryan Corbett (@rjcorb)

This module calculates normalized splice junction counts from PBTA rMATS data, and compares with normal brain reference cohort junction counts to identify tumor-enriched junctions 

## Usage
### script to run analysis

```
bash run_module.sh
```

## Folder content
* `run_module.sh` shell script to run analysis
* `01-get-junction-counts.R` extract junction counts from PBTA rMATS and normalize. 
* `02-calculate-tumor-enriched-splicing.R` identify tumor-enriched junctions by comparing expression in PBTA versus normal control cohorts. 
* `util/rmats-processing-functions.R` scripts containing functions to process raw rMATS data, normal junction and target counts, calculate mean CPMs by subgroups, and generate matrices. 

## Input files
* `pbta-rna-high-intron-samples.tsv` samples to be removed from analyses due to high intronic read fraction

## Directory structure
```
.
├── 01-get-junction-counts.R
├── 02-calculate-tumor-enriched-splicing.R
├── README.md
├── input
│   └── pbta-rna-high-intron-samples.tsv
├── results
│   ├── junction-annot.tsv.gz
│   ├── pbta-merged-norm-junction-cts.qs2
│   ├── tumor-enriched-oncofetal-splice-junctions.bed
│   └── tumor-enriched-oncofetal-splice-junctions.tsv.gz
├── run_module.sh
└── util
    └── rmats-processing-functions.R
```