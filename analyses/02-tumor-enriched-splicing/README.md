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
* `03-uniprot-domain-annotation.sh` annotate junction to uniprot topological domains. 
* `04-pfam-annotation.R` annotate junctions to Pfam functional domains. 
* `05-domain-expression-filtering.R` merge tumor-enriched splice junction annotations, filter for gene TPM > 10. 
* `util/rmats-processing-functions.R` scripts containing functions to process raw rMATS data, normal junction and target counts, calculate mean CPMs by subgroups, and generate matrices. 

## Input files
* `pbta-rna-high-intron-samples.tsv` samples to be removed from analyses due to high intronic read fraction

## Directory structure
```
.
├── 01-get-junction-counts.R
├── 02-calculate-tumor-enriched-splicing.R
├── 03-uniprot-domain-annotation.sh
├── 04-pfam-annotation.R
├── 05-domain-expression-filtering.Rs
├── README.md
├── input
│   ├── add-tpm-values.R
│   └── pbta-rna-high-intron-samples.tsv
├── results
│   ├── junction-annot.tsv.gz
│   ├── pbta-merged-norm-junction-cts.qs2
│   ├── tumor-enriched-oncofetal-splice-junctions-pfam-annotated.tsv.gz
│   ├── tumor-enriched-oncofetal-splice-junctions-domain-expr-annotated.tsv.gz
│   ├── tumor-enriched-oncofetal-splice-junctions.bed
│   ├── tumor-enriched-oncofetal-splice-junctions.tsv.gz
│   ├── tumor-enriched-oncofetal-splice-events-domain-anno.uniq.tsv
│   ├── tumor-enriched-oncofetal-splice-junctions-cds.bed
│   ├── tumor-enriched-oncofetal-splice-junctions.unipLocCytopl.bed
│   ├── tumor-enriched-oncofetal-splice-junctions.unipLocExtra.bed
│   └── tumor-enriched-oncofetal-splice-junctions.unipLocTransMemb.bed
├── run_module.sh
└── util
    └── rmats-processing-functions.R
```