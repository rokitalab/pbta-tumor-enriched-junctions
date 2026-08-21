# Tumor-enriched splice junction assessment

Module authors: Ryan Corbett (@rjcorb)

This module calculates normalized splice junction counts from PBTA rMATS data, and compares with normal brain reference cohort junction counts to identify tumor-enriched junctions 

## Usage
### script to run analysis

```
bash run_module.sh
```

## Tumor-enriched junction classification strategy

1. Compare PBTA junction CPMs with normal-brain reference matrices from GTEx, postnatal Evo-Devo brain, and normal pediatric brain.

   - A junction is eligible for this comparison only when its mean CPM is below 10 in every reference group.
   - For each PBTA sample and junction, calculate:
     - fold change relative to the highest control-group mean CPM; and
     - the minimum signal-to-noise ratio (SNR) across control groups, where SNR is `(tumor CPM − control mean CPM) / control CPM SD`.

2. Classify a junction as **tumor-enriched** when it is more than fivefold higher than every control group and has an SNR greater than 5 for every control group.

   - Retain a call only when no junction sharing either splice boundary is classified as non-specific in the same sample.
   - Junctions absent from all control matrices are also classified as tumor-enriched when they use at least one unannotated splice site and occur in fewer than 250 PBTA samples (less than 10% of the cohort).

3. Classify tumor-enriched junctions as **oncofetal** using Evo-Devo expression.

   - Compare each prenatal region/week-bin group with every postnatal region/stage group.
   - A junction is oncofetal when at least one prenatal group has a minimum prenatal-to-postnatal fold change greater than 2 and a minimum prenatal-versus-postnatal SNR greater than 2 across all postnatal groups.
   - Record the prenatal group or groups meeting these criteria in the output.

## Folder content
* `run_module.sh` shell script to run analysis
* `01-get-junction-counts.R` extract junction counts from PBTA rMATS and normalize. 
* `02-classify-tumor-enriched-junctions.R` identify tumor-enriched junctions by comparing expression in PBTA versus normal control cohorts.
* `03-classify-oncofetal-junctions.R` classify tumor-enriched junctions as oncofetal using prenatal versus postnatal Evo-Devo expression.
* `04-uniprot-domain-annotation.sh` annotate junctions to uniprot topological domains. To retain junctions resulting in alterations or gains of topological domains, we filter junctions to those that either 1) partially overlap a topological domain interval or 2) are competely within a topological domain interval
* `05-pfam-annotation.R` annotate junctions to Pfam functional domains. To assess domain gain, domain loss, and/or domain alteration due to tumor-enriched splicing, we retain all junction-domain overlap types for downstream analyses.
* `06-domain-expression-filtering.R` merge tumor-enriched splice junction annotations, filter for gene TPM > 10.
* `util/rmats-processing-functions.R` scripts containing functions to process raw rMATS data, normal junction and target counts, calculate mean CPMs by subgroups, and generate matrices. 

## Input files
* `pbta-rna-high-intron-samples.tsv` samples to be removed from analyses due to high intronic read fraction

## Directory structure
```
.
├── 01-get-junction-counts.R
├── 02-classify-tumor-enriched-junctions.R
├── 03-classify-oncofetal-junctions.R
├── 04-uniprot-domain-annotation.sh
├── 05-pfam-annotation.R
├── 06-domain-expression-filtering.R
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
│   ├── tumor-enriched-oncofetal-splice-junctions.tsv.gz
│   ├── tumor-enriched-oncofetal-splice-junctions-cds.bed
│   ├── tumor-enriched-oncofetal-splice-junctions.unipLocCytopl.bed
│   ├── tumor-enriched-oncofetal-splice-junctions.unipLocExtra.bed
│   └── tumor-enriched-oncofetal-splice-junctions.unipLocTransMemb.bed
├── run_module.sh
└── util
    └── rmats-processing-functions.R
```
