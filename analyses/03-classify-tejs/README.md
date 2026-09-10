# Classify tumor-enriched splice junctions

Module authors: Ryan Corbett (@rjcorb)

This module classifies PBTA splice junctions as tumor-enriched and oncofetal, then annotates them with protein-domain and gene-expression information.

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
   - Junctions absent from all control matrices are also classified as tumor-enriched when they use at least one unannotated splice site and occur in fewer than 250 PBTA samples (less than 10% of the cohort). Before this classification, these candidates are cross-checked against the control STAR junction file (`SJ.merged.control-cohort.tsv.gz`) and excluded when chromosome and intron coordinates match exactly: the STAR intron start is the upstream exon boundary plus one, and the intron end is the downstream exon boundary minus one.

3. Classify tumor-enriched junctions as **oncofetal** using Evo-Devo expression.

   - Compare each prenatal region/week-bin group with every postnatal region/stage group.
   - A junction is oncofetal when at least one prenatal group has a minimum prenatal-to-postnatal fold change greater than 2 and a minimum prenatal-versus-postnatal SNR greater than 2 across all postnatal groups.
   - Record the prenatal group or groups meeting these criteria in the output.

## Folder content
* `run_module.sh` shell script to run analysis
* `01-classify-tumor-enriched-junctions.R` identify tumor-enriched junctions by comparing PBTA expression with normal control cohorts.
* `02-classify-oncofetal-junctions.R` classify tumor-enriched junctions as oncofetal using prenatal versus postnatal Evo-Devo expression.
* `03-uniprot-domain-annotation.sh` annotate junctions to UniProt topological domains.
* `04-pfam-annotation.R` annotate junctions to Pfam functional domains.
* `05-domain-expression-filtering.R` merge junction annotations and retain junctions in genes with TPM >= 10.

## Input files
* `analyses/02-pbta-junction-processing/results/pbta-merged-norm-batch-corrected-junction-cts.qs2` batch-corrected PBTA junction CPMs.
* `analyses/02-pbta-junction-processing/results/junction-annot.tsv.gz` PBTA splice-junction annotations.

## Directory structure
```
.
├── 01-classify-tumor-enriched-junctions.R
├── 02-classify-oncofetal-junctions.R
├── 03-uniprot-domain-annotation.sh
├── 04-pfam-annotation.R
├── 05-domain-expression-filtering.R
├── README.md
├── results
│   ├── tumor-enriched-junctions-atrt.qs2
│   ├── tumor-enriched-oncofetal-splice-junctions-atrt.bed
│   ├── tumor-enriched-oncofetal-splice-junctions-cds-atrt.bed
│   ├── tumor-enriched-oncofetal-splice-junctions-domain-anno-atrt.uniq.tsv
│   ├── tumor-enriched-oncofetal-splice-junctions-domain-expr-annotated-atrt.tsv.gz
│   ├── tumor-enriched-oncofetal-splice-junctions-pfam-annotated-atrt.tsv.gz
│   └── tumor-enriched-oncofetal-splice-junctions-atrt.tsv.gz
├── run_module.sh
└── util
    └── add-tpm-values.R
```
