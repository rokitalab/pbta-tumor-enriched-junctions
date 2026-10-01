# PBTA junction processing

Module authors: Ryan Corbett (@rjcorb)

This module extracts and normalizes PBTA rMATS junction counts, creates an all-PBTA splice-junction CPM matrix, and corrects junction CPMs for RNA library-preparation batch effects.

## Usage
### script to run analysis

```
bash run_module.sh
```

## Folder content
* `run_module.sh` shell script to run analysis
* `01-get-junction-counts.R` extract PBTA junction counts from rMATS files and normalize them by sample input-read count.
* `02-create-all-pbta-junction-cpm-matrix.R` create an all-PBTA junction CPM matrix from normalized PBTA junction counts, excluding junctions expressed above the CPM cutoff in normal control cohorts.
* `03-batch-correct-all-pbta-junction-cpm.R` correct all-PBTA junction CPMs for RNA library-preparation batch effects using ComBat.
* `util/rmats-processing-functions.R` functions for processing rMATS data and normalizing junction counts.
* `util/junction-cpm-matrix.functions.R` functions for constructing junction CPM matrices.

## Input files
* `data/pbta-rmats_merged_raw_{SE,RI,A3SS,A5SS}.qs2` PBTA rMATS splice-event data.
* `data/pbta_input_read_counts.tsv` PBTA input-read counts used to normalize junction counts.
* `analyses/00-create-cohort-histologies/input/pbta-rna-high-intron-samples.tsv` PBTA samples excluded because of high intronic read fraction.
* `analyses/00-create-cohort-histologies/results/cohort-histologies.tsv` PBTA sample identifiers used to select samples for the CPM matrix.
* `analyses/01-ctrl-rmats-processing/results/*-merged-norm-junction-ct-mat.qs2` normal control-cohort junction CPM matrices used to exclude expressed junctions.

## Directory structure
```
.
├── 01-get-junction-counts.R
├── 02-create-all-pbta-junction-cpm-matrix.R
├── 03-batch-correct-all-pbta-junction-cpm.R
├── README.md
├── results
│   ├── all-pbta-splice-junction-cpm.qs2
│   ├── all-pbta-splice-junction-log2-cpm-combat-corrected.qs2
│   ├── junction-annot.tsv.gz
│   ├── pbta-merged-norm-batch-corrected-junction-cts.qs2
│   └── pbta-merged-norm-junction-cts.qs2
├── run_module.sh
└── util
    ├── junction-cpm-matrix.functions.R
    └── rmats-processing-functions.R
```
