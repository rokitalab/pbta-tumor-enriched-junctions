# Generate normal brain splice event expression matrices

Module authors: Ryan Corbett (@rjcorb)

This module generates splice event junction counts per million (CPM) matrices in normal brain control cohorts to be used as references in PBTA tumor-enriched splicing analyses. 

## Usage
### script to run analysis

```
bash run_module.sh
```

## Folder content
* `run_module.sh` shell script to run analysis
* `01-create-gtex-matrices.R` generate splice event matrices from GTEx <40 samples
* `02-create-evodevo-matrices.R` generate splice event matrices from Evo-Devo brain samples (excluding middle adult, elderly participants); calculate separate junction count stats for 1) narrow evo-devo subgroups (region + age) and 2) broad subgroups (region + broad stage, fetal or postnatal)
* `03-create-pedbrain-matrices.R` generated splice event matrices from normal pediatric brain samples
* `04-create-brain-celltype-matrices.R` generated splice event matrices from normal pediatric brain cell types
* `util/rmats-processing-functions.R` scripts containing functions to process raw rMATS data, normal junction and target counts, calculate mean CPMs by subgroups, and generate matrices. 

## Input files
* `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` file to obtain participant age ranges
* `gtex-v10-removed-samples.tsv` low-quality samples removed from GTEx v10

## Directory structure
```
.
├── 01-create-gtex-matrices.R
├── 02-create-evodevo-matrices.R
├── 03-create-pedbrain-matrices.R
├── 04-create-brain-celltype-matrices.R
├── README.md
├── input
│   ├── GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
│   └── gtex-v10-removed-samples.tsv
├── results
│   ├── evodevo-merged-broadgroup-psi-mat.qs2
│   ├── evodevo-merged-broadgroup-norm-junction-ct-mat.qs2
│   ├── evodevo-merged-broadgroup-norm-junction-sd-mat.qs2
│   ├── evodevo-merged-subgroup-norm-junction-ct-mat.qs2
│   ├── evodevo-merged-subgroup-norm-junction-sd-mat.qs2
│   ├── evodevo-merged-subgroup-psi-mat.qs2
│   ├── gtex-merged-norm-junction-ct-mat.qs2
│   ├── gtex-merged-norm-junction-sd-mat.qs2
│   ├── gtex-merged-psi-mat.qs2
│   ├── normal-brain-celltype-merged-norm-junction-ct-mat.qs2
│   ├── normal-brain-celltype-merged-norm-junction-sd-mat.qs2
│   ├── normal-brain-celltype-merged-psi-mat.qs2
│   ├── normal-pedbrain-merged-norm-junction-ct-mat.qs2
│   ├── normal-pedbrain-merged-norm-junction-sd-mat.qs2
│   └── normal-pedbrain-merged-psi-mat.qs2
├── run_module.sh
└── util
    └── rmats-processing-functions.R
```