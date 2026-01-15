# Generate normal brain splice event expression matrices

Module authors: Ryan Corbett (@rjcorb)

This module generates splice event junction and target counts per million (CPM) matrices in normal brain control cohorts to be used as references in PBTA tumor-enriched splicing analyses. 

## Usage
### script to run analysis

```
bash run_module.sh  --dpsi <dpsi> --n_na <n_na>
```

## Folder content
* `run_module.sh` shell script to run analysis
* `01-create-gtex-matrices.R` generate splice event matrices from GTEx <40 samples
* `02-create-evodevo-matrices.R` generate splice event matrices from Evo-Devo brain samples (excluding middle adult, elderly participants)
* `03-create-pedbrain-matrices.R` generated splice event matrices from normal pediatric brain samples
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
├── README.md
├── input
│   ├── GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
│   └── gtex-v10-removed-samples.tsv
├── results
│   ├── evodevo-a3ss-norm-junction-ct-mat.qs2
│   ├── evodevo-a5ss-norm-junction-ct-mat.qs2
│   ├── evodevo-merged-norm-junction-ct-mat.qs2
│   ├── evodevo-ri-norm-junction-ct-mat.qs2
│   ├── evodevo-ri-norm-target-ct-mat.qs2
│   ├── evodevo-se-norm-junction-ct-mat.qs2
│   ├── evodevo-se-norm-target-ct-mat.qs2
│   ├── gtex-a3ss-norm-junction-ct-mat.qs2
│   ├── gtex-a5ss-norm-junction-ct-mat.qs2
│   ├── gtex-merged-norm-junction-ct-mat.qs2
│   ├── gtex-ri-norm-junction-ct-mat.qs2
│   ├── gtex-ri-norm-target-ct-mat.qs2
│   ├── gtex-se-norm-junction-ct-mat.qs2
│   ├── gtex-se-norm-target-ct-mat.qs2
│   ├── normal-pedbrain-a3ss-norm-junction-ct-mat.qs2
│   ├── normal-pedbrain-a5ss-norm-junction-ct-mat.qs2
│   ├── normal-pedbrain-merged-norm-junction-ct-mat.qs2
│   ├── normal-pedbrain-ri-norm-junction-ct-mat.qs2
│   ├── normal-pedbrain-ri-norm-target-ct-mat.qs2
│   ├── normal-pedbrain-se-norm-junction-ct-mat.qs2
│   └── normal-pedbrain-se-norm-target-ct-mat.qs2
├── run_module.sh
└── util
    └── rmats-processing-functions.R
```