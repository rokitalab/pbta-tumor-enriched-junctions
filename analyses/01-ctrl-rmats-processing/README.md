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
* `02-create-evodevo-matrices.R` generates splice-event matrices from Evo-Devo forebrain and hindbrain samples, excluding middle-adult and elderly participants. Postnatal matrices are grouped by region and recorded stage and use mean PSI/CPM values; prenatal matrices are grouped by region and gestational-age bins (4–5, 6–7, 8–12, 13–16, and 17–19 post-conception weeks) and use median PSI/CPM values. The script also writes the derived prenatal-bin metadata.
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
│   ├── evodevo-brain-prenatal-week-binned-metadata.tsv
│   ├── evodevo-merged-postnatal-norm-junction-ct-mat.qs2
│   ├── evodevo-merged-postnatal-norm-junction-sd-mat.qs2
│   ├── evodevo-merged-postnatal-psi-mat.qs2
│   ├── evodevo-merged-prenatal-week-binned-norm-junction-ct-mat.qs2
│   ├── evodevo-merged-prenatal-week-binned-norm-junction-sd-mat.qs2
│   ├── evodevo-merged-prenatal-week-binned-psi-mat.qs2
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
