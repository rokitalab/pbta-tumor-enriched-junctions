# Create cohort histologies file

Module authors: Ryan Corbett (@rjcorb)

This module creates a cohort histologies file for the tumor-enriched splicing project, including all sample- and patient-level metadata to conduct downstream analyses. 

## Usage
### script to run analysis

```
bash run_module.sh
```

## Folder content
* `run_module.sh` shell script to run analysis
* `01-create-histologies.R` create cohort histologies file that includes data from OPC v15 histologies, somalier genetic ancestry prediction, and updated opc survival data from Sep 2025. 

## Input files
* `plot-mapping.tsv` curated cancer groups for plotting based on `broad_histology` and `cancer_group` field values. 
* `somalier-ancestry-prediction-superpopulation.tsv` somalier genetic ancestry predicition results with BS_IDs appended. From [germline-preprocessing repo](https://github.com/rokitalab/germline-preprocessing/blob/f1fdeced5ad59a29bc0dc88fcba4f4366ab7fa8e/analyses/collapse-tumor-histologies/results/somalier-ancestry-prediction-superpopulation.tsv)
* `openpedcan_histologies_0311.csv` opc histologies warehouse pull dated Mar 11 2026

## Directory structure
```
.
├── 01-create-histologies.R
├── README.md
├── input
│   ├── openpedcan_histologies_20250924.csv
│   ├── plot-mapping.tsv
│   └── somalier-ancestry-prediction-superpopulation.tsv
├── results
│   ├── cohort-histologies.tsv
│   └── cohort-somalier-genetic-ancestry-prediction.tsv
└── run_module.sh
```