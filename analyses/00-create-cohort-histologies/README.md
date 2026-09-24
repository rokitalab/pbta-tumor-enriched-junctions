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
* `01-create-histologies.R` create cohort histologies file that includes data from OPC v15 histologies, somalier genetic ancestry prediction, and updated opc survival data from Mar 2026. 
* `02-circos-plot.Rmd` generate cohort circos plot split by plot group 

## Input files
* `plot-mapping.tsv` curated cancer groups for plotting based on `broad_histology` and `cancer_group` field values. 
* `somalier-ancestry-prediction-superpopulation.tsv` somalier genetic ancestry predicition results with BS_IDs appended. From [germline-preprocessing repo](https://github.com/rokitalab/germline-preprocessing/blob/f1fdeced5ad59a29bc0dc88fcba4f4366ab7fa8e/analyses/collapse-tumor-histologies/results/somalier-ancestry-prediction-superpopulation.tsv)
* `openpedcan_histologies_0311.csv` opc histologies warehouse pull dated Mar 11 2026
* `histologies-rare-cns.tsv` histologies file with rare subtyping information pulled from [OpenPedCan](https://github.com/rokitalab/OpenPedCan-Project-CNH/blob/dev/analyses/molecular-subtyping-integrate/results/histologies.tsv) 

## Directory structure
```
.
├── 01-create-histologies.R
├── 02-circos-plot.nb.html
├── 02-circos-plot.Rmd
├── input
│   ├── histologies-rare-cns.tsv
│   ├── openpedcan_histologies0311.csv
│   ├── pbta-rna-high-intron-samples.tsv
│   ├── plot-mapping.tsv
│   └── somalier-ancestry-prediction-superpopulation.tsv
├── plots
│   └── cohort-circos-plot.pdf
├── README.md
├── results
│   ├── cohort-histologies.tsv
│   └── cohort-somalier-genetic-ancestry-prediction.tsv
└── run_module.sh
```