# Tumor-enriched junction summary analyses

Module authors: Ryan Corbett (@rjcorb)

This analysis module performs the following:

* Generates tumor-enriched splice junction (TEJ) summary plots and identifies recurrent TEJs in primary tumors by histology.
* Annotates TEJs to transcript isoforms and exons. 
* Generates CPM expression matrices for recurrent primary TEJs in all PBTA and control cohort samples. 

## Usage
### script to run analysis

```
bash run_module.sh
```

## Folder content
* `run_module.sh` shell script to run analysis
* `01-summary.Rmd`: Identify primary & recurrent TEJs in PBTA, and generate summary plots 
* `02-splice-site-annotation.Rmd`: Annotate all TEJs to transcript isoforms, exons, and introns when annotation exists

## Directory structure
```
.
├── 01-summary.Rmd
├── 01-summary.nb.html
├── 02-splice-site-annotation.Rmd
├── 02-splice-site-annotation.nb.html
├── README.md
├── plots
│   ├── recurrent-tej-n-by-sample-hist.pdf
│   ├── tej-by-specificity-hist-barplot.pdf
│   ├── tej-n-by-hist-barplot.pdf
│   ├── tejs-by-event-type-hist-barplot.pdf
│   ├── tejs-by-novelss-hist-barplot.pdf
│   └── total-tej-n-by-sample-hist.pdf
├── results
│   ├── all-tumor-enriched-oncofetal-splice-junctions-annotated.tsv.gz
│   ├── recurrent-primary-tumor-enriched-oncofetal-splice-junctions-annotated.tsv.gz
│   ├── recurrent-primary-tumor-enriched-oncofetal-splice-junctions.tsv.gz
│   └── tumor-enriched-oncofetal-splice-junctions-cohort-filtered.tsv.gz
├── run_module.sh
└── util
    └── annotation_functions.R
```