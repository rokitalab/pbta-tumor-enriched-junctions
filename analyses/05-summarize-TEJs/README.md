# Tumor-enriched junction summary analyses

Module authors: Ryan Corbett (@rjcorb)

This analysis module performs the following:

* Generates tumor-enriched splice junction (TEJ) summary plots and identifies recurrent TEJs in primary tumors by histology.
* Annotates TEJs to transcript isoforms and exons. 
* Generates CPM expression matrices for recurrent primary TEJs in all PBTA and control cohort samples. 
* Performs broad functional characterization of TEJs

## Usage
### script to run analysis

```
bash run_module.sh
```

## Folder content
* `run_module.sh` shell script to run analysis
* `01-summary.Rmd`: Identify primary & recurrent TEJs in PBTA, and generate summary plots 
* `02-splice-site-annotation.Rmd`: Annotate all TEJs to transcript isoforms, exons, and introns when annotation exists
* `03-create-recurrent-tej-cpm-matrix.R`: create PBTA TEJ CPM matrix
* `04-create-control-tej-cpm-matrix.R`: create PBTA TEJ CPM matrix
* `05-tej-functional-summary.Rmd`: assess number of TEJs associated with protein pfam and extracellular domains. 
* `06-generate-tej-psi-matrix.Rmd`: extract PSIs for all TEJ-associated splice events and save as matrix. 


## Directory structure
```
.
├── 01-summary.Rmd
├── 01-summary.nb.html
├── 02-splice-site-annotation.Rmd
├── 02-splice-site-annotation.nb.html
├── 03-create-recurrent-tej-cpm-matrix.R
├── 04-create-control-tej-cpm-matrix.R
├── 05-tej-functional-summary.Rmd
├── 05-tej-functional-summary.nb.html
├── 06-generate-tej-psi-matrix.Rmd
├── 06-generate-tej-psi-matrix.nb.html
├── README.md
├── plots
│   ├── ATRT-recurrent-ec-tejs-by-gene.pdf
│   ├── CNSemb-recurrent-ec-tejs-by-gene.pdf
│   ├── CPT-recurrent-ec-tejs-by-gene.pdf
│   ├── CRANIO-recurrent-ec-tejs-by-gene.pdf
│   ├── DMG-recurrent-ec-tejs-by-gene.pdf
│   ├── EPN-recurrent-ec-tejs-by-gene.pdf
│   ├── GCT-recurrent-ec-tejs-by-gene.pdf
│   ├── GNT-recurrent-ec-tejs-by-gene.pdf
│   ├── HGG-recurrent-ec-tejs-by-gene.pdf
│   ├── LGG-recurrent-ec-tejs-by-gene.pdf
│   ├── MB-recurrent-ec-tejs-by-gene.pdf
│   ├── MES-recurrent-ec-tejs-by-gene.pdf
│   ├── MNG-recurrent-ec-tejs-by-gene.pdf
│   ├── NFP-recurrent-ec-tejs-by-gene.pdf
│   ├── NNT-recurrent-ec-tejs-by-gene.pdf
│   ├── OLIGO-recurrent-ec-tejs-by-gene.pdf
│   ├── Other-recurrent-ec-tejs-by-gene.pdf
│   ├── SWN-recurrent-ec-tejs-by-gene.pdf
│   ├── atrt-hgg-recurrent-ec-tejs-by-gene.pdf
│   ├── recurrent-tej-n-by-sample-hist.pdf
│   ├── tej-by-specificity-hist-barplot.pdf
│   ├── tej-n-by-consequence.pdf
│   ├── tej-n-by-hist-barplot.pdf
│   ├── tejs-by-event-type-hist-barplot.pdf
│   ├── tejs-by-junction-status-hist-barplot.pdf
│   └── total-tej-n-by-sample-hist.pdf
├── results
│   ├── all-tumor-enriched-oncofetal-splice-junctions-annotated.tsv.gz
│   ├── control-cohort-histologies.tsv
│   ├── recurrent-primary-tej-associated-splice-events.tsv
│   ├── recurrent-primary-tumor-enriched-oncofetal-splice-event-psis.rds
│   ├── recurrent-primary-tumor-enriched-oncofetal-splice-junctions-annotated-domain-updated.tsv.gz
│   ├── recurrent-primary-tumor-enriched-oncofetal-splice-junctions-annotated.tsv.gz
│   ├── recurrent-primary-tumor-enriched-oncofetal-splice-junctions.tsv.gz
│   ├── tumor-enriched-oncofetal-splice-junction-cpm-ctrls.rds
│   ├── tumor-enriched-oncofetal-splice-junction-cpm.rds
│   └── tumor-enriched-oncofetal-splice-junctions-cohort-filtered.tsv.gz
├── run_module.sh
└── util
    ├── annotation_functions.R
    └── rmats-processing-functions.R
```