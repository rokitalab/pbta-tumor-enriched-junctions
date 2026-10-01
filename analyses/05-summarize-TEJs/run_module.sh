#!/bin/sh

# Identify primary and recurrent TEJs in PBTA and generate summary plots
R -e "rmarkdown::render('01-summary.Rmd')"

# Identify primary and recurrent TEJs in PBTA and generate summary plots
R -e "rmarkdown::render('02-splice-site-annotation.Rmd')"

# Create PBTA cpm matrix
Rscript --vanilla 03-create-recurrent-tej-cpm-matrix.R

# Create PBTA cpm matrix
Rscript --vanilla 04-create-control-tej-cpm-matrix.R

# summarize TEJ functional consequences
R -e "rmarkdown::render('05-tej-functional-summary.Rmd')"

# generate TEJ splice event PSI matrix
R -e "rmarkdown::render('06-generate-tej-psi-matrix.Rmd')"