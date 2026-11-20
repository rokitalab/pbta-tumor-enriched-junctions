#!/bin/sh

# create cohort histologies file
Rscript --vanilla 01-create-histologies.R

# create control cohort histologies file
Rscript --vanilla 03-create-control-cohort-histologies.R

# generate cohort circos plot
R -e "rmarkdown::render('02-circos-plot.Rmd')"