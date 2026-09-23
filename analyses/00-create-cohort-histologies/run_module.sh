#!/bin/sh

# create cohort histologies file
Rscript --vanilla 01-create-histologies.R

# generate cohort circos plot
R -e "rmarkdown::render('02-circos-plot.Rmd')"