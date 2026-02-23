#!/bin/sh

# get TESJ-associated splice events
Rscript --vanilla 01-get-junction-splice-events.R

# calculate TESJ-associated splice event differential splicing
Rscript --vanilla 02-get-junction-diff-splicing.R