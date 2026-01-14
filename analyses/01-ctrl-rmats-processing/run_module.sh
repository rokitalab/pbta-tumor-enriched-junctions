#!/bin/sh

# Generate GTEx matrices
Rscript --vanilla 01-create-gtex-matrices.R

# Generate Evo-devo matrices
Rscript --vanilla 02-create-evodevo-matrices.R

# Generate Evo-devo matrices
Rscript --vanilla 03-create-pedbrain-matrices.R
