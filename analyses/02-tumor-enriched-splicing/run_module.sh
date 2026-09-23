#!/bin/sh

# if junction counts file not present, run first script
if [ -f "results/pbta-merged-norm-junction-cts.qs2" ]; then
    echo "Found junction counts file. Proceeding..."
else
    echo "Junction counts file does not exist. Running 01-get-junction-counts.R..."
    Rscript 01-get-junction-counts.R
fi

# identify tumor-enriched junctions
Rscript --vanilla 02-calculate-tumor-enriched-splicing.R

# uniprot domain annotation
bash 03-uniprot-domain-annotation.sh

# pfam domain annotation
Rscript --vanilla 04-pfam-annotation.R

# domain and expression filtering
Rscript --vanilla 05-domain-expression-filtering.R