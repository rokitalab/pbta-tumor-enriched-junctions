#!/bin/sh

# if junction counts file not present, run first script
if [ -f "results/pbta-merged-norm-junction-cts.qs2" ]; then
    echo "Found junction counts file. Proceeding..."
else
    echo "Junction counts file does not exist. Running 01-get-junction-counts.R..."
    Rscript 01-get-junction-counts.R
fi

# identify tumor-enriched junctions and classify oncofetal junctions
Rscript --vanilla 02-classify-tumor-enriched-junctions.R
Rscript --vanilla 03-classify-oncofetal-junctions.R

# uniprot domain annotation
bash 04-uniprot-domain-annotation.sh

# pfam domain annotation
Rscript --vanilla 05-pfam-annotation.R

# domain and expression filtering
Rscript --vanilla 06-domain-expression-filtering.R
