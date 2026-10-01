#!/bin/sh

# identify tumor-enriched junctions and classify oncofetal junctions
Rscript --vanilla 01-classify-tumor-enriched-junctions.R
Rscript --vanilla 02-classify-oncofetal-junctions.R

# uniprot domain annotation
bash 03-uniprot-domain-annotation.sh

# pfam domain annotation
Rscript --vanilla 04-pfam-annotation.R

# domain and expression filtering
Rscript --vanilla 05-domain-expression-filtering.R
