#!/bin/sh

# if junction counts file not present, run first script
if [ -f "results/pbta-merged-norm-junction-cts.qs2" ]; then
    echo "Found junction counts file. Proceeding..."
else
    echo "Junction counts file does not exist. Running 01-get-junction-counts.R..."
    Rscript 01-get-junction-counts.R
fi

