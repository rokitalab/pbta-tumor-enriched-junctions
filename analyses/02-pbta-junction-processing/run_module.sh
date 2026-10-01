#!/bin/sh

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root_dir=$(CDPATH= cd -- "$script_dir/../.." && pwd)
results_dir="$root_dir/analyses/02-pbta-junction-processing/results"

cd "$script_dir" || exit 1

# If normalized junction counts are not available, generate them first.
if [ -f "$results_dir/pbta-merged-norm-junction-cts.qs2" ]; then
    echo "Found normalized PBTA junction counts. Proceeding..."
else
    echo "Normalized PBTA junction counts do not exist. Running 01-get-junction-counts.R..."
    Rscript --vanilla 01-get-junction-counts.R
fi

# Create the all-PBTA junction CPM matrix.
Rscript --vanilla 02-create-all-pbta-junction-cpm-matrix.R

# Correct all-PBTA junction CPMs for RNA library-preparation batch effects.
Rscript --vanilla 03-batch-correct-all-pbta-junction-cpm.R
