#!/bin/sh

# Run tumor-enriched junction (TEJ) analysis modules, in order

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Get base directory of project
cd ..
BASEDIR="$(pwd)"
cd -

analyses_dir="$BASEDIR/analyses"

# Create cohort histologies
echo "Create cohort histologies..."
cd ${analyses_dir}/00-create-cohort-histologies
bash run_module.sh

# Process control rMATS splice event matrices
echo "Process control rMATS splice event matrices..."
cd ${analyses_dir}/01-ctrl-rmats-processing
bash run_module.sh

# Process PBTA junction counts and CPM matrices
echo "Process PBTA junction counts and CPM matrices..."
cd ${analyses_dir}/02-pbta-junction-processing
bash run_module.sh

# Classify tumor-enriched junctions
echo "Classify tumor-enriched junctions..."
cd ${analyses_dir}/03-classify-tejs
bash run_module.sh

# Assess TEJ differential splicing
echo "Assess TEJ differential splicing..."
cd ${analyses_dir}/04-tej-differential-splicing
bash run_module.sh

# Summarize tumor-enriched junctions
echo "Summarize tumor-enriched junctions..."
cd ${analyses_dir}/05-summarize-TEJs
bash run_module.sh
