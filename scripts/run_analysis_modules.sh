#!/bin/sh

# Run tumor-enriched splicing analysis modules

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Get base directory of project
cd ..
BASEDIR="$(pwd)"
cd -

analyses_dir="$BASEDIR/analyses"

# Run process-ctrls-variants analysis module
# echo "Process control SE events..."
# cd ${analyses_dir}/process-ctrls-variants
# bash run-module.sh
# 
# # Run tumor-enriched SE events detection analysis module
# echo "Run tumor-enriched SE event detection..."
# cd ${analyses_dir}/01-tumor-specific-variants
# bash run_module.sh
# 
# # Run tumor-enriched AltSS events detection analysis module
# echo "Run tumor-enriched AltSS event detection..."
# cd ${analyses_dir}/tumor-specific-alt-splice-sites
# bash run-module.sh
# 
# # Run tumor-enriched RI events detection analysis module
# echo "Run tumor-enriched RI event detection..."
# cd ${analyses_dir}/tumor-specific-retained-introns
# bash run-module.sh
# 
# # Run TESE summary analysis module
# echo "Run tumor-enriched splice event summary..."
# cd ${analyses_dir}/summarize-tumor-enriched-splicing
# bash run-module.sh
# 
# # Run TESE clustering analysis module
# echo "Run tumor-enriched splice event clustering..."
# cd ${analyses_dir}/sample-clustering
# bash run-module.sh
# 
# # Run SF correlation analysis module
# echo "Run SF correlation analyses..."
# cd ${analyses_dir}/splicing-rna-binding-protein-correlations
# bash run_module.sh
# 
# # Run translation analysis module
# echo "Run translation..."
# cd ${analyses_dir}/translation
# bash run_module.sh

# Run atrt analysis module
echo "Run ATRT-specific analyses..."
cd ${analyses_dir}/atrt-enriched-splicing
bash run_module.sh

# Run hgg analysis module
echo "Run HGG-specific analyses..."
cd ${analyses_dir}/hgg-enriched-splicing
bash run_module.sh

# Run pepMLM pre-processing module
echo "Run pepMLM pre-processing..."
cd ${analyses_dir}/pepMLM-preprocessing
bash scripts/run_module.sh