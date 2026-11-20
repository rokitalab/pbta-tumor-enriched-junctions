## GTEx splice event junction & target count matrix generation 
##
## Ryan Corbett
##
## Jan 2026

# This script performs the following:
# - loads GTEx rMATS files and filters for GTEx samples in the control cohort histologies
#   (brain regions, participants <40 years old, excluding samples removed in v10)
# - pulls splice event junction and target coordinates and normalizes against rMATS input reads
# - calculates mean junction/target normalized counts per GTEx subgroup

library(tidyverse)
library(qs2)
library(stringr)
library(data.table)

### Set directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "01-ctrl-rmats-processing")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")

source(file.path(analysis_dir, "util", "rmats-processing-functions.R"))

# Set file paths
gtex_se_file <- file.path(data_dir,
                          "gtex-rmats_merged_raw_SE.qs2")
gtex_ri_file <- file.path(data_dir,
                          "gtex-rmats_merged_raw_RI.qs2")
gtex_a3ss_file <- file.path(data_dir,
                            "gtex-rmats_merged_raw_A3SS.qs2")
gtex_a5ss_file <- file.path(data_dir,
                            "gtex-rmats_merged_raw_A5SS.qs2")
gtex_read_file <- file.path(data_dir,
                            "gtex_input_read_counts.tsv")

# GTEx metadata
hist_file <- file.path(data_dir,
                       "histologies.tsv")

# Control cohort histologies; defines the samples included in each control cohort
ctrl_hist_file <- file.path(root_dir, "analyses",
                            "00-create-cohort-histologies",
                            "results",
                            "control-cohort-histologies.tsv")

# Load control cohort histologies and get GTEx sample IDs
gtex_ids <- read_tsv(ctrl_hist_file) %>%
  dplyr::filter(cohort == "GTEx") %>%
  pull(sample_id)

# Load OPC hist and filter for GTEx control samples
gtex_brain_under40_hist <- read_tsv(hist_file) %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% gtex_ids) %>%
  # only need BS ID, subgroup columns
  dplyr::select(Kids_First_Biospecimen_ID, gtex_subgroup)

# Load rMATS SE results
se_df <- qs2::qs_read(gtex_se_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id),
                exonStart_0base = exonStart_0base + 1,
                upstreamES = upstreamES + 1,
                downstreamES = downstreamES + 1) %>%
  dplyr::filter(sample_id %in% gtex_brain_under40_hist$Kids_First_Biospecimen_ID)

# create average psi matrix (splice event ids x subgroups)
se_psi_mat <- create_psi_matrix(se_df,
                                event_type = "SE",
                                hist = gtex_brain_under40_hist,
                                group_col = "gtex_subgroup",
                                id_col = "Kids_First_Biospecimen_ID")

# Define SE junction and target IDs, and select relevant columns 
se_df <- define_junctions_targets(se_df,
                                  event_type = "SE")

# Extract all junctions and pivot longer
se_junction_df <- create_junction_df(se_df, 
                                     event_type = "SE")

# Load retained intron (RI) rMATS results, update sample ID, and filter for brain under40
ri_df <- qs2::qs_read(gtex_ri_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id),
                upstreamES = upstreamES + 1,
                downstreamES = downstreamES + 1) %>%
  dplyr::filter(sample_id %in% gtex_brain_under40_hist$Kids_First_Biospecimen_ID)

# create average psi matrix (splice event ids x subgroups)
ri_psi_mat <- create_psi_matrix(ri_df,
                                event_type = "RI",
                                hist = gtex_brain_under40_hist,
                                group_col = "gtex_subgroup",
                                id_col = "Kids_First_Biospecimen_ID")

# define columns specifying junction coordinates and retain only relevant columns
ri_df <- define_junctions_targets(ri_df,
                                  event_type = "RI")

# Build the long-form junction table, merging data from:
# upstream-intron 
# intron-downstream
# upstream-downstream junctions
ri_junction_df <- create_junction_df(ri_df,
                                     event_type = "RI")

# A3SS events, modify sample ID and filter out cell samples 
a3ss_df <- qs2::qs_read(gtex_a3ss_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id), 
                flankingES = flankingES + 1,
                longExonStart_0base = longExonStart_0base + 1,
                shortES = shortES + 1) %>%
  dplyr::filter(sample_id %in% gtex_brain_under40_hist$Kids_First_Biospecimen_ID)

# create average psi matrix (splice event ids x subgroups)
a3ss_psi_mat <- create_psi_matrix(a3ss_df,
                                event_type = "A3SS",
                                hist = gtex_brain_under40_hist,
                                group_col = "gtex_subgroup",
                                id_col = "Kids_First_Biospecimen_ID")

# define junction coordinates and filter for relevant columns
a3ss_df <- define_junctions_targets(a3ss_df,
                                    event_type = "A3SS")

# pivot longer for single row per unique sample & junction 
a3ss_junction_df <- create_junction_df(a3ss_df,
                                       event_type = "A3SS")

# A5SS events
a5ss_df <- qs2::qs_read(gtex_a5ss_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id), 
                flankingES = flankingES + 1,
                longExonStart_0base = longExonStart_0base + 1,
                shortES = shortES + 1) %>%
  dplyr::filter(sample_id %in% gtex_brain_under40_hist$Kids_First_Biospecimen_ID)

# create average psi matrix (splice event ids x subgroups)
a5ss_psi_mat <- create_psi_matrix(a5ss_df,
                                  event_type = "A5SS",
                                  hist = gtex_brain_under40_hist,
                                  group_col = "gtex_subgroup",
                                  id_col = "Kids_First_Biospecimen_ID")

# define junction coordinates and filter for relevant columns 
a5ss_df <- define_junctions_targets(a5ss_df,
                                    event_type = "A5SS")

# pivot longer for single row per unique sample & junction 
a5ss_junction_df <- create_junction_df(a5ss_df,
                                       event_type = "A5SS")

junction_list <- list("se" = se_junction_df,
                      "ri" = ri_junction_df,
                      "a3ss" = a3ss_junction_df,
                      "a5ss" = a5ss_junction_df)

# convert read cts and hist to data tables
gtex_hist_dt <- as.data.table(gtex_brain_under40_hist)

gtex_read_cts <- read_tsv(gtex_read_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id))

# create empty list to scores junction cpm matrices
junction_mat_list <- list()
junction_sd_list <- list()

# loop through junction dfs to generate matrices
for (event in names(junction_list)){
  
  # define current junction df
  junction_df <- junction_list[[event]]
  
  # run `generate_norm_junction_mat` to obtain desired normalized mean cpm matrix
  junction_mat_list[[event]] <- generate_norm_junction_mat(junction_df, 
                                                         gtex_read_cts,
                                                         gtex_hist_dt,
                                                         group_col = "gtex_subgroup",
                                                         id_col = "Kids_First_Biospecimen_ID")

  # run `generate_junction_sd_mat` to obtain desired normalized sd cpm matrix
  junction_sd_list[[event]] <- generate_junction_sd_mat(junction_df, 
                                                        gtex_read_cts,
                                                        gtex_hist_dt,
                                                        group_col = "gtex_subgroup",
                                                        id_col = "Kids_First_Biospecimen_ID")
  
}

# Merge junction matrices and filter for unique junction IDs
merged_norm_junction_mat <- junction_mat_list[["se"]] %>%
  bind_rows(junction_mat_list[["ri"]],
            junction_mat_list[["a3ss"]],
            junction_mat_list[["a5ss"]]) %>%
  distinct(junction, .keep_all = TRUE)
  
# Save merged junction output
arrow::write_parquet(merged_norm_junction_mat,
             file.path(results_dir,
                       "gtex-merged-norm-junction-ct-mat.parquet"),
             compression = "zstd")

# Merge junction sd matrices and filter for unique junction IDs
merged_junction_sd_mat <- junction_sd_list[["se"]] %>%
  bind_rows(junction_sd_list[["ri"]],
            junction_sd_list[["a3ss"]],
            junction_sd_list[["a5ss"]]) %>%
  distinct(junction, .keep_all = TRUE)

# Save merged junction output
arrow::write_parquet(merged_junction_sd_mat,
             file.path(results_dir,
                       "gtex-merged-norm-junction-sd-mat.parquet"),
             compression = "zstd")

# Merge PSI matrices 
merged_psi_mat <- se_psi_mat %>%
  dplyr::mutate(splicing_case = "SE") %>%
  bind_rows(ri_psi_mat %>%
              dplyr::mutate(splicing_case = "RI"),
            a3ss_psi_mat %>%
              dplyr::mutate(splicing_case = "A3SS"),
            a5ss_psi_mat %>%
              dplyr::mutate(splicing_case = "A5SS")) %>%
  dplyr::select(splice_id, splicing_case, everything())

arrow::write_parquet(merged_psi_mat,
             file.path(results_dir,
                       "gtex-merged-psi-mat.parquet"),
             compression = "zstd")

# Print session info
sessionInfo()
