## Pediatric brain splice event junction & target count matrix generation 
##
## Ryan Corbett
##
## Jan 2026

# This script performs the following:
# - loads pediatric brain rMATS files 
# - pulls splice event junction and target coordinates and normalizes against rMATS input reads
# - calculates mean junction/target normalized counts per brain region

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
se_file <- file.path(data_dir,
                     "normal_ped_brain-rmats_merged_raw_SE.qs2")
ri_file <- file.path(data_dir,
                     "normal_ped_brain-rmats_merged_raw_RI.qs2")
a3ss_file <- file.path(data_dir,
                       "normal_ped_brain-rmats_merged_raw_A3SS.qs2")
a5ss_file <- file.path(data_dir,
                       "normal_ped_brain-rmats_merged_raw_A5SS.qs2")
read_file <- file.path(data_dir,
                       "normal_ped_brain_input_read_counts.tsv")

hist_file <- file.path(data_dir,
                       "ped-normal-brain-histologies.tsv")

# Load hist
hist <- read_tsv(hist_file) %>%
  # filter out tumor-infiltrated pons
  dplyr::filter(sample_id != "7316-7585") %>%
  # only need sample ID, primary_site columns
  dplyr::select(sample_id, primary_site)

# Load rMATS SE results
se_df <- qs2::qs_read(se_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id),
                exonStart_0base = exonStart_0base + 1,
                upstreamES = upstreamES + 1,
                downstreamES = downstreamES + 1) %>%
  dplyr::filter(sample_id %in% hist$sample_id)

# create average psi matrix (splice event ids x subgroups)
se_psi_mat <- create_psi_matrix(se_df,
                                event_type = "SE",
                                hist = hist,
                                group_col = "primary_site",
                                id_col = "sample_id")

# Define SE junction and target IDs, and select relevant columns 
se_df <- define_junctions_targets(se_df,
                                  event_type = "SE")

# Extract all junctions and pivot longer
se_junction_df <- create_junction_df(se_df, 
                                     event_type = "SE")

# Load retained intron (RI) rMATS results, update sample ID
ri_df <- qs2::qs_read(ri_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id),
                upstreamES = upstreamES + 1,
                downstreamES = downstreamES + 1) %>%
  dplyr::filter(sample_id %in% hist$sample_id)

# create average psi matrix (splice event ids x subgroups)
ri_psi_mat <- create_psi_matrix(ri_df,
                                event_type = "RI",
                                hist = hist,
                                group_col = "primary_site",
                                id_col = "sample_id")

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
a3ss_df <- qs2::qs_read(a3ss_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id), 
                flankingES = flankingES + 1,
                longExonStart_0base = longExonStart_0base + 1,
                shortES = shortES + 1) %>%
  dplyr::filter(sample_id %in% hist$sample_id)

# create average psi matrix (splice event ids x subgroups)
a3ss_psi_mat <- create_psi_matrix(a3ss_df,
                                event_type = "A3SS",
                                hist = hist,
                                group_col = "primary_site",
                                id_col = "sample_id")

# define junction coordinates and filter for relevant columns
a3ss_df <- define_junctions_targets(a3ss_df,
                                    event_type = "A3SS")

# pivot longer for single row per unique sample & junction 
a3ss_junction_df <- create_junction_df(a3ss_df,
                                       event_type = "A3SS")

# A5SS events
a5ss_df <- qs2::qs_read(a5ss_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id), 
                flankingES = flankingES + 1,
                longExonStart_0base = longExonStart_0base + 1,
                shortES = shortES + 1) %>%
  dplyr::filter(sample_id %in% hist$sample_id)

# create average psi matrix (splice event ids x subgroups)
a5ss_psi_mat <- create_psi_matrix(a5ss_df,
                                  event_type = "A5SS",
                                  hist = hist,
                                  group_col = "primary_site",
                                  id_col = "sample_id")

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
hist_dt <- as.data.table(hist)

read_cts <- read_tsv(read_file) %>%
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
                                                           read_cts,
                                                           hist_dt,
                                                           group_col = "primary_site",
                                                           id_col = "sample_id")
  
  # run `generate_junction_sd_mat` to obtain desired normalized sd cpm matrix
  junction_sd_list[[event]] <- generate_junction_sd_mat(junction_df, 
                                                        read_cts,
                                                        hist_dt,
                                                        group_col = "primary_site",
                                                        id_col = "sample_id")
  
}

# Merge junction matrices and filter for unique junction IDs
merged_norm_junction_mat <- junction_mat_list[["se"]] %>%
  bind_rows(junction_mat_list[["ri"]],
            junction_mat_list[["a3ss"]],
            junction_mat_list[["a5ss"]]) %>%
  distinct(junction, .keep_all = TRUE)

# Save merged junction output
qs2::qs_save(merged_norm_junction_mat,
             file.path(results_dir,
                       "normal-pedbrain-merged-norm-junction-ct-mat.qs2"))

# Merge junction sd matrices and filter for unique junction IDs
merged_junction_sd_mat <- junction_sd_list[["se"]] %>%
  bind_rows(junction_sd_list[["ri"]],
            junction_sd_list[["a3ss"]],
            junction_sd_list[["a5ss"]]) %>%
  distinct(junction, .keep_all = TRUE)

# Save merged junction output
qs2::qs_save(merged_junction_sd_mat,
             file.path(results_dir,
                       "normal-pedbrain-merged-norm-junction-sd-mat.qs2"))

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

qs2::qs_save(merged_psi_mat,
             file.path(results_dir,
                       "normal-pedbrain-merged-psi-mat.qs2"))

# print session info
sessionInfo()
