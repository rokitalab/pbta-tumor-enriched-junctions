## Evo-Devo splice event junction & target count matrix generation 
##
## Ryan Corbett
##
## Jan 2026

# This script performs the following:
# - loads Evo-Devo rMATS files and filters for brain regions (all ages)
# - pulls splice event junction and target coordinates and normalizes against rMATS input reads
# - calculates mean junction/target normalized counts per Evo-Devo subgroup

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
                     "evodevo-rmats_merged_raw_SE.qs2")
ri_file <- file.path(data_dir,
                     "evodevo-rmats_merged_raw_RI.qs2")
a3ss_file <- file.path(data_dir,
                       "evodevo-rmats_merged_raw_A3SS.qs2")
a5ss_file <- file.path(data_dir,
                       "evodevo-rmats_merged_raw_A5SS.qs2")
read_file <- file.path(data_dir,
                       "evodevo_input_read_counts.tsv")

hist_file <- file.path(data_dir,
                       "evodevo-histologies.tsv")

# Load evodevo hist and filter for brain regions
evodevo_hist <- read_tsv(hist_file) %>%
  # filter for brain samples, filter out older life stages
  dplyr::filter(primary_site %in% c("Forebrain", "Hindbrain"),
                !pathology_free_text_diagnosis %in% c("Middle Adult",
                                                      "Elderly")) %>%
  dplyr::mutate(evodevo_subgroup = glue::glue("{primary_site}-{pathology_free_text_diagnosis}")) %>%
  # only need BS ID, subgroup columns
  dplyr::select(Kids_First_Biospecimen_ID, evodevo_subgroup)

# Load rMATS SE results
se_df <- qs2::qs_read(se_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id),
                exonStart_0base = exonStart_0base + 1,
                upstreamES = upstreamES + 1,
                downstreamES = downstreamES + 1) %>%
  dplyr::filter(sample_id %in% evodevo_hist$Kids_First_Biospecimen_ID)

# Define SE junction and target IDs, and select relevant columns 
se_df <- define_junctions_targets(se_df,
                                  event_type = "SE")

# get median target counts 
se_target_df <- create_target_df(se_df)

# Extract all junctions and pivot longer
se_junction_df <- create_junction_df(se_df, 
                                     event_type = "SE")

# Load retained intron (RI) rMATS results, update sample ID, and filter for brain under40
ri_df <- qs2::qs_read(ri_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id),
                upstreamES = upstreamES + 1,
                downstreamES = downstreamES + 1) %>%
  dplyr::filter(sample_id %in% evodevo_hist$Kids_First_Biospecimen_ID)

# define columns specifying junction coordinates and retain only relevant columns
ri_df <- define_junctions_targets(ri_df,
                                  event_type = "RI")

# get median intron counts 
ri_target_df <- create_target_df(ri_df)

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
  dplyr::filter(sample_id %in% evodevo_hist$Kids_First_Biospecimen_ID)

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
  dplyr::filter(sample_id %in% evodevo_hist$Kids_First_Biospecimen_ID)

# define junction coordinates and filter for relevant columns 
a5ss_df <- define_junctions_targets(a5ss_df,
                                    event_type = "A5SS")

# pivot longer for single row per unique sample & junction 
a5ss_junction_df <- create_junction_df(a5ss_df,
                                       event_type = "A5SS")

# create lists of junction and target dfs
target_list <- list("se"= se_target_df,
                    "ri" = ri_target_df)

junction_list <- list("se" = se_junction_df,
                      "ri" = ri_junction_df,
                      "a3ss" = a3ss_junction_df,
                      "a5ss" = a5ss_junction_df)

# convert read cts and hist to data tables
evodevo_hist_dt <- as.data.table(evodevo_hist)

evodevo_read_cts <- read_tsv(read_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id))

# loops through target dfs to generate matrices
for (event in names(target_list)){
  
  # define current target df
  target_df <- target_list[[event]]
  
  # run `generate_norm_target_mat` funtion to obtain desired matrix
  norm_target_mat <- generate_norm_target_mat(target_df, 
                                              evodevo_read_cts,
                                              evodevo_hist_dt,
                                              group_col = "evodevo_subgroup",
                                              id_col = "Kids_First_Biospecimen_ID")
  
  # save mat
  qs2::qs_save(norm_target_mat,
               file.path(results_dir, 
                         glue::glue("evodevo-{event}-norm-target-ct-mat.qs2")))
  
}

# create empty list to scores junction cpm matrices
junction_mat_list <- list()

# loop through junction dfs to generate matrices
for (event in names(junction_list)){
  
  # define current junction df
  junction_df <- junction_list[[event]]
  
  # run `generate_norm_junction_mat` to obtain desired normalized mean cpm matrix
  junction_mat_list[[event]] <- generate_norm_junction_mat(junction_df, 
                                                           evodevo_read_cts,
                                                           evodevo_hist_dt,
                                                           group_col = "evodevo_subgroup",
                                                           id_col = "Kids_First_Biospecimen_ID")
  
  # save mat
  qs2::qs_save(junction_mat_list[[event]],
               file.path(results_dir, 
                         glue::glue("evodevo-{event}-norm-junction-ct-mat.qs2")))
  
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
                       "evodevo-merged-norm-junction-ct-mat.qs2"))

# print session info
sessionInfo()
