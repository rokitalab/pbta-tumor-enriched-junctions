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
  # define subgroups (region + stage)
  dplyr::mutate(evodevo_subgroup = glue::glue("{primary_site}-{pathology_free_text_diagnosis}")) %>%
  # define broad groups (region + fetal OR postnatal)
  dplyr::mutate(evodevo_broadgroup = case_when(
    primary_site == "Forebrain" & grepl("Conception", pathology_free_text_diagnosis) ~ "Forebrain, fetal",
    primary_site == "Forebrain" & !grepl("Conception", pathology_free_text_diagnosis) ~ "Forebrain, postnatal",
    primary_site == "Hindbrain" & grepl("Conception", pathology_free_text_diagnosis) ~ "Hindbrain, fetal",
    primary_site == "Hindbrain" & !grepl("Conception", pathology_free_text_diagnosis) ~ "Hindbrain, postnatal"
  )) %>%
  # only need BS ID, subgroup columns
  dplyr::select(Kids_First_Biospecimen_ID, evodevo_subgroup, evodevo_broadgroup)

# Load rMATS SE results
se_df <- qs2::qs_read(se_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id),
                exonStart_0base = exonStart_0base + 1,
                upstreamES = upstreamES + 1,
                downstreamES = downstreamES + 1) %>%
  dplyr::filter(sample_id %in% evodevo_hist$Kids_First_Biospecimen_ID)

# create average psi matrix (splice event ids x subgroups)
se_subgroup_psi_mat <- create_psi_matrix(se_df,
                                event_type = "SE",
                                hist = evodevo_hist,
                                group_col = "evodevo_subgroup",
                                id_col = "Kids_First_Biospecimen_ID")

# create average psi matrix (splice event ids x broadgroups)
se_broadgroup_psi_mat <- create_psi_matrix(se_df,
                                         event_type = "SE",
                                         hist = evodevo_hist,
                                         group_col = "evodevo_broadgroup",
                                         id_col = "Kids_First_Biospecimen_ID")

# Define SE junction and target IDs, and select relevant columns 
se_df <- define_junctions_targets(se_df,
                                  event_type = "SE")

# Extract all junctions and pivot longer
se_junction_df <- create_junction_df(se_df, 
                                     event_type = "SE")

# Load retained intron (RI) rMATS results, update sample ID, and filter for brain under40
ri_df <- qs2::qs_read(ri_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id),
                upstreamES = upstreamES + 1,
                downstreamES = downstreamES + 1) %>%
  dplyr::filter(sample_id %in% evodevo_hist$Kids_First_Biospecimen_ID)

# create average psi matrix (splice event ids x subgroups)
ri_subgroup_psi_mat <- create_psi_matrix(ri_df,
                                        event_type = "RI",
                                        hist = evodevo_hist,
                                        group_col = "evodevo_subgroup",
                                        id_col = "Kids_First_Biospecimen_ID")

# create average psi matrix (splice event ids x subgroups)
ri_broadgroup_psi_mat <- create_psi_matrix(ri_df,
                                          event_type = "RI",
                                          hist = evodevo_hist,
                                          group_col = "evodevo_broadgroup",
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
a3ss_df <- qs2::qs_read(a3ss_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id), 
                flankingES = flankingES + 1,
                longExonStart_0base = longExonStart_0base + 1,
                shortES = shortES + 1) %>%
  dplyr::filter(sample_id %in% evodevo_hist$Kids_First_Biospecimen_ID)

# create average psi matrix (splice event ids x subgroups)
a3ss_subgroup_psi_mat <- create_psi_matrix(a3ss_df,
                                event_type = "A3SS",
                                hist = evodevo_hist,
                                group_col = "evodevo_subgroup",
                                id_col = "Kids_First_Biospecimen_ID")

# create average psi matrix (splice event ids x broadgroups)
a3ss_broadgroup_psi_mat <- create_psi_matrix(a3ss_df,
                                  event_type = "A3SS",
                                  hist = evodevo_hist,
                                  group_col = "evodevo_broadgroup",
                                  id_col = "Kids_First_Biospecimen_ID")

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

# create average psi matrix (splice event ids x subgroups)
a5ss_subgroup_psi_mat <- create_psi_matrix(a5ss_df,
                                  event_type = "A5SS",
                                  hist = evodevo_hist,
                                  group_col = "evodevo_subgroup",
                                  id_col = "Kids_First_Biospecimen_ID")

# create average psi matrix (splice event ids x subgroups)
a5ss_broadgroup_psi_mat <- create_psi_matrix(a5ss_df,
                                  event_type = "A5SS",
                                  hist = evodevo_hist,
                                  group_col = "evodevo_broadgroup",
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
evodevo_hist_dt <- as.data.table(evodevo_hist)

evodevo_read_cts <- read_tsv(read_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id))

# create empty list to scores junction cpm matrices
junction_subgroup_mat_list <- list()
junction_broadgroup_mat_list <- list()
junction_subgroup_sd_list <- list()
junction_broadgroup_sd_list <- list()

# loop through junction dfs to generate matrices
for (event in names(junction_list)){
  
  # define current junction df
  junction_df <- junction_list[[event]]
  
  # run `generate_norm_junction_mat` to obtain desired normalized mean cpm matrix
  
  # by subgroup:
  junction_subgroup_mat_list[[event]] <- generate_norm_junction_mat(junction_df, 
                                                                     evodevo_read_cts,
                                                                     evodevo_hist_dt,
                                                                     group_col = "evodevo_subgroup",
                                                                     id_col = "Kids_First_Biospecimen_ID")
  
  # by broadgroup:
  junction_broadgroup_mat_list[[event]] <- generate_norm_junction_mat(junction_df, 
                                                                    evodevo_read_cts,
                                                                    evodevo_hist_dt,
                                                                    group_col = "evodevo_broadgroup",
                                                                    id_col = "Kids_First_Biospecimen_ID")
  
  # run `generate_junction_sd_mat` to obtain desired normalized sd cpm matrix
  
  # by subgroup:
  junction_subgroup_sd_list[[event]] <- generate_junction_sd_mat(junction_df, 
                                                         evodevo_read_cts,
                                                         evodevo_hist_dt,
                                                         group_col = "evodevo_subgroup",
                                                         id_col = "Kids_First_Biospecimen_ID")
  
  # by broadgroup:
  junction_broadgroup_sd_list[[event]] <- generate_junction_sd_mat(junction_df, 
                                                                 evodevo_read_cts,
                                                                 evodevo_hist_dt,
                                                                 group_col = "evodevo_broadgroup",
                                                                 id_col = "Kids_First_Biospecimen_ID")
  
}

# Merge junction matrices and filter for unique junction IDs
merged_norm_subgroup_junction_mat <- junction_subgroup_mat_list[["se"]] %>%
  bind_rows(junction_subgroup_mat_list[["ri"]],
            junction_subgroup_mat_list[["a3ss"]],
            junction_subgroup_mat_list[["a5ss"]]) %>%
  distinct(junction, .keep_all = TRUE)

# Save merged junction output
qs2::qs_save(merged_norm_subgroup_junction_mat,
             file.path(results_dir,
                       "evodevo-merged-subgroup-norm-junction-ct-mat.qs2"))

merged_norm_broadgroup_junction_mat <- junction_broadgroup_mat_list[["se"]] %>%
  bind_rows(junction_broadgroup_mat_list[["ri"]],
            junction_broadgroup_mat_list[["a3ss"]],
            junction_broadgroup_mat_list[["a5ss"]]) %>%
  distinct(junction, .keep_all = TRUE)

# Save merged junction output
qs2::qs_save(merged_norm_broadgroup_junction_mat,
             file.path(results_dir,
                       "evodevo-merged-broadgroup-norm-junction-ct-mat.qs2"))

# Merge junction sd matrices and filter for unique junction IDs
merged_subgroup_junction_sd_mat <- junction_subgroup_sd_list[["se"]] %>%
  bind_rows(junction_subgroup_sd_list[["ri"]],
            junction_subgroup_sd_list[["a3ss"]],
            junction_subgroup_sd_list[["a5ss"]]) %>%
  distinct(junction, .keep_all = TRUE)

# Save merged junction output
qs2::qs_save(merged_subgroup_junction_sd_mat,
             file.path(results_dir,
                       "evodevo-merged-subgroup-norm-junction-sd-mat.qs2"))

merged_broadgroup_junction_sd_mat <- junction_broadgroup_sd_list[["se"]] %>%
  bind_rows(junction_broadgroup_sd_list[["ri"]],
            junction_broadgroup_sd_list[["a3ss"]],
            junction_broadgroup_sd_list[["a5ss"]]) %>%
  distinct(junction, .keep_all = TRUE)

# Save merged junction output
qs2::qs_save(merged_broadgroup_junction_sd_mat,
             file.path(results_dir,
                       "evodevo-merged-broadgroup-norm-junction-sd-mat.qs2"))

# Merge PSI matrices 
merged_subgroup_psi_mat <- se_subgroup_psi_mat %>%
  dplyr::mutate(splicing_case = "SE") %>%
  bind_rows(ri_subgroup_psi_mat %>%
              dplyr::mutate(splicing_case = "RI"),
            a3ss_subgroup_psi_mat %>%
              dplyr::mutate(splicing_case = "A3SS"),
            a5ss_subgroup_psi_mat %>%
              dplyr::mutate(splicing_case = "A5SS")) %>%
  dplyr::select(splice_id, splicing_case, everything())

qs2::qs_save(merged_subgroup_psi_mat,
             file.path(results_dir,
                       "evodevo-merged-subgroup-psi-mat.qs2"))


merged_broadgroup_psi_mat <- se_broadgroup_psi_mat %>%
  dplyr::mutate(splicing_case = "SE") %>%
  bind_rows(ri_broadgroup_psi_mat %>%
              dplyr::mutate(splicing_case = "RI"),
            a3ss_broadgroup_psi_mat %>%
              dplyr::mutate(splicing_case = "A3SS"),
            a5ss_broadgroup_psi_mat %>%
              dplyr::mutate(splicing_case = "A5SS")) %>%
  dplyr::select(splice_id, splicing_case, everything())

qs2::qs_save(merged_broadgroup_psi_mat,
             file.path(results_dir,
                       "evodevo-merged-broadgroup-psi-mat.qs2"))

# print session info
sessionInfo()
