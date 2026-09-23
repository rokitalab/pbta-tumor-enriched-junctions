## Evo-Devo splice event junction & target count matrix generation 
##
## Ryan Corbett
##
## Jan 2026

# This script performs the following:
# - loads Evo-Devo rMATS files and filters for brain regions (all ages)
# - pulls splice event junction and target coordinates and normalizes against rMATS input reads
# - calculates matrices for postnatal brain region/stage groups only
# - calculates prenatal brain matrices grouped by region and developmental stage:
#   4-5 PCW early embryonic, 6-7 PCW late embryonic, 8-12 PCW early fetal,
#   13-16 PCW early mid-fetal, and 17-19 PCW late mid-fetal

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

# Load evodevo hist, filter for brain regions, and define sample groupings
evodevo_hist <- read_tsv(hist_file) %>%
  # filter for brain samples, filter out older life stages
  dplyr::filter(primary_site %in% c("Forebrain", "Hindbrain"),
                !pathology_free_text_diagnosis %in% c("Middle Adult",
                                                      "Elderly")) %>%
  # define subgroups (region + stage)
  dplyr::mutate(evodevo_subgroup = glue::glue("{primary_site}-{pathology_free_text_diagnosis}")) %>%
  # Extract gestational age only for prenatal samples and assign developmental-stage bins.
  dplyr::mutate(
    weeks_post_conception = case_when(
      grepl("Week Post Conception", pathology_free_text_diagnosis) ~
        readr::parse_number(pathology_free_text_diagnosis),
      TRUE ~ NA_real_
    ),
    prenatal_week_bin = case_when(
      dplyr::between(weeks_post_conception, 4, 5) ~ "4-5pcw_early_embryonic",
      dplyr::between(weeks_post_conception, 6, 7) ~ "6-7pcw_late_embryonic",
      dplyr::between(weeks_post_conception, 8, 12) ~ "8-12pcw_early_fetal",
      dplyr::between(weeks_post_conception, 13, 16) ~ "13-16pcw_early_midfetal",
      dplyr::between(weeks_post_conception, 17, 19) ~ "17-19pcw_late_midfetal",
      TRUE ~ NA_character_
    ),
    evodevo_prenatal_week_group = if_else(
      !is.na(prenatal_week_bin),
      str_c(primary_site, prenatal_week_bin, sep = "-"),
      NA_character_
    ),
    evodevo_postnatal_group = if_else(
      !grepl("Conception", pathology_free_text_diagnosis),
      evodevo_subgroup,
      NA_character_
    )
  ) %>%
  # Retain a derived metadata file without altering the source metadata.
  {readr::write_tsv(., file.path(results_dir,
                                 "evodevo-brain-prenatal-week-binned-metadata.tsv")); .} %>%
  # Retain the sample ID and the two requested matrix groupings.
  dplyr::select(Kids_First_Biospecimen_ID, evodevo_postnatal_group,
                evodevo_prenatal_week_group)

# Prenatal samples only; used for matrices grouped by brain region and developmental stage.
evodevo_prenatal_hist <- evodevo_hist %>%
  dplyr::filter(!is.na(evodevo_prenatal_week_group))
evodevo_postnatal_hist <- evodevo_hist %>%
  dplyr::filter(!is.na(evodevo_postnatal_group))

# Load rMATS SE results
se_df <- qs2::qs_read(se_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id),
                exonStart_0base = exonStart_0base + 1,
                upstreamES = upstreamES + 1,
                downstreamES = downstreamES + 1) %>%
  dplyr::filter(sample_id %in% evodevo_hist$Kids_First_Biospecimen_ID)

# create average PSI matrix for postnatal samples by region and stage
se_postnatal_psi_mat <- create_psi_matrix(se_df %>%
                                            dplyr::filter(sample_id %in% evodevo_postnatal_hist$Kids_First_Biospecimen_ID),
                                          event_type = "SE",
                                          hist = evodevo_postnatal_hist,
                                          group_col = "evodevo_postnatal_group",
                                          id_col = "Kids_First_Biospecimen_ID")

# create median PSI matrix for prenatal samples by region and developmental stage
se_prenatal_week_psi_mat <- create_psi_matrix(se_df %>%
                                                dplyr::filter(sample_id %in% evodevo_prenatal_hist$Kids_First_Biospecimen_ID),
                                              event_type = "SE",
                                              hist = evodevo_prenatal_hist,
                                              group_col = "evodevo_prenatal_week_group",
                                              id_col = "Kids_First_Biospecimen_ID",
                                              aggregation_fun = median)

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

# create average PSI matrix for postnatal samples by region and stage
ri_postnatal_psi_mat <- create_psi_matrix(ri_df %>%
                                            dplyr::filter(sample_id %in% evodevo_postnatal_hist$Kids_First_Biospecimen_ID),
                                          event_type = "RI",
                                          hist = evodevo_postnatal_hist,
                                          group_col = "evodevo_postnatal_group",
                                          id_col = "Kids_First_Biospecimen_ID")

# create median PSI matrix for prenatal samples by region and developmental stage
ri_prenatal_week_psi_mat <- create_psi_matrix(ri_df %>%
                                                dplyr::filter(sample_id %in% evodevo_prenatal_hist$Kids_First_Biospecimen_ID),
                                              event_type = "RI",
                                              hist = evodevo_prenatal_hist,
                                              group_col = "evodevo_prenatal_week_group",
                                              id_col = "Kids_First_Biospecimen_ID",
                                              aggregation_fun = median)

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

# create average PSI matrix for postnatal samples by region and stage
a3ss_postnatal_psi_mat <- create_psi_matrix(a3ss_df %>%
                                              dplyr::filter(sample_id %in% evodevo_postnatal_hist$Kids_First_Biospecimen_ID),
                                            event_type = "A3SS",
                                            hist = evodevo_postnatal_hist,
                                            group_col = "evodevo_postnatal_group",
                                            id_col = "Kids_First_Biospecimen_ID")

# create median PSI matrix for prenatal samples by region and developmental stage
a3ss_prenatal_week_psi_mat <- create_psi_matrix(a3ss_df %>%
                                                  dplyr::filter(sample_id %in% evodevo_prenatal_hist$Kids_First_Biospecimen_ID),
                                                event_type = "A3SS",
                                                hist = evodevo_prenatal_hist,
                                                group_col = "evodevo_prenatal_week_group",
                                                id_col = "Kids_First_Biospecimen_ID",
                                                aggregation_fun = median)

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

# create average PSI matrix for postnatal samples by region and stage
a5ss_postnatal_psi_mat <- create_psi_matrix(a5ss_df %>%
                                              dplyr::filter(sample_id %in% evodevo_postnatal_hist$Kids_First_Biospecimen_ID),
                                            event_type = "A5SS",
                                            hist = evodevo_postnatal_hist,
                                            group_col = "evodevo_postnatal_group",
                                            id_col = "Kids_First_Biospecimen_ID")

# create median PSI matrix for prenatal samples by region and developmental stage
a5ss_prenatal_week_psi_mat <- create_psi_matrix(a5ss_df %>%
                                                  dplyr::filter(sample_id %in% evodevo_prenatal_hist$Kids_First_Biospecimen_ID),
                                                event_type = "A5SS",
                                                hist = evodevo_prenatal_hist,
                                                group_col = "evodevo_prenatal_week_group",
                                                id_col = "Kids_First_Biospecimen_ID",
                                                aggregation_fun = median)

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
evodevo_prenatal_hist_dt <- as.data.table(evodevo_prenatal_hist)
evodevo_postnatal_hist_dt <- as.data.table(evodevo_postnatal_hist)

evodevo_read_cts <- read_tsv(read_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id))

# create empty list to scores junction cpm matrices
junction_postnatal_mat_list <- list()
junction_prenatal_week_mat_list <- list()
junction_postnatal_sd_list <- list()
junction_prenatal_week_sd_list <- list()

# loop through junction dfs to generate matrices
for (event in names(junction_list)){
  
  # define current junction df
  junction_df <- junction_list[[event]]
  prenatal_junction_df <- junction_df[
    sample_id %in% evodevo_prenatal_hist_dt$Kids_First_Biospecimen_ID
  ]
  postnatal_junction_df <- junction_df[
    sample_id %in% evodevo_postnatal_hist_dt$Kids_First_Biospecimen_ID
  ]
  
  # run `generate_norm_junction_mat` to obtain normalized CPM matrices
  
  # by postnatal region and stage:
  junction_postnatal_mat_list[[event]] <- generate_norm_junction_mat(postnatal_junction_df,
                                                                      evodevo_read_cts,
                                                                      evodevo_postnatal_hist_dt,
                                                                      group_col = "evodevo_postnatal_group",
                                                                      id_col = "Kids_First_Biospecimen_ID")

  # by prenatal region and developmental stage:
  junction_prenatal_week_mat_list[[event]] <- generate_norm_junction_mat(prenatal_junction_df,
                                                                           evodevo_read_cts,
                                                                           evodevo_prenatal_hist_dt,
                                                                           group_col = "evodevo_prenatal_week_group",
                                                                           id_col = "Kids_First_Biospecimen_ID",
                                                                           aggregation_fun = median)
  
  # run `generate_junction_sd_mat` to obtain desired normalized sd cpm matrix
  
  # by postnatal region and stage:
  junction_postnatal_sd_list[[event]] <- generate_junction_sd_mat(postnatal_junction_df,
                                                                    evodevo_read_cts,
                                                                    evodevo_postnatal_hist_dt,
                                                                    group_col = "evodevo_postnatal_group",
                                                                    id_col = "Kids_First_Biospecimen_ID")

  # by prenatal region and developmental stage:
  junction_prenatal_week_sd_list[[event]] <- generate_junction_sd_mat(prenatal_junction_df,
                                                                        evodevo_read_cts,
                                                                        evodevo_prenatal_hist_dt,
                                                                        group_col = "evodevo_prenatal_week_group",
                                                                        id_col = "Kids_First_Biospecimen_ID")
  
}

# Merge and save mean junction CPM matrix for postnatal region/stage groups.
merged_norm_postnatal_junction_mat <- junction_postnatal_mat_list[["se"]] %>%
  bind_rows(junction_postnatal_mat_list[["ri"]],
            junction_postnatal_mat_list[["a3ss"]],
            junction_postnatal_mat_list[["a5ss"]]) %>%
  distinct(junction, .keep_all = TRUE)

qs2::qs_save(merged_norm_postnatal_junction_mat,
             file.path(results_dir,
                       "evodevo-merged-postnatal-norm-junction-ct-mat.qs2"))

# Merge and save median junction CPM matrix for prenatal region/developmental-stage groups.
merged_norm_prenatal_week_junction_mat <- junction_prenatal_week_mat_list[["se"]] %>%
  bind_rows(junction_prenatal_week_mat_list[["ri"]],
            junction_prenatal_week_mat_list[["a3ss"]],
            junction_prenatal_week_mat_list[["a5ss"]]) %>%
  distinct(junction, .keep_all = TRUE)

qs2::qs_save(merged_norm_prenatal_week_junction_mat,
             file.path(results_dir,
                       "evodevo-merged-prenatal-week-binned-norm-junction-ct-mat.qs2"))

# Merge and save junction CPM SD matrix for postnatal region/stage groups.
merged_postnatal_junction_sd_mat <- junction_postnatal_sd_list[["se"]] %>%
  bind_rows(junction_postnatal_sd_list[["ri"]],
            junction_postnatal_sd_list[["a3ss"]],
            junction_postnatal_sd_list[["a5ss"]]) %>%
  distinct(junction, .keep_all = TRUE)

qs2::qs_save(merged_postnatal_junction_sd_mat,
             file.path(results_dir,
                       "evodevo-merged-postnatal-norm-junction-sd-mat.qs2"))

# Merge and save junction CPM SD matrix for prenatal region/developmental-stage groups.
merged_prenatal_week_junction_sd_mat <- junction_prenatal_week_sd_list[["se"]] %>%
  bind_rows(junction_prenatal_week_sd_list[["ri"]],
            junction_prenatal_week_sd_list[["a3ss"]],
            junction_prenatal_week_sd_list[["a5ss"]]) %>%
  distinct(junction, .keep_all = TRUE)

qs2::qs_save(merged_prenatal_week_junction_sd_mat,
             file.path(results_dir,
                       "evodevo-merged-prenatal-week-binned-norm-junction-sd-mat.qs2"))

# Merge and save PSI matrix for postnatal region/stage groups.
merged_postnatal_psi_mat <- se_postnatal_psi_mat %>%
  dplyr::mutate(splicing_case = "SE") %>%
  bind_rows(ri_postnatal_psi_mat %>%
              dplyr::mutate(splicing_case = "RI"),
            a3ss_postnatal_psi_mat %>%
              dplyr::mutate(splicing_case = "A3SS"),
            a5ss_postnatal_psi_mat %>%
              dplyr::mutate(splicing_case = "A5SS")) %>%
  dplyr::select(splice_id, splicing_case, everything())

qs2::qs_save(merged_postnatal_psi_mat,
             file.path(results_dir,
                       "evodevo-merged-postnatal-psi-mat.qs2"))

merged_prenatal_week_psi_mat <- se_prenatal_week_psi_mat %>%
  dplyr::mutate(splicing_case = "SE") %>%
  bind_rows(ri_prenatal_week_psi_mat %>%
              dplyr::mutate(splicing_case = "RI"),
            a3ss_prenatal_week_psi_mat %>%
              dplyr::mutate(splicing_case = "A3SS"),
            a5ss_prenatal_week_psi_mat %>%
              dplyr::mutate(splicing_case = "A5SS")) %>%
  dplyr::select(splice_id, splicing_case, everything())

qs2::qs_save(merged_prenatal_week_psi_mat,
             file.path(results_dir,
                       "evodevo-merged-prenatal-week-binned-psi-mat.qs2"))

# print session info
sessionInfo()
