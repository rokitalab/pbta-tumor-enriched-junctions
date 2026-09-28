## Calculate TESJ differential splicing
##
## Ryan Corbett
##
## Feb 2026

# This script performs the following:
# - Loads PBTA splice event data frame and control cohort PSI matrices
# - calculates dPSIs between tumor event PSIs and control cohort PSIs
# - Determine consensus differential splicing event associated with each junction, when available
# - Filter TESJs for those associated with differential splicing event

library(tidyverse)
library(qs2)
library(stringr)
library(data.table)

### Set directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "04-tej-differential-splicing")
results_dir <- file.path(analysis_dir, "results")

# Set file paths
enr_jc_splice_events_file <- file.path(results_dir, 
                                       "tumor-enriched-oncofetal-junction-splice-events.qs2")

gtex_psi_file <- file.path(root_dir, "analyses",
                           "01-ctrl-rmats-processing",
                           "results",
                           "gtex-merged-psi-mat.qs2")

evodevo_psi_file <- file.path(root_dir, "analyses",
                              "01-ctrl-rmats-processing",
                              "results",
                              "evodevo-merged-postnatal-collapsed-psi-mat.qs2")

pedbrain_psi_file <- file.path(root_dir, "analyses",
                               "01-ctrl-rmats-processing",
                               "results",
                               "normal-pedbrain-merged-psi-mat.qs2")

enr_jcs_annotated_file <- file.path(root_dir, "analyses",
                                    "03-classify-tejs",
                                    "results",
                                    "tumor-enriched-oncofetal-splice-junctions-domain-expr-annotated.tsv.gz")

# Wrangle data
enr_jc_splice_event_df <- qs2::qs_read(enr_jc_splice_events_file)

# Retain the original tumor-enrichment classification in each module output.
# The splice-event table is generated from rMATS records, so it does not
# itself contain this junction-level annotation.
enr_jcs_annotated_df <- read_tsv(enr_jcs_annotated_file)
criteria_by_id <- enr_jcs_annotated_df %>%
  dplyr::transmute(id = glue::glue("{sample_id}-{junction}"), criteria) %>%
  distinct()

gtex_psi_mat <- qs2::qs_read(gtex_psi_file)

evodevo_psi_mat <- qs2::qs_read(evodevo_psi_file)

pedbrain_psi_mat <- qs2::qs_read(pedbrain_psi_file)

# calculate mean dPSIs between tumor and normals
enr_jc_dpsi_df <- enr_jc_splice_event_df %>%
  # join control cohort PSIs
  left_join(gtex_psi_mat,
            by = c("splice_id", "splicing_case")) %>%
  left_join(evodevo_psi_mat,
            by = c("splice_id", "splicing_case")) %>%
  left_join(pedbrain_psi_mat,
            by = c("splice_id", "splicing_case")) %>%
  # calculate dPSIs for each tumor-control cohort pair
  dplyr::mutate(across(contains(c("Brain", "brain", "Cerebellum", "Cortex", "Pituitary", "Pons")), ~ psi - ., .names = "dPSI_{.col}")) %>%
  # calculate mean psi across all control cohorts
  dplyr::mutate(mean_psi = rowMeans(across(contains(c("Brain", "brain", "Cerebellum", "Cortex", "Pituitary", "Pons")) & !contains("dPSI")), na.rm = TRUE)) %>%
  # calculate mean dPSI
  dplyr::mutate(mean_dpsi = rowMeans(across(contains(c("dPSI"))), na.rm = TRUE)) %>%
  # create sample-junction
  dplyr::mutate(id = glue::glue("{sample_id}-{junction}")) %>%
  # define junction differential splicing classifications
  dplyr::mutate(type_resolved = case_when(
    # If junction type is exon inclusion and dpsi and/or psi support exon inclusion 
    type == "exon inclusion" & (mean_dpsi > 0.1 | (mean_psi < 0.025 & psi > 0.1)) ~ "exon inclusion",
    # If junction type is exon skipping and dpsi and/or psi support exon skipping 
    type == "exon skipping" & (mean_dpsi < -0.1 | (mean_psi > 0.975 & psi < 0.9)) ~ "exon skipping",
    # If junction type is intron retention and dpsi and/or psi support intron retention
    type == "intron retention" & (mean_dpsi > 0.1 | (mean_psi < 0.025 & psi > 0.1)) ~ "intron retention",
    # If junction type is Alt SS and dpsi and/or psi support event
    type == "A3SS+" & (mean_dpsi > 0.1 | (mean_psi < 0.025 & psi > 0.1)) ~ "A3SS+",
    type == "A3SS-" & (mean_dpsi < -0.1 | (mean_psi > 0.975 & psi < 0.9)) ~ "A3SS-",
    type == "A5SS+" & (mean_dpsi > 0.1 | (mean_psi < 0.025 & psi > 0.1)) ~ "A5SS+",
    type == "A5SS-" & (mean_dpsi < -0.1 | (mean_psi > 0.975 & psi < 0.9)) ~ "A5SS-",
    is.nan(mean_dpsi) ~ "Unknown",
    TRUE ~ "No diff splice event"
  ))

### Get sample consensus diff splice events

# First, identify all sample-junction pairs with single diff splice event
ids_1_n <- enr_jc_dpsi_df %>%
  dplyr::filter(!is.nan(mean_dpsi),
                type_resolved != "No diff splice event") %>%
  dplyr::count(id) %>%
  dplyr::filter(n == 1)  %>%
  pull(id)

# extract from full df
resolved_1_n <- enr_jc_dpsi_df %>%
  dplyr::filter(id %in% ids_1_n,
                !is.nan(mean_dpsi),
                type_resolved != "No diff splice event") %>%
  dplyr::select(id, sample_id, 
                junction, splice_id, 
                psi, mean_dpsi,
                type_resolved) %>%
  dplyr::rename(event_type_sample = type_resolved)

# Next, identify sample-junction pairs associated with diff splice events of the same type
ids_1_n_collapsed <- enr_jc_dpsi_df %>% 
  dplyr::filter(!type_resolved %in% c("No diff splice event",
                                      "Unknown"),
                !id %in% ids_1_n) %>%
  distinct(id, type_resolved) %>%
  count(id) %>%
  dplyr::filter(n == 1) %>%
  pull(id)
  
# extract ids from full data frame
resolved_1_n_collapsed <- enr_jc_dpsi_df %>%
  dplyr::filter(id %in% ids_1_n_collapsed,
                !is.nan(mean_dpsi),
                type_resolved != "No diff splice event") %>%
  dplyr::select(id, sample_id, 
                junction, splice_id, 
                psi, mean_dpsi,
                type_resolved) %>%
  # order by descending |dPSI|
  dplyr::arrange(desc(abs(mean_dpsi))) %>%
  # only retain splice event with max |dPSI|
  distinct(id, .keep_all = TRUE) %>%
  dplyr::rename(event_type_sample = type_resolved)

# Next, take splice event with greatest mean_dPSI and assign as resolved splice event
resolved_top_dpsi <- enr_jc_dpsi_df %>% 
  dplyr::filter(!type_resolved %in% c("No diff splice event",
                                      "Unknown"),
                !id %in% c(ids_1_n,
                           ids_1_n_collapsed)) %>%
  dplyr::arrange(desc(abs(mean_dpsi))) %>%
  distinct(id, .keep_all = TRUE) %>%
  dplyr::select(id, sample_id, 
                junction, splice_id, 
                psi, mean_dpsi,
                type_resolved) %>%
  dplyr::rename(event_type_sample = type_resolved)

# Next, take all splice events with no records in any control cohort
resolved_na_dpsi <- enr_jc_dpsi_df %>% 
  dplyr::filter(type_resolved == "Unknown",
                !id %in% c(ids_1_n,
                           ids_1_n_collapsed,
                           resolved_top_dpsi$id)) %>%
  # determine if PSI supports greater than 10% reads supporting event (psi > 0.1 for inclusion events, psi < 0.9 for loss events)
  dplyr::mutate(type_resolved = case_when(
    type %in% c("A3SS+", "A5SS+", "exon inclusion") & psi > 0.1 ~ type,
    type == "intron retention" & psi > 0.1 ~ "intron retention",
    type %in% c("A3SS-", "A5SS-", "exon skipping") & psi < 0.9 ~ type,
    TRUE ~ type_resolved
  )) %>%
  dplyr::filter(type_resolved != "Unknown") %>%
  dplyr::select(id, sample_id, 
                junction, splice_id, 
                psi, mean_dpsi,
                type,
                type_resolved) %>%
  # calculate deviation to rank splice event psis
  dplyr::mutate(deviation = case_when(
    type_resolved %in% c("A3SS+", "A5SS+", 
                         "exon inclusion", "intron retention") ~ psi,
    type_resolved %in% c("A3SS-", "A5SS-", "exon skipping") ~ 1 - psi
  )) %>%
  # arrange by decreasing deviation
  dplyr::arrange(desc(deviation)) %>%
  distinct(id, .keep_all = TRUE) %>%
  dplyr::rename(event_type_sample = type_resolved)

# merge sample junction consensus diff splice events
merged_resolved_df <- resolved_1_n %>%
  bind_rows(resolved_1_n_collapsed,
            resolved_top_dpsi,
            resolved_na_dpsi)

### define junction-level consensus diff splice event

# Get consensus diff splice event type for each junction
consensus_jc_df <- merged_resolved_df %>%
  group_by(junction, event_type_sample) %>%
  # get N diff splice event types and median mean dPSI per junction
  summarise(type_n = n(),
            median_dpsi = median(mean_dpsi)) %>%
  # Order junctions by desc freq of diff splice event types, and median dPSI
  dplyr::arrange(junction,
                 desc(type_n),
                 desc(abs(median_dpsi))) %>%
  # keep diff splice event with highest N (or higher median dPSI) as consensus call
  distinct(junction, .keep_all = TRUE) %>%
  dplyr::rename(consensus_jc_event_type = event_type_sample)

merged_resolved_df <- merged_resolved_df %>%
  left_join(consensus_jc_df %>%
              dplyr::select(junction, 
                            consensus_jc_event_type)) %>%
  left_join(criteria_by_id, by = "id")

# write tsv
write_tsv(merged_resolved_df,
          file.path(results_dir,
                    "tumor-enriched-oncofetal-junction-diff-splice-event-annotation.tsv.gz"))

# load domain-annotated, expression filtered junction df
enr_jcs_annotated_filtered_df <- enr_jcs_annotated_df %>%
  dplyr::mutate(id = glue::glue("{sample_id}-{junction}")) %>%
  dplyr::filter(id %in% merged_resolved_df$id) %>%
  left_join(merged_resolved_df %>%
              dplyr::select(id, event_type_sample, consensus_jc_event_type)) %>%
  # rename sample_id col to BS_ID for downstream analyses
  dplyr::rename(Kids_First_Biospecimen_ID = sample_id)

# write to output
qs2::qs_save(enr_jcs_annotated_filtered_df,
              file.path(results_dir,
                        "tumor-enriched-oncofetal-diff-splice-junctions.qs2"))
