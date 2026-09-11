## Identify tumor-enriched junctions in PBTA
##
## Ryan Corbett
##
## Jan 2025
## Revised Aug 2026

# This script loads PBTA splice junction counts and compares to normalized GTEx, Evo-devo and pediatric brain counts to identify tumor-enriched junctions 

# Load libraries
library(tidyverse)
library(qs2)

### Set directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "03-classify-tejs")
results_dir <- file.path(analysis_dir, "results")

source(file.path(analysis_dir, "util", "other-functions.R"))

## file paths

# PBTA
pbta_junction_file <- file.path(root_dir, "analyses",
                                "02-pbta-junction-processing",
                                "results",
                                "pbta-merged-norm-batch-corrected-junction-cts.qs2")

# Control
cohort_histologies_file <- file.path(root_dir, "analyses",
                                     "00-create-cohort-histologies",
                                     "results",
                                     "cohort-histologies.tsv")

control_sj_file <- file.path(data_dir,
                             "SJ.merged.control-cohort.tsv.gz")

# GTEx
gtex_junction_mat_file <- file.path(root_dir, "analyses",
                                    "01-ctrl-rmats-processing",
                                    "results",
                                    "gtex-merged-norm-junction-ct-mat.qs2")
gtex_junction_sd_file <- file.path(root_dir, "analyses",
                                   "01-ctrl-rmats-processing",
                                   "results",
                                   "gtex-merged-norm-junction-sd-mat.qs2")

# Evo-devo postnatal brain (collapsed to forebrain and hindbrain reference groups)
evodevo_junction_mat_file <- file.path(root_dir, "analyses",
                                       "01-ctrl-rmats-processing",
                                       "results",
                                       "evodevo-merged-postnatal-norm-junction-ct-mat.qs2")
evodevo_junction_sd_file <- file.path(root_dir, "analyses",
                                      "01-ctrl-rmats-processing",
                                      "results",
                                      "evodevo-merged-postnatal-norm-junction-sd-mat.qs2")
evodevo_metadata_file <- file.path(root_dir, "analyses",
                                   "01-ctrl-rmats-processing",
                                   "results",
                                   "evodevo-brain-prenatal-week-binned-metadata.tsv")

# pediatric normal brain
pedbrain_junction_mat_file <- file.path(root_dir, "analyses",
                                        "01-ctrl-rmats-processing",
                                        "results",
                                        "normal-pedbrain-merged-norm-junction-ct-mat.qs2")
pedbrain_junction_sd_file <- file.path(root_dir, "analyses",
                                       "01-ctrl-rmats-processing",
                                       "results",
                                       "normal-pedbrain-merged-norm-junction-sd-mat.qs2")

## Wrangle data

# Restrict this test run to ATRT samples.
atrt_sample_ids <- read_tsv(cohort_histologies_file,
                            show_col_types = FALSE) %>%
  dplyr::filter(plot_group == "Atypical Teratoid Rhabdoid Tumor") %>%
  dplyr::pull(Kids_First_Biospecimen_ID) %>%
  unique()

print("Loading PBTA junctions...")
pbta_junction_df <- qs2::qs_read(pbta_junction_file) %>%
  dplyr::filter(sample_id %in% atrt_sample_ids)
print(glue::glue("Retained {length(atrt_sample_ids)} ATRT samples and ",
                 "{nrow(pbta_junction_df)} PBTA junction rows."))

# Define the novel-junction prevalence cutoff for the cohort currently being
# surveyed. A junction must occur in fewer than 25% of retained samples.
n_samples_surveyed <- dplyr::n_distinct(pbta_junction_df$sample_id)
novel_junction_sample_cutoff <- 0.5 * n_samples_surveyed

# Load ctrl matrices 

# GTEx mean cpms
print("Loading control matrices...")
gtex_junction_mat <- qs2::qs_read(gtex_junction_mat_file) %>%
  # add "mean_" prefix to distinguish from cpm sd 
  rename_with(
    ~ paste0("mean_cpm_", .x),
    -junction
  )

# GTEx sd cpms
gtex_sd_mat <- qs2::qs_read(gtex_junction_sd_file) %>%
  # add "sd_ prefix
  rename_with(
    ~ paste0("sd_cpm_", .x),
    -junction
  )

evodevo_group_sizes <- readr::read_tsv(evodevo_metadata_file,
                                       show_col_types = FALSE) %>%
  dplyr::filter(!is.na(evodevo_postnatal_group)) %>%
  dplyr::transmute(
    evodevo_postnatal_group,
    region = sub("-.*", "", evodevo_postnatal_group)
  ) %>%
  dplyr::count(region, evodevo_postnatal_group, name = "n")

evodevo_postnatal_ref <- collapse_evodevo_postnatal_groups(
  qs2::qs_read(evodevo_junction_mat_file),
  qs2::qs_read(evodevo_junction_sd_file),
  evodevo_group_sizes
)
evodevo_junction_mat <- evodevo_postnatal_ref %>%
  dplyr::select(junction, starts_with("mean_cpm_"))
evodevo_sd_mat <- evodevo_postnatal_ref %>%
  dplyr::select(junction, starts_with("sd_cpm_"))

# Normal ped brain
pedbrain_junction_mat <- qs2::qs_read(pedbrain_junction_mat_file) %>%
  rename_with(
    ~ paste0("mean_cpm_", .x),
    -junction
  )

pedbrain_sd_mat <- qs2::qs_read(pedbrain_junction_sd_file) %>%
  rename_with(
    ~ paste0("sd_cpm_", .x),
    -junction
  )

# merge control mean cpm matrices
ctrl_junction_mat <- gtex_junction_mat %>%
  left_join(evodevo_junction_mat) %>%
  left_join(pedbrain_junction_mat) 

# Filter out junctions with mean CPM >= 10 in any control group.
keep_cols <- !grepl("junction", colnames(ctrl_junction_mat))

ctrl_junction_mat <- ctrl_junction_mat[
  rowSums(as.data.frame(ctrl_junction_mat)[, keep_cols, drop = FALSE] >= 10) == 0,
]

# merge control sd cpm matrices
ctrl_sd_mat <- gtex_sd_mat %>%
  left_join(evodevo_sd_mat) %>%
  left_join(pedbrain_sd_mat)

print("Filtering PBTA junctions...")
pbta_junction_sub_df <- pbta_junction_df %>%
  dplyr::filter(junction %in% ctrl_junction_mat$junction)

# tumor enriched calculations will be run on 10M rows at a time
n <- nrow(pbta_junction_sub_df)
chunk_size = 1e7

starts <- seq(1, n, by = chunk_size)

# create empty list to store chunk results
enr_jc_list <- list()

#loop through chunks to assess for junction enrichment
print("Classifying tumor-enriched junctions...")
for (i in 1:length(starts)) {
  
  # define current index
  start <- starts[i]
  
  idx <- start:min(start + chunk_size - 1, n)
  print(glue::glue("Processing chunk {i}..."))
  
  # filter full df 
  jc_chunk <- pbta_junction_sub_df[idx,]
  
  # join tumor and normal junction counts
  enr_jc_chunk <- jc_chunk %>%
    # define junction boundaries
    dplyr::mutate(boundary = sub(".*-(.*)-.*", "\\1", junction)) %>%
    dplyr::mutate(sample_boundary = glue::glue("{sample_id}-{boundary}")) %>%
    # join reference cohort mean and sd cpm matrices
    left_join(ctrl_junction_mat) %>%
    left_join(ctrl_sd_mat)
  
  # Calculate summary metrics without retaining a separate SNR column for
  # every control group.
  cpm_cols <- grep("^mean_cpm_", names(enr_jc_chunk), value = TRUE)
  n_chunk_rows <- nrow(enr_jc_chunk)
  max_control_cpm <- rep(-Inf, n_chunk_rows)
  min_control_snr <- rep(Inf, n_chunk_rows)

  for (j in seq_along(cpm_cols)) {
    cpm <- cpm_cols[j]
    sd <- sub("^mean_cpm_", "sd_cpm_", cpm)
    control_cpm <- enr_jc_chunk[[cpm]]
    snr <- (enr_jc_chunk$junction_cpm - control_cpm) / enr_jc_chunk[[sd]]

    max_control_cpm <- pmax(max_control_cpm, control_cpm, na.rm = TRUE)
    min_control_snr <- pmin(min_control_snr, snr, na.rm = TRUE)
  }
  
  # calculate summary stats & define tumor-enriched junctions
  enr_jc_chunk <- enr_jc_chunk %>%
    dplyr::mutate(
      # Minimum junction CPM fold-change across all control groups.
      max_mean_cpm_all = max_control_cpm,
      min_cpm_snr_all = min_control_snr,
      min_cpm_fc_all = junction_cpm/(max_control_cpm +  1e-5),
      ) %>%
    # determine if junction is tumor-enriched
    dplyr::mutate(junction_preference = case_when(
      # Junction enriched relative to all ref groups (FC & SNR > 5), ref groups do not have mean cpm > 10
      min_cpm_fc_all > 5 & min_cpm_snr_all > 5 & max_mean_cpm_all < 10 ~ "Tumor-enriched",
      TRUE ~ "Non-specific"
    ))
  
  # pull all boundaries called non-specific 
  nonspecific_boundaries <- enr_jc_chunk %>%
    dplyr::filter(junction_preference == "Non-specific") %>%
    pull(sample_boundary)
  
  # Retain only tumor-specific boundaries and the columns used downstream.
  enr_jc_list[[i]] <- enr_jc_chunk %>%
    dplyr::filter(junction_preference != "Non-specific" & !sample_boundary %in% nonspecific_boundaries) %>%
    dplyr::select(junction, sample_id, junction_count, junction_cpm, boundary,
                  junction_preference, min_cpm_fc_all, min_cpm_snr_all,
                  max_mean_cpm_all)

}

# merge into single df
merged_enr_jc_df <- bind_rows(enr_jc_list)

#### extraction of junctions not found in controls

# get junctions in control matrices
junctions_in_ctrl_df <- pbta_junction_df %>%
  distinct(junction) %>%
  dplyr::filter(junction %in% ctrl_sd_mat$junction) %>%
  dplyr::mutate(boundary = sub(".*-(.*)-.*", "\\1", junction))

# get junctions not in controls
ts_junctions <- pbta_junction_df %>%
  distinct(junction) %>%
  dplyr::filter(!junction %in% ctrl_sd_mat$junction) %>%
  # Exclude retained introns, represented as zero-length junctions whose two
  # splice boundaries are identical (fields 3 and 4).
  dplyr::mutate(
    junction_fields = strsplit(junction, ":|-|_"),
    up_boundary_coord = vapply(junction_fields, `[[`, character(1), 3),
    down_boundary_coord = vapply(junction_fields, `[[`, character(1), 4)
  ) %>%
  dplyr::filter(up_boundary_coord != down_boundary_coord) %>%
  dplyr::select(-junction_fields, -up_boundary_coord, -down_boundary_coord) %>%
  dplyr::mutate(boundary = sub(".*-(.*)-.*", "\\1", junction)) %>%
  dplyr::filter(!boundary %in% junctions_in_ctrl_df$boundary)

# Remove junctions that were absent from the rMATS control matrices but are
# observed in the independent control STAR SJ data.  Save the matched pairs
# for auditability before excluding them from the tumor-specific set.
control_sj_matches <- find_control_sj_matches(ts_junctions$junction,
                                              control_sj_file)

print(glue::glue(
  "Removing {dplyr::n_distinct(control_sj_matches$junction)} tumor-specific ",
  "junction(s) observed in the control SJ file."
))
ts_junctions <- ts_junctions %>%
  dplyr::filter(!junction %in% control_sj_matches$junction)

# get ts junction counts
ts_junction_ct_df <- pbta_junction_df %>%
  dplyr::filter(junction %in% ts_junctions$junction) %>%
  dplyr::group_by(junction) %>%
  dplyr::summarise(n_samples = dplyr::n_distinct(sample_id), .groups = "drop") %>%
  dplyr::arrange(desc(n_samples))

# import gtf, convert to df, filter for exons
gtf <- rtracklayer::import(file.path(data_dir,
                                     "gencode.v39.primary_assembly.annotation.gtf.gz"))

gtf_df <- as.data.frame(gtf) %>%
  dplyr::filter(type == "exon") %>%
  distinct(seqnames, start, end, .keep_all = TRUE)

# define exon starts and ends
starts <- unique(glue::glue("{gtf_df$seqnames}:{gtf_df$start}"))
ends <- unique(glue::glue("{gtf_df$seqnames}:{gtf_df$end}"))

# annotate junctions to exon boundaries
ts_junction_annot <- ts_junctions %>%
  # get unique junctions
  distinct(junction) %>%
  # split junction boundaries coords
  dplyr::mutate(chr = unlist(lapply(strsplit(junction, ":|-|_"), function(x) x[[1]])),
                up_bound = unlist(lapply(strsplit(junction, ":|-|_"), function(x) x[[3]])),
                down_bound = unlist(lapply(strsplit(junction, ":|-|_"), function(x) x[[4]]))) %>%
  dplyr::mutate(up_bound = glue::glue("{chr}:{up_bound}"),
                down_bound = glue::glue("{chr}:{down_bound}")) %>%
  # determine if junction bounds are annotated
  dplyr::mutate(up_bound_annotated = case_when(
    up_bound %in% ends ~ "Yes",
    TRUE ~ "No"
  )) %>%
  dplyr::mutate(down_bound_annotated = case_when(
    down_bound %in% starts ~ "Yes",
    TRUE ~ "No"
  ))

# Filter novel splice-site junctions to those present in <25% of the
# currently surveyed cohort.
ts_junction_ct_df <- ts_junction_ct_df %>%
  left_join(ts_junction_annot) %>%
  dplyr::mutate(novel_ss_usage = case_when(
    (up_bound_annotated == "No" | down_bound_annotated == "No") & up_bound != down_bound ~ "Yes",
    up_bound_annotated == "No" & down_bound_annotated == "No" ~ "Yes",
    TRUE ~ "No"
  )) %>%
  dplyr::filter(n_samples < novel_junction_sample_cutoff,
                novel_ss_usage == "Yes")

# get all ts junctions meeting criteria above
tesjs_na_ctrl <- pbta_junction_df %>%
  dplyr::filter(junction %in% ts_junction_ct_df$junction) %>%
  left_join(ts_junction_ct_df %>% dplyr::select(junction,
                                                novel_ss_usage)) %>%
  # Define the shared splice boundary here too, so it is retained after these
  # no-control junctions are appended to the control-compared calls.
  dplyr::mutate(boundary = sub(".*-(.*)-.*", "\\1", junction),
                criteria = case_when(
                  novel_ss_usage == "Yes" ~ "No ctrl expr, novel SS usage",
                  TRUE ~ "No ctrl expr"
                  ),
                junction_preference = "Tumor-enriched")

# append above junctions to merged enr jc df 
merged_enr_jc_df <- merged_enr_jc_df %>%
  dplyr::mutate(criteria = "enriched expr vs. ctrls") %>%
  bind_rows(tesjs_na_ctrl)
 
# Save the intermediate tumor-enriched calls for the next pipeline stage.
qs2::qs_save(merged_enr_jc_df,
             file.path(results_dir, "tumor-enriched-junctions-atrt.qs2"))

# print session info
sessionInfo()
