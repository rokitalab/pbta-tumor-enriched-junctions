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
analysis_dir <- file.path(root_dir, "analyses", "02-tumor-enriched-splicing")
results_dir <- file.path(analysis_dir, "results")

## file paths

# PBTA
pbta_junction_file <- file.path(root_dir, "analyses",
                                "02-tumor-enriched-splicing",
                                "results",
                                "pbta-merged-norm-junction-cts.qs2")

# GTEx
gtex_junction_mat_file <- file.path(root_dir, "analyses",
                                    "01-ctrl-rmats-processing",
                                    "results",
                                    "gtex-merged-norm-junction-ct-mat.qs2")
gtex_junction_sd_file <- file.path(root_dir, "analyses",
                                   "01-ctrl-rmats-processing",
                                   "results",
                                   "gtex-merged-norm-junction-sd-mat.qs2")

# Evo-devo postnatal brain (region/stage-specific groups)
evodevo_junction_mat_file <- file.path(root_dir, "analyses",
                                       "01-ctrl-rmats-processing",
                                       "results",
                                       "evodevo-merged-postnatal-norm-junction-ct-mat.qs2")
evodevo_junction_sd_file <- file.path(root_dir, "analyses",
                                      "01-ctrl-rmats-processing",
                                      "results",
                                      "evodevo-merged-postnatal-norm-junction-sd-mat.qs2")

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
print("Loading PBTA junctions...")
pbta_junction_df <- qs2::qs_read(pbta_junction_file)

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

# Evo-devo
evodevo_junction_mat <- qs2::qs_read(evodevo_junction_mat_file) %>%
  rename_with(
    ~ paste0("mean_cpm_", .x),
    -junction
  )

evodevo_sd_mat <- qs2::qs_read(evodevo_junction_sd_file) %>%
  rename_with(
    ~ paste0("sd_cpm_", .x),
    -junction
  )

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
  rowSums(ctrl_junction_mat[, keep_cols, drop = FALSE] >= 10) == 0,
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
  dplyr::mutate(boundary = sub(".*-(.*)-.*", "\\1", junction)) %>%
  dplyr::filter(!boundary %in% junctions_in_ctrl_df$boundary)

# get ts junction counts
ts_junction_ct_df <- pbta_junction_df %>%
  dplyr::filter(junction %in% ts_junctions$junction) %>%
  dplyr::count(junction) %>% 
  dplyr::arrange(desc(n))

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

# filter ts junctions for those using novel SSs and in <10% of cohort
ts_junction_ct_df <- ts_junction_ct_df %>%
  left_join(ts_junction_annot) %>%
  dplyr::mutate(novel_ss_usage = case_when(
    (up_bound_annotated == "No" | down_bound_annotated == "No") & up_bound != down_bound ~ "Yes",
    up_bound_annotated == "No" & down_bound_annotated == "No" ~ "Yes",
    TRUE ~ "No"
  )) %>%
  # filter for junctions in <10% of cohort
  dplyr::filter(n < 250 & novel_ss_usage == "Yes") 

# get all ts junctions meeting criteria above
tesjs_na_ctrl <- pbta_junction_df %>%
  dplyr::filter(junction %in% ts_junction_ct_df$junction) %>%
  dplyr::mutate(criteria = "No ctrl expr, novel SS usage",
                junction_preference = "Tumor-enriched")

# append above junctions to merged enr jc df 
merged_enr_jc_df <- merged_enr_jc_df %>%
  dplyr::mutate(criteria = "enriched expr vs. ctrls") %>%
  bind_rows(tesjs_na_ctrl)
 
# Save the intermediate tumor-enriched calls for the next pipeline stage.
qs2::qs_save(merged_enr_jc_df,
             file.path(results_dir, "tumor-enriched-junctions.qs2"))

# print session info
sessionInfo()
