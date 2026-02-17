## Identify tumor-enriched junctions in PBTA 
##
## Ryan Corbett
##
## Jan 2025

# This script loads PBTA splice junction counts and compares to normalized GTEx, Evo-devo and pediatric brain counts to identify tumor-enriched junctions 

# Load libraries
library(tidyverse)
library(qs2)

### Set directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "02-tumor-enriched-splicing")
input_dir <- file.path(analysis_dir, "input")
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

# Evo-devo
evodevo_junction_mat_file <- file.path(root_dir, "analyses",
                                       "01-ctrl-rmats-processing",
                                       "results",
                                       "evodevo-merged-broadgroup-norm-junction-ct-mat.qs2")
evodevo_junction_sd_file <- file.path(root_dir, "analyses",
                                      "01-ctrl-rmats-processing",
                                      "results",
                                      "evodevo-merged-broadgroup-norm-junction-sd-mat.qs2")

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

# filter out any junctions with mean cpm >= 10 in any postnatal cohort (will not be considered tumor-enriched) 
keep_cols <- !grepl("junction|fetal", colnames(ctrl_junction_mat))

ctrl_junction_mat <- ctrl_junction_mat[
  rowSums(ctrl_junction_mat[, keep_cols, drop = FALSE] >= 10) == 0,
]

# merge control sd cpm matrices
ctrl_sd_mat <- gtex_sd_mat %>%
  left_join(evodevo_sd_mat) %>%
  left_join(pedbrain_sd_mat)

# calculate number of control groups
n_ctrl_grps <- ncol(ctrl_junction_mat) - 1

n_fetal_grps <- sum(grepl("fetal", colnames(ctrl_junction_mat)))
n_postnatal_grps <- n_ctrl_grps - n_fetal_grps

print("Filtering PBTA junctions...")
pbta_junction_sub_df <- pbta_junction_df %>%
  dplyr::filter(junction %in% ctrl_junction_mat$junction)

# tumor enriched calculations will be run on 10M rows at a time
n <- nrow(pbta_junction_sub_df)
chunk_size = 1e7

starts <- seq(1, n, by = chunk_size)

# create empty list to store chunk results
enr_jc_list <- list()

# define cpm mean and sd column names
mean_cpm_cols <- grep("^mean_cpm_", colnames(ctrl_junction_mat), value = TRUE)
sd_cpm_cols  <- grep("^sd_cpm_", colnames(ctrl_sd_mat), value = TRUE)

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
    # join reference cohort mean and sd cpm matrices
    left_join(ctrl_junction_mat) %>%
    left_join(ctrl_sd_mat)
  
  # Get all "mean_cpm_" columns and corresponding "mean_sd_" columns
  cpm_cols <- grep("^mean_cpm_", names(enr_jc_chunk), value = TRUE)
  sd_cols  <- sub("^mean_cpm_", "sd_cpm_", cpm_cols)
  
  # Loop over each column pair and compute signal-to-noise ratio
  for (j in seq_along(cpm_cols)) {
    
    cpm <- cpm_cols[j]
    sd  <- sd_cols[j]
    
    # define column name
    snr_name <- paste0("cpm_SNR_", sub("^mean_cpm_", "", cpm))
    
    # calculate SNR as (tumor junction cpm - control cpm mean)/(control cpm sd)
    enr_jc_chunk[[snr_name]] <- (enr_jc_chunk$junction_cpm - enr_jc_chunk[[cpm]]) / enr_jc_chunk[[sd]]
    
  }
  
  # calculate summary stats & define tumor-enriched junctions
  enr_jc_chunk <- enr_jc_chunk %>%
    dplyr::mutate(
      # maximum junction mean CPM in fetal and postnatal reference cohorts
      max_mean_cpm_all = do.call(pmax, c(select(., contains("mean_cpm")), na.rm = TRUE)),
      max_mean_cpm_postnatal = do.call(pmax, c(select(., contains("mean_cpm") & !contains("fetal")), na.rm = TRUE)),
      # minimum junction CPM fold-change across all reference cohorts and postnatal cohorts
      min_cpm_fc_all = junction_cpm/(max_mean_cpm_all +  1e-5),
      min_cpm_fc_postnatal = junction_cpm/(max_mean_cpm_postnatal + 1e-5),
      # minimum junction CPM signal-to-noise ratio across all reference cohorts and postnatal cohorts
      min_cpm_snr_all = do.call(pmin, c(select(., contains("cpm_SNR")), na.rm = TRUE)),
      min_cpm_snr_postnatal = do.call(pmin, c(select(., contains("cpm_SNR") & !contains("fetal")), na.rm = TRUE)),
      # Number of ref groups with NA junction counts
      # n_na_ctrl_all = rowSums(across(contains("mean_cpm"), ~ is.na(.), .names = NULL), na.rm = TRUE),
      # n_na_ctrl_postnatal = rowSums(across(!contains("fetal") & contains("mean_cpm"), ~ is.na(.), .names = NULL), na.rm = TRUE)
      ) %>%
    # determine if junction is tumor-enriched
    dplyr::mutate(junction_preference = case_when(
      # Junction enriched relative to all ref groups (FC & SNR > 5), ref groups do not have mean cpm > 10
      min_cpm_fc_all > 5 & min_cpm_snr_all > 5 & max_mean_cpm_all < 10 ~ "Tumor-enriched",
      # same criteria as above but for postnatal ref groups only
      min_cpm_fc_postnatal > 5 & min_cpm_snr_postnatal > 5 & max_mean_cpm_postnatal < 10 ~ "Oncofetal",
      TRUE ~ "Non-specific"
    ))
  
  # pull all boundaries called non-specific 
  nonspecific_boundaries <- enr_jc_chunk %>%
    dplyr::filter(junction_preference == "Non-specific") %>%
    dplyr::mutate(sample_boundary = glue::glue("{sample_id}-{boundary}")) %>%
    pull(sample_boundary)
  
  # retain only tumor-specific boundaries
  enr_jc_list[[i]] <- enr_jc_chunk %>%
    dplyr::mutate(sample_boundary = glue::glue("{sample_id}-{boundary}")) %>%
    dplyr::filter(junction_preference != "Non-specific" & !sample_boundary %in% nonspecific_boundaries)

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
  # filter for novel ss usage junctions in <10% of cohort
  dplyr::filter(novel_ss_usage == "Yes",
                n < 250) 

# get all ts junctions meeting criteria above
tesjs_na_ctrl <- pbta_junction_df %>%
  dplyr::filter(junction %in% ts_junction_ct_df$junction) %>%
  dplyr::mutate(criteria = "No ctrl expr, novel SS usage")

# append above junctions to merged enr jc df 
merged_enr_jc_df <- merged_enr_jc_df %>%
  dplyr::mutate(criteria = "enriched expr vs. ctrls") %>%
  bind_rows(tesjs_na_ctrl)
 
# Load junction annotation file
jc_annot <- read_tsv(file.path(results_dir, "junction-annot.tsv.gz")) %>% 
  dplyr::filter(junction %in% merged_enr_jc_df$junction) %>%
  distinct(junction, geneSymbol, event_type, .keep_all = TRUE) %>%
  group_by(junction, strand, chr, 
           up_jc_start, up_jc_end, down_jc_start,
           down_jc_end) %>%
  summarise(geneSymbol = str_c(geneSymbol, collapse = ",")) %>%
  ungroup()

# merge annotation to tumor-enriched junctions
merged_enr_jc_annot_df <- merged_enr_jc_df %>%
  dplyr::select(junction, sample_id, 
               junction_count, junction_cpm,
               boundary,
               junction_preference,
               min_cpm_fc_all, min_cpm_fc_postnatal,
               min_cpm_snr_all, min_cpm_snr_postnatal,
               max_mean_cpm_all, max_mean_cpm_postnatal) %>%
  left_join(jc_annot)

# save to output
write_tsv(merged_enr_jc_annot_df,
          file.path(results_dir, "tumor-enriched-oncofetal-splice-junctions.tsv.gz"))

# create bed files for uniprot annotation

# upstream exon 
jc_bed_df <- merged_enr_jc_annot_df %>%
  dplyr::mutate(
    up_jc_end     = as.integer(up_jc_end),
    down_jc_start = as.integer(down_jc_start)
  ) %>%
  dplyr::select(chr, up_jc_end, down_jc_start,
                junction, junction_preference, 
                strand) %>%
  dplyr::rename(start = up_jc_end,
                end = down_jc_start,
                preference = junction_preference,
                id = junction) %>%
  distinct()

# write to output
write.table(jc_bed_df,
            file.path(results_dir,
                      "tumor-enriched-oncofetal-splice-junctions.bed"),
            col.names = FALSE, row.names = FALSE,
            quote = FALSE, sep = "\t")

# print session info
sessionInfo()

