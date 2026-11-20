## Brain cell type splice event junction & target count matrix generation 
##
## Ryan Corbett
##
## Jan 2026

# This script performs the following:
# - loads brain cell type rMATS files 
# - pulls splice event junction and target coordinates and normalizes against rMATS input reads
# - calculates mean junction/target normalized counts per cell type

library(tidyverse)
library(qs2)
library(stringr)
library(data.table)

### Set directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "05-summarize-TEJs")
results_dir <- file.path(analysis_dir, "results")

source(file.path(analysis_dir,
                 "util", "rmats-processing-functions.R"))

# Set file paths

# GTEx
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

# Evo-devo
evodevo_se_file <- file.path(data_dir,
                             "evodevo-rmats_merged_raw_SE.qs2")
evodevo_ri_file <- file.path(data_dir,
                             "evodevo-rmats_merged_raw_RI.qs2")
evodevo_a3ss_file <- file.path(data_dir,
                               "evodevo-rmats_merged_raw_A3SS.qs2")
evodevo_a5ss_file <- file.path(data_dir,
                               "evodevo-rmats_merged_raw_A5SS.qs2")
evodevo_read_file <- file.path(data_dir,
                               "evodevo_input_read_counts.tsv")

# Pediatric brain 
pedbrain_se_file <- file.path(data_dir,
                              "normal_ped_brain-rmats_merged_raw_SE.qs2")
pedbrain_ri_file <- file.path(data_dir,
                              "normal_ped_brain-rmats_merged_raw_RI.qs2")
pedbrain_a3ss_file <- file.path(data_dir,
                                "normal_ped_brain-rmats_merged_raw_A3SS.qs2")
pedbrain_a5ss_file <- file.path(data_dir,
                                "normal_ped_brain-rmats_merged_raw_A5SS.qs2")
pedbrain_read_file <- file.path(data_dir,
                                "normal_ped_brain_input_read_counts.tsv")

# cell type
celltype_se_file <- file.path(data_dir,
                              "brain_cell_type-rmats_merged_raw_SE.qs2")
celltype_ri_file <- file.path(data_dir,
                              "brain_cell_type-rmats_merged_raw_RI.qs2")
celltype_a3ss_file <- file.path(data_dir,
                                "brain_cell_type-rmats_merged_raw_A3SS.qs2")
celltype_a5ss_file <- file.path(data_dir,
                                "brain_cell_type-rmats_merged_raw_A5SS.qs2")
celltype_read_file <- file.path(data_dir,
                                "brain_cell_type_input_read_counts.tsv")

ctrl_hist_file <- file.path(root_dir, "analyses",
                            "00-create-cohort-histologies",
                            "results",
                            "control-cohort-histologies.tsv")

recur_tej_file <- file.path(results_dir,
                             "recurrent-primary-tumor-enriched-oncofetal-splice-junctions.tsv.gz")

# Wrangle data
recur_tej_df <- read_tsv(recur_tej_file)

# Load read counts
gtex_read_cts <- read_tsv(gtex_read_file) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"))


# Load control cohort histologies (generated in 00-create-cohort-histologies)
ctrl_hist <- read_tsv(ctrl_hist_file)

evodevo_read_cts <- read_tsv(evodevo_read_file) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"))

pedbrain_read_cts <- read_tsv(pedbrain_read_file) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"))

celltype_read_cts <- read_tsv(celltype_read_file) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"))

gtex_list <- list("A3SS" = gtex_a3ss_file,
                  "A5SS" = gtex_a5ss_file,
                  "SE" = gtex_se_file,
                  "RI" = gtex_ri_file)

evodevo_list <- list("A3SS" = evodevo_a3ss_file,
                    "A5SS" = evodevo_a5ss_file,
                    "SE" = evodevo_se_file,
                    "RI" = evodevo_ri_file)

pedbrain_list <- list("A3SS" = pedbrain_a3ss_file,
                      "A5SS" = pedbrain_a5ss_file,
                      "SE" = pedbrain_se_file,
                      "RI" = pedbrain_ri_file)

celltype_list <- list("A3SS" = celltype_a3ss_file,
                      "A5SS" = celltype_a5ss_file,
                      "SE" = celltype_se_file,
                      "RI" = celltype_ri_file)

cohort_list <- list("GTEx" = gtex_list,
                    "Evo-devo" = evodevo_list,
                    "Pedbrain" = pedbrain_list,
                    "Celltype" = celltype_list)

read_cts <- gtex_read_cts %>%
  bind_rows(evodevo_read_cts,
            pedbrain_read_cts,
            celltype_read_cts)

# map cohort list names to histologies `cohort` values
hist_cohort_names <- c("GTEx" = "GTEx",
                       "Evo-devo" = "Evo-devo",
                       "Pedbrain" = "Pediatric brain",
                       "Celltype" = "Pediatric brain cell type")

cohort_jc_list <- list()

for (cohort in names(cohort_list)){
  
  print(glue::glue("Processing {cohort}..."))
  
  # get cohort rmats file list
  file_list <- cohort_list[[cohort]]
  
  # define cohort histologies
  hist <- ctrl_hist %>%
    dplyr::filter(cohort == hist_cohort_names[[.env$cohort]])
  
  # create empty list to store junction counts
  jc_list <- list()
  
  # loop through event types
  for (event in names(file_list)){
    
    print(glue::glue("loading {cohort} {event} rmats..."))
    rmats <- qs2::qs_read(file_list[[event]]) %>%
      dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$")) %>%
      dplyr::filter(geneSymbol %in% recur_tej_df$gene_symbol,
                    sample_id %in% hist$sample_id)
    
    n <- nrow(rmats)
    chunk_size = 1e7
    
    starts <- seq(1, n, by = chunk_size)
    
    chunk_df_list <-list()
    
    print(glue::glue("processing {cohort} {event} rmats chunks..."))
    for (i in 1:length(starts)) {
      
      print(glue::glue("processing chunk {i}..."))
      
      start <- starts[i]
      
      idx <- start:min(start + chunk_size - 1, n)
      
      rmats_chunk <- rmats[idx,]
      
      chunk_df_list[[i]] <- create_junction_df(rmats_chunk,
                                               event_type = event) %>%
        dplyr::filter(junction %in% recur_tej_df$junction)
      
    }
    
    # merge chunks
    jc_list[[event]] <- bind_rows(chunk_df_list)
    
  }
  
  # merge jc counts across events
  cohort_jc_list[[cohort]] <- bind_rows(jc_list)
  
}

print("Merging cohorts...")
merged_ctrl_tej_df <- cohort_jc_list[["GTEx"]] %>%
  bind_rows(cohort_jc_list[["Evo-devo"]],
            cohort_jc_list[["Pedbrain"]],
            cohort_jc_list[["Celltype"]])

# Generate tej cpm matrix
print("Generating TEJ cpm matrix...")
ctrl_tej_cpm_mat <- generate_norm_junction_mat(merged_ctrl_tej_df,
                                           read_cts,
                                           "sample_id")

# write to output
saveRDS(ctrl_tej_cpm_mat,
        file.path(results_dir,
                  "tumor-enriched-oncofetal-splice-junction-cpm-ctrls.rds"))


# print session info
sessionInfo()

