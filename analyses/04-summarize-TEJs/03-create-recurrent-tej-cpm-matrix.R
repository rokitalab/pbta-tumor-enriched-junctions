## Create recurrent TEJ expression matrix (PBTA)
##
## Ryan Corbett
##
## Jan 2026

# This script performs the following:
# - loads PBTA rMATS files 
# - pulls splice event junction coordinates, filters for TEJs, and normalizes against rMATS input reads
# - Creates matrix of CPM values for all PBTA samples

# Load libraries
library(tidyverse)
library(qs2)
library(stringr)
library(data.table)

### Set directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "04-summarize-TEJs")
results_dir <- file.path(analysis_dir, "results")

source(file.path(analysis_dir,
                 "util", "rmats-processing-functions.R"))

# Set file paths
se_file <- file.path(data_dir,
                     "pbta-rmats_merged_raw_SE.qs2")
ri_file <- file.path(data_dir,
                     "pbta-rmats_merged_raw_RI.qs2")
a3ss_file <- file.path(data_dir,
                       "pbta-rmats_merged_raw_A3SS.qs2")
a5ss_file <- file.path(data_dir,
                       "pbta-rmats_merged_raw_A5SS.qs2")
read_file <- file.path(data_dir,
                       "pbta_input_read_counts.tsv")

hist_file <- file.path(root_dir, "analyses",
                       "00-create-cohort-histologies",
                       "results",
                       "cohort-histologies.tsv")

recur_tej_file <- file.path(results_dir,
                         "recurrent-primary-tumor-enriched-oncofetal-splice-junctions.tsv.gz")

recur_tej_df <- read_tsv(recur_tej_file)

file_list <- list("A3SS" = a3ss_file,
                  "A5SS" = a5ss_file,
                  "SE" = se_file,
                  "RI" = ri_file)

# Load read counts, cohort histologies
read_cts <- read_tsv(read_file) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"))

hist <- read_tsv(hist_file)

# create empty list to store junctions
jc_list <- list()

# loop through rmats files
for (event in names(file_list)){
  
  print(glue::glue("loading {event} rmats..."))
  # filter rmats for samples in hist, genes in TEJ df
  rmats <- qs2::qs_read(file_list[[event]]) %>%
    dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$")) %>%
    dplyr::filter(geneSymbol %in% recur_tej_df$gene_symbol,
                  sample_id %in% hist$Kids_First_Biospecimen_ID)
  
  # define chunk sizes
  n <- nrow(rmats)
  chunk_size = 1e7
  
  starts <- seq(1, n, by = chunk_size)
  
  chunk_df_list <-list()
  
  # loop through chunks
  print(glue::glue("processing {event} rmats chunks..."))
  for (i in 1:length(starts)) {
    
    print(glue::glue("processing chunk {i}..."))
    
    start <- starts[i]
    
    # define index
    idx <- start:min(start + chunk_size - 1, n)
    
    rmats_chunk <- rmats[idx,]
    
    # extract junctions
    chunk_df_list[[i]] <- create_junction_df(rmats_chunk,
                                             event_type = event) %>%
      # filter for TEJ junctions
      dplyr::filter(junction %in% recur_tej_df$junction)
    
  }
  
  # merge junctions
  jc_list[[event]] <- bind_rows(chunk_df_list)
  
}

# merge junctions across splice event types
print("Merging splice events...")
merged_tej_jc_df <- jc_list[["SE"]] %>%
  bind_rows(jc_list[["RI"]],
            jc_list[["A3SS"]],
            jc_list[["A5SS"]])

# Generate tej cpm matrix
print("Generating tej cpm matrix...")
tej_cpm_mat <- generate_norm_junction_mat(merged_tej_jc_df,
                                           read_cts,
                                           "sample_id")

# write to output
saveRDS(tej_cpm_mat,
        file.path(results_dir,
                  "tumor-enriched-oncofetal-splice-junction-cpm.rds"))

# print session info
sessionInfo()
