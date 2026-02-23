## Extract PBTA splice events associated with TESJs
##
## Ryan Corbett
##
## Feb 2026

# This script performs the following:
# - Loads PBTA rMATS files and generates splice event IDs
# - Filter splice events for those utilizing any TESJ
# - Saves to output

library(tidyverse)
library(qs2)
library(stringr)
library(data.table)

### Set directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "03-tesj-differential-splicing")
results_dir <- file.path(analysis_dir, "results")

source(file.path(analysis_dir, "util", "rmats-processing-functions.R"))

# Set file paths
se_file <- file.path(data_dir,
                     "pbta-rmats_merged_raw_SE.qs2")
ri_file <- file.path(data_dir,
                     "pbta-rmats_merged_raw_RI.qs2")
a3ss_file <- file.path(data_dir,
                       "pbta-rmats_merged_raw_A3SS.qs2")
a5ss_file <- file.path(data_dir,
                       "pbta-rmats_merged_raw_A5SS.qs2")

enr_jc_file <- file.path(root_dir, "analyses",
                         "02-tumor-enriched-splicing",
                         "results",
                         "tumor-enriched-oncofetal-splice-junctions-domain-expr-annotated.tsv.gz")

# Wrangle data
enr_jc_df <- read_tsv(enr_jc_file)

# define list of rmats file
file_list <- list("A3SS" = a3ss_file,
                  "A5SS" = a5ss_file,
                  "SE" = se_file,
                  "RI" = ri_file)

# create empty list to store splice events
event_list <- list()

# loop through rmats files
for (event in names(file_list)){

  print(glue::glue("loading {event} rmats..."))
  # Load rMATS
  rmats <- qs2::qs_read(file_list[[event]]) %>%
    # filter for gene symbols in TESJ df, rm events with NA PSI
    dplyr::filter(geneSymbol %in% enr_jc_df$gene_symbol,
                  !is.na(IncLevel1))
  
  # Run splice event extraction in chunks 
  # define chunk sizes and number of chunks
  n <- nrow(rmats)
  chunk_size = 1e7
  
  starts <- seq(1, n, by = chunk_size)
  
  # create empty list to store chunk splice events
  chunk_df_list <-list()
  
  # loop through chunks
  print(glue::glue("processing {event} rmats chunks..."))
  for (i in 1:length(starts)) {
    
    print(glue::glue("processing chunk {i}..."))
    
    # define interval for subsetting rmats
    start <- starts[i]
    
    idx <- start:min(start + chunk_size - 1, n)
    
    rmats_chunk <- rmats[idx,]
  
    # define junction and splice event coordinates from chunk
    rmats_chunk_df <- define_junctions_targets(rmats_chunk,
                                        event_type = event)
  
    # Create data frame of all junctions + associated splice events
    chunk_df_list[[i]] <- create_junction_df(rmats_chunk_df,
                                           event_type = event) %>%
      # Filter for splice events associated with TESJs
      dplyr::filter(glue::glue("{sample_id}:{junction}") %in% glue::glue("{enr_jc_df$sample_id}:{enr_jc_df$junction}")) %>%
      dplyr::rename(psi = IncLevel1) %>%
      dplyr::mutate(splicing_case = event)
    
  }
    # merge chunk results
    event_list[[event]] <- bind_rows(chunk_df_list)
  
}

# merge results across splice events
print("Merging splice events...")
merged_enr_jc_event_df <- event_list[["SE"]] %>%
  bind_rows(event_list[["RI"]],
            event_list[["A3SS"]],
            event_list[["A5SS"]])

# write to output
qs2::qs_save(merged_enr_jc_event_df,
             file.path(results_dir,
                       "tumor-enriched-oncofetal-junction-splice-events.qs2"))

# print session info
sessionInfo()

