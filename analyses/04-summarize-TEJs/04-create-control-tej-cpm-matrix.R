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
analysis_dir <- file.path(root_dir, "analyses", "04-summarize-TEJs")
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

hist_file <- file.path(data_dir,
                       "histologies.tsv")

gtex_meta_file <- file.path(root_dir, "analyses",
                            "01-ctrl-rmats-processing", 
                            "input",
                            "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
gtex_v10_rm_file <- file.path(root_dir, "analyses",
                              "01-ctrl-rmats-processing", 
                              "input",
                              "gtex-v10-removed-samples.tsv")

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

evodevo_hist_file <- file.path(data_dir,
                       "evodevo-histologies.tsv")

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

pedbrain_hist_file <- file.path(data_dir,
                       "ped-normal-brain-histologies.tsv")

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

celltype_hist_file <- file.path(data_dir,
                                "GSE73721-normal-histologies.tsv")

recur_tej_file <- file.path(results_dir,
                             "recurrent-primary-tumor-enriched-oncofetal-splice-junctions.tsv.gz")

# Wrangle data
recur_tej_df <- read_tsv(recur_tej_file)

# Load read counts
gtex_read_cts <- read_tsv(gtex_read_file) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"))


### Prepare cohort histologies files 

# GTEx
gtex_meta_df <- read_tsv(gtex_meta_file)

pts_under40 <- gtex_meta_df %>%
  dplyr::filter(AGE %in% c("20-29",
                           "30-39")) %>%
  pull(SUBJID)

# Load file containing samples removed in v10, and get IDs
gtex_v10_rm_samples <- read_tsv(gtex_v10_rm_file) %>%
  pull(SAMPID)

# Load OPC hist and filter for GTEx RNA-seq
gtex_brain_under40_hist <- read_tsv(hist_file) %>%
  dplyr::filter(cohort == "GTEx",
                experimental_strategy == "RNA-Seq",
                grepl("Brain", gtex_subgroup)) %>%
  dplyr::mutate(id = unlist(lapply(strsplit(Kids_First_Biospecimen_ID, "-"), function(x) x[2]))) %>%
  dplyr::mutate(id = glue::glue("GTEX-{id}")) %>%
  # filter for pts under 40yo, filter out samples removed in v10
  dplyr::filter(id %in% pts_under40,
                !Kids_First_Biospecimen_ID %in% gtex_v10_rm_samples) %>%
  # only need BS ID, subgroup columns
  dplyr::select(Kids_First_Biospecimen_ID, gtex_subgroup) %>%
  dplyr::rename(sample_id = Kids_First_Biospecimen_ID,
                subgroup = gtex_subgroup) %>%
  dplyr::mutate(cohort = "GTEx")

# Evo-devo

evodevo_hist <- read_tsv(evodevo_hist_file) %>%
  # filter for brain samples, filter out older life stages
  dplyr::filter(primary_site %in% c("Forebrain", "Hindbrain"),
                !pathology_free_text_diagnosis %in% c("Middle Adult",
                                                      "Elderly")) %>%
  # define subgroups (region + stage)
  dplyr::mutate(evodevo_subgroup = glue::glue("{primary_site}-{pathology_free_text_diagnosis}")) %>%
  # only need BS ID, subgroup columns
  dplyr::select(Kids_First_Biospecimen_ID, evodevo_subgroup) %>%
  dplyr::rename(sample_id = Kids_First_Biospecimen_ID,
                subgroup = evodevo_subgroup) %>%
  dplyr::mutate(cohort = "Evo-devo")

evodevo_read_cts <- read_tsv(evodevo_read_file) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"))

# ped brain

pedbrain_hist <- read_tsv(pedbrain_hist_file) %>%
  # filter out tumor-infiltrated pons
  dplyr::filter(sample_id != "7316-7585") %>%
  # only need sample ID, primary_site columns
  dplyr::select(sample_id, primary_site) %>%
  dplyr::rename(subgroup = primary_site) %>%
  dplyr::mutate(cohort = "Pediatric brain")

pedbrain_read_cts <- read_tsv(pedbrain_read_file) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"))

# cell type 

celltype_hist <- read_csv(celltype_hist_file) %>%
  # only need run ID, cell type columns
  distinct(Run, cell_type) %>%
  dplyr::rename(sample_id = Run,
                group = cell_type) %>%
  dplyr::mutate(cohort = "Pediatric brain cell type")

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

cohort_jc_list <- list()

for (cohort in names(cohort_list)){
  
  print(glue::glue("Processing {cohort}..."))
  
  # get cohort rmats file list
  file_list <- cohort_list[[cohort]]
  
  # define cohort histologies
  if (cohort == "GTEx"){
    
    hist <- gtex_brain_under40_hist
    
  } else if (cohort == "Evo-devo"){
    
    hist <- evodevo_hist
    
  } else if (cohort == "Pedbrain"){
    
    hist <- pedbrain_hist
    
  } else {
    
    hist <- celltype_hist
    
  }
  
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


# merge histologies 
ctrl_hist <- gtex_brain_under40_hist %>%
  bind_rows(evodevo_hist, pedbrain_hist, celltype_hist) %>%
  dplyr::select(sample_id, cohort, everything())

# write to output
write_tsv(ctrl_hist,
          file.path(results_dir, "control-cohort-histologies.tsv"))

# print session info
sessionInfo()

