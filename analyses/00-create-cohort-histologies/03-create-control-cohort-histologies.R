# Create control cohort histologies file
#
# Ryan Corbett
#
# Jan 2026

# This script performs the following:
# - filters GTEx brain RNA-seq samples to donors <40yo, excluding samples removed in GTEx v10
# - filters Evo-devo samples to forebrain/hindbrain samples, excluding middle adult and elderly life stages
# - filters pediatric normal brain samples, excluding tumor-infiltrated pons and cerebellum samples
# - loads brain cell type (GSE73721) samples
# - merges all control cohorts into a single histologies file

# load packages
library(tidyverse)

### Set directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "00-create-cohort-histologies")
results_dir <- file.path(analysis_dir, "results")

# Set file paths
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

evodevo_hist_file <- file.path(data_dir,
                               "evodevo-histologies.tsv")

pedbrain_hist_file <- file.path(data_dir,
                                "ped-normal-brain-histologies.tsv")

celltype_hist_file <- file.path(data_dir,
                                "GSE73721-normal-histologies.tsv")

### Prepare cohort histologies

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
  dplyr::select(Kids_First_Biospecimen_ID, gtex_subgroup,
                RNA_library) %>%
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
  dplyr::select(Kids_First_Biospecimen_ID, evodevo_subgroup,
                RNA_library) %>%
  dplyr::rename(sample_id = Kids_First_Biospecimen_ID,
                subgroup = evodevo_subgroup) %>%
  dplyr::mutate(cohort = "Evo-devo")

# ped brain
pedbrain_hist <- read_tsv(pedbrain_hist_file) %>%
  # filter out tumor-infiltrated samples (pons, cerebellum)
  dplyr::filter(!sample_id %in% c("7316-7585", "SRR26129084")) %>%
  # only need sample ID, primary_site columns
  dplyr::select(sample_id, primary_site, RNA_library) %>%
  dplyr::rename(subgroup = primary_site) %>%
  dplyr::mutate(cohort = "Pediatric brain")

# cell type
celltype_hist <- read_csv(celltype_hist_file) %>%
  # only need run ID, cell type columns
  distinct(Run, cell_type) %>%
  dplyr::rename(sample_id = Run,
                subgroup = cell_type) %>%
  dplyr::mutate(cohort = "Pediatric brain cell type") %>%
  dplyr::mutate(RNA_library = "stranded")

# merge histologies
ctrl_hist <- gtex_brain_under40_hist %>%
  bind_rows(evodevo_hist, pedbrain_hist, celltype_hist) %>%
  dplyr::select(sample_id, cohort, everything())

# write to output
write_tsv(ctrl_hist,
          file.path(results_dir, "control-cohort-histologies.tsv"))

# print session info
sessionInfo()
