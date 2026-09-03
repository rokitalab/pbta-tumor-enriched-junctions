## Merge tumor-enriched splice junction annotations and filter for high expression
##
## Ryan Corbett
##
## 2025

# This script performs the following
# 
# 1) Merges tumor-enriched junction coordinates, uniprot domain, and pfam domain annotations
# 2) Filters TEJs for those in genes and samples with TPM >= 10

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
  library(rprojroot)
  library(R.utils)
})

# create directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "02-tumor-enriched-splicing")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# source function
source(file.path(analysis_dir, "util", "add-tpm-values.R"))

# Set file paths
enr_jc_file <- file.path(results_dir, 
                         "tumor-enriched-oncofetal-splice-junctions-atrt.tsv.gz")

cds_file <- file.path(results_dir, 
                      "tumor-enriched-oncofetal-splice-junctions-cds-atrt.bed")

pfam_file <- file.path(results_dir,
                       "tumor-enriched-oncofetal-splice-junctions-pfam-annotated-atrt.tsv.gz")

domain_file <- file.path(results_dir,
                         "tumor-enriched-oncofetal-splice-junctions-domain-anno-atrt.uniq.tsv")

exp_file <- file.path(data_dir, 
                      "pbta_gene-expression-rsem-tpm-collapsed.rds")

# Define tpm_cutff
tpm_cutoff <- 10

### Wrangle data
enr_jc_df <- read_tsv(enr_jc_file)

# Load cds, uniprot topological domain, and pfam domain annotation files
cds_df <- read_tsv(cds_file,
                   col_names = c("chr", "start", "end",
                                 "junction", "junction_preference",
                                 "strand", "type"))

# format uniprot annotation df
uniprot_domain_df <- read_tsv(domain_file) %>%
  filter(domain_type %in% c("EC", "TM", "IC")) %>%
  dplyr::filter(junction %in% enr_jc_df$junction) %>%
  dplyr::arrange(domain_type) %>%
  distinct(junction, .keep_all = TRUE) %>%
  # rename columns
  dplyr::rename(uniprot_domain_start = domain_start,
                uniprot_domain_end = domain_end,
                uniprot_domain_type = domain_type) %>% 
  # deselect columns to not include in final output
  dplyr::select(-overlap_start,
                -overlap_end,
                -junction_preference)

# format pfam annotation df
pfam_domain_df <- read_tsv(pfam_file) %>%
  # rename columns
  dplyr::rename(pfam_domain_start = domain_start,
                pfam_domain_end = domain_end,
                pfam_domain_overlap_type = domain_overlap_status,
                junction_overlaps_pfam_domain = junction_overlaps_domain) %>%
  # select relevant columns to include in final output
  dplyr::select(junction, pfam_id, pfam_name,
                pfam_description, pfam_domain_start,
                pfam_domain_end, pfam_domain_overlap_type,
                junction_overlaps_pfam_domain)

### Merge annotations

# append cds, uniprot, and pfam domain annotations to tej df
enr_jc_domain_df <- enr_jc_df %>%
  # define if splice event exon is in cds
  dplyr::mutate(in_cds = case_when(
    junction %in% cds_df$junction ~ "Yes",
    TRUE ~ "No"
  )) %>%
  # Join uniprot annotation
  left_join(uniprot_domain_df,
            by = "junction") %>%
  # Join Pfam annotation
  left_join(pfam_domain_df) %>%
  dplyr::rename(gene_symbol = geneSymbol)

### Expression filtering 

samples_of_int <- unique(enr_jc_domain_df$sample_id)

# read in exp file and filter for samples in tej df
exp <- readRDS(exp_file) %>%
  select(any_of(samples_of_int))

# append TPM and only retain junctions in genes with TPM >= defined cutoff
enr_jc_domain_expr_df <- add_TPM_values(enr_jc_domain_df, exp) %>%
  filter(gene_tpm >= tpm_cutoff) 

# remove unnecessary columns, rename and reorder remaining columns
enr_jc_domain_expr_df <- enr_jc_domain_expr_df %>%
  dplyr::select(sample_id, junction, chr, strand, 
                gene_symbol, junction_cpm,
                junction_preference,
                gene_tpm,
                everything()) %>%
  # Define consequence based on uniprot and pfam domain annotations
  dplyr::mutate(consequence = case_when(
    uniprot_domain_type == "EC" & !is.na(pfam_id) ~ "EC domain alteration, Other functional",
    uniprot_domain_type == "EC" ~ "EC domain alteration",
    !is.na(pfam_id) ~ "Other functional",
    TRUE ~ "Unknown"
  ))

# write to output
write_tsv(enr_jc_domain_expr_df,
          file.path(results_dir,
                    "tumor-enriched-oncofetal-splice-junctions-domain-expr-annotated-atrt.tsv.gz"))

# Print session info
sessionInfo()
