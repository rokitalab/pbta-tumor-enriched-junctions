## Filter tumor-specific splice variants
##
## Ryan Corbett
##
## 2025

# This script filters tumor-specific splice variants as follows:
# 
# 1)
# 2)

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
                         "tumor-enriched-oncofetal-splice-junctions.tsv.gz")

cds_file <- file.path(results_dir, 
                      "tumor-enriched-oncofetal-splice-junctions-cds.bed")

pfam_file <- file.path(results_dir,
                       "tumor-enriched-oncofetal-splice-junctions-pfam-annotated.tsv.gz")

domain_file <- file.path(results_dir,
                         "tumor-enriched-oncofetal-splice-junctions-domain-anno.uniq.tsv")

exp_file <- file.path(data_dir, 
                      "pbta_gene-expression-rsem-tpm-collapsed.rds")

# Assign parsed arguments to variables
tpm_cutoff <- 10

enr_jc_df <- read_tsv(enr_jc_file)

# Load cds, uniprot topological domain, and pfam domain annotation files
cds_df <- read_tsv(cds_file,
                   col_names = c("chr", "start", "end",
                                 "junction", "junction_preference",
                                 "strand", "type"))

# read in results file, get bs ids for that file only
uniprot_domain_df <- read_tsv(domain_file) %>%
  dplyr::mutate(coverage = as.double(coverage)) %>%
  filter(domain_type %in% c("EC", "TM", "IC")) %>%
  dplyr::filter(junction %in% enr_jc_df$junction) %>%
  dplyr::arrange(domain_type) %>%
  distinct(junction, .keep_all = TRUE)

pfam_domain_df <- read_tsv(pfam_file)

# append cds, uniprot, and pfam domain annotations to ts events df
enr_jc_domain_df <- enr_jc_df %>%
  # define if splice event exon is in cds
  dplyr::mutate(in_cds = case_when(
    junction %in% cds_df$junction ~ "Yes",
    TRUE ~ "No"
  )) %>%
  # Join uniprot annotation
  left_join(uniprot_domain_df %>% 
              dplyr::select(-overlap_start,
                            -overlap_end,
                            -junction_preference),
            by = "junction") %>%
  dplyr::rename("uniprot_domain" = domain_type) %>%
  # Join Pfam annotation
  left_join(pfam_domain_df %>% 
              dplyr::select(junction, pfam_id, pfam_name,
                            pfam_description, domain_start,
                            domain_end, junction_overlaps_domain)) %>%
  dplyr::rename(pfam_domain_start = domain_start,
                pfam_domain_end = domain_end,
                gene_symbol = geneSymbol)

### Expression filtering 

samples_of_int <- unique(enr_jc_domain_df$sample_id)

# read in exp file and filter for samples in ts events df
exp <- readRDS(exp_file) %>%
  select(any_of(samples_of_int))

# append TPM and only retain splice events in genes with TPM >= defined cutoff
enr_jc_domain_expr_df <- add_TPM_values(enr_jc_domain_df, exp) %>%
  filter(gene_tpm >= tpm_cutoff) 

# remove unnecessary columns, rename and reorder remaining columns
enr_jc_domain_expr_df <- enr_jc_domain_expr_df %>%
  dplyr::select(sample_id, junction, chr, strand, 
                gene_symbol, junction_cpm,
                junction_preference,
              #  criteria,
                gene_tpm,
                everything()) %>%
  # Define consequence based on uniprot and pfam domain annotations
  dplyr::mutate(consequence = case_when(
    uniprot_domain == "EC" & !is.na(pfam_id) ~ "EC domain alteration, Other functional",
    uniprot_domain == "EC" ~ "EC domain alteration",
    !is.na(pfam_id) ~ "Other functional",
    TRUE ~ "Unknown"
  ))

write_tsv(enr_jc_domain_expr_df,
          file.path(results_dir,
                    "tumor-enriched-oncofetal-splice-junctions-domain-expr-annotated.tsv.gz"))
