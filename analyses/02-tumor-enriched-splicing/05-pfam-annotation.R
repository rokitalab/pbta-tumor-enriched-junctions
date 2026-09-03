## Pfam domain annotation
##
## Ryan Corbett
##
## 2025
#
# This scripts annotates tumor-specific junctions to Pfam protein domains


# Load libraries
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

enr_jc_events_file <- file.path(results_dir,
                                "tumor-enriched-oncofetal-splice-junctions-atrt.tsv.gz")

enr_jc_df <- read_tsv(enr_jc_events_file) %>%
  dplyr::rename(gene_symbol = geneSymbol) %>%
  distinct(junction, .keep_all = TRUE)

## Pfam annotation

# Load Pfam domain database from annoFuseData package
bioMartDataPfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuseData")) %>%
  dplyr::rename(gene_symbol = hgnc_symbol,
                pfam_name = NAME,
                pfam_description = DESC)

# extract relevant columns from merged res and left join with pfam df by gene symbol
enr_jc_pfam_df <- enr_jc_df %>%
  dplyr::select(junction, gene_symbol, 
                up_jc_start, up_jc_end,
                down_jc_start, down_jc_end) %>%
  distinct() %>%
  left_join(bioMartDataPfam %>%
              dplyr::select(gene_symbol, pfam_id,
                            pfam_name, pfam_description,
                            domain_start, domain_end), by = "gene_symbol")

# calculate junction overlap with domains
enr_jc_pfam_df <- enr_jc_pfam_df %>%
  # Determine if junction overlaps domain
  dplyr::mutate(domain_overlap_status = case_when(
    # domain completely within junction interval
    domain_start >= up_jc_end & domain_end <= down_jc_start ~ "Domain completely in junction interval",
    # junction interval completely within domain
    up_jc_end >= domain_start & down_jc_start <= domain_end ~ "Junction interval completely in domain",
    # partial overlap b/w junction interval and domain
    (up_jc_end >= domain_start & up_jc_end <= domain_end) | (down_jc_start >= domain_start & down_jc_start <= domain_end) ~ "Partial overlap",
    TRUE ~ "No overlap"
  )) %>%
  dplyr::mutate(junction_overlaps_domain = case_when(
    domain_overlap_status != "No overlap" ~ "Yes",
    TRUE ~ "No"
  )) %>%
  # only retain junctions that overlap domains
  dplyr::filter(junction_overlaps_domain == "Yes")

# add pfam annotation to merged ts events df
enr_jc_events_w_pfam <- enr_jc_df %>%
  left_join(enr_jc_pfam_df) %>%
  distinct(junction, .keep_all = TRUE)

# print number of junctions overlapping pfam domain
print("Junction overlaps Pfam domain:")
table(!is.na(enr_jc_events_w_pfam$junction_overlaps_domain))

# write to output
write_tsv(enr_jc_events_w_pfam,
          file.path(results_dir,
                    "tumor-enriched-oncofetal-splice-junctions-pfam-annotated-atrt.tsv.gz"))

# print session info
sessionInfo()
