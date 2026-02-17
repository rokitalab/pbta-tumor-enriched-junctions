## Extract PBTA normalized junction counts
##
## Ryan Corbett
##
## Jan 2026

# This script performs the following:
# - loads PBTA rMATS files 
# - pulls splice event junction coordinates and counts
# - normalizes against rMATS input reads
# - saves merged, normalized junction expression file

# Load libraries 
library(tidyverse)
library(qs2)
library(stringr)
library(data.table)

### Set directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "02-tumor-enriched-splicing")
input_dir <- file.path(analysis_dir, "input")
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
read_file <- file.path(data_dir,
                       "pbta_input_read_counts.tsv")

samples_to_rm_file <- file.path(input_dir,
                                "pbta-rna-high-intron-samples.tsv")

# Load samples to remove (high intronic read count)
samples_to_rm <- read_tsv(samples_to_rm_file) %>%
  pull(Kids_First_Biospecimen_ID)

### Load rMATS files and prepare junction dfs

# Load rMATS SE results
print("Loading SE rMATS...")

se_df <- qs2::qs_read(se_file) %>%
  # retain relevant columns
  dplyr::select(
    sample_id, geneSymbol, strand,
    chr,
    upstreamES, upstreamEE,
    exonStart_0base, exonEnd,
    downstreamES, downstreamEE,
    upstream_to_target_count,
    target_to_downstream_count,
    upstream_to_downstream_count
  )
  # remove suffix from sample ID, convert coordinates to one-based
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"),
                exonStart_0base = exonStart_0base + 1,
                upstreamES = upstreamES + 1,
                downstreamES = downstreamES + 1) %>%
  # rm low-quality RNA samples, events with low junction counts
  dplyr::filter(!sample_id %in% samples_to_rm,
                (upstream_to_target_count >= 10 | target_to_downstream_count >= 10 | upstream_to_downstream_count >= 10))

# Define SE junction and target IDs, and select relevant columns
print("Extracting SE event junctions...")

# do this in chunks to reduce memory
n <- nrow(se_df)
chunk_size = 1e7

starts <- seq(1, n, by = chunk_size)

# create empty list to store junction coordinates & counts
se_jc_list <- list()

# loop through chunks
  for (i in 1:length(starts)) {

    print(glue::glue("processing chunk {i}..."))

    start <- starts[i]

    # get current index
    idx <- start:min(start + chunk_size - 1, n)

    # subset full df
    se_chunk <- se_df[idx,]

    # Define junction coordinates
    se_jc_list[[i]] <- define_junctions_targets(se_chunk,
                                         event_type = "SE")

  }

# merge results from chunks
se_df <- bind_rows(se_jc_list)

# Build the long-form junction table, merging data from:
# upstream-exon
# exon-downstream
# upstream-downstream junctions
print("Creating SE junction df...")

se_junction_df <- create_junction_df(se_df,
                                     event_type = "SE") %>%
  # filter out low junction counts
  dplyr::filter(junction_ct >= 10)

rm(se_df)

# Load retained intron (RI) rMATS results
print("Loading RI rMATS...")
ri_df <- qs2::qs_read(ri_file) %>%
  #retain relevant columns
  dplyr::select(
    sample_id, geneSymbol, strand,
    chr,
    upstreamES, upstreamEE,
    downstreamES, downstreamEE,
    upstream_to_intron_count,
    intron_to_downstream_count,
    upstream_to_downstream_count
  ) %>%
  # remove suffix from sample ID, convert coordinates to one-based
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"),
                upstreamES = upstreamES + 1,
                downstreamES = downstreamES + 1) %>%
  dplyr::filter(!sample_id %in% samples_to_rm,
                (upstream_to_intron_count >= 10 | intron_to_downstream_count >= 10 | upstream_to_downstream_count >= 10))

# define columns specifying junction coordinates and retain only relevant columns
print("Extracting RI event junctions...")

ri_df <- define_junctions_targets(ri_df,
                                  event_type = "RI")

# Build the long-form junction table, merging data from:
# upstream-intron
# intron-downstream
# upstream-downstream junctions
print("Creating RI junction df...")

ri_junction_df <- create_junction_df(ri_df,
                                     event_type = "RI") %>%
  dplyr::filter(junction_ct >= 10)

rm(ri_df)

# A3SS events, modify sample ID and filter out cell samples
print("Loading A3SS rMATS...")
a3ss_df <- qs2::qs_read(a3ss_file) %>%
  # retain relevenat columns
  dplyr::select(
    sample_id, geneSymbol, strand,
    chr, shortES, shortEE,
    longExonStart_0base, longExonEnd,
    flankingES, flankingEE,
    long_to_flanking_count,
    short_to_flanking_count
  ) %>%
  # remove suffix from sample ID, convert coordinates to one-based
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"),
                flankingES = flankingES + 1,
                longExonStart_0base = longExonStart_0base + 1,
                shortES = shortES + 1) %>%
  dplyr::filter(!sample_id %in% samples_to_rm,
                (long_to_flanking_count >= 10 | short_to_flanking_count >= 10))

# define junction coordinates and filter for relevant columns
print("Extracting A3SS junctions...")
a3ss_df <- define_junctions_targets(a3ss_df,
                                    event_type = "A3SS")

# pivot longer for single row per unique sample & junction
print("Creating A3SS junction df...")
a3ss_junction_df <- create_junction_df(a3ss_df,
                                       event_type = "A3SS") %>%
  dplyr::filter(junction_ct >= 10)

rm(a3ss_df)

# A5SS events
print("Loading A5SS rMATS...")
a5ss_df <- qs2::qs_read(a5ss_file) %>%
  # select relevant columns
  dplyr::select(
    sample_id, geneSymbol, strand,
    chr,
    shortES, shortEE,
    longExonStart_0base, longExonEnd,
    flankingES, flankingEE,
    long_to_flanking_count,
    short_to_flanking_count
  ) %>%
  # remove suffix from sample ID, convert coordinates to one-based
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"),
                flankingES = flankingES + 1,
                longExonStart_0base = longExonStart_0base + 1,
                shortES = shortES + 1) %>%
  dplyr::filter(!sample_id %in% samples_to_rm,
                (long_to_flanking_count >= 10 | short_to_flanking_count >= 10))


# define junction coordinates and filter for relevant columns
print("Extracting A5SS junctions...")
a5ss_df <- define_junctions_targets(a5ss_df,
                                    event_type = "A5SS")

# pivot longer for single row per unique sample & junction
print("Creating A5SS junction df...")
a5ss_junction_df <- create_junction_df(a5ss_df,
                                       event_type = "A5SS") %>%
  dplyr::filter(junction_ct >= 10)

### Normalize junction counts

# merge junction dfs across event types
merged_junction_df <- se_junction_df %>%
  bind_rows(ri_junction_df,
            a3ss_junction_df,
            a5ss_junction_df) %>%
  dplyr::arrange(sample_id, junction)

# Load rMATS input read counts to normalize junction counts
read_cts <- read_tsv(read_file) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"))

# operate on merged junction df in chunks
n <- nrow(merged_junction_df)
chunk_size = 1e7

starts <- seq(1, n, by = chunk_size)

# create empty list to store junction cpms by chunk
norm_jc_chunk_list <- list()

# loop through chunks
print("Normalizing junction counts...")
system.time({
  for (i in 1:length(starts)) {

    print(glue::glue("processing chunk {i}..."))

    start <- starts[i]

    idx <- start:min(start + chunk_size - 1, n)

    jc_chunk <- merged_junction_df[idx,]

    # run `generate_norm_junction_mat` to obtain desired normalized mean cpm matrix
    norm_jc_chunk_list[[i]] <- generate_norm_junction_mat(jc_chunk,
                                                           read_cts,
                                                           group_col = "sample_id")

  }
})

# merge df list
merged_norm_junction_df <- bind_rows(norm_jc_chunk_list)

# Save merged junction output
qs2::qs_save(merged_norm_junction_df,
             file.path(results_dir,
                       "pbta-merged-norm-junction-cts.qs2"))

# create junction annotation df

print("creating junction annotation df...")
system.time({
  merged_junction_df <- se_junction_df %>%
    distinct(junction, geneSymbol, strand) %>%
    bind_rows(ri_junction_df %>%
                distinct(junction, geneSymbol, strand)) %>%
    bind_rows(a3ss_junction_df %>%
                distinct(junction, geneSymbol, strand)) %>%
    bind_rows(a5ss_junction_df %>%
                distinct(junction, geneSymbol, strand)) %>%
    distinct(junction, geneSymbol, strand) %>%
    # separate coordinates into individual columns
    dplyr::mutate(
      chr = unlist(lapply(strsplit(junction, ":|-|_"), function(x) x[[1]])),
      up_jc_start = unlist(lapply(strsplit(junction, ":|-|_"), function(x) x[[2]])),
      up_jc_end = unlist(lapply(strsplit(junction, ":|-|_"), function(x) x[[3]])),
      down_jc_start = unlist(lapply(strsplit(junction, ":|-|_"), function(x) x[[4]])),
      down_jc_end = unlist(lapply(strsplit(junction, ":|-|_"), function(x) x[[5]])),
                 )
})

# save junction annotation
write_tsv(merged_junction_df,
          file.path(results_dir,
                    "junction-annot.tsv.gz"))

# print session info
sessionInfo()
