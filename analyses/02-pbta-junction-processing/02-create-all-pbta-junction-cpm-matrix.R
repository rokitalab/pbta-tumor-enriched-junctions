## Create all PBTA rMATS junction CPM matrix
##
## Ryan Corbett
##
## Sep 2026

# This script uses the junction identities in pbta-merged-norm-junction-cts.qs2
# to select junctions from the four raw PBTA rMATS event files and calculate
# per-sample CPMs. Unlike 03-create-recurrent-tej-cpm-matrix.R, it does not
# restrict junctions to TEJs or apply a junction read-count cutoff.

library(tidyverse)
library(qs2)
library(stringr)
library(data.table)

### Set directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "02-pbta-junction-processing")
results_dir <- file.path(analysis_dir, "results")

source(file.path(
  analysis_dir,
  "util",
  "junction-cpm-matrix.functions.R"
))

# Set file paths
file_list <- list(
  "A3SS" = file.path(data_dir, "pbta-rmats_merged_raw_A3SS.qs2"),
  "A5SS" = file.path(data_dir, "pbta-rmats_merged_raw_A5SS.qs2"),
  "SE" = file.path(data_dir, "pbta-rmats_merged_raw_SE.qs2"),
  "RI" = file.path(data_dir, "pbta-rmats_merged_raw_RI.qs2")
)
read_file <- file.path(data_dir, "pbta_input_read_counts.tsv")
hist_file <- file.path(
  root_dir,
  "analyses",
  "00-create-cohort-histologies",
  "results",
  "cohort-histologies.tsv"
)
pbta_junction_file <- file.path(
  results_dir,
  "pbta-merged-norm-junction-cts.qs2"
)
output_file <- file.path(
  results_dir,
  "all-pbta-splice-junction-cpm.qs2"
)
control_junction_files <- list(
  "GTEx" = file.path(
    root_dir,
    "analyses",
    "01-ctrl-rmats-processing",
    "results",
    "gtex-merged-norm-junction-ct-mat.qs2"
  ),
  "normal pediatric brain" = file.path(
    root_dir,
    "analyses",
    "01-ctrl-rmats-processing",
    "results",
    "normal-pedbrain-merged-norm-junction-ct-mat.qs2"
  )
)
control_cpm_cutoff <- 10

# Load normalization read counts, PBTA samples, and the complete set of
# junction identities retained in the normalized PBTA junction table.
read_cts <- read_tsv(read_file, show_col_types = FALSE) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"))

pbta_sample_ids <- read_tsv(hist_file, show_col_types = FALSE) %>%
  dplyr::pull(Kids_First_Biospecimen_ID)

print("Loading PBTA normalized junction counts...")
pbta_junctions <- qs2::qs_read(pbta_junction_file) %>%
  dplyr::distinct(junction) %>%
  dplyr::pull(junction)

# Exclude junctions expressed above the CPM cutoff in any GTEx sample or
# normal pediatric brain group before extracting rMATS data. Scan columns one
# at a time so this filter does not create a second full dense copy of either
# control CPM matrix.
control_junctions_above_cutoff <- unlist(lapply(
  names(control_junction_files),
  function(cohort) {
    print(glue::glue("Loading {cohort} control junction CPM matrix..."))
    control_mat <- as.data.table(qs2::qs_read(control_junction_files[[cohort]]))
    control_sample_cols <- setdiff(names(control_mat), "junction")
    control_above_cutoff <- rep(FALSE, nrow(control_mat))

    for (sample_col in control_sample_cols) {
      control_cpm <- control_mat[[sample_col]]
      control_above_cutoff <- control_above_cutoff |
        (!is.na(control_cpm) & control_cpm > control_cpm_cutoff)
    }

    control_mat[control_above_cutoff, junction]
  }
)) %>%
  unique()
gc()

print(glue::glue(
  "Removing {length(control_junctions_above_cutoff)} junctions with CPM > ",
  "{control_cpm_cutoff} in at least one control sample/group..."
))
pbta_junctions <- setdiff(pbta_junctions, control_junctions_above_cutoff)
rm(control_junctions_above_cutoff)
gc()

print(glue::glue(
  "Retaining {length(pbta_junctions)} PBTA junction identities for extraction."
))

# Create an empty list to store junctions from all event types.
jc_list <- list()

for (event in names(file_list)) {
  print(glue::glue("Loading {event} rMATS..."))
  
  rmats <- qs2::qs_read(file_list[[event]])
  
  rmats <- rmats %>%
    dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$")) %>%
    dplyr::filter(sample_id %in% pbta_sample_ids)
  
  n <- nrow(rmats)
  chunk_size <- 1e7
  starts <- seq(1, n, by = chunk_size)
  chunk_df_list <- list()
  
  print(glue::glue("Extracting {event} junctions in {length(starts)} chunk(s)..."))
  for (i in seq_along(starts)) {
    print(glue::glue("Processing {event} chunk {i}..."))
    idx <- starts[[i]]:min(starts[[i]] + chunk_size - 1, n)
    
    chunk_df_list[[i]] <- create_junction_df(
      rmats[idx, ],
      event_type = event
    ) %>%
      dplyr::filter(junction %in% pbta_junctions)
  }
  
  jc_list[[event]] <- bind_rows(chunk_df_list)
  rm(rmats, chunk_df_list)
  gc()
}

print("Merging splice-event junctions...")
merged_junction_df <- bind_rows(jc_list)
rm(jc_list)
gc()

# At this point all event types have been combined, so junction support can be
# counted across the complete PBTA data set without first creating a dense CPM
# matrix. Remove junctions with non-missing counts in fewer than three samples.
merged_junction_dt <- as.data.table(merged_junction_df)
rm(merged_junction_df)

junction_sample_counts <- merged_junction_dt[
  , .(n_non_na_samples = uniqueN(sample_id[!is.na(junction_ct)])),
  by = junction
]
low_support_junctions <- junction_sample_counts[
  n_non_na_samples < 3L,
  junction
]

print(glue::glue(
  "Removing {length(low_support_junctions)} junctions with non-NA counts ",
  "in fewer than three PBTA samples..."
))
merged_junction_dt <- merged_junction_dt[
  !junction %in% low_support_junctions
]
rm(junction_sample_counts, low_support_junctions)
gc()

# Generate the CPM matrix in junction chunks.  A full junction-by-sample cast
# would require materializing all >3M junctions at once; keeping each junction
# in exactly one chunk also ensures duplicate junction/sample records are
# collapsed correctly by generate_norm_junction_mat().
junction_chunk_size <- 100000L
setkey(merged_junction_dt, junction)

junction_ids <- unique(merged_junction_dt$junction)
chunk_starts <- seq.int(1L, length(junction_ids), by = junction_chunk_size)
all_junction_cpm_chunks <- vector("list", length(chunk_starts))

print(glue::glue(
  "Calculating all-junction CPMs in {length(chunk_starts)} chunks of up to ",
  "{junction_chunk_size} junctions..."
))
for (i in seq_along(chunk_starts)) {
  start <- chunk_starts[[i]]
  end <- min(start + junction_chunk_size - 1L, length(junction_ids))
  junction_chunk_ids <- junction_ids[start:end]
  
  print(glue::glue(
    "Normalizing junction chunk {i}/{length(chunk_starts)} ",
    "({length(junction_chunk_ids)} junctions)..."
  ))
  
  # A keyed data.table subset avoids repeatedly scanning all junction records.
  junction_chunk_df <- merged_junction_dt[
    junction_chunk_ids,
    on = .(junction),
    nomatch = 0L
  ]
  
  all_junction_cpm_chunks[[i]] <- generate_norm_junction_mat(
    junction_chunk_df,
    read_cts,
    group_col = "sample_id"
  )
}

all_junction_cpm_mat <- rbindlist(all_junction_cpm_chunks, use.names = TRUE)
rm(merged_junction_dt, all_junction_cpm_chunks)
gc()

print("Saving all-junction CPM matrix...")
qs2::qs_save(all_junction_cpm_mat, output_file)

sessionInfo()
