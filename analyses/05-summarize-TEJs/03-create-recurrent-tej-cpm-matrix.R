## Create recurrent TEJ expression matrices (PBTA)
##
## Ryan Corbett
## Jan 2026
##
## This script subsets the all-PBTA junction CPM matrices to recurrent TEJs.
## The ComBat-corrected source matrix is stored as log2(CPM + 1), so it is
## converted back to CPM before saving.

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(rprojroot)
  library(qs2)
})

## Set directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "05-summarize-TEJs")
results_dir <- file.path(analysis_dir, "results")
pbta_results_dir <- file.path(
  root_dir,
  "analyses",
  "02-pbta-junction-processing",
  "results"
)

## Set file paths
recur_tej_file <- file.path(
  results_dir,
  "recurrent-primary-tumor-enriched-oncofetal-splice-junctions.tsv.gz"
)
all_cpm_file <- file.path(
  pbta_results_dir,
  "all-pbta-splice-junction-cpm.qs2"
)
all_combat_log2_cpm_file <- file.path(
  pbta_results_dir,
  "all-pbta-splice-junction-log2-cpm-combat-corrected.qs2"
)
recur_cpm_file <- file.path(
  results_dir,
  "tumor-enriched-oncofetal-splice-junction-cpm.qs2"
)
recur_combat_cpm_file <- file.path(
  results_dir,
  "tumor-enriched-oncofetal-splice-junction-cpm-combat-corrected.qs2"
)

## Load recurrent TEJ identities
recurrent_junctions <- readr::read_tsv(
  recur_tej_file,
  show_col_types = FALSE
)[["junction"]] |> unique()

if (!length(recurrent_junctions)) {
  stop("No recurrent TEJ junctions were found in: ", recur_tej_file)
}

subset_recurrent_junctions <- function(input_file, output_file,
                                       back_transform_log2 = FALSE) {
  message("Loading: ", input_file)
  junction_mat <- data.table::as.data.table(qs2::qs_read(input_file))

  if (!"junction" %in% names(junction_mat)) {
    stop("The input matrix must contain a 'junction' column: ", input_file)
  }

  sample_ids <- setdiff(names(junction_mat), "junction")
  if (!length(sample_ids)) {
    stop("The input matrix contains no sample columns: ", input_file)
  }

  data.table::setkey(junction_mat, junction)
  recurrent_mat <- junction_mat[
    recurrent_junctions,
    on = .(junction),
    nomatch = 0L
  ]

  missing_junctions <- setdiff(recurrent_junctions, recurrent_mat$junction)
  if (length(missing_junctions)) {
    warning(
      length(missing_junctions),
      " recurrent TEJ(s) were absent from ",
      basename(input_file),
      "."
    )
  }

  if (back_transform_log2) {
    message("Converting ComBat-corrected log2(CPM + 1) values back to CPM...")
    recurrent_mat[, (sample_ids) := lapply(
      .SD,
      function(x) round(pmax(2^x - 1, 0), 5)
    ), .SDcols = sample_ids]
  }

  message("Saving: ", output_file)
  qs2::qs_save(recurrent_mat, output_file)
}

## Subset uncorrected CPMs and ComBat-corrected CPMs.
subset_recurrent_junctions(all_cpm_file, recur_cpm_file)
subset_recurrent_junctions(
  all_combat_log2_cpm_file,
  recur_combat_cpm_file,
  back_transform_log2 = TRUE
)

sessionInfo()
