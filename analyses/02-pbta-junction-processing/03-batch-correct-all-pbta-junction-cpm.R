################################################################################
# 03-batch-correct-all-pbta-junction-cpm.R
# Batch correction using ComBat for RNA library prep effects
#
# Author: Ryan Corbett
# usage: Rscript --vanilla 03-batch-correct-all-pbta-junction-cpm.R
################################################################################

## Load packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(rprojroot)
  library(sva)
  library(qs2)
})

## Set up paths
root_dir <- find_root(has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "02-pbta-junction-processing")
results_dir <- file.path(analysis_dir, "results")

## Test settings. Set test_mode to FALSE for the complete matrix.
test_mode <- FALSE
test_n_junctions <- 1000000L

## Exclude junctions observed in fewer than this many samples.
min_non_na_values <- 10L
min_non_na_per_batch <- 2L

## Use a fixed seed so the randomized junction order is reproducible.
random_seed <- 1L

## Input and output files
cohort_hist_file <- file.path(
  root_dir,
  "analyses",
  "00-create-cohort-histologies",
  "results",
  "cohort-histologies.tsv"
)
input_file <- file.path(results_dir, "all-pbta-splice-junction-cpm.qs2")
pbta_merged_file <- file.path(
  results_dir,
  "pbta-merged-norm-junction-cts.qs2"
)
output_file <- file.path(
  results_dir,
  if (test_mode) {
    "all-pbta-splice-junction-log2-cpm-combat-corrected-test.qs2"
  } else {
    "all-pbta-splice-junction-log2-cpm-combat-corrected.qs2"
  }
)
batch_corrected_pbta_file <- file.path(
  results_dir,
  "pbta-merged-norm-batch-corrected-junction-cts.qs2"
)

## Process 50,000 junctions in each ComBat chunk.
chunk_size <- 50000L

## Load the junction CPM matrix and cohort annotations
message("Loading all-PBTA junction CPM matrix...")
junction_cpm <- as.data.table(qs2::qs_read(input_file))

if (!"junction" %in% names(junction_cpm)) {
  stop("The input matrix must contain a 'junction' column.")
}

sample_ids <- setdiff(names(junction_cpm), "junction")
if (!length(sample_ids)) {
  stop("The input matrix contains no sample columns.")
}

## Randomize row indices without copying the full CPM matrix. The randomized
## order is used for both test-row selection and ComBat chunk formation.
set.seed(random_seed)
junction_indices <- sample.int(nrow(junction_cpm))

if (test_mode) {
  n_test_junctions <- min(test_n_junctions, length(junction_indices))
  message(
    "Test mode: retaining ", n_test_junctions,
    " randomly selected junctions."
  )
  junction_indices <- junction_indices[seq_len(n_test_junctions)]
}

cohort_hist <- readr::read_tsv(cohort_hist_file, show_col_types = FALSE)
required_columns <- c("Kids_First_Biospecimen_ID", "RNA_library")
if (!all(required_columns %in% names(cohort_hist))) {
  stop(
    "The cohort-histologies file must contain: ",
    paste(required_columns, collapse = ", "),
    "."
  )
}

## Prepare the RNA-library batch variable in the exact order of matrix columns.
batch_lookup <- cohort_hist %>%
  dplyr::select(Kids_First_Biospecimen_ID, RNA_library) %>%
  dplyr::distinct(Kids_First_Biospecimen_ID, .keep_all = TRUE)

batch <- batch_lookup$RNA_library[
  match(sample_ids, batch_lookup$Kids_First_Biospecimen_ID)
]
names(batch) <- sample_ids

if (anyNA(batch)) {
  missing_samples <- names(batch)[is.na(batch)]
  stop(
    "RNA library type is missing for ",
    length(missing_samples),
    " matrix sample(s): ",
    paste(missing_samples, collapse = ", "),
    "."
  )
}

batch <- as.factor(batch)
batches <- lapply(levels(batch), function(level) which(batch == level))

## ComBat cannot estimate a scale adjustment for a junction with zero
## within-batch variance. Rows with this condition are retained in the output
## but left as uncorrected log2(CPM + 1), consistent with ComBat's treatment
## of uniform-expression rows.
combat_eligible_rows <- function(cpm_mat, batches) {
  eligible <- rep(TRUE, nrow(cpm_mat))

  for (sample_indices in batches) {
    batch_variance <- apply(
      cpm_mat[, sample_indices, drop = FALSE],
      1L,
      stats::var,
      na.rm = TRUE
    )
    eligible <- eligible & is.finite(batch_variance) & batch_variance > 0
  }

  eligible
}

## Filter sparsely observed junctions before batch correction. Calculate
## non-missing counts in blocks to avoid allocating a full logical matrix.
filter_chunk_size <- 10000L
filter_starts <- seq.int(1L, length(junction_indices), by = filter_chunk_size)
keep_rows <- logical(length(junction_indices))

message(
  "Removing junctions with fewer than ", min_non_na_values,
  " non-NA CPM values or fewer than ", min_non_na_per_batch,
  " non-NA values in any RNA-library batch..."
)
for (i in seq_along(filter_starts)) {
  start <- filter_starts[[i]]
  end <- min(start + filter_chunk_size - 1L, length(junction_indices))
  row_indices <- junction_indices[start:end]
  cpm_filter_chunk <- as.matrix(junction_cpm[row_indices, ..sample_ids])
  total_non_na <- rowSums(!is.na(cpm_filter_chunk))
  batch_non_na <- vapply(
    batches,
    function(sample_indices) {
      rowSums(!is.na(cpm_filter_chunk[, sample_indices, drop = FALSE]))
    },
    numeric(nrow(cpm_filter_chunk))
  )
  keep_rows[start:end] <- total_non_na >= min_non_na_values &
    rowSums(batch_non_na >= min_non_na_per_batch) == length(batches)
}

n_removed <- sum(!keep_rows)
junction_indices <- junction_indices[keep_rows]
rm(keep_rows, cpm_filter_chunk)
gc()

message(
  "Retained ", length(junction_indices), " junctions; removed ", n_removed,
  " failing the non-NA observation thresholds."
)
if (!length(junction_indices)) {
  stop("No junctions remain after filtering sparse rows.")
}

## Run ComBat separately for junction chunks to limit peak memory use.
## Randomizing junction_indices makes each chunk broadly representative.
n_junctions <- length(junction_indices)
chunk_starts <- seq.int(1L, n_junctions, by = chunk_size)
combat_chunks <- vector("list", length(chunk_starts))

message(
  "Batch-correcting ", n_junctions, " junctions in ",
  length(chunk_starts), " chunk(s) of up to ", chunk_size, " rows..."
)

for (i in seq_along(chunk_starts)) {
  start <- chunk_starts[[i]]
  end <- min(start + chunk_size - 1L, n_junctions)
  row_indices <- junction_indices[start:end]
  message(
    "Processing chunk ", i, "/", length(chunk_starts),
    " (", length(row_indices), " randomized junctions)..."
  )

  cpm_chunk <- as.matrix(junction_cpm[row_indices, ..sample_ids])
  storage.mode(cpm_chunk) <- "double"
  log2_cpm_chunk <- log2(cpm_chunk + 1)
  eligible_rows <- combat_eligible_rows(cpm_chunk, batches)
  combat_chunk <- log2_cpm_chunk

  if (any(eligible_rows)) {
    combat_chunk[eligible_rows, ] <- sva::ComBat(
      dat = log2_cpm_chunk[eligible_rows, , drop = FALSE],
      batch = batch,
      par.prior = TRUE
    )
  }
  if (any(!eligible_rows)) {
    message(
      "Leaving ", sum(!eligible_rows),
      " junctions uncorrected because at least one library batch has ",
      "zero within-batch variance."
    )
  }

  combat_chunks[[i]] <- data.table::as.data.table(round(combat_chunk, 5))
  data.table::setnames(combat_chunks[[i]], sample_ids)
  combat_chunks[[i]][, junction := junction_cpm$junction[row_indices]]
  data.table::setcolorder(combat_chunks[[i]], c("junction", sample_ids))

  rm(cpm_chunk, log2_cpm_chunk, eligible_rows, combat_chunk)
  gc()
}

combat_mat <- data.table::rbindlist(combat_chunks, use.names = TRUE)
rm(junction_cpm, combat_chunks)
gc()

message("Saving ComBat-corrected log2 CPM matrix...")
qs2::qs_save(combat_mat, output_file)

## Replace the long-form PBTA junction CPM values with their ComBat-corrected
## counterparts. ComBat operates on log2(CPM + 1), so transform its output
## back to CPM here. Bound small negative values introduced by the adjustment
## at zero, as CPMs cannot be negative.
message("Loading long-form PBTA junction CPM table for replacement...")
pbta_merged <- as.data.table(qs2::qs_read(pbta_merged_file))
required_pbta_columns <- c("junction", "sample_id", "junction_cpm")
if (!all(required_pbta_columns %in% names(pbta_merged))) {
  stop(
    "The PBTA merged junction table must contain: ",
    paste(required_pbta_columns, collapse = ", "),
    "."
  )
}

message("Converting ComBat-corrected log2 CPMs back to CPMs...")
combat_mat[, (sample_ids) := lapply(
  .SD,
  function(x) round(pmax(2^x - 1, 0), 5)
), .SDcols = sample_ids]

## The ComBat matrix includes only rows that passed the observation filters.
## Rows not present in it (for example, sparse junctions) retain their
## original CPM values in the final table.
message("Replacing PBTA junction CPMs with batch-corrected CPMs...")
corrected_cpm_long <- data.table::melt(
  combat_mat,
  id.vars = "junction",
  variable.name = "sample_id",
  value.name = "batch_corrected_junction_cpm",
  na.rm = TRUE,
  variable.factor = FALSE
)
data.table::setkey(pbta_merged, junction, sample_id)
data.table::setkey(corrected_cpm_long, junction, sample_id)
pbta_merged[corrected_cpm_long,
  junction_cpm := i.batch_corrected_junction_cpm,
  on = .(junction, sample_id)
]

message("Saving PBTA junction table with batch-corrected CPMs...")
qs2::qs_save(pbta_merged, batch_corrected_pbta_file)

sessionInfo()
