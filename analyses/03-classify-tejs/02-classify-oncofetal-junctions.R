## Classify tumor-enriched PBTA junctions as oncofetal
##
## Ryan Corbett
## Aug 2026

# Load packages used for data manipulation and serialized matrices.
library(tidyverse)
library(qs2)

# Define analysis directories and input/output file paths.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "03-classify-tejs")
results_dir <- file.path(analysis_dir, "results")
junction_annotation_file <- file.path(
  root_dir,
  "analyses",
  "02-pbta-junction-processing",
  "results",
  "junction-annot.tsv.gz"
)

tumor_enriched_file <- file.path(results_dir, "tumor-enriched-junctions.qs2")
prenatal_junction_mat_file <- file.path(root_dir, "analyses",
                                        "01-ctrl-rmats-processing", "results",
                                        "evodevo-merged-prenatal-week-binned-norm-junction-ct-mat.qs2")
postnatal_junction_mat_file <- file.path(root_dir, "analyses",
                                         "01-ctrl-rmats-processing", "results",
                                         "evodevo-merged-postnatal-stage-norm-junction-ct-mat.qs2")
postnatal_junction_sd_file <- file.path(root_dir, "analyses",
                                        "01-ctrl-rmats-processing", "results",
                                        "evodevo-merged-postnatal-stage-norm-junction-sd-mat.qs2")

print("Loading tumor-enriched junctions and Evo-Devo matrices...")
# Load stage-01 calls and postnatal Evo-Devo CPM/SD matrices by brain region and stage.
merged_enr_jc_df <- qs2::qs_read(tumor_enriched_file)
evodevo_junction_mat <- qs2::qs_read(postnatal_junction_mat_file) %>%
  rename_with(~ paste0("mean_cpm_", .x), -junction)
evodevo_sd_mat <- qs2::qs_read(postnatal_junction_sd_file) %>%
  rename_with(~ paste0("sd_cpm_", .x), -junction)

# Restrict prenatal/postnatal comparisons to tumor-enriched junctions.
junctions_to_classify <- merged_enr_jc_df %>%
  dplyr::filter(junction_preference == "Tumor-enriched") %>%
  dplyr::distinct(junction)

# Load prenatal CPMs for the candidate junctions.
prenatal_junction_mat <- qs2::qs_read(prenatal_junction_mat_file) %>%
  dplyr::filter(junction %in% junctions_to_classify$junction) %>%
  rename_with(~ paste0("prenatal_mean_cpm_", .x), -junction)

# Identify the corresponding prenatal CPM and postnatal CPM/SD columns.
postnatal_mean_cols <- grep("^mean_cpm_", names(evodevo_junction_mat), value = TRUE)
postnatal_sd_cols <- sub("^mean_cpm_", "sd_cpm_", postnatal_mean_cols)
prenatal_mean_cols <- grep("^prenatal_mean_cpm_", names(prenatal_junction_mat), value = TRUE)
prenatal_groups <- sub("^prenatal_mean_cpm_", "", prenatal_mean_cols)

# Assemble one comparison table containing prenatal and postnatal measurements.
prenatal_vs_postnatal_df <- prenatal_junction_mat %>%
  left_join(evodevo_junction_mat %>%
              dplyr::filter(junction %in% junctions_to_classify$junction) %>%
              dplyr::select(junction, all_of(postnatal_mean_cols)),
            by = "junction") %>%
  left_join(evodevo_sd_mat %>%
              dplyr::filter(junction %in% junctions_to_classify$junction) %>%
              dplyr::select(junction, all_of(postnatal_sd_cols)),
            by = "junction") %>%
  as.data.frame()

# A junction must have an observed postnatal CPM in at least one group to be
# eligible for oncofetal classification. Keep postnatal NAs as missing.
prenatal_vs_postnatal_df$has_postnatal_cpm <- rowSums(
  !is.na(prenatal_vs_postnatal_df[, postnatal_mean_cols, drop = FALSE])
) > 0

# For each prenatal group, calculate the minimum FC and SNR across every
# postnatal group without retaining the intermediate comparison columns.
prenatal_min_fc_cols <- character(length(prenatal_groups))
prenatal_min_snr_cols <- character(length(prenatal_groups))
oncofetal_call_cols <- character(length(prenatal_groups))

for (k in seq_along(prenatal_groups)) {
  # Compare one prenatal group with every postnatal group.
  prenatal_group <- prenatal_groups[k]
  prenatal_mean_col <- prenatal_mean_cols[k]
  prenatal_cpm <- prenatal_vs_postnatal_df[[prenatal_mean_col]]
  min_fc <- rep(Inf, nrow(prenatal_vs_postnatal_df))
  min_snr <- rep(Inf, nrow(prenatal_vs_postnatal_df))
  n_snr <- integer(nrow(prenatal_vs_postnatal_df))
  
  for (j in seq_along(postnatal_mean_cols)) {
    postnatal_cpm <- prenatal_vs_postnatal_df[[postnatal_mean_cols[j]]]
    fc <- prenatal_cpm / (postnatal_cpm + 1e-5)
    snr <- (prenatal_cpm - postnatal_cpm) /
      prenatal_vs_postnatal_df[[postnatal_sd_cols[j]]]
    min_fc <- pmin(min_fc, fc, na.rm = TRUE)
    min_snr <- pmin(min_snr, snr, na.rm = TRUE)
    n_snr <- n_snr + !is.na(snr)
  }

  # Missing prenatal CPM is not evidence of prenatal expression. Likewise,
  # junctions absent from every postnatal group are not eligible for an
  # oncofetal call. Without these guards, pmin(..., na.rm = TRUE) retains the
  # Inf initializer for rows without valid comparisons.
  no_comparison <- is.na(prenatal_cpm) |
      !prenatal_vs_postnatal_df$has_postnatal_cpm
  min_fc[no_comparison] <- NA_real_
  
  # A postnatal mean can exist without any postnatal SD (e.g. every observed
  # group has a single sample), so the SNR criterion cannot be evaluated.
  min_snr[no_comparison | n_snr == 0] <- NA_real_

  prenatal_min_fc_cols[k] <- paste0("min_prenatal_fc_", prenatal_group)
  prenatal_min_snr_cols[k] <- paste0("min_prenatal_snr_", prenatal_group)
  oncofetal_call_cols[k] <- paste0("oncofetal_", prenatal_group)

  # Store the per-prenatal-group minima and threshold call.
  prenatal_vs_postnatal_df[[prenatal_min_fc_cols[k]]] <- min_fc
  prenatal_vs_postnatal_df[[prenatal_min_snr_cols[k]]] <- min_snr
  prenatal_vs_postnatal_df[[oncofetal_call_cols[k]]] <-
    prenatal_vs_postnatal_df[[prenatal_min_fc_cols[k]]] > 2 &
    prenatal_vs_postnatal_df[[prenatal_min_snr_cols[k]]] > 2
}

# Summarize the strongest qualifying prenatal group for each junction.
oncofetal_calls <- prenatal_vs_postnatal_df %>%
  dplyr::select(junction, all_of(prenatal_min_fc_cols),
                all_of(prenatal_min_snr_cols), all_of(oncofetal_call_cols))
prenatal_min_fc_df <- as.data.frame(
  oncofetal_calls[, prenatal_min_fc_cols, drop = FALSE]
)
prenatal_min_snr_df <- as.data.frame(
  oncofetal_calls[, prenatal_min_snr_cols, drop = FALSE]
)
oncofetal_call_df <- as.data.frame(
  oncofetal_calls[, oncofetal_call_cols, drop = FALSE]
)
oncofetal_calls$max_prenatal_min_cpm_fc <-
  do.call(pmax, c(as.list(prenatal_min_fc_df), na.rm = TRUE))
oncofetal_calls$max_prenatal_min_cpm_snr <-
  do.call(pmax, c(as.list(prenatal_min_snr_df), na.rm = TRUE))
oncofetal_calls$oncofetal_prenatal_group <- apply(
  oncofetal_call_df, 1,
  function(x) str_c(prenatal_groups[which(as.logical(x) %in% TRUE)], collapse = ";")
)
oncofetal_calls <- oncofetal_calls %>%
  dplyr::mutate(oncofetal_prenatal_group = na_if(oncofetal_prenatal_group, "")) %>%
  dplyr::select(junction, max_prenatal_min_cpm_fc,
                max_prenatal_min_cpm_snr, oncofetal_prenatal_group)

# Add prenatal summary statistics and relabel qualifying junctions as oncofetal.
merged_enr_jc_df <- merged_enr_jc_df %>%
  left_join(oncofetal_calls, by = "junction") %>%
  dplyr::mutate(junction_preference = case_when(
    junction_preference == "Tumor-enriched" & !is.na(oncofetal_prenatal_group) ~ "Oncofetal",
    TRUE ~ junction_preference
  ))

# Load genomic junction annotations and merge them with the classified calls.
jc_annot <- read_tsv(junction_annotation_file) %>%
  dplyr::filter(junction %in% merged_enr_jc_df$junction) %>%
  distinct(junction, geneSymbol, .keep_all = TRUE) %>%
  group_by(junction, strand, chr, up_jc_start, up_jc_end, down_jc_start,
           down_jc_end) %>%
  summarise(geneSymbol = str_c(geneSymbol, collapse = ","), .groups = "drop")

merged_enr_jc_annot_df <- merged_enr_jc_df %>%
  dplyr::select(junction, sample_id, junction_count, junction_cpm, boundary,
                criteria, junction_preference, min_cpm_fc_all, min_cpm_snr_all,
                max_mean_cpm_all, max_prenatal_min_cpm_fc,
                max_prenatal_min_cpm_snr, oncofetal_prenatal_group) %>%
  left_join(jc_annot, by = "junction")

# Write the annotated junction table for downstream analyses.
write_tsv(merged_enr_jc_annot_df,
          file.path(results_dir, "tumor-enriched-oncofetal-splice-junctions.tsv.gz"))

# Create and write a BED file for downstream protein-domain annotation.
jc_bed_df <- merged_enr_jc_annot_df %>%
  dplyr::mutate(up_jc_end = as.integer(up_jc_end),
                down_jc_start = as.integer(down_jc_start)) %>%
  dplyr::select(chr, up_jc_end, down_jc_start, junction, junction_preference,
                strand) %>%
  dplyr::rename(start = up_jc_end, end = down_jc_start,
                preference = junction_preference, id = junction) %>%
  distinct()

write.table(jc_bed_df,
            file.path(results_dir, "tumor-enriched-oncofetal-splice-junctions.bed"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# Report unique junction classification counts and infinite prenatal/postnatal
# comparison statistics among oncofetal junctions.
unique_classified_junctions <- merged_enr_jc_annot_df %>%
  dplyr::distinct(junction, junction_preference,
                  max_prenatal_min_cpm_fc, max_prenatal_min_cpm_snr)

classification_counts <- unique_classified_junctions %>%
  dplyr::filter(junction_preference %in% c("Tumor-enriched", "Oncofetal")) %>%
  dplyr::count(junction_preference, name = "n_unique_junctions") %>%
  tidyr::complete(junction_preference = c("Tumor-enriched", "Oncofetal"),
                  fill = list(n_unique_junctions = 0L))

oncofetal_infinite_counts <- unique_classified_junctions %>%
  dplyr::filter(junction_preference == "Oncofetal") %>%
  dplyr::summarise(
    n_oncofetal_junctions = dplyr::n(),
    n_with_infinite_fc = sum(is.infinite(max_prenatal_min_cpm_fc)),
    n_with_infinite_snr = sum(is.infinite(max_prenatal_min_cpm_snr)),
    n_with_infinite_fc_or_snr = sum(
      is.infinite(max_prenatal_min_cpm_fc) |
        is.infinite(max_prenatal_min_cpm_snr)
    )
  )

print(classification_counts)
print(oncofetal_infinite_counts)

# Record package and R versions used for this run.
sessionInfo()
