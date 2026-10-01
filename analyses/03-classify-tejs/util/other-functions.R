# Helper functions for tumor-enriched junction classification.

# Identify rMATS junctions observed in the control STAR SJ file without
# loading that multi-gigabyte file into R. rMATS stores the flanking exon
# coordinates, whereas STAR stores the intervening intron coordinates. The
# STAR intron starts one base after the upstream exon boundary and ends one
# base before the downstream exon boundary.
find_control_sj_matches <- function(junctions, sj_file) {
  if (!file.exists(sj_file)) {
    stop(glue::glue("Control SJ file does not exist: {sj_file}"))
  }

  sj_header <- readr::read_tsv(sj_file, n_max = 0, show_col_types = FALSE)
  required_sj_cols <- c("chr", "intron_start", "intron_end")
  if (!all(required_sj_cols %in% names(sj_header))) {
    stop("Control SJ file must contain chr, intron_start, and intron_end columns.")
  }

  candidate_coords <- tibble::tibble(junction = junctions) %>%
    dplyr::distinct(junction) %>%
    dplyr::transmute(
      junction,
      chr = sub(":.*", "", junction),
      upstream_exon_end = as.integer(sub(".*:[0-9]+-([0-9]+)_.*", "\\1", junction)),
      downstream_exon_start = as.integer(sub(".*_([0-9]+)-[0-9]+$", "\\1", junction))
    )

  if (nrow(candidate_coords) == 0) {
    return(tibble::tibble(
      junction = character(), chr = character(),
      upstream_exon_end = integer(), downstream_exon_start = integer(),
      control_intron_start = integer(), control_intron_end = integer()
    ))
  }

  control_lookup <- candidate_coords %>%
    dplyr::transmute(
      junction, chr, upstream_exon_end, downstream_exon_start,
      control_intron_start = upstream_exon_end + 1L,
      control_intron_end = downstream_exon_start - 1L
    ) %>%
    dplyr::distinct()

  lookup_file <- tempfile("control-sj-lookup-", fileext = ".tsv")
  matched_keys_file <- tempfile("control-sj-matches-", fileext = ".tsv")
  on.exit(unlink(c(lookup_file, matched_keys_file)), add = TRUE)

  readr::write_tsv(
    control_lookup %>%
      dplyr::distinct(chr, control_intron_start, control_intron_end),
    lookup_file,
    col_names = FALSE
  )

  # Print each matched coordinate once. The first input is the small lookup
  # table; the second is streamed line-by-line from the control SJ file.
  # Use gzip when appropriate so the multi-gigabyte input never needs to be
  # expanded on disk or loaded into R.
  awk_program <- paste(
    "BEGIN { FS=OFS=\"\\t\" }",
    "NR==FNR { wanted[$1 FS $2 FS $3]=1; next }",
    "FNR==1 { next }",
    "{ key=$2 FS $3 FS $4; if (key in wanted && !(key in seen)) { print key; seen[key]=1 } }"
  )
  sj_reader <- if (grepl("\\.gz$", sj_file, ignore.case = TRUE)) {
    paste("gzip -cd --", shQuote(sj_file))
  } else {
    paste("cat --", shQuote(sj_file))
  }
  scan_command <- paste(
    sj_reader, "| awk", shQuote(awk_program), shQuote(lookup_file), "-"
  )
  awk_status <- system2(
    "bash",
    args = c("-c", shQuote(scan_command)),
    stdout = matched_keys_file
  )
  if (awk_status != 0) {
    stop("Failed to scan the control SJ file for matching junctions.")
  }

  if (file.info(matched_keys_file)$size == 0) {
    return(tibble::tibble(
      junction = character(), chr = character(),
      upstream_exon_end = integer(), downstream_exon_start = integer(),
      control_intron_start = integer(), control_intron_end = integer()
    ))
  }

  matched_keys <- readr::read_tsv(
    matched_keys_file,
    col_names = c("chr", "control_intron_start", "control_intron_end"),
    col_types = readr::cols(
      chr = readr::col_character(),
      control_intron_start = readr::col_integer(),
      control_intron_end = readr::col_integer()
    )
  )

  control_lookup %>%
    dplyr::inner_join(
      matched_keys,
      by = c("chr", "control_intron_start", "control_intron_end")
    ) %>%
    dplyr::distinct()
}

# Collapse stage-specific postnatal Evo-Devo summaries into one reference
# group per brain region, using pooled group means and sample SDs.
collapse_evodevo_postnatal_groups <- function(mean_mat, sd_mat, group_sizes) {
  if (!setequal(mean_mat$junction, sd_mat$junction)) {
    stop("Evo-Devo mean and SD matrices do not contain the same junctions.")
  }

  sd_mat <- sd_mat[match(mean_mat$junction, sd_mat$junction), ]
  reference_groups <- c("Forebrain" = "postnatal-forebrain",
                        "Hindbrain" = "postnatal-hindbrain")
  collapsed_mat <- tibble(junction = mean_mat$junction)

  for (region in names(reference_groups)) {
    stage_groups <- group_sizes %>%
      dplyr::filter(region == .env$region) %>%
      dplyr::pull(evodevo_postnatal_group)

    if (!all(stage_groups %in% names(mean_mat)) ||
        !all(stage_groups %in% names(sd_mat))) {
      stop(glue::glue("Missing {region} stage group(s) in an Evo-Devo matrix."))
    }

    stage_n <- group_sizes$n[match(stage_groups,
                                   group_sizes$evodevo_postnatal_group)]
    # qs2 may preserve a data.table class without loading the corresponding
    # indexing method. Convert to a base data frame before selecting columns
    # so this works for either serialized object variant.
    select_stage_columns <- function(x, columns) {
      as.matrix(as.data.frame(x)[, columns, drop = FALSE])
    }
    stage_means <- select_stage_columns(mean_mat, stage_groups)
    stage_sds <- select_stage_columns(sd_mat, stage_groups)
    observed <- !is.na(stage_means)
    total_n <- rowSums(sweep(observed, 2, stage_n, `*`))

    # Pooled mean and sample SD across the original stage groups.
    stage_means[!observed] <- 0
    pooled_mean <- rowSums(sweep(stage_means, 2, stage_n, `*`)) / total_n
    stage_sds[is.na(stage_sds)] <- 0
    within_ss <- rowSums(sweep(stage_sds^2, 2, stage_n - 1, `*`) * observed)
    mean_deviation <- sweep(stage_means, 1, pooled_mean, `-`)
    mean_deviation[!observed] <- 0
    between_ss <- rowSums(sweep(mean_deviation^2, 2, stage_n, `*`))
    pooled_sd <- sqrt((within_ss + between_ss) / (total_n - 1))
    pooled_sd[total_n <= 1] <- NA_real_

    collapsed_mat[[paste0("mean_cpm_", reference_groups[[region]])]] <- pooled_mean
    collapsed_mat[[paste0("sd_cpm_", reference_groups[[region]])]] <- pooled_sd
  }

  collapsed_mat
}
