## Functions for creating PBTA junction CPM matrices

# Create a junction-by-group CPM matrix from rMATS-derived junction counts.
# Multiple rMATS records can represent the same junction in a sample; retain
# the maximum count before normalization, as in the PBTA junction-counting
# workflow.
generate_norm_junction_mat <- function(junction_df,
                                       read_cts,
                                       group_col = "sample_id") {
  junction_df <- as.data.table(junction_df)
  read_cts <- as.data.table(read_cts)

  if (!all(c("sample_id", "junction", "junction_ct") %in% names(junction_df))) {
    stop("junction_df must contain sample_id, junction, and junction_ct.")
  }
  if (!all(c("sample_id", "used_read_count") %in% names(read_cts))) {
    stop("read_cts must contain sample_id and used_read_count.")
  }
  if (!group_col %in% names(junction_df)) {
    stop(sprintf("group_col '%s' is not present in junction_df.", group_col))
  }

  # Collapse duplicate junction/sample records before the cast. Including
  # geneSymbol in this grouping would leave duplicate junction/sample pairs,
  # causing data.table::dcast() to default to a record-count aggregation.
  junction_df <- junction_df[
    ,
    .(junction_count = max(junction_ct, na.rm = TRUE)),
    by = eval(unique(c("sample_id", "junction", group_col)))
  ]

  junction_df <- merge(
    junction_df,
    read_cts[, .(sample_id, used_read_count)],
    by = "sample_id",
    all.x = TRUE
  )

  if (anyNA(junction_df$used_read_count)) {
    warning("Some junction records have no matching used_read_count.")
  }

  junction_df[, junction_cpm := junction_count / used_read_count * 1e6]

  dcast(
    junction_df,
    as.formula(paste("junction ~", group_col)),
    value.var = "junction_cpm",
    fun.aggregate = max,
    fill = NA_real_
  )
}
