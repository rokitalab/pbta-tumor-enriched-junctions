suppressPackageStartupMessages({
  library(tidyverse)
})


# Function to add TPM values to a dataframe from a matrix
add_TPM_values <- function(df, matrix) {
  # Ensure gene_symbol and sample_id columns exist in df
  if (!("gene_symbol" %in% names(df)) | !("sample_id" %in% names(df))) {
    stop("DataFrame must contain 'gene_symbol' and 'sample_id' as columns")
  }
  
  # Ensure rownames and colnames are properly set in the matrix
  if (is.null(rownames(matrix)) | is.null(colnames(matrix))) {
    stop("Matrix must have rownames and colnames set for gene_symbol and Sample_id, respectively")
  }
  
  expr_df <- matrix %>%
    as.data.frame() %>%
    rownames_to_column("gene_symbol") %>%
    dplyr::filter(gene_symbol %in% df$gene_symbol) %>%
    pivot_longer(-gene_symbol, names_to = "sample_id",
                 values_to = "gene_tpm")
  
  # Convert to lowercase and trim whitespace for consistency
  df <- df %>%
    # mutate(gene = gene_symbol,
    #        sample = sample_id) %>%
    left_join(expr_df, 
              by = c("sample_id", "gene_symbol"))

  # # Extract the matching TPM values
  # df$gene_tpm <- mapply(function(gene, sample) {
  #   if (gene %in% rownames(matrix) & sample %in% colnames(matrix)) {
  #     matrix[gene, sample]
  #   } else {
  #     NA
  #   }
  # }, df$gene_symbol, df$sample_id)
  
  return(df)
}