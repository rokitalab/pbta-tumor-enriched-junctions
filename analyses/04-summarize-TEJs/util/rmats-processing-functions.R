## Functions for processing PBTA rMATS data
#
## Author: Ryan Corbett
#
## January 2026

# reformat rMATS data frames to include junction coordinates and filters for relevant columns
define_junctions_targets <- function(df, event_type){
  
  # processing will differ based on event_type
  if (event_type == "SE"){
    
    reformatted_df <- df %>%
      # define inclusion junction, skip junction, and target coordinates
      dplyr::mutate(
        up_incl_jc = str_c(chr, ":", upstreamES, "-", upstreamEE, "_",
                           exonStart_0base, "-", exonEnd),
        down_incl_jc = str_c(chr, ":", exonStart_0base, "-", exonEnd, "_",
                             downstreamES, "-", downstreamEE),
        skip_jc  = str_c(chr, ":", upstreamES, "-", upstreamEE, "_",
                         downstreamES, "-", downstreamEE)
      ) %>%
      # select sample, gene, coordinate, and count columns
      dplyr::select(sample_id, geneSymbol, up_incl_jc,
                    down_incl_jc,
                    skip_jc,
                    strand,
                    upstream_to_target_count,
                    target_to_downstream_count,
                    upstream_to_downstream_count)
    
  } else if (event_type == "RI"){
    
    # define inclusion junction, skip junction, and target coordinates
    reformatted_df <- df %>%
      dplyr::mutate(
        up_incl_jc = paste0(chr, ":", upstreamES, "-", upstreamEE, "_",
                            upstreamEE, "-", downstreamES),
        down_incl_jc = paste0(chr, ":", upstreamEE, "-", downstreamES, "_",
                              downstreamES, "-", downstreamEE),
        skip_jc = paste0(chr, ":", upstreamES, "-", upstreamEE, "_",
                         downstreamES, "-", downstreamEE)
      ) %>%
      # rename `intron_count` as `target_count` 
      dplyr::select(sample_id, geneSymbol, up_incl_jc,
                    down_incl_jc, 
                    skip_jc,
                    strand,
                    upstream_to_intron_count,
                    intron_to_downstream_count,
                    upstream_to_downstream_count)
    
  } else if (event_type == "A3SS"){
    
    # define long and short inclusion junction coordinates
    reformatted_df <- df %>%
      dplyr::mutate(
        long_incl_jc = case_when(
          strand == "+" ~ paste0(chr, ":", flankingES, "-", flankingEE, "_",
                                 longExonStart_0base, "-", longExonEnd),
          strand == "-" ~ paste0(chr, ":", longExonStart_0base, "-", longExonEnd, "_",
                                 flankingES, "-", flankingEE))) %>%
      dplyr::mutate(
        short_incl_jc = case_when(
          strand == "+" ~ paste0(chr, ":", flankingES, "-", flankingEE, "_",
                                 shortES, "-", shortEE),
          strand == "-" ~ paste0(chr, ":", shortES, "-", shortEE, "_",
                                 flankingES, "-", flankingEE))
      ) %>%
      # retain sample, gene, coordinates, and count columns
      dplyr::select(sample_id, geneSymbol, long_incl_jc,
                    short_incl_jc,
                    strand,
                    long_to_flanking_count,
                    short_to_flanking_count)
    
  } else if (event_type == "A5SS"){
    
    # define long and short inclusion junction coordinates
    reformatted_df <- df %>%
      dplyr::mutate(
        long_incl_jc = case_when(
          strand == "+" ~ paste0(chr, ":", longExonStart_0base, "-", longExonEnd, "_",
                                 flankingES, "-", flankingEE),
          strand == "-" ~ paste0(chr, ":", flankingES, "-", flankingEE, "_",
                                 longExonStart_0base, "-", longExonEnd))
      ) %>%
      dplyr::mutate(
        short_incl_jc = case_when(
          strand == "+" ~ paste0(chr, ":", shortES, "-", shortEE, "_",
                                 flankingES, "-", flankingEE),
          strand == "-" ~ paste0(chr, ":", flankingES, "-", flankingEE, "_",
                                 shortES, "-", shortEE))
      ) %>%
      # retain sample, gene, coordinates, and count columns
      dplyr::select(sample_id, geneSymbol, long_incl_jc,
                    short_incl_jc,
                    strand,
                    long_to_flanking_count,
                    short_to_flanking_count)
    
  }
  
  # convert to data.table
  reformatted_df <- as.data.table(reformatted_df)
  
  # return reformatted df
  return(reformatted_df)
  
}

# create data frame of rmats-derived target counts by sample
create_target_df <- function(df){
  
  target_df <- df[
    # in rare cases that target coordinates are duplicated, calculate max
    , .(target_count = max(target_count, na.rm = TRUE)),
    by = .(sample_id, geneSymbol, strand, target)
  ]
  
  return(target_df)
  
}

# create data frame of rmats-derived junction counts by sample
create_junction_df <- function(df, event_type){
  
  # processing will differ based on event_type
  if (event_type == "SE"){
    
    df <- df %>%
      dplyr::mutate(exonStart_0base = exonStart_0base + 1,
                    upstreamES = upstreamES + 1,
                    downstreamES = downstreamES + 1)
    
    # create unique rows for each inclusion and skipping junction
    junction_df <- rbindlist(list(
      df[, .(sample_id, geneSymbol, chr, strand,
             up_start = upstreamES,
             up_end = upstreamEE,
             down_start = exonStart_0base,
             down_end = exonEnd,
             junction_ct = upstream_to_target_count
      )],
      df[, .(sample_id, geneSymbol, chr, strand,
             up_start = exonStart_0base,
             up_end = exonEnd,
             down_start = downstreamES,
             down_end = downstreamEE,
             junction_ct = target_to_downstream_count)],
      df[, .(sample_id, geneSymbol, chr, strand,
             up_start = upstreamES,
             up_end = upstreamEE,
             down_start = downstreamES,
             down_end = downstreamEE,
             junction_ct = upstream_to_downstream_count)]
    ))
    
    jc_map <- junction_df %>%
      distinct(chr, up_start, up_end,
               down_start, down_end) %>%
      dplyr::mutate(up_start = as.integer(up_start),
                    up_end = as.integer(up_end),
                    down_start = as.integer(down_start),
                    down_end = as.integer(down_end)) %>%
      dplyr::mutate(junction = paste0(chr, ":", up_start, "-", up_end, "_",
                                      down_start, "-", down_end))
    
    junction_df <- junction_df %>%
      left_join(jc_map) %>%
      dplyr::select(-chr, -up_start,
                    -up_end, -down_start,
                    -down_end) %>%
      dplyr::select(sample_id, geneSymbol, strand, junction, junction_ct)
    
  } else if (event_type == "RI"){
    
    df <- df %>%
      dplyr::mutate(riExonStart_0base = riExonStart_0base + 1,
                    upstreamES = upstreamES + 1,
                    downstreamES = downstreamES + 1)
    
    # create unique rows for each inclusion and skipping junction
    junction_df <- rbindlist(list(
      df[, .(sample_id, geneSymbol, chr, strand,
             up_start = upstreamES,
             up_end = upstreamEE,
             down_start = upstreamEE,
             down_end = downstreamES,
             junction_ct = upstream_to_intron_count
      )],
      df[, .(sample_id, geneSymbol, chr, strand,
             up_start = upstreamEE,
             up_end = downstreamES,
             down_start = downstreamES,
             down_end = downstreamEE,
             junction_ct = intron_to_downstream_count)],
      df[, .(sample_id, geneSymbol, chr, strand,
             up_start = upstreamES,
             up_end = upstreamEE,
             down_start = downstreamES,
             down_end = downstreamEE,
             junction_ct = upstream_to_downstream_count)]
    ))
    
    jc_map <- junction_df %>%
      distinct(chr, up_start, up_end,
               down_start, down_end) %>%
      dplyr::mutate(up_start = as.integer(up_start),
                    up_end = as.integer(up_end),
                    down_start = as.integer(down_start),
                    down_end = as.integer(down_end)) %>%
      dplyr::mutate(junction = paste0(chr, ":", up_start, "-", up_end, "_",
                                      down_start, "-", down_end))
    
    junction_df <- junction_df %>%
      left_join(jc_map) %>%
      dplyr::select(-chr, -up_start,
                    -up_end, -down_start,
                    -down_end) %>%
      dplyr::select(sample_id, geneSymbol, strand, junction, junction_ct)
    
  } else if (event_type == "A3SS"){
    
    df <- df %>%
      dplyr::mutate(longExonStart_0base = longExonStart_0base + 1,
                    shortES = shortES + 1,
                    flankingES = flankingES + 1)
    
    junction_df <- rbindlist(list(
      df[, .(sample_id, geneSymbol, chr, strand,
             up_start = fifelse(strand == "+", flankingES, longExonStart_0base),
             up_end = fifelse(strand == "+", flankingEE, longExonEnd),
             down_start = fifelse(strand == "+", longExonStart_0base, flankingES),
             down_end = fifelse(strand == "+", longExonEnd, flankingEE),
             junction_ct = long_to_flanking_count
      )],
      df[, .(sample_id, geneSymbol, chr, strand,
             up_start = fifelse(strand == "+", flankingES, shortES),
             up_end = fifelse(strand == "+", flankingEE, shortEE),
             down_start = fifelse(strand == "+", shortES, flankingES),
             down_end = fifelse(strand == "+", shortEE, flankingEE),
             junction_ct = short_to_flanking_count)]
    ))
    
    jc_map <- junction_df %>%
      distinct(chr, up_start, up_end,
               down_start, down_end) %>%
      dplyr::mutate(up_start = as.integer(up_start),
                    up_end = as.integer(up_end),
                    down_start = as.integer(down_start),
                    down_end = as.integer(down_end)) %>%
      dplyr::mutate(junction = paste0(chr, ":", up_start, "-", up_end, "_",
                                      down_start, "-", down_end))
    
    junction_df <- junction_df %>%
      left_join(jc_map) %>%
      dplyr::select(-chr, -up_start,
                    -up_end, -down_start,
                    -down_end) %>%
      dplyr::select(sample_id, geneSymbol, strand, junction, junction_ct)
    
    
  } else if (event_type == "A5SS"){
    
    df <- df %>%
      dplyr::mutate(longExonStart_0base = longExonStart_0base + 1,
                    shortES = shortES + 1,
                    flankingES = flankingES + 1)
    
    junction_df <- rbindlist(list(
      df[, .(sample_id, geneSymbol, chr, strand,
             up_start = fifelse(strand == "+", longExonStart_0base, flankingES),
             up_end = fifelse(strand == "+", longExonEnd, flankingEE),
             down_start = fifelse(strand == "+", flankingES, longExonStart_0base),
             down_end = fifelse(strand == "+", flankingEE, longExonEnd),
             junction_ct = long_to_flanking_count
      )],
      df[, .(sample_id, geneSymbol, chr, strand,
             up_start = fifelse(strand == "+", shortES, flankingES),
             up_end = fifelse(strand == "+", shortEE, flankingEE),
             down_start = fifelse(strand == "+", flankingES, shortES),
             down_end = fifelse(strand == "+", flankingEE, shortEE),
             junction_ct = short_to_flanking_count)]
    ))
    
    jc_map <- junction_df %>%
      distinct(chr, up_start, up_end,
               down_start, down_end) %>%
      dplyr::mutate(up_start = as.integer(up_start),
                    up_end = as.integer(up_end),
                    down_start = as.integer(down_start),
                    down_end = as.integer(down_end)) %>%
      dplyr::mutate(junction = paste0(chr, ":", up_start, "-", up_end, "_",
                                      down_start, "-", down_end))
    
    junction_df <- junction_df %>%
      left_join(jc_map) %>%
      dplyr::select(-chr, -up_start,
                    -up_end, -down_start,
                    -down_end) %>%
      dplyr::select(sample_id, geneSymbol, strand, junction, junction_ct)
    
  }
  
  return(junction_df)
  
}

# generate mean normalized target cpm matrices 
generate_norm_target_mat <- function(target_df, 
                                     read_cts,
                                     group_col){
  
  # # Merge histology info
  # target_df <- merge(
  #   target_df,
  #   hist[, c(id_col, group_col), with = FALSE],
  #   by.x = "sample_id",
  #   by.y = id_col,
  #   all.x = TRUE
  # )
  # 
  # Merge read count info
  target_df <- merge(
    target_df,
    read_cts,
    by.x = "sample_id",
    by.y = "sample_id",
    all.x = TRUE
  )
  
  # calculate cpm
  target_df[, target_cpm :=
              target_count / used_read_count * 1000000
  ]
  
  final_target_df <- target_df %>%
    distinct(target, !!sym(group_col), .keep_all = TRUE) %>%
    dplyr::select(target, !!sym(group_col),
                  target_count, target_cpm)
  
  # Group and compute means
  # agg_target_df <- target_df[
  #   , .(mean_target_cpm = mean(target_cpm, na.rm = TRUE)),
  #   by = c(group_col, "target")
  # ]
  
  # Pivot wider
  # norm_target_mat <- dcast(
  #   target_df,
  #   as.formula(paste("target ~", group_col)),
  #   value.var = "target_cpm"
  # )
  
  # return mat
  return(final_target_df)
  
}


# generate mean normalized junction cpm matrices 
generate_norm_junction_mat <- function(junction_df, 
                                       read_cts,
                                       group_col){
  
  # in rare cases where junctions are duplicated, calculate max junction counts
  junction_df <- junction_df[
    , .(junction_count = max(junction_ct)),
    by = .(sample_id, geneSymbol, junction)
  ]
  
  # # Merge histology info
  # junction_df <- merge(
  #   junction_df,
  #   hist[, c(id_col, group_col), with = FALSE],
  #   by.x = "sample_id",
  #   by.y = id_col,
  #   all.x = TRUE
  # )
  
  # Merge read count info
  junction_df <- merge(
    junction_df,
    read_cts,
    by.x = "sample_id",
    by.y = "sample_id",
    all.x = TRUE
  )
  
  # calculate cpm
  junction_df[, junction_cpm :=
                junction_count / used_read_count * 1000000
  ]
  
  # Group and compute means
  # agg_junction_df <- junction_df[
  #   , .(mean_junction_cpm = mean(junction_cpm, na.rm = TRUE)),
  #   by = c(group_col, "junction")
  # ]
  
  # final_junction_df <- junction_df %>%
  #   #  distinct(junction, !!sym(group_col), .keep_all = TRUE) %>%
  #   dplyr::select(junction, !!sym(group_col),
  #                 junction_count, junction_cpm)
  
  # Wide format
  norm_junction_mat <- dcast(
    junction_df,
    as.formula(paste("junction ~", group_col)),
    value.var = "junction_cpm"
  )
  
  # norm_junction_mat <- junction_df_unique %>%
  #   distinct(junction, !!sym(group_col), .keep_all = TRUE) %>%
  #   pivot_wider(
  #     id_cols = -c(geneSymbol, junction_count, used_read_count),
  #     names_from = !!sym(group_col),
  #     values_from = junction_cpm
  #     )
  
  # return mat
  return(norm_junction_mat)
  
}

