# reformat rMATS data frames to include junction and target coordinate intervals, and filter for relevant columns
define_junctions_targets <- function(df, event_type){
  
  # processing will differ based on event_type
  if (event_type == "SE"){
    
    reformatted_df <- df  %>%
      dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"), 
                    exonStart_0base = exonStart_0base + 1,
                    upstreamES = upstreamES + 1,
                    downstreamES = downstreamES + 1) %>%
      # define inclusion junction, skip junction, and target coordinates
      dplyr::mutate(
        up_incl_jc = str_c(chr, ":", upstreamES, "-", upstreamEE, "_",
                           exonStart_0base, "-", exonEnd),
        down_incl_jc = str_c(chr, ":", exonStart_0base, "-", exonEnd, "_",
                             downstreamES, "-", downstreamEE),
        skip_jc  = str_c(chr, ":", upstreamES, "-", upstreamEE, "_",
                         downstreamES, "-", downstreamEE)
      ) %>%
      dplyr::filter((up_incl_jc %in% enr_jc_df$junction | down_incl_jc %in% enr_jc_df$junction | skip_jc %in% enr_jc_df$junction),
                    sample_id %in% enr_jc_df$sample_id) %>%
      dplyr::mutate(splice_id = glue::glue("{chr}:{exonStart_0base}-{exonEnd}_{upstreamES}-{upstreamEE}_{downstreamES}-{downstreamEE}_{strand}")) %>%
      # select sample, gene, coordinate, and count columns
      dplyr::select(sample_id, geneSymbol, up_incl_jc,
                    down_incl_jc, skip_jc,
                    strand,
                    splice_id,
                    IncLevel1)
    
  } else if (event_type == "RI"){
    
    # define inclusion junction, skip junction, and target coordinates
    reformatted_df <- df %>%
      dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"), 
                    upstreamES = upstreamES + 1,
                    downstreamES = downstreamES + 1) %>%
      dplyr::mutate(
        up_incl_jc = str_c(chr, ":", upstreamES, "-", upstreamEE, "_",
                           upstreamEE, "-", downstreamES),
        down_incl_jc = str_c(chr, ":", upstreamEE, "-", downstreamES, "_",
                             downstreamES, "-", downstreamEE),
        skip_jc = str_c(chr, ":", upstreamES, "-", upstreamEE, "_",
                        downstreamES, "-", downstreamEE)
      ) %>%
      dplyr::filter(up_incl_jc %in% enr_jc_df$junction | down_incl_jc %in% enr_jc_df$junction | skip_jc %in% enr_jc_df$junction) %>%
      dplyr::mutate(splice_id = glue::glue("{chr}:{riExonStart_0base}-{riExonEnd}_{upstreamES}-{upstreamEE}_{downstreamES}-{downstreamEE}_{strand}")) %>%
      # rename `intron_count` as `target_count` 
      dplyr::rename(target_count = intron_count) %>%
      dplyr::select(sample_id, geneSymbol, up_incl_jc,
                    down_incl_jc, skip_jc,
                    strand,
                    splice_id,
                    IncLevel1)
    
  } else if (event_type == "A3SS"){
    
    # define long and short inclusion junction coordinates
    reformatted_df <- df %>%
      dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"), 
                    longExonStart_0base = longExonStart_0base + 1,
                    shortES = shortES + 1,
                    flankingES = flankingES + 1) %>%
      dplyr::mutate(
        long_incl_jc = case_when(
          strand == "+" ~ str_c(chr, ":", flankingES, "-", flankingEE, "_",
                                longExonStart_0base, "-", longExonEnd),
          strand == "-" ~ str_c(chr, ":", longExonStart_0base, "-", longExonEnd, "_",
                                flankingES, "-", flankingEE))) %>%
      dplyr::mutate(
        short_incl_jc = case_when(
          strand == "+" ~ str_c(chr, ":", flankingES, "-", flankingEE, "_",
                                shortES, "-", shortEE),
          strand == "-" ~ str_c(chr, ":", shortES, "-", shortEE, "_",
                                flankingES, "-", flankingEE))
      ) %>%
      dplyr::filter(long_incl_jc %in% enr_jc_df$junction | short_incl_jc %in% enr_jc_df$junction) %>%
      dplyr::mutate(splice_id = glue::glue("{chr}:{longExonStart_0base}-{longExonEnd}_{shortES}-{shortEE}_{flankingES}-{flankingEE}_{strand}")) %>%
      # retain sample, gene, coordinates, and count columns
      dplyr::select(sample_id, geneSymbol, long_incl_jc,
                    short_incl_jc,
                    strand,
                    splice_id,
                    IncLevel1)
    
  } else if (event_type == "A5SS"){
    
    # define long and short inclusion junction coordinates
    reformatted_df <- df %>%
      dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$"), 
                    longExonStart_0base = longExonStart_0base + 1,
                    shortES = shortES + 1,
                    flankingES = flankingES + 1) %>%
      dplyr::mutate(
        long_incl_jc = case_when(
          strand == "+" ~ str_c(chr, ":", longExonStart_0base, "-", longExonEnd, "_",
                                flankingES, "-", flankingEE),
          strand == "-" ~ str_c(chr, ":", flankingES, "-", flankingEE, "_",
                                longExonStart_0base, "-", longExonEnd))
      ) %>%
      dplyr::mutate(
        short_incl_jc = case_when(
          strand == "+" ~ str_c(chr, ":", shortES, "-", shortEE, "_",
                                flankingES, "-", flankingEE),
          strand == "-" ~ str_c(chr, ":", flankingES, "-", flankingEE, "_",
                                shortES, "-", shortEE))
      ) %>%
      dplyr::filter(long_incl_jc %in% enr_jc_df$junction | short_incl_jc %in% enr_jc_df$junction) %>%
      dplyr::mutate(splice_id = glue::glue("{chr}:{longExonStart_0base}-{longExonEnd}_{shortES}-{shortEE}_{flankingES}-{flankingEE}_{strand}")) %>%
      # retain sample, gene, coordinates, and count columns
      dplyr::select(sample_id, geneSymbol, long_incl_jc,
                    short_incl_jc,
                    strand,
                    # long_to_flanking_count,
                    # short_to_flanking_count,
                    splice_id,
                    IncLevel1)
    
  }
  
  # convert to data.table
  reformatted_df <- as.data.table(reformatted_df)
  
  # return reformatted df
  return(reformatted_df)
  
}




# create data frame of rmats-dervied junction counts by sample
create_junction_df <- function(df, event_type){
  
  # processing will differ based on event_type
  if (event_type == "SE"){
    
    # create unique rows for each inclusion and skipping junction
    junction_df <- rbindlist(list(
      df[, .(sample_id, strand,
             junction = up_incl_jc,
             type = "exon inclusion",
             splice_id,
             IncLevel1
      )],
      df[, .(sample_id, strand,
             junction = down_incl_jc,
             type = "exon inclusion",
             splice_id,
             IncLevel1)],
      df[, .(sample_id, strand,
             junction = skip_jc,
             type = "exon skipping",
             splice_id,
             IncLevel1)]
    ))
    
  } else if (event_type == "RI"){
    
    # create unique rows for each inclusion and skipping junction
    junction_df <- rbindlist(list(
      df[, .(sample_id, strand,
             junction = up_incl_jc,
             type = "intron retention",
             splice_id,
             IncLevel1)],
      df[, .(sample_id, strand,
             junction = down_incl_jc,
             type = "intron retention",
             splice_id,
             IncLevel1)],
      df[, .(sample_id, strand,
             junction = skip_jc,
             type = "intron exclusion",
             splice_id,
             IncLevel1)]
    ))
    
  } else if (event_type %in% c("A5SS", "A3SS")){
    
    # create unique rows for each long and short inclusion junction
    junction_df <- rbindlist(list(
      df[, .(sample_id, strand,
             junction = long_incl_jc,
             type = ifelse(event_type == "A5SS", "A5SS+",
                           "A3SS+"),
             splice_id,
             IncLevel1)],
      df[, .(sample_id, strand,
             junction = short_incl_jc,
             type = ifelse(event_type == "A5SS", "A5SS-",
                           "A3SS-"),
             splice_id,
             IncLevel1)]
    ))
    
  }
  
  return(junction_df)
  
}
