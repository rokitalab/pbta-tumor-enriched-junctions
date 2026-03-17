# Function to annotate splice junctions to transcripts, exons, introns

annotate_junctions <- function(chr = chr,
                               up_exon_start = up_exon_start,
                               up_exon_end = up_exon_end,
                               down_exon_start = down_exon_start,
                               down_exon_end = down_exon_end,
                               gene_symbol = gene_symbol,
                               preference = preference,
                               gtf = gtf) {
  
  # convert gtf df to data table for easier filtering
  gtf_dt <- as.data.table(gtf)
  
  setkey(gtf_dt, seqnames, start, end, gene_name)
  
  # Get GTF entries for upstream and downstream exons
  up_exon_gtf <- gtf_dt[.(chr, up_exon_start, up_exon_end, gene_symbol)]
  up_exon_gtf <- up_exon_gtf[!is.na(up_exon_gtf$width)]
  
  down_exon_gtf <- gtf_dt[.(chr, down_exon_start, down_exon_end, gene_symbol)]
  down_exon_gtf <- down_exon_gtf[!is.na(down_exon_gtf$width)]
  
  # Scenario 1: junction is RI
  if (preference == "RI") {
    
    # if adjacent exon is not annotated:
    if (nrow(up_exon_gtf) == 0 & nrow(down_exon_gtf) == 0) {
      
      return("RI event adjacent to novel exon")
      
      # if adjacent exon is upstream and annotated:
    } else if (nrow(up_exon_gtf) > 0 & nrow(down_exon_gtf) == 0) {
      
      # choose top TSL transcript as annotation if available
      exon_best_TSL <- up_exon_gtf %>%
        dplyr::arrange(transcript_support_level,
                       transcript_name) %>%
        head(n = 1) %>%
        dplyr::select(strand, transcript_name,
                      exon_number)
      
      # get intron number from exon number and strand
      intron_number <- ifelse(exon_best_TSL$strand == "+", 
                              exon_best_TSL$exon_number,
                              as.numeric(exon_best_TSL$exon_number) - 1)
      
      junction_annot <- glue::glue("{exon_best_TSL$transcript_name}, exon{exon_best_TSL$exon_number}:intron{intron_number}")
      return(junction_annot)
      
      # if adjacent exon is downstream and annotated:
    } else if (nrow(up_exon_gtf) == 0 & nrow(down_exon_gtf) > 0) {
      
      exon_best_TSL <- down_exon_gtf %>%
        dplyr::arrange(transcript_support_level,
                       transcript_name) %>%
        head(n = 1) %>%
        dplyr::select(strand, transcript_name,
                      exon_number)
      
      intron_number <- ifelse(exon_best_TSL$strand == "+", 
                              exon_best_TSL$exon_number,
                              as.numeric(exon_best_TSL$exon_number) - 1)
      
      junction_annot <- glue::glue("{exon_best_TSL$transcript_name}, intron{intron_number}:exon{exon_best_TSL$exon_number}")
      return(junction_annot)

    }
    
  }
  
  # For non-RI events: determine if exons are part of same transcript
  common_transcripts <- intersect(up_exon_gtf$transcript_name,
                                  down_exon_gtf$transcript_name)
  
  # Scenario 2: exons share single common transcript
  if (length(common_transcripts) == 1) {
    
    # filter exon for annotation in common transcript isoform
    up_exon_number <- up_exon_gtf$exon_number[
      up_exon_gtf$transcript_name == common_transcripts
    ]
    
    down_exon_number <- down_exon_gtf$exon_number[
      down_exon_gtf$transcript_name == common_transcripts
    ]
    
    # assign junction annotation
    junction_annot <- glue::glue("{common_transcripts}, exon{up_exon_number}:exon{down_exon_number}")
    return(junction_annot)

  } 
  
  if (length(common_transcripts) > 1){
    
    # Scenario 3: exons each map to single protein-coding transcript
    common_pc_transcripts <- intersect(up_exon_gtf$transcript_name[up_exon_gtf$transcript_type == "protein_coding"],
                                       down_exon_gtf$transcript_name[down_exon_gtf$transcript_type == "protein_coding"])
    
    if(length(common_pc_transcripts) == 1){
      
      # filter exon numbers for those in single protein-coding transcript
      up_exon_number <- up_exon_gtf$exon_number[
        up_exon_gtf$transcript_name == common_pc_transcripts
      ]
      
      down_exon_number <- down_exon_gtf$exon_number[
        down_exon_gtf$transcript_name == common_pc_transcripts
      ]
      
      junction_annot <- glue::glue("{common_pc_transcripts}, exon{up_exon_number}:exon{down_exon_number}")
      return(junction_annot)

    }
    
    # Scenario 4: junction is annotated as exon skipping --> take transcript that supports skipping with highest TSL
    if (preference == "ES"){
      
      # identify transcripts with at least on exon between junction exons
      common_transcripts_es_df <- up_exon_gtf %>%
        dplyr::filter(transcript_name %in% common_transcripts) %>%
        dplyr::select(transcript_name, transcript_support_level,
                      exon_number) %>%
        dplyr::rename(up_exon_number = exon_number) %>%
        left_join(down_exon_gtf %>%
                    dplyr::select(transcript_name, 
                                  exon_number) %>%
                    dplyr::rename(down_exon_number = exon_number)) %>%
        dplyr::filter(abs(as.numeric(up_exon_number) - as.numeric(down_exon_number)) > 1) 
      
      # take transcript with best TSL, if multiple
      if (nrow(common_transcripts_es_df) > 0){
        
        common_transcripts_es_best_tsl <- common_transcripts_es_df %>%
          dplyr::arrange(transcript_support_level,
                         transcript_name) %>%
          head(n = 1) %>%
          pull(transcript_name)
        
        up_exon_number <- up_exon_gtf$exon_number[
          up_exon_gtf$transcript_name == common_transcripts_es_best_tsl
        ]
        
        down_exon_number <- down_exon_gtf$exon_number[
          down_exon_gtf$transcript_name == common_transcripts_es_best_tsl
        ]
        
        junction_annot <- glue::glue("{common_transcripts_es_best_tsl}, exon{up_exon_number}:exon{down_exon_number}")
        return(junction_annot)

      }
      
    }
    
    # Scenario 5: exons map to multiple common transcripts --> take transcript with highest TSL (1 = most supported)
    common_transcripts_best_tsl <- up_exon_gtf %>%
      dplyr::filter(transcript_name %in% common_transcripts) %>%
      dplyr::arrange(transcript_support_level,
                     transcript_name) %>% 
      dplyr::filter(transcript_support_level == min(transcript_support_level, na.rm = TRUE)) %>%
      head(n = 1) %>%
      pull(transcript_name)
    
    if (length(common_transcripts_best_tsl) == 1){

      up_exon_number <- up_exon_gtf$exon_number[
        up_exon_gtf$transcript_name == common_transcripts_best_tsl
      ]
      
      down_exon_number <- down_exon_gtf$exon_number[
        down_exon_gtf$transcript_name == common_transcripts_best_tsl
      ]
      
      junction_annot <- glue::glue("{common_transcripts_best_tsl}, exon{up_exon_number}:exon{down_exon_number}")
      return(junction_annot)

    }
    
    # Scenario 6: exons map to multiple common transcripts with no TSL --> take first transcripts
    common_transcripts_no_tsl_first <- up_exon_gtf %>%
      dplyr::filter(transcript_name %in% common_transcripts) %>%
      dplyr::arrange(transcript_name) %>% 
      head(n = 1) %>%
      pull(transcript_name)
    
    if (length(common_transcripts_no_tsl_first) == 1){
      
      up_exon_number <- up_exon_gtf$exon_number[
        up_exon_gtf$transcript_name == common_transcripts_no_tsl_first
      ]
      
      down_exon_number <- down_exon_gtf$exon_number[
        down_exon_gtf$transcript_name == common_transcripts_no_tsl_first
      ]
      
      junction_annot <- glue::glue("{common_transcripts_no_tsl_first}, exon{up_exon_number}:exon{down_exon_number}")
      return(junction_annot)

    }
    
  }
  
  # Scenario 7: exons do not map to any common transcripts --> take transcripts with best TSL for each exon
  if (length(common_transcripts) == 0 & nrow(up_exon_gtf) > 0 & nrow(down_exon_gtf) > 0){
    
    up_exon <- up_exon_gtf %>%
      dplyr::arrange(transcript_support_level,
                     transcript_name) %>% 
      head(n = 1) %>%
      dplyr::mutate(exon_id = glue::glue("{transcript_name}, exon{exon_number}")) %>%
      pull(exon_id)
    
    down_exon <- down_exon_gtf %>%
      dplyr::arrange(transcript_support_level,
                     transcript_name) %>% 
      head(n = 1) %>%
      dplyr::mutate(exon_id = glue::glue("{transcript_name}, exon{exon_number}")) %>%
      pull(exon_id)
    
    junction_annot <- glue::glue("{up_exon}:{down_exon}")
    return(junction_annot)

  }
  
  # Scenario 8: one exon does not map to gtf exon --> annotate as novel donor or acceptor
  if ((nrow(up_exon_gtf) == 0 & nrow(down_exon_gtf) != 0) | (nrow(up_exon_gtf) != 0 & nrow(down_exon_gtf) == 0)) {
    
    if (nrow(up_exon_gtf) == 0){
      
      down_exon <- down_exon_gtf %>%
        dplyr::arrange(transcript_support_level,
                       transcript_name) %>%
        head(n = 1) %>%
        dplyr::select(strand, transcript_name,
                      exon_number)
      
      novel_type <- ifelse(down_exon$strand == "+", 
                           "novel_donor",
                           "novel_acceptor")
      
      junction_annot <- glue::glue("{novel_type}:{down_exon$transcript_name}, exon{down_exon$exon_number}")
      return(junction_annot)

    } else if (nrow(down_exon_gtf) == 0){
      
      up_exon <- up_exon_gtf %>%
        dplyr::arrange(transcript_support_level,
                       transcript_name) %>%
        head(n = 1) %>%
        dplyr::select(strand, transcript_name,
                      exon_number)
      
      novel_type <- ifelse(up_exon$strand == "+", 
                           "novel_acceptor",
                           "novel_donor")
      
      junction_annot <- glue::glue("{novel_type}:{up_exon$transcript_name}, exon{up_exon$exon_number}")
      return(junction_annot)

    }
    
  }
  
  # Scenario 9: Neither exon maps to gtf exon
  if (nrow(up_exon_gtf) == 0 & nrow(down_exon_gtf) == 0){
    
    return("Novel exon:exon junction")
    
  }
  
}
