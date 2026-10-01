# Function to annotate splice junctions to transcripts, exons, introns

annotate_junctions <- function(chr = chr,
                               up_exon_start = up_exon_start,
                               up_exon_end = up_exon_end,
                               down_exon_start = down_exon_start,
                               down_exon_end = down_exon_end,
                               gene_symbol = gene_symbol,
                               preference = preference,
                               gtf = gtf) {
  
  # Reuse the indexed GTF data table supplied by the calling script. This
  # avoids copying and sorting the same GTF subset for every junction.
  gtf_dt <- if (data.table::is.data.table(gtf)) {
    gtf
  } else {
    data.table::as.data.table(gtf)
  }
  
  # Match the junction-facing splice boundaries, rather than requiring each
  # complete rMATS interval to equal an annotated exon. The distal boundary of
  # an rMATS interval can differ from the GTF exon boundary without creating a
  # novel splice site at this junction.
  up_exon_gtf <- gtf_dt[
    .(chr, gene_symbol, up_exon_end),
    on = .(seqnames, gene_name, end),
    nomatch = 0L
  ]
  
  down_exon_gtf <- gtf_dt[
    .(chr, gene_symbol, down_exon_start),
    on = .(seqnames, gene_name, start),
    nomatch = 0L
  ]
  
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
                              as.numeric(exon_best_TSL$exon_number) - 1,
                              as.numeric(exon_best_TSL$exon_number))
      
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
    
    #   # Scenario 3: exons each map to single protein-coding transcript
    #   common_pc_transcripts <- intersect(up_exon_gtf$transcript_name[up_exon_gtf$transcript_type == "protein_coding"],
    #                                      down_exon_gtf$transcript_name[down_exon_gtf$transcript_type == "protein_coding"])
    #   
    #   if(length(common_pc_transcripts) == 1){
    #     
    #     # filter exon numbers for those in single protein-coding transcript
    #     up_exon_number <- up_exon_gtf$exon_number[
    #       up_exon_gtf$transcript_name == common_pc_transcripts
    #     ]
    #     
    #     down_exon_number <- down_exon_gtf$exon_number[
    #       down_exon_gtf$transcript_name == common_pc_transcripts
    #     ]
    #     
    #     junction_annot <- glue::glue("{common_pc_transcripts}, exon{up_exon_number}:exon{down_exon_number}")
    #     return(junction_annot)
    #     
    #   }
    #   
    #   # Scenario 5: exons map to multiple common transcripts --> take transcript with highest TSL (1 = most supported)
    #   common_transcripts_best_tsl <- up_exon_gtf %>%
    #     dplyr::filter(transcript_name %in% common_transcripts) %>%
    #     dplyr::arrange(transcript_support_level,
    #                    transcript_name) %>% 
    #     dplyr::filter(transcript_support_level == min(transcript_support_level, na.rm = TRUE)) %>%
    #     head(n = 1) %>%
    #     pull(transcript_name)
    #   
    #   if (length(common_transcripts_best_tsl) == 1){
    #     
    #     up_exon_number <- up_exon_gtf$exon_number[
    #       up_exon_gtf$transcript_name == common_transcripts_best_tsl
    #     ]
    #     
    #     down_exon_number <- down_exon_gtf$exon_number[
    #       down_exon_gtf$transcript_name == common_transcripts_best_tsl
    #     ]
    #     
    #     junction_annot <- glue::glue("{common_transcripts_best_tsl}, exon{up_exon_number}:exon{down_exon_number}")
    #     return(junction_annot)
    #     
    #   }
    #   
    #   # Scenario 6: exons map to multiple common transcripts with no TSL --> take first transcripts
    #   common_transcripts_no_tsl_first <- up_exon_gtf %>%
    #     dplyr::filter(transcript_name %in% common_transcripts) %>%
    #     dplyr::arrange(transcript_name) %>% 
    #     head(n = 1) %>%
    #     pull(transcript_name)
    #   
    #   if (length(common_transcripts_no_tsl_first) == 1){
    #     
    #     up_exon_number <- up_exon_gtf$exon_number[
    #       up_exon_gtf$transcript_name == common_transcripts_no_tsl_first
    #     ]
    #     
    #     down_exon_number <- down_exon_gtf$exon_number[
    #       down_exon_gtf$transcript_name == common_transcripts_no_tsl_first
    #     ]
    #     
    #     junction_annot <- glue::glue("{common_transcripts_no_tsl_first}, exon{up_exon_number}:exon{down_exon_number}")
    #     return(junction_annot)
    #     
    #   }
    #   
    # }
    
    
    candidate_tbl <- tibble(transcript_name = common_transcripts) %>%
      rowwise() %>%
      mutate(
        up_exon_number = up_exon_gtf$exon_number[
          up_exon_gtf$transcript_name == transcript_name
        ][1],
        down_exon_number = down_exon_gtf$exon_number[
          down_exon_gtf$transcript_name == transcript_name
        ][1],
        strand = up_exon_gtf$strand[
          up_exon_gtf$transcript_name == transcript_name
        ][1],
        tsl = up_exon_gtf$transcript_support_level[
          up_exon_gtf$transcript_name == transcript_name
        ][1]
      ) %>%
      ungroup() %>%
      mutate(
        up_exon_number = as.numeric(up_exon_number),
        down_exon_number = as.numeric(down_exon_number),
        tsl = suppressWarnings(as.numeric(tsl)),
        
        # sequential logic
        sequential = case_when(
          strand == "+" ~ down_exon_number == up_exon_number + 1,
          strand == "-" ~ down_exon_number == up_exon_number - 1,
          TRUE ~ FALSE
        ),
        
        # protein coding (intersection of both exons)
        is_protein_coding = transcript_name %in% intersect(
          up_exon_gtf$transcript_name[up_exon_gtf$transcript_type == "protein_coding"],
          down_exon_gtf$transcript_name[down_exon_gtf$transcript_type == "protein_coding"]
        ),
        
        # ranking fields
        seq_rank = ifelse(sequential, 1, 0),
        pc_rank = ifelse(is_protein_coding, 1, 0)
      )
    
    # unified ranking
    best_tx <- candidate_tbl %>%
      arrange(
        desc(seq_rank),   # sequential first
        desc(pc_rank),    # protein-coding next
        tsl,              # lowest TSL best
        transcript_name   # deterministic fallback
      ) %>%
      dplyr::slice_head(n = 1)
    
    return(glue::glue(
      "{best_tx$transcript_name}, exon{best_tx$up_exon_number}:exon{best_tx$down_exon_number}"
    ))
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
      
      junction_annot <- if (novel_type == "novel_acceptor") {
        glue::glue("{down_exon$transcript_name}, exon{down_exon$exon_number}:{novel_type}")
      } else {
        glue::glue("{novel_type}:{down_exon$transcript_name}, exon{down_exon$exon_number}")
      }
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
      
      junction_annot <- if (novel_type == "novel_acceptor") {
        glue::glue("{up_exon$transcript_name}, exon{up_exon$exon_number}:{novel_type}")
      } else {
        glue::glue("{novel_type}:{up_exon$transcript_name}, exon{up_exon$exon_number}")
      }
      return(junction_annot)
      
    }
    
  }
  
  # Scenario 9: Neither exon maps to gtf exon
  if (nrow(up_exon_gtf) == 0 & nrow(down_exon_gtf) == 0){
    
    return("Novel exon:exon junction")
    
  }
  
}
