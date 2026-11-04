# Function to annotate splice junctions to transcripts, exons, introns

# Helper to package junction annotation with Ensembl transcript ID(s). When
# up- and downstream exons map to different transcripts, IDs are collapsed
# with ";". Missing IDs (novel exons) are dropped. `cassette_status` records
# whether the rMATS skipped exon of an ES junction is in the chosen transcript.
make_annot <- function(junction_name, transcript_ids = NA_character_,
                       cassette_status = NA_character_) {

  transcript_ids <- unique(transcript_ids[!is.na(transcript_ids)])

  tibble::tibble(
    junction_name = as.character(junction_name),
    transcript_id = if (length(transcript_ids) == 0) NA_character_ else paste(transcript_ids, collapse = ";"),
    cassette_status = as.character(cassette_status)
  )
}

# Convert GTF transcript support levels ("1"-"5", "NA") to numbers so that
# 1 (best supported) sorts first and NA sorts last.
parse_tsl <- function(tsl) {
  suppressWarnings(readr::parse_number(as.character(tsl)))
}

# MANE Select flag for GTF rows. The calling script adds a logical `is_mane`
# column to the GTF; without it no transcript is treated as MANE.
mane_flag <- function(exon_gtf) {
  if ("is_mane" %in% names(exon_gtf)) {
    exon_gtf$is_mane %in% TRUE
  } else {
    rep(FALSE, nrow(exon_gtf))
  }
}

# Order GTF exon rows for a single junction-facing exon: protein-coding
# transcripts first, then MANE Select, then best transcript support level,
# then transcript name as a deterministic fallback.
rank_exon_rows <- function(exon_gtf) {
  tibble::as_tibble(exon_gtf) %>%
    dplyr::mutate(
      pc_rank = transcript_type %in% "protein_coding",
      is_mane = mane_flag(exon_gtf),
      tsl = parse_tsl(transcript_support_level)
    ) %>%
    dplyr::arrange(dplyr::desc(pc_rank),
                   dplyr::desc(is_mane),
                   tsl,
                   transcript_name)
}

annotate_junctions <- function(chr = chr,
                               up_exon_start = up_exon_start,
                               up_exon_end = up_exon_end,
                               down_exon_start = down_exon_start,
                               down_exon_end = down_exon_end,
                               gene_symbol = gene_symbol,
                               preference = preference,
                               gtf = gtf,
                               skipped_exon_start = NA,
                               skipped_exon_end = NA) {
  
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

  # For ES junctions, find transcripts containing the rMATS skipped (cassette)
  # exon, matched on both of its boundaries within the same gene.
  has_cassette <- isTRUE(preference == "ES") &&
    !is.na(skipped_exon_start) && !is.na(skipped_exon_end)

  cassette_tx_ids <- if (has_cassette) {
    gtf_dt[
      .(chr, gene_symbol, skipped_exon_start),
      on = .(seqnames, gene_name, start),
      nomatch = 0L
    ][end == skipped_exon_end, unique(transcript_id)]
  } else {
    character(0)
  }

  # Package an annotation together with its cassette exon status
  annot <- function(junction_name, transcript_ids = NA_character_) {
    cassette_status <- if (!has_cassette) {
      NA_character_
    } else if (any(transcript_ids %in% cassette_tx_ids)) {
      "in chosen transcript"
    } else if (length(cassette_tx_ids) > 0) {
      "annotated elsewhere"
    } else {
      "novel"
    }
    make_annot(junction_name, transcript_ids, cassette_status)
  }

  # Scenario 1: junction is RI
  if (preference == "RI") {
    
    # if adjacent exon is not annotated:
    if (nrow(up_exon_gtf) == 0 & nrow(down_exon_gtf) == 0) {
      
      return(annot("RI event adjacent to novel exon"))

      # if adjacent exon is upstream and annotated:
    } else if (nrow(up_exon_gtf) > 0 & nrow(down_exon_gtf) == 0) {

      # choose protein-coding, then MANE, then best-TSL transcript
      exon_best_TSL <- rank_exon_rows(up_exon_gtf) %>%
        head(n = 1) %>%
        dplyr::select(strand, transcript_name, transcript_id,
                      exon_number)

      # get intron number from exon number and strand
      intron_number <- ifelse(exon_best_TSL$strand == "+",
                              exon_best_TSL$exon_number,
                              as.numeric(exon_best_TSL$exon_number) - 1)

      junction_annot <- glue::glue("{exon_best_TSL$transcript_name}, exon{exon_best_TSL$exon_number}:intron{intron_number}")
      return(annot(junction_annot, exon_best_TSL$transcript_id))

      # if adjacent exon is downstream and annotated:
    } else if (nrow(up_exon_gtf) == 0 & nrow(down_exon_gtf) > 0) {

      exon_best_TSL <- rank_exon_rows(down_exon_gtf) %>%
        head(n = 1) %>%
        dplyr::select(strand, transcript_name, transcript_id,
                      exon_number)

      intron_number <- ifelse(exon_best_TSL$strand == "+",
                              as.numeric(exon_best_TSL$exon_number) - 1,
                              as.numeric(exon_best_TSL$exon_number))

      junction_annot <- glue::glue("{exon_best_TSL$transcript_name}, intron{intron_number}:exon{exon_best_TSL$exon_number}")
      return(annot(junction_annot, exon_best_TSL$transcript_id))
      
    }
    
  }
  
  # For non-RI events: determine if exons are part of same transcript
  common_transcripts <- intersect(up_exon_gtf$transcript_name,
                                  down_exon_gtf$transcript_name)
  
  # Scenario 2: exons share one or more common transcripts --> rank candidates
  if (length(common_transcripts) >= 1){
    
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
        transcript_id = up_exon_gtf$transcript_id[
          up_exon_gtf$transcript_name == transcript_name
        ][1],
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
        
        # For ES calls, prefer transcript isoforms with exon(s) between the
        # two junction-facing exons.  This selects an isoform that contains
        # the skipped exon rather than one in which the exons are
        # naturally adjacent.
        n_skipped_exons = pmax(abs(down_exon_number - up_exon_number) - 1, 0),
        contains_skipped_exon = preference == "ES" & n_skipped_exons > 0,

        # For ES calls with a consensus rMATS event, prefer isoforms that
        # contain that event's skipped exon over ones skipping other exons.
        contains_cassette = transcript_id %in% cassette_tx_ids,

        # MANE Select transcript
        is_mane = transcript_name %in%
          up_exon_gtf$transcript_name[mane_flag(up_exon_gtf)],

        # ranking fields
        skip_rank = ifelse(contains_skipped_exon, 1, 0),
        seq_rank = ifelse(sequential, 1, 0),
        pc_rank = ifelse(is_protein_coding, 1, 0)
      )

    # unified ranking
    best_tx <- candidate_tbl %>%
      arrange(
        desc(contains_cassette), # ES: transcript contains the rMATS skipped exon
        desc(skip_rank),         # ES: otherwise, transcript skips some exon(s)
        desc(pc_rank),           # protein-coding first
        desc(is_mane),           # then MANE Select
        tsl,                     # then the best-supported transcript (NA last)
        desc(seq_rank),          # then flanking exons adjacent in the transcript
        transcript_name          # deterministic fallback
      ) %>%
      dplyr::slice_head(n = 1)

    return(annot(
      glue::glue(
        "{best_tx$transcript_name}, exon{best_tx$up_exon_number}:exon{best_tx$down_exon_number}"
      ),
      best_tx$transcript_id
    ))
  }
  
  # Scenario 7: exons do not map to any common transcripts --> take best
  # protein-coding/MANE/TSL-ranked transcript for each exon
  if (length(common_transcripts) == 0 & nrow(up_exon_gtf) > 0 & nrow(down_exon_gtf) > 0){

    up_exon <- rank_exon_rows(up_exon_gtf) %>%
      head(n = 1) %>%
      dplyr::mutate(exon_id = glue::glue("{transcript_name}, exon{exon_number}"))

    down_exon <- rank_exon_rows(down_exon_gtf) %>%
      head(n = 1) %>%
      dplyr::mutate(exon_id = glue::glue("{transcript_name}, exon{exon_number}"))

    junction_annot <- glue::glue("{up_exon$exon_id}:{down_exon$exon_id}")
    return(annot(junction_annot,
                 c(up_exon$transcript_id, down_exon$transcript_id)))
    
  }
  
  # Scenario 8: one exon does not map to gtf exon --> annotate as novel donor or acceptor
  if ((nrow(up_exon_gtf) == 0 & nrow(down_exon_gtf) != 0) | (nrow(up_exon_gtf) != 0 & nrow(down_exon_gtf) == 0)) {
    
    if (nrow(up_exon_gtf) == 0){
      
      down_exon <- rank_exon_rows(down_exon_gtf) %>%
        head(n = 1) %>%
        dplyr::select(strand, transcript_name, transcript_id,
                      exon_number)
      
      novel_type <- ifelse(down_exon$strand == "+", 
                           "novel_donor",
                           "novel_acceptor")
      
      junction_annot <- if (novel_type == "novel_acceptor") {
        glue::glue("{down_exon$transcript_name}, exon{down_exon$exon_number}:{novel_type}")
      } else {
        glue::glue("{novel_type}:{down_exon$transcript_name}, exon{down_exon$exon_number}")
      }
      return(annot(junction_annot, down_exon$transcript_id))

    } else if (nrow(down_exon_gtf) == 0){

      up_exon <- rank_exon_rows(up_exon_gtf) %>%
        head(n = 1) %>%
        dplyr::select(strand, transcript_name, transcript_id,
                      exon_number)
      
      novel_type <- ifelse(up_exon$strand == "+", 
                           "novel_acceptor",
                           "novel_donor")
      
      junction_annot <- if (novel_type == "novel_acceptor") {
        glue::glue("{up_exon$transcript_name}, exon{up_exon$exon_number}:{novel_type}")
      } else {
        glue::glue("{novel_type}:{up_exon$transcript_name}, exon{up_exon$exon_number}")
      }
      return(annot(junction_annot, up_exon$transcript_id))
      
    }
    
  }
  
  # Scenario 9: Neither exon maps to gtf exon
  if (nrow(up_exon_gtf) == 0 & nrow(down_exon_gtf) == 0){
    
    return(annot("Novel exon:exon junction"))
    
  }
  
}
