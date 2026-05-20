# Create cohort histologies file
#
# Ryan Corbett
#
# Feb 2026

# load packages
library(tidyverse)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "00-create-cohort-histologies")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

# Set file paths
pbta_a5ss_rmats_file <- file.path(data_dir,
                                "pbta-rmats_merged_raw_A5SS.qs2")

independent_specimens_file <- file.path(data_dir,
                                        "independent-specimens.rnaseq.primary-plus-pre-release.tsv")

plot_mapping_file <- file.path(input_dir,
                               "plot-mapping.tsv")

hist_file <- file.path(data_dir,
                       "histologies.tsv")

samples_to_rm_file <- file.path(root_dir, "analyses",
                                "02-tumor-enriched-splicing",
                                "input",
                                "pbta-rna-high-intron-samples.tsv")

ancestry_file <- file.path(input_dir, 
                           "somalier-ancestry-prediction-superpopulation.tsv")

survival_file <- file.path(input_dir,
                           "openpedcan_histologies0311.csv")

# Wrangle data 

# pull PBTA samples from one of the rMATS files
pbta_samples <- qs2::qs_read(pbta_a5ss_rmats_file) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "_[^_]+$")) %>%
  distinct(sample_id) %>%
  pull(sample_id)

samples_to_rm <- read_tsv(samples_to_rm_file) %>%
  pull(Kids_First_Biospecimen_ID)

hist <- read_tsv(hist_file)

# filter for initial/primary tumor independent specimens
independent_specimens <- read_tsv(independent_specimens_file) %>%
  dplyr::filter(tumor_descriptor %in% c("Initial CNS Tumor",
                                        "Primary Tumor")) %>%
  pull(Kids_First_Biospecimen_ID)

# Load plot mapping df 
plot_mapping_df <- read_tsv(plot_mapping_file)

### Merge OPC histologies data, plot group assignments

# Create cohort histologies 
cohort_hist <- hist %>% 
  # filter for PBTA
  dplyr::filter(Kids_First_Biospecimen_ID %in% pbta_samples) %>%
  # select relevant columns
  dplyr::select(Kids_First_Biospecimen_ID,
                Kids_First_Participant_ID,
                cohort,
                sub_cohort,
                cohort_participant_id,
                match_id,
                sample_id,
                sample_type,
                tumor_descriptor,
                composition,
                cell_line_composition,
                broad_histology,
                cancer_group,
                molecular_subtype,
                molecular_subtype_methyl,
                pathology_diagnosis,
                pathology_free_text_diagnosis,
                primary_site,
                CNS_region,
                extent_of_tumor_resection,
                RNA_library,
                age_at_diagnosis_days,
                race,
                ethnicity,
                reported_gender,
                germline_sex_estimate,
                cancer_predispositions) %>%
  # add plot group & hex codes
  left_join(plot_mapping_df %>%
              dplyr::select(broad_histology,
                            cancer_group,
                            plot_group,
                            plot_group_hex)) %>%
  # resolve NA plot group assignments
  dplyr::mutate(plot_group = case_when(
    is.na(plot_group) ~ "Other tumor",
    TRUE ~ plot_group
  )) %>%
  dplyr::mutate(plot_group = case_when(
    Kids_First_Biospecimen_ID == "BS_N7VQ1GQB" ~ "Low-grade glioma",
    TRUE ~ plot_group
  )) %>%
  dplyr::mutate(molecular_subtype = case_when(
    Kids_First_Biospecimen_ID == "BS_N7VQ1GQB" ~ "LGG, KIAA1549-BRAF",
    TRUE ~ molecular_subtype
  )) %>%
  # update plot group hex codes
  dplyr::mutate(plot_group_hex = case_when(
    plot_group == "Oligodendroglioma" ~ "tan",
    Kids_First_Biospecimen_ID == "BS_N7VQ1GQB" ~ "#8f8fbf",
    is.na(plot_group_hex) ~ "#b5b5b5",
    TRUE ~ plot_group_hex
  )) %>%
  # add missing age dx for following patients
  mutate(age_at_diagnosis_days = case_when(Kids_First_Participant_ID == "PT_AEDWCP8Z" ~
                                             as.integer(365.25*17),
                                           Kids_First_Participant_ID == "PT_6PYNBA9C" ~ as.integer(365.25*5/12),
                                           TRUE ~ age_at_diagnosis_days)
         )%>%
  # filter out exome capture libraries, patients > 40, normal samples
  dplyr::filter(!RNA_library %in% c("exome capture",
                                    "exome_capture"),
                (age_at_diagnosis_days < 365.25*40 | sub_cohort == "PNOC"),
                sample_type != "Normal",
                # rm nbl, metastatic tumors
                !cancer_group %in% c("Neuroblastoma", "Metastatic secondary tumors"),
                # rm low quality PNOC samples
                !Kids_First_Biospecimen_ID %in% samples_to_rm
                ) %>%
  # indicate if independent primary tumor 
  dplyr::mutate(is_independent_primary = case_when(
    Kids_First_Biospecimen_ID %in% independent_specimens ~ "Yes",
    TRUE ~ "No"
  ))

### Append genetic ancestry data

# Load genetic ancestry df, filter for pts in cohort, and select distinct samples per pt
ancestry_df <- read_tsv(ancestry_file) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_normal) %>%
  left_join(hist %>% dplyr::select(Kids_First_Biospecimen_ID,
                                   Kids_First_Participant_ID)) %>%
  dplyr::filter(Kids_First_Participant_ID %in% cohort_hist$Kids_First_Participant_ID,
                # one pt has two normals with different superpop calls; we will take the sample used in germline project and remove the other
                Kids_First_Biospecimen_ID != "BS_5MNX1W83") %>%
  distinct(Kids_First_Participant_ID, .keep_all = TRUE)

# Save filtered, comprehensive ancestry results to output
write_tsv(ancestry_df,
          file.path(results_dir, 
                    "cohort-somalier-genetic-ancestry-prediction.tsv"))

# append genetic ancestry superpopulation to cohort hist
cohort_hist <- cohort_hist %>%
  left_join(ancestry_df %>%
              dplyr::select(Kids_First_Participant_ID,
                            predicted_ancestry)) %>%
  distinct()


### Append updated survival data

# Load survival df and select relevant columns
survival_df <- read_csv(survival_file) %>%
  dplyr::select(Kids_First_Participant_ID,
                contains(c("EFS", "OS_days", "OS_status"))) %>%
  distinct()

# append survival data to cohort hist and define EFS_status
cohort_hist <- cohort_hist %>%
  left_join(survival_df) %>%
  # when disease-related event reported, report status as EVENT
  dplyr::mutate(EFS_status = case_when(
    EFS_event_type %in% c("Deceased-due to disease",
                          "Progressive",
                          "Progressive - Metastatic",
                          "Recurrence",
                          "Recurrence - Metastatic",
                          "Second Malignancy",
                          "Second Malignancy - Metastatic") ~ "EVENT",
    TRUE ~ "NO EVENT"
  ))
              
# write cohort hist to output
write_tsv(cohort_hist,
          file.path(results_dir,
                    "cohort-histologies.tsv"))

# print session info
sessionInfo()
