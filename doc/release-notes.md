# release notes


## current release (v13)
- Data release date: 2026-09-04
- OpenPedCan data release date: 2024-03-01 (v15)
- status: available

New files: 
* `depmap-splice-events-rmats.tsv.gz`
* `depmap-gene-expression-rsem-tpm-collapsed.rds`
* `depmap-rna-isoform-expression-rsem-tpm.rds`
* `depmap-rna-metadata.tsv`; depmap RNA-seq sample metadata with Bioassay IDs
* `SJ.merged.control-cohort.tsv.gz`; STAR junction counts, filtered for normal brain control cohort samples

```
v13
.
├── AlphaFold-tumor-enriched-splice-junction-predictions.zip
├── GSE73721-normal-histologies.tsv
├── brain_cell_type-rmats_merged_raw_A3SS.qs2
├── brain_cell_type-rmats_merged_raw_A5SS.qs2
├── brain_cell_type-rmats_merged_raw_RI.qs2
├── brain_cell_type-rmats_merged_raw_SE.qs2
├── brain_cell_type_input_read_counts.tsv
├── consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only_autosomes.tsv.gz
├── consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only_x_and_y.tsv.gz
├── cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz
├── cptac-protein-imputed-prot-expression-abundance.tsv.gz
├── depmap-gene-expression-rsem-tpm-collapsed.rds
├── depmap-rna-isoform-expression-rsem-tpm.rds
├── depmap-rna-metadata.tsv
├── depmap-splice-events-rmats.tsv.gz
├── evodevo-histologies.tsv
├── evodevo-rmats_merged_raw_A3SS.qs2
├── evodevo-rmats_merged_raw_A5SS.qs2
├── evodevo-rmats_merged_raw_RI.qs2
├── evodevo-rmats_merged_raw_SE.qs2
├── evodevo_gene-expression-rsem-tpm-collapsed.all.rds
├── evodevo_input_read_counts.tsv
├── evodevo_rna-isoform-expression-rsem-tpm.rds
├── gbm-protein-imputed-phospho-expression-abundance.tsv.gz
├── gbm-protein-imputed-prot-expression-abundance.tsv.gz
├── gtex-harmonized-gene-expression-rsem-tpm-collapsed.brain-under40.rds
├── gtex-harmonized-isoform-expression-rsem-tpm.rds
├── gtex-rmats_merged_raw_A3SS.qs2
├── gtex-rmats_merged_raw_A5SS.qs2
├── gtex-rmats_merged_raw_RI.qs2
├── gtex-rmats_merged_raw_SE.qs2
├── gtex_input_read_counts.tsv
├── histologies.tsv
├── hope-protein-imputed-phospho-expression-abundance.tsv.gz
├── hope-protein-imputed-prot-expression-abundance.tsv.gz
├── independent-specimens.methyl.primary-plus.eachcohort.tsv
├── independent-specimens.methyl.primary-plus.tsv
├── independent-specimens.methyl.primary.eachcohort.tsv
├── independent-specimens.methyl.primary.tsv
├── independent-specimens.methyl.relapse.eachcohort.tsv
├── independent-specimens.methyl.relapse.tsv
├── independent-specimens.rnaseq.primary-plus-pre-release.tsv
├── independent-specimens.rnaseq.primary-pre-release.tsv
├── independent-specimens.rnaseq.relapse-pre-release.tsv
├── independent-specimens.rnaseqpanel.primary-plus.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary-plus.tsv
├── independent-specimens.rnaseqpanel.primary.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary.tsv
├── independent-specimens.rnaseqpanel.relapse.eachcohort.tsv
├── independent-specimens.rnaseqpanel.relapse.tsv
├── independent-specimens.wgs.primary-plus.eachcohort.tsv
├── independent-specimens.wgs.primary-plus.tsv
├── independent-specimens.wgs.primary.eachcohort.tsv
├── independent-specimens.wgs.primary.tsv
├── independent-specimens.wgs.relapse.eachcohort.tsv
├── independent-specimens.wgs.relapse.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv
├── long-read-RNAseq.zip
├── md5sum.txt
├── normal-brain-isoform-expression-rsem-tpm.rds
├── normal-brain.gene-expression-rsem-tpm-collapsed.all.rds
├── normal_ped_brain-rmats_merged_raw_A3SS.qs2
├── normal_ped_brain-rmats_merged_raw_A5SS.qs2
├── normal_ped_brain-rmats_merged_raw_RI.qs2
├── normal_ped_brain-rmats_merged_raw_SE.qs2
├── normal_ped_brain_input_read_counts.tsv
├── pbta-gene-expr-log2-tpm-combat-corrected.qs2
├── pbta-isoform-expr-log2-tpm-combat-corrected.qs2
├── pbta-rmats_merged_raw_A3SS.qs2
├── pbta-rmats_merged_raw_A5SS.qs2
├── pbta-rmats_merged_raw_RI.qs2
├── pbta-rmats_merged_raw_SE.qs2
├── pbta_gene-expression-rsem-tpm-collapsed.rds
├── pbta_input_read_counts.tsv
├── ped-normal-brain-histologies.tsv
├── release-notes.md
├── rna-isoform-expression-rsem-tpm.rds
├── SJ.merged.control-cohort.tsv.gz
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── snv-mutation-tmb-all.tsv
├── snv-mutation-tmb-coding.tsv
├── snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz
├── tumor-enriched-oncofetal-diff-splice-junctions.qs2
├── tumor-enriched-oncofetal-junction-diff-splice-event-annotation.tsv.gz
├── unipLocCytopl.hg38.col.txt
├── unipLocExtra.hg38.col.txt
└── unipLocTransMemb.hg38.col.txt
```


## archived release (v12)
- Data release date: 2026-07-28
- OpenPedCan data release date: 2024-03-01 (v15)
- status: available

New files: 
* `AlphaFold-tumor-enriched-splice-junction-predictions.zip`; zipped file containing AlphaFold2-multimer results from proteoform & pepMLM-predicted peptide binder pairs
* `long-read-RNAseq.zip`; zipped file containing unfiltered long-read RNA-seq workflow abundance, junction, and gff files. Replaces `atrt-hgg-longread-rna-seq.zip` from v11. 
* `consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only_autosomes.tsv.gz` and `consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only_x_and_y.tsv.gz`; OpenPedCan-derived CNV files. 

```
v12
.
├── AlphaFold-tumor-enriched-splice-junction-predictions.zip
├── GSE73721-normal-histologies.tsv
├── brain_cell_type-rmats_merged_raw_A3SS.qs2
├── brain_cell_type-rmats_merged_raw_A5SS.qs2
├── brain_cell_type-rmats_merged_raw_RI.qs2
├── brain_cell_type-rmats_merged_raw_SE.qs2
├── brain_cell_type_input_read_counts.tsv
├── consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only_autosomes.tsv.gz
├── consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only_x_and_y.tsv.gz
├── cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz
├── cptac-protein-imputed-prot-expression-abundance.tsv.gz
├── evodevo-histologies.tsv
├── evodevo-rmats_merged_raw_A3SS.qs2
├── evodevo-rmats_merged_raw_A5SS.qs2
├── evodevo-rmats_merged_raw_RI.qs2
├── evodevo-rmats_merged_raw_SE.qs2
├── evodevo_gene-expression-rsem-tpm-collapsed.all.rds
├── evodevo_input_read_counts.tsv
├── evodevo_rna-isoform-expression-rsem-tpm.rds
├── gbm-protein-imputed-phospho-expression-abundance.tsv.gz
├── gbm-protein-imputed-prot-expression-abundance.tsv.gz
├── gtex-harmonized-gene-expression-rsem-tpm-collapsed.brain-under40.rds
├── gtex-harmonized-isoform-expression-rsem-tpm.rds
├── gtex-rmats_merged_raw_A3SS.qs2
├── gtex-rmats_merged_raw_A5SS.qs2
├── gtex-rmats_merged_raw_RI.qs2
├── gtex-rmats_merged_raw_SE.qs2
├── gtex_input_read_counts.tsv
├── histologies.tsv
├── hope-protein-imputed-phospho-expression-abundance.tsv.gz
├── hope-protein-imputed-prot-expression-abundance.tsv.gz
├── independent-specimens.methyl.primary-plus.eachcohort.tsv
├── independent-specimens.methyl.primary-plus.tsv
├── independent-specimens.methyl.primary.eachcohort.tsv
├── independent-specimens.methyl.primary.tsv
├── independent-specimens.methyl.relapse.eachcohort.tsv
├── independent-specimens.methyl.relapse.tsv
├── independent-specimens.rnaseq.primary-plus-pre-release.tsv
├── independent-specimens.rnaseq.primary-pre-release.tsv
├── independent-specimens.rnaseq.relapse-pre-release.tsv
├── independent-specimens.rnaseqpanel.primary-plus.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary-plus.tsv
├── independent-specimens.rnaseqpanel.primary.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary.tsv
├── independent-specimens.rnaseqpanel.relapse.eachcohort.tsv
├── independent-specimens.rnaseqpanel.relapse.tsv
├── independent-specimens.wgs.primary-plus.eachcohort.tsv
├── independent-specimens.wgs.primary-plus.tsv
├── independent-specimens.wgs.primary.eachcohort.tsv
├── independent-specimens.wgs.primary.tsv
├── independent-specimens.wgs.relapse.eachcohort.tsv
├── independent-specimens.wgs.relapse.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv
├── long-read-RNAseq.zip
├── md5sum.txt
├── normal-brain-isoform-expression-rsem-tpm.rds
├── normal-brain.gene-expression-rsem-tpm-collapsed.all.rds
├── normal_ped_brain-rmats_merged_raw_A3SS.qs2
├── normal_ped_brain-rmats_merged_raw_A5SS.qs2
├── normal_ped_brain-rmats_merged_raw_RI.qs2
├── normal_ped_brain-rmats_merged_raw_SE.qs2
├── normal_ped_brain_input_read_counts.tsv
├── pbta-gene-expr-log2-tpm-combat-corrected.qs2
├── pbta-isoform-expr-log2-tpm-combat-corrected.qs2
├── pbta-rmats_merged_raw_A3SS.qs2
├── pbta-rmats_merged_raw_A5SS.qs2
├── pbta-rmats_merged_raw_RI.qs2
├── pbta-rmats_merged_raw_SE.qs2
├── pbta_gene-expression-rsem-tpm-collapsed.rds
├── pbta_input_read_counts.tsv
├── ped-normal-brain-histologies.tsv
├── release-notes.md
├── rna-isoform-expression-rsem-tpm.rds
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── snv-mutation-tmb-all.tsv
├── snv-mutation-tmb-coding.tsv
├── snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz
├── tumor-enriched-oncofetal-diff-splice-junctions.qs2
├── tumor-enriched-oncofetal-junction-diff-splice-event-annotation.tsv.gz
├── unipLocCytopl.hg38.col.txt
├── unipLocExtra.hg38.col.txt
└── unipLocTransMemb.hg38.col.txt
```

## archived release (v11)
- Data release date: 2026-06-11
- OpenPedCan data release date: 2024-03-01 (v15)
- status: available

New files: 
* `atrt-hgg-longread-rna-seq.zip`; zipped file containing abundance, gff, isoform classification, and junction files outputted by Kinnex/MAS-Iso-Seq Complete Long-Read Pipeline. Contains subdirectories for ATRT and HGG sequencing projects along with sample metadata files. 

```
v11
.
├── atrt-hgg-longread-rna-seq.zip
├── GSE73721-normal-histologies.tsv
├── brain_cell_type-rmats_merged_raw_A3SS.qs2
├── brain_cell_type-rmats_merged_raw_A5SS.qs2
├── brain_cell_type-rmats_merged_raw_RI.qs2
├── brain_cell_type-rmats_merged_raw_SE.qs2
├── brain_cell_type_input_read_counts.tsv
├── cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz
├── cptac-protein-imputed-prot-expression-abundance.tsv.gz
├── evodevo-histologies.tsv
├── evodevo-rmats_merged_raw_A3SS.qs2
├── evodevo-rmats_merged_raw_A5SS.qs2
├── evodevo-rmats_merged_raw_RI.qs2
├── evodevo-rmats_merged_raw_SE.qs2
├── evodevo_gene-expression-rsem-tpm-collapsed.all.rds
├── evodevo_input_read_counts.tsv
├── evodevo_rna-isoform-expression-rsem-tpm.rds
├── gbm-protein-imputed-phospho-expression-abundance.tsv.gz
├── gbm-protein-imputed-prot-expression-abundance.tsv.gz
├── gtex-harmonized-gene-expression-rsem-tpm-collapsed.brain-under40.rds
├── gtex-harmonized-isoform-expression-rsem-tpm.rds
├── gtex-rmats_merged_raw_A3SS.qs2
├── gtex-rmats_merged_raw_A5SS.qs2
├── gtex-rmats_merged_raw_RI.qs2
├── gtex-rmats_merged_raw_SE.qs2
├── gtex_input_read_counts.tsv
├── histologies.tsv
├── hope-protein-imputed-phospho-expression-abundance.tsv.gz
├── hope-protein-imputed-prot-expression-abundance.tsv.gz
├── independent-specimens.methyl.primary-plus.eachcohort.tsv
├── independent-specimens.methyl.primary-plus.tsv
├── independent-specimens.methyl.primary.eachcohort.tsv
├── independent-specimens.methyl.primary.tsv
├── independent-specimens.methyl.relapse.eachcohort.tsv
├── independent-specimens.methyl.relapse.tsv
├── independent-specimens.rnaseq.primary-plus-pre-release.tsv
├── independent-specimens.rnaseq.primary-pre-release.tsv
├── independent-specimens.rnaseq.relapse-pre-release.tsv
├── independent-specimens.rnaseqpanel.primary-plus.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary-plus.tsv
├── independent-specimens.rnaseqpanel.primary.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary.tsv
├── independent-specimens.rnaseqpanel.relapse.eachcohort.tsv
├── independent-specimens.rnaseqpanel.relapse.tsv
├── independent-specimens.wgs.primary-plus.eachcohort.tsv
├── independent-specimens.wgs.primary-plus.tsv
├── independent-specimens.wgs.primary.eachcohort.tsv
├── independent-specimens.wgs.primary.tsv
├── independent-specimens.wgs.relapse.eachcohort.tsv
├── independent-specimens.wgs.relapse.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv
├── md5sum.txt
├── normal-brain-isoform-expression-rsem-tpm.rds
├── normal-brain.gene-expression-rsem-tpm-collapsed.all.rds
├── normal_ped_brain-rmats_merged_raw_A3SS.qs2
├── normal_ped_brain-rmats_merged_raw_A5SS.qs2
├── normal_ped_brain-rmats_merged_raw_RI.qs2
├── normal_ped_brain-rmats_merged_raw_SE.qs2
├── normal_ped_brain_input_read_counts.tsv
├── pbta-gene-expr-log2-tpm-combat-corrected.qs2
├── pbta-isoform-expr-log2-tpm-combat-corrected.qs2
├── pbta-rmats_merged_raw_A3SS.qs2
├── pbta-rmats_merged_raw_A5SS.qs2
├── pbta-rmats_merged_raw_RI.qs2
├── pbta-rmats_merged_raw_SE.qs2
├── pbta_gene-expression-rsem-tpm-collapsed.rds
├── pbta_input_read_counts.tsv
├── ped-normal-brain-histologies.tsv
├── release-notes.md
├── rna-isoform-expression-rsem-tpm.rds
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── snv-mutation-tmb-all.tsv
├── snv-mutation-tmb-coding.tsv
├── snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz
├── tumor-enriched-oncofetal-diff-splice-junctions.qs2
├── tumor-enriched-oncofetal-junction-diff-splice-event-annotation.tsv.gz
├── unipLocCytopl.hg38.col.txt
├── unipLocExtra.hg38.col.txt
└── unipLocTransMemb.hg38.col.txt
```

## archived release (v10)
- Data release date: 2026-04-20
- OpenPedCan data release date: 2024-03-01 (v15)
- status: available

New files: 
* `normal-brain-isoform-expression-rsem-tpm.rds`: TPM matrix including normal ped brain and cell type samples. Generated in https://cavatica.sbgenomics.com/u/harenzaj/pediatric-normal-brain-harmonization, copied from s3://bti-openaccess-us-east-1-prd-rokita-lab/data-merges/v4/
* `normal-brain.gene-expression-rsem-tpm-collapsed.all.rds`: TPM matrix including normal ped brain and cell type samples. Generated in https://cavatica.sbgenomics.com/u/harenzaj/pediatric-normal-brain-harmonization, copied from s3://bti-openaccess-us-east-1-prd-rokita-lab/data-merges/v4/
* `pbta-gene-expr-log2-tpm-combat-corrected.qs2`: generated in tumor-enriched-splicing repository, commit 1cefce467a56197f3d10a79c0897848e6a872ba8
* `pbta-isoform-expr-log2-tpm-combat-corrected.qs2`: generated in tumor-enriched-splicing repository, commit 1cefce467a56197f3d10a79c0897848e6a872ba8
* `tumor-enriched-oncofetal-diff-splice-junctions.qs2`: generated in tumor-enriched-splicing repository, commit 9464fd106ea06cbabe3d9efd0a18252d2ec0b9db
* `tumor-enriched-oncofetal-junction-diff-splice-event-annotation.tsv.gz`: generated in tumor-enriched-splicing repository, commit 9464fd106ea06cbabe3d9efd0a18252d2ec0b9db

```
v10
.
├── GSE73721-normal-histologies.tsv
├── brain_cell_type-rmats_merged_raw_A3SS.qs2
├── brain_cell_type-rmats_merged_raw_A5SS.qs2
├── brain_cell_type-rmats_merged_raw_RI.qs2
├── brain_cell_type-rmats_merged_raw_SE.qs2
├── brain_cell_type_input_read_counts.tsv
├── cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz
├── cptac-protein-imputed-prot-expression-abundance.tsv.gz
├── evodevo-histologies.tsv
├── evodevo-rmats_merged_raw_A3SS.qs2
├── evodevo-rmats_merged_raw_A5SS.qs2
├── evodevo-rmats_merged_raw_RI.qs2
├── evodevo-rmats_merged_raw_SE.qs2
├── evodevo_gene-expression-rsem-tpm-collapsed.all.rds
├── evodevo_input_read_counts.tsv
├── evodevo_rna-isoform-expression-rsem-tpm.rds
├── gbm-protein-imputed-phospho-expression-abundance.tsv.gz
├── gbm-protein-imputed-prot-expression-abundance.tsv.gz
├── gtex-harmonized-gene-expression-rsem-tpm-collapsed.brain-under40.rds
├── gtex-harmonized-isoform-expression-rsem-tpm.rds
├── gtex-rmats_merged_raw_A3SS.qs2
├── gtex-rmats_merged_raw_A5SS.qs2
├── gtex-rmats_merged_raw_RI.qs2
├── gtex-rmats_merged_raw_SE.qs2
├── gtex_input_read_counts.tsv
├── histologies.tsv
├── hope-protein-imputed-phospho-expression-abundance.tsv.gz
├── hope-protein-imputed-prot-expression-abundance.tsv.gz
├── independent-specimens.methyl.primary-plus.eachcohort.tsv
├── independent-specimens.methyl.primary-plus.tsv
├── independent-specimens.methyl.primary.eachcohort.tsv
├── independent-specimens.methyl.primary.tsv
├── independent-specimens.methyl.relapse.eachcohort.tsv
├── independent-specimens.methyl.relapse.tsv
├── independent-specimens.rnaseq.primary-plus-pre-release.tsv
├── independent-specimens.rnaseq.primary-pre-release.tsv
├── independent-specimens.rnaseq.relapse-pre-release.tsv
├── independent-specimens.rnaseqpanel.primary-plus.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary-plus.tsv
├── independent-specimens.rnaseqpanel.primary.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary.tsv
├── independent-specimens.rnaseqpanel.relapse.eachcohort.tsv
├── independent-specimens.rnaseqpanel.relapse.tsv
├── independent-specimens.wgs.primary-plus.eachcohort.tsv
├── independent-specimens.wgs.primary-plus.tsv
├── independent-specimens.wgs.primary.eachcohort.tsv
├── independent-specimens.wgs.primary.tsv
├── independent-specimens.wgs.relapse.eachcohort.tsv
├── independent-specimens.wgs.relapse.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv
├── md5sum.txt
├── normal-brain-isoform-expression-rsem-tpm.rds
├── normal-brain.gene-expression-rsem-tpm-collapsed.all.rds
├── normal_ped_brain-rmats_merged_raw_A3SS.qs2
├── normal_ped_brain-rmats_merged_raw_A5SS.qs2
├── normal_ped_brain-rmats_merged_raw_RI.qs2
├── normal_ped_brain-rmats_merged_raw_SE.qs2
├── normal_ped_brain_input_read_counts.tsv
├── pbta-gene-expr-log2-tpm-combat-corrected.qs2
├── pbta-isoform-expr-log2-tpm-combat-corrected.qs2
├── pbta-rmats_merged_raw_A3SS.qs2
├── pbta-rmats_merged_raw_A5SS.qs2
├── pbta-rmats_merged_raw_RI.qs2
├── pbta-rmats_merged_raw_SE.qs2
├── pbta_gene-expression-rsem-tpm-collapsed.rds
├── pbta_input_read_counts.tsv
├── ped-normal-brain-histologies.tsv
├── release-notes.md
├── rna-isoform-expression-rsem-tpm.rds
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── snv-mutation-tmb-all.tsv
├── snv-mutation-tmb-coding.tsv
├── snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz
├── tumor-enriched-oncofetal-diff-splice-junctions.qs2
├── tumor-enriched-oncofetal-junction-diff-splice-event-annotation.tsv.gz
├── unipLocCytopl.hg38.col.txt
├── unipLocExtra.hg38.col.txt
└── unipLocTransMemb.hg38.col.txt
```


## archived release (v9)
- Data release date: 2026-01-09
- OpenPedCan data release date: 2024-03-01 (v15)
- status: available

This data release adds rMATS-turbo v4.3 output by splicing case for PBTA, GTEx, Evo-Devo, and pediatric brain and brain cell type cohorts. Each cohort contains an `*input_read_counts.tsv` indicating the number of input reads used by rMATS-turbo.  
New files: 
`brain_cell_type-rmats_merged_raw_A3SS.qs2`
`brain_cell_type-rmats_merged_raw_A5SS.qs2`
`brain_cell_type-rmats_merged_raw_RI.qs2`
`brain_cell_type-rmats_merged_raw_SE.qs2`
`brain_cell_type_input_read_counts.tsv`
`evodevo-rmats_merged_raw_A3SS.qs2`
`evodevo-rmats_merged_raw_A5SS.qs2`
`evodevo-rmats_merged_raw_RI.qs2`
`evodevo-rmats_merged_raw_SE.qs2`
`evodevo_input_read_counts.tsv`
`gtex-rmats_merged_raw_A3SS.qs2`
`gtex-rmats_merged_raw_A5SS.qs2`
`gtex-rmats_merged_raw_RI.qs2`
`gtex-rmats_merged_raw_SE.qs2`
`gtex_input_read_counts.tsv`
`normal_ped_brain-rmats_merged_raw_A3SS.qs2`
`normal_ped_brain-rmats_merged_raw_A5SS.qs2`
`normal_ped_brain-rmats_merged_raw_RI.qs2`
`normal_ped_brain-rmats_merged_raw_SE.qs2`
`normal_ped_brain_input_read_counts.tsv`
`pbta-rmats_merged_raw_A3SS.qs2`
`pbta-rmats_merged_raw_A5SS.qs2`
`pbta-rmats_merged_raw_RI.qs2`
`pbta-rmats_merged_raw_SE.qs2`
`pbta_input_read_counts.tsv`

```
v9
.
├── GSE73721-normal-histologies.tsv
├── GSE73721-normal-rna-isoform-expression-rsem-tpm.rds
├── brain_cell_type-rmats_merged_raw_A3SS.qs2
├── brain_cell_type-rmats_merged_raw_A5SS.qs2
├── brain_cell_type-rmats_merged_raw_RI.qs2
├── brain_cell_type-rmats_merged_raw_SE.qs2
├── brain_cell_type_input_read_counts.tsv
├── cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz
├── cptac-protein-imputed-prot-expression-abundance.tsv.gz
├── evodevo-histologies.tsv
├── evodevo-rmats_merged_raw_A3SS.qs2
├── evodevo-rmats_merged_raw_A5SS.qs2
├── evodevo-rmats_merged_raw_RI.qs2
├── evodevo-rmats_merged_raw_SE.qs2
├── evodevo_gene-expression-rsem-tpm-collapsed.all.rds
├── evodevo_input_read_counts.tsv
├── evodevo_rna-isoform-expression-rsem-tpm.rds
├── gbm-protein-imputed-phospho-expression-abundance.tsv.gz
├── gbm-protein-imputed-prot-expression-abundance.tsv.gz
├── gtex-harmonized-gene-expression-rsem-tpm-collapsed.brain-under40.rds
├── gtex-harmonized-isoform-expression-rsem-tpm.rds
├── gtex-rmats_merged_raw_A3SS.qs2
├── gtex-rmats_merged_raw_A5SS.qs2
├── gtex-rmats_merged_raw_RI.qs2
├── gtex-rmats_merged_raw_SE.qs2
├── gtex_input_read_counts.tsv
├── histologies.tsv
├── hope-protein-imputed-phospho-expression-abundance.tsv.gz
├── hope-protein-imputed-prot-expression-abundance.tsv.gz
├── independent-specimens.methyl.primary-plus.eachcohort.tsv
├── independent-specimens.methyl.primary-plus.tsv
├── independent-specimens.methyl.primary.eachcohort.tsv
├── independent-specimens.methyl.primary.tsv
├── independent-specimens.methyl.relapse.eachcohort.tsv
├── independent-specimens.methyl.relapse.tsv
├── independent-specimens.rnaseq.primary-plus-pre-release.tsv
├── independent-specimens.rnaseq.primary-pre-release.tsv
├── independent-specimens.rnaseq.relapse-pre-release.tsv
├── independent-specimens.rnaseqpanel.primary-plus.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary-plus.tsv
├── independent-specimens.rnaseqpanel.primary.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary.tsv
├── independent-specimens.rnaseqpanel.relapse.eachcohort.tsv
├── independent-specimens.rnaseqpanel.relapse.tsv
├── independent-specimens.wgs.primary-plus.eachcohort.tsv
├── independent-specimens.wgs.primary-plus.tsv
├── independent-specimens.wgs.primary.eachcohort.tsv
├── independent-specimens.wgs.primary.tsv
├── independent-specimens.wgs.relapse.eachcohort.tsv
├── independent-specimens.wgs.relapse.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv
├── normal_ped_brain-rmats_merged_raw_A3SS.qs2
├── normal_ped_brain-rmats_merged_raw_A5SS.qs2
├── normal_ped_brain-rmats_merged_raw_RI.qs2
├── normal_ped_brain-rmats_merged_raw_SE.qs2
├── normal_ped_brain_input_read_counts.tsv
├── pbta-rmats_merged_raw_A3SS.qs2
├── pbta-rmats_merged_raw_A5SS.qs2
├── pbta-rmats_merged_raw_RI.qs2
├── pbta-rmats_merged_raw_SE.qs2
├── pbta_gene-expression-rsem-tpm-collapsed.rds
├── pbta_input_read_counts.tsv
├── ped-normal-brain-gene-expression-rsem-tpm.all.rds
├── ped-normal-brain-histologies.tsv
├── ped-normal-brain-isoform-expression-rsem-tpm.rds
├── release-notes.md
├── rna-isoform-expression-rsem-tpm.rds
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── snv-mutation-tmb-all.tsv
├── snv-mutation-tmb-coding.tsv
├── snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz
├── unipLocCytopl.hg38.col.txt
├── unipLocExtra.hg38.col.txt
└── unipLocTransMemb.hg38.col.txt
```

## archived release (v8)
- Data release date: 2025-04-03
- OpenPedCan data release date: 2024-03-01 (v15)
- status: available

This data release adds PBTA and control cohort rMATs file separated by splicing case (single exon [SE], retained intron [RI], and alternative splice site [AltSS]), and several additional OPC data files outlined below.

rMATs files:
- `GSE73721-normal-splice-events-rmats.AltSS.tsv.gz`
- `GSE73721-normal-splice-events-rmats.RI.tsv.gz`
- `GSE73721-normal-splice-events-rmats.SE.tsv.gz`
- `evodevo-harmonized-splice-events-rmats.AltSS.tsv.gz`
- `evodevo-harmonized-splice-events-rmats.RI.tsv.gz`
- `evodevo-harmonized-splice-events-rmats.SE.tsv.gz`
- `gtex-brain-under40-harmonized-splice-events-rmats.AltSS.tsv.gz`
- `gtex-brain-under40-harmonized-splice-events-rmats.RI.tsv.gz`
- `gtex-brain-under40-harmonized-splice-events-rmats.SE.tsv.gz`
- `ped-normal-brain-splice-events-rmats.AltSS.tsv.gz`
- `ped-normal-brain-splice-events-rmats.RI.tsv.gz`
- `ped-normal-brain-splice-events-rmats.SE.tsv.gz`
- `pbta-splice-events-rmats.AltSS.tsv.gz`
- `pbta-splice-events-rmats.RI.tsv.gz`
- `pbta-splice-events-rmats.SE.tsv.gz`

To reduce rMATs file sizes, splicing case-specific files were subsetted for the following columns:
- All rMATs files: `splicing_case`, `sample_id`, `GeneID`, `geneSymbol`, `chr`, `strand`, `IJC_SAMPLE_1`, `SJC_SAMPLE_1`, `IncFormLen`, `SkipFormLen`, `IncLevel1`
- `*SE.tsv.gz` files: `exonStart_0base`, `exonEnd`, `upstreamES`, `upstreamEE`, `downstreamES`, `downstreamEE`,
- `*RI.tsv.gz` files: `riExonStart_0base`, `riExonEnd`, `upstreamES`, `upstreamEE`, `downstreamES`, `downstreamEE`
- `*AltSS.tsv.gz` files: `longExonStart_0base`, `longExonEnd`, `shortES`, `shortEE`, `flankingES`, `flankingEE`

OpenPedCan files:
- `cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz`
- `cptac-protein-imputed-prot-expression-abundance.tsv.gz`
- `gbm-protein-imputed-phospho-expression-abundance.tsv.gz`
- `gbm-protein-imputed-prot-expression-abundance.tsv.gz`
- `hope-protein-imputed-phospho-expression-abundance.tsv.gz`
- `hope-protein-imputed-prot-expression-abundance.tsv.gz`
- `independent-specimens.methyl.primary-plus.eachcohort.tsv`
- `independent-specimens.methyl.primary-plus.tsv`
- `independent-specimens.methyl.primary.eachcohort.tsv`
- `independent-specimens.methyl.primary.tsv`
- `independent-specimens.methyl.relapse.eachcohort.tsv`
- `independent-specimens.methyl.relapse.tsv`
- `independent-specimens.rnaseq.primary-plus-pre-release.tsv`
- `independent-specimens.rnaseq.primary-pre-release.tsv`
- `independent-specimens.rnaseq.relapse-pre-release.tsv`
- `independent-specimens.rnaseqpanel.primary-plus.eachcohort.tsv`
- `independent-specimens.rnaseqpanel.primary-plus.tsv`
- `independent-specimens.rnaseqpanel.primary.eachcohort.tsv`
- `independent-specimens.rnaseqpanel.primary.tsv`
- `independent-specimens.rnaseqpanel.relapse.eachcohort.tsv`
- `independent-specimens.rnaseqpanel.relapse.tsv`
- `independent-specimens.wgs.primary-plus.eachcohort.tsv`
- `independent-specimens.wgs.primary-plus.tsv`
- `independent-specimens.wgs.primary.eachcohort.tsv`
- `independent-specimens.wgs.primary.tsv`
- `independent-specimens.wgs.relapse.eachcohort.tsv`
- `independent-specimens.wgs.relapse.tsv`
- `independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wgs.tsv`
- `independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv`
- `independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv`
- `independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv`
- `independent-specimens.wgswxspanel.primary.eachcohort.prefer.wgs.tsv`
- `independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv`
- `independent-specimens.wgswxspanel.primary.prefer.wgs.tsv`
- `independent-specimens.wgswxspanel.primary.prefer.wxs.tsv`
- `independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wgs.tsv`
- `independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv`
- `independent-specimens.wgswxspanel.relapse.prefer.wgs.tsv`
- `independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv`
- `snv-consensus-plus-hotspots.maf.tsv.gz` (subsetted for PBTA)
- `snv-mutation-tmb-all.tsv`
- `snv-mutation-tmb-coding.tsv`
- `snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz`


```
v8
.
├── GSE73721-normal-histologies.tsv
├── GSE73721-normal-rna-isoform-expression-rsem-tpm.rds
├── GSE73721-normal-splice-events-rmats.AltSS.tsv.gz
├── GSE73721-normal-splice-events-rmats.RI.tsv.gz
├── GSE73721-normal-splice-events-rmats.SE.tsv.gz
├── cptac-protein-imputed-phospho-expression-log2-ratio.tsv.gz
├── cptac-protein-imputed-prot-expression-abundance.tsv.gz
├── evodevo-harmonized-splice-events-rmats.AltSS.tsv.gz
├── evodevo-harmonized-splice-events-rmats.RI.tsv.gz
├── evodevo-harmonized-splice-events-rmats.SE.tsv.gz
├── evodevo-histologies.tsv
├── evodevo_gene-expression-rsem-tpm-collapsed.all.rds
├── evodevo_rna-isoform-expression-rsem-tpm.rds
├── gbm-protein-imputed-phospho-expression-abundance.tsv.gz
├── gbm-protein-imputed-prot-expression-abundance.tsv.gz
├── gtex-brain-under40-harmonized-splice-events-rmats.AltSS.tsv.gz
├── gtex-brain-under40-harmonized-splice-events-rmats.RI.tsv.gz
├── gtex-brain-under40-harmonized-splice-events-rmats.SE.tsv.gz
├── gtex-harmonized-gene-expression-rsem-tpm-collapsed.brain-under40.rds
├── gtex-harmonized-isoform-expression-rsem-tpm.rds
├── histologies.tsv
├── hope-protein-imputed-phospho-expression-abundance.tsv.gz
├── hope-protein-imputed-prot-expression-abundance.tsv.gz
├── independent-specimens.methyl.primary-plus.eachcohort.tsv
├── independent-specimens.methyl.primary-plus.tsv
├── independent-specimens.methyl.primary.eachcohort.tsv
├── independent-specimens.methyl.primary.tsv
├── independent-specimens.methyl.relapse.eachcohort.tsv
├── independent-specimens.methyl.relapse.tsv
├── independent-specimens.rnaseq.primary-plus-pre-release.tsv
├── independent-specimens.rnaseq.primary-pre-release.tsv
├── independent-specimens.rnaseq.relapse-pre-release.tsv
├── independent-specimens.rnaseqpanel.primary-plus.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary-plus.tsv
├── independent-specimens.rnaseqpanel.primary.eachcohort.tsv
├── independent-specimens.rnaseqpanel.primary.tsv
├── independent-specimens.rnaseqpanel.relapse.eachcohort.tsv
├── independent-specimens.rnaseqpanel.relapse.tsv
├── independent-specimens.wgs.primary-plus.eachcohort.tsv
├── independent-specimens.wgs.primary-plus.tsv
├── independent-specimens.wgs.primary.eachcohort.tsv
├── independent-specimens.wgs.primary.tsv
├── independent-specimens.wgs.relapse.eachcohort.tsv
├── independent-specimens.wgs.relapse.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary-plus.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.primary.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wgs.tsv
├── independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv
├── md5sum.txt
├── pbta-splice-events-rmats.AltSS.tsv.gz
├── pbta-splice-events-rmats.RI.tsv.gz
├── pbta-splice-events-rmats.SE.tsv.gz
├── pbta_gene-expression-rsem-tpm-collapsed.rds
├── ped-normal-brain-gene-expression-rsem-tpm.all.rds
├── ped-normal-brain-histologies.tsv
├── ped-normal-brain-isoform-expression-rsem-tpm.rds
├── ped-normal-brain-splice-events-rmats.AltSS.tsv.gz
├── ped-normal-brain-splice-events-rmats.RI.tsv.gz
├── ped-normal-brain-splice-events-rmats.SE.tsv.gz
├── release-notes.md
├── rna-isoform-expression-rsem-tpm.rds
├── snv-consensus-plus-hotspots.maf.tsv.gz
├── snv-mutation-tmb-all.tsv
├── snv-mutation-tmb-coding.tsv
├── snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz
├── unipLocCytopl.hg38.col.txt
├── unipLocExtra.hg38.col.txt
└── unipLocTransMemb.hg38.col.txt
```


## archived release (v7)
- Data release date: 2025-02-04
- OpenPedCan data release date: 2024-03-01 (v15)
- status: available

This data release adds an additional normal pons sample to all ped brain normal data files, and renames `GSE243682*` files with `ped-normal-brain*` prefix.

Additional files: 
- ped-normal-brain-gene-expression-rsem-tpm.all.rds
- ped-normal-brain-histologies.tsv
- ped-normal-brain-isoform-expression-rsem-tpm.rds
- ped-normal-brain-splice-events-rmats.tsv

```
v7
├── GSE73721-normal-histologies.tsv
├── GSE73721-normal-rna-isoform-expression-rsem-tpm.rds
├── GSE73721-normal-splice-events-rmats.tsv.gz
├── evodevo-harmonized.splice-events-rmats.SE.tsv.gz
├── evodevo-histologies.tsv
├── evodevo_gene-expression-rsem-tpm-collapsed.all.rds
├── evodevo_rna-isoform-expression-rsem-tpm.rds
├── gtex-brain-under40-harmonized-splice-events-rmats.SE.tsv.gz
├── gtex-harmonized-gene-expression-rsem-tpm-collapsed.brain-under40.rds
├── gtex-harmonized-isoform-expression-rsem-tpm.rds
├── gtex-harmonized-splice-events-rmats.SE.tsv.gz
├── histologies.tsv
├── md5sum.txt
├── pbta-splice-events-rmats.SE.tsv.gz
├── pbta_gene-expression-rsem-tpm-collapsed.rds
├── ped-normal-brain-gene-expression-rsem-tpm.all.rds
├── ped-normal-brain-isoform-expression-rsem-tpm.rds
├── ped-normal-brain-splice-events-rmats.tsv
├── ped-normal-brain-histologies.tsv
├── release-notes.md
├── rna-isoform-expression-rsem-tpm.rds
├── unipLocCytopl.hg38.col.txt
├── unipLocExtra.hg38.col.txt
└── unipLocTransMemb.hg38.col.txt
```

## archived release (v6)
- Data release data: 2025-01-29
- OpenPedCan data release date: 2024-03-01 (v15)
- status: available

Additional files: 
- GSE73721-normal-histologies.tsv
- GSE73721-normal-rna-isoform-expression-rsem-tpm.rds
- GSE73721-normal-splice-events-rmats.tsv.gz

```
v6
├── GSE243682_normal_gene-expression-rsem-tpm-collapsed.rds
├── GSE243682_normal_rna-isoform-expression-rsem-tpm.rds
├── GSE243682_normal_splice-events-rmats.tsv.gz
├── GSE73721-normal-histologies.tsv
├── GSE73721-normal-rna-isoform-expression-rsem-tpm.rds
├── GSE73721-normal-splice-events-rmats.tsv.gz
├── evodevo-harmonized.splice-events-rmats.SE.tsv.gz
├── evodevo-histologies.tsv
├── evodevo_gene-expression-rsem-tpm-collapsed.all.rds
├── evodevo_rna-isoform-expression-rsem-tpm.rds
├── gtex-brain-under40-harmonized-splice-events-rmats.SE.tsv.gz
├── gtex-harmonized-gene-expression-rsem-tpm-collapsed.brain-under40.rds
├── gtex-harmonized-isoform-expression-rsem-tpm.rds
├── gtex-harmonized-splice-events-rmats.SE.tsv.gz
├── histologies.tsv
├── md5sum.txt
├── pbta-splice-events-rmats.SE.tsv.gz
├── pbta_gene-expression-rsem-tpm-collapsed.rds
├── ped-normal-brain-histologies.tsv
├── release-notes.md
├── rna-isoform-expression-rsem-tpm.rds
├── unipLocCytopl.hg38.col.txt
├── unipLocExtra.hg38.col.txt
└── unipLocTransMemb.hg38.col.txt
```

## archived release (v5)
- Data release data: 2024-09-19
- OpenPedCan data release date: 2024-03-01 (v15)
- status: available

Additional files:
- evodevo-harmonized.splice-events-rmats.SE.tsv.gz
- gtex-brain-under40-harmonized-splice-events-rmats.SE.tsv.gz
- gtex-harmonized-gene-expression-rsem-tpm-collapsed.brain-under40.rds
- gtex-harmonized-splice-events-rmats.SE.tsv.gz
- md5sum-subset.txt
- pbta-splice-events-rmats.SE.tsv.gz
- pbta_gene-expression-rsem-tpm-collapsed.rds

```
v5
├── GSE243682_normal_gene-expression-rsem-tpm-collapsed.rds
├── GSE243682_normal_rna-isoform-expression-rsem-tpm.rds
├── GSE243682_normal_splice-events-rmats.tsv.gz
├── evodevo-harmonized.splice-events-rmats.SE.tsv.gz
├── evodevo-histologies.tsv
├── evodevo_gene-expression-rsem-tpm-collapsed.all.rds
├── evodevo_rna-isoform-expression-rsem-tpm.rds
├── gtex-brain-under40-harmonized-splice-events-rmats.SE.tsv.gz
├── gtex-harmonized-gene-expression-rsem-tpm-collapsed.brain-under40.rds
├── gtex-harmonized-isoform-expression-rsem-tpm.rds
├── gtex-harmonized-splice-events-rmats.SE.tsv.gz
├── histologies.tsv
├── md5sum.txt
├── pbta-splice-events-rmats.SE.tsv.gz
├── pbta_gene-expression-rsem-tpm-collapsed.rds
├── ped-normal-brain-histologies.tsv
├── release-notes.md
├── rna-isoform-expression-rsem-tpm.rds
├── unipLocCytopl.hg38.col.txt
├── unipLocExtra.hg38.col.txt
└── unipLocTransMemb.hg38.col.txt
```

## archived release (v4)
- Data release data: 2024-09-11
- OpenPedCan data release date: 2024-03-01 (v15)
- status: available

Additional files:
- `GSE243682_normal_gene-counts-rsem-expected_count-collapsed.rds`
- `GSE243682_normal_gene-expression-rsem-tpm-collapsed.rds`
- `GSE243682_normal_rna-isoform-expression-rsem-expected-counts.rds`
- `GSE243682_normal_rna-isoform-expression-rsem-tpm.rds`
- `GSE243682_normal_splice-events-rmats.tsv.gz`
- `ped-normal-brain-histologies.tsv`

```
v4
├── GSE243682_normal_gene-counts-rsem-expected_count-collapsed.rds
├── GSE243682_normal_gene-expression-rsem-tpm-collapsed.rds
├── GSE243682_normal_rna-isoform-expression-rsem-expected-counts.rds
├── GSE243682_normal_rna-isoform-expression-rsem-tpm.rds
├── GSE243682_normal_splice-events-rmats.tsv.gz
├── evodevo_gene-counts-rsem-expected_count-collapsed.all.rds
├── evodevo_gene-expression-rsem-tpm-collapsed.all.rds
├── evodevo_rna-isoform-expression-rsem-expected_count.rds
├── evodevo_rna-isoform-expression-rsem-tpm.rds
├── gtex-harmonized-gene-counts-rsem-expected_count-collapsed.all.rds
├── gtex-harmonized-gene-expression-rsem-tpm-collapsed.all.rds
├── gtex-harmonized-isoform-expression-rsem-expected_count.rds
├── gtex-harmonized-isoform-expression-rsem-tpm.rds
├── histologies.tsv
├── md5sum.txt
├── ped-normal-brain-histologies.tsv
├── release-notes.md
├── rna-isoform-expression-rsem-tpm.rds
├── unipLocCytopl.hg38.col.txt
├── unipLocExtra.hg38.col.txt
└── unipLocTransMemb.hg38.col.txt
```

## archived release (v3)
- Data release data: 2024-06-24
- OpenPedCan data release date: 2024-03-01 (v15)
- status: available

Files removed and replaced below (see [ticket](https://github.com/d3b-center/splicing-neoepitopes/issues/31))
- `pedbrain_gene-counts-rsem-expected_count-collapsed.all.rds`
- `pedbrain_gene-expression-rsem-tpm-collapsed.all.rds`
- `pedbrain_rna-isoform-expression-rsem-expected_count.rds`
- `pedbrain_rna-isoform-expression-rsem-tpm.rds`
- `pedbrain_rna-isoform-expression-rsem-tpm.rds`


```
v3
├── ctrls.filtered.SE.MATS.JC.txt
├── evodevo_gene-counts-rsem-expected_count-collapsed.all.rds
├── evodevo_gene-expression-rsem-tpm-collapsed.all.rds
├── evodevo_rna-isoform-expression-rsem-expected_count.rds
├── evodevo_rna-isoform-expression-rsem-tpm.rds
├── evodevo_splice-events-rmats.SE.tsv.gz
├── gene-expression-rsem-tpm-collapsed.rds
├── gtex-brain-under40-harmonized-splice-events-rmats.SE.tsv.gz
├── gtex-harmonized-gene-counts-rsem-expected_count-collapsed.all.rds
├── gtex-harmonized-gene-expression-rsem-tpm-collapsed.all.rds
├── gtex-harmonized-isoform-expression-rsem-expected_count.rds
├── gtex-harmonized-isoform-expression-rsem-tpm.rds
├── gtex-harmonized-splice-events-rmats.SE.tsv.gz
├── histologies.tsv
├── md5sum.txt
├── release-notes.md
├── rna-isoform-expression-rsem-tpm.rds
├── splice-events-rmats.tsv.gz
├── unipLocCytopl.hg38.col.txt
├── unipLocExtra.hg38.col.txt
└── unipLocTransMemb.hg38.col.txt
```

## archived release (v2)
- Data release data: 2024-06-24
- OpenPedCan data release date: 2024-03-01 (v15)
- status: available

Files removed and replaced below (see [ticket](https://github.com/d3b-center/splicing-neoepitopes/issues/31))
- `control-rna-isoform-expression-rsem-counts-tpm.rds`
- `ctrls.filtered.SE.MATS.JC.txt`

```
v2
├── ctrls.filtered.SE.MATS.JC.txt
├── evodevo_gene-counts-rsem-expected_count-collapsed.all.rds
├── evodevo_gene-expression-rsem-tpm-collapsed.all.rds
├── evodevo_rna-isoform-expression-rsem-expected_count.rds
├── evodevo_rna-isoform-expression-rsem-tpm.rds
├── evodevo_splice-events-rmats.tsv.gz
├── gene-expression-rsem-tpm-collapsed.rds
├── histologies.tsv
├── md5sum.txt
├── pedbrain_gene-counts-rsem-expected_count-collapsed.all.rds
├── pedbrain_gene-expression-rsem-tpm-collapsed.all.rds
├── pedbrain_rna-isoform-expression-rsem-expected_count.rds
├── pedbrain_rna-isoform-expression-rsem-tpm.rds
├── pedbrain_splice-events-rmats.tsv.gz
├── release-notes.md
├── rna-isoform-expression-rsem-tpm.rds
├── splice-events-rmats.tsv.gz
├── unipLocCytopl.hg38.col.txt
├── unipLocExtra.hg38.col.txt
└── unipLocTransMemb.hg38.col.txt
```

## archived release (v1)
- Data release data: 2024-05-02
- OpenPedCan data release date: 2024-03-01 (v15)
- status: available

Additional files:
- `control-rna-isoform-expression-rsem-counts-tpm.rds` GENCODE v39 processed controls
- `gene-expression-rsem-tpm-collapsed.rds` subsetted for PBTA samples
- `rna-isoform-expression-rsem-tpm.rds` subsetted for PBTA samples
- `splice-events-rmats-pbta.tsv.gz`: rmats subsetted by PBTA cohort
- `unipLocCytopl.hg38.col.txt`: uniprot cytoplasm protein domain annotation in bed format
- `unipLocExtra.hg38.col.txt`: uniprot extracellular protein domain annotation in bed format
- `unipLocTransMemb.hg38.col.txt`: uniprot transmembrane protein domain annotation in bed format

```
v1
    ├── control-rna-isoform-expression-rsem-counts-tpm.rds
    ├── ctrls.filtered.SE.MATS.JC.txt
    ├── gene-expression-rsem-tpm-collapsed.rds
    ├── histologies.tsv
    ├── rna-isoform-expression-rsem-tpm.rds
    ├── splice-events-rmats.tsv.gz
    ├── unipLocCytopl.hg38.col.txt
    ├── unipLocExtra.hg38.col.txt
    └── unipLocTransMemb.hg38.col.txt
```
