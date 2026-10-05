# Classify tumor-enriched splice junctions

Module authors: Ryan Corbett (@rjcorb)

This module classifies PBTA splice junctions as tumor-enriched and oncofetal, then annotates them with protein-domain and gene-expression information.

## Usage
### script to run analysis

```
bash run_module.sh
```

Scripts are run in numerical order; each uses results from the one before it. Results are written to `results/` (not tracked in git).

## Tumor-enriched junction classification strategy

Junction names have the form `chr:upstreamExonStart-upstreamExonEnd_downstreamExonStart-downstreamExonEnd`. The **boundary** of a junction is the intron it spans, `upstreamExonEnd_downstreamExonStart`; junctions in the same sample that share a boundary are evaluated together below.

### 1. Compare PBTA junctions with normal-brain references (`01-classify-tumor-enriched-junctions.R`)

Reference groups are GTEx brain, postnatal Evo-Devo brain (collapsed to forebrain and hindbrain groups), and normal pediatric brain.

- A junction is eligible for this comparison only when its mean CPM is below 10 in every reference group (junctions with mean CPM >= 10 in any group are dropped).
- For each PBTA sample and junction, calculate:
  - fold change (`min_cpm_fc_all`) relative to the highest reference-group mean CPM (a pseudocount of 1e-5 is added to the denominator);
  - the minimum signal-to-noise ratio (`min_cpm_snr_all`) across reference groups, where SNR is `(tumor CPM − reference mean CPM) / reference CPM SD`. Reference groups with an undefined SD (e.g. a single sample) are skipped; if no SNR can be computed, the SNR is missing and the junction is not called.
- A junction is called **tumor-enriched** in a sample when its fold change is greater than 5 and its minimum SNR is greater than 5.
- A call is retained only when no junction sharing the same boundary (intron) is classified as non-specific in the same sample.

Sample-level calls from this comparison are labeled `criteria` = `enriched expr vs. ctrls`.

### 2. Junctions absent from the reference matrices

Junctions that are not found in any reference matrix are also classified as tumor-enriched (`criteria` = `No ctrl expr, novel SS usage`) when all of the following hold:

- the junction is not a zero-length retained-intron junction (identical upstream and downstream boundary coordinates);
- its boundary is not shared with any junction present in the reference matrices;
- it uses at least one unannotated splice site, compared with GENCODE v39 exon boundaries (upstream boundary matched to exon ends, downstream boundary to exon starts); and
- it occurs in fewer than 25% of the PBTA samples surveyed.

Before this classification, these candidates are cross-checked against the control STAR junction file (`SJ.merged.control-cohort.tsv.gz`) and excluded when chromosome and intron coordinates match exactly: the STAR intron start is the upstream exon boundary plus one, and the intron end is the downstream exon boundary minus one (see `util/other-functions.R`).

### 3. Classify tumor-enriched junctions as oncofetal (`02-classify-oncofetal-junctions.R`)

All tumor-enriched junctions from both routes above are compared using Evo-Devo expression, with prenatal groups (by brain region and prenatal week bin) compared against every postnatal group (by brain region and stage).

- Fold change is `prenatal CPM / postnatal CPM` (pseudocount of 1e-5), and SNR is `(prenatal CPM − postnatal CPM) / postnatal CPM SD`.
- For each prenatal group, take the minimum fold change and minimum SNR across all postnatal groups.
- A junction is **oncofetal** when at least one prenatal group has a minimum fold change greater than 2 and a minimum SNR greater than 2 relative to all postnatal groups.
- A junction is not eligible if it has no observed postnatal CPM in any group, or if its prenatal CPM is missing in a given group (missing values are never treated as evidence of prenatal expression). If no postnatal SD is available, SNR is missing and the criterion is not met.
- The prenatal group(s) meeting the criteria are recorded in `oncofetal_prenatal_group` (semicolon-separated). The best minimum fold change and SNR across prenatal groups are recorded in `max_prenatal_min_cpm_fc` and `max_prenatal_min_cpm_snr`.
- Tumor-enriched junctions that do not meet the oncofetal criteria keep the `Tumor-enriched` label (called "tumor-specific" in the TAPESTRY web app).

The oncofetal call is made per junction using Evo-Devo only, so every sample carrying a junction receives the same oncofetal/tumor-enriched label; the tumor-enriched call itself is made per sample.

### 4. Domain annotation and expression filtering

- `03-uniprot-domain-annotation.sh` restricts junctions to protein-coding CDS (GENCODE v39) and annotates overlap with UniProt extracellular, transmembrane and cytoplasmic topological domains. Domains lying completely within the junction interval are excluded.
- `04-pfam-annotation.R` annotates junctions to Pfam domains (`annoFuseData`) by gene symbol and classifies the overlap as domain completely in junction interval, junction interval completely in domain, or partial overlap.
- `05-domain-expression-filtering.R` merges the annotations, adds `in_cds`, retains junctions in samples where the gene's TPM is >= 10 (`gene_tpm`), and defines `consequence` (`EC domain alteration`, `Other functional`, both, or `Unknown`) from the UniProt and Pfam annotations.

## Folder content
* `run_module.sh` shell script to run analysis
* `01-classify-tumor-enriched-junctions.R` identify tumor-enriched junctions by comparing PBTA expression with normal control cohorts.
* `02-classify-oncofetal-junctions.R` classify tumor-enriched junctions as oncofetal using prenatal versus postnatal Evo-Devo expression, and merge with junction annotations.
* `03-uniprot-domain-annotation.sh` annotate junctions to UniProt topological domains.
* `04-pfam-annotation.R` annotate junctions to Pfam functional domains.
* `05-domain-expression-filtering.R` merge junction annotations and retain junctions in genes with TPM >= 10.
* `util/other-functions.R` helper functions, including the control STAR junction cross-check.
* `util/add-tpm-values.R` helper function to append gene TPM values.

## Input files
From other modules (`analyses/`):
* `02-pbta-junction-processing/results/pbta-merged-norm-batch-corrected-junction-cts.qs2` batch-corrected PBTA junction CPMs.
* `02-pbta-junction-processing/results/junction-annot.tsv.gz` PBTA splice-junction annotations.
* `01-ctrl-rmats-processing/results/` mean and SD junction CPM matrices (`-ct-mat.qs2` / `-sd-mat.qs2`) for:
  * GTEx (`gtex-merged-norm-junction-*`);
  * postnatal Evo-Devo collapsed to forebrain/hindbrain (`evodevo-merged-postnatal-norm-junction-*`);
  * normal pediatric brain (`normal-pedbrain-merged-norm-junction-*`);
  * postnatal Evo-Devo by region and stage (`evodevo-merged-postnatal-stage-norm-junction-*`); and
  * prenatal Evo-Devo by region and week bin, mean only (`evodevo-merged-prenatal-week-binned-norm-junction-ct-mat.qs2`).

From `data/` (see `download_data.sh`):
* `SJ.merged.control-cohort.tsv.gz` normal brain STAR junction counts, used to exclude junctions observed in controls.
* `gencode.v39.primary_assembly.annotation.gtf.gz` GENCODE v39 annotation.
* `unipLocExtra.hg38.col.txt`, `unipLocTransMemb.hg38.col.txt`, `unipLocCytopl.hg38.col.txt` UniProt topological domain coordinates.
* `pbta_gene-expression-rsem-tpm-collapsed.rds` PBTA gene expression (TPM).

## Output files
Columns of `tumor-enriched-oncofetal-splice-junctions.tsv.gz` include `junction`, `sample_id`, `junction_count`, `junction_cpm`, `boundary`, `criteria`, `junction_preference` (`Tumor-enriched` or `Oncofetal`), `min_cpm_fc_all`, `min_cpm_snr_all`, `max_mean_cpm_all` (missing for junctions absent from the references), `max_prenatal_min_cpm_fc`, `max_prenatal_min_cpm_snr`, `oncofetal_prenatal_group`, and junction annotation columns (`strand`, `chr`, `up_jc_start`, `up_jc_end`, `down_jc_start`, `down_jc_end`, `geneSymbol`).

## Directory structure
```
.
├── 01-classify-tumor-enriched-junctions.R
├── 02-classify-oncofetal-junctions.R
├── 03-uniprot-domain-annotation.sh
├── 04-pfam-annotation.R
├── 05-domain-expression-filtering.R
├── README.md
├── results
│   ├── tumor-enriched-junctions.qs2
│   ├── tumor-enriched-oncofetal-splice-junctions.bed
│   ├── tumor-enriched-oncofetal-splice-junctions-cds.bed
│   ├── tumor-enriched-oncofetal-splice-junctions.unipLocCytopl.bed
│   ├── tumor-enriched-oncofetal-splice-junctions.unipLocExtra.bed
│   ├── tumor-enriched-oncofetal-splice-junctions.unipLocTransMemb.bed
│   ├── tumor-enriched-oncofetal-splice-junctions-domain-anno.uniq.tsv
│   ├── tumor-enriched-oncofetal-splice-junctions-domain-expr-annotated.tsv.gz
│   ├── tumor-enriched-oncofetal-splice-junctions-pfam-annotated.tsv.gz
│   └── tumor-enriched-oncofetal-splice-junctions.tsv.gz
├── run_module.sh
└── util
    ├── add-tpm-values.R
    └── other-functions.R
```
