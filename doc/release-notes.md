# release notes

## current release (v4)
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

## current release (v3)
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


## current release (v2)
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
