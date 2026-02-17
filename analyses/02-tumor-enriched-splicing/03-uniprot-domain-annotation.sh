#!/bin/sh

# set results dir
output_dir="results"

# set file names
# input
input_filename="$output_dir/tumor-enriched-oncofetal-splice-junctions.bed"

# output files
output_cds_filename="$output_dir/tumor-enriched-oncofetal-splice-junctions-cds.bed"
output_cds_ec_filename="${output_dir%.*}/tumor-enriched-oncofetal-splice-junctions.unipLocExtra.bed"
output_cds_tmem_filename="${output_dir%.*}/tumor-enriched-oncofetal-splice-junctions.unipLocTransMemb.bed"
output_cds_cyto_filename="${output_dir%.*}/tumor-enriched-oncofetal-splice-junctions.unipLocCytopl.bed"
output_target_list_filename="${output_dir%.*}/tumor-enriched-oncofetal-splice-junctions-domain-anno.uniq.tsv"

## Filter tumor-enriched junction bed for coding exons only
bedtools intersect -wo -a "$input_filename" -b <(zgrep 'CDS.*transcript_type "protein_coding"' ../../data/gencode.v39.primary_assembly.annotation.gtf.gz) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | sort -u > "$output_cds_filename"

## Run bedtools intersect command against extracellular, transmembrane, and intra-cellular/cytosolic domains to get domain coordinates
bedtools intersect -wo -a "$output_cds_filename" -b ../../data/unipLocExtra.hg38.col.txt | \
awk '{
    overlap_start = ($2 > $9) ? $2 : $9;
    overlap_end = ($3 < $10) ? $3 : $10;
    sample_id = $5;
    print $0"\t"overlap_start"\t"overlap_end}' > "$output_cds_ec_filename"

bedtools intersect -wo -a "$output_cds_filename" -b ../../data/unipLocTransMemb.hg38.col.txt |\
awk '{
    overlap_start = ($2 > $8) ? $2 : $9;
    overlap_end = ($3 < $9) ? $3 : $10;
    sample_id = $5;
    print $0"\t"overlap_start"\t"overlap_end }' > "$output_cds_tmem_filename"

bedtools intersect -wo -a "$output_cds_filename" -b ../../data/unipLocCytopl.hg38.col.txt |\
awk '{
    overlap_start = ($2 > $8) ? $2 : $9;
    overlap_end = ($3 < $9) ? $3 : $10;
    sample_id = $5;
    print $0"\t"overlap_start"\t"overlap_end }' > "$output_cds_cyto_filename"

## Extract overlapping junctions for each domain type
echo -e  "junction\tjunction_preference\tcoverage\tdomain_type\toverlap_start\toverlap_end" > "$output_target_list_filename"
cat "$output_cds_ec_filename" | awk '{print $4"\t"$5"\t"($11/($3-$2))*100"\tEC\t"$12"\t"$13}' | sort -u  >> "$output_target_list_filename"
cat "$output_cds_tmem_filename" | awk '{print $4"\t"$5"\t"($11/($3-$2))*100"\tTM\t"$12"\t"$13}' | sort -u >> "$output_target_list_filename"
cat "$output_cds_cyto_filename" | awk '{print $4"\t"$5"\t"($11/($3-$2))*100"\tIC\t"$12"\t"$13}' | sort -u  >> "$output_target_list_filename"
