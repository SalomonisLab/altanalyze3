#!/bin/bash

WORKING_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
DATA_DIR="${WORKING_DIR}/data"

mkdir -p $DATA_DIR
cd $DATA_DIR
echo "Fetching data for tests into ${DATA_DIR} directory. Previously downloaded files will be overwritten."

wget -q --show-progress \
    https://github.com/nsalomonis/BAM-to-Junction-BED/raw/master/ReferenceExonCoordinates/Mm_Ensembl_exon_Mm10.txt -O Mm_Ensembl_exon_Mm10.tsv
cat Mm_Ensembl_exon_Mm10.tsv | grep -v "gene" | awk '$2 ~ /^I/ {print $3"\t"$5"\t"$6"\t"$2"-"$1"\t"0"\t"$4}' | sort -k1,1 -k2,2n -k3,3n | bgzip > Mm_Ensembl_exon_Mm10_filtered.bed.gz
tabix -p bed Mm_Ensembl_exon_Mm10_filtered.bed.gz