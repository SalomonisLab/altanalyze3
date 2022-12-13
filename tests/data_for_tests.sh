#!/bin/bash

WORKING_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
DATA_DIR="${WORKING_DIR}/data"

mkdir -p $DATA_DIR
cd $DATA_DIR
echo "Fetching data for tests into ${DATA_DIR} directory. Previously downloaded files will be overwritten."

wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/hg19_ref.tsv.gz -O hg19_ref.tsv.gz
gzip -d -f hg19_ref.tsv.gz
cat hg19_ref.tsv | grep -v "gene" | awk '$2 ~ /^I/ {print $3"\t"$5"\t"$6"\t"$2"-"$1"\t"0"\t"$4}' | sort -k1,1 -k2,2n -k3,3n | bgzip > hg19_ref_introns.bed.gz
tabix -p bed hg19_ref_introns.bed.gz

wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/hg19_pe.bam -O hg19_pe.bam
wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/hg19_pe.bam.bai -O hg19_pe.bam.bai
cp hg19_pe.bam hg19_pe_not_indexed.bam
cp hg19_pe.bam hg19_pe_copy.bam
cp hg19_pe.bam.bai hg19_pe_copy.bam.bai