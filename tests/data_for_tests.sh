#!/bin/bash

WORKING_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
DATA_DIR="${WORKING_DIR}/data"

mkdir -p $DATA_DIR
cd $DATA_DIR
echo "Fetching data for tests into ${DATA_DIR} directory. Previously downloaded files will be overwritten."

wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/gene_model_all.tsv.gz -O gene_model_all.tsv.gz
gzip -d -f gene_model_all.tsv.gz

wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/Cal27P5-1.bam -O Cal27P5-1.bam
wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/Cal27P5-1.bam.bai -O Cal27P5-1.bam.bai

wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/Cal27P5-2.bam -O Cal27P5-2.bam
wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/Cal27P5-2.bam.bai -O Cal27P5-2.bam.bai

wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/Cal27P5-3.bam -O Cal27P5-3.bam
wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/Cal27P5-3.bam.bai -O Cal27P5-3.bam.bai

cp Cal27P5-1.bam Cal27P5-1-copy.bam

mkdir -p controls && cd controls
wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/controls/Cal27P5_1_aggregated_counts.bed.gz -O Cal27P5_1_aggregated_counts.bed.gz
wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/controls/Cal27P5_1_intcounts.bed -O Cal27P5_1_intcounts.bed
wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/controls/Cal27P5_1_juncounts.bed -O Cal27P5_1_juncounts.bed
wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/controls/Cal27P5_2_intcounts.bed -O Cal27P5_2_intcounts.bed
wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/controls/Cal27P5_2_juncounts.bed -O Cal27P5_2_juncounts.bed
wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/controls/Cal27P5_3_intcounts.bed -O Cal27P5_3_intcounts.bed
wget -q --show-progress https://github.com/michael-kotliar/altanalyze3_data/raw/main/data/controls/Cal27P5_3_juncounts.bed -O Cal27P5_3_juncounts.bed
cd ..