#!/bin/bash

WORKING_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
DATA_DIR="${WORKING_DIR}/data"

mkdir -p $DATA_DIR
cd $DATA_DIR
echo "Fetching data for tests into ${DATA_DIR} directory. Previously downloaded files will be overwritten."

wget -q --show-progress \
    https://github.com/nsalomonis/BAM-to-Junction-BED/raw/master/ReferenceExonCoordinates/Mm_Ensembl_exon_Mm10.txt -O Mm_Ensembl_exon_Mm10.tsv