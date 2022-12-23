#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

$DIR/data_for_tests.sh

TMP=$DIR/data/e2e_tmp

mkdir -p $TMP
cd $TMP

echo "Count junctions reads from Cal27P5-1.bam"
altanalyze3 juncount --bam ../Cal27P5-1.bam \
                     --output Cal27P5_1_juncounts

echo "Count introns reads from Cal27P5-1.bam"
altanalyze3 intcount --bam ../Cal27P5-1.bam \
                     --ref ../gene_model_all.tsv \
                     --output Cal27P5_1_intcounts \

echo "Aggregating counts from Cal27P5_juncounts.bed and Cal27P5_intcounts.bed"
altanalyze3 aggregate --juncounts ./Cal27P5_1_juncounts.bed \
                      --intcounts ./Cal27P5_1_intcounts.bed \
                      --ref ../gene_model_all.tsv \
                      --bed \
                      --loglevel info \
                      --output Cal27P5_1_aggregated_counts