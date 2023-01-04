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

diff Cal27P5_1_aggregated_counts.bed.gz ../controls/Cal27P5_1_aggregated_counts.bed.gz
DIFF_ERROR=$?
if [ $DIFF_ERROR -eq 2 ]
then
   echo "something wrong with the diff command"
   exit 1
elif [ $DIFF_ERROR -eq 1 ]
then
   echo "Cal27P5_1_aggregated_counts.bed.gz file is not equal to the control one"
   exit 1
else
   echo "success"
fi


echo "Count junctions reads from Cal27P5-2.bam"
altanalyze3 juncount --bam ../Cal27P5-2.bam \
                     --output Cal27P5_2_juncounts

echo "Count introns reads from Cal27P5-2.bam"
altanalyze3 intcount --bam ../Cal27P5-2.bam \
                     --ref ../gene_model_all.tsv \
                     --output Cal27P5_2_intcounts \

echo "Count junctions reads from Cal27P5-3.bam"
altanalyze3 juncount --bam ../Cal27P5-3.bam \
                     --output Cal27P5_3_juncounts

echo "Count introns reads from Cal27P5-3.bam"
altanalyze3 intcount --bam ../Cal27P5-3.bam \
                     --ref ../gene_model_all.tsv \
                     --output Cal27P5_3_intcounts \

echo "Aggregating counts from Cal27P5_[1/2/3]_juncounts.bed and Cal27P5_[1/2/3]_intcounts.bed"
altanalyze3 aggregate --juncounts ./Cal27P5_1_juncounts.bed ./Cal27P5_2_juncounts.bed ./Cal27P5_3_juncounts.bed \
                      --intcounts ./Cal27P5_1_intcounts.bed ./Cal27P5_2_intcounts.bed ./Cal27P5_3_intcounts.bed \
                      --ref ../gene_model_all.tsv \
                      --bed \
                      --loglevel info \
                      --output Cal27P5_1_2_3_aggregated_counts