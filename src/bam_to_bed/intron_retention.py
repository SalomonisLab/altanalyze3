import os
import sys
import time
import pysam
import string
import random
import argparse
from functools import partial
from multiprocessing import Pool


###############################################################################################################################################################
# CURRENT INSTRUCTIONS
#
# To obtain properly formatted coordinate sorted indexed gene model reference BED file
# that includes only introns, do the following actions with your AltAnalyze formatted
# Hs_Ensembl_exon.txt (for example) file.
#
# 1. Make sure you have these columns in you input Hs_Ensembl_exon.txt file
#    head -n 1 Hs_Ensembl_exon.txt
#    gene	exon-id	chromosome	strand	exon-region-start(s)	exon-region-stop(s)	constitutive_call	ens_exon_ids	splice_events	splice_junctions
#
# 2. Convert to BED format excluding all not intron records
#    cat Hs_Ensembl_exon.txt | grep -v "gene" | awk '$2 ~ /^I/ {print $3"\t"$5"\t"$6"\t"$2"-"$1"\t"0"\t"$4}' | sort -k1,1 -k2,2n -k3,3n | bgzip > hs_ref.bed.gz
#    tabix -p bed hs_ref.bed.gz
###############################################################################################################################################################


#################################################################################################
# DEPRECATED INSTRUCTIONS
# keep it in case we decide to support GTF input in future
#
# (grep ^"#" input.gtf; grep -v ^"#" input.gtf | sort -k1,1 -k4,4n) | bgzip > input.sorted.gtf.gz
# tabix -p gff input.sorted.gtf.gz
#################################################################################################


def is_paired(args):
    with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as bam_handler:
        bam_iter = bam_handler.fetch()
        for i in range(20):                       # check max 20 reads
            try:
                bam_handler.mate(next(bam_iter))  # Raises ValueError – if the read is unpaired or the mate is unmapped
            except ValueError:
                return False
        return True


def get_default_chr(args):
    with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as bam_handler:
        return [_.contig for _ in bam_handler.get_index_statistics()]


def mate_mapped_to_intron(pileups):  # TODO
    return pileups


def count_alignments(column_iter, filter_by_mate=None):
    filter_by_mate = False if filter_by_mate is None else filter_by_mate
    collected_alignments = []
    for column in column_iter:
        if filter_by_mate:
            collected_alignments.extend(mate_mapped_to_intron(column.pileups))
        else:
            collected_alignments.extend(column.pileups)
    return len(collected_alignments)


def process_contig(args, job):
    contig = job[0]
    location = job[1]
    print(f"""Process {contig} into {location}""")
    with open(location + "__" + contig, "w") as output_stream:
        with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as bam_handler:
            with pysam.TabixFile(args.ref, mode="r", parser=pysam.asBed(), threads=args.threads) as ref_handler:
                for intron in ref_handler.fetch(contig):
                    five_prime_counts = count_alignments(
                        column_iter=bam_handler.pileup(       # should return exactly one item
                            contig=contig,
                            start=intron.start-args.span,
                            end=intron.start-args.span+1,
                            truncate=True                     # only columns in the exact region are returned
                        ),
                        filter_by_mate=args.paired
                    )
                    three_prime_counts = count_alignments(
                        column_iter=bam_handler.pileup(      # should return exactly one item
                            contig=contig,
                            start=intron.end+args.span,
                            end=intron.end+args.span+1,
                            truncate=True                    # only columns in the exact region are returned
                        ),
                        filter_by_mate=args.paired
                    )
                    if five_prime_counts > 0 or three_prime_counts > 0:
                        output_stream.write(
                            f"""{contig}\t{intron.start}\t{intron.end}\t{intron.name}\t0\t{five_prime_counts}\t{three_prime_counts}\n"""
                        )


def collect_results(args, jobs):
    with open(args.output + ".bed", "w") as output_stream:
        for contig, location in jobs:
            chunk = location + "__" + contig
            with open(chunk, "r") as input_stream:
                output_stream.write(input_stream.read())
                os.remove(chunk)


def get_jobs(args):
    marker = "".join(random.choices(string.ascii_uppercase + string.digits, k=10))
    with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as bam_handler:
        with pysam.TabixFile(args.ref, mode="r", parser=pysam.asBed(), threads=args.threads) as ref_handler:
            return [
                (_.contig, args.output + "__" + marker)
                    for _ in bam_handler.get_index_statistics()
                        if _.contig in ref_handler.contigs and _.contig in args.chr       # safety measure to include only chromosomes present in BAM, BED, and --chr
            ]


def normalize_args(args, skip_list=[]):
    normalized_args = {}
    for key,value in args.__dict__.items():
        if key not in skip_list:
            normalized_args[key] = value if not value or os.path.isabs(value) else os.path.normpath(os.path.join(os.getcwd(), value))
        else:
            normalized_args[key]=value
    return argparse.Namespace (**normalized_args)


def assert_args(args):
    args.paired = is_paired(args)
    args.chr = get_default_chr(args) if len(args.chr) == 0 else args.chr
    return args


def arg_parser():
    general_parser = argparse.ArgumentParser()
    general_parser.add_argument("--bam",     help="Path to the coordinate-sorted indexed BAM file (should include only mapped reads)", type=str, required=True)
    general_parser.add_argument("--ref",     help="Path to the coordinate-sorted indexed gene model reference BED file", type=str, required=True)
    general_parser.add_argument("--span",    help="5' and 3' overlap that read should have over a splice-site to be counted", type=int, default=10)
    general_parser.add_argument("--threads", help="Number of threads to decompress BAM file", type=int, default=1)
    general_parser.add_argument("--cpus",    help="Number of processes to run in parallel", type=int, default=1)
    general_parser.add_argument("--chr",     help="Select chromosomes to process. Default: all available", type=str, nargs="*", default=[])
    general_parser.add_argument("--output",  help="Output file prefix", type=str, default="intron")
    return general_parser


def main(argsl=None):
    if argsl is None:
        argsl = sys.argv[1:]
    args, _ = arg_parser().parse_known_args(argsl)
    args = assert_args(normalize_args(args, ["span", "threads", "cpus", "chr"]))
    start = time.time()
    jobs = get_jobs(args)
    with Pool(args.cpus) as pool:
        pool.map(partial(process_contig, args), jobs)
    collect_results(args, jobs)
    print (f"""Elapsed time: {time.time() - start}""")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
