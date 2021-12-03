import os
import sys
import gzip
import time
import pysam
import argparse


def process_introns(args):
    with gzip.open(args.output, "wt") as output_stream:
        with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as handler:
            for chr_stats in handler.get_index_statistics():
                print (f"""Processing {chr_stats.contig}""")
                chr_introns = handler.find_introns((r for r in handler.fetch(chr_stats.contig)))
                for position, score in chr_introns.items():
                    output_stream.write(
                        f"""{chr_stats.contig}\t{position[0]}\t{position[1]}\tJUNC:{chr_stats.contig}-{position[0]}-{position[1]}\t{score}\n"""
                    )


def normalize_args(args, skip_list=[]):
    normalized_args = {}
    for key,value in args.__dict__.items():
        if key not in skip_list:
            normalized_args[key] = value if not value or os.path.isabs(value) else os.path.normpath(os.path.join(os.getcwd(), value))
        else:
            normalized_args[key]=value
    return argparse.Namespace (**normalized_args)


def arg_parser():
    general_parser = argparse.ArgumentParser()
    general_parser.add_argument("--bam",     help="Path to the coordinate-sorted indexed BAM file", type=str, required=True)
    general_parser.add_argument("--threads", help="Number of threads to decompress BAM file", type=int, default=1)
    general_parser.add_argument("--output",  help="Output file name", type=str, default="spliced.bed.gz")
    return general_parser


def main(argsl=None):
    if argsl is None:
        argsl = sys.argv[1:]
    args, _ = arg_parser().parse_known_args(argsl)
    args = normalize_args(args, ["threads"])
    start = time.time()
    process_introns(args)
    print (f"""Elapsed time: {time.time() - start}""")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

