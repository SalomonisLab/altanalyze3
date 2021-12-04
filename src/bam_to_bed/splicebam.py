import os
import sys
import gzip
import time
import pysam
import string
import random
import argparse
from multiprocessing import Pool
from functools import partial


def process_contig(args, job):
    contig = job[0]
    location = job[1]
    with gzip.open(location + "__" + contig, "wt") as output_stream:
        with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as handler:
            introns = handler.find_introns((r for r in handler.fetch(contig)))           # need () instead of [] to use as iterator
            for position, score in introns.items():
                output_stream.write(
                    f"""{contig}\t{position[0]}\t{position[1]}\tJUNC:{contig}-{position[0]}-{position[1]}\t{score}\n"""
                )


def collect_results(args, jobs):
    with open(args.output + ".bed.gz", "wb") as output_stream:
        for contig, location in jobs:
            chunk = location + "__" + contig
            with open(chunk, "rb") as input_stream:
                output_stream.write(input_stream.read())
                os.remove(chunk)


def get_jobs(args):
    marker = "".join(random.choices(string.ascii_uppercase + string.digits, k=10))
    with pysam.AlignmentFile(args.bam, mode="rb") as handler:
        return [(_.contig, args.output + "__" + marker) for _ in handler.get_index_statistics()]


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
    general_parser.add_argument("--cpus",    help="Number of processes to run in parallel", type=int, default=1)
    general_parser.add_argument("--output",  help="Output file prefix", type=str, default="spliced")
    return general_parser


def main(argsl=None):
    if argsl is None:
        argsl = sys.argv[1:]
    args, _ = arg_parser().parse_known_args(argsl)
    args = normalize_args(args, ["threads", "cpus"])
    start = time.time()
    jobs = get_jobs(args)
    with Pool(args.cpus) as pool:
        pool.map(partial(process_contig, args), jobs)
    collect_results(args, jobs)
    print (f"""Elapsed time: {time.time() - start}""")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

