import os
import sys
import time
import enum
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
#
# If input BAM file is paired-end, it should include only primary alignment (only two reads should have identical names)
###############################################################################################################################################################


class Category (enum.Enum):

    PRIME_5 = 1
    PRIME_3 = 2
    INTRON = 3
    DISCARD = 4


class Overlaps:

    def __init__(self):
        self.reset()

    def reset(self):
        self.overlaps = {}

    def __getitem__(self, contig, start, end, name):
        return self.overlaps.setdefault(
            (contig, start, end, name),
            {
                "p5": 0,
                "p3": 0
            }
        )

    def __iter__(self):
        for (contig, start, end, name), value in self.overlaps.items():
            yield contig, start, end, name, value["p5"], value["p3"]

    def increment(self, contig, start, end, name, category, step=None):
        step = 1 if step is None else step
        current = self.__getitem__(contig, start, end, name)
        if category == Category.PRIME_5:
            self.overlaps[(contig, start, end, name)] = {
                "p5": current["p5"] + step,
                "p3": current["p3"]
            }
        elif category == Category.PRIME_3:
            self.overlaps[(contig, start, end, name)] = {
                "p5": current["p5"],
                "p3": current["p3"] + step
            }


class Counter:

    def __init__(self, bam, ref, threads=None):
        self.bam = bam
        self.ref = ref
        self.threads = 1 if threads is None else threads
        self.paired = self.__is_paired()
        self.reset()

    def reset(self):
        self.cache = {}
        self.overlaps = Overlaps()

    def get_category(self, read, intron, span=None):
        span = 0 if span is None else span
        if read.reference_end - intron.start >= span and intron.start - read.reference_start >= span:
            return Category.PRIME_5
        elif read.reference_start - intron.start >= 0 and intron.end - read.reference_end >= 0:
            return Category.INTRON
        elif read.reference_end - intron.end >= span and intron.end - read.reference_start >= span:
            return Category.PRIME_3
        else:
            return Category.DISCARD

    def __is_paired(self, n_reads=None):
        n_reads = 20 if n_reads is None else n_reads
        with pysam.AlignmentFile(self.bam, mode="rb", threads=self.threads) as bam_handler:
            bam_iter = bam_handler.fetch()
            for i in range(n_reads):
                try:
                    bam_handler.mate(next(bam_iter))
                except ValueError:
                    return False
            return True

    def process_overlaps(self, intron_data, cached_data):
        contig, start, end, name, category = intron_data
        if cached_data is None:
            self.overlaps.increment(contig, start, end, name, category)
        else:
            cached_contig, cached_start, cached_end, cached_name, cached_category = cached_data
            if cached_category is Category.INTRON:
                self.overlaps.increment(contig, start, end, name, category)
            elif category is Category.INTRON:
                self.overlaps.increment(cached_contig, cached_start, cached_end, cached_name, cached_category)

    def calculate(self, contig, span):
        with pysam.AlignmentFile(self.bam, mode="rb", threads=self.threads) as bam_handler:
            with pysam.TabixFile(self.ref, mode="r", parser=pysam.asBed(), threads=self.threads) as ref_handler:
                intron_iter = ref_handler.fetch(contig)
                intron = next(intron_iter)                         # get initial value from intron iterator
                for read in bam_handler.fetch(contig):
                    if read.reference_start - intron.end >= 0:
                        try:
                            intron = next(intron_iter)
                        except StopIteration:
                            break
                    else:
                        intron_data = (
                            contig,
                            intron.start,
                            intron.end,
                            intron.name,
                            self.get_category(read, intron, span)
                        )
                        if self.paired and read.query_name not in self.cache:
                            self.cache[read.query_name] = intron_data
                        else:
                            self.process_overlaps(
                                intron_data,
                                self.cache.get(read.query_name, None)
                            )
                            try:
                                del self.cache[read.query_name]
                            except KeyError:
                                pass

    def export(self, out):
        print(f"Export temporary results to {out}")
        with open(out, "w") as out_handler:
            for contig, start, end, name, p5, p3 in self.overlaps:
                out_handler.write(f"""{contig}\t{start}\t{end}\t{name}\t0\t{p5}\t{p3}\n""")


def get_jobs(args):
    marker = "".join(random.choices(string.ascii_uppercase + string.digits, k=10))
    with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as bam_handler:
        with pysam.TabixFile(args.ref, mode="r", parser=pysam.asBed(), threads=args.threads) as ref_handler:
            return [
                (_.contig, args.output + "__" + marker)
                    for _ in bam_handler.get_index_statistics()
                        if _.contig in ref_handler.contigs and _.contig in args.chr       # safety measure to include only chromosomes present in BAM, BED, and --chr
            ]


def process_contig(args, job):
    counter = Counter(args.bam, args.ref, args.threads)
    counter.calculate(contig=job[0], span=args.span)
    counter.export(out=job[1]+"__"+job[0])


def collect_results(args, jobs):
    with open(args.output + ".bed", "w") as output_stream:
        for contig, location in jobs:
            chunk = location + "__" + contig
            print(f"Collect results from {chunk}")
            with open(chunk, "r") as input_stream:
                output_stream.write(input_stream.read())
                os.remove(chunk)


def get_default_chr(args):
    with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as bam_handler:
        return [_.contig for _ in bam_handler.get_index_statistics()]


def assert_args(args):
    args.chr = get_default_chr(args) if len(args.chr) == 0 else args.chr
    return args


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
