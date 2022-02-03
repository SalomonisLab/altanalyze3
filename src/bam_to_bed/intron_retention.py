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
###############################################################################################################################################################


class Cat (enum.Enum):

    PRIME_5 = 1
    PRIME_3 = 2
    INTRON = 3
    DISCARD = 4

class Overlaps:

    def __init__(self):
        self.reset()

    def reset(self):
        self.overlaps = {}

    def __getitem__(self, key):                                          # key is a tuple (contig, start, end, name)
        return self.overlaps.setdefault(key, {"p5": 0, "p3": 0})

    def __iter__(self):
        for (contig, start, end, name, strand), value in self.overlaps.items():
            yield contig, start, end, name, strand, value["p5"], value["p3"]

    def increment(self, contig, start, end, name, strand, category, step=None):
        step = 1 if step is None else step
        if category is Cat.PRIME_5:
            self[(contig, start, end, name, strand)]["p5"] += step
        elif category is Cat.PRIME_3:
            self[(contig, start, end, name, strand)]["p3"] += step


class Counter:

    def __init__(self, bam, ref, span, strandness, threads=None):
        self.bam = bam
        self.ref = ref
        self.span = span
        self.strandness = strandness
        self.threads = 1 if threads is None else threads
        self.paired = self.__is_paired()
        self.reset()

    def reset(self):
        self.cache = {}
        self.overlaps = Overlaps()

    def get_category(self, read, intron, span=None, strandness=None):
        span = 0 if span is None else span
        if strandness=="forward":
            if self.__is_paired:
                if read.is_read1:
                    if intron.strand=="+" and read.is_reverse:
                        return Cat.DISCARD
                    if intron.strand=="-" and not read.is_reverse:
                        return Cat.DISCARD
                elif read.is_read2:
                    if intron.strand=="+" and not read.is_reverse:
                        return Cat.DISCARD
                    if intron.strand=="-" and read.is_reverse:
                        return Cat.DISCARD
            else:
                if intron.strand=="+" and read.is_reverse:
                    return Cat.DISCARD
        elif strandness=="reverse":
            if self.__is_paired:
                if read.is_read1:
                    if intron.strand=="+" and not read.is_reverse:
                        return Cat.DISCARD
                    if intron.strand=="-" and read.is_reverse:
                        return Cat.DISCARD
                elif read.is_read2:
                    if intron.strand=="+" and read.is_reverse:
                        return Cat.DISCARD
                    if intron.strand=="-" and not read.is_reverse:
                        return Cat.DISCARD
            else:
                if intron.strand=="-" and not read.is_reverse:
                    return Cat.DISCARD
        if read.reference_end - intron.start >= span and intron.start - read.reference_start >= span:
            return Cat.PRIME_5
        elif read.reference_start - intron.start >= 0 and intron.end - read.reference_end >= 0:
            return Cat.INTRON
        elif read.reference_end - intron.end >= span and intron.end - read.reference_start >= span:
            return Cat.PRIME_3
        else:
            return Cat.DISCARD

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
        contig, start, end, name, strand, category = intron_data
        if cached_data is None:
            self.overlaps.increment(contig, start, end, name, strand, category)
        else:
            cached_contig, cached_start, cached_end, cached_name, cached_strand, cached_category = cached_data
            if cached_category is Cat.INTRON:
                self.overlaps.increment(contig, start, end, name, strand, category)
            elif category is Cat.INTRON:
                self.overlaps.increment(cached_contig, cached_start, cached_end, cached_name, cached_strand, cached_category)

    def skip_read(self, read):
        return read.is_secondary or \
               read.is_duplicate or \
               read.is_unmapped or \
               read.is_supplementary or \
               (read.is_paired and read.mate_is_unmapped)

    def correct_contig(self, contig, handler):
        try:
            handler.fetch(contig)
            return contig
        except ValueError:
            return contig.lstrip("chr")

    def calculate(self, contig):
        with pysam.AlignmentFile(self.bam, mode="rb", threads=self.threads) as bam_handler:
            with pysam.TabixFile(self.ref, mode="r", parser=pysam.asBed(), threads=self.threads) as ref_handler:
                contig_ref = self.correct_contig(contig, ref_handler)  # contig can be with or without 'chr' prefix so we need to try to fetch both and see which one works
                contig_bam = self.correct_contig(contig, bam_handler)
                intron_iter = ref_handler.fetch(contig_ref)
                intron = next(intron_iter)                             # get initial value from intron iterator
                for read in bam_handler.fetch(contig_bam):
                    if self.skip_read(read):                           # skip all "bad" reads
                        continue
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
                            intron.strand,
                            self.get_category(read, intron, self.span, self.strandness)
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

    def export(self, location):
        print(f"Export temporary results to {location}")
        with open(location, "w") as out_handler:
            for contig, start, end, name, strand, p5, p3 in self.overlaps:
                # out_handler.write(f"{contig}\t{start}\t{end}\t{name}\t0\t{p5}\t{p3}\n")
                out_handler.write(f"{contig}\t{start-self.span}\t{start+self.span}\t{name}-{start}\t{p5}\t{strand}\n")
                out_handler.write(f"{contig}\t{end-self.span}\t{end+self.span}\t{name}-{end}\t{p3}\t{strand}\n")


def chr_decorator(function):
    """
    Decorator to prepend any str or [str] with 'chr' prefix if it was missing
    """
    def __prefix(c):
        return c if c.startswith("chr") else f"chr{c}"
    def __wrapper(contig):
        raw_res = function(contig)
        if isinstance(raw_res, list):
            return [__prefix(c) for c in raw_res]
        else:
            return __prefix(raw_res)
    return __wrapper


def get_jobs(args):
    marker = "".join(random.choices(string.ascii_uppercase + string.digits, k=10))
    return [
        (contig, args.output + "__" + marker)                                # contig is always prepended with 'chr'
            for contig in get_all_bam_chr(args)
                if contig in get_all_ref_chr(args) and contig in args.chr    # safety measure to include only chromosomes present in BAM, BED, and --chr
    ]


def process_contig(args, job):
    print(f"Process {job[0]}")
    counter = Counter(
        bam=args.bam,
        ref=args.ref,
        span=args.span,
        strandness=args.strandness,
        threads=args.threads
    )
    counter.calculate(job[0])
    counter.export(job[1]+"__"+job[0])


def collect_results(args, jobs):
    with open(args.output + ".bed", "w") as output_stream:
        for contig, location in jobs:
            chunk = location + "__" + contig
            print(f"Collect results from {chunk}")
            with open(chunk, "r") as input_stream:
                output_stream.write(input_stream.read())
                os.remove(chunk)


@chr_decorator
def get_all_bam_chr(args):
    with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as bam_handler:
        return [_.contig for _ in bam_handler.get_index_statistics()]


@chr_decorator
def get_all_ref_chr(args):
    with pysam.TabixFile(args.ref, mode="r", parser=pysam.asBed(), threads=args.threads) as ref_handler:
        return ref_handler.contigs


def assert_args(args):
    args.chr = get_all_bam_chr(args) if len(args.chr) == 0 else [c if c.startswith("chr") else f"chr{c}" for c in args.chr]
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
    general_parser.add_argument("--strandness",
        help=" ".join(
            [
                "Takes strand information into account when counting overlaps.",
                "For 'forward' and single-end reads, the read has to be mapped",
                "to the same strand as intron. For paired-end reads, the first",
                "read has to be on the same strand as intron and the second read",
                "on the opposite strand. For 'reverse', these rules are reversed.",
                "Default: strand information is ignored"
            ]
        ),
        choices=["forward", "reverse"]
    )
    general_parser.add_argument("--threads",  help="Number of threads to decompress BAM file", type=int, default=1)
    general_parser.add_argument("--cpus",     help="Number of processes to run in parallel", type=int, default=1)
    general_parser.add_argument("--chr",      help="Select chromosomes to process. Default: all available", type=str, nargs="*", default=[])
    general_parser.add_argument("--output",   help="Output file prefix", type=str, default="intron")
    return general_parser


def main(argsl=None):
    if argsl is None:
        argsl = sys.argv[1:]
    args, _ = arg_parser().parse_known_args(argsl)
    args = assert_args(normalize_args(args, ["span", "threads", "cpus", "chr", "strandness"]))
    start = time.time()
    jobs = get_jobs(args)
    with Pool(args.cpus) as pool:
        pool.map(partial(process_contig, args), jobs)
    collect_results(args, jobs)
    print (f"""Elapsed time: {time.time() - start}""")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
