import os
import sys
import time
import enum
import pysam
import string
import random
import logging
import argparse
import multiprocessing
from functools import partial


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


def setup_logger(logger, log_level, log_format=None):
    log_format = "%(processName)12s (%(asctime)s): %(message)s" if log_format is None else log_format
    for log_handler in logger.handlers:
        logger.removeHandler(log_handler)
    for log_filter in logger.filters:
        logger.removeFilter(log_filter)
    logging.basicConfig(level=log_level, format=log_format)


class Cat (enum.Enum):
    """
    A class to include any categorical infromation to be used intead of
    string or other data types comparison.
    """
    
    # to define categories of overlap
    PRIME_5 = 1
    PRIME_3 = 2
    INTRON = 3
    DISCARD = 4
    
    # to define strand specificity of RNA library
    AUTO = 5
    FORWARD = 6
    REVERSE = 7
    UNSTRANDED = 8

class IntronOverlaps:
    """
    Class to store counters for each inron overlapped with at least one read.
    The instance of this class can be used as Iterator.
    """

    def __init__(self):
        self.reset()

    def reset(self):
        self.overlaps = {}

    def __getitem__(self, key):                                                    # key is a tuple (contig, start, end, name, strand)
        return self.overlaps.setdefault(key, {"p5": 0, "p3": 0})

    def __iter__(self):
        for (contig, start, end, name, strand), value in self.overlaps.items():
            yield contig, start, end, name, strand, value["p5"], value["p3"]

    def update(self, key, category, step=None):                                    # key is a tuple (contig, start, end, name, strand)
        """
        Increments 'p5' or 'p3' overlap counter by 'step' value.
        Any other categories except PRIME_5 and PRIME_3 are ignored
        """

        step = 1 if step is None else step
        if category is Cat.PRIME_5:
            self[key]["p5"] += step
        elif category is Cat.PRIME_3:
            self[key]["p3"] += step


class Counter:

    def __init__(self, bam, ref, span, strandness, threads=None):
        self.bam = bam
        self.ref = ref
        self.span = span
        self.strandness = strandness
        self.threads = 1 if threads is None else threads
        self.paired = self.is_paired()
        self.reset()

    def reset(self):
        self.cache = {}
        self.overlaps = IntronOverlaps()

    def get_overlap_category(self, r_start, r_end, i_start, i_end, span=None):
        span = self.span if span is None else span
        if r_end - i_start >= span and i_start - r_start >= span:
            return Cat.PRIME_5
        elif r_start - i_start >= span and i_end - r_end >= span:
            return Cat.INTRON
        elif r_end - i_end >= span and i_end - r_start >= span:
            return Cat.PRIME_3
        else:
            return Cat.DISCARD

    def is_paired(self, n_reads=None):
        n_reads = 20 if n_reads is None else n_reads
        with pysam.AlignmentFile(self.bam, mode="rb", threads=self.threads) as bam_handler:
            bam_iter = bam_handler.fetch()
            for i in range(n_reads):
                try:
                    bam_handler.mate(next(bam_iter))
                except ValueError:
                    return False
            return True

    def guard_strandness(function):
        def check(self, current_data, cached_data=None):
            if self.strandness is Cat.AUTO:
                if current_data[4] == current_data[8]:
                    if cached_data is None or cached_data[4] == cached_data[8]:
                        function(self, current_data, cached_data)
            elif self.strandness is Cat.UNSTRANDED:
                function(self, current_data, cached_data)
        return check

    def guard_distance(function):
        def __check(self, current_data, cached_data=None):
            if cached_data is None or current_data[0:5] == cached_data[0:5]:            # make sure we didn't accidentally jumped to the next intron
                function(self, current_data, cached_data)
        return __check
    
    @guard_distance
    @guard_strandness
    def update_overlaps(self, current_data, cached_data=None):
        current_category = self.get_overlap_category(
            r_start=current_data[5],
            r_end=current_data[6],
            i_start=current_data[1],
            i_end=current_data[2]
        )
        if cached_data is None:                                                         # working with single read data as we didn't store any cache
            self.overlaps.update(current_data[0:5], current_category)
        else:
            cached_category = self.get_overlap_category(
                r_start=cached_data[5],
                r_end=cached_data[6],
                i_start=cached_data[1],
                i_end=cached_data[2]
            )
            if cached_category is Cat.INTRON:
                self.overlaps.update(
                    current_data[0:5],
                    current_category
                )
            elif current_category is Cat.INTRON:
                self.overlaps.update(
                    cached_data[0:5],
                    cached_category
                )

    def skip_read(self, read):
        return read.is_secondary or \
               read.is_duplicate or \
               read.is_supplementary or \
               (read.is_paired and read.mate_is_unmapped)

    def get_correct_contig(self, contig, handler):
        try:
            handler.fetch(contig)
            return contig
        except ValueError:
            return contig.lstrip("chr")

    def calculate(self, contig):
        with pysam.AlignmentFile(self.bam, mode="rb", threads=self.threads) as bam_handler:
            with pysam.TabixFile(self.ref, mode="r", parser=pysam.asBed(), threads=self.threads) as ref_handler:
                contig_ref = self.get_correct_contig(contig, ref_handler)                                           # contig in the file can be both with or without 'chr' prefix
                contig_bam = self.get_correct_contig(contig, bam_handler)                                           # the same as above
                intron_iter = ref_handler.fetch(contig_ref)
                intron = next(intron_iter)                                                                          # get initial value from intron iterator
                no_introns = False                                                                                  # to break the outer fetching reads loop
                for read in bam_handler.fetch(contig_bam):                                                          # fetches only mapped reads
                    if self.skip_read(read):                                                                        # gate to skip all "bad" reads
                        continue
                    while read.reference_start - intron.end >= 0:
                        try:
                            intron = next(intron_iter)
                        except StopIteration:
                            no_introns = True
                            break
                    if no_introns:                                                                              # no need to iterate over the reads if no introns left
                        break
                    transcript_strand = None
                    try:
                        transcript_strand = read.get_tag("XS")
                    except KeyError:
                        pass
                    collected_data = (                                                                          # tuple takes the smallest amount of memory
                        contig,
                        intron.start,
                        intron.end,
                        intron.name,
                        intron.strand,
                        read.reference_start,
                        read.reference_end,
                        "-" if read.is_reverse else "+",                                                        # to which strand the read was mapped
                        transcript_strand
                    )
                    try:
                        cached_data = self.cache[read.query_name]
                        self.update_overlaps(collected_data, cached_data)
                        del self.cache[read.query_name]
                    except KeyError:
                        if self.paired:
                            self.cache[read.query_name] = collected_data
                        else:
                            self.update_overlaps(collected_data)


    def export(self, location):
        logging.info(f"Export results to {location}")
        with open(location, "w") as out_handler:
            for contig, start, end, name, strand, p5, p3 in self.overlaps:
                out_handler.write(f"{contig}\t{start-self.span-1}\t{start+self.span-1}\t{name}-{start}\t{p5}\t{strand}\n")
                out_handler.write(f"{contig}\t{end-self.span}\t{end+self.span}\t{name}-{end}\t{p3}\t{strand}\n")


def guard_chr(function):
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
    multiprocessing.current_process().name = job[0]                          # mostly for logging purposes
    setup_logger(multiprocessing.get_logger(), args.loglevel)
    logging.info(f"Process chromosome {job[0]}")
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
            logging.info(f"Collect results from {chunk}")
            with open(chunk, "r") as input_stream:
                output_stream.write(input_stream.read())
                os.remove(chunk)


@guard_chr
def get_all_bam_chr(args):
    with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as bam_handler:
        return [_.contig for _ in bam_handler.get_index_statistics()]


@guard_chr
def get_all_ref_chr(args):
    with pysam.TabixFile(args.ref, mode="r", parser=pysam.asBed(), threads=args.threads) as ref_handler:
        return ref_handler.contigs


def assert_args(args):
    args.chr = get_all_bam_chr(args) if len(args.chr) == 0 else [c if c.startswith("chr") else f"chr{c}" for c in args.chr]
    args.strandness = {
        "auto": Cat.AUTO,
        "forward": Cat.FORWARD,
        "reverse": Cat.REVERSE,
        "unstranded": Cat.UNSTRANDED
    }[args.strandness]
    args.loglevel = {
        "fatal": logging.FATAL,
        "error": logging.ERROR,
        "warning": logging.WARNING,
        "info": logging.INFO,
        "debug": logging.DEBUG
    }[args.loglevel]
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
                "Strand specificty of the RNA library."
                "Default: first 'auto' (try to detect strand from the XS tag",
                "of the read), then downgrade to 'unstranded'"
            ]
        ),
        type=str,
        default="auto",
        choices=["auto", "forward", "reverse", "unstranded"]
    )
    general_parser.add_argument("--threads",  help="Number of threads to decompress BAM file", type=int, default=1)
    general_parser.add_argument("--cpus",     help="Number of processes to run in parallel", type=int, default=1)
    general_parser.add_argument("--chr",      help="Select chromosomes to process. Default: all available", type=str, nargs="*", default=[])
    general_parser.add_argument("--loglevel",
            help="Logging level. Default: info",
            type=str,
            default="info",
            choices=["fatal", "error", "warning", "info", "debug"]
    )
    general_parser.add_argument("--output",   help="Output file prefix", type=str, default="intron")
    return general_parser


def main(argsl=None):
    if argsl is None:
        argsl = sys.argv[1:]
    args, _ = arg_parser().parse_known_args(argsl)
    args = assert_args(normalize_args(args, ["span", "threads", "cpus", "chr", "strandness", "loglevel"]))
    setup_logger(logging.root, args.loglevel)
    start = time.time()
    jobs = get_jobs(args)
    logging.info(f"Span {len(jobs)} job(s) among a pool of size {args.cpus}")
    with multiprocessing.Pool(args.cpus) as pool:
        pool.map(partial(process_contig, args), jobs)
    collect_results(args, jobs)
    logging.info (f"""Elapsed time: {round(time.time() - start)} sec""")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
