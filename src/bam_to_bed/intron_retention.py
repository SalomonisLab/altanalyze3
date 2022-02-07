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
from collections import namedtuple


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


Job = namedtuple("Job", "contig location")
RawData = namedtuple(
    "RawData",
    "contig intron_start intron_end intron_name intron_strand read_start read_end read_strand xs_strand read_name"
)


class Cat (enum.IntEnum):
    """
    A class to include any categorical infromation to be used intead of
    string or other data types comparison.
    """
    
    # to define categories of overlap
    PRIME_5 = enum.auto()
    PRIME_3 = enum.auto()
    INTRON = enum.auto()
    DISCARD = enum.auto()
    
    # to define strand specificity of RNA library
    AUTO = enum.auto()
    FORWARD = enum.auto()
    REVERSE = enum.auto()
    UNSTRANDED = enum.auto()

    def __str__(self):
        return self.name


class IntronOverlaps:
    """
    Class to store counters for each inron overlapped with at least one read.
    The instance of this class can be used as Iterator.
    """

    def __init__(self):
        self.reset()

    def reset(self):
        self.overlaps = {}

    def __getitem__(self, key):                                                # key is a tuple (contig, start, end, name, strand)
        return self.overlaps.setdefault(key, {"p5": 0, "p3": 0})

    def __iter__(self):
        for (contig, start, end, name, strand), counters in self.overlaps.items():
            yield contig, start, end, name, strand, counters["p5"], counters["p3"]

    def increment_p5(self, key, step=None):                                    # key is a tuple (contig, start, end, name, strand)
        step = 1 if step is None else step
        self[key]["p5"] += step
        logging.debug(f"""Increment p5 counter on {step} for {key}""")

    def increment_p3(self, key, step=None):                                    # key is a tuple (contig, start, end, name, strand)
        step = 1 if step is None else step
        self[key]["p3"] += step
        logging.debug(f"""Increment p3 counter on {step} for {key}""")


class Counter:

    def __init__(self, bam, ref, span, strandness, location, threads=None):
        self.bam = bam
        self.ref = ref
        self.span = span
        self.strandness = strandness
        self.location = location
        self.threads = 1 if threads is None else threads
        self.paired = self.is_paired()
        self.reset()

    def reset(self):
        self.cache = {}
        self.overlaps = IntronOverlaps()
        self.used_reads = {
            Cat.PRIME_5: [],
            Cat.PRIME_3: [],
            Cat.DISCARD: []
        }

    def get_overlap_category(self, raw_data, span=None):
        span = self.span if span is None else span
        if raw_data.read_end - raw_data.intron_start >= span and \
                raw_data.intron_start - raw_data.read_start >= span:
            return Cat.PRIME_5
        elif raw_data.read_start - raw_data.intron_start >= 0 and \
                raw_data.intron_end - raw_data.read_end >= 0:
            return Cat.INTRON
        elif raw_data.read_end - raw_data.intron_end >= span and \
                raw_data.intron_end - raw_data.read_start >= span:
            return Cat.PRIME_3
        else:
            return Cat.DISCARD

    def is_paired(self):
        with pysam.AlignmentFile(self.bam, mode="rb", threads=self.threads) as bam_handler:
            for read in bam_handler.fetch():                                                 # this fetches only mapped reads
                if self.skip_read(read):
                    continue
                logging.info(f"""Alignments are in {"paired-end" if read.is_paired else "single read"} format""")
                return read.is_paired

    def guard_strandness(function):
        def check(self, current_data, cached_data=None):
            if self.strandness is Cat.AUTO:
                if current_data[4] == current_data[8]:
                    if cached_data is None or cached_data[4] == cached_data[8]:
                        return function(self, current_data, cached_data)
                    else:
                        logging.debug(f"""Strandness guard blocked the overlap for""")
                        logging.debug(f"""{current_data}""")
                        logging.debug(f"""{cached_data}""")
                        return Cat.DISCARD
                else:
                    logging.debug(f"""Strandness guard blocked the overlap for""")
                    logging.debug(f"""{current_data}""")
                    return Cat.DISCARD
            elif self.strandness is Cat.UNSTRANDED:
                return function(self, current_data, cached_data)
        return check

    def guard_distance(function):
        def check(self, current_data, cached_data=None):
            if cached_data is None or current_data[0:5] == cached_data[0:5]:            # make sure we didn't accidentally jumped to the next intron
                return function(self, current_data, cached_data)
            else:
                logging.debug(f"""Distance guard blocked the overlap for""")
                logging.debug(f"""{current_data}""")
                logging.debug(f"""{cached_data}""")
                return Cat.DISCARD
        return check
    
    @guard_distance
    @guard_strandness
    def update_overlaps(self, current_data, cached_data=None):
        current_category = self.get_overlap_category(current_data)
        intron_key = (                                                  # will be the same for both current_data and cached_data because of guard_distance
            current_data.contig,
            current_data.intron_start,
            current_data.intron_end,
            current_data.intron_name,
            current_data.intron_strand
        )
        if cached_data is None:
            logging.debug("Check overlap for single read")
            logging.debug(f"""{current_data}, {current_category}""")
            if current_category is Cat.PRIME_5:
                self.overlaps.increment_p5(intron_key)
            elif current_category is Cat.PRIME_3:
                self.overlaps.increment_p3(intron_key)
            return current_category
        else:
            logging.debug("Check overlap for paired-end")
            cached_category = self.get_overlap_category(cached_data)
            logging.debug(f"""{current_data}, {current_category}""")
            logging.debug(f"""{cached_data}, {cached_category}""")
            if Cat.DISCARD in [current_category, cached_category] or \
                    current_category == cached_category == Cat.INTRON:
                return Cat.DISCARD
            elif current_category in [Cat.PRIME_5, Cat.INTRON] and cached_category in [Cat.PRIME_5, Cat.INTRON]:
                self.overlaps.increment_p5(intron_key)
                return Cat.PRIME_5
            elif current_category in [Cat.PRIME_3, Cat.INTRON] and cached_category in [Cat.PRIME_3, Cat.INTRON]:
                self.overlaps.increment_p3(intron_key)
                return Cat.PRIME_3
            else:
                logging.debug(f"""Not implemented combination of {current_category} and {cached_category} categories""")
                return Cat.DISCARD

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
                    logging.debug(
                        f"""Fetch a read {read.query_name} {contig_bam}:{read.reference_start}-{read.reference_end} {"-" if read.is_reverse else "+"}"""
                    )
                    if self.skip_read(read):                                                                        # gate to skip all "bad" reads
                        logging.debug(f"""Skip a read {read.query_name}""")
                        logging.debug(f"""is_secondary: {read.is_secondary}""")
                        logging.debug(f"""is_duplicate: {read.is_duplicate}""")
                        logging.debug(f"""is_supplementary: {read.is_supplementary}""")
                        logging.debug(f"""is_paired and mate_is_unmapped: {read.is_paired and read.mate_is_unmapped}""")
                        continue
                    while read.reference_start - intron.end >= 0:
                        try:
                            intron = next(intron_iter)
                            logging.debug(
                                f"""Switched to a new intron {intron.name} {contig_ref}:{intron.start}-{intron.end} {intron.strand}"""
                            )
                        except StopIteration:
                            no_introns = True
                            break
                    if no_introns:                                                                              # no need to iterate over the reads if no introns left
                        logging.debug("Halt read iteration - run out of introns")
                        break
                    xs_strand = None
                    try:
                        xs_strand = read.get_tag("XS")
                    except KeyError:
                        pass
                    current_data = RawData(
                        contig = contig,
                        intron_start = intron.start,
                        intron_end = intron.end,
                        intron_name = intron.name,
                        intron_strand = intron.strand,
                        read_start = read.reference_start,
                        read_end = read.reference_end,
                        read_strand = "-" if read.is_reverse else "+",
                        xs_strand = xs_strand,
                        read_name = read.query_name
                    )
                    if self.paired:
                        if read.query_name in self.cache:
                            cached_data = self.cache[read.query_name]
                            overlapped_as = self.update_overlaps(current_data, cached_data)
                            self.used_reads[overlapped_as].append(read.query_name)
                            del self.cache[read.query_name]
                            logging.debug(f"""Remove cached data for {read.query_name}""")
                        else:
                            self.cache[read.query_name] = current_data
                            logging.debug(f"""Add cached data for {read.query_name}""")
                    else:
                        overlapped_as = self.update_overlaps(current_data)
                        self.used_reads[overlapped_as].append(read.query_name)

    def export_counts(self):
        logging.info(f"""Save counts to {self.location}""")
        with open(self.location, "w") as out_handler:
            for contig, start, end, name, strand, p5, p3 in self.overlaps:
                out_handler.write(f"{contig}\t{start-self.span-1}\t{start+self.span-1}\t{name}-{start}\t{p5}\t{strand}\n")
                out_handler.write(f"{contig}\t{end-self.span}\t{end+self.span}\t{name}-{end}\t{p3}\t{strand}\n")

    def export_reads(self):
        bam_location = os.path.splitext(self.location)[0] + ".bam"
        logging.info(f"""Save reads to {bam_location}""")
        with pysam.AlignmentFile(self.bam, mode="rb", threads=self.threads) as in_bam_handler:
            with pysam.AlignmentFile(bam_location, "wb", threads=self.threads, template=in_bam_handler) as out_bam_handler:
                for read in in_bam_handler.fetch():
                    logging.debug(f"""Fetch read {read.query_name}""")
                    if self.skip_read(read):
                        logging.debug("Skip")
                        continue
                    if read.query_name in self.used_reads[Cat.PRIME_5]:
                        read.set_tag("XI", "P5")
                        logging.debug("Assign XI=P5")
                    elif read.query_name in self.used_reads[Cat.PRIME_3]:
                        read.set_tag("XI", "P3")
                        logging.debug("Assign XI=P3")
                    elif read.query_name in self.used_reads[Cat.DISCARD]:
                        read.set_tag("XI", "D")
                        logging.debug("Assign XI=D")
                    else:
                        read.set_tag("XI", "U")
                        logging.debug("Assign XI=U")
                    out_bam_handler.write(read)


class ArgsParser():

    def __init__(self, args):
        self.args, _ = self.get_parser().parse_known_args(args)
        self.normalize_args(["span", "threads", "cpus", "chr", "strandness", "loglevel", "savereads"])
        self.assert_args()
        self.set_args_as_attributes()
        logging.debug("Use argument")
        logging.debug(self.args)

    def set_args_as_attributes(self):
        for arg, value in vars(self.args).items():
            setattr(self, arg, value)

    def get_parser(self):
        general_parser = argparse.ArgumentParser()
        general_parser.add_argument("--bam",      help="Path to the coordinate-sorted indexed BAM file (should include only mapped reads)", type=str, required=True)
        general_parser.add_argument("--ref",      help="Path to the coordinate-sorted indexed gene model reference BED file", type=str, required=True)
        general_parser.add_argument("--span",     help="5' and 3' overlap that read should have over a splice-site to be counted", type=int, default=10)
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
        general_parser.add_argument("--threads",   help="Number of threads to decompress BAM file", type=int, default=1)
        general_parser.add_argument("--cpus",      help="Number of processes to run in parallel", type=int, default=1)
        general_parser.add_argument("--chr",       help="Select chromosomes to process. Default: all available", type=str, nargs="*", default=[])
        general_parser.add_argument("--loglevel",
                help="Logging level. Default: info",
                type=str,
                default="info",
                choices=["fatal", "error", "warning", "info", "debug"]
        )
        general_parser.add_argument("--savereads", help="Export processed reads into the BAM file. Default: False", action="store_true")
        general_parser.add_argument("--output",    help="Output file prefix", type=str, default="intron")
        return general_parser

    def normalize_args(self, skip=None):
        skip = [] if skip is None else skip
        normalized_args = {}
        for key,value in self.args.__dict__.items():
            if key not in skip:
                normalized_args[key] = value if not value or os.path.isabs(value) else os.path.normpath(os.path.join(os.getcwd(), value))
            else:
                normalized_args[key]=value
        self.args = argparse.Namespace (**normalized_args)

    def assert_args(self):
        self.args.chr = get_all_bam_chr(self.args) if len(self.args.chr) == 0 else [c if c.startswith("chr") else f"chr{c}" for c in self.args.chr]
        self.args.strandness = {
            "auto": Cat.AUTO,
            "forward": Cat.FORWARD,
            "reverse": Cat.REVERSE,
            "unstranded": Cat.UNSTRANDED
        }[self.args.strandness]
        self.args.loglevel = {
            "fatal": logging.FATAL,
            "error": logging.ERROR,
            "warning": logging.WARNING,
            "info": logging.INFO,
            "debug": logging.DEBUG
        }[self.args.loglevel]


def setup_logger(logger, log_level, log_format=None):
    log_format = "%(processName)12s (%(asctime)s): %(message)s" if log_format is None else log_format
    for log_handler in logger.handlers:
        logger.removeHandler(log_handler)
    for log_filter in logger.filters:
        logger.removeFilter(log_filter)
    logging.basicConfig(level=log_level, format=log_format)


def get_jobs(args):
    tmp_marker = "".join(random.choices(string.ascii_uppercase + string.digits, k=10))
    return [
        Job(
            contig=contig,                                                      # contig is always prepended with 'chr'
            location=args.output + "__" + contig + "__" + tmp_marker + ".bed"
        )
        for contig in get_all_bam_chr(args)
            if contig in get_all_ref_chr(args) and contig in args.chr           # safety measure to include only chromosomes present in BAM, BED, and --chr
    ]


def process_contig(args, job):
    setup_logger(
        multiprocessing.get_logger(),
        args.loglevel
    )
    multiprocessing.current_process().name = job.contig
    logging.info(f"""Process chromosome {job.contig} to {job.location}""")
    counter = Counter(
        bam=args.bam,
        ref=args.ref,
        span=args.span,
        strandness=args.strandness,
        location=job.location,
        threads=args.threads
    )
    counter.calculate(job[0])
    counter.export_counts()
    if args.savereads:
        counter.export_reads()


def collect_results(args, jobs):
    with open(args.output + ".bed", "w") as output_stream:
        for job in jobs:
            logging.info(f"""Collect counts from {job.location}""")
            with open(job.location, "r") as input_stream:
                output_stream.write(input_stream.read())
                logging.debug(f"""Remove {job.location}""")
                os.remove(job.location)
    if args.savereads:
        tmp_bam = args.output + "".join(random.choices(string.ascii_uppercase + string.digits, k=10)) + ".bam"
        with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as template_handler:
            with pysam.AlignmentFile(tmp_bam, "wb", threads=args.threads, template=template_handler) as out_bam_handler:
                for job in jobs:
                    bam_location = os.path.splitext(job.location)[0] + ".bam"
                    logging.info(f"""Collect reads from {bam_location}""")
                    with pysam.AlignmentFile(bam_location, mode="rb", threads=args.threads) as in_bam_handler:
                        for read in in_bam_handler.fetch(until_eof=True):                                         # we don't need index because of until_eof
                            out_bam_handler.write(read)
                        logging.debug(f"""Remove {bam_location}""")
                        os.remove(bam_location)
        pysam.sort("-o", args.output + ".bam", tmp_bam)
        pysam.index(args.output + ".bam")
        os.remove(tmp_bam)


def guard_chr(function):
    def prefix(c):
        return c if c.startswith("chr") else f"chr{c}"
    def wrapper(contig):
        raw_res = function(contig)
        if isinstance(raw_res, list):
            return [prefix(c) for c in raw_res]
        else:
            return prefix(raw_res)
    return wrapper


@guard_chr
def get_all_bam_chr(args):
    with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as bam_handler:
        return [_.contig for _ in bam_handler.get_index_statistics()]


@guard_chr
def get_all_ref_chr(args):
    with pysam.TabixFile(args.ref, mode="r", parser=pysam.asBed(), threads=args.threads) as ref_handler:
        return ref_handler.contigs


def main(args=None):
    args = ArgsParser(sys.argv[1:] if args is None else args)
    setup_logger(logging.root, args.loglevel)
    start = time.time()
    jobs = get_jobs(args)
    logging.info(f"""Span {len(jobs)} job(s) between {args.cpus} CPU(s)""")
    with multiprocessing.Pool(args.cpus) as pool:
        pool.map(partial(process_contig, args), jobs)
    collect_results(args, jobs)
    logging.info (f"""Elapsed time: {round(time.time() - start)} sec""")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
