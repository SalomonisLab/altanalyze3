import pysam
import shutil
import logging
import multiprocessing
from functools import partial

from altanalyze3.utilities.logger import setup_logger
from altanalyze3.utilities.constants import (
    IntRetCat,
    IntRetRawData,
    Job
)
from altanalyze3.utilities.helpers import get_tmp_suffix
from altanalyze3.utilities.io import (
    get_all_bam_chr,
    get_all_ref_chr,
    get_correct_contig,
    skip_bam_read,
    is_bam_paired
)


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
        self.paired = is_bam_paired(self.bam, self.threads)
        self.reset()

    def reset(self):
        self.cache = {}
        self.overlaps = IntronOverlaps()
        self.used_reads = {
            IntRetCat.PRIME_5: [],
            IntRetCat.PRIME_3: [],
            IntRetCat.DISCARD: []
        }

    def get_overlap_category(self, raw_data, span=None):
        span = self.span if span is None else span
        if raw_data.read_end - raw_data.intron_start >= span and \
                raw_data.intron_start - raw_data.read_start >= span:
            return IntRetCat.PRIME_5
        elif raw_data.read_start - raw_data.intron_start >= 0 and \
                raw_data.intron_end - raw_data.read_end >= 0:
            return IntRetCat.INTRON
        elif raw_data.read_end - raw_data.intron_end >= span and \
                raw_data.intron_end - raw_data.read_start >= span:
            return IntRetCat.PRIME_3
        else:
            return IntRetCat.DISCARD

    def guard_strandness(function):
        def check(self, current_data, cached_data=None):
            if self.strandness is IntRetCat.AUTO:
                if self.paired and \
                    (
                        current_data.xs_strand == cached_data.xs_strand == None or                                                    # downgrade to "unstranded"
                        current_data.xs_strand == cached_data.xs_strand == current_data.intron_strand == cached_data.intron_strand    # pass strandness check
                    ):
                    return function(self, current_data, cached_data)
                elif not self.paired and \
                    (
                        current_data.xs_strand is None or                                                                             # downgrade to "unstranded"
                        current_data.xs_strand == current_data.intron_strand                                                          # pass strandness check
                    ):
                    return function(self, current_data)
                else:
                    logging.debug("Strandness guard blocked the overlap for")
                    logging.debug(f"""{current_data}""")
                    return IntRetCat.DISCARD
            elif self.strandness is IntRetCat.UNSTRANDED:
                if self.paired:
                    return function(self, current_data, cached_data)
                else:
                    return function(self, current_data)
            elif self.strandness is IntRetCat.FORWARD:
                if self.paired and current_data.intron_strand == "+" and \
                    current_data.read_1 and current_data.read_strand == "+" and \
                        cached_data.read_2 and cached_data.read_strand == "-":
                    return function(self, current_data, cached_data)
                elif not self.paired and current_data.intron_strand == "+" and \
                    current_data.read_strand == "+":                                                           # all of the reads came from the "+" strand
                    return function(self, current_data)
                else:
                    logging.debug("Strandness guard blocked the overlap for")
                    logging.debug(f"""{current_data}""")
                    return IntRetCat.DISCARD
            elif self.strandness is IntRetCat.REVERSE:
                if self.paired and current_data.intron_strand == "-" and \
                    current_data.read_1 and current_data.read_strand == "-" and \
                        cached_data.read_2 and cached_data.read_strand == "+":
                    return function(self, current_data, cached_data)
                elif not self.paired and current_data.intron_strand == "-" and \
                    current_data.read_strand == "-":                                                           # all of the reads came from the "-" strand
                    return function(self, current_data)
                else:
                    logging.debug("Strandness guard blocked the overlap for")
                    logging.debug(f"""{current_data}""")
                    return IntRetCat.DISCARD
        return check

    def guard_distance(function):
        def check(self, current_data, cached_data=None):
            if not self.paired:
                return function(self, current_data)
            elif self.paired and current_data[0:5] == cached_data[0:5]:
                return function(self, current_data, cached_data)
            else:
                logging.debug("Distance guard blocked the overlap for")
                logging.debug(f"""{current_data}""")
                logging.debug(f"""{cached_data}""")
                return IntRetCat.DISCARD
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
        if not self.paired:
            logging.debug("Check overlap for single read")
            logging.debug(f"""{current_data}, {current_category}""")
            if current_category is IntRetCat.PRIME_5:
                self.overlaps.increment_p5(intron_key)
            elif current_category is IntRetCat.PRIME_3:
                self.overlaps.increment_p3(intron_key)
            return current_category
        else:
            assert(cached_data != None)
            logging.debug("Check overlap for paired-end")
            cached_category = self.get_overlap_category(cached_data)
            logging.debug(f"""{current_data}, {current_category}""")
            logging.debug(f"""{cached_data}, {cached_category}""")
            if IntRetCat.DISCARD in [current_category, cached_category] or \
                    current_category == cached_category:                                                           # both introns, both 5', or both 3'
                return IntRetCat.DISCARD
            elif current_category in [IntRetCat.PRIME_5, IntRetCat.INTRON] and cached_category in [IntRetCat.PRIME_5, IntRetCat.INTRON]:   # 5' and intron or intron and 5'
                self.overlaps.increment_p5(intron_key)
                return IntRetCat.PRIME_5
            elif current_category in [IntRetCat.PRIME_3, IntRetCat.INTRON] and cached_category in [IntRetCat.PRIME_3, IntRetCat.INTRON]:   # 3' and intron or intron and 3'
                self.overlaps.increment_p3(intron_key)
                return IntRetCat.PRIME_3
            else:
                logging.debug(f"""Not implemented combination of {current_category} and {cached_category} categories""")
                return IntRetCat.DISCARD

    def calculate(self, contig):
        with pysam.AlignmentFile(self.bam, mode="rb", threads=self.threads) as bam_handler:
            with pysam.TabixFile(str(self.ref), mode="r", parser=pysam.asBed(), threads=self.threads) as ref_handler:
                contig_ref = get_correct_contig(contig, ref_handler)                                                # contig in the file can be both with or without 'chr' prefix
                contig_bam = get_correct_contig(contig, bam_handler)                                                # the same as above
                intron_iter = ref_handler.fetch(contig_ref)
                intron = next(intron_iter)                                                                          # get initial value from intron iterator
                no_introns = False                                                                                  # to break the outer fetching reads loop
                for read in bam_handler.fetch(contig_bam):                                                          # fetches only mapped reads
                    logging.debug(
                        f"""Fetch a read {read.query_name} {contig_bam}:{read.reference_start}-{read.reference_end} {"-" if read.is_reverse else "+"}"""
                    )
                    if skip_bam_read(read):                                                                        # gate to skip all "bad" reads
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
                    current_data = IntRetRawData(
                        contig = contig,
                        intron_start = intron.start,
                        intron_end = intron.end,
                        intron_name = intron.name,
                        intron_strand = intron.strand,
                        read_start = read.reference_start,
                        read_end = read.reference_end,
                        read_strand = "-" if read.is_reverse else "+",
                        xs_strand = xs_strand,
                        read_name = read.query_name,
                        read_1 = read.is_read1,
                        read_2 = read.is_read2,
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
        with self.location.open("w") as out_handler:
            for contig, start, end, name, strand, p5, p3 in self.overlaps:
                out_handler.write(f"{contig}\t{start-self.span-1}\t{start+self.span-1}\t{name}_{start}\t{p5}\t{strand}\n")
                out_handler.write(f"{contig}\t{end-self.span}\t{end+self.span}\t{name}_{end}\t{p3}\t{strand}\n")

    def export_reads(self):
        bam_location = self.location.with_suffix(".bam")
        logging.info(f"""Save reads to {bam_location}""")
        with pysam.AlignmentFile(self.bam, mode="rb", threads=self.threads) as in_bam_handler:
            with pysam.AlignmentFile(bam_location, "wb", threads=self.threads, template=in_bam_handler) as out_bam_handler:
                for read in in_bam_handler.fetch():
                    logging.debug(f"""Fetch read {read.query_name}""")
                    if skip_bam_read(read):
                        logging.debug("Skip")
                        continue
                    if read.query_name in self.used_reads[IntRetCat.PRIME_5]:
                        read.set_tag("XI", "P5")
                        logging.debug("Assign XI=P5")
                    elif read.query_name in self.used_reads[IntRetCat.PRIME_3]:
                        read.set_tag("XI", "P3")
                        logging.debug("Assign XI=P3")
                    elif read.query_name in self.used_reads[IntRetCat.DISCARD]:
                        read.set_tag("XI", "D")
                        logging.debug("Assign XI=D")
                    else:
                        read.set_tag("XI", "U")
                        logging.debug("Assign XI=U")
                    out_bam_handler.write(read)


def get_jobs(args):
    return [
        Job(
            contig=contig,                                                            # contig is always prepended with 'chr'
            location=args.tmp.joinpath(get_tmp_suffix())
        )
        for contig in get_all_bam_chr(args.bam, args.threads)
            if contig in get_all_ref_chr(args.ref, args.threads) and contig in args.chr                 # safety measure to include only chromosomes present in BAM, BED, and --chr
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
    with args.output.with_suffix(".bed").open("w") as output_stream:
        for job in jobs:
            logging.info(f"""Collect counts from {job.location}""")
            with job.location.open("r") as input_stream:
                output_stream.write(input_stream.read())
                logging.debug(f"""Remove {job.location}""")
                job.location.unlink()
    if args.savereads:
        tmp_bam = args.output.with_suffix(get_tmp_suffix() + ".bam")
        with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as template_handler:
            with pysam.AlignmentFile(tmp_bam, "wb", threads=args.threads, template=template_handler) as out_bam_handler:
                for job in jobs:
                    bam_location = job.location.with_suffix(".bam")
                    logging.info(f"""Collect reads from {bam_location}""")
                    with pysam.AlignmentFile(bam_location, mode="rb", threads=args.threads) as in_bam_handler:
                        for read in in_bam_handler.fetch(until_eof=True):                                         # we don't need index because of until_eof
                            out_bam_handler.write(read)
                        logging.debug(f"""Remove {bam_location}""")
                        bam_location.unlink()
        pysam.sort("-o", str(args.output.with_suffix(".bam")), str(tmp_bam))
        pysam.index(str(args.output.with_suffix(".bam")))
        tmp_bam.unlink()


def count_introns(args):
    jobs = get_jobs(args)
    logging.info(f"""Span {len(jobs)} job(s) between {args.cpus} CPU(s)""")
    with multiprocessing.Pool(args.cpus) as pool:
        pool.map(partial(process_contig, args), jobs)
    collect_results(args, jobs)

    logging.debug("Removing temporary directory")
    shutil.rmtree(args.tmp)