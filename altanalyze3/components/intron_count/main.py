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


### Reference-consuming CIGAR operations: M, D, N, =, X
REF_CONSUMING = frozenset({0, 2, 3, 7, 8})


def aligned_blocks(read):
    """Return the contiguous reference intervals the read actually aligns over.

    Intron retention asks where a read has aligned bases, so the alignment is split at every
    skipped span. A read that splices an intron out has no aligned bases inside it and simply
    does not overlap it.
    """
    blocks = []
    pos = read.reference_start
    block_start = pos
    for op, length in read.cigartuples:
        if op == 3:                                   # the read skips this reference span
            if pos > block_start:
                blocks.append((block_start, pos))
            pos += length
            block_start = pos
        elif op in REF_CONSUMING:
            pos += length
    if pos > block_start:
        blocks.append((block_start, pos))
    return blocks


def covering_extent(blocks, intron_start, intron_end):
    """Return the aligned block that overlaps this intron, as ``(start, end)``, or None.

    The alignment envelope of a long read spans many genes, so the envelope cannot say whether
    one intron is retained. The block that touches the intron can.
    """
    best = None
    for block_start, block_end in blocks:
        if block_end > intron_start and block_start < intron_end:
            overlap = min(block_end, intron_end) - max(block_start, intron_start)
            if best is None or overlap > best[0]:
                best = (overlap, block_start, block_end)
    if best is None:
        return None
    return best[1], best[2]


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
            ### A read lying wholly inside an intron. update_overlaps returns this category but
            ### counts no boundary evidence for it, matching the prior behaviour. The key was
            ### missing, which only stayed hidden while one intron was examined per read.
            IntRetCat.INTRON: [],
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

    def _raw_data(self, contig, intron, read_info, extent):
        """Build one IntRetRawData for a single read against a single intron."""
        return IntRetRawData(
            contig = contig,
            intron_start = intron[0],
            intron_end = intron[1],
            intron_name = intron[2],
            intron_strand = intron[3],
            read_start = extent[0],
            read_end = extent[1],
            read_strand = read_info["read_strand"],
            xs_strand = read_info["xs_strand"],
            read_name = read_info["read_name"],
            read_1 = read_info["read_1"],
            read_2 = read_info["read_2"],
        )

    def calculate(self, contig):
        """Count intron-retention evidence for EVERY intron each read touches.

        The previous version advanced one sequential intron iterator, so exactly one intron was
        examined per read, and it judged overlap from the alignment envelope. A long read spans
        many introns and skips most of them, so both had to change. Each read is now scored
        against every intron it overlaps, independently, from its aligned blocks.
        """
        reads_seen = 0
        pairs_scored = 0
        with pysam.AlignmentFile(self.bam, mode="rb", threads=self.threads) as bam_handler:
            with pysam.TabixFile(str(self.ref), mode="r", parser=pysam.asBed(), threads=self.threads) as ref_handler:
                contig_ref = get_correct_contig(contig, ref_handler)                                                # contig in the file can be both with or without 'chr' prefix
                contig_bam = get_correct_contig(contig, bam_handler)                                                # the same as above
                for read in bam_handler.fetch(contig_bam):                                                          # fetches only mapped reads
                    if skip_bam_read(read):                                                                        # gate to skip all "bad" reads
                        continue
                    if read.cigartuples is None:
                        continue
                    reads_seen += 1
                    blocks = aligned_blocks(read)
                    xs_strand = None
                    try:
                        xs_strand = read.get_tag("XS")
                    except KeyError:
                        pass
                    read_info = {
                        "read_strand": "-" if read.is_reverse else "+",
                        "xs_strand": xs_strand,
                        "read_name": read.query_name,
                        "read_1": read.is_read1,
                        "read_2": read.is_read2,
                    }
                    ### Every intron this alignment reaches, not just the next one in the file.
                    try:
                        candidates = list(ref_handler.fetch(contig_ref, read.reference_start, read.reference_end))
                    except ValueError:                                                                              # contig absent from the reference
                        continue
                    for entry in candidates:
                        intron = (entry.start, entry.end, entry.name, entry.strand)
                        extent = covering_extent(blocks, entry.start, entry.end)
                        if extent is None:
                            continue
                        pairs_scored += 1
                        current_data = self._raw_data(contig, intron, read_info, extent)
                        if self.paired:
                            key = (read.query_name, entry.start, entry.end, entry.name)
                            if key in self.cache:
                                cached_data = self.cache.pop(key)
                                overlapped_as = self.update_overlaps(current_data, cached_data)
                                self.used_reads[overlapped_as].append(read.query_name)
                            else:
                                self.cache[key] = current_data
                        else:
                            overlapped_as = self.update_overlaps(current_data)
                            self.used_reads[overlapped_as].append(read.query_name)
        logging.info(f"{contig}: {reads_seen} reads; {pairs_scored} read-intron pairs scored")

    def export_counts(self):
        logging.info(f"""Save counts to {self.location}""")
        with self.location.open("w") as out_handler:
            for contig, start, end, name, strand, p5, p3 in self.overlaps:
                if p5>0 and p3>0: # PE reads exist at the 5' and 3' end of the intron
                    if strand == "+":
                        out_handler.write(f"{contig}\t{start-1}\t{start}\t{name}_{start}\t{p5}\t{strand}\n")
                        out_handler.write(f"{contig}\t{end}\t{end+1}\t{name}_{end}\t{p3}\t{strand}\n")
                    elif strand == "-":
                        out_handler.write(f"{contig}\t{end+1}\t{end}\t{name}_{end}\t{p5}\t{strand}\n")
                        out_handler.write(f"{contig}\t{start}\t{start+1}\t{name}_{start+1}\t{p3}\t{strand}\n")


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
            location=args.tmp.joinpath(args.bam.stem, get_tmp_suffix())
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
            logging.info(f"Collect counts from {job.location}")
            with job.location.open("r") as input_stream:
                output_stream.write(input_stream.read())
            # Only delete now if NOT saving reads
            if not args.savereads:
                logging.debug(f"Remove {job.location}")
                job.location.unlink()

    if args.savereads:
        tmp_bam = args.output.with_suffix(get_tmp_suffix() + ".bam")
        with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as template_handler:
            with pysam.AlignmentFile(tmp_bam, "wb", threads=args.threads, template=template_handler) as out_bam_handler:
                for job in jobs:
                    bam_location = job.location.with_suffix(".bam")
                    logging.info(f"Collect reads from {bam_location}")
                    with pysam.AlignmentFile(bam_location, mode="rb", threads=args.threads) as in_bam_handler:
                        for read in in_bam_handler.fetch(until_eof=True):
                            out_bam_handler.write(read)
                        logging.debug(f"Remove {bam_location}")
                        bam_location.unlink()  # delete bam parts here
        pysam.sort("-o", str(args.output.with_suffix(".bam")), str(tmp_bam))
        pysam.index(str(args.output.with_suffix(".bam")))
        tmp_bam.unlink()


def count_introns(args):
    sample_path = args.tmp.joinpath(args.bam.stem)
    sample_path.mkdir(parents=True, exist_ok=True)
    jobs = get_jobs(args)
    logging.info(f"""Span {len(jobs)} job(s) between {args.cpus} CPU(s)""")
    with multiprocessing.Pool(args.cpus) as pool:
        pool.map(partial(process_contig, args), jobs)
    collect_results(args, jobs)

    #"""
    import pathlib
    # Ensure you are getting the .bed.gz file path (this is args.ref after indexing)
    ref_gz = pathlib.Path(args.ref)
    # Delete .bed.gz (if you want) and .bed.gz.tbi
    files_to_delete = [ref_gz, ref_gz.with_name(ref_gz.name + ".tbi")]
    for f in files_to_delete:
        if f.exists():
            logging.info(f"Removing temporary file: {f}")
            f.unlink()
    #"""
    logging.debug(f"Removing temporary directory and all contents for sample {args.bam.stem} at {sample_path}")
    shutil.rmtree(sample_path)

if __name__ == "__main__":
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description="Count intron overlaps from BAM using GFF/GTF annotation.")

    parser.add_argument("--bam", type=Path, required=False, default=Path("/path/to/input.bam"), help="Input BAM file")
    parser.add_argument("--ref", type=Path, required=False, default=Path("/path/to/annotation.bed.gz"), help="Input intron BED file (.gz)")
    parser.add_argument("--output", type=Path, required=False, default=Path("/path/to/output_prefix"), help="Prefix path for output files (no extension)")
    parser.add_argument("--tmp", type=Path, required=False, default=Path("/tmp/altanalyze_tmp"), help="Temporary directory for intermediate files")
    parser.add_argument("--span", type=int, default=10, help="Span threshold for 5'/3' prime assignment")
    parser.add_argument("--strandness", type=int, default=0, help="0: Auto, 1: Unstranded, 2: Forward, 3: Reverse")
    parser.add_argument("--cpus", type=int, default=4, help="Number of CPU processes to use")
    parser.add_argument("--threads", type=int, default=2, help="Number of BAM/Tabix threads per process")
    parser.add_argument("--savereads", action="store_true", help="Save classified reads in output BAM")
    parser.add_argument("--chr", nargs='+', default=[], help="Specific chromosomes to process (e.g., --chr chr1 chr2 chrX)")
    parser.add_argument("--loglevel", type=str, default="INFO", help="Logging level: DEBUG, INFO, WARNING, ERROR")

    args = parser.parse_args()
    """
    # --- DIRECT OVERRIDE: you want the script runnable without needing manual input ---
    args.bam = Path("/Users/your_username/path_to_input.bam")
    args.ref = Path("/Users/your_username/path_to_introns.bed.gz")
    args.output = Path("/Users/your_username/path_to_output_prefix")
    args.tmp = Path("/Users/your_username/tmp_folder")
    args.chr = ["chr1", "chr2", "chrX"]  # Example: change to [] if you want all chromosomes
    args.savereads = True  # Set to False if you don't want output BAM
    """
    args.chr =["chr1", "chr2", "chr6"] 
    # Ensure tmp directory exists
    args.tmp.mkdir(parents=True, exist_ok=True)
    print (args)
    count_introns(args)
