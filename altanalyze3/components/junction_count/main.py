import sys
import pysam
import shutil
import logging
import multiprocessing
from functools import partial

from altanalyze3.utilities.io import (
    get_all_bam_chr,
    get_correct_contig
)
from altanalyze3.utilities.logger import setup_logger
from altanalyze3.utilities.constants import Job
from altanalyze3.utilities.helpers import get_tmp_suffix


def process_contig_fast_no_strand(args, job):
    setup_logger(
        multiprocessing.get_logger(),
        args.loglevel
    )
    multiprocessing.current_process().name = job.contig
    logging.info(f"Process chromosome {job.contig} to {job.location}")
    with job.location.open("wt") as output_stream:
        with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as handler:
            introns = handler.find_introns(
                (r for r in handler.fetch(get_correct_contig(job.contig, handler)))
            )
            for position, score in introns.items():
                output_stream.write(
                    f"{job.contig}\t{position[0]}\t{position[1]}\tJUNC:{job.contig}-{position[0]}-{position[1]}\t{score}\t.\n"
                )

### Reference-consuming CIGAR operations: M, D, N, =, X
REF_CONSUMING = frozenset({0, 2, 3, 7, 8})


def resolve_read_strand(read, transcript_oriented):
    """Return the transcript strand for a read, and the source it came from.

    A read is never dropped for lacking a strand call. Order of preference:

    ``XS``  the aligner's own splice-strand call (HISAT2, TopHat, STAR run with
            ``--outSAMstrandField intronMotif``). Authoritative when present.
    ``ts``  minimap2 and pbmm2 transcript-strand tag. Relative to the read, so it is combined
            with the read orientation.
    read    only when the library is transcript-oriented, which holds for PacBio isoseq after
            ``isoseq refine``. Meaningless for unstranded short-read data, so it is gated.
    ``.``   nothing to go on. The junction is still counted and written.
    """
    try:
        return read.get_tag("XS"), "XS"
    except KeyError:
        pass
    try:
        ts = read.get_tag("ts")                     # '+' means the read matches the transcript
        if ts in ("+", "-"):
            flip = read.is_reverse
            fwd = (ts == "+") != flip
            return ("+" if fwd else "-"), "ts"
    except KeyError:
        pass
    if transcript_oriented:
        return ("-" if read.is_reverse else "+"), "read"
    return ".", "none"


def is_transcript_oriented(bam_path, threads):
    """True when each read's orientation carries the transcript direction (PacBio isoseq)."""
    with pysam.AlignmentFile(bam_path, mode="rb", threads=threads) as handler:
        header = handler.header.to_dict()
    platforms = {rg.get("PL", "").upper() for rg in header.get("RG", [])}
    programs = {pg.get("PN", pg.get("ID", "")).lower() for pg in header.get("PG", [])}
    if "PACBIO" in platforms:
        return True
    return any("isoseq" in p or "pbmm2" in p for p in programs)


def process_contig(args, job):
    from collections import defaultdict
    setup_logger(multiprocessing.get_logger(), args.loglevel)
    multiprocessing.current_process().name = job.contig
    logging.info(f"Process chromosome {job.contig} to {job.location}")
    transcript_oriented = getattr(args, "transcript_oriented", None)
    if transcript_oriented is None:
        transcript_oriented = is_transcript_oriented(args.bam, args.threads)
    junctions = defaultdict(int)
    strand_sources = defaultdict(int)
    reads_with_junctions = 0
    with job.location.open("wt") as output_stream:
        with pysam.AlignmentFile(args.bam, mode="rb", threads=args.threads) as bam:
            contig = get_correct_contig(job.contig, bam)
            for read in bam.fetch(contig):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                cigar = read.cigartuples
                if cigar is None:
                    continue
                if not any(op == 3 for op, _ in cigar):
                    continue
                strand, source = resolve_read_strand(read, transcript_oriented)
                strand_sources[source] += 1
                reads_with_junctions += 1
                ### Count EVERY N operation. A long read routinely carries many junctions, and
                ### each one is an independent observation.
                ref_pos = read.reference_start
                for op, length in cigar:
                    if op == 3:
                        junctions[(ref_pos, ref_pos + length, strand)] += 1
                    if op in REF_CONSUMING:
                        ref_pos += length
        for (start, end, strand), count in junctions.items():
            if strand == '-':
                start, end = end, start  # reverse the coordinates
                start += 1 # to be compliant with the exon annotations
            else:
                end += 1
            output_stream.write(
                f"{job.contig}\t{start}\t{end}\tJUNC:{job.contig}-{start}-{end}\t{count}\t{strand}\n"
            )
    total = sum(junctions.values())
    logging.info(
        f"{job.contig}: {reads_with_junctions} reads with junctions, {total} junction "
        f"observations, {len(junctions)} distinct junctions; strand source "
        + ", ".join(f"{k}={v}" for k, v in sorted(strand_sources.items()))
    )

def collect_results(args, jobs):
    output_file = args.output.with_suffix(".bed")
    rows = 0
    with output_file.open("w") as output_stream:
        for job in jobs:
            logging.info(f"Collect counts from {job.location}")
            with job.location.open("r") as input_stream:
                for line in input_stream:
                    output_stream.write(line)
                    rows += 1
    logging.info(f"Wrote {rows} junctions to {output_file}")
    ### An empty count file is a failure, not a success. Reporting it as one hid a total loss of
    ### junctions on every BAM whose aligner writes no XS tag.
    if rows == 0:
        raise RuntimeError(
            f"No junctions were counted from {args.bam}. The output {output_file} is empty. "
            f"Check that the BAM contains spliced alignments on the requested chromosomes "
            f"({', '.join(map(str, args.chr))})."
        )
    return rows


def get_jobs(args):
    return [
        Job(
            contig=contig,
            location=args.tmp.joinpath(args.bam.stem, get_tmp_suffix())
        )
        for contig in get_all_bam_chr(args.bam, args.threads)
        if contig in args.chr
    ]


def count_junctions(args):
    sample_path = args.tmp.joinpath(args.bam.stem)
    sample_path.mkdir(parents=True, exist_ok=True)  # <<< ADD THIS LINE

    jobs = get_jobs(args)
    ### Decide once in the parent so each worker does not re-open the BAM header.
    args.transcript_oriented = is_transcript_oriented(args.bam, args.threads)
    logging.info(
        f"Span {len(jobs)} job(s) between {args.cpus} CPU(s); "
        f"transcript-oriented reads: {args.transcript_oriented}"
    )
    with multiprocessing.Pool(args.cpus) as pool:
        pool.map(partial(process_contig, args), jobs)
    collect_results(args, jobs)

    logging.debug(f"Removing temporary directory for sample {args.bam.stem} at {sample_path}")
    shutil.rmtree(sample_path)
