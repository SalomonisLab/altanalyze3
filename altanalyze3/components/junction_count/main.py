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


def process_contig(args, job):
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


def collect_results(args, jobs):
    output_file = args.output.with_suffix(".bed")
    with output_file.open("w") as output_stream:
        for job in jobs:
            logging.info(f"Collect counts from {job.location}")
            with job.location.open("r") as input_stream:
                output_stream.write(input_stream.read())
    # No deletion of individual files here


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
    logging.info(f"Span {len(jobs)} job(s) between {args.cpus} CPU(s)")
    with multiprocessing.Pool(args.cpus) as pool:
        pool.map(partial(process_contig, args), jobs)
    collect_results(args, jobs)

    logging.debug(f"Removing temporary directory for sample {args.bam.stem} at {sample_path}")
    shutil.rmtree(sample_path)
