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
###############################################################################################################################################################

class Segment (enum.Enum):

    PRIME_5 = 1
    PRIME_3 = 2
    INTRON = 3
    UKNOWN = 4

class Counter:

    def __init__(self, bam, ref, threads=None):
        self.bam = bam
        self.ref = ref
        self.threads = 1 if threads is None else threads
        self.paired = self.__is_paired()
        self.reset()

    def reset(self):
        self.overlaps = {}
        self.cached_reads = {}

    def get_segment(self, read, intron, span=None):
        span = 0 if span is None else span
        if read.reference_end - intron.start >= span and intron.start - read.reference_start >= span:
            return Segment.PRIME_5
        elif read.reference_start - intron.start >= 0 and intron.end - read.reference_end >= 0:
            return Segment.INTRON
        elif read.reference_end - intron.end >= span and intron.end - read.reference_start >= span:
            return Segment.PRIME_3
        else:
            return Segment.UKNOWN

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

    def update_overlaps(self, current_data, cached_data):
        current_contig, current_intron_start, current_intron_end, current_intron_name, current_segment = current_data
        if cached_data is not None:
            cached_contig, cached_intron_start, cached_intron_end, cached_intron_name, cached_segment = cached_data
            assert(current_contig == cached_contig)
            if current_segment != cached_segment and current_segment in [Segment.PRIME_5, Segment.INTRON] and cached_segment in [Segment.PRIME_5, Segment.INTRON]:
                contig = current_contig if cached_segment == Segment.INTRON else cached_contig
                intron_start = current_intron_start if cached_segment == Segment.INTRON else cached_intron_start
                intron_end = current_intron_end if cached_segment == Segment.INTRON else cached_intron_end
                intron_name = current_intron_name if cached_segment == Segment.INTRON else cached_intron_name
                default_or_current_value = self.overlaps.get(
                    (contig, intron_start, intron_end),
                    {
                        "name": intron_name,
                        "p5_count": 0,
                        "p3_count": 0
                    }
                )
                self.overlaps[(contig, intron_start, intron_end)] = {
                    "name": intron_name,
                    "p5_count": default_or_current_value["p5_count"] + 1,
                    "p3_count": default_or_current_value["p3_count"]
                }
            elif current_segment != cached_segment and current_segment in [Segment.PRIME_3, Segment.INTRON] and cached_segment in [Segment.PRIME_3, Segment.INTRON]:
                contig = current_contig if cached_segment == Segment.INTRON else cached_contig
                intron_start = current_intron_start if cached_segment == Segment.INTRON else cached_intron_start
                intron_end = current_intron_end if cached_segment == Segment.INTRON else cached_intron_end
                intron_name = current_intron_name if cached_segment == Segment.INTRON else cached_intron_name
                default_or_current_value = self.overlaps.get(
                    (contig, intron_start, intron_end),
                    {
                        "name": intron_name,
                        "p5_count": 0,
                        "p3_count": 0
                    }
                )
                self.overlaps[(contig, intron_start, intron_end)] = {
                    "name": intron_name,
                    "p5_count": default_or_current_value["p5_count"],
                    "p3_count": default_or_current_value["p3_count"] + 1
                }
            else:
                pass
        else:
            default_or_current_value = self.overlaps.get(
                (current_contig, current_intron_start, current_intron_end),
                {
                    "name": current_intron_name,
                    "p5_count": 0,
                    "p3_count": 0
                }
            )
            if current_segment == Segment.PRIME_5:
                self.overlaps[(current_contig, current_intron_start, current_intron_end)] = {
                    "name": current_intron_name,
                    "p5_count": default_or_current_value["p5_count"] + 1,
                    "p3_count": default_or_current_value["p3_count"]
                }
            elif current_segment == Segment.PRIME_3:
                self.overlaps[(current_contig, current_intron_start, current_intron_end)] = {
                    "name": current_intron_name,
                    "p5_count": default_or_current_value["p5_count"],
                    "p3_count": default_or_current_value["p3_count"] + 1
                }

    def calculate(self, contig, span):
        with pysam.AlignmentFile(self.bam, mode="rb", threads=self.threads) as bam_handler:
            with pysam.TabixFile(self.ref, mode="r", parser=pysam.asBed(), threads=self.threads) as ref_handler:
                intron_iter = ref_handler.fetch(contig)
                intron = next(intron_iter)                                                                           # get initial value from intron iterator
                for read in bam_handler.fetch(contig):
                    if read.reference_start - intron.end >= 0:
                        try:
                            intron = next(intron_iter)
                        except StopIteration:
                            break
                    else:
                        current_data = (
                            contig,
                            intron.start,
                            intron.end,
                            intron.name,
                            self.get_segment(read, intron, span)
                        )
                        if self.paired and read.query_name not in self.cached_reads:
                            self.cached_reads[read.query_name] = current_data
                        else:
                            self.update_overlaps(
                                current_data,
                                self.cached_reads.get(read.query_name, None)
                            )
                            try:
                                del self.cached_reads[read.query_name]
                            except KeyError:
                                pass

    def calculate_se_old(self, contig, span):
        with pysam.AlignmentFile(self.bam, mode="rb", threads=self.threads) as bam_handler:
            with pysam.TabixFile(self.ref, mode="r", parser=pysam.asBed(), threads=self.threads) as ref_handler:
                for intron in ref_handler.fetch(contig):
                    fiv_p, thr_p = [
                        bam_handler.count(
                            contig=contig,
                            start=pos,
                            end=pos+1,
                            read_callback=lambda r: r.reference_end - pos >= span or pos - r.reference_start >= span
                        )
                        for pos in [intron.start, intron.end]
                    ]
                    if fiv_p > 0 or thr_p > 0:
                        self.overlaps[(contig, intron.start, intron.end)] = {
                            "name": intron.name,
                            "p5_count": fiv_p,
                            "p3_count": thr_p
                        }

    def export(self, out):
        print(f"Export temporary results to {out}")
        with open(out, "w") as out_handler:
            for (contig, start, end), v in self.overlaps.items():
                out_handler.write(f"""{contig}\t{start}\t{end}\t{v["name"]}\t0\t{v["p5_count"]}\t{v["p3_count"]}\n""")


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
