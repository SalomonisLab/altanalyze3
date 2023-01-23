import pysam
import pandas
import shutil
import logging
import itertools
import multiprocessing
from functools import partial
from collections import deque
from altanalyze3.utilities.logger import setup_logger
from altanalyze3.utilities.helpers import get_tmp_suffix
from altanalyze3.utilities.constants import (
    Job,
    Annotation,
    JunctionsParams,
    IntronsParams,
    AnnotationsParams,
    AnnMatchCat
)
from altanalyze3.utilities.io import (
    get_all_ref_chr,
    get_correct_contig,
    get_indexed_bed,
    export_counts_to_anndata
)


# https://pythonspeed.com/articles/pandas-string-dtype-memory/
# https://h2oai.github.io/db-benchmark/
# https://stackoverflow.com/questions/20459536/convert-pandas-dataframe-to-sparse-numpy-matrix-directly
# https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.write.html


def get_jun_annotation_jobs(args):
    # we include only chromosomes present both in jun_coords_location and references,
    # and user provided list. Contig is always prepended with 'chr'
    return [
        Job(
            contig=contig,
            location=args.tmp.joinpath(get_tmp_suffix())
        )
        for contig in get_all_ref_chr(args.jun_coords_location, args.threads)
            if contig in get_all_ref_chr(args.ref, args.threads) and contig in args.chr
    ]


class GroupedAnnotations():

    def __init__(self):
        self.__exact_match = []
        self.__partial_match = []
        self.__distant_match = []

    def add(self, sa, ea):                                                                                           # sa - start annotation, ea - end annotation

        sa_shift = "" if sa.position == 0 else f"""_{sa.position}"""                                                 # junction start coord if not exact exon match
        ea_shift = "" if ea.position == 0 else f"""_{ea.position + 1}"""                                             # junction end coorf if not exact exon match

        if sa.match == AnnMatchCat.CLOSEST or ea.match == AnnMatchCat.CLOSEST:                                       # not enough data to make a conclusion
            pass
        elif sa.gene == ea.gene and sa.strand == ea.strand:                                                          # the same gene and strand
            if sa.match == AnnMatchCat.EXON_END and ea.match == AnnMatchCat.EXON_START:                              #    exact match
                if sa.strand == "+":
                    self.__exact_match.append((f"""{sa.gene}:{sa.exon}-{ea.exon}""", sa.strand))                     #        start - end exons order
                elif sa.strand == "-":
                    self.__exact_match.append((f"""{sa.gene}:{ea.exon}-{sa.exon}""", sa.strand))                     #        end - start exons order
                else:
                    assert(False), "Not implemented logic"
            else:                                                                                                    #    partial match
                if sa.strand == "+":
                    self.__partial_match.append((f"""{sa.gene}:{sa.exon}{sa_shift}-{ea.exon}{ea_shift}""", sa.strand))  #     start - end exons order
                elif sa.strand == "-":
                    self.__partial_match.append((f"""{sa.gene}:{ea.exon}{ea_shift}-{sa.exon}{sa_shift}""", sa.strand))  #     end - start exons order
                else:
                    assert(False), "Not implemented logic"
        else:                                                                                                        # different gene and/or strand
            self.__distant_match.append((f"""{sa.gene}:{sa.exon}{sa_shift}-{ea.gene}:{ea.exon}{ea_shift}""", "."))   #    distant match

    def best(self):
        all_annotations = self.__exact_match + self.__partial_match + self.__distant_match
        if len(all_annotations) > 0:
            return all_annotations[0]
        else:
            return ("undefined", ".")


class RefDeque():

    def __init__(self, reference_iter):
        self.__reference_iter = reference_iter
        self.__exhausted = False
        self.__deque = deque()
        self.__buffer = []
        self.__last_ref = None
        self.contig = None
        self.start = None
        self.end = None
        self.name = None
        self.strand = None
        self.gene = None
        self.exon = None
        self.type = None

    def __set_state(self, reference, safe=None):
        safe = False if safe is None else safe
        self.contig = reference.contig
        self.start = reference.start
        self.end = reference.end
        self.name = reference.name
        self.strand = reference.strand
        self.gene, self.exon = self.name.split(":")
        self.type = self.exon[0].upper()
        logging.debug(f"""===> {"safely " if safe else ""}get new reference {self.contig}:{self.start:,}-{self.end:,}, {self.name}, {self.strand}""")
        if safe:
            self.__buffer.append(reference)

    def pop_left(self, safe=None):                 # if safe is True, we will be able to restore used items
        safe = False if safe is None else safe
        if safe:
            if len(self.__buffer) == 0:            # this is the first time we called pop_left with safe=True
                self.__last_ref = (
                    self.contig,
                    self.start,
                    self.end,
                    self.name,
                    self.strand,
                    self.gene,
                    self.exon,
                    self.type
                )
        else:                                      # we called unsafe pop_left, so keepeing buffer doesn't make sense anymore
            self.__buffer = []
            self.__last_ref = None
        try:
            self.__set_state(self.__deque.popleft(), safe)
        except IndexError:
            try:
                self.__deque.append(next(self.__reference_iter))
                self.__set_state(self.__deque.popleft(), safe)
            except StopIteration:
                self.__exhausted = True

    def restore(self):
        if len(self.__buffer) > 0 and self.__last_ref is not None:        # restore only when we have what to restore
            logging.debug("Restore references deque")
            self.__buffer.reverse()
            self.__deque.extendleft(self.__buffer)
            self.contig, self.start, self.end, self.name, self.strand, self.gene, self.exon, self.type = self.__last_ref
            self.__buffer = []
            self.__last_ref = None

    def is_empty(self):
        return self.__exhausted and len(self.__deque) == 0


def get_annotation(job, query_location, references_location, threads=None):
    threads = 1 if threads is None else threads
    all_annotations = []
    with pysam.TabixFile(str(query_location), mode="r", parser=pysam.asBed(), threads=threads) as query_handler:
        with pysam.TabixFile(str(references_location), mode="r", parser=pysam.asBed(), threads=threads) as references_handler:

            reference_iter = references_handler.fetch(get_correct_contig(job.contig, references_handler))
            ref_deque = RefDeque(reference_iter)                                                            # we need it because we can't roll back pysam iterator
            ref_deque.pop_left()                                                                            # getting a new reference from the start of the left

            for current_query in query_handler.fetch(get_correct_contig(job.contig, query_handler)):        # fetches all jucntions from the specific chromosome
                logging.debug(f"""---> get new region {current_query.contig}:{current_query.start:,}-{current_query.end:,}""")
                check_position = current_query.start                                                        # we always use start
                current_annotations = []

                while not ref_deque.is_empty():
                    if check_position < ref_deque.start:
                        logging.debug(f"""{check_position:,} is behind {ref_deque.contig}:{ref_deque.start:,}-{ref_deque.end:,}, {ref_deque.name}, {ref_deque.strand}""")
                        logging.debug("Closest match")
                        current_annotations.append(Annotation(ref_deque.gene, ref_deque.exon, ref_deque.strand, check_position, AnnMatchCat.CLOSEST))
                        break
                    elif check_position > ref_deque.end:
                        logging.debug(f"{check_position:,} is after {ref_deque.contig}:{ref_deque.start:,}-{ref_deque.end:,}, {ref_deque.name}, {ref_deque.strand}")
                        ref_deque.pop_left()
                    else:
                        while not ref_deque.is_empty() and check_position >= ref_deque.start:
                            if ref_deque.type == "E" and check_position == ref_deque.start:
                                logging.debug("Exact match - exon start")
                                current_annotations.append(Annotation(ref_deque.gene, ref_deque.exon, ref_deque.strand, 0, AnnMatchCat.EXON_START))
                            elif ref_deque.type == "E" and check_position == ref_deque.end:
                                logging.debug("Exact match - exon end")
                                current_annotations.append(Annotation(ref_deque.gene, ref_deque.exon, ref_deque.strand, 0, AnnMatchCat.EXON_END))
                            elif check_position > ref_deque.start and check_position < ref_deque.end:
                                if ref_deque.type == "E":
                                    logging.debug("Not exact exon match")
                                    current_annotations.append(Annotation(ref_deque.gene, ref_deque.exon, ref_deque.strand, check_position, AnnMatchCat.EXON_MID))
                                elif ref_deque.type == "I":
                                    logging.debug("Not exact intron match")
                                    current_annotations.append(Annotation(ref_deque.gene, ref_deque.exon, ref_deque.strand, check_position, AnnMatchCat.INTRON_MID))
                                else:
                                    assert(False), "Not implemented logic"
                            ref_deque.pop_left(safe=True)
                        ref_deque.restore()
                        break

                if len(current_annotations) == 0:                                            # in case we failed to annotate junction
                    logging.debug(f"""Failed to annotate {check_position:,}""")
                    current_annotations.append(None)

                all_annotations.append((current_annotations, int(current_query.score)))      # score keeps the original rows order

    all_annotations.sort(key = lambda i:i[1])    # need to sort the results so the order corresponds to the jun_coords_location
    return all_annotations


def process_jun_annotation(args, job):
    setup_logger(
        multiprocessing.get_logger(),
        args.loglevel
    )
    multiprocessing.current_process().name = job.contig

    logging.info(f"""Annotating junctions coordinates from {args.jun_coords_location}, subset to {job.contig} chromosome""")
    logging.info(f"""Temporary results will be saved to {job.location}""")

    logging.debug("Annotating junctions start coordinates")
    sorted_start_annotations = get_annotation(
        job=job,
        query_location=args.jun_starts_location,
        references_location=args.ref,
        threads=args.threads
    )

    logging.debug("Annotating junctions end coordinates")
    sorted_end_annotations = get_annotation(
        job=job,
        query_location=args.jun_ends_location,
        references_location=args.ref,
        threads=args.threads
    )

    with job.location.open("wt") as output_stream:
        with pysam.TabixFile(str(args.jun_coords_location), mode="r", parser=pysam.asBed(), threads=args.threads) as jun_coords_handler:
            for (current_coords, (start_annotations, _), (end_annotations, _)) in zip(
                jun_coords_handler.fetch(get_correct_contig(job.contig, jun_coords_handler)), sorted_start_annotations, sorted_end_annotations
            ):
                logging.debug(f"""Assigning annotation for {job.contig}:{current_coords.start:,}-{current_coords.end:,}, {current_coords.name}""")
                start_annotations = [_ for _ in start_annotations if _ is not None]
                end_annotations = [_ for _ in end_annotations if _ is not None]
                grouped_annotation = GroupedAnnotations()
                for annotations_pair in list(itertools.product(start_annotations, end_annotations)):   # if any of [start_]/[end_]annotations is [], returns []
                    grouped_annotation.add(annotations_pair[0], annotations_pair[1])
                annotation, strand = grouped_annotation.best()
                output_stream.write(
                    f"""{job.contig}\t{current_coords.start}\t{current_coords.end}\t{current_coords.name}\t{annotation}\t{strand}\n"""   # no reason to keep it BED-formatted
                )


def load_counts_data(query_locations, query_aliases, selected_chr, tmp_location, as_junctions=None, save_bed=None):
    save_bed = False if save_bed is None else save_bed
    loading_params = JunctionsParams if as_junctions is not None and as_junctions else IntronsParams

    counts_df = None
    for query_location, query_alias in zip(query_locations, query_aliases):
        logging.info(f"""Loading counts from {query_location} as {query_alias}""")
        current_df = pandas.read_csv(query_location, **loading_params)
        current_df.rename(
            columns = {"score": query_alias},
            inplace = True
        )
        current_df = current_df[current_df.index.get_level_values("chr").isin(selected_chr)]                         # subset only to those chromosomes that are provided in --chr
        counts_df = current_df if counts_df is None else counts_df.join(current_df, how="outer").fillna(0)
        counts_df = counts_df.astype(pandas.SparseDtype("int32", 0))                                                # saves memory and is required before exporting to h5ad

    logging.info("Sorting counts by coordinates in ascending order")
    counts_df.sort_index(ascending=True, inplace=True)

    counts_location = tmp_location.joinpath(get_tmp_suffix()).with_suffix(".pickle")
    logging.info(f"""Saving counts with coordinates as a temporary pickled file to {counts_location}""")
    counts_df.to_pickle(counts_location)

    coords_location, starts_location, ends_location = None, None, None
    if save_bed:                                                                                                     # no need to allow to save BED when we alreade have intron annotations
        coords_df = counts_df.index.to_frame(index=False)
        coords_df["score"] = range(0, len(coords_df.index))
        coords_location = tmp_location.joinpath(get_tmp_suffix()).with_suffix(".bed")
        logging.info(f"""Saving only coordinates as a temporary BED file to {coords_location}""")
        coords_df.to_csv(
            coords_location,
            sep="\t",
            header=False,
            index=False,
            columns=loading_params["names"]
        )
        coords_location = get_indexed_bed(coords_location, keep_original=False, force=True)

        start_coords_df = coords_df.copy()
        start_coords_df["end"] = start_coords_df["start"]
        start_coords_df.sort_values(
            by=loading_params["index_col"],
            ascending=True,
            inplace=True
        )
        starts_location = tmp_location.joinpath(get_tmp_suffix()).with_suffix(".bed")
        logging.info(f"""Saving only start coordinates as a temporary BED file to {starts_location}""")
        start_coords_df.to_csv(
            starts_location,
            sep="\t",
            header=False,
            index=False,
            columns=loading_params["names"]
        )
        starts_location = get_indexed_bed(starts_location, keep_original=False, force=True)

        end_coords_df = coords_df.copy()
        end_coords_df["start"] = end_coords_df["end"]
        end_coords_df.sort_values(
            by=loading_params["index_col"],
            ascending=True,
            inplace=True
        )
        ends_location = tmp_location.joinpath(get_tmp_suffix()).with_suffix(".bed")
        logging.info(f"""Saving only end coordinates as a temporary BED file to {ends_location}""")
        end_coords_df.to_csv(
            ends_location,
            sep="\t",
            header=False,
            index=False,
            columns=loading_params["names"]
        )
        ends_location = get_indexed_bed(ends_location, keep_original=False, force=True)

    return (counts_location, coords_location, starts_location, ends_location)


def collect_results(args):

    counts_df = None
    if args.juncounts is not None:
        collected_annotations_df = None
        for job in args.jun_annotation_jobs:
            logging.info(f"""Loading annotated junctions coordinates from {job.location}""")
            annotations_df = pandas.read_csv(job.location, **AnnotationsParams)
            collected_annotations_df = annotations_df if collected_annotations_df is None else pandas.concat([collected_annotations_df, annotations_df])
            logging.info(f"""Removing {job.location}""")
            job.location.unlink()
        logging.info(f"""Loading pickled junctions counts from {args.jun_counts_location}""")
        counts_df = pandas.read_pickle(args.jun_counts_location)
        logging.info("Joining pickled junctions counts with annotations")
        counts_df = counts_df.join(collected_annotations_df, how="left")

    int_counts_df = None
    if args.intcounts is not None:
        logging.info(f"""Loading pickled introns counts from {args.int_counts_location}""")
        int_counts_df = pandas.read_pickle(args.int_counts_location)
        int_counts_df["annotation"] = int_counts_df.index.to_frame(index=False).name.values
        int_counts_df.reset_index(level="strand", inplace=True)                                     # need to move "strand" from the index to a column

    if counts_df is None:
        counts_df = int_counts_df
    elif int_counts_df is not None:
        logging.info("Concatenating annotated junctions and introns counts from")
        counts_df = pandas.concat([counts_df, int_counts_df])

    counts_df.sort_index(ascending=True, inplace=True)

    adata_location = args.output.with_suffix(".h5ad")
    logging.info(f"""Exporting aggregated counts to {adata_location}""")
    metadata_columns = ["annotation", "strand"]
    export_counts_to_anndata(
        counts_df=counts_df,
        location=adata_location,
        counts_columns=[c for c in counts_df.columns.values if c not in metadata_columns],
        strand_coords=True,
        metadata_columns=metadata_columns
    )

    if args.bed:
        bed_location = args.output.with_suffix(".bed")
        logging.info(f"""Exporting annotated coordinates to {bed_location}""")
        counts_df.reset_index(inplace=True)
        counts_df.drop(columns=["name"], inplace=True)
        counts_df.rename(columns={"annotation": "name"}, inplace=True)
        counts_df["score"] = 0
        counts_df.to_csv(
            bed_location,
            sep="\t",
            header=False,
            index=False,
            columns=["chr", "start", "end", "name", "score", "strand"]
        )
        bed_location = get_indexed_bed(bed_location, keep_original=False, force=True)


def aggregate(args):

    if args.juncounts is not None:
        logging.info("Processing junctions counts")
        args.jun_counts_location, args.jun_coords_location, args.jun_starts_location, args.jun_ends_location = load_counts_data(
            query_locations=args.juncounts,
            query_aliases=args.aliases,
            selected_chr=args.chr,
            tmp_location=args.tmp,
            as_junctions=True,                                                       # "strand" column will be ignored
            save_bed=True
        )
        args.jun_annotation_jobs = get_jun_annotation_jobs(args)
        logging.info(f"""Span {len(args.jun_annotation_jobs)} junctions annotation job(s) between {args.cpus} CPU(s)""")
        with multiprocessing.Pool(args.cpus) as pool:
            pool.map(partial(process_jun_annotation, args), args.jun_annotation_jobs)

    if args.intcounts is not None:
        logging.info("Processing introns counts")
        args.int_counts_location, _, _, _ = load_counts_data(
            query_locations=args.intcounts,
            query_aliases=args.aliases,
            selected_chr=args.chr,
            tmp_location=args.tmp,
            save_bed=False
        )

    logging.info("Collecting results")
    collect_results(args)

    logging.info("Removing temporary directory")
    shutil.rmtree(args.tmp)