import sys
import pysam
import numpy
import pandas
import shutil
import logging
import multiprocessing
from functools import partial
from altanalyze3.utilities.logger import setup_logger
from altanalyze3.utilities.helpers import get_tmp_suffix
from altanalyze3.utilities.constants import (
    Job,
    Annotation,
    JunctionsParams,
    IntronsParams,
    AnnotationsParams
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


def get_annotation(job, query_location, references_location, threads=None):
    threads = 1 if threads is None else threads
    collected_annotations = []
    with pysam.TabixFile(str(query_location), mode="r", parser=pysam.asBed(), threads=threads) as query_handler:
        with pysam.TabixFile(str(references_location), mode="r", parser=pysam.asBed(), threads=threads) as references_handler:
            references_iter = references_handler.fetch(get_correct_contig(job.contig, references_handler))         # to iterate within the specific chromosome
            current_reference = next(references_iter)                                                              # get initial value from reference iterator
            run_out_of_references = False                                                                          # to stop iteration over references when they run out
            for current_query in query_handler.fetch(get_correct_contig(job.contig, query_handler)):               # fetches all jucntions from the specific chromosome
                logging.debug(f"""Current query [{current_query.start:,}, {current_query.end:,})""")

                while not run_out_of_references:
                    logging.debug(f"""Current reference [{current_reference.start:,}, {current_reference.end:,}),  {current_reference.name}, {current_reference.strand}""")
                    current_gene, current_exon = current_reference.name.split(":")
                    current_type = current_exon[0].upper()
                    current_position = current_query.start             # we always use start, because query_location will be either start or end sorted with identical columns

                    if current_position < current_reference.start:
                        logging.debug(f"""{current_position:,} is behind [{current_reference.start:,}, {current_reference.end:,})""")
                        current_annotation = Annotation(
                            gene=current_gene,
                            exon=current_exon,
                            strand=current_reference.strand,
                            position=current_position,
                            order=int(current_query.score)
                        )
                        break
                    elif current_position > current_reference.end:
                        try:
                            logging.debug(f"{current_position:,} is after [{current_reference.start:,}, {current_reference.end:,})")
                            logging.debug("Attempting to switch to the next reference")
                            current_reference = next(references_iter)
                        except StopIteration:
                            logging.debug("Run out of references")
                            current_annotation = Annotation(
                                gene=numpy.nan,
                                exon=numpy.nan,
                                strand=numpy.nan,
                                position=numpy.nan,
                                order=int(current_query.score)
                            )
                            run_out_of_references = True
                    elif current_type == "E":
                        if current_position >= current_reference.start and current_position <= current_reference.end:
                            logging.debug(f"""{current_position:,} is within or adjacent to exon [{current_reference.start:,}, {current_reference.end:,})""")
                            if current_position == current_reference.start or current_position == current_reference.end:
                                logging.debug("Exact match")
                                current_annotation = Annotation(
                                    gene=current_gene,
                                    exon=current_exon,
                                    strand=current_reference.strand,
                                    position=0,
                                    order=int(current_query.score)
                                )
                            else:
                                logging.debug("Not exact match")
                                current_annotation = Annotation(
                                    gene=current_gene,
                                    exon=current_exon,
                                    strand=current_reference.strand,
                                    position=current_position,
                                    order=int(current_query.score)
                                )
                            break
                    elif current_type == "I":
                        if current_position > current_reference.start and current_position < current_reference.end:
                            logging.debug(f"""{current_position:,} is within intron [{current_reference.start:,}, {current_reference.end:,})""")
                            current_annotation = Annotation(
                                gene=current_gene,
                                exon=current_exon,
                                strand=current_reference.strand,
                                position=current_position,
                                order=int(current_query.score)
                            )
                            break
                        elif current_position == current_reference.start or current_position == current_reference.end:
                            try:
                                logging.debug(f"{current_position:,} is outside of intron [{current_reference.start:,}, {current_reference.end:,})")
                                logging.debug("Attempting to switch to the next reference")
                                current_reference = next(references_iter)
                            except StopIteration:
                                logging.debug("Run out of references")
                                current_annotation = Annotation(
                                    gene=numpy.nan,
                                    exon=numpy.nan,
                                    strand=numpy.nan,
                                    position=numpy.nan,
                                    order=int(current_query.score)
                                )
                                run_out_of_references = True
                        else:
                            sys.exit(-1)  # exit, not correct logic
                    else:
                        sys.exit(-1)      # exit, not correct logic

                logging.debug(f"""Assigning {current_annotation}""")
                collected_annotations.append(current_annotation)

    collected_annotations.sort(key = lambda i:i.order)    # need to sort the results so the order corresponds to the jun_coords_location
    return collected_annotations


def process_jun_annotation(args, job):
    setup_logger(
        multiprocessing.get_logger(),
        args.loglevel
    )
    multiprocessing.current_process().name = job.contig

    logging.info(f"""Annotating junctions coordinates from {args.jun_coords_location}, subset to {job.contig} chromosome""")
    logging.info(f"""Temporary results will be saved to {job.location}""")

    logging.debug("Annotating junctions start coordinates")
    start_annotations = get_annotation(
        job=job,
        query_location=args.jun_starts_location,
        references_location=args.ref,
        threads=args.threads
    )

    logging.debug("Annotating junctions end coordinates")
    end_annotations = get_annotation(
        job=job,
        query_location=args.jun_ends_location,
        references_location=args.ref,
        threads=args.threads
    )

    with job.location.open("wt") as output_stream:
        with pysam.TabixFile(str(args.jun_coords_location), mode="r", parser=pysam.asBed(), threads=args.threads) as jun_coords_handler:
            for (current_coords, start_annotation, end_annotation) in zip(
                jun_coords_handler.fetch(get_correct_contig(job.contig, jun_coords_handler)), start_annotations, end_annotations
            ):
                start_shift = "" if start_annotation.position == 0 else f"""_{start_annotation.position}"""
                end_shift = "" if end_annotation.position == 0 else f"""_{end_annotation.position}"""
                gene = (f"""{start_annotation.gene}:""", "") if start_annotation.gene == end_annotation.gene else (f"""{start_annotation.gene}:""", f"""{end_annotation.gene}:""")
                annotation = f"""{gene[0]}{start_annotation.exon}{start_shift}-{gene[1]}{end_annotation.exon}{end_shift}"""
                strand = start_annotation.strand if start_annotation.strand == start_annotation.strand else "."
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
        counts_df = counts_df.astype(pandas.SparseDtype("uint32", 0))                                                # saves memory and is required before exporting to h5ad

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
            logging.debug(f"""Removing {job.location}""")
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

    logging.debug("Collecting results")
    collect_results(args)

    logging.debug("Removing temporary directory")
    shutil.rmtree(args.tmp)