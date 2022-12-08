import sys
import pysam
import numpy
import pandas
import anndata
import logging
import multiprocessing
from functools import partial
from altanalyze3.utilities.logger import setup_logger
from altanalyze3.utilities.helpers import (
    lambda_chr_converter,
    get_tmp_suffix,
    TimeIt
)
from altanalyze3.utilities.constants import (
    Job,
    Annotation
)
from altanalyze3.utilities.io import (
    get_all_ref_chr,
    get_correct_contig,
    get_indexed_bed
)


# https://pythonspeed.com/articles/pandas-string-dtype-memory/
# https://h2oai.github.io/db-benchmark/
# https://stackoverflow.com/questions/20459536/convert-pandas-dataframe-to-sparse-numpy-matrix-directly
# https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.write.html


def get_jobs(args):
    # we include only chromosomes present both in tmp_coord and tmp_reference,
    # and user provided list. Contig is always prepended with 'chr'
    return [
        Job(
            contig=contig,
            location=args.tmp.joinpath(get_tmp_suffix())
        )
        for contig in get_all_ref_chr(args.tmp_coord, args.threads)
            if contig in get_all_ref_chr(args.tmp_reference, args.threads) and contig in args.chr
    ]


def get_annotation(job, query_location, reference_location, threads=None):
    threads = 1 if threads is None else threads
    collected_annotations = []
    with pysam.TabixFile(str(query_location), mode="r", parser=pysam.asBed(), threads=threads) as query_handler:
        with pysam.TabixFile(str(reference_location), mode="r", parser=pysam.asBed(), threads=threads) as reference_handler:
            reference_iter = reference_handler.fetch(get_correct_contig(job.contig, reference_handler))            # to iterate within the specific chromosome
            current_reference = next(reference_iter)                                                               # get initial value from reference iterator
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
                            name=int(current_query.name)
                        )
                        break
                    elif current_position > current_reference.end:
                        try:
                            logging.debug(f"{current_position:,} is after [{current_reference.start:,}, {current_reference.end:,})")
                            logging.debug("Attempting to switch to the next reference")
                            current_reference = next(reference_iter)
                        except StopIteration:
                            logging.debug("Run out of references")
                            current_annotation = Annotation(
                                gene=numpy.nan,
                                exon=numpy.nan,
                                strand=numpy.nan,
                                position=numpy.nan,
                                name=int(current_query.name)
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
                                    name=int(current_query.name)
                                )
                            else:
                                logging.debug("Not exact match")
                                current_annotation = Annotation(
                                    gene=current_gene,
                                    exon=current_exon,
                                    strand=current_reference.strand,
                                    position=current_position,
                                    name=int(current_query.name)
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
                                name=int(current_query.name)
                            )
                            break
                        elif current_position == current_reference.start or current_position == current_reference.end:
                            try:
                                logging.debug(f"{current_position:,} is outside of intron [{current_reference.start:,}, {current_reference.end:,})")
                                logging.debug("Attempting to switch to the next reference")
                                current_reference = next(reference_iter)
                            except StopIteration:
                                logging.debug("Run out of references")
                                current_annotation = Annotation(
                                    gene=numpy.nan,
                                    exon=numpy.nan,
                                    strand=numpy.nan,
                                    position=numpy.nan,
                                    name=int(current_query.name)
                                )
                                run_out_of_references = True
                        else:
                            sys.exit(-1)  # exit, not correct logic
                    else:
                        sys.exit(-1)  # exit, not correct logic

                logging.debug(f"""Assigning {current_annotation}""")
                collected_annotations.append(current_annotation)

    collected_annotations.sort(key = lambda i:i.name)    # need to sort the results so the order corresponds to the tmp_coord
    return collected_annotations


def process_contig(args, job):
    setup_logger(
        multiprocessing.get_logger(),
        args.loglevel
    )
    multiprocessing.current_process().name = job.contig

    logging.info(f"""Annotating aggregated junctions coordinates from {args.tmp_coord}, subset to {job.contig} chromosome""")
    logging.info(f"""Temporary results will be saved to {job.location}""")

    logging.debug("Annotating junctions start positions")
    start_annotations = get_annotation(
        job=job,
        query_location=args.tmp_start,
        reference_location=args.tmp_reference,
        threads=args.threads
    )
    logging.debug("Annotating junctions end positions")
    end_annotations = get_annotation(
        job=job,
        query_location=args.tmp_end,
        reference_location=args.tmp_reference,
        threads=args.threads
    )

    with job.location.open("wt") as output_stream:
        with pysam.TabixFile(str(args.tmp_coord), mode="r", parser=pysam.asBed(), threads=args.threads) as coord_handler:
            for (current_query, start_annotation, end_annotation) in zip(
                coord_handler.fetch(get_correct_contig(job.contig, coord_handler)), start_annotations, end_annotations
            ):
                start_position = "" if start_annotation.position == 0 else f"""_{start_annotation.position}"""
                end_position = "" if end_annotation.position == 0 else f"""_{end_annotation.position}"""
                gene = (f"""{start_annotation.gene}:""", "") if start_annotation.gene == end_annotation.gene else (f"""{start_annotation.gene}:""", f"""{end_annotation.gene}:""")
                name = f"""{gene[0]}{start_annotation.exon}{start_position}-{gene[1]}{end_annotation.exon}{end_position}"""
                strand = start_annotation.strand if start_annotation.strand == start_annotation.strand else "."
                output_stream.write(
                    f"""{job.contig}\t{current_query.start}\t{current_query.end}\t{name}\t0\t{strand}\n"""
                )


def load_query_data(args):
    query_df = None
    for jun_location, alias in zip(args.juncounts, args.aliases):
        logging.info(f"""Loading junction counts from {jun_location} as {alias}""")
        jun_df = pandas.read_csv(
            jun_location,
            usecols=[0, 1, 2, 4],
            names=["chr", "start", "end", alias],
            converters={"chr": lambda_chr_converter},
            sep="\t",
        )
        jun_df.set_index(["chr", "start", "end"], inplace=True)
        jun_df = jun_df[jun_df.index.get_level_values("chr").isin(args.chr)]                                       # subset only to those chromosomes that are provided in --chr
        query_df = jun_df if query_df is None else query_df.join(jun_df, how="outer").fillna(0)
        query_df = query_df.astype(pandas.SparseDtype("uint32", 0))                                                # saves memory and is required before exporting to h5ad

    logging.info("Sorting aggregated junctions counts by coordinates (chr, start, end) in ascending order")
    query_df.sort_index(ascending=True, inplace=True)

    tmp_counts_location = args.tmp.joinpath(get_tmp_suffix()).with_suffix(".pickle")                               # we keep index in file for safety reasons
    logging.info(f"""Saving aggregated junctions counts as a temporary pickled file to {tmp_counts_location}""")
    query_df.to_pickle(tmp_counts_location)

    query_df_coord = query_df.index.to_frame(index=False)
    tmp_coord_location = args.tmp.joinpath(get_tmp_suffix()).with_suffix(".bed")
    logging.info(f"""Saving aggregated junctions coordinates as a temporary BED file to {tmp_coord_location}""")
    query_df_coord.to_csv(
        tmp_coord_location,
        sep="\t",
        header=False,
        index=False
    )
    tmp_coord_location = get_indexed_bed(tmp_coord_location)

    query_df_start_coord = query_df_coord.copy()
    query_df_start_coord["name"] = range(0, len(query_df_start_coord.index))
    query_df_start_coord["end"] = query_df_start_coord["start"]
    query_df_start_coord.sort_values(by=["chr", "start", "end"], ascending=True, inplace=True)
    tmp_coord_start_location = args.tmp.joinpath(get_tmp_suffix()).with_suffix(".bed")
    logging.info(f"""Saving aggregated junctions start coordinates as a temporary BED file to {tmp_coord_start_location}""")
    query_df_start_coord.to_csv(
        tmp_coord_start_location,
        sep="\t",
        header=False,
        index=False
    )
    tmp_coord_start_location = get_indexed_bed(tmp_coord_start_location)

    query_df_end_coord = query_df_coord.copy()
    query_df_end_coord["name"] = range(0, len(query_df_end_coord.index))
    query_df_end_coord["start"] = query_df_end_coord["end"]
    query_df_end_coord.sort_values(by=["chr", "start", "end"], ascending=True, inplace=True)
    tmp_coord_end_location = args.tmp.joinpath(get_tmp_suffix()).with_suffix(".bed")
    logging.info(f"""Saving aggregated junctions end coordinates as a temporary BED file to {tmp_coord_end_location}""")
    query_df_end_coord.to_csv(
        tmp_coord_end_location,
        sep="\t",
        header=False,
        index=False
    )
    tmp_coord_end_location = get_indexed_bed(tmp_coord_end_location)

    return (tmp_coord_location, tmp_coord_start_location, tmp_coord_end_location, tmp_counts_location)


def load_reference_data(args, shift_start_by=None):
    logging.info(f"""Loading reference annotation from {args.ref}""")
    reference_df = pandas.read_csv(
        args.ref,
        usecols=[0, 1, 2, 3, 4, 5],
        names=["gene", "exon", "chr", "strand", "start", "end"],
        converters={"chr": lambda_chr_converter},
        skiprows=1,
        sep="\t",
    )
    if shift_start_by is not None:
        reference_df["start"] = reference_df["start"] + shift_start_by                        # to correct 1-based coordinates
    reference_df.set_index(["chr", "start", "end"], inplace=True)
    reference_df = reference_df[reference_df.index.get_level_values("chr").isin(args.chr)]    # subset only to those chromosomes that are provided in --chr
    reference_df.sort_index(ascending=True, inplace=True)                                     # this may potentially mix overlapping genes from different strands
    reference_df["name"] = reference_df["gene"] + ":" + reference_df["exon"]
    reference_df["score"] = 0                                                                 # dummy column to correspond to BED format
    reference_df.drop(["gene", "exon"], axis=1, inplace=True)                                 # droping unused columns

    tmp_reference_location = args.tmp.joinpath(get_tmp_suffix()).with_suffix(".bed")
    logging.info(f"""Saving sorted reference annotation as a temporary BED file to {tmp_reference_location}""")
    reference_df.to_csv(
        tmp_reference_location,
        sep="\t",
        columns=["name", "score", "strand"],                                                  # we have "chr", "start", "end" in the index
        header=False,
        index=True
    )

    logging.info("Compressing and indexing saved temporary reference annotation BED file")
    tmp_reference_location = get_indexed_bed(tmp_reference_location)

    return tmp_reference_location


def collect_results(args, jobs):
    collected_annotations_df = None
    for job in jobs:
        logging.info(f"""Loading annotated junction coordinates from {job.location}""")
        annotations_df = pandas.read_csv(
            job.location,
            usecols=[0, 1, 2, 3],
            names=["chr", "start", "end", "name"],
            converters={"chr": lambda_chr_converter},
            sep="\t",
        )
        annotations_df.set_index(["chr", "start", "end"], inplace=True)
        collected_annotations_df = annotations_df if collected_annotations_df is None else pandas.concat([collected_annotations_df, annotations_df])
        logging.debug(f"""Removing {job.location}""")
        job.location.unlink()

    collected_annotations_df.sort_index(ascending=True, inplace=True)
    query_df = pandas.read_pickle(args.tmp_counts)
    query_df = query_df.join(collected_annotations_df, how="left")

    var_index = query_df.index.to_frame(index=False).astype(str)
    var_index.insert(0, "name", query_df["name"].values)
    var_index = var_index["name"] + " " + var_index["chr"] + ":" + var_index["start"] + "-" + var_index["end"]
    query_df.drop(columns=["name"], inplace=True)

    logging.info("Converting aggregated junction counts to csr matrix")
    query_adata = anndata.AnnData(
        query_df.T.sparse.to_coo().tocsr(),                                        # csr matrix with intervals as columns, and samples as rows
        obs=pandas.DataFrame(index=query_df.columns),                              # how is it different from creating a columns with the name sample?
        var=pandas.DataFrame(index=var_index),                                     # the same question as above
        dtype="uint32"
    )

    h5ad_location = args.output.with_suffix(".h5ad")
    logging.info(f"""Exporting aggregated counts to {h5ad_location}""")
    query_adata.write(h5ad_location)   # no need to compress temporary file (should be faster)


def aggregate(args):
    with TimeIt():
        args.tmp_coord, args.tmp_start, args.tmp_end, args.tmp_counts = load_query_data(args)
        args.tmp_reference = load_reference_data(args, shift_start_by=-1)                     # need to shift start position to make 0-based half-open coordinates
        jobs = get_jobs(args)
        logging.info(f"""Span {len(jobs)} job(s) between {args.cpus} CPU(s)""")
        with multiprocessing.Pool(args.cpus) as pool:
            pool.map(partial(process_contig, args), jobs)
        collect_results(args, jobs)
