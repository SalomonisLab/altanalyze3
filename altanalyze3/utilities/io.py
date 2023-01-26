import pysam
import numpy
import pandas
import anndata
import logging
from altanalyze3.utilities.constants import (
    ChrConverter,
    ReferencesParams
)
from altanalyze3.utilities.helpers import get_tmp_suffix


def guard_chr(function):
    def wrapper(location, threads):
        raw_res = function(location, threads)
        if isinstance(raw_res, list):
            return list(map(ChrConverter, raw_res))
        else:
            return ChrConverter(raw_res)
    return wrapper


@guard_chr
def get_all_bam_chr(location, threads):
    with pysam.AlignmentFile(location, mode="rb", threads=threads) as bam_handler:
        return [_.contig for _ in bam_handler.get_index_statistics()]


@guard_chr
def get_all_ref_chr(location, threads):
    with pysam.TabixFile(str(location), mode="r", parser=pysam.asBed(), threads=threads) as ref_handler:
        return ref_handler.contigs


def get_indexed_bed(location, keep_original=None, force=None):
    keep_original = True if keep_original is None else keep_original
    force = False if force is None else force
    return pysam.tabix_index(
        str(location),
        preset="bed",
        keep_original=keep_original,
        force=force
    )


def get_correct_contig(contig, handler):   # Attempting to fetch both types of choromosome names
    try:
        handler.fetch(contig)
        return contig
    except ValueError:
        return contig.lstrip("chr")


def skip_bam_read(read):
    """
    Returns true of read should be skipped based on the
    specified conditions
    """
    return read.is_secondary or \
           read.is_duplicate or \
           read.is_supplementary or \
           (read.is_paired and read.mate_is_unmapped) or \
           (read.is_paired and not read.is_proper_pair)


def is_bam_paired(location, threads=None):
    """
    Returns true if alignments in the BAM file are paired-end
    """
    threads = 1 if threads is None else threads
    with pysam.AlignmentFile(location, mode="rb", threads=threads) as bam_handler:
        for read in bam_handler.fetch():                                                 # this fetches only mapped reads
            if skip_bam_read(read):
                continue
            return read.is_paired


def is_bam_indexed(location, threads=None):
    """
    Returns true if BAM file is indexed
    """
    threads = 1 if threads is None else threads
    with pysam.AlignmentFile(location, mode="rb", threads=threads) as bam_handler:
        return bam_handler.has_index()


def get_indexed_references(location, tmp_location, selected_chr=None, only_introns=None):
    only_introns = False if only_introns is None else only_introns

    logging.info(f"""Loading references from {location}""")
    references_df = pandas.read_csv(location, **ReferencesParams)
    logging.debug(f"""Loaded {len(references_df.index)} lines""")

    selected_chr = references_df.index.get_level_values("chr").unique().tolist()
    if selected_chr is not None:
        if isinstance(selected_chr, list):
            selected_chr = list(map(ChrConverter, selected_chr))
        else:
            selected_chr = ChrConverter(selected_chr)

    references_df = references_df[references_df.index.get_level_values("chr").isin(selected_chr)]         # subset to the selected chromosomes

    if only_introns:
        logging.info("Filtering references to include only introns")
        references_df = references_df[references_df["exon"].str.contains("^I")]
        logging.debug(f"""After filtering {len(references_df.index)} lines remained""")

    logging.info("Sorting references by coordinates in ascending order")
    references_df.sort_index(ascending=True, inplace=True)                                                # this may potentially mix overlapping genes from different strands
    references_df["name"] = references_df["gene"] + ":" + references_df["exon"]
    references_df["score"] = 0                                                                            # dummy column to correspond to BED format
    coords_df = references_df.index.to_frame(index=False)                                                 # we need it only for thickStart and thickEnd
    references_df["thickStart"] = coords_df.start.values
    references_df["thickEnd"] = numpy.where(
        references_df["exon"].str.upper().str.startswith("E"),
        coords_df["end"],
        coords_df["start"]
    )
    references_df.drop(["gene", "exon"], axis=1, inplace=True)                                            # droping unused columns

    target_location = tmp_location.joinpath(get_tmp_suffix()).with_suffix(".bed")
    logging.info(f"""Saving references as a BED file to {target_location}""")
    references_df.to_csv(
        target_location,
        sep="\t",
        columns=["name", "score", "strand", "thickStart", "thickEnd"],                                    # we have "chr", "start", "end" in the index
        header=False,
        index=True
    )

    return get_indexed_bed(target_location, keep_original=False, force=True)


def export_counts_to_anndata(counts_df, location, counts_columns=None, metadata_columns=None, sparse_dtype=None, fill_value=None, strand_coords=None):
    counts_columns = counts_df.columns.values if counts_columns is None else counts_columns
    metadata_columns = [] if metadata_columns is None else metadata_columns
    sparse_dtype = "int32" if sparse_dtype is None else sparse_dtype
    fill_value = 0 if fill_value is None else fill_value
    strand_coords = False if strand_coords is None else strand_coords

    def __get_name(series, strand_coords):
        if strand_coords and "strand" in series and series.at["strand"] == "-":
            return f"""{series.at["chr"]}:{series.at["end"] + 1}-{series.at["start"]}"""
        else:
            return f"""{series.at["chr"]}:{series.at["start"]}-{series.at["end"] + 1}"""

    csr_matrix = counts_df.loc[:, counts_columns].astype(pandas.SparseDtype(sparse_dtype, fill_value)).T.sparse.to_coo().tocsr()
    adata = anndata.AnnData(csr_matrix, dtype=sparse_dtype)
    adata.obs_names = counts_columns
    adata_var = counts_df.copy().loc[:, metadata_columns].astype(str)    # can be empty df if metadata_columns is []
    adata_var.index = adata_var.reset_index().agg(__get_name, axis="columns", strand_coords=strand_coords)
    adata.var = adata_var
    adata.write(location)
