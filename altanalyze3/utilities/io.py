import pysam
import pandas
import anndata
import logging
from altanalyze3.utilities.helpers import lambda_chr_converter


def guard_chr(function):
    def prefix(c):
        return c if c.startswith("chr") else f"chr{c}"
    def wrapper(location, threads):
        raw_res = function(location, threads)
        if isinstance(raw_res, list):
            return [prefix(c) for c in raw_res]
        else:
            return prefix(raw_res)
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


def is_bam_paired(location, threads):
    """
    Returns true of alignments in the BAM file are paired-end
    """
    with pysam.AlignmentFile(location, mode="rb", threads=threads) as bam_handler:
        for read in bam_handler.fetch():                                                 # this fetches only mapped reads
            if skip_bam_read(read):
                continue
            return read.is_paired


def get_indexed_reference(location, selected_chr_list=None, shift_start_by=None, only_introns=None):
    only_introns = False if only_introns is None else only_introns
    logging.info(f"""Loading references from {location}""")
    references_df = pandas.read_csv(
        location,
        usecols=[0, 1, 2, 3, 4, 5],
        names=["gene", "chr", "strand", "exon", "start", "end"],
        converters={"chr": lambda_chr_converter},
        skiprows=1,
        sep="\t",
    )
    logging.debug(f"""Loaded {len(references_df.index)} lines""")
    
    selected_chr_list = references_df.chr.unique().tolist()
    if selected_chr_list is not None:
        selected_chr_list = list(map(lambda_chr_converter, selected_chr_list))

    if only_introns:
        logging.info("Filtering references to include only introns")
        references_df = references_df[references_df["exon"].str.contains("^I")]
        logging.debug(f"""After filtering {len(references_df.index)} lines remained""")

    if shift_start_by is not None:
        logging.debug(f"""Shifting start coordinates by {shift_start_by}""")
        references_df["start"] = references_df["start"] + shift_start_by                                  # to correct 1-based coordinates

    references_df.set_index(["chr", "start", "end"], inplace=True)
    references_df = references_df[references_df.index.get_level_values("chr").isin(selected_chr_list)]    # subset to the selected chromosomes

    logging.info("Sorting references by coordinates in ascending order")
    references_df.sort_index(ascending=True, inplace=True)                                                # this may potentially mix overlapping genes from different strands
    references_df["name"] = references_df["gene"] + ":" + references_df["exon"]
    references_df["score"] = 0                                                                            # dummy column to correspond to BED format
    references_df.drop(["gene", "exon"], axis=1, inplace=True)                                            # droping unused columns

    target_location = location.with_suffix(".bed")
    logging.info(f"""Saving references as a BED file to {target_location}""")
    references_df.to_csv(
        target_location,
        sep="\t",
        columns=["name", "score", "strand"],                                                              # we have "chr", "start", "end" in the index
        header=False,
        index=True
    )

    return get_indexed_bed(target_location, keep_original=False, force=True)


def export_counts_to_anndata(counts_df, location, counts_columns=None, metadata_columns=None, sparse_dtype=None, fill_value=None):
    counts_columns = counts_df.columns.values if counts_columns is None else counts_columns
    sparse_dtype = "uint32" if sparse_dtype is None else sparse_dtype
    fill_value = 0 if fill_value is None else fill_value

    csr_matrix = counts_df.loc[:, counts_columns].astype(pandas.SparseDtype(sparse_dtype, fill_value)).T.sparse.to_coo().tocsr()
    adata = anndata.AnnData(csr_matrix, dtype=sparse_dtype)
    adata.obs_names = counts_columns
    adata.var_names = counts_df.index.to_frame(index=False).astype(str)[counts_df.index.names].agg("-".join, axis=1)
    if metadata_columns is not None:
        adata_var = counts_df.copy().loc[:, metadata_columns]
        adata_var.index = adata.var.index.copy()
        adata.var = adata_var
    adata.write(location)
