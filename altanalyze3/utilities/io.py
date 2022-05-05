import pysam


def guard_chr(function):
    def prefix(c):
        return c if c.startswith("chr") else f"chr{c}"
    def wrapper(contig):
        raw_res = function(contig)
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
    with pysam.TabixFile(location, mode="r", parser=pysam.asBed(), threads=threads) as ref_handler:
        return ref_handler.contigs


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