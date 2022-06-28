import pandas as pd


class ImportAnnotations:
    """
    Class to get all intron and exon junction annotations from Gene model
    """


def __init__(self, genemodel):
    self.genemodel = genemodel
    self.chr,
    self.position


def splicesite_annotations(self):
    """
    Import splice-site to exon annotations from the gene model
    """
    data = pd.read_csv(self.genemodel, header=None, names=[
                       "gene_id", "chr", "exon_region_id", "exon_start", "exon_stop", "exon_annotations"])
    self.chr = data.chr
    self.position = data.exon_start
    return {self.chr, self.position}


def intron_annotations(self):
    """
    import all start, end position and chromosome for each intron with its name (e.g. , ENSG00000123:I2.1)
    """
    data = pd.read_csv(self.genemodel, header=None, names=[
        "gene_id", "chr", "exon_region_id", "exon_start", "exon_stop", "exon_annotations"])

    if (data.exon_region_id.startswith("I")):
        intron_start = data.exon_start
        intron_stop = data.exon_stop
        intron_chr = data.chr
    return {intron_start, intron_stop, intron_chr}
