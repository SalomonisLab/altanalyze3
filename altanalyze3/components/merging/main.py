import pandas as pd
import os


class ImportAnnotations:
    """
    Class to get all intron and exon junction annotations from Gene model
    """

    def __init__(self, genemodel):
        self.genemodel = genemodel
        self.chr,
        self.position,
        self.intron_start,
        self.intron_stop

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


class JunctionCoordinates:
    """
    This class is to get all the junction coordinates from exon-exon and exon-intron junction files.
    """

    def __init__(self, intronexonfile, exonexonfile):
        self.exonexon = exonexonfile
        self.intronexon = intronexonfile
        self.chromosome,
        self.strand,
        self.exon_region_start,
        self.exon_region_stop

    @property
    def strand(self):
        return self.strand

    @property
    def chromosome(self):
        return self.chromosome

    def getUniqueExonJunctionKey(self, df):
        return df.gene + df.chr + df.strand

    def findUniqueExons(self, junction_data):
        dataframe = pd.read_csv(junction_data)
        unique_junction_object = {}
        if(dataframe['exon-id'].startswith('E')):
            if (dataframe['exon-region-start'].unique() and dataframe['exon-region-end'].unique()):
                unique_junction_object[self.getUniqueExonJunctionKey(dataframe)] = [
                    dataframe['exon-region-start'], dataframe['exon-region-end']]
        return unique_junction_object

    def findUniqueExonExonJunctions(self, all_samples_dir):
        '''
        Iterate over all samples and find unique exon-exon junctions
        '''
        file_list = os.listdir(all_samples_dir)
        all_samples_count = len(file_list)
        dataframes_list = []
        unique_junction_objects = []

        for i in range(all_samples_count):
            temp_df = pd.read_csv("./csv/" + all_samples_dires[i])
            dataframes_list.append(temp_df)

        for eachSample in dataframes_list:
            unique_junction_objects.app(self.findUniqueExons(eachSample))
