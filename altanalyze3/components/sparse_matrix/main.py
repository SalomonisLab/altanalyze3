import pandas as pd
import os


class SparseMatrix:

    '''

    This class goes through all the exon-exon juction count's file for all samples
    and generates a sparse matrix. In this sparse matrix the y-axis is all the Unique junctions
    and x-axis is corresponding samples and their counts inside the matrix
    Eventually this sparse matrix will be saved in a h5ad format

    '''

    def __int__(self, junctioncounts):
        self.
        self.junctioncounts,
        self.chr,
        self.sample,
        self.exon_region_start,
        self.exon_region_stop

    def getUniqueExonJunctionCoordinates(self):
        '''
        Iterate over a exon junction count file(outputed from Misha's juncount)
        And save all the unique coordinates as a map
        '''
        

