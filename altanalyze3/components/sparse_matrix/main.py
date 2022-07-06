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
        self.junctionfiles
        self.junctioncounts,
        self.chr,
        self.sample,
        self.exon_region_start,
        self.exon_region_stop

    def getAllUniqueExonJunctionCoordinates(self):
        '''
        Iterate over a exon junction count file(outputed from Misha's juncount)
        And concat all the unique coordinates in one dataframe
        '''
        junctionfiles = ['intcount_results.bed']
        main_dataframe = pd.DataFrame(pd.read_csv(junctionfiles[0]))

        for i in range(1, len(junctionfiles)):
            # get unique juncitons
            data = pd.read_csv(junctionfiles[i])
            df = pd.DataFrame(data)
            df.drop_duplicates(subset=["junction"])
            #concat all the junctions in one file
            main_dataframe = pd.concat([main_dataframe,df],axis =1)
        
        return main_dataframe

 
