import pandas as pd
import os


class SparseMatrix:

    '''

    This class goes through all the exon-exon juction count's file for all samples
    and generates a sparse matrix. In this sparse matrix the y-axis is all the Unique junctions
    and x-axis is corresponding samples and their counts inside the matrix
    Eventually this sparse matrix will be saved in a h5ad format

    '''
    
    def __int__(self):
        self.junctionfiles,
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
        junctionfiles = ['/Users/sin9gp/altanalyze3/tests/data/intcount_results.bed']
        main_dataframe = pd.DataFrame(pd.read_csv(junctionfiles[0],sep='\t'))
        print(main_dataframe.iloc[:,3:])

        for i in range(0, len(self.junctionfiles)):
            # get unique juncitons
            data = pd.read_csv(self.junctionfiles[i])
            df = pd.DataFrame(data)
            df.drop_duplicates(subset=[3])
            #concat all the junctions in one file
        return main_dataframe

    def calculateTotalJunctionCounts(self):
        '''
        This function will iterate over all the sample files and get the junction counts for the junctions
        stored in "unique junction coordinate" file (from above function).
        '''
        junctionfiles = ['/Users/sin9gp/altanalyze3/tests/data/intcount_results.bed']
        junctionCoordinateMap = {}
        #uniqueCoordinates = self.getAllUniqueExonJunctionCoordinates()
        for i in range( 0,len(junctionfiles)):
            data = pd.read_csv(junctionfiles[i])
            df = pd.DataFrame(data)
            uniqueCoordinates = df.drop_duplicates(subset=[3])
            # main_dataframe = pd.concat([main_dataframe,df],axis =1)
            for j in range(0,len(uniqueCoordinates)):
                #check if the coordinate match, then add the counts
                if(df['junction'] == uniqueCoordinates[j]):
                    junctionCoordinateMap[uniqueCoordinates[j]] += df['counts']
        return junctionCoordinateMap

    
    def createSparseMatrix(self):
        '''
        This function will take 
        '''
        print("i am going to make sparse matrix")
            
            


        

        
sparseM = SparseMatrix()
sparseM.getAllUniqueExonJunctionCoordinates()



 
