from xml.etree.ElementTree import tostring
import pandas as pd
import os
import numpy as np
from scipy.sparse import csr_matrix, dok_matrix
import time
from timeit import default_timer as timer
from multiprocessing import Pool, cpu_count

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


class JunctionCoordinateCalculation:

    '''
    This class is a structure to divide a bed (junction) file into chunks (by chromosome)
    and multi-process them (multi-CPU cores). 

    chunk-size = based on chromosome
    '''
    
    def __init__(self):
        self.chunksize = 50
    
    def read_chunk(self):
        
        with open("/Users/sin9gp/altanalyze3/tests/data/intcount_results.bed") as f:
            count = 0
            lines = []
            for line in f:
                count +=1
                lines.append(line)
                if count % self.chunksize == 0:
                    self.createSparseMatrix(lines) 


    def createSparseMatrix(self,dataarray, sampleFileArray):
        '''
        This function will take each sample file (junction coordinates) and divide them by chromosome 
        csr_matrix(row,col)
        '''
        rowCount = self.uniqueJunctionCoordinatecount(dataarray)
        colCount = self.calculateTotalSampleFiles(sampleFileArray)

        # creating sparse matrix
        sparseMatrix = csr_matrix((dataarray, (rowCount, colCount)), 
                          shape = (3, 3)).toarray()
        print(sparseMatrix)

    def totalJunctionCoordinatecount(self, dataarray):
        '''
        this function will take the chunks of files and count them
        This will be used to determine the number of rows of the sparse matrix.
        '''
        return len(dataarray)

    def totalSampleFiles(self,sampleFileArray):
        '''
        this funnction calculates total number of sample files and returns it.
        This will determine the number of columns of sparse matrix
        '''
        return len(sampleFileArray)
    
    def getJunctionAndCounts(self,sampleFile):
        '''
        args: sample file
        output: all unique junctions 
        '''
        bamdata = pd.read_csv(sampleFile, sep='\t')
        junctions = pd.DataFrame(bamdata).iloc[:,[3]].values.tolist()
        counts = pd.DataFrame(bamdata).iloc[:,[4]].values.tolist()
        return [{'junction1':1},{'junction2':2},{'junction3':1},{'junction4':0},{'junction5':0},{'junction6':0},
        {'junction7':0},{'junction8':0},{'junction9':0},{'junction10':0},{'junction11':0}]

    def reader(self,filelist):
        lines = []
        chunksize = 50
        
        for file in filelist:
            print(file)
            with open('/Users/sin9gp/altanalyze3/tests/data/bed/' + file) as f:
                count = 0
                for line in f:
                    count +=1
                    lines.append(line)
                    if count % chunksize == 0:
                        print("read first 50 lines")
                    
            

        
    
def readBamFiles(self):
    '''to d0 - dont forget there re two sbed files per sample - one for intron and other for exon-exon junciotn'''
    bamfiles = ['/Users/sin9gp/altanalyze3/tests/data/sample1.bed','/Users/sin9gp/altanalyze3/tests/data/sample2.bed']
    sampleids = []
    for file in bamfiles:
        # print(sampleids)
        bamdata = pd.read_csv(file, sep='\t')
        df = pd.DataFrame(bamdata)
        sampleids = [os.path.basename(file)]
        junctionCoordinates = self.getJunctionAndCounts(file)
        # calculate matrix dimensions
        M = self.totalJunctionCoordinatecount(df.iloc[:,[3]]) #number of rows
        N = self.totalSampleFiles(bamfiles) #number of columns
    dok_sparse = dok_matrix((M,N))
    data = []
    for col in sampleids:
        for junction in junctionCoordinates:
            for val in junction.values():
                data.append(val)
                dok_sparse[junction,col] = val
    
    csr = dok_sparse.tocsr()
    print(csr)

  

# def main():
#     '''
#     Use Python filesystem to load the sample files
#     Use Multi-processing to read and divide them by chromosome
#     call sparse matrix function on each process
#     '''
#     pool = Pool(4) # number of cores you want to use
#     file_list = os.listdir("/Users/sin9gp/altanalyze3/tests/data/bed")
#     df_list = pool.map(reader, file_list) #creates a list of the loaded df's
#     df = pd.concat(df_list) # concatenates all the df's into a single df

#     for file in file_list:
#         if file.endswith(".bed"):
#             print("hello", file)

   
def square(n):

    time.sleep(2)

    return n * n


def main():

    start = timer()

    print(f'starting computations on {cpu_count()} cores')

    values = (2, 4, 6, 8)

    with Pool() as pool:
        res = pool.map(square, values)
        print(res)

    end = timer()
    print(f'elapsed time: {end - start}')
            
            

        
    
if __name__ == '__main__':
    '''
    Pool object which offers a convenient means of parallelizing the execution of a function across 
    multiple input values, distributing the input data across processes

    Map() - This method chops the iterable into a number of chunks which it submits to the process pool as separate tasks.
    The (approximate) size of these chunks can be specified by setting chunksize to a positive integer.
    '''
    
    
    # junctionCoordinates = JunctionCoordinateCalculation()
    # file_list = os.listdir("/Users/sin9gp/altanalyze3/tests/data/bed")
    # print(file_list)
    # junctionCoordinates.reader(file_list)
    main()

        
   
# benchmark - time the matrix conversion 
# convert to csr_matrix?
# test this with 6-8 samples (run on the cluster)
            
        





 
