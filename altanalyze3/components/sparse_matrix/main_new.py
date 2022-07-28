from xml.etree.ElementTree import tostring
import pandas as pd
import os
import numpy as np
from scipy.sparse import dok_matrix, dok_array, coo_matrix
import time
from timeit import default_timer as timer
from multiprocessing import Pool, cpu_count
from csv import DictReader
from csv import reader
import argparse



def grouper():
    '''
    Reads one sample file line by line and yeilds the data along with the sample name 
    '''
    return 1

def createSparseMatrix(count, samplefile):
    '''
    Input: takes one file at a time and puts it into sparse matrix

    dok_matrix - subclass of Python dict
        keys are (row, column) index tuples (no duplicate entries allowed)
        values are corresponding non-zero values
    
    in our case row - junction id, column - sample id
    '''
    bedfolder = '/Users/sin9gp/altanalyze3/tests/data/bed'
    #challenges to solve - figure out how to get the dimensions of matrix
    bamdata = pd.read_csv(file,sep='\t')
    df = pd.DataFrame(bamdata)
    M = len(df.iloc[:,[3]])
    N = len(os.listdir(bedfolder))
    print("M,N", M, N)
    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build Gene Model from GTF')
    parser.add_argument('--samplebam',type=str,default=None,help='the path to the bed/txt folder')
    args = parser.parse_args()
    p = Pool(4)
    '''
    # Use 'sampleFileReader' to split test data into
	# groups you can process without using a
	# ton of RAM. You'll probably want to 
	# increase the chunk size considerably
	# to something like 1000 lines per core.
    '''
    sampleFiles = os.listdir("/Users/sin9gp/altanalyze3/tests/data/bed/")
    totaljunctions = 0
    samplefileticker = 0
    for file in sampleFiles:
        samplefileticker += 1
        with open('/Users/sin9gp/altanalyze3/tests/data/bed/' + file) as eachsamplefile:
            r = pd.read_csv(eachsamplefile, sep='\t')
            df = pd.DataFrame(r)
            junctioncoordinates = df.iloc[:,[3]]
            junctioncoordinatesCount = len(junctioncoordinates)
            totaljunctions+=junctioncoordinatesCount
            #print(file, junctioncoordinatesCount, totaljunctions)
    M = totaljunctions
    N = samplefileticker
    #print(M,N)
    dok_sparse = dok_array((M,N))
    concatenateddfs = []
    row, col, data = [],[],[]
    giant_dict = {}
    for file in sampleFiles:
        with open('/Users/sin9gp/altanalyze3/tests/data/bed/' + file) as eachsamplefile:
            for line in eachsamplefile:
                print(line.split('\t')[3])
                giant_dict[line.split('\t')[3]] = {'sampleid':file,'spliceevents':line.split('\t')[4].rstrip()}
                # temp = map(int,line.split('\t'))
                
    print(giant_dict)
            # df = pd.DataFrame(r)
            # junctioncoordinates = df.iloc[:,[3]]
            # spliceEvents = df.iloc[:,[4]]
            # concatenateddfs.append(df)
            # for i in len(junctioncoordinates):
            #     dok_sparse




    




            #results = p.map(createSparseMatrix, 1, eachsamplefile)







    

