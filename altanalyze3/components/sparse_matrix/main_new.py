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


data_path = '/Users/sin9gp/altanalyze3/tests/data/junction_dir/'

def grouper():
    '''
    Reads one sample file line by line and yeilds the data along with the sample name 
    '''
    return 1

def createSparseMatrix(giant_dictionary):
    '''
    Input: takes giant dictionary where key is the junction and 
    value is a map of sample name and number of spliceevents


    dok_matrix - subclass of Python dict
        keys are (row, column) index tuples (no duplicate entries allowed)
        values are corresponding non-zero values
    
    in our case row - junction id, column - sample id and value would be the number of splice events
    '''
    M = len(giant_dictionary.keys()) #number of rows
    N = len(os.listdir(data_path))
    S = dok_matrix((M,N), dtype=np.float32)
    junctions = giant_dict.keys()
    print(junctions)
        
            
           
    print(str(M), str(N) + ' dimensions of the sparse matrix')
    
    
    
    

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
    
    sampleFiles = os.listdir(data_path)
    totaljunctions = 0
    samplefileticker = 0
    for file in sampleFiles:
        samplefileticker += 1
        with open(data_path + file) as eachsamplefile:
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
        with open(data_path + file) as eachsamplefile:
            for line in eachsamplefile:
                junctionCoordinateKey = line.split('\t')[3]
                if(junctionCoordinateKey in giant_dict.keys()):
                    #if the junction already exists in the dictionary           
                    updatedspliceeventcount = giant_dict[junctionCoordinateKey]['spliceevents'] + line.split('\t')[4].rstrip()
                    giant_dict[junctionCoordinateKey] = {'sampleid':file,'spliceevents':line.split('\t')[4].rstrip()}
                giant_dict[junctionCoordinateKey] = {'sampleid':file,'spliceevents':line.split('\t')[4].rstrip()}
    createSparseMatrix(giant_dict)

           








    

