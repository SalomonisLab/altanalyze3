from xml.etree.ElementTree import tostring
import pandas as pd
import os
import numpy as np
from scipy.sparse import csr_matrix, dok_matrix
import time
from timeit import default_timer as timer
from multiprocessing import Pool, cpu_count
from csv import DictReader
from csv import reader

def createMatrix(row):
    #print('/Users/sin9gp/altanalyze3/tests/data/bed/' + lines)
    M = len(row) #number of rows
    N = 5 #number of columns
    dok_sparse = dok_matrix((M,N))
    data = []
    file_list = os.listdir("/Users/sin9gp/altanalyze3/tests/data/bed")
    for i in range(0,len(file_list)):
        for junction in row[3]:
            data.append(junction)
            dok_sparse[junction,i] = row[3]
    return dok_sparse

    

def main():
    '''
    Pool object which offers a means of parallelizing the execution of a function across 
    multiple input values, distributing the input data across processes

    Map() - This method chops the iterable into a number of chunks which it submits to the process pool as separate tasks.
    The (approximate) size of these chunks can be specified by setting chunksize to a positive integer.
    '''
    start = timer()
    print(f'starting computations on {cpu_count()} cores')
    file_list = os.listdir("/Users/sin9gp/altanalyze3/tests/data/bed")
    
    p = Pool(processes=cpu_count() - 7)
    print(p)
    for file in file_list:
        with open('/Users/sin9gp/altanalyze3/tests/data/bed/' + file) as f:
            csv_dict_reader = reader(f, delimiter='\t')
            #print(row[0])
            with Pool() as pool:
                res = pool.map(createMatrix, csv_dict_reader)
                print(res)

    end = timer()
    print(f'elapsed time: {end - start}')


#fix this divide process - divide by sample (each sample has two files (intron and junction files))
            

        
    
if __name__ == '__main__':
    main()
