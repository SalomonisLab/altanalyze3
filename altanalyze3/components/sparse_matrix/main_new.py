from xml.etree.ElementTree import tostring
import pandas as pd
import os
import numpy as np
from scipy.sparse import csr_matrix, dok_matrix
import time
from timeit import default_timer as timer
from multiprocessing import Pool, cpu_count

def square(n):

    time.sleep(2)

    return n * n


def main():

    start = timer()
    print(f'starting computations on {cpu_count()} cores')
    file_list = os.listdir("/Users/sin9gp/altanalyze3/tests/data/bed")
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
