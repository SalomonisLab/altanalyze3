from fileinput import filename
import timeit
import logging
import pysam
import multiprocessing as mp
import subprocess
import shlex
import pandas as pd

filename = "/Users/sin9gp/altanalyze3/tests/data/junction_dir/Cal27P5-2__junction.bed"
def bgzip(filename):
    """Call bgzip to compress a file."""
    print("i am bgzippping")
    subprocess.Popen(['bgzip', '-f', filename])

def tabix_index_junctions(filename,
        preset="gff", chrom=0, start=1, end=2, skip=0, comment="#"):
    """Call tabix to create an index for a bgzip-compressed file."""
    print("I am indexing")

    df = pd.read_csv(filename, sep = '\t', header = None, names=["chr", "start", "stop", "annotation", "splice_count"])
    df.sort_values(by=["chr", "start","stop"], inplace=True)
                                                                                                                                                                                  
    # pysam.tabix_index(df, force = False, preset = "bed")


def tabix_query(filename, chrom, start, end):
    """Call tabix and generate an array of strings for each line it returns."""
    query = '{}:{}-{}'.format(chrom, start, end)
    process = subprocess.Popen(['tabix', '-f', filename, query], stdout=subprocess.PIPE)
    for line in process.stdout:
        yield line.strip().split()

def sort_index_junctions(filename):
    sort_command = ['sort','-k1V','-k2n','-k3n',filename]
    command = shlex.join(sort_command)
    p1 = subprocess.Popen(sort_command)
    #p2 = subprocess.Popen(['bgzip', '-f'],stdin=p1.stdout)
    #p1.wait()
    print("p1 return: ", p1.returncode)
    #print("p2 return: ", p2.returncode)

#sort_index_junctions(filename)

#bgzip(filename=filename)
# tabix_index(zippedfilename,preset="gff", chrom=0, start=1, end=2, skip=0, comment="#")

# out_sorted = 'myfile.sorted'
# out_zipped= out_sorted + ".gz"

# with open(out_zipped,'w') as sort_zip_out :
#     #-k1,1 -k2,2n -k3,3n
#     #cmd="sort -k1V -k2n -k3n /Users/sin9gp/altanalyze3/tests/data/junction_dir/Cal27P5-2__junction.bed"
#     remote_command = ['sort','-k1V','-k2n','-k3n',filename]
#     command = shlex.join(remote_command)
#     p1 = subprocess.Popen(command,stdout=subprocess.PIPE)
#     p2 = subprocess.Popen(['bgzip', '-f'],stdin=p1.stdout,stdout = sort_zip_out)
#     #working
#     # p1 = Popen(['sort','-k1V','-k2n','-k3n', '-f',filename])
    
#     #p2 = Popen(['bgzip','-f'], stdin=p1.stdout, stdout= sort_zip_out)
#     p1.stdout.close()  #finish first subprocess before starting second
#     print("i am waiting")
#     p1.wait()  #wait for results to be written

#when these two subprocesses are finished, 
#tabix_index(out_zipped) 

tabix_index_junctions("/Users/sin9gp/altanalyze3/tests/data/junction_dir/subset/Cal27P5-1__junction.bed")