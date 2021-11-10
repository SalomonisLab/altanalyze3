# Development guide

## Input

The users can either have a single bam file or a bunch of bam file (i.e. 1000 bam files). The format of the bam file should follow certain rules {`FILL THE BLANK`}.

## Running the program 

The current program conducted the following steps:

1. For each alignment record in bam file, if its CIGAR string contains "N" (potentially supportive reads for splicing junctions), we store its start and end position in the reference genome (GRCh38). After this step, we should have a dictionary-like object {`FILL THE BLANK`} storing the start and end position of all junction reads.

2. Merge reads that are overlapping in a certain {`FILL THE BLANK`} region, write it to a bed file as a feature.

The above description should verbose enough, and programming language-agnostic, a nice [example](https://timoast.github.io/sinto/basic_usage.html#how-the-fragment-file-is-generated) I liked a lot.

Potential ideas for making it faster:

1. Algorithm level, iterating each alignment record may not the most efficient way, think about leveraging data structure like suffix array or index the genome fasta. (I am not an expert at all.)

2. Parallelization, So far `multiprocessing` is the module that nevel fails me, refer to a [notebook](https://github.com/frankligy/ScPyT/blob/main/tricks/8_parallelization.ipynb) I wrote a while ago for the live example. Otherthan that, I know the options including `pathos` which used `dill` instead of `pickle`, and `dask` for scheduling (already beyond simple parallelization, it's more for workflow).

## Output

Each bam file will have a corresponding bed file (or two bed files) {`FILL THE BLANK`}, the schema of the bed file can be found in this [altanalyze link](https://altanalyze.readthedocs.io/en/latest/RNASeq_sample_data/).