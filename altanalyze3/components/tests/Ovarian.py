import os
import sys
import csv
import anndata as ad
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from scipy.io import mmread
import collections
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import long_read.isoform_junctions as junc
import long_read.isoform_matrix as iso
import long_read.isoform_ratios as isor
import long_read.isoform_translation as isot
import long_read.gff_process as gff_process
from Bio import SeqIO

query_gff_file = '/Users/saljh8/Dropbox/LongRead/PacBio/Ovarian/correct_merge.collapsed.gff'
ensembl_exon_dir = '/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt'
#gff_output_dir = gff_process.consolidateLongReadGFFs(query_gff_file, ensembl_exon_dir)

ref_gff_file = "/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/gencode.v45.annotation.gff3"
genome_fasta = "/Users/saljh8/Dropbox/Revio/Other/Variants/SNV/genome.fa"
transcript_associations_file = "/Users/saljh8/Dropbox/LongRead/PacBio/Ovarian/gff-output/transcript_associations.txt"
cds_records, transcript_records, protein_records = isot.gff_translate(query_gff_file,genome_fasta,ref_gff_file,transcript_associations_file)

with open("transcript_sequences.fasta", "w") as cds_file:
    SeqIO.write(transcript_records, cds_file, "fasta")
    
with open("protein_sequences.fasta", "w") as protein_file:
    SeqIO.write(protein_records, protein_file, "fasta")

with open("cds_sequences.fasta", "w") as cds_file:
    SeqIO.write(cds_records, cds_file, "fasta")
