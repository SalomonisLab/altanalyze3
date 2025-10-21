import os,sys
import pandas as pd
import glob
from tqdm import tqdm

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from ..long_read import gff_process as gff_process
from ..long_read import isoform_translation as isot

def main(gff_source, genome_fasta, ensembl_exon_dir, ref_gff_file):
    """ This module assigns isoform annotations to the statistical output from comparisons.py """

    #### Create a combined GFF and find the longest consensus isoforms
    output = gff_process.consolidateLongReadGFFs(gff_source, ensembl_exon_dir)
    if '.txt' in output:
        transcript_associations_file = output
        gff_output_dir = os.path.dirname(transcript_associations_file)
    else:
        gff_output_dir = output
        transcript_associations_file = gff_output_dir + "/transcript_associations.txt"
    
    #### Translate isoform into protein sequence and check for NMD
    print (transcript_associations_file)
    if os.path.isfile(os.path.join(gff_output_dir, "combined.gff")):
        query_gff_file = gff_output_dir+"/combined.gff"
    else:
        query_gff_file = gff_source[0] # Use the input gff - not a collapsed version
    
    print ('Using ref_gff:',ref_gff_file)
    cds_records, transcript_records, protein_records = isot.gff_translate(query_gff_file,genome_fasta,ref_gff_file,transcript_associations_file)

    from Bio import SeqIO
    with open("protein_sequences.fasta", "w") as protein_file:
        SeqIO.write(protein_records, protein_file, "fasta")

    with open("transcript_sequences.fasta", "w") as cds_file:
        SeqIO.write(transcript_records, cds_file, "fasta")

    with open("cds_sequences.fasta", "w") as cds_file:
        SeqIO.write(cds_records, cds_file, "fasta")
    
    return query_gff_file, transcript_associations_file
