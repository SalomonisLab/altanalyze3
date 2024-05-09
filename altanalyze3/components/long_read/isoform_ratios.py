import os
import sys
import csv
import anndata as ad
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from scipy.io import mmread
from scipy.sparse import lil_matrix
import collections
import argparse
from tqdm import tqdm # progress bar

sys.path.insert(1, os.path.join(sys.path[0], '..'))
import long_read.gff_process as gff_process

def collapse_isoforms(gff_source,ensembl_exon_dir):
    """ 
    Derive common isoform IDs across study samples for downstream comparison
    """
    
    gff_output_dir = gff_process.consolidateLongReadGFFs(gff_source, ensembl_exon_dir)


if __name__ == '__main__':
    # python iso_matrix_to_junctions.py /path/to/matrix /path/to/ensembl_exon /path/to/gff --barcode_cluster /path/to/barcode_cluster

    parser = argparse.ArgumentParser(description='Process single-cell data to generate junction matrices.')
    parser.add_argument('matrix_path', type=str, help='Path to the matrix directory.')
    parser.add_argument('ensembl_exon_dir', type=str, help='Path to the Ensembl exon directory.')
    parser.add_argument('gff_source', type=str, help='Path to the GFF source file.')
    parser.add_argument('--barcode_cluster', type=str, default=None, help='Optional: Path to the barcode cluster file.')
    parser.add_argument('--reverse_complement', type=bool, default=False, help='Optional: Reverse complement the sparse matrix cell barcodes.')
    args = parser.parse_args()

    barcode_clusters = pd.read_csv(args.barcode_cluster, sep='\t', index_col=0)
    print(barcode_clusters);sys.exit()
    exportJunctionMatrix(args.matrix_path, args.ensembl_exon_dir, args.gff_source, barcode_clusters=barcode_clusters, rev=args.reverse_complement)

    """
    gff = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/HSC_3k/ND167_HSC-3k.gtf'
    matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/HSC_3k/filtered_matrices'
    barcode_cluster = 'HSC.txt'
    ensembl_exon_dir = '/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt'
    exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff, barcode_cluster=barcode_cluster)
    """
 