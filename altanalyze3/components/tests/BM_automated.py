import os
import sys
import csv
import anndata as ad
import pandas as pd
import collections
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import long_read.isoform_matrix as iso
import long_read.isoform_automate as isoa
import long_read.comparisons as comp

#import multiprocessing
# Create multiprocessing context using 'spawn'
#mp_context = multiprocessing.get_context('spawn')
#num_cores = mp_context.cpu_count()

if __name__ == '__main__':
    # Load local database files
    ensembl_exon_dir = '/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt'
    gene_symbol_file = '/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl-annotations.txt'
    gencode_gff = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/gencode.v45.annotation.gff3'
    genome_fasta = "/Users/saljh8/Dropbox/Revio/Other/Variants/SNV/genome.fa"

    # Load cell barcode to cluster associations (order in terms of similar clusters - if possible)
    bc_dir1 = '/Users/saljh8/Dropbox/Revio/Young-Healthy-Revio/Illumina/SoupX-hg38/cellHarmony-labels-young-old.txt'
    bc_dir2 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/QueryGroups.cellHarmony.txt'
    bc_dir3 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/Illumina/SoupX-hg38/cellHarmony/QueryGroups.cellHarmony.txt'
    bc_dir4 = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/Illumina/SoupX-hg38/cellHarmony/QueryGroups.cellHarmony.txt'
    bc_dir5 = '/Users/saljh8/Dropbox/Revio/MDS-Revio/Aggregate/QueryGroups.cellHarmony.txt'
    barcode_cluster_dirs = [bc_dir1,bc_dir2,bc_dir3,bc_dir4,bc_dir5]
    # Union of all imported clusters (in identified order)
    cluster_order = iso.return_cluster_order(barcode_cluster_dirs)

    # Load sample metadata
    metadata_file = '/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/tests/BM_metadata.txt'
    sample_dict = isoa.import_metadata(metadata_file)

    # Pre-process all samples in the metadata file
    isoa.pre_process_samples(metadata_file, barcode_cluster_dirs, ensembl_exon_dir)

    # Integrate all sample processed results into combined junction/splicing, isoform, and isoform ratio files
    isoa.combine_processed_samples(metadata_file, barcode_cluster_dirs, ensembl_exon_dir, gencode_gff, genome_fasta)

    # Perform pairwise comparisons between groups and cell-types
    condition1 = 'young'
    condition2 = 'aged'
    condition3 = 'MDS-pre'
    conditions = [(condition1,condition2),(condition2,condition3)]
    #cluster_order = ['HSC-1','MEP'] ### Only consider specific cell type rather than all
    analyses = ['junction','isoform','isoform-ratio']
    comp.compute_differentials(sample_dict,conditions,cluster_order,gene_symbol_file,analyses=analyses)
