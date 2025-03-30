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

def exportRatios(ob,out):
    pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_cluster_counts(ob,cell_threshold=5,count_threshold=0,compute_tpm=True)
    pseudo_pdf.to_csv(out+'_pseudo_cluster_counts.txt', sep='\t')
    tpm.to_csv(out+'_pseudo_cluster-tpm.txt', sep='\t')
    isoform_to_gene_ratio.to_csv(out+'_pseudo_cluster-ratio.txt', sep='\t')

ensembl_exon_dir = '/Users/saljh8/Desktop/Code/AltAnalyze/AltDatabase/EnsMart112/ensembl/Mm/Mm_Ensembl_exon.txt'

build_healthy_iso_ref = False
if build_healthy_iso_ref:
    #### Create a combined GFF and find the longest consensus isoforms
    ref_gff_file = "gencode.vM36.annotation.gff3"
    gff1 = '/Users/saljh8/Dropbox/Revio/TC/flnc-1.collapsed.sorted.filtered_lite.gff'
    gff2 = '/Users/saljh8/Dropbox/Revio/TC/flnc-2.collapsed.sorted.filtered_lite.gff'

    gff_source = [ref_gff_file,gff1,gff2]
    gff_output_dir = gff_process.consolidateLongReadGFFs(gff_source, ensembl_exon_dir)

    #### Translate isoform into protein sequence and check for NMD
    query_gff_file = gff_output_dir+"/combined.gff"
    genome_fasta = "/Users/saljh8/Dropbox/Revio/Other/Variants/SNV/genome.fa"
    transcript_associations_file = gff_output_dir+"/transcript_associations.txt"
    cds_records, transcript_records, protein_records = isot.gff_translate(query_gff_file,genome_fasta,ref_gff_file,transcript_associations_file)

    with open("protein_sequences.fasta", "w") as protein_file:
        SeqIO.write(protein_records, protein_file, "fasta")

    with open("cds_sequences.fasta", "w") as cds_file:
        SeqIO.write(cds_records, cds_file, "fasta")

process_bulk_samples = True
if process_bulk_samples:

    isoform_association_path = "/Users/saljh8/Dropbox/Revio/TC/gff-output/isoform_links.txt"

    # Produce an isoform counts h5ad for the K562 combined captures
    matrix_path = '/Users/saljh8/Dropbox/Revio/TC/flnc-1.collapsed.abundance.txt'
    gff_source = 'flnc-1.collapsed.sorted.filtered_lite'
    flnc1_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

    matrix_path = '/Users/saljh8/Dropbox/Revio/TC/flnc-2.collapsed.abundance.txt'
    gff_source = 'flnc-2.collapsed.sorted.filtered_lite'
    flnc2_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

    flnc1_adata = ad.read_h5ad('flnc-1.collapsed.sorted.filtered_lite-isoform.h5ad')
    pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_counts_no_cluster(flnc1_adata,count_threshold=0,compute_tpm=True)
    pseudo_pdf.to_csv('flnc1-isoform_counts.txt', sep='\t')
    tpm.to_csv('flnc1-isoform-tpm.txt', sep='\t')
    isoform_to_gene_ratio.to_csv('flnc1-isoform-ratio.txt', sep='\t')

    flnc2_adata = ad.read_h5ad('flnc-2.collapsed.sorted.filtered_lite-isoform.h5ad')
    pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_counts_no_cluster(flnc2_adata,count_threshold=0,compute_tpm=True)
    pseudo_pdf.to_csv('flnc2-isoform_counts.txt', sep='\t')
    tpm.to_csv('flnc2-isoform-tpm.txt', sep='\t')
    isoform_to_gene_ratio.to_csv('flnc2-isoform-ratio.txt', sep='\t')
