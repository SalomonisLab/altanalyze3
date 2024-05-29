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

ensembl_exon_dir ='/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt'

#"""
################################### JUNCTION LEVEL ANALYSIS ###################################
######## Process data from sample ND167 ########
# Import the cellHarmony results and save cell type pd frame for each cell barcode for each separate capture
barcode_cluster_dir = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/QueryGroups.cellHarmony.txt'
barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dir)

# Produce an junction counts h5ad for the WM35_ND20-167__HSC_3k capture
matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/HSC_3k/filtered_matrices'
gff1 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/HSC_3k/ND167_HSC-3k.gtf'
hsc_nd167_adata = junc.exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff1, barcode_clusters=barcode_sample_dict['WM35_ND20-167__HSC_3k'])

# Produce an junction counts h5ad for the WM35_ND20-167__MPP_3k capture
matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/MPP_3k/filtered_matrices'
gff2 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/MPP_3k/ND167_MPP_3k.gtf'
mpp_nd167_adata = junc.exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff2, barcode_clusters=barcode_sample_dict['WM35_ND20-167__MPP_3k'])

# Produce an junction counts h5ad for the WM35_ND20-167__LMPP_3k capture
matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/LMPP_3k/filtered_matrices'
gff3 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/LMPP_3k/ND167_LMPP_3k.gtf'
lmpp_nd167_adata = junc.exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff3, barcode_clusters=barcode_sample_dict['WM35_ND20-167__LMPP_3k'])

# Reimport adata objects (automatically named after the GTF/GFF in the analysis directory)
hsc_nd167_adata = ad.read_h5ad('ND167_HSC-3k.h5ad')
mpp_nd167_adata = ad.read_h5ad('ND167_MPP_3k.h5ad')
lmpp_nd167_adata = ad.read_h5ad('ND167_LMPP_3k.h5ad')

# Ensure each has adata has capture specific cell barcodes and combine
hsc_nd167_adata = iso.append_sample_name(hsc_nd167_adata, 'WM35_ND20-167__HSC_3k')
mpp_nd167_adata = iso.append_sample_name(mpp_nd167_adata, 'WM35_ND20-167__MPP_3k')
lmpp_nd167_adata = iso.append_sample_name(lmpp_nd167_adata, 'WM35_ND20-167__LMPP_3k')

nd167_adata = ad.concat([hsc_nd167_adata, mpp_nd167_adata, lmpp_nd167_adata], axis=0, join='outer')
nd167_adata.write_h5ad('ND167.h5ad')

nd167_adata = ad.read_h5ad('ND167.h5ad')
pseudo_pdf = iso.pseudo_cluster_counts(nd167_adata,cell_threshold=5,count_threshold=0)
pseudo_pdf.to_csv('nd167-combined_pseudo_cluster_counts.txt', sep='\t')


######## Process data from sample ND251 ########

# Import the cellHarmony results and save cell type pd frame for each cell barcode for each separate capture
barcode_cluster_dir = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/Illumina/SoupX-hg38/cellHarmony/QueryGroups.cellHarmony.txt'
barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dir)

# Produce an junction counts h5ad for the WF40_ND21-251_HSC_3k capture (first run)
matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/PacBio/first_run/scisoseq.seurat_info/isoforms_seurat'
gff4 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/PacBio/first_run/ND251-1.gff'
nd251_1_adata = junc.exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff4, barcode_clusters=barcode_sample_dict['WF40_ND21-251_HSC_3k'],rev=True)

# Produce an junction counts h5ad for the WF40_ND21-251_HSC_3k capture (second run)
matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/PacBio/size_selected/scisoseq.seurat_info/isoforms_seurat'
gff5 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/PacBio/size_selected/ND251-2.gff'
nd251_2_adata = junc.exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff5, barcode_clusters=barcode_sample_dict['WF40_ND21-251_HSC_3k'],rev=True)

# Reimport adata objects (automatically named after the GTF/GFF in the analysis directory)
nd251_1_adata = ad.read_h5ad('ND251-1.h5ad')
nd251_2_adata = ad.read_h5ad('ND251-2.h5ad')

# Ensure each has adata has capture specific cell barcodes and combine
nd251_1_adata = iso.append_sample_name(nd251_1_adata, 'WF40_ND21-251_HSC_3k-1')
nd251_2_adata = iso.append_sample_name(nd251_2_adata, 'WF40_ND21-251_HSC_3k-2')

nd251_adata = ad.concat([nd251_1_adata, nd251_2_adata], axis=0, join='outer')
nd251_adata.write_h5ad('ND251-apharesis.h5ad')

pseudo_pdf = iso.pseudo_cluster_counts(nd251_adata,cell_threshold=5,count_threshold=0)
pseudo_pdf.to_csv('nd251-combined_pseudo_cluster_counts.txt', sep='\t')

######## Process data from elderly Bone Marrow samples ########
# Import the cellHarmony results and save cell type pd frame for each cell barcode for each separate capture
barcode_cluster_dir = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/Illumina/SoupX-hg38/cellHarmony/QueryGroups.cellHarmony.txt'
barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dir)

# Produce an junction counts h5ad for the AM72-N-0206_CD34
matrix_path = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/AM72/scisoseq.seurat_info/isoforms_seurat'
gff6 = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/AM72/AM72.gff'
am72_adata = junc.exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff6, barcode_clusters=barcode_sample_dict['AM72-N-0206_CD34'],rev=True)
am72_adata = ad.read_h5ad('AM72.h5ad')
pseudo_pdf = iso.pseudo_cluster_counts(am72_adata,cell_threshold=5,count_threshold=0)
pseudo_pdf.to_csv('am72-combined_pseudo_cluster_counts.txt', sep='\t')

# Produce an junction counts h5ad for the WF83-N-0215_CD34 capture
matrix_path = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/WF83/scisoseq.seurat_info/isoforms_seurat'
gff7 = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/WF83/WF83.gff'
wf83_adata = junc.exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff7, barcode_clusters=barcode_sample_dict['WF83-N-0215_CD34'],rev=True)
wf83_adata = ad.read_h5ad('WF83.h5ad')

pseudo_pdf = iso.pseudo_cluster_counts(wf83_adata,cell_threshold=5,count_threshold=0)
pseudo_pdf.to_csv('wf83-combined_pseudo_cluster_counts.txt', sep='\t')

#### Run AltAnalyze2 with min_exp=4;med_exp=9, as
# python ExpressionBuilder.py --i /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/ --a unbiased --species Hs --platform RNASeq --additional /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/ExpressionInput/counts.ND-junctions.txt
# python import_scripts/AugmentEventAnnotations.py --i /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI.txt --s Hs --platform RNASeq
# python stats_scripts/metaDataAnalysis.py --i /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt --s Hs --p PSI --dPSI 0.1 --percentExp 40 --m /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/ExpressionInput/groups.ND-junctions.txt

"""
"""
################################### ISOFORM ANNOTATION ###################################
#### Create a combined GFF and find the longest consensus isoforms
gff1 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/HSC_3k/ND167_HSC-3k.gtf'
gff2 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/MPP_3k/ND167_MPP_3k.gtf'
gff3 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/LMPP_3k/ND167_LMPP_3k.gtf'
gff4 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/PacBio/first_run/ND251-1.gff'
gff5 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/PacBio/size_selected/ND251-2.gff'
gff6 = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/AM72/AM72.gff'
gff7 = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/WF83/WF83.gff'
gencode_gff = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/gencode.v45.annotation.gff3'
gff_RefSeq = '/Users/saljh8/Dropbox/Revio/ENCODE/Aggregate/hg38.ncbiRefSeq.gtf'
gff_GTEx = '/Users/saljh8/Dropbox/Revio/GTEx/flair_filter_transcripts.gtf'
gff_hep1 = '/Users/saljh8/Dropbox/Revio/ENCODE/HEPG2/GFF/ENCFF231EVD.gtf'
gff_hep2 = '/Users/saljh8/Dropbox/Revio/ENCODE/HEPG2/GFF/ENCFF303OEI.gtf'
gff_hep3 = '/Users/saljh8/Dropbox/Revio/ENCODE/HEPG2/GFF/ENCFF333KIC.gtf'
gff_hep4 = '/Users/saljh8/Dropbox/Revio/ENCODE/HEPG2/GFF/ENCFF382WAU.gtf'
gff_hep5 = '/Users/saljh8/Dropbox/Revio/ENCODE/HEPG2/GFF/ENCFF930AKE.gtf'
gff_k1 = '/Users/saljh8/Dropbox/Revio/ENCODE/K562/GTF/ENCFF207OOQ.gtf'
gff_k2 = '/Users/saljh8/Dropbox/Revio/ENCODE/K562/GTF/ENCFF298JGJ.gtf'
gff_k3 = '/Users/saljh8/Dropbox/Revio/ENCODE/K562/GTF/ENCFF391VFS.gtf'
gff_k4 = '/Users/saljh8/Dropbox/Revio/ENCODE/K562/GTF/ENCFF584GRG.gtf'
gff_k5 = '/Users/saljh8/Dropbox/Revio/ENCODE/K562/GTF/ENCFF588BMG.gtf'
gff_k6 = '/Users/saljh8/Dropbox/Revio/ENCODE/K562/GTF/ENCFF688NNQ.gtf'
gff_k7 = '/Users/saljh8/Dropbox/Revio/ENCODE/K562/GTF/ENCFF691BVH.gtf'
gff_tfiso = '/Users/saljh8/Dropbox/Revio/ENCODE/TFiso/c_6k_unique_acc_aligns.gtf'
gff_tfiso2 = '/Users/saljh8/Dropbox/Revio/ENCODE/TFiso/c_6k_unique_acc_aligns.gff'
#gff_process.reformat_uncoventional_gtf(gff_tfiso,gff_tfiso2)

gff_source = [gff_tfiso2,gff1,gff2,gff3,gff4,gff5,gff6,gff7,gencode_gff,gff_RefSeq,gff_GTEx,gff_hep1,gff_hep2,gff_hep3,gff_hep4,gff_hep5,gff_k1,gff_k2,gff_k3,gff_k4,gff_k5,gff_k6,gff_k7]
#gff_source = [gff2,gff1,gencode_gff]
gff_output_dir = gff_process.consolidateLongReadGFFs(gff_source, ensembl_exon_dir)
#"""

#"""
#### Translate isoform into protein sequence and check for NMD
ref_gff_file = "/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/gencode.v45.annotation.gff3"
query_gff_file = "/Users/saljh8/Dropbox/Revio/Aggregate2/gff-output/combined.gff"
genome_fasta = "/Users/saljh8/Dropbox/Revio/Other/Variants/SNV/genome.fa"
transcript_associations_file = "/Users/saljh8/Dropbox/Revio/Aggregate2/gff-output/transcript_associations.txt"
cds_records, protein_records = isot.gff_translate(query_gff_file,genome_fasta,ref_gff_file,transcript_associations_file)

"""
with open("protein_sequences.fasta", "w") as protein_file:
    SeqIO.write(protein_records, protein_file, "fasta")"""

with open("cds_sequences.fasta", "w") as cds_file:
    SeqIO.write(cds_records, cds_file, "fasta")
#"""
#"""
################################### ISOFORM LEVEL ANALYSIS ###################################
#### Further collapse the isoform matrix for each sample using multi-sample consensus IDs

# Indicates the "best" isoform across smaples
#isoform_association_path = "/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/gff-output/isoform_links.txt"
isoform_association_path = "/Users/saljh8/Dropbox/Revio/Aggregate/gff-output/isoform_links.txt"

# Import the cellHarmony results and save cell type pd frame for each cell barcode for each separate capture
barcode_cluster_dir = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/QueryGroups.cellHarmony.txt'
barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dir)

# Produce an isoform counts h5ad for the WM35_ND20-167__HSC_3k capture
matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/HSC_3k/filtered_matrices'
gff_source = 'ND167_HSC-3k'
hsc_nd167_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=barcode_sample_dict['WM35_ND20-167__HSC_3k'])

# Produce an junction counts h5ad for the WM35_ND20-167__MPP_3k capture
matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/MPP_3k/filtered_matrices'
gff_source = 'ND167_MPP_3k'
mpp_nd167_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=barcode_sample_dict['WM35_ND20-167__MPP_3k'])

# Produce an junction counts h5ad for the WM35_ND20-167__LMPP_3k capture
matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/MPP_3k/filtered_matrices'
gff_source = 'ND167_LMPP_3k'
lmpp_nd167_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=barcode_sample_dict['WM35_ND20-167__LMPP_3k'])

# Combine and export pseudobulk results
#hsc_nd167_adata = ad.read_h5ad('ND167_HSC-3k-isoform.h5ad')
#mpp_nd167_adata = ad.read_h5ad('ND167_MPP-3k-isoform.h5ad')
#lmpp_nd167_adata = ad.read_h5ad('ND167_LMPP-3k-isoform.h5ad')

# Ensure each has adata has capture specific cell barcodes and combine
hsc_nd167_adata = iso.append_sample_name(hsc_nd167_adata, 'WM35_ND20-167__HSC_3k')
mpp_nd167_adata = iso.append_sample_name(mpp_nd167_adata, 'WM35_ND20-167__MPP_3k')
lmpp_nd167_adata = iso.append_sample_name(lmpp_nd167_adata, 'WM35_ND20-167__LMPP_3k')
nd167_adata = ad.concat([hsc_nd167_adata, mpp_nd167_adata, lmpp_nd167_adata], axis=0, join='outer')
nd167_adata.write_h5ad('ND167-isoform.h5ad')
pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_cluster_counts(hsc_nd167_adata,cell_threshold=5,count_threshold=0,compute_tpm=True)
pseudo_pdf.to_csv('nd167-combined_pseudo_cluster_counts.txt', sep='\t')
tpm.to_csv('nd167-combined_pseudo_cluster-tpm.txt', sep='\t')
isoform_to_gene_ratio.to_csv('nd167-combined_pseudo_cluster-ratio.txt', sep='\t')

######## Process data from sample ND251 ########
# Import the cellHarmony results and save cell type pd frame for each cell barcode for each separate capture
barcode_cluster_dir = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/Illumina/SoupX-hg38/cellHarmony/QueryGroups.cellHarmony.txt'
barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dir)

# Produce an junction counts h5ad for the WF40_ND21-251_HSC_3k capture (first run)
matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/PacBio/first_run/scisoseq.seurat_info/isoforms_seurat'
gff_source = 'ND251-1'
nd251_1_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=barcode_sample_dict['WF40_ND21-251_HSC_3k'],rev=True)

# Produce an junction counts h5ad for the WF40_ND21-251_HSC_3k capture (second run)
matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/PacBio/size_selected/scisoseq.seurat_info/isoforms_seurat'
gff_source = 'ND251-2'
nd251_2_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=barcode_sample_dict['WF40_ND21-251_HSC_3k'],rev=True)

# Reimport adata objects (automatically named after the GTF/GFF in the analysis directory)
#nd251_1_adata = ad.read_h5ad('ND251-1-isoform.h5ad')
#nd251_2_adata = ad.read_h5ad('ND251-2-isoform.h5ad')
# Ensure each has adata has capture specific cell barcodes and combine
nd251_1_adata = iso.append_sample_name(nd251_1_adata, 'WF40_ND21-251_HSC_3k-1')
nd251_2_adata = iso.append_sample_name(nd251_2_adata, 'WF40_ND21-251_HSC_3k-2')
nd251_adata = ad.concat([nd251_1_adata, nd251_2_adata], axis=0, join='outer')
nd251_adata.write_h5ad('ND251-apharesis-isoform.h5ad')
pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_cluster_counts(nd251_adata,cell_threshold=5,count_threshold=0,compute_tpm=True)
pseudo_pdf.to_csv('nd251-combined_pseudo_cluster_counts.txt', sep='\t')
tpm.to_csv('nd251-combined_pseudo_cluster-tpm.txt', sep='\t')
isoform_to_gene_ratio.to_csv('nd251-combined_pseudo_cluster-ratio.txt', sep='\t')

######## Process data from elderly Bone Marrow samples ########
# Import the cellHarmony results and save cell type pd frame for each cell barcode for each separate capture
barcode_cluster_dir = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/Illumina/SoupX-hg38/cellHarmony/QueryGroups.cellHarmony.txt'
barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dir)

# Produce an junction counts h5ad for the AM72-N-0206_CD34
matrix_path = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/AM72/scisoseq.seurat_info/isoforms_seurat'
gff_source = 'AM72'
am72_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=barcode_sample_dict['AM72-N-0206_CD34'],rev=True)
#am72_adata = ad.read_h5ad('AM72-isoform-isoform.h5ad')
pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_cluster_counts(am72_adata,cell_threshold=5,count_threshold=0,compute_tpm=True)
pseudo_pdf.to_csv('AM72-combined_pseudo_cluster_counts.txt', sep='\t')
tpm.to_csv('AM72-combined_pseudo_cluster-tpm.txt', sep='\t')
isoform_to_gene_ratio.to_csv('AM72-combined_pseudo_cluster-ratio.txt', sep='\t')

# Produce an junction counts h5ad for the WF83-N-0215_CD34 capture
matrix_path = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/WF83/scisoseq.seurat_info/isoforms_seurat'
gff_source = 'WF83'
wf83_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=barcode_sample_dict['WF83-N-0215_CD34'],rev=True)
#wf83_adata = ad.read_h5ad('WF83-isoform-isoform.h5ad')
pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_cluster_counts(wf83_adata,cell_threshold=5,count_threshold=0,compute_tpm=True)
pseudo_pdf.to_csv('WF83-combined_pseudo_cluster_counts.txt', sep='\t')
tpm.to_csv('WF83-combined_pseudo_cluster-tpm.txt', sep='\t')
isoform_to_gene_ratio.to_csv('WF83-combined_pseudo_cluster-ratio.txt', sep='\t')

#### Run AltAnalyze2
#AltAnalyze.py --accessoryAnalysis mergeFiles --input /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/Isoforms/nd167-combined_pseudo_cluster-ratio.txt --input /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/Isoforms/nd251-combined_pseudo_cluster-ratio.txt --input /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/Isoforms/AM72-combined_pseudo_cluster-ratio.txt --input /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/Isoforms/WF83-combined_pseudo_cluster-ratio.txt --replaceAsZero False --join union --output /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/Isoforms
#python AltAnalyze.py --species Hs --platform "3'array" --expdir /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/Isoforms/AltAnalyze/ExpressionInput/exp.ND-isoforms-sparse.txt  --groupdir /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/Isoforms/AltAnalyze/ExpressionInput/groups.ND-isoforms-sparse.txt --compdir /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/Isoforms/AltAnalyze/ExpressionInput/comps.ND-isoforms-sparse.txt --output /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/Isoforms/AltAnalyze --vendor Ensembl
#"""
##### Process ENCODE - K562
#"""
#isoform_association_path = "/Users/saljh8/Dropbox/Revio/Aggregate/gff-output/isoform_links.txt

# Produce an isoform counts h5ad for the K562 combined captures
matrix_path = '/Users/saljh8/Dropbox/Revio/ENCODE/K562/ENCFF078DIN.tsv'
gff_source = 'ENCFF207OOQ'
ENCFF207OOQ_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

matrix_path = '/Users/saljh8/Dropbox/Revio/ENCODE/K562/ENCFF991VSH.tsv'
gff_source = 'ENCFF298JGJ'
ENCFF298JGJ_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

matrix_path = '/Users/saljh8/Dropbox/Revio/ENCODE/K562/ENCFF516IEG.tsv'
gff_source = 'ENCFF391VFS'
ENCFF391VFS_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

matrix_path = '/Users/saljh8/Dropbox/Revio/ENCODE/K562/ENCFF668BLB.tsv'
gff_source = 'ENCFF584GRG'
ENCFF584GRG_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

matrix_path = '/Users/saljh8/Dropbox/Revio/ENCODE/K562/ENCFF448BDX.tsv'
gff_source = 'ENCFF588BMG'
ENCFF588BMG_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

matrix_path = '/Users/saljh8/Dropbox/Revio/ENCODE/K562/ENCFF946RLA.tsv'
gff_source = 'ENCFF688NNQ'
ENCFF688NNQ_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

matrix_path = '/Users/saljh8/Dropbox/Revio/ENCODE/K562/ENCFF302ZXZ.tsv'
gff_source = 'ENCFF691BVH'
ENCFF691BVH_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

ENCFF207OOQ_adata = iso.append_sample_name(ENCFF207OOQ_adata, 'ENCFF207OOQ')
ENCFF298JGJ_adata = iso.append_sample_name(ENCFF298JGJ_adata, 'ENCFF298JGJ')
ENCFF391VFS_adata = iso.append_sample_name(ENCFF391VFS_adata, 'ENCFF391VFS')
ENCFF584GRG_adata = iso.append_sample_name(ENCFF584GRG_adata, 'ENCFF584GRG')
ENCFF588BMG_adata = iso.append_sample_name(ENCFF588BMG_adata, 'ENCFF588BMG')
ENCFF688NNQ_adata = iso.append_sample_name(ENCFF688NNQ_adata, 'ENCFF688NNQ')
ENCFF691BVH_adata = iso.append_sample_name(ENCFF691BVH_adata, 'ENCFF691BVH')

K562_adata = ad.concat([ENCFF207OOQ_adata, ENCFF298JGJ_adata, ENCFF391VFS_adata, ENCFF584GRG_adata, ENCFF588BMG_adata, ENCFF688NNQ_adata, ENCFF691BVH_adata], axis=0, join='outer')
K562_adata.write_h5ad('K562-isoform.h5ad')

K562_adata = ad.read_h5ad('K562-isoform.h5ad')
pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_counts_no_cluster(K562_adata,count_threshold=0,compute_tpm=True)
pseudo_pdf.to_csv('K562-isoform_counts.txt', sep='\t')
tpm.to_csv('K562-isoform-tpm.txt', sep='\t')
isoform_to_gene_ratio.to_csv('K562-isoform-ratio.txt', sep='\t')
#"""
#"""
##### Process ENCODE - HEPG2

#isoform_association_path = "/Users/saljh8/Dropbox/Revio/Aggregate/gff-output/isoform_links.txt

# Produce an isoform counts h5ad for the K562 combined captures
matrix_path = '/Users/saljh8/Dropbox/Revio/ENCODE/HEPG2/ENCFF382KCL.tsv'
gff_source = 'ENCFF231EVD'
ENCFF231EVD_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

matrix_path = '/Users/saljh8/Dropbox/Revio/ENCODE/HEPG2/ENCFF382KVB.tsv'
gff_source = 'ENCFF333KIC'
ENCFF333KIC_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

matrix_path = '/Users/saljh8/Dropbox/Revio/ENCODE/HEPG2/ENCFF859LEL.tsv'
gff_source = 'ENCFF382WAU'
ENCFF382WAU_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

matrix_path = '/Users/saljh8/Dropbox/Revio/ENCODE/HEPG2/ENCFF151KTU.tsv'
gff_source = 'ENCFF930AKE'
ENCFF930AKE_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

matrix_path = '/Users/saljh8/Dropbox/Revio/ENCODE/HEPG2/ENCFF833JVB.tsv'
gff_source = 'ENCFF303OEI'
ENCFF303OEI_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

ENCFF231EVD_adata = iso.append_sample_name(ENCFF231EVD_adata, 'ENCFF231EVD')
ENCFF333KIC_adata = iso.append_sample_name(ENCFF333KIC_adata, 'ENCFF333KIC')
ENCFF382WAU_adata = iso.append_sample_name(ENCFF382WAU_adata, 'ENCFF382WAU')
ENCFF930AKE_adata = iso.append_sample_name(ENCFF930AKE_adata, 'ENCFF930AKE')
ENCFF303OEI_adata = iso.append_sample_name(ENCFF303OEI_adata, 'ENCFF303OEI')

HEPG2_adata = ad.concat([ENCFF231EVD_adata, ENCFF333KIC_adata, ENCFF382WAU_adata, ENCFF930AKE_adata, ENCFF303OEI_adata], axis=0, join='outer')
HEPG2_adata.write_h5ad('HEPG2-isoform.h5ad')

HEPG2_adata = ad.read_h5ad('HEPG2-isoform.h5ad')
pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_counts_no_cluster(HEPG2_adata,count_threshold=0,compute_tpm=True)
pseudo_pdf.to_csv('HEPG2-isoform_counts.txt', sep='\t')
tpm.to_csv('HEPG2-isoform-tpm.txt', sep='\t')
isoform_to_gene_ratio.to_csv('HEPG2-isoform-ratio.txt', sep='\t')

##### Process GTEx
matrix_path = '/Users/saljh8/Dropbox/Revio/GTEx/quantification_flair_filter.counts.clean.txt'
gff_source = 'flair_filter_transcripts'
GTEx_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=None,mtx=False)

pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_counts_no_cluster(GTEx_adata,count_threshold=0,compute_tpm=True)
pseudo_pdf.to_csv('GTEx-isoform_counts.txt', sep='\t')
tpm.to_csv('GTEx-isoform-tpm.txt', sep='\t')
isoform_to_gene_ratio.to_csv('GTEx-isoform-ratio.txt', sep='\t')

#"""
"""
gff_marrow = '/Users/saljh8/Dropbox/Revio/ENCODE/Aggregate/combined-PacBio.gff'
gff_gencode = '/Users/saljh8/Dropbox/Revio/ENCODE/Aggregate/Gencode.gff'
gff_k562 = '/Users/saljh8/Dropbox/Revio/ENCODE/Aggregate/K562.gff'
gff_HEPG2 = '/Users/saljh8/Dropbox/Revio/ENCODE/Aggregate/HEPG2.gff'
gff_RefSeq = '/Users/saljh8/Dropbox/Revio/ENCODE/Aggregate/hg38.ncbiRefSeq.gtf'
gff_GTEx = '/Users/saljh8/Dropbox/Revio/GTEx/flair_filter_transcripts.gtf'
gff_source = [gff_marrow,gff_gencode,gff_RefSeq]
gff_output_dir = gff_process.consolidateLongReadGFFs(gff_source, ensembl_exon_dir)
"""
