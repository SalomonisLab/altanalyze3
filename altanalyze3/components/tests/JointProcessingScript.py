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
import long_read.isoform_ratio as isor

ensembl_exon_dir ='/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt'
"""
######## Process data from sample ND167 ########
# Import the cellHarmony results and save cell type pd frame for each cell barcode for each separate capture
barcode_cluster_dir = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/QueryGroups.cellHarmony.txt'
barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dir)

# Produce an junction counts h5ad for the WM35_ND20-167__HSC_3k capture
matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/HSC_3k/filtered_matrices'
gff1 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/HSC_3k/ND167_HSC-3k.gtf'
hsc_nd167_adata = junc.exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff_source, barcode_clusters=barcode_sample_dict['WM35_ND20-167__HSC_3k'])

# Produce an junction counts h5ad for the WM35_ND20-167__MPP_3k capture
matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/MPP_3k/filtered_matrices'
gff2 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/MPP_3k/ND167_MPP_3k.gtf'
mpp_nd167_adata = junc.exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff_source, barcode_clusters=barcode_sample_dict['WM35_ND20-167__MPP_3k'])

# Produce an junction counts h5ad for the WM35_ND20-167__HSC_3k capture
matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/LMPP_3k/filtered_matrices'
gff3 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/LMPP_3k/ND167_LMPP_3k.gtf'
lmpp_nd167_adata = junc.exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff_source, barcode_clusters=barcode_sample_dict['WM35_ND20-167__LMPP_3k'])

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
nd251_1_adata = junc.exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff_source, barcode_clusters=barcode_sample_dict['WF40_ND21-251_HSC_3k'],rev=True)

# Produce an junction counts h5ad for the WF40_ND21-251_HSC_3k capture (second run)
matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/PacBio/size_selected/scisoseq.seurat_info/isoforms_seurat'
gff5 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/PacBio/size_selected/ND251-2.gff'
nd251_2_adata = junc.exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff_source, barcode_clusters=barcode_sample_dict['WF40_ND21-251_HSC_3k'],rev=True)

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
am72_adata = junc.exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff_source, barcode_clusters=barcode_sample_dict['AM72-N-0206_CD34'],rev=True)
am72_adata = ad.read_h5ad('AM72.h5ad')
pseudo_pdf = iso.pseudo_cluster_counts(am72_adata,cell_threshold=5,count_threshold=0)
pseudo_pdf.to_csv('am72-combined_pseudo_cluster_counts.txt', sep='\t')

# Produce an junction counts h5ad for the WF83-N-0215_CD34 capture
matrix_path = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/WF83/scisoseq.seurat_info/isoforms_seurat'
gff7 = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/WF83/WF83.gff'
wf83_adata = junc.exportJunctionMatrix(matrix_path, ensembl_exon_dir, gff_source, barcode_clusters=barcode_sample_dict['WF83-N-0215_CD34'],rev=True)
wf83_adata = ad.read_h5ad('WF83.h5ad')

pseudo_pdf = iso.pseudo_cluster_counts(wf83_adata,cell_threshold=5,count_threshold=0)
pseudo_pdf.to_csv('wf83-combined_pseudo_cluster_counts.txt', sep='\t')

#### Run AltAnalyze2 with min_exp=4;med_exp=9, as
# python ExpressionBuilder.py --i /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/ --a unbiased --species Hs --platform RNASeq --additional /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/ExpressionInput/counts.ND-junctions.txt
# python import_scripts/AugmentEventAnnotations.py --i /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI.txt --s Hs --platform RNASeq
# python stats_scripts/metaDataAnalysis.py --i /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt --s Hs --p PSI --dPSI 0.1 --percentExp 40 --m /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/ExpressionInput/groups.ND-junctions.txt
"""

isor.collapse_isoforms(gff_source,ensembl_exon_dir)