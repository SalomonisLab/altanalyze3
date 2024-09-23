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
from oncosplice.metadataAnalysis import *

def exportRatios(ob,out):
    pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_cluster_counts(ob,cell_threshold=5,count_threshold=0,compute_tpm=True)
    pseudo_pdf.to_csv(out+'_pseudo_cluster_counts.txt', sep='\t')
    tpm.to_csv(out+'_pseudo_cluster-tpm.txt', sep='\t')
    isoform_to_gene_ratio.to_csv(out+'_pseudo_cluster-ratio.txt', sep='\t')

# Reference transcriptome
ensembl_exon_dir = '/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt'

############# GTF Files ###############
# Reference annotations (note - the program saves results automatically witht the GTF name)
gencode_gff = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/gencode.v45.annotation.gff3'
gff_RefSeq = '/Users/saljh8/Dropbox/Revio/ENCODE/Aggregate/hg38.ncbiRefSeq.gtf'
gff_GTEx = '/Users/saljh8/Dropbox/Revio/GTEx/flair_filter_transcripts.gtf'
# Young
gff_wm34 = '/Users/saljh8/Dropbox/Revio/Young-Healthy-Revio/PacBio/WM34/WM34.gff'
gff_bm27 = '/Users/saljh8/Dropbox/Revio/Young-Healthy-Revio/PacBio/BM27/BM27.gff'
gff_bf21 = '/Users/saljh8/Dropbox/Revio/Young-Healthy-Revio/PacBio/BF21/BF21.gff'
gff_167h = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/HSC_3k/ND167_HSC-3k.gtf'
gff_167m = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/MPP_3k/ND167_MPP_3k.gtf'
gff_167l = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/LMPP_3k/ND167_LMPP_3k.gtf'
gff_251_1 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/PacBio/first_run/ND251-1.gff'
gff_251_2 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/PacBio/size_selected/ND251-2.gff'
# Aging
gff_am72 = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/AM72/AM72.gff'
gff_wf83 = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/WF83/WF83.gff'
gff_bf71 = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/BF71/BF71.gff'
# MDS
gff_5801M_pre ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/MDS-pre/5801-pre.gff'
gff_5801M_post ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/MDS-post/5801-post.gff'
gff_5801A_pre1 = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML/run1/5801-AML-1.gff'
gff_5801A_pre2 ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML/run2/5801-AML-2.gff'
gff_5801A_therapy ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML-therapy/5801IR.gff'
gff_5801A_post ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML-post/5801-9_341.gff'
gff_am75M_pre ='/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/pre/Y1258.gff'
gff_am75M_therapy ='/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/on/Y1607.gff'
gff_am75M_post ='/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/post/Y5116.gff'
gff_wf78M_pre ='/Users/saljh8/Dropbox/Revio/MDS-Revio/WF78-MDS/pre/Y1493.gff'
gff_wf78M_post ='/Users/saljh8/Dropbox/Revio/MDS-Revio/WF78-MDS/post/Y1737.gff'

############# Matrix Directories ###############
# Young
matrix_wm34 = '/Users/saljh8/Dropbox/Revio/Young-Healthy-Revio/PacBio/WM34/isoforms_seurat'
matrix_bm27 = '/Users/saljh8/Dropbox/Revio/Young-Healthy-Revio/PacBio/BM27/isoforms_seurat'
matrix_bf21 = '/Users/saljh8/Dropbox/Revio/Young-Healthy-Revio/PacBio/BF21/isoforms_seurat'
matrix_167h = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/HSC_3k/filtered_matrices'
matrix_167m = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/MPP_3k/filtered_matrices'
matrix_167l = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/LMPP_3k/filtered_matrices'
matrix_251_1 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/PacBio/first_run/scisoseq.seurat_info/isoforms_seurat'
matrix_251_2 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/PacBio/size_selected/scisoseq.seurat_info/isoforms_seurat'
# Aging
matrix_am72 = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/AM72/scisoseq.seurat_info/isoforms_seurat'
matrix_wf83 = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/WF83/scisoseq.seurat_info/isoforms_seurat'
matrix_bf71 = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/BF71/isoforms_seurat'
# MDS
matrix_5801M_pre ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/MDS-pre/isoforms_seurat'
matrix_5801M_post ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/MDS-post/isoforms_seurat'
matrix_5801A_pre1 = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML/run1/isoforms_seurat'
matrix_5801A_pre2 ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML/run2/isoforms_seurat'
matrix_5801A_therapy ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML-therapy/isoforms_seurat'
matrix_5801A_post ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML-post/isoforms_seurat'
matrix_am75M_pre ='/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/pre/isoforms_seurat'
matrix_am75M_therapy ='/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/on/isoforms_seurat'
matrix_am75M_post ='/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/post/isoforms_seurat'
matrix_wf78M_pre ='/Users/saljh8/Dropbox/Revio/MDS-Revio/WF78-MDS/pre/isoforms_seurat'
matrix_wf78M_post ='/Users/saljh8/Dropbox/Revio/MDS-Revio/WF78-MDS/post/isoforms_seurat'

bc_dir1 = '/Users/saljh8/Dropbox/Revio/Young-Healthy-Revio/Illumina/SoupX-hg38/cellHarmony-labels-young-old.txt'
bc_dir2 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/QueryGroups.cellHarmony.txt'
bc_dir3 = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND251_MAS-Seq/Illumina/SoupX-hg38/cellHarmony/QueryGroups.cellHarmony.txt'
bc_dir4 = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/Illumina/SoupX-hg38/cellHarmony/QueryGroups.cellHarmony.txt'
bc_dir5 = '/Users/saljh8/Dropbox/Revio/MDS-Revio/Aggregate/QueryGroups.cellHarmony.txt'
barcode_cluster_dirs = [bc_dir1,bc_dir2,bc_dir3,bc_dir4,bc_dir5]
barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dirs)

healthy_junction_counts = False
if healthy_junction_counts:

    ######## Process data from Young Bone Marrow samples ########
    ### Requires 256GB of RAM
    adata = junc.exportJunctionMatrix(matrix_wm34, ensembl_exon_dir, gff_wm34, barcode_clusters=barcode_sample_dict['WM34_CD34'],rev=True)
    adata = junc.exportJunctionMatrix(matrix_bf21, ensembl_exon_dir, gff_bf21, barcode_clusters=barcode_sample_dict['BF21-CD34'],rev=True)
    adata = junc.exportJunctionMatrix(matrix_bm27, ensembl_exon_dir, gff_bm27, barcode_clusters=barcode_sample_dict['BM27_CD34'],rev=True)
    
    hsc_nd167_adata = junc.exportJunctionMatrix(matrix_167h, ensembl_exon_dir, gff_167h, barcode_clusters=barcode_sample_dict['WM35_ND20-167__HSC_3k'])
    mpp_nd167_adata = junc.exportJunctionMatrix(matrix_167m, ensembl_exon_dir, gff_167m, barcode_clusters=barcode_sample_dict['WM35_ND20-167__MPP_3k'])
    lmpp_nd167_adata = junc.exportJunctionMatrix(matrix_167l, ensembl_exon_dir, gff_167l, barcode_clusters=barcode_sample_dict['WM35_ND20-167__LMPP_3k'])
    
    # Ensure each has adata has capture specific cell barcodes and combine
    hsc_nd167_adata = iso.append_sample_name(hsc_nd167_adata, 'WM35_ND20-167__HSC_3k')
    mpp_nd167_adata = iso.append_sample_name(mpp_nd167_adata, 'WM35_ND20-167__MPP_3k')
    lmpp_nd167_adata = iso.append_sample_name(lmpp_nd167_adata, 'WM35_ND20-167__LMPP_3k')
    nd167_adata = ad.concat([hsc_nd167_adata, mpp_nd167_adata, lmpp_nd167_adata], axis=0, join='outer')
    nd167_adata.write_h5ad('ND167.h5ad')
    
    nd251_1_adata = junc.exportJunctionMatrix(matrix_251_1, ensembl_exon_dir, gff_251_1, barcode_clusters=barcode_sample_dict['WF40_ND21-251_HSC_3k'],rev=True)
    nd251_2_adata = junc.exportJunctionMatrix(matrix_251_2, ensembl_exon_dir, gff_251_2, barcode_clusters=barcode_sample_dict['WF40_ND21-251_HSC_3k'],rev=True)
    nd251_1_adata = iso.append_sample_name(nd251_1_adata, 'WF40_ND21-251_HSC_3k-1')
    nd251_2_adata = iso.append_sample_name(nd251_2_adata, 'WF40_ND21-251_HSC_3k-2')
    nd251_adata = ad.concat([nd251_1_adata, nd251_2_adata], axis=0, join='outer')
    nd251_adata.write_h5ad('ND251-apharesis.h5ad')
    
    ######## Process data from elderly Bone Marrow samples ########

    # Produce an junction counts h5ad for the AM72-N-0206_CD34
    junc.exportJunctionMatrix(matrix_am72, ensembl_exon_dir, gff_am72, barcode_clusters=barcode_sample_dict['AM72-N-0206_CD34'],rev=True)
    junc.exportJunctionMatrix(matrix_wf83, ensembl_exon_dir, gff_wf83, barcode_clusters=barcode_sample_dict['WF83-N-0215_CD34'],rev=True)
    junc.exportJunctionMatrix(matrix_bf71, ensembl_exon_dir, gff_bf71, barcode_clusters=barcode_sample_dict['BF71-N-0421_CD34'],rev=True)
    
    ######## Process data from MDS Bone Marrow samples ########

    adata = junc.exportJunctionMatrix(matrix_5801M_pre,ensembl_exon_dir,gff_5801M_pre,barcode_clusters=barcode_sample_dict['WM71-0121-5801_CD34'],rev=True)
    adata = junc.exportJunctionMatrix(matrix_5801M_post,ensembl_exon_dir,gff_5801M_post,barcode_clusters=barcode_sample_dict['MDS5801-postHMA_CD34'],rev=True)

    p5801_aml_adata_1 = junc.exportJunctionMatrix(matrix_5801A_pre1,ensembl_exon_dir,gff_5801A_pre1,barcode_clusters=barcode_sample_dict['MDS5801-preVEN_342'],rev=True)
    p5801_aml_adata_2 = junc.exportJunctionMatrix(matrix_5801A_pre2,ensembl_exon_dir,gff_5801A_pre2,barcode_clusters=barcode_sample_dict['MDS5801-preVEN_342'],rev=True,)
    
    nd251_1_adata = ad.read_h5ad('ND251-1.h5ad')
    nd251_2_adata = ad.read_h5ad('ND251-2.h5ad')

    # Ensure each has adata has capture specific cell barcodes and combine
    nd251_1_adata = iso.append_sample_name(nd251_1_adata, 'WF40_ND21-251_HSC_3k-1')
    nd251_2_adata = iso.append_sample_name(nd251_2_adata, 'WF40_ND21-251_HSC_3k-2')
    nd251_adata = ad.concat([nd251_1_adata, nd251_2_adata], axis=0, join='outer')
    nd251_adata.write_h5ad('5801-AML.h5ad')

    p5801_aml_therapy_adata = junc.exportJunctionMatrix(matrix_5801A_therapy,ensembl_exon_dir,gff_5801A_therapy,barcode_sample_dict['MDS5801-1229IR_CD34'],rev=True)
    p5801_aml_post_adata = junc.exportJunctionMatrix(matrix_5801A_post,ensembl_exon_dir,gff_5801A_post,barcode_clusters=barcode_sample_dict['MDS5801R-341'],rev=True)

    am75_post_adata = junc.exportJunctionMatrix(matrix_am75M_post,ensembl_exon_dir,gff_source=gff_am75M_post,barcode_clusters=barcode_sample_dict['Y5116'],rev=True)
    am75_pre_adata = junc.exportJunctionMatrix(matrix_am75M_pre,ensembl_exon_dir,gff_source=gff_am75M_pre,barcode_clusters=barcode_sample_dict['Y1258'],rev=True)
    am75_on_adata = junc.exportJunctionMatrix(matrix_am75M_therapy,ensembl_exon_dir,gff_source=gff_am75M_therapy,barcode_clusters=barcode_sample_dict['Y1607'],rev=True)
    wf78_post_adata = junc.exportJunctionMatrix(matrix_wf78M_post,ensembl_exon_dir,gff_source=gff_wf78M_post,barcode_clusters=barcode_sample_dict['Y1737'],rev=True)
    wf78_pre_adata = junc.exportJunctionMatrix(matrix_wf78M_pre,ensembl_exon_dir,gff_source=gff_wf78M_pre,barcode_clusters=barcode_sample_dict['Y1493'],rev=True)

aggregate_junction_data = True
if aggregate_junction_data:
    # Define the list of sample files
    sample_files = [
        'WM34.h5ad','BM27.h5ad','BF21.h5ad','ND251-apharesis.h5ad','ND167.h5ad',
        'BF71.h5ad','WF83.h5ad','AM72.h5ad',
        '5801-pre','5801-post', '5801-AML', '5801IR','5801-9_341',
        'Y1258','Y1607','Y5116','Y1493','Y1737'
    ]

    pseudo_counts,pseudo_tpm,pseudo_ratios = iso.concatenate_h5ad_and_compute_pseudobulks(sample_files,collection_name='Old-Young')

    # Filter and re-order exported pseudobulks
    input_file = pseudo_counts
    output_file = 'filtered_combined_junction_counts.txt'
    groups_file = 'groups_file2.txt'
    cell_type_order_file = 'cell_type_order.txt'  # Optional
    iso.export_and_filter_pseudobulks(input_file, output_file, groups_file, cell_type_order_file)


    #### Run AltAnalyze2 with min_exp=4;med_exp=9, as
    # python ExpressionBuilder.py --i /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/ --a unbiased --species Hs --platform RNASeq --additional /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/ExpressionInput/counts.ND-junctions.txt
    # python import_scripts/AugmentEventAnnotations.py --i /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI.txt --s Hs --platform RNASeq
    # python stats_scripts/metaDataAnalysis.py --i /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt --s Hs --p PSI --dPSI 0.1 --percentExp 40 --m /Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/pseudocounts/AltAnalyze2/ExpressionInput/groups.ND-junctions.txt

################################### ISOFORM ANNOTATION ###################################

build_healthy_iso_ref = False
if build_healthy_iso_ref:
    #### Create a combined GFF and find the longest consensus isoforms
    gff_source =  [gencode_gff,gff_RefSeq,gff_GTEx]
    gff_source += [gff_wm34,gff_bm27,gff_bf21,gff_167h,gff_167m,gff_167l,gff_251_1,gff_251_2]
    gff_source += [gff_am72,gff_wf83,gff_bf71]
    gff_source += [gff_5801M_pre,gff_5801M_post,gff_5801A_pre1,gff_5801A_pre2,gff_5801A_therapy,gff_5801A_post]
    gff_source += [gff_am75M_pre,gff_am75M_therapy,gff_am75M_post,gff_am75M_therapy,gff_am75M_post,gff_wf78M_pre,gff_wf78M_post]
    gff_output_dir = gff_process.consolidateLongReadGFFs(gff_source, ensembl_exon_dir,mode='Ensembl')

    #### Translate isoform into protein sequence and check for NMD
    ref_gff_file = "/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/Accessory/AltAnalyze3/gencode.v45.annotation.gff3"
    query_gff_file = gff_output_dir+"/combined.gff"
    genome_fasta = "/Users/saljh8/Dropbox/Revio/Other/Variants/SNV/genome.fa"
    transcript_associations_file = gff_output_dir+"/transcript_associations.txt"
    cds_records, transcript_records, protein_records = isot.gff_translate(query_gff_file,genome_fasta,ref_gff_file,transcript_associations_file)

    with open("protein_sequences.fasta", "w") as protein_file:
        SeqIO.write(protein_records, protein_file, "fasta")

    with open("cds_sequences.fasta", "w") as cds_file:
        SeqIO.write(cds_records, cds_file, "fasta")
    sys.exit()

integrate_healthy_malignant_isoforms = False
if integrate_healthy_malignant_isoforms:
    print("Integrating malignant and healthy for pseudobulk cell populations...")

    ######## Process data from Young Bone Marrow samples ########
    isoform_association_path = "/Users/saljh8/Dropbox/Revio/BoneMarrow/gff-output/isoform_links.txt"

    wm34_adata = isor.exportConsensusIsoformMatrix(matrix_wm34,isoform_association_path,gff_source=gff_wm34,barcode_clusters=barcode_sample_dict['WM34_CD34'],rev=True)
    bf21_adata = isor.exportConsensusIsoformMatrix(matrix_bf21,isoform_association_path,gff_source=gff_bf21,barcode_clusters=barcode_sample_dict['BF21-CD34'],rev=True)
    bm27_adata = isor.exportConsensusIsoformMatrix(matrix_bm27,isoform_association_path,gff_source=gff_bm27,barcode_clusters=barcode_sample_dict['BM27_CD34'],rev=True)
    hsc_nd167_adata = isor.exportConsensusIsoformMatrix(matrix_167h,isoform_association_path,gff_source=gff_167h,barcode_clusters=barcode_sample_dict['WM35_ND20-167__HSC_3k'],sample_id='ND20-167')
    mpp_nd167_adata = isor.exportConsensusIsoformMatrix(matrix_167m,isoform_association_path,gff_source=gff_167m,barcode_clusters=barcode_sample_dict['WM35_ND20-167__MPP_3k'],sample_id='ND20-167')
    lmpp_nd167_adata = isor.exportConsensusIsoformMatrix(matrix_167l,isoform_association_path,gff_source=gff_167l,barcode_clusters=barcode_sample_dict['WM35_ND20-167__LMPP_3k'],sample_id='ND20-167')
    hsc_nd167_adata = iso.append_sample_name(hsc_nd167_adata, 'WM35_ND20-167__HSC_3k')
    mpp_nd167_adata = iso.append_sample_name(mpp_nd167_adata, 'WM35_ND20-167__MPP_3k')
    lmpp_nd167_adata = iso.append_sample_name(lmpp_nd167_adata, 'WM35_ND20-167__LMPP_3k')
    nd167_adata = ad.concat([hsc_nd167_adata, mpp_nd167_adata, lmpp_nd167_adata], axis=0, join='outer')
    nd167_adata.write_h5ad('ND167-isoform.h5ad')
    nd251_1_adata = isor.exportConsensusIsoformMatrix(matrix_251_1,isoform_association_path,gff_source=gff_251_1,rev=True,barcode_clusters=barcode_sample_dict['WF40_ND21-251_HSC_3k'],sample_id='ND251')
    nd251_2_adata = isor.exportConsensusIsoformMatrix(matrix_251_1,isoform_association_path,gff_source=gff_251_2,rev=True,barcode_clusters=barcode_sample_dict['WF40_ND21-251_HSC_3k'],sample_id='ND251')
    nd251_1_adata = iso.append_sample_name(nd251_1_adata, 'WF40_ND21-251_HSC_3k-1')
    nd251_2_adata = iso.append_sample_name(nd251_2_adata, 'WF40_ND21-251_HSC_3k-2')
    nd251_adata = ad.concat([nd251_1_adata, nd251_2_adata], axis=0, join='outer')
    nd251_adata.write_h5ad('ND251-apharesis-isoform.h5ad')

    ######## Process data from Elderly Bone Marrow samples ########

    bf71_adata = isor.exportConsensusIsoformMatrix(matrix_bf71,isoform_association_path,gff_source=gff_bf71,barcode_clusters=barcode_sample_dict['BF71-N-0421_CD34'],rev=True)



    ######## Process data from elderly Bone Marrow samples ########
    barcode_cluster_dir = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/Illumina/SoupX-hg38/cellHarmony/QueryGroups.cellHarmony.txt'
    barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dir)
    
    # Produce an junction counts h5ad for the AM72-N-0206_CD34
    matrix_path = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/AM72/scisoseq.seurat_info/isoforms_seurat'
    gff_source = 'AM72'
    am72_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=barcode_sample_dict['AM72-N-0206_CD34'],rev=True)

    # Produce an junction counts h5ad for the WF83-N-0215_CD34 capture
    matrix_path = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/PacBio/WF83/scisoseq.seurat_info/isoforms_seurat'
    gff_source = 'WF83'
    wf83_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=barcode_sample_dict['WF83-N-0215_CD34'],rev=True)

    barcode_cluster_dir = '/Users/saljh8/Dropbox/Revio/MDS-Revio/Aggregate/QueryGroups.cellHarmony.txt'
    barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dir)
    
    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML/run1/isoforms_seurat'
    gff_source = '5801-AML-1'
    p5801_aml_adata_1 = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=barcode_sample_dict['MDS5801-preVEN_342'],rev=True,sample_id='5801')
    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML/run2/isoforms_seurat'
    gff_source = '5801-AML-2'
    p5801_aml_adata_2 = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=barcode_sample_dict['MDS5801-preVEN_342'],rev=True,sample_id='5801')
    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML-therapy/isoforms_seurat'
    gff_source = '5801IR'
    p5801_aml_therapy_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,barcode_clusters=barcode_sample_dict['MDS5801-1229IR_CD34'],rev=True,sample_id='5801')
    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML-post/isoforms_seurat'
    gff_source = '5801-9_341'
    p5801_aml_post_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['MDS5801R-341'],sample_id='5801')
    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/MDS-post/isoforms_seurat'
    gff_source = '5801-post'
    p5801_mds_post_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['MDS5801-postHMA_CD34'],sample_id='5801')
    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/MDS-pre/isoforms_seurat'
    gff_source = '5801-pre'
    p5801_mds_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['WM71-0121-5801_CD34'],sample_id='5801')
    p5801_aml_adata_1 = iso.append_sample_name(p5801_aml_adata_1, '5801-AML-1')
    p5801_aml_adata_2 = iso.append_sample_name(p5801_aml_adata_2, '5801-AML-2')
    p5801_aml_therapy_adata = iso.append_sample_name(p5801_aml_therapy_adata, '5801IR')
    p5801_aml_post_adata = iso.append_sample_name(p5801_aml_post_adata, '5801-9_341')
    p5801_mds_post_adata = iso.append_sample_name(p5801_mds_post_adata, '5801-post')
    p5801_mds_adata = iso.append_sample_name(p5801_mds_adata, '5801-pre')
    p5801_adata = ad.concat([p5801_aml_adata_1, p5801_aml_adata_2,p5801_aml_therapy_adata,p5801_aml_post_adata,p5801_mds_post_adata,p5801_mds_adata], axis=0, join='outer')
    p5801_adata.write_h5ad('5801.h5ad')

    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/post/isoforms_seurat'
    gff_source = 'Y5116' #AM75-MDS-post
    am75_post_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['Y5116'],sample_id='AM75')
    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/pre/isoforms_seurat'
    gff_source = 'Y1258' # AM75-MDSpre
    am75_pre_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['Y1258'],sample_id='AM75')
    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/on/isoforms_seurat'
    gff_source = 'Y1607' #AM75-MDSpre
    am75_on_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['Y1607'],sample_id='AM75')
    am75_post_adata = iso.append_sample_name(am75_post_adata, 'Y5116')
    am75_pre_adata = iso.append_sample_name(am75_pre_adata, 'Y1258')
    am75_on_adata = iso.append_sample_name(am75_on_adata, 'Y1607')
    am75_adata = ad.concat([am75_post_adata, am75_pre_adata, am75_on_adata], axis=0, join='outer')
    am75_adata.write_h5ad('am75.h5ad')

    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/WF78-MDS/post/isoforms_seurat'
    gff_source = 'Y1737' #WF78-MDS-post
    wf78_post_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['Y1737'],sample_id='WF78')
    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/WF78-MDS/pre/isoforms_seurat'
    gff_source = 'Y1493' #WF78-MDS-pre
    wf78_pre_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['Y1493'],sample_id='WF78')
    wf78_post_adata = iso.append_sample_name(wf78_post_adata, 'Y1737')
    wf78_pre_adata = iso.append_sample_name(wf78_pre_adata, 'Y1493')
    wf78_adata = ad.concat([wf78_post_adata, wf78_pre_adata], axis=0, join='outer')
    wf78_adata.write_h5ad('wf78.h5ad')

    barcode_cluster_dir = '/Users/saljh8/Dropbox/Revio/AML-Revio/QueryGroups.cellHarmony.txt'
    barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dir)
    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML14/run1/isoforms_seurat'
    gff_source = 'AML14'
    aml14_1_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['AML-14_CITE_GEX'],sample_id='AML14')
    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML14/run2/isoforms_seurat'
    gff_source = 'AML14-2'
    aml14_2_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['AML-14_CITE_GEX'],sample_id='AML14')
    aml14_1_adata = iso.append_sample_name(aml14_1_adata, 'AML14')
    aml14_2_adata = iso.append_sample_name(aml14_2_adata, 'AML14-2')
    aml14_adata = ad.concat([aml14_1_adata, aml14_2_adata], axis=0, join='outer')
    aml14_adata.write_h5ad('aml14.h5ad')

    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML12/12-run1/isoforms_seurat'
    gff_source = 'AML12'
    aml12_1_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['AML-12_1_CITE_GEX'],sample_id='AML12')
    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML12/12-run2/isoforms_seurat'
    gff_source = 'AML12-2'
    aml12_2_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['AML-12_1_CITE_GEX'],sample_id='AML12')
   
    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML12/13-run1/isoforms_seurat'
    gff_source = 'AML13'
    aml12_3_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['AML-13_CITE_GEX'],sample_id='AML12')
    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML12/13-run2/isoforms_seurat'
    gff_source = 'AML13-2'
    aml12_4_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['AML-13_CITE_GEX'],sample_id='AML12')
    aml12_1_adata = iso.append_sample_name(aml12_1_adata, 'AML12')
    aml12_2_adata = iso.append_sample_name(aml12_2_adata, 'AML12-2')
    aml12_3_adata = iso.append_sample_name(aml12_3_adata, 'AML13')
    aml12_4_adata = iso.append_sample_name(aml12_4_adata, 'AML13-2')
    aml12_adata = ad.concat([aml12_1_adata, aml12_2_adata,aml12_3_adata,aml12_4_adata], axis=0, join='outer')
    aml12_adata.write_h5ad('aml12.h5ad')

    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML7/run1/isoforms_seurat'
    gff_source = 'AML7'
    aml7_1_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['AML-7_1_CITE_GEX'],sample_id='AML7')
    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML7/run2/isoforms_seurat'
    gff_source = 'AML7-2'
    aml7_2_adata = isor.exportConsensusIsoformMatrix(matrix_path,isoform_association_path,gff_source=gff_source,rev=True,barcode_clusters=barcode_sample_dict['AML-7_1_CITE_GEX'],sample_id='AML7')
    aml7_1_adata = iso.append_sample_name(aml7_1_adata, 'AML7')
    aml7_2_adata = iso.append_sample_name(aml7_2_adata, 'AML7-2')
    aml7_adata = ad.concat([aml7_1_adata, aml7_2_adata], axis=0, join='outer')
    aml7_adata.write_h5ad('aml7.h5ad')
    
    am75_pre_adata = ad.read_h5ad('Y1258-isoform.h5ad')
    am75_post_adata = ad.read_h5ad('Y5116-isoform.h5ad')
    am75_post_adata = iso.append_sample_name(am75_post_adata, 'Y5116')
    am75_pre_adata = iso.append_sample_name(am75_pre_adata, 'Y1258')
    am75_adata = ad.concat([am75_pre_adata, am75_post_adata], axis=0, join='outer')
    am75_adata.write_h5ad('am75-pre-post.h5ad')
    exportRatios(am75_adata,'am75-pre-post')
    
    
    p5801_aml_post_adata = ad.read_h5ad('5801-9_341-isoform.h5ad')
    p5801_mds_post_adata = ad.read_h5ad('5801-post-isoform.h5ad')
    p5801_mds_adata = ad.read_h5ad('5801-pre-isoform.h5ad')
    p5801_aml_therapy_adata = ad.read_h5ad('5801IR-isoform.h5ad')
    p5801_aml_adata_2 = ad.read_h5ad('5801-AML-2-isoform.h5ad')
    p5801_aml_adata_1 = ad.read_h5ad('5801-AML-1-isoform.h5ad')

    p5801_aml_adata_1 = iso.append_sample_name(p5801_aml_adata_1, '5801-AML-1')
    p5801_aml_adata_2 = iso.append_sample_name(p5801_aml_adata_2, '5801-AML-2')
    p5801_aml_therapy_adata = iso.append_sample_name(p5801_aml_therapy_adata, '5801IR')
    p5801_aml_post_adata = iso.append_sample_name(p5801_aml_post_adata, '5801-9_341')
    p5801_mds_post_adata = iso.append_sample_name(p5801_mds_post_adata, '5801-post')
    p5801_mds_adata = iso.append_sample_name(p5801_mds_adata, '5801-pre')
    p5801_adata = ad.concat([p5801_aml_therapy_adata,p5801_aml_post_adata], axis=0, join='outer')
    p5801_adata.write_h5ad('5801-AML-therapy.h5ad')
    exportRatios(p5801_adata,'5801-AML-therapy')

    #nd167_adata = ad.read_h5ad('ND167-isoform.h5ad');print('ND167-isoform.h5ad')
    #nd251_adata = ad.read_h5ad('ND251-apharesis-isoform.h5ad')

    """

    am72_adata = ad.read_h5ad('AM72-isoform.h5ad'); print('AM72-isoform.h5ad')
    wf83_adata = ad.read_h5ad('WF83-isoform.h5ad'); print('WF83-isoform.h5ad')
    wm34_adata = ad.read_h5ad('WM34-isoform.h5ad')
    bf71_adata = ad.read_h5ad('BF71-isoform.h5ad')
    bf21_adata = ad.read_h5ad('BF21-isoform.h5ad')
    bm27_adata = ad.read_h5ad('BM27-isoform.h5ad')

    # Make sure the indices are strings
    wm34_adata.obs.index = wm34_adata.obs.index.astype(str)
    bf21_adata.obs.index = bf21_adata.obs.index.astype(str)
    bf71_adata.obs.index = bf71_adata.obs.index.astype(str)
    bm27_adata.obs.index = bm27_adata.obs.index.astype(str)
    wf83_adata.obs.index = wf83_adata.obs.index.astype(str)
    am72_adata.obs.index = am72_adata.obs.index.astype(str)

    # Make observation names unique to avoid duplication issues
    wm34_adata.obs_names_make_unique()
    bf21_adata.obs_names_make_unique()
    bf71_adata.obs_names_make_unique()
    bm27_adata.obs_names_make_unique()
    wf83_adata.obs_names_make_unique()
    am72_adata.obs_names_make_unique()
    
    #### WARNING - MEMORY INTENSIVE - NOT RECOMMENDED - USE BELOW
    combined_adata = ad.concat([wm34_adata,bf21_adata,bm27_adata,bf71_adata,am72_adata, wf83_adata], axis=0, join='outer')

    pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_cluster_counts(combined_adata,cell_threshold=5,count_threshold=0,compute_tpm=True)
    pseudo_pdf.to_csv('ND-MDS-AML-combined_pseudo_cluster_counts.txt', sep='\t')
    tpm.to_csv('ND-MDS-AM-combined_pseudo_cluster-tpm.txt', sep='\t')
    isoform_to_gene_ratio.to_csv('ND-MDS-AM-combined_pseudo_cluster-ratio.txt', sep='\t')
    """

    # Define the list of sample files
    sample_files = [
        'WM34-isoform.h5ad',
        'BF21-isoform.h5ad',
        'BF71-isoform.h5ad',
        'BM27-isoform.h5ad',
        'WF83-isoform.h5ad',
        'AM72-isoform.h5ad',
        'ND251-apharesis-isoform.h5ad',
        'ND167-isoform.h5ad'
    ]

    # Memory efficient combine and pseudobulk of many h5ad files
    pseudo_counts,pseudo_tpm,pseudo_ratios = concatenate_h5ad_and_compute_pseudobulks(sample_files,collection_name='ND-MDS-AML')

    # Filter and re-order exported pseudobulks
    input_file = pseudo_tpm
    input_file = pseudo_ratios
    input_file = 'filtered_ratio-filtered.txt'
    output_file = 'filtered_ratio2.txt'
    groups_file = 'groups_file2.txt'
    cell_type_order_file = 'cell_type_order.txt'  # Optional

    iso.export_and_filter_pseudobulks(input_file, output_file, groups_file, cell_type_order_file)




run_healthy_marrow_iso = False
if run_healthy_marrow_iso:

    ################################### ISOFORM LEVEL ANALYSIS ###################################
    #### Further collapse the isoform matrix for each sample using multi-sample consensus IDs

    # Indicates the "best" isoform across smaples
    isoform_association_path = "/Users/saljh8/Dropbox/Revio/BoneMarrow/gff-output/isoform_links.txt"

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
    matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/LMPP_3k/filtered_matrices'
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
    pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_cluster_counts(nd167_adata,cell_threshold=5,count_threshold=0,compute_tpm=True)
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






integrate_healthy_malignant_junction = False
if integrate_healthy_malignant_junction:
    isoform_association_path = "/Users/saljh8/Dropbox/Revio/BoneMarrow/gff-output/isoform_links.txt"

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


    #### AML and MDS

    gff_a1 = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML/run1/5801-AML-1.gff'
    gff_a2 ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML/run2/5801-AML-2.gff'
    gff_a3 ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML-post/5801-9_341.gff'
    gff_a4 ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML-therapy/5801IR.gff'
    gff_a5 ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/MDS-post/5801-post.gff'
    gff_a6 ='/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/MDS-pre/5801-pre.gff'
    gff_a7 ='/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/on/Y1607.gff'
    gff_a8 ='/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/pre/Y1258.gff'
    gff_a9 ='/Users/saljh8/Dropbox/Revio/MDS-Revio/WF78-MDS/post/Y1737.gff'
    gff_a10 ='/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/post/Y5116.gff'
    gff_a11 ='/Users/saljh8/Dropbox/Revio/MDS-Revio/WF78-MDS/pre/Y1493.gff'
    gff_a12 ='/Users/saljh8/Dropbox/Revio/AML-Revio/AML7/run1/AML7.gff'
    gff_a13 ='/Users/saljh8/Dropbox/Revio/AML-Revio/AML7/run2/AML7-2.gff'
    gff_a14 ='/Users/saljh8/Dropbox/Revio/AML-Revio/AML12/12-run1/AML12.gff'
    gff_a15 ='/Users/saljh8/Dropbox/Revio/AML-Revio/AML12/12-run1/AML12.gff'
    gff_a16 ='/Users/saljh8/Dropbox/Revio/AML-Revio/AML12/12-run2/AML12-2.gff'
    gff_a17 ='/Users/saljh8/Dropbox/Revio/AML-Revio/AML12/13-run1/AML13.gff'
    gff_a18 ='/Users/saljh8/Dropbox/Revio/AML-Revio/AML12/13-run2/AML13-2.gff'
    gff_a19 ='/Users/saljh8/Dropbox/Revio/AML-Revio/AML14/run1/AML14.gff'
    gff_a20 ='/Users/saljh8/Dropbox/Revio/AML-Revio/AML14/run2/AML14-2.gff'

    barcode_cluster_dir = '/Users/saljh8/Dropbox/Revio/Aged-Healthy-Revio/Illumina/SoupX-hg38/cellHarmony/QueryGroups.cellHarmony.txt'
    barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dir)

    # Produce an junction counts h5ad for the 5801 captures
    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML/run1/isoforms_seurat'
    p5801_aml_adata_1 = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_a1,rev=True,barcode_clusters=barcode_sample_dict['MDS5801-preVEN_342'])
    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML/run2/isoforms_seurat'
    p5801_aml_adata_2 = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_a2,rev=True,barcode_clusters=barcode_sample_dict['MDS5801-preVEN_342'])

    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML-therapy/isoforms_seurat'
    p5801_aml_therapy_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_a4,rev=True,sample_id='5801')

    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/AML-post/isoforms_seurat'
    p5801_aml_post_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_a3,rev=True,sample_id='5801')

    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/MDS-post/isoforms_seurat'
    p5801_mds_post_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_a5,rev=True,sample_id='5801')

    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/5801-MDS/MDS-pre/isoforms_seurat'
    p5801_mds_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_a6,rev=True,sample_id='5801')

    p5801_aml_adata_1 = iso.append_sample_name(p5801_aml_adata_1, '5801-AML-1')
    p5801_aml_adata_2 = iso.append_sample_name(p5801_aml_adata_2, '5801-AML-2')
    p5801_aml_therapy_adata = iso.append_sample_name(p5801_aml_therapy_adata, '5801IR')
    p5801_aml_post_adata = iso.append_sample_name(p5801_aml_post_adata, '5801-9_341')
    p5801_mds_post_adata = iso.append_sample_name(p5801_mds_post_adata, '5801-post')
    p5801_mds_adata = iso.append_sample_name(p5801_mds_adata, '5801-pre')
    p5801_adata = ad.concat([p5801_aml_adata_1, p5801_aml_adata_2,p5801_aml_therapy_adata,p5801_aml_post_adata,p5801_mds_post_adata,p5801_mds_adata], axis=0, join='outer')
    p5801_adata.write_h5ad('5801_junction.h5ad')

    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/post/isoforms_seurat'
    gff_source = 'Y5116'
    am75_post_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_source=gff_source,rev=True,sample_id='AM75')
    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/pre/isoforms_seurat'
    gff_source = 'Y1258'
    am75_pre_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_source=gff_source,rev=True,sample_id='AM75')
    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/AM75-MDS/on/isoforms_seurat'
    gff_source = 'Y1607'
    am75_on_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_source=gff_source,rev=True,sample_id='AM75')
    am75_post_adata = iso.append_sample_name(am75_post_adata, 'Y5116')
    am75_pre_adata = iso.append_sample_name(am75_pre_adata, 'Y1258')
    am75_on_adata = iso.append_sample_name(am75_on_adata, 'Y1607')
    am75_adata = ad.concat([am75_post_adata, am75_pre_adata, am75_on_adata], axis=0, join='outer')
    am75_adata.write_h5ad('am75.h5ad')

    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/WF78-MDS/post/isoforms_seurat'
    gff_source = 'Y1737'
    wf78_post_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_source=gff_source,rev=True,sample_id='WF78')
    matrix_path = '/Users/saljh8/Dropbox/Revio/MDS-Revio/WF78-MDS/pre/isoforms_seurat'
    gff_source = 'Y1493'
    wf78_pre_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_source=gff_source,rev=True,sample_id='WF78')
    wf78_post_adata = iso.append_sample_name(wf78_post_adata, 'Y1737')
    wf78_pre_adata = iso.append_sample_name(wf78_pre_adata, 'Y1493')
    wf78_adata = ad.concat([wf78_post_adata, wf78_pre_adata], axis=0, join='outer')
    wf78_adata.write_h5ad('wf78.h5ad')

    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML14/run1/isoforms_seurat'
    gff_source = 'AML14'
    aml14_1_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_source=gff_source,rev=True,sample_id='AML14')
    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML14/run2/isoforms_seurat'
    gff_source = 'AML14-2'
    aml14_2_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_source=gff_source,rev=True,sample_id='AML14')
    aml14_1_adata = iso.append_sample_name(aml14_1_adata, 'AML14')
    aml14_2_adata = iso.append_sample_name(aml14_2_adata, 'AML14-2')
    aml14_adata = ad.concat([aml14_1_adata, aml14_2_adata], axis=0, join='outer')
    aml14_adata.write_h5ad('aml14.h5ad')

    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML12/12-run1/isoforms_seurat'
    gff_source = 'AML12'
    aml12_1_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_source=gff_source,rev=True,sample_id='AML12')
    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML12/12-run2/isoforms_seurat'
    gff_source = 'AML12-2'
    aml12_2_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_source=gff_source,rev=True,sample_id='AML12')
    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML12/13-run1/isoforms_seurat'
    gff_source = 'AML13'
    aml12_3_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_source=gff_source,rev=True,sample_id='AML12')
    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML12/13-run2/isoforms_seurat'
    gff_source = 'AML13-2'
    aml12_4_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_source=gff_source,rev=True,sample_id='AML12')
    aml12_1_adata = iso.append_sample_name(aml12_1_adata, 'AML12')
    aml12_2_adata = iso.append_sample_name(aml12_2_adata, 'AML12-2')
    aml12_3_adata = iso.append_sample_name(aml12_3_adata, 'AML13')
    aml12_4_adata = iso.append_sample_name(aml12_4_adata, 'AML13-2')
    aml12_adata = ad.concat([aml12_1_adata, aml12_2_adata,aml12_3_adata,aml12_4_adata], axis=0, join='outer')
    aml12_adata.write_h5ad('aml12.h5ad')

    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML7/run1/isoforms_seurat'
    gff_source = 'AML7'
    aml7_1_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_source=gff_source,rev=True,sample_id='AML7')
    matrix_path = '/Users/saljh8/Dropbox/Revio/AML-Revio/AML7/run2/isoforms_seurat'
    gff_source = 'AML7-2'
    aml7_2_adata = junc.exportJunctionMatrix(matrix_path,ensembl_exon_dir,gff_source=gff_source,rev=True,sample_id='AML7')
    aml7_1_adata = iso.append_sample_name(aml7_1_adata, 'AML7')
    aml7_2_adata = iso.append_sample_name(aml7_2_adata, 'AML7-2')
    aml7_adata = ad.concat([aml7_1_adata, aml7_2_adata], axis=0, join='outer')
    aml7_adata.write_h5ad('aml7.h5ad')

    am75_pre_adata = ad.read_h5ad('Y1258-isoform.h5ad')
    am75_post_adata = ad.read_h5ad('Y5116-isoform.h5ad')
    am75_post_adata = iso.append_sample_name(am75_post_adata, 'Y5116')
    am75_pre_adata = iso.append_sample_name(am75_pre_adata, 'Y1258')
    am75_adata = ad.concat([am75_pre_adata, am75_post_adata], axis=0, join='outer')
    am75_adata.write_h5ad('am75-pre-post.h5ad')
    exportRatios(am75_adata,'am75-pre-post')

    
    p5801_aml_post_adata = ad.read_h5ad('5801-9_341-isoform.h5ad')
    p5801_mds_post_adata = ad.read_h5ad('5801-post-isoform.h5ad')
    p5801_mds_adata = ad.read_h5ad('5801-pre-isoform.h5ad')
    p5801_aml_therapy_adata = ad.read_h5ad('5801IR-isoform.h5ad')
    p5801_aml_adata_2 = ad.read_h5ad('5801-AML-2-isoform.h5ad')
    p5801_aml_adata_1 = ad.read_h5ad('5801-AML-1-isoform.h5ad')

    p5801_aml_adata_1 = iso.append_sample_name(p5801_aml_adata_1, '5801-AML-1')
    p5801_aml_adata_2 = iso.append_sample_name(p5801_aml_adata_2, '5801-AML-2')
    p5801_aml_therapy_adata = iso.append_sample_name(p5801_aml_therapy_adata, '5801IR')
    p5801_aml_post_adata = iso.append_sample_name(p5801_aml_post_adata, '5801-9_341')
    p5801_mds_post_adata = iso.append_sample_name(p5801_mds_post_adata, '5801-post')
    p5801_mds_adata = iso.append_sample_name(p5801_mds_adata, '5801-pre')
    p5801_adata = ad.concat([p5801_aml_therapy_adata,p5801_aml_post_adata], axis=0, join='outer')
    p5801_adata.write_h5ad('5801-AML-therapy.h5ad')
    exportRatios(p5801_adata,'5801-AML-therapy')

    nd167_adata = ad.read_h5ad('ND167-isoform.h5ad')
    nd251_adata = ad.read_h5ad('ND251-apharesis-isoform.h5ad')
    am72_adata = ad.read_h5ad('AM72-isoform.h5ad')
    wf83_adata = ad.read_h5ad('WF83-isoform.h5ad')
    p5801_adata = ad.read_h5ad('5801.h5ad')
    am75_adata = ad.read_h5ad('am75.h5ad')
    wf78_adata = ad.read_h5ad('wf78.h5ad')
    aml14_adata = ad.read_h5ad('aml14.h5ad')
    aml12_adata = ad.read_h5ad('aml12.h5ad')
    aml7_adata = ad.read_h5ad('aml12.h5ad')
    
    a_adata = ad.read_h5ad('5801-9_341-isoform.h5ad')
    b_adata = ad.read_h5ad('5801-post-isoform.h5ad')
    c_adata = ad.read_h5ad('5801-pre-isoform.h5ad')
    d_adata = ad.read_h5ad('5801IR-isoform.h5ad')
    e_adata = ad.read_h5ad('5801-AML-2-isoform.h5ad')
    f_adata = ad.read_h5ad('5801-AML-1-isoform.h5ad')

    exportRatios(a_adata,'5801-9_341')
    exportRatios(b_adata,'5801-post')
    exportRatios(c_adata,'5801-pre')
    exportRatios(d_adata,'5801IR')
    exportRatios(e_adata,'5801-AML-2')
    exportRatios(f_adata,'5801-AML-1')

    exportRatios(nd167_adata,'ND167')
    exportRatios(nd251_adata,'ND251')
    exportRatios(am72_adata,'AM72')
    exportRatios(wf83_adata,'WF83')
    exportRatios(p5801_adata,'5801')
    exportRatios(am75_adata,'am75')
    exportRatios(wf78_adata,'wf78')
    exportRatios(aml14_adata,'aml14')
    exportRatios(aml12_adata,'aml12')
    exportRatios(aml7_adata,'aml12')
    
    combined_adata = ad.concat([nd167_adata, nd251_adata, am72_adata, wf83_adata, p5801_adata, am75_adata, wf78_adata, aml14_adata, aml12_adata, aml7_adata], axis=0, join='outer')
    pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_cluster_counts(combined_adata,cell_threshold=5,count_threshold=0,compute_tpm=True)
    pseudo_pdf.to_csv('ND-MDS-AML-combined_pseudo_cluster_counts.txt', sep='\t')
    tpm.to_csv('ND-MDS-AM-combined_pseudo_cluster-tpm.txt', sep='\t')
    isoform_to_gene_ratio.to_csv('ND-MDS-AM-combined_pseudo_cluster-ratio.txt', sep='\t')

old_young_junction_integrate = False
if old_young_junction_integrate:
    ################################### JUNCTION LEVEL ANALYSIS ###################################
    ######## Process data from sample ND167 ########
    

    start_time = time.time()
    print(("Finding top " + str(top_n_differential_events) + " differential splicing events per NMF cluster..."))
    diff_events, remove_clusters, results_all = find_top_differential_events(exp_file=full_psi_file, groups_file=groups_file,
                                                                             metadata=metadata,
                                                                             output_loc=os.path.join(cd, filename),
                                                                             min_group_size=min_group_size, dPSI=dPSI,
                                                                             dPSI_p_val=dPSI_p_val,
                                                                             min_differential_events=min_differential_events,
                                                                             top_n=top_n_differential_events)
    elapsed_time = (time.time() - start_time) / 60  # Convert seconds to minutes
    print("--- %.3f minutes ---" % elapsed_time)