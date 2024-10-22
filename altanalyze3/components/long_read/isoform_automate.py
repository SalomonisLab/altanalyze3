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
from . import isoform_junctions_chunk as junc
from . import isoform_matrix as iso
from . import isoform_ratios as isor
from . import isoform_automate as isoa
from . import isoform_translation as isot
from . import gff_process as gff_process
from ..psi import psi_single as psi
from Bio import SeqIO
import asyncio

def exportRatios(ob,out):
    pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_cluster_counts(ob,cell_threshold=5,count_threshold=0,compute_tpm=True)
    pseudo_pdf.to_csv(out+'_pseudo_cluster_counts.txt', sep='\t')
    tpm.to_csv(out+'_pseudo_cluster-tpm.txt', sep='\t')
    isoform_to_gene_ratio.to_csv(out+'_pseudo_cluster-ratio.txt', sep='\t')

def import_metadata(metadata_file):
    """Reads metadata file and groups samples by uid."""
    metadata = pd.read_csv(metadata_file, sep='\t')
    sample_dict = {}
    for _, row in metadata.iterrows():
        uid = row['uid']
        if uid not in sample_dict:
            sample_dict[uid] = []
        sample_dict[uid].append({
            'gff': row['gff'], 
            'gff_name': os.path.basename(row['gff']).split('.')[0],
            'matrix': row['matrix'], 
            'library': row['library'], 
            'reverse': row['reverse'],
            'groups': row['groups']
        })
    return sample_dict

def export_junction_matrix(matrix, gff, library, reverse, ensembl_exon_dir, barcode_sample_dict):
    """Exports junction matrix for a single sample."""
    adata = junc.exportJunctionMatrix(matrix, ensembl_exon_dir, gff, barcode_clusters=barcode_sample_dict[library], rev=reverse)
    return iso.append_sample_name(adata, library)

def export_isoform_matrix(matrix, gff, library, isoform_associations_path, reverse, barcode_sample_dict):
    """Exports an integrated isoform matrix for a single sample."""
    adata = isor.exportConsensusIsoformMatrix(matrix,isoform_associations_path,gff_source=gff,barcode_clusters=barcode_sample_dict[library],rev=reverse)
    return iso.append_sample_name(adata, library)

def export_junction_h5ad(sample_dict, ensembl_exon_dir, barcode_sample_dict):
    """Processes samples, concatenates if necessary, and writes h5ad files."""
    for uid, samples in sample_dict.items():
        num_samples = len(samples)
        print(f'Processing {num_samples} library for sample: {uid}')
        adata_list = [export_junction_matrix(s['matrix'], s['gff'], s['library'], s['reverse'], ensembl_exon_dir, barcode_sample_dict) for s in samples]
        if num_samples>1:
            combined_adata = ad.concat(adata_list, axis=0, join='outer') if len(adata_list) > 1 else adata_list[0]
            combined_adata.write_h5ad(f'{uid}.h5ad', compression='gzip')

def export_isoform_h5ad(sample_dict, ensembl_exon_dir, barcode_sample_dict, reference_gff, genome_fasta):
    """Processes samples, concatenates if necessary, and writes h5ad files."""

    current_dir = os.getcwd()
    gff_output = "/gff-output/"
    query_gff_file = current_dir+gff_output+"combined.gff"
    transcript_associations_path = current_dir+gff_output+"transcript_associations.txt"
    isoform_associations_path = current_dir+gff_output+"isoform_links.txt"
    if os.path.exists(current_dir+gff_output):
        pass
    else:
        # Create a combined GFF
        gff_files = [reference_gff]
        for uid, samples in sample_dict.items():
            for s in samples:
                gff_files.append(s['gff'])
        gff_output_dir = gff_process.consolidateLongReadGFFs(gff_files, ensembl_exon_dir, mode='Ensembl')
        
        # Export isoform translations
        cds_records, transcript_records, protein_records = isot.gff_translate(query_gff_file,genome_fasta,reference_gff,transcript_associations_path)
        with open("protein_sequences.fasta", "w") as protein_file:
            SeqIO.write(protein_records, protein_file, "fasta")
        with open("transcript_sequences.fasta", "w") as cds_file:
            SeqIO.write(transcript_records, cds_file, "fasta")
        with open("cds_sequences.fasta", "w") as cds_file:
            SeqIO.write(cds_records, cds_file, "fasta")

    # Export individual library/sample isoform h5ad
    for uid, samples in sample_dict.items():
        num_samples = len(samples)
        print(f'Processing {num_samples} library for sample: {uid}')
        adata_list = [export_isoform_matrix(s['matrix'], s['gff_name'], s['library'], isoform_associations_path, s['reverse'], barcode_sample_dict) for s in samples]
        if num_samples>1:
            combined_adata = ad.concat(adata_list, axis=0, join='outer') if len(adata_list) > 1 else adata_list[0]
            combined_adata.write_h5ad(f'{uid}_isoform.h5ad')

def get_valid_h5ad(sample_dict,dataType):
    sample_files=[]
    uids=[]
    current_dir = os.getcwd()
    for uid, samples in sample_dict.items():
        h5ad = uid + ('.h5ad' if dataType == 'junction' else '-isoform.h5ad')
        if not os.path.exists(f'{current_dir}/{h5ad}'):
            for s in samples:
                h5ad = s['gff_name'] + ('.h5ad' if dataType == 'junction' else '-isoform.h5ad')
        sample_files.append(h5ad)
        uids.append(uid)
    return sample_files,uids

def get_sample_to_group(sample_dict,dataType):
    sample_group={}
    current_dir = os.getcwd()
    for uid, samples in sample_dict.items():
        sample = uid + ('' if dataType == 'junction' else '-isoform')
        group = samples[0]['groups']
        sample_group[sample]=group
        if not os.path.exists(f'{current_dir}/{sample}'):
            for s in samples:
                sample = s['gff_name'] + ('' if dataType == 'junction' else '-isoform')
                sample_group[sample]=group
    return sample_group

import os

def export_sorted(filename, sort_col, exclude_header=True):
    """
    Efficiently sort a large file without loading everything into memory.
    The file is sorted based on the specified column (0-indexed).
    """
    output_file = filename[:-4] + '-sorted'  # Temporary sorted file
    index = []  # Store tuples of (sort key, offset, line length)

    with open(filename, 'r', encoding='utf-8') as f:
        header = None
        if exclude_header:
            header = f.readline()  # Store header separately
        while True:
            offset = f.tell()  # Track the byte offset of each line
            line = f.readline()
            if not line:
                break  # Stop if end of file is reached
            length = len(line)
            col_value = line.split('\t')[sort_col].strip()  # Extract the sort key
            index.append((col_value, offset, length))  # Add to index
    index.sort()

    with open(output_file, 'w', encoding='utf-8') as o:
        if exclude_header and header:
            o.write(header)  # Write header to the output file
        with open(filename, 'r', encoding='utf-8') as f:
            for col_value, offset, length in index:
                f.seek(offset)  # Move to the line's position
                o.write(f.read(length))  # Write the line to the output
    try:
        os.remove(filename)  # Delete the original file
        os.rename(output_file, filename)  # Rename sorted file to original name
        return filename
    except Exception as e:
        print(f"Error replacing the original file: {e}")
        return output_file

def export_pseudo_counts(metadata_file,barcode_cluster_dirs,dataType='junction',compute_tpm=False):
    sample_dict = isoa.import_metadata(metadata_file)
    sample_files,uids = get_valid_h5ad(sample_dict,dataType)
    cluster_order = iso.return_cluster_order(barcode_cluster_dirs)

    # Memory efficient combine and pseudobulk of many h5ad files
    pseudo_counts = dataType+'_combined_pseudo_cluster_counts.txt'
    pseudo_tpm = dataType+'_combined_pseudo_cluster_tpm.txt'
    pseudo_ratios = dataType+'_combined_pseudo_cluster_ratio.txt'
    
    iso.concatenate_h5ad_and_compute_pseudobulks_optimized(sample_files,collection_name=dataType,compute_tpm=compute_tpm)

    # Organize and filter
    iso.export_and_filter_pseudobulks(pseudo_counts, pseudo_counts[:-4]+'-filtered.txt', cluster_order)
    if dataType == 'isoform':
        iso.export_and_filter_pseudobulks(pseudo_ratios, pseudo_ratios[:-4]+'-filtered.txt', cluster_order)
        iso.export_and_filter_pseudobulks(pseudo_tpm, pseudo_tpm[:-4]+'-filtered.txt', cluster_order)

def pre_process_samples(metadata_file, barcode_cluster_dirs, ensembl_exon_dir):
    # Perform junction quantification across samples
    sample_dict = import_metadata(metadata_file)
    barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dirs)
    export_junction_h5ad(sample_dict, ensembl_exon_dir, barcode_sample_dict)
    export_pseudo_counts(metadata_file,barcode_cluster_dirs,'junction')

def combine_processed_samples(metadata_file, barcode_cluster_dirs, ensembl_exon_dir, gencode_gff, genome_fasta):
    sample_dict = import_metadata(metadata_file)
    barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dirs)

    # Perform isoform quantification across samples
    export_isoform_h5ad(sample_dict, ensembl_exon_dir, barcode_sample_dict, gencode_gff, genome_fasta)
    export_pseudo_counts(metadata_file,barcode_cluster_dirs,'isoform',compute_tpm=True)
    
    junction_coords_file = 'junction_combined_pseudo_cluster_counts-filtered.txt'
    outdir = 'psi_combined_pseudo_cluster_counts.txt'
    junction_coords_file = export_sorted(junction_coords_file, 0) ### Sort the expression file
    #psi.main(junction_path=junction_coords_file, query_gene=None, outdir=outdir, use_multiprocessing=False, mp_context=mp_context, num_cores=num_cores)
    run_psi_analysis(junction_coords_file, outdir)

# Wrap in a function to avoid event loop conflicts
def run_psi_analysis(junction_path, outdir):
    # Use asyncio.run to ensure event loop is properly created
    asyncio.run(psi.main(junction_path=junction_path, query_gene=None, outdir=outdir))
