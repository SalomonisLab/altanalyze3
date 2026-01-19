import os, shutil
import sys
import csv
import importlib
from tqdm import tqdm
import anndata as ad
import pandas as pd
import numpy as np
from scipy.sparse import vstack, csr_matrix
from scipy.io import mmread
import collections
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from . import isoform_junctions_chunk as junc
from . import isoform_matrix as iso
from . import isoform_ratios as isor
from . import isoform_translation as isot
from . import gff_process as gff_process
from ..psi import psi_single as psi
importlib.reload(psi)
from Bio import SeqIO
import asyncio
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

def exportRatios(ob,out):
    pseudo_pdf, tpm, isoform_to_gene_ratio = iso.pseudo_cluster_counts(ob,cell_threshold=5,count_threshold=0,compute_tpm=True)
    pseudo_pdf.to_csv(out+'_pseudo_cluster_counts.txt', sep='\t')
    tpm.to_csv(out+'_pseudo_cluster-tpm.txt', sep='\t')
    isoform_to_gene_ratio.to_csv(out+'_pseudo_cluster-ratio.txt', sep='\t')

def import_metadata(metadata_file, return_size = False, include_hashed_samples = False):
    """Reads metadata file and groups samples by uid."""
    metadata = pd.read_csv(metadata_file, sep='\t')
    sample_dict = {}
    gff_name_tracker = {}
    hashed_path = None
    gff_names = []
    comparison_groups = {}

    for _, row in metadata.iterrows():
        uid = row['uid']
        library = row['library']
        gff_path = row['gff']
        if 'hashed_barcodes' in row:
            hashed_path = row['hashed_barcodes']
            if isinstance(hashed_path, str) and '.' in hashed_path:
                pass
            else:
                hashed_path = None
        gff_name = os.path.basename(gff_path)  # Full GFF filename
        gff_basename = gff_name[:-4]  # Name without extension

        # Check for duplicate GFF filenames or "." in the filename before the extension
        if gff_name in gff_name_tracker or '.' in gff_basename:
            new_gff_name = f"{library}.gff"
            new_gff_path = os.path.join(os.path.dirname(gff_path), new_gff_name)

            # Copy the original GFF to the new path
            shutil.copy2(gff_path, new_gff_path)
            #print(f"Incompatible gff name... copying to:\n {gff_path} to {new_gff_path}")

            # Update the GFF path to point to the new file
            gff_path = new_gff_path
            gff_basename = library

        # Track the GFF filename (to check for duplicates -  not allowed)
        gff_name_tracker[gff_name] = gff_path
        
        if include_hashed_samples == False:
            # Only include hashed libraries as separate when exporting sample-level h5ad files
            if gff_basename in gff_names:
                continue

        # Populate the sample dictionary
        if uid not in sample_dict:
            sample_dict[uid] = []

        sample_dict[uid].append({
            'gff': gff_path, 
            'gff_name': gff_basename,
            'matrix': row['matrix'], 
            'library': library, 
            'reverse': row['reverse'],
            'groups': row['groups'],
            'hashed': hashed_path
        })
        if row['groups'] in comparison_groups:
            if uid not in comparison_groups[row['groups']]:
                comparison_groups[row['groups']].append(uid)
        else:
            comparison_groups[row['groups']]=[uid]
        
        gff_names.append(gff_basename)

    group_sizes = list(set(len(v) for v in comparison_groups.values()))
    two_groups_with_min_3 = sum(x > 2 for x in group_sizes) >= 2
    if two_groups_with_min_3:
        min_group_size = 3
    else:
        min_group_size = min(group_sizes) if group_sizes else 0

    if return_size:
        print(f"Minimum group size: {min_group_size}")
        return sample_dict, min_group_size
    else:
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

def export_isoform_h5ad(sample_dict, ensembl_exon_dir, barcode_sample_dict, reference_gff, genome_fasta, deleteGFF=False):
    """Processes samples, concatenates if necessary, and writes h5ad files."""

    current_dir = os.getcwd()
    gff_output = "/gff-output/"
    query_gff_file = current_dir+gff_output+"combined.gff"
    transcript_associations_path = current_dir+gff_output+"transcript_associations.txt"
    isoform_associations_path = current_dir+gff_output+"isoform_links.txt"

    if os.path.exists(isoform_associations_path) and deleteGFF == False:
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
        with open("orf_sequences.fasta", "w") as cds_file:
            SeqIO.write(cds_records, cds_file, "fasta")
            
    # Export individual library/sample isoform h5ad
    for uid, samples in sample_dict.items():
        num_samples = len(samples)
        print(f'Processing {num_samples} libraries for sample: {uid}')
        
        adata_list = []
        for s in samples:
            adata = export_isoform_matrix(
                s['matrix'], s['gff_name'], s['library'],
                isoform_associations_path, s['reverse'], barcode_sample_dict)
            # Resolve the index conflict issue:
            if adata.obs.index.name in adata.obs.columns:
                print(f"Conflict found in {s['library']}: index name '{adata.obs.index.name}' is also a column.")
                adata.obs.rename_axis('barcode', inplace=True)  # Rename index to 'barcode' to avoid conflict
            adata.obs_names_make_unique()
            adata_list.append(adata)
            
            # Check the overlaping labels
            from collections import Counter
            labels = [s['library'] for s in samples]
            label_counts = Counter(labels)
            duplicates = {label: count for label, count in label_counts.items() if count > 1}
            if duplicates:
                print("Duplicated labels found:")
                for label, count in duplicates.items():
                    print(f"{label}: {count} occurrences")
            else:
                #print("No duplicated labels found.")
                pass

        """The below code combines h5ads. Instead of combining counts (see alternative code below)
        the same cell barcode is referenced twice with a sample-level suffix for the barcode ID. 
        Combined reads are the union not intersection."""
        sum_reads_per_cell = False

        if num_samples > 1:
            if sum_reads_per_cell:
                # Combine without changing cell names
                combined = ad.concat(adata_list, axis=0, join='outer', index_unique=None)
                # Group by cell barcode and sum counts
                print("Summing counts across identical barcodes...")
                counts_df = combined.to_df()
                summed_counts = counts_df.groupby(combined.obs_names).sum()
                # Create new AnnData with summed counts
                new_obs = combined.obs.groupby(combined.obs_names).first()
                new_adata = ad.AnnData(X=csr_matrix(summed_counts.values), obs=new_obs, var=combined.var.loc[summed_counts.columns])
                print(f"Saving merged and summed AnnData for {uid}")
                new_adata.write_h5ad(f'{uid}_isoform.h5ad')
            else:
                # Combine without the 'keys' parameter to avoid duplicate label issues
                combined_adata = ad.concat(adata_list, axis=0, join='outer', label="sample")
                combined_adata.obs_names_make_unique()
                combined_adata.write_h5ad(f'{uid}_isoform.h5ad')

def get_valid_h5ad(sample_dict,dataType):
    sample_index={}
    sample_files=[]
    uids=[]
    current_dir = os.getcwd()
    for uid, samples in sample_dict.items():
        h5ad = uid + ('.h5ad' if dataType == 'junction' else '-isoform.h5ad')
        if not os.path.exists(f'{current_dir}/{h5ad}'):
            for s in samples:
                h5ad = s['gff_name'] + ('.h5ad' if dataType == 'junction' else '-isoform.h5ad')
                if h5ad not in sample_files:
                    sample_files.append(h5ad)
                sample_index[s['gff_name']] = h5ad
        else:
            for s in samples:
                sample_index[s['gff_name']] = h5ad
                if h5ad not in sample_files:
                    sample_files.append(h5ad)
        uids.append(uid)

    return sample_files,sample_index,uids

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

def deconvolute_libraries_by_hashed_index(sample_dict,sample_h5ads,dataType):
    sample_metadata={}
    update_sample_files=[]
    added_samples=[] 
    for uid, samples in sample_dict.items():
        for s in samples:
            if s['hashed'] is not None:
                sample_metadata.setdefault((s['hashed'], s['gff_name']), []).append(uid)
                if s['gff_name'] not in added_samples:
                    added_samples.append(s['gff_name'])
            else:
                name = sample_h5ads[s['gff_name']]
                if name not in update_sample_files:
                    update_sample_files.append(name)

    for info, samples in sample_metadata.items():
        hashed_path, gff_name = info
        # Read barcode-hashed_sample mapping
        df = pd.read_csv(hashed_path, sep='\t', header=None, names=['barcode', 'hashed_sample'])
        #print("First 5 barcodes in df['barcode']:", df['barcode'].head().tolist())
        df['barcode'] = df['barcode'].str.replace(r'\..*$', '', regex=True).str.strip()
        hashed_sample_dict = df.groupby('hashed_sample')['barcode'].apply(list).to_dict()
    
        # Process each sample
        if gff_name in sample_h5ads:
            print (sample_h5ads[gff_name])
            adata = ad.read_h5ad(sample_h5ads[gff_name])  # Load the corresponding h5ad file - P1-isoform.h5ad
            #print (adata.obs_names)
            #print("adata.obs first 5 rows:\n", adata.obs.head())
            #print("adata.obs columns:\n", adata.obs.columns.tolist())

            for sample in samples:
                barcodes = hashed_sample_dict[sample]
                # Subset and save h5ad files for each hashed sample
                sample_adata = adata[adata.obs_names.isin(barcodes)]  # Subset the data
                # Generate output filename
                if dataType == 'junction':
                    #output_filename = f"{sample}_{gff_name}.h5ad"
                    output_filename = f"{sample}.h5ad"
                else:
                    #output_filename = f"{sample}_{gff_name}-isoform.h5ad"
                    output_filename = f"{sample}-isoform.h5ad"
                output_path = os.path.join(os.getcwd(), output_filename)
                # Save subset h5ad file
                sample_adata.write(output_path,compression='gzip')
                print(f"Saved {output_filename}")
                if output_filename not in update_sample_files:
                    update_sample_files.append(output_filename)

    """
    for file in update_sample_files:
        adata = ad.read_h5ad(file)
        um_barcodes = adata.n_obs
        print (file, um_barcodes)
    """
    return update_sample_files

def export_pseudo_counts(metadata_file,barcode_cluster_dirs,dataType='junction',compute_tpm=False):
    sample_dict, min_group_size = import_metadata(metadata_file, return_size = True, include_hashed_samples = True)
    cluster_order = iso.return_cluster_order(barcode_cluster_dirs)
    # Memory efficient combine and pseudobulk of many h5ad files
    pseudo_counts = dataType+'_combined_pseudo_cluster_counts.txt'
    pseudo_tpm = dataType+'_combined_pseudo_cluster_tpm.txt'
    pseudo_ratios = dataType+'_combined_pseudo_cluster_ratio.txt'
    
    # Deconvolute multi-sample libraries (hashed) - if present
    sample_files,sample_h5ads,uids = get_valid_h5ad(sample_dict,dataType)
    print (sample_files)
    sample_files = deconvolute_libraries_by_hashed_index(sample_dict,sample_h5ads,dataType)
    print (sample_files)

    iso.concatenate_h5ad_and_compute_pseudobulks_optimized(sample_files,collection_name=dataType,compute_tpm=compute_tpm)

    # Organize and filter
    iso.export_and_filter_pseudobulk_chunks(pseudo_counts, pseudo_counts[:-4]+'-filtered.txt', cluster_order, min_group_size=min_group_size)
    if dataType == 'isoform':
        iso.export_and_filter_pseudobulk_chunks(pseudo_ratios, pseudo_ratios[:-4]+'-filtered.txt', cluster_order, min_group_size=min_group_size)
        iso.export_and_filter_pseudobulk_chunks(pseudo_tpm, pseudo_tpm[:-4]+'-filtered.txt', cluster_order, min_group_size=min_group_size)

def pre_process_samples(metadata_file, barcode_cluster_dirs, ensembl_exon_dir):
    # Perform junction quantification across samples
    sample_dict = import_metadata(metadata_file)
    barcode_sample_dict = iso.import_barcode_clusters(barcode_cluster_dirs)
    export_junction_h5ad(sample_dict, ensembl_exon_dir, barcode_sample_dict)
    export_pseudo_counts(metadata_file,barcode_cluster_dirs,'junction')

def combine_processed_samples(metadata_file, barcode_cluster_dirs, ensembl_exon_dir, gencode_gff, genome_fasta):
    sample_dict, min_group_size = import_metadata(metadata_file, return_size = True)
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


def combine_h5ad_files(sample_h5ads, output_file='combined.h5ad'):
    """
    Combine multiple h5ad files into a single file with a union of features.
    
    Parameters:
    -----------
    sample_h5ads : dictionary
        Dictionary of samples (keys) and h5ad filenames (values)
    output_file : str
        Name of the output combined h5ad file
        
    Returns:
    --------
    anndata.AnnData
        Combined AnnData object with union of features
    """

    # Get unique h5ad files
    sample_h5ads = sorted(set(sample_h5ads.values()))
    print(f"Combining {len(sample_h5ads)} h5ad files...")
    
    # Initialize lists to store data
    adatas = []
    all_features = set()
    
    # First pass: read all files and collect unique features
    for h5ad_file in tqdm(sample_h5ads, desc="Reading files"):
        adata = ad.read_h5ad(h5ad_file)
        adatas.append(adata)
        all_features.update(adata.var_names)
    
    # Convert to sorted list for consistent ordering
    all_features = sorted(list(all_features))
    print(f"Total unique features: {len(all_features)}")
    
    # Create a mapping of features to their indices in the combined matrix
    feature_to_combined_idx = {feat: idx for idx, feat in enumerate(all_features)}
    
    # Second pass: align features and combine
    aligned_adatas = []
    for adata in tqdm(adatas, desc="Aligning features"):
        # Create a mapping of current features to their indices
        feature_to_idx = {feat: idx for idx, feat in enumerate(adata.var_names)}
        
        # Create a mapping from old indices to new indices
        old_to_new_idx = {}
        for feat, old_idx in feature_to_idx.items():
            if feat in feature_to_combined_idx:
                old_to_new_idx[old_idx] = feature_to_combined_idx[feat]
        
        # Create a new sparse matrix with all features
        # Use a more efficient approach for large matrices
        rows, cols, data = [], [], []
        
        # Get the sparse matrix in COO format for efficient iteration
        X_coo = adata.X.tocoo()
        
        # Only process non-zero elements
        for i, j, v in zip(X_coo.row, X_coo.col, X_coo.data):
            if j in old_to_new_idx:
                rows.append(i)
                cols.append(old_to_new_idx[j])
                data.append(v)
        
        # Create the new sparse matrix
        new_X = csr_matrix((data, (rows, cols)), shape=(adata.n_obs, len(all_features)))
        
        # Create new var DataFrame with all features
        new_var = pd.DataFrame(index=all_features)
        
        # Copy over any existing var annotations
        for col in adata.var.columns:
            # Initialize with appropriate dtype
            if pd.api.types.is_numeric_dtype(adata.var[col]):
                new_var[col] = pd.Series(index=all_features, dtype=float)
            else:
                new_var[col] = pd.Series(index=all_features, dtype='object')
            
            # Fill in values for existing features
            new_var.loc[adata.var.index, col] = adata.var[col]
        
        # Create new AnnData with aligned features
        new_adata = ad.AnnData(
            X=new_X,
            obs=adata.obs,
            var=new_var,
            uns=adata.uns
        )
        aligned_adatas.append(new_adata)
    
    # Combine all aligned AnnData objects
    print("Combining aligned AnnData objects...")
    combined_adata = ad.concat(aligned_adatas, axis=0, join='outer')
    
    # Write combined file
    print(f"Writing combined file to {output_file}")
    combined_adata.write_h5ad(output_file, compression='gzip')
    
    return combined_adata