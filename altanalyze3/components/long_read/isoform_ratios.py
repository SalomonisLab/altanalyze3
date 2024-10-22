import os
import sys
import csv
import anndata as ad
import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy.sparse import lil_matrix
import collections
import argparse
from tqdm import tqdm # progress bar

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from . import isoform_matrix as iso
from . import gff_process as gff_process
    
def exportConsensusIsoformMatrix(matrix_dir, isoform_association_path, gff_source=None, barcode_clusters=None, rev=False, mtx=True, sample_id=None):
    """ 
    Decompose isoforms into exon-exon and exon-junctions counts at the single-cell and pseudobulk level
    """

    # Load the 10x matrix keyed by isoforms in the first column and cell barcode cluster annotations

    if mtx:
        adata = iso.mtx_to_adata(int_folder=matrix_dir, gene_is_index=True, feature='genes.tsv', feature_col=0, barcode='barcodes.tsv', barcode_col=0, matrix='matrix.mtx', rev=rev)
    else:
        adata = iso.tsv_to_adata(matrix_dir)    

    # Ensure obs index and column names are string type
    adata.obs.index = adata.obs.index.astype(str)
    adata.obs.columns = adata.obs.columns.astype(str)

    # If existing cell clusters are provided already (e.g., supervised classification from gene-level analyses)
    if barcode_clusters is not None and not barcode_clusters.empty:
        adata = iso.calculate_barcode_match_percentage(adata, barcode_clusters)
    else:
        if sample_id is not None:
            adata.obs['cluster'] = sample_id
        else:
            adata.obs['cluster'] = gff_source

    
    # Need to filter cell barcodes to those with a sufficient number of expressed genes
    adata.obs['total_counts'] = np.sum(adata.X, axis=1)
    adata = adata[adata.obs['total_counts'] >= 100, :]
    remaining_cells = adata.n_obs
    print(f"Number of remaining cells with >=100 reads: {remaining_cells}")

    isoform_names = adata.var.index.tolist()
    print('isoform features |',isoform_names[:5])

    def parse_isoform_mapping(isoform_association_path, gff_source = None):
        """ Parse the isoform mapping file and return dictionaries of isoform-to-collapsed-isoforms """
        isoform_to_meta_isoform = collections.OrderedDict()

        print ('Mapping isoforms to meta-isoforms...',)
        unique_isoforms={}
        unique_ref_isoforms={}
        isoform_to_gene = {}
        with open(isoform_association_path, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                #ENSG00000241860	PB.10.25	ND167_MPP_3k	ENST00000655252.1	gencode.v45.annotation
                gene, ref_isoform, ref_gff, source_isoform, source_gff = row
                if source_gff == gff_source:
                    isoform_to_meta_isoform[source_isoform] = gene+':'+ref_isoform
                    unique_isoforms[source_isoform,source_gff]=[]
                    unique_ref_isoforms[ref_isoform,ref_gff]=[]
                if ref_gff == gff_source:
                    isoform_to_meta_isoform[ref_isoform] = gene+':'+ref_isoform
                    unique_isoforms[ref_isoform,ref_gff]=[]
                    unique_ref_isoforms[ref_isoform,ref_gff]=[]
                isoform_to_gene[ref_isoform] = gene

        print('unique isoforms |',len(unique_ref_isoforms),'out of',len(unique_isoforms))
        return isoform_to_meta_isoform, isoform_to_gene

    # Load isoform mapping information
    isoform_to_meta_isoform, isoform_to_gene = parse_isoform_mapping(isoform_association_path, gff_source=gff_source)

    isoform_counts = {}
    # Iterate through each cell barcode and isoform
    print ('Computing junctions counts for cells')
    cell_data = adata.X.toarray()
    isoform_ids = adata.var_names

    for isoform_idx, isoform in tqdm(enumerate(isoform_ids), total=len(isoform_ids), desc="Processing isoforms"):
        counts = cell_data[:, isoform_idx]
        if np.any(counts > 0):
            if isoform in isoform_to_meta_isoform: # Exclude isoforms without gene annotations
                ref_isoform = isoform_to_meta_isoform[isoform]
                for cell_index, count in zip(np.where(counts > 0)[0], counts[counts > 0]):
                    if ref_isoform not in isoform_counts:
                        isoform_counts[ref_isoform] = np.zeros(len(adata.obs_names), dtype=int)
                    isoform_counts[ref_isoform][cell_index] += count

    # Convert the counts to a sparse matrix (and report progress)
    isoform_counts_df = pd.DataFrame({k: pd.Series(v, index=adata.obs_names) for k, v in tqdm(isoform_counts.items(), desc="Creating DataFrame")})
    sparse_junction_matrix = lil_matrix(isoform_counts_df.shape)
    for i, (index, row) in enumerate(tqdm(isoform_counts_df.iterrows(), total=isoform_counts_df.shape[0], desc="Converting to sparse matrix")):
        sparse_junction_matrix[i, :] = row.values
    sparse_junction_matrix = sparse_junction_matrix.tocsr()

    # Create a new AnnData object for the junction counts
    isoform_adata = ad.AnnData(X=sparse_junction_matrix, obs=adata.obs, var=pd.DataFrame(index=isoform_counts.keys()))

    # Add gene annotations
    isoform_adata.var['gene'] = [isoform_to_gene.get(isoform, '') for isoform in isoform_adata.var_names]
    
    # Some weird key's get added when excluding barcode to cluster associations
    isoform_adata.obs = isoform_adata.obs.dropna().reset_index(drop=True)


    # Write the junction counts to an h5ad file
    h5ad_dir = os.path.basename(gff_source)
    
    #isoform_adata.write_h5ad("filtered_junction.h5ad")
    isoform_adata.write_h5ad(f"{h5ad_dir.split('.g')[0]}-isoform.h5ad", compression='gzip')
    print ('h5ad exported')

    export_pseudobulk = False
    if export_pseudobulk:
        # Compute pseudo-cluster counts and write them to a file
        grouped = isoform_counts_df.groupby(adata.obs['cluster'])
        
        # Remove clusters with less than 10 cells
        filtered_summed_groups = grouped.filter(lambda x: len(x) >= 10).groupby(adata.obs['cluster']).sum()
        
        # Remove junctions with less than 10 reads total
        filtered_summed_groups = filtered_summed_groups.loc[:, filtered_summed_groups.sum(axis=0) >= 10]
        filtered_summed_groups_transposed = filtered_summed_groups.transpose()
        filtered_summed_groups_transposed.to_csv(f"{h5ad_dir.split('.g')[0]}_isocounts.txt", sep='\t')
        print ('pseudobulk cluster junction counts exported')

    return isoform_adata

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
    exportConsensusIsoformMatrix(args.matrix_path,isoform_association_path,gff_source=args.gff_source,barcode_clusters=barcode_clusters,rev=args.reverse_complement)

