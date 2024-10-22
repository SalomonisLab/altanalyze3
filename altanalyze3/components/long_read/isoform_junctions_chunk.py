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
import time

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from . import isoform_matrix as iso
from . import gff_process as gff_process
    
def parse_exon_file(ensembl_exon_dir):
    """ Import Ensembl exon genomic information """
    exon_dict = {}
    gene_dict = {}
    with open(ensembl_exon_dir, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            key = (row['gene'], row['exon-id'])
            if row['strand'] == '+':
                value = row['chromosome'], row['exon-region-start(s)'],row['exon-region-stop(s)']
            else:
                value = row['chromosome'], row['exon-region-stop(s)'],row['exon-region-start(s)']
            exon_dict[key] = value
            gene_dict[row['gene']] = row['chromosome']
    return exon_dict, gene_dict

def exportJunctionMatrix(matrix_dir, ensembl_exon_dir, gff_source, barcode_clusters=None, rev=False):
    """ 
    Decompose isoforms into exon-exon and exon-junctions counts at the single-cell and pseudobulk level
    """

    # Annotate isoforms relative to ensembl exons
    transcript_associations = gff_process.consolidateLongReadGFFs(gff_source, ensembl_exon_dir)

    # Load the prior computed exon annotations and coordinates for junction mapping
    exon_dict, gene_dict = parse_exon_file(ensembl_exon_dir)

    # Load the 10x matrix keyed by isoforms in the first column and cell barcode cluster annotations
    adata = iso.mtx_to_adata(int_folder=matrix_dir, gene_is_index=True, feature='genes.tsv', feature_col=0, barcode='barcodes.tsv', barcode_col=0, matrix='matrix.mtx', rev=rev)
    
    # If existing cell clusters are provided already (e.g., supervised classification from gene-level analyses)
    if barcode_clusters is not None and not barcode_clusters.empty:
        adata = iso.calculate_barcode_match_percentage(adata, barcode_clusters)

    # Need to filter cell barcodes to those with a sufficient number of expressed genes
    adata.obs['total_counts'] = np.sum(adata.X, axis=1)
    adata = adata[adata.obs['total_counts'] >= 100, :]
    remaining_cells = adata.n_obs
    print(f"Number of remaining cells with >=100 reads: {remaining_cells}")

    isoform_names = adata.var.index.tolist()
    print('isoform features |',isoform_names[:5])

    def parse_isoform_mapping(isoform_mapping_file):
        """ Parse the isoform mapping file and return dictionaries of isoform-to-junctions and isoform-to-genes """
        isoform_to_junctions = collections.OrderedDict()
        isoform_to_gene = {}

        print ('Mapping isoforms to junctions...',)
        with open(isoform_mapping_file, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                gene, strand, exon_structure, isoform, gff_name = row
                # If multiple gff files used - restrict to isoforms to the appropriate gff
                if 'UNK' not in gene: # (gff_source==gff_name or gff_source==None) and
                    exons = exon_structure.split('|')[1:] #[1:-1] - consider APA
                    junctions = []
                    for exon in exons:
                        geneID = gene
                        if ':' in exon:
                            geneID,exon = exon.split(':')
                        if '_' in exon: # Novel exon
                            exon_prefix,start = exon.split('_')
                            stop = start
                            chr = gene_dict[geneID]
                        else:
                            try: chr, start, stop = exon_dict[geneID,exon]
                            except:
                                print ('GFF Annotation Construction ERROR!!!',gene,geneID,exon,isoform)
                                print (exon_structure, exons);sys.exit() 
                        junctions.append((exon,(start, stop),geneID))

                    # Extract unique junction position and name per AltAnalyze ID conventions (exclude non-real junctions - exon start, stop)
                    # [('I10.1', ('16606', '15948')), ('E11.1', ('15947', '15796')), ('E12.1', ('15038', '15005')), ('E12.1_14970', ('14970', '14970')), ('I13.1_14829', ('14829', '14829'))] gives ['ENSG00000227232:I10.1-E11.1=chr1:15948-15947', 'ENSG00000227232:E11.1-E12.1=chr1:15796-15038', 'ENSG00000227232:E12.1_14970-I13.1_14829=chr1:14970-14829']
                    splice_sites = []
                    previous_exon = None
                    for i in range(len(junctions) - 1):
                        current_exon, current_range, current_gene = junctions[i]
                        next_exon, next_range, next_gene = junctions[i + 1]
                        base_current_exon = current_exon.split('_')[0]
                        base_next_exon = next_exon.split('_')[0]
                        #if base_current_exon == base_next_exon and '_' in next_exon:
                        #    continue
                        if current_gene == next_gene:
                            splice_site = f"{current_gene}:{current_exon}-{next_exon}={chr}:{current_range[1]}-{next_range[0]}"
                        else:
                            splice_site = f"{current_gene}:{current_exon}-{next_gene}:{next_exon}={chr}:{current_range[1]}-{next_range[0]}"
                        splice_sites.append(splice_site)
                    """
                    if 'PB.145.87' == isoform:
                        print(junctions)
                        print (exon_structure, exons, splice_sites);sys.exit() 
                    """
                    isoform_to_junctions[isoform] = splice_sites
                    isoform_to_gene[isoform] = gene
        return isoform_to_junctions, isoform_to_gene

    # Load isoform mapping information
    isoform_to_junctions, isoform_to_gene = parse_isoform_mapping(transcript_associations)

    junction_counts = {}
    # Iterate through each cell barcode and isoform
    print ('Computing junctions counts for cells')

    isoform_ids = adata.var_names
    # Define the chunk size
    chunk_size = 2000  # Adjust this number based on memory usage and processing speed

    # Iterate over the isoforms in chunks
    for start in tqdm(range(0, len(isoform_ids), chunk_size), desc="Processing isoforms in chunks"):
        end = min(start + chunk_size, len(isoform_ids))
        chunk_isoform_ids = isoform_ids[start:end]
        
        # Convert the current chunk of isoforms to dense format
        chunk_data = adata.X[:, start:end].toarray()  # Only this chunk of isoforms is converted to dense
        
        # Process isoforms in this chunk
        for local_idx, isoform in enumerate(chunk_isoform_ids):
            counts = chunk_data[:, local_idx]
            if np.any(counts > 0):
                junctions = isoform_to_junctions.get(isoform, [])
                for cell_index, count in zip(np.where(counts > 0)[0], counts[counts > 0]):
                    for junction in junctions:
                        if junction not in junction_counts:
                            junction_counts[junction] = np.zeros(len(adata.obs_names), dtype=int)
                        junction_counts[junction][cell_index] += count


    # Convert the counts to a dense DataFrame and then to a sparse matrix
    junction_counts_df = pd.DataFrame({k: pd.Series(v, index=adata.obs_names) for k, v in tqdm(junction_counts.items(), desc="Creating DataFrame")})
    def lil_matrix_method():
        from scipy.sparse import lil_matrix
        junction_counts_df = pd.DataFrame({k: pd.Series(v, index=adata.obs_names) for k, v in tqdm(junction_counts.items(), desc="Creating DataFrame")})
        sparse_junction_matrix = lil_matrix(junction_counts_df.shape)
        for i, (index, row) in enumerate(tqdm(junction_counts_df.iterrows(), total=junction_counts_df.shape[0], desc="Converting to sparse matrix")):
            sparse_junction_matrix[i, :] = row.values
        return sparse_junction_matrix.tocsr()

    def coo_matrix_method():
        from scipy.sparse import coo_matrix
        rows, cols, data = [], [], []
        for col_index, (k, values) in enumerate(tqdm(junction_counts.items(), desc="Building sparse matrix (COO)")):
            for row_index, count in enumerate(values):
                if count != 0:
                    rows.append(row_index)
                    cols.append(col_index)
                    data.append(count)
        return coo_matrix((data, (rows, cols)), shape=(len(adata.obs_names), len(junction_counts))).tocsr()

    def dok_matrix_method():
        from scipy.sparse import dok_matrix
        sparse_junction_matrix = dok_matrix((len(adata.obs_names), len(junction_counts)))
        for col_index, (k, values) in enumerate(tqdm(junction_counts.items(), desc="Building sparse matrix (DOK)")):
            for row_index, count in enumerate(values):
                if count != 0:
                    sparse_junction_matrix[row_index, col_index] = count
        return sparse_junction_matrix.tocsr()

    sparse_junction_matrix = coo_matrix_method()
    
    # Create a new AnnData object for the junction counts
    junction_adata = ad.AnnData(X=sparse_junction_matrix, obs=adata.obs, var=pd.DataFrame(index=junction_counts.keys()))

    # Add gene annotations
    junction_adata.var['gene'] = [isoform_to_gene.get(isoform, '') for isoform in junction_adata.var_names]

    # Write the junction counts to an h5ad file
    h5ad_dir = os.path.basename(gff_source)
    #junction_adata.write_h5ad("filtered_junction.h5ad")
    junction_adata.write_h5ad(f"{h5ad_dir.split('.g')[0]}.h5ad", compression='gzip')
    print ('h5ad exported')
    
    export_pseudobulk = False
    if export_pseudobulk:
        if barcode_clusters is not None and not barcode_clusters.empty:
            # Compute pseudo-cluster counts and write them to a file
            grouped = junction_counts_df.groupby(adata.obs['cluster'])
            
            # Remove clusters with less than 10 cells
            filtered_summed_groups = grouped.filter(lambda x: len(x) >= 10).groupby(adata.obs['cluster']).sum()
            
            # Remove junctions with less than 10 reads total
            filtered_summed_groups = filtered_summed_groups.loc[:, filtered_summed_groups.sum(axis=0) >= 10]
            filtered_summed_groups_transposed = filtered_summed_groups.transpose()
            filtered_summed_groups_transposed.to_csv(f"{h5ad_dir.split('.g')[0]}_counts.txt", sep='\t')
            print ('pseudobulk cluster junction counts exported')

    return junction_adata

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
    exportJunctionMatrix(args.matrix_path, args.ensembl_exon_dir, args.gff_source, barcode_clusters=barcode_clusters, rev=args.reverse_complement)

    """
    Example paths:
    gff = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/HSC_3k/ND167_HSC-3k.gtf'
    matrix_path = '/Users/saljh8/Dropbox/Revio/Apharesis-Sequel2/ND167_Iso-Seq/PacBio/HSC_3k/filtered_matrices'
    barcode_cluster = 'HSC.txt'
    ensembl_exon_dir = '/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl_exon.txt'
    """
 