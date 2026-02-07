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
import gc
import ctypes
import array
from pathlib import Path
import subprocess
import resource

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from . import isoform_matrix as iso
from . import gff_process as gff_process

_MEMTRACE_TAG = "MEMTRACE_TEMP"  # MEMTRACE_TEMP
_MEMTRACE_ENABLED = os.environ.get("ISOFORM_JUNC_MEMTRACE") == "1"


def _get_rss_stats():
    rss_mb = None
    rss_max_mb = None
    try:
        rss_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform == 'darwin':
            rss_max_mb = rss_kb / (1024 * 1024)
        else:
            rss_max_mb = rss_kb / 1024
    except Exception:
        pass
    try:
        output = subprocess.check_output(
            ["ps", "-o", "rss=", "-p", str(os.getpid())],
            text=True,
        ).strip()
        if output:
            rss_mb = float(output) / 1024
    except Exception:
        pass
    return rss_mb, rss_max_mb


def _memtrace(label):
    if not _MEMTRACE_ENABLED:
        return
    rss_mb, rss_max_mb = _get_rss_stats()
    parts = []
    if rss_mb is not None:
        parts.append(f"rss_mb={rss_mb:.1f}")
    if rss_max_mb is not None:
        parts.append(f"rss_max_mb={rss_max_mb:.1f}")
    status = " ".join(parts) if parts else "rss_mb=n/a"
    print(f"[{_MEMTRACE_TAG}] {label} {status}")


def _trim_memory():
    gc.collect()
    if sys.platform.startswith('linux'):
        try:
            libc = ctypes.CDLL("libc.so.6")
            libc.malloc_trim(0)
        except OSError:
            pass
        return
    if sys.platform == 'darwin':
        try:
            libmalloc = ctypes.CDLL("libmalloc.dylib")
            libmalloc.malloc_default_zone.restype = ctypes.c_void_p
            zone = libmalloc.malloc_default_zone()
            libmalloc.malloc_zone_pressure_relief.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
            libmalloc.malloc_zone_pressure_relief(zone, 0)
        except OSError:
            pass
        return
    if sys.platform.startswith('win'):
        try:
            process = ctypes.windll.kernel32.GetCurrentProcess()
            ctypes.windll.psapi.EmptyWorkingSet(process)
        except Exception:
            pass
    
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

def exportJunctionMatrix(matrix_dir, ensembl_exon_dir, gff_source, barcode_clusters=None, rev=False,
                         return_adata=True):
    """ 
    Decompose isoforms into exon-exon and exon-junctions counts at the single-cell and pseudobulk level
    """

    _memtrace("start")
    print(f"Reading GFF for junction-quantification: {gff_source}")
    
    # Annotate isoforms relative to ensembl exons
    transcript_associations = gff_process.consolidateLongReadGFFs(gff_source, ensembl_exon_dir, mode="collapse")
    _memtrace("after_consolidate_gff")

    # Load the prior computed exon annotations and coordinates for junction mapping
    exon_dict, gene_dict = parse_exon_file(ensembl_exon_dir)

    # Load the 10x matrix keyed by isoforms in the first column and cell barcode cluster annotations
    adata = iso.matrix_dir_to_adata(int_folder=matrix_dir, gene_is_index=True, feature='genes.tsv', feature_col=0, barcode='barcodes.tsv', barcode_col=0, matrix='matrix.mtx', rev=rev)
    
    # If existing cell clusters are provided already (e.g., supervised classification from gene-level analyses)
    if barcode_clusters is not None and not barcode_clusters.empty:
        adata = iso.calculate_barcode_match_percentage(adata, barcode_clusters)

    # Need to filter cell barcodes to those with a sufficient number of expressed genes
    counts_threshold = 0 ### Previously set to 100 but this cuts the pseudobulk counts in half
    adata.obs['total_counts'] = np.sum(adata.X, axis=1)
    try:
        adata = adata[adata.obs['total_counts'] >= counts_threshold, :]
    except Exception as e:
        print(f"An error occurred: {e}")
        print(f"AnnData: {anndata.__version__}")
        print ("WARNING -- SOFTWARE EXIT - likely outdated version of anndata (recommend version 0.10.9).");sys.exit()

    remaining_cells = adata.n_obs
    print(f"Number of remaining cells with >={counts_threshold} reads: {remaining_cells}")

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
    del transcript_associations

    junction_index = {}
    junction_list = []
    triplet_count = 0
    buffer = array.array('i')
    buffer_triplets = 100000
    output_prefix = gff_source.split('.g')[0]
    output_path = Path(f"{output_prefix}-junction.h5ad")
    tmp_triplets_path = output_path.parent / f".{output_path.name}.triplets.bin"
    memmap_dir = output_path.parent / f".{output_path.stem}.memmap"
    # Iterate through each cell barcode and isoform
    print ('Computing junctions counts for cells')

    isoform_ids = adata.var_names
    matrix = adata.X
    if not hasattr(matrix, "tocsc"):
        from scipy import sparse
        matrix = sparse.csr_matrix(matrix)
    matrix = matrix.tocsc(copy=False)

    with open(tmp_triplets_path, 'wb') as tmp_handle:
        for isoform_idx, isoform in tqdm(enumerate(isoform_ids), total=len(isoform_ids), desc="Processing isoforms"):
            junctions = isoform_to_junctions.get(isoform, [])
            if not junctions:
                continue
            start = matrix.indptr[isoform_idx]
            end = matrix.indptr[isoform_idx + 1]
            if start == end:
                continue
            row_indices = matrix.indices[start:end]
            values = matrix.data[start:end]
            junction_cols = []
            for junction in junctions:
                col_index = junction_index.get(junction)
                if col_index is None:
                    col_index = len(junction_list)
                    junction_index[junction] = col_index
                    junction_list.append(junction)
                junction_cols.append(col_index)
            for row_index, count in zip(row_indices, values):
                if count == 0:
                    continue
                count_value = int(count)
                for col_index in junction_cols:
                    buffer.extend([int(row_index), int(col_index), count_value])
                    triplet_count += 1
                if len(buffer) >= buffer_triplets * 3:
                    buffer.tofile(tmp_handle)
                    buffer = array.array('i')
        if buffer:
            buffer.tofile(tmp_handle)
            buffer = array.array('i')
    del matrix

    from scipy.sparse import coo_matrix
    if triplet_count:
        memmap_dir.mkdir(parents=True, exist_ok=True)
        rows_path = memmap_dir / "rows.dat"
        cols_path = memmap_dir / "cols.dat"
        data_path = memmap_dir / "data.dat"
        rows = np.memmap(rows_path, dtype=np.int32, mode='w+', shape=(triplet_count,))
        cols = np.memmap(cols_path, dtype=np.int32, mode='w+', shape=(triplet_count,))
        data = np.memmap(data_path, dtype=np.int32, mode='w+', shape=(triplet_count,))
        idx = 0
        chunk_triplets = 500000
        with open(tmp_triplets_path, 'rb') as tmp_handle:
            while True:
                chunk = np.fromfile(tmp_handle, dtype=np.int32, count=chunk_triplets * 3)
                if chunk.size == 0:
                    break
                triplets = chunk.reshape(-1, 3)
                end = idx + triplets.shape[0]
                rows[idx:end] = triplets[:, 0]
                cols[idx:end] = triplets[:, 1]
                data[idx:end] = triplets[:, 2]
                idx = end
        sparse_junction_matrix = coo_matrix(
            (data, (rows, cols)),
            shape=(len(adata.obs_names), len(junction_list))
        ).tocsr()
        try:
            rows.flush()
            cols.flush()
            data.flush()
        except AttributeError:
            pass
        for path in (rows_path, cols_path, data_path):
            try:
                path.unlink()
            except FileNotFoundError:
                continue
        try:
            memmap_dir.rmdir()
        except OSError:
            pass
    else:
        sparse_junction_matrix = coo_matrix(
            (len(adata.obs_names), len(junction_list))
        ).tocsr()
    try:
        tmp_triplets_path.unlink()
    except FileNotFoundError:
        pass
    
    # Create a new AnnData object for the junction counts
    #print("Original adata.obs_names:", adata.obs_names[:5].tolist())
    print(f"Junctions retained for h5ad: {len(junction_list)}")
    junction_adata = ad.AnnData(X=sparse_junction_matrix, obs=adata.obs, var=pd.DataFrame(index=junction_list))
    #junction_adata = ad.AnnData(X=sparse_junction_matrix, obs=adata.obs.copy(), var=pd.DataFrame(index=junction_counts.keys()))
    junction_adata.obs_names = adata.obs_names
    #print("Final junction_adata.obs_names:", junction_adata.obs_names[:5].tolist())
    del sparse_junction_matrix
    del adata
    del isoform_to_junctions
    del junction_index
    del junction_list
    _memtrace("after_cleanup_pre_h5ad")

    # Add gene annotations
    junction_adata.var['gene'] = [isoform_to_gene.get(isoform, '') for isoform in junction_adata.var_names]
    del isoform_to_gene

    # Write the junction counts to an h5ad file
    #h5ad_dir = os.path.basename(gff_source)
    #junction_adata.write_h5ad("filtered_junction.h5ad")
    junction_adata.write_h5ad(output_path, compression='gzip')
    print(f"Saved junction h5ad: {output_path}")
    _memtrace("after_write_h5ad")
    #loaded_adata = ad.read_h5ad(f"{h5ad_dir.split('.g')[0]}.h5ad")
    #print("Loaded barcodes after saving:", loaded_adata.obs_names[:5].tolist())

    print ('h5ad exported')
    
    export_pseudobulk = False
    if export_pseudobulk:
        if barcode_clusters is not None and not barcode_clusters.empty:
            junction_counts_df = pd.DataFrame.sparse.from_spmatrix(
                sparse_junction_matrix, index=adata.obs_names, columns=junction_list
            )
            # Compute pseudo-cluster counts and write them to a file
            grouped = junction_counts_df.groupby(adata.obs['cluster'])
            
            # Remove clusters with less than 10 cells
            filtered_summed_groups = grouped.filter(lambda x: len(x) >= 10).groupby(adata.obs['cluster']).sum()
            
            # Remove junctions with less than 10 reads total
            filtered_summed_groups = filtered_summed_groups.loc[:, filtered_summed_groups.sum(axis=0) >= 10]
            filtered_summed_groups_transposed = filtered_summed_groups.transpose()
            filtered_summed_groups_transposed.to_csv(f"{h5ad_dir.split('.g')[0]}_counts.txt", sep='\t')
            print ('pseudobulk cluster junction counts exported')

    gff_process.clearEnsemblCache()
    _memtrace("after_clear_gff_cache")

    _trim_memory()
    _memtrace("after_trim")
    if not return_adata:
        # RETURN_ADATA_TEMP
        del junction_adata
        _trim_memory()
        _memtrace("after_drop_return")
        return None
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
 
