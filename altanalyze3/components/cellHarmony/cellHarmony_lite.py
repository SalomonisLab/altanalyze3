import os, time, sys
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
from glob import glob
from tqdm import tqdm
from scipy.spatial.distance import cdist
import warnings
warnings.filterwarnings("ignore", message="Variable names are not unique. To make them unique, call `.var_names_make_unique`.")

# Main function
def combine_and_align_h5(
    h5_files, 
    cellharmony_ref,
    h5ad_file=None,
    output_dir="output",
    export_cptt=False,
    export_h5ad=False,
    min_genes=500,
    min_cells=5,
    min_counts=1000,
    mit_percent=10,
):

    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)

    # Load reference centroid gene expression matrix
    reference_df = pd.read_csv(cellharmony_ref, sep='\t', index_col=0)
    cell_populations = reference_df.columns.tolist()

    # If supplied with a precomputed h5ad - proceed to filtering and alignment
    if h5ad_file!=None:
        adata_combined = sc.read_h5ad(h5ad_file)
        print(f"reimported adata shape: {adata_combined.shape} (cells x genes)")

    else:
        # Load and combine all 10x .h5 files with sample/group metadata
        adata_list = []
        for h5_file in tqdm(h5_files, desc="Loading input files"):
            sample_name = os.path.basename(h5_file).replace(".h5", "").replace(".mtx", "").replace(".gz", "")   
            if '__' in sample_name:
                group_name = sample_name.split('__')[0]  # If groups are separated by __
            else:
                if '_' in sample_name:
                    group_name = '_'.join(sample_name.split('_')[:-1])
                else:
                    group_name = sample_name
            if h5_file.endswith(".h5"):
                adata = sc.read_10x_h5(h5_file)
            else:
                mtx_dir = os.path.dirname(h5_file)
                adata = sc.read_10x_mtx(mtx_dir, var_names='gene_symbols', cache=True)
            adata.var_names_make_unique()
            adata.obs_names = [f"{bc}.{sample_name}" for bc in adata.obs_names]
            adata.obs["sample"] = sample_name
            adata.obs["group"] = group_name
            adata_list.append(adata)

        adata_combined = ad.concat(adata_list, label="sample", join="outer", fill_value=0)
        print(f"adata shape: {adata_combined.shape} (cells x genes)")

        # Scanpy QC filtering
        sc.pp.filter_cells(adata_combined, min_genes=min_genes)
        print(f"Cells remaining after min_genes {min_genes} filtering: {adata_combined.n_obs}")
        sc.pp.filter_genes(adata_combined, min_cells=min_cells)
        sc.pp.filter_cells(adata_combined, min_counts=min_counts)
        print(f"Cells remaining after min_counts {min_counts} filtering: {adata_combined.n_obs}")

        mito_genes = adata_combined.var_names.str.upper().str.startswith("Mt-")
        adata_combined.obs["pct_counts_mt"] = (
            np.sum(adata_combined[:, mito_genes].X, axis=1).A1 /
            np.sum(adata_combined.X, axis=1).A1
        ) * 100

        # Filter cells with >10% mitochondrial gene expression
        adata_combined = adata_combined[adata_combined.obs["pct_counts_mt"] < mit_percent].copy()
        print(f"Cells remaining after gene-count filtering: {adata_combined.n_obs}")

    # Log2 normalization to counts per ten thousand (CPTT)
    with tqdm(total=3, desc="Normalization steps") as pbar:
        sc.pp.normalize_total(adata_combined, target_sum=1e4); pbar.update(1)
        sc.pp.log1p(adata_combined); pbar.update(1)
        adata_combined.X = adata_combined.X / np.log(2); pbar.update(1)  # Convert natural log to log2
        
    if export_h5ad:
        adata_combined.write(output_dir+"/combined_qc_normalized.h5ad", compression="gzip")

    # Filter to marker genes
    marker_genes = reference_df.index
    genes_present = reference_df.index.intersection(adata_combined.var_names)
    adata_filtered = adata_combined[:, genes_present].copy()

    missing_genes = set(reference_df.index) - set(genes_present)
    if missing_genes:
        print(f"Warning: {len(missing_genes)} marker genes not found in dataset and will be excluded out of {len(marker_genes)}.")
        if len(missing_genes) <= 10:
            print(f"Missing genes: {sorted(missing_genes)}")
        else:
            print(f"(First 10) Missing genes: {sorted(list(missing_genes))[:10]}")

    if export_cptt:
        # Step 7: Export CPTT dense matrix (gene x cells)
        cptt_df = pd.DataFrame(
            adata_filtered.X,
            index=adata_filtered.obs_names,
            columns=adata_filtered.var_names
        ).T
        cptt_df.insert(0, "UID", cptt_df.index)
        cptt_df.to_csv(output_dir+"/CPTT_matrix.txt", sep="\t", index=False)

    # Match each cell to the best reference profile using cosine similarity
    print ('Aligning cells to reference...')
    align_start_time = time.time()
    ref_matrix = reference_df.loc[genes_present].T
    ref_norm = ref_matrix.div(np.linalg.norm(ref_matrix, axis=1), axis=0)

    query_matrix = pd.DataFrame(
        adata_filtered.X.toarray() if sp.issparse(adata_filtered.X) else adata_filtered.X,
        index=adata_filtered.obs_names,
        columns=adata_filtered.var_names
    )

    query_norm = query_matrix.div(np.linalg.norm(query_matrix, axis=1), axis=0)

    # Compute cosine similarities (1 - cosine distance)
    similarities = 1 - cdist(query_norm.values, ref_norm.values, metric="cosine")

    # Get assignments and alignment scores
    best_matches = np.argmax(similarities, axis=1)
    assignments = ref_norm.index[best_matches]
    alignment_scores = similarities[np.arange(len(best_matches)), best_matches]  # Alignment score per cell

    # Create output DataFrame
    cell_populations
    match_df = pd.DataFrame({
        "CellBarcode": query_matrix.index,
        "AssignedCellState": assignments,
        "AlignmentScore": alignment_scores
    })

    ordered_match_df = pd.concat([
        match_df[match_df["AssignedCellState"] == pop].sort_values("AlignmentScore", ascending=False)
        for pop in reference_df.columns if pop in match_df["AssignedCellState"].values
    ], ignore_index=True)

    ordered_match_df.to_csv(output_dir + "/cellHarmony_lite_assignments.txt", sep="\t", index=False)
    elapsed_time = time.time() - align_start_time
    print(f"Cosine similarity computation completed in {elapsed_time:.2f} seconds.")
    
    assignment_summary = match_df["AssignedCellState"].value_counts()
    print("Assignment summary (cells per reference state):")
    print(assignment_summary.to_string())

    elapsed_time = time.time() - start_time
    print(f"Analysis completed in {elapsed_time:.2f} seconds.")
    return ordered_match_df

# Usage example
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Compute cellHarmony from h5 files')
    parser.add_argument('--h5dir', type=str, default=None, help='the path to folder of h5 files')
    parser.add_argument('--refdir', type=str, default=None, help='cellHarmony reference states and genes')
    parser.add_argument('--outdir', type=str, default='output', help='output dir for the output file')

    # optional args
    parser.add_argument('--h5ad', type=str, default=None, help='existing h5ad')
    parser.add_argument('--cptt', action='store_true', help='export a dense tsv for ref gene normalized exp')
    parser.add_argument('--export_h5ad', action='store_true', help='export an h5ad with all counts and normalized exp')
    parser.add_argument('--min_genes', type=int, default=500, help='min_genes for scanpy QC')
    parser.add_argument('--min_cells', type=int, default=0, help='min_cells for scanpy QC')
    parser.add_argument('--min_counts', type=int, default=1000, help='min_counts for scanpy QC')
    parser.add_argument('--mit_percent', type=int, default=10, help='mit_percent for scanpy QC')
    args = parser.parse_args()

    h5_directory = args.h5dir
    cellharmony_ref = args.refdir
    output_dir = args.outdir
    h5ad_file=args.h5ad
    export_cptt=args.cptt
    export_h5ad=args.export_h5ad
    min_genes=args.min_genes
    min_cells=args.min_cells
    min_counts=args.min_counts
    mit_percent=args.mit_percent

    if '.h5' not in h5_directory:
        h5_files = glob(os.path.join(h5_directory, "*.h5"))
    else:
        h5_files = [h5_directory]

    combine_and_align_h5(
        h5_files=h5_files,
        h5ad_file=h5ad_file,
        cellharmony_ref=cellharmony_ref,
        output_dir=output_dir,
        export_cptt=export_cptt,
        export_h5ad=export_h5ad,
        min_genes=min_genes,
        min_cells=min_cells,
        min_counts=min_counts,
        mit_percent=mit_percent,
    )
