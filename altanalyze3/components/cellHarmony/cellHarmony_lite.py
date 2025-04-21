import os, time, sys
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
from glob import glob
from tqdm import tqdm
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore", message="Variable names are not unique. To make them unique, call `.var_names_make_unique`.")
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*invalid value encountered in log2.*")

plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['figure.facecolor'] = 'white'

def save_marker_genes(adata, groupby, output_file):
    deg = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    deg.to_csv(output_file, sep='\t', index=False)

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
    generate_umap=False,
    save_adata=False,
    unsupervised_cluster=False
):
    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)

    reference_df = pd.read_csv(cellharmony_ref, sep='\t', index_col=0)
    cell_populations = reference_df.columns.tolist()

    if h5ad_file is not None:
        adata_combined = sc.read_h5ad(h5ad_file)
        print(f"reimported adata shape: {adata_combined.shape} (cells x genes)")
    else:
        adata_list = []
        for h5_file in tqdm(h5_files, desc="Loading input files"):
            sample_name = os.path.basename(h5_file).replace(".h5", "").replace(".mtx", "").replace(".gz", "")   
            group_name = sample_name.split('__')[0] if '__' in sample_name else '_'.join(sample_name.split('_')[:-1])
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

        sc.pp.filter_cells(adata_combined, min_genes=min_genes)
        print(f"Cells remaining after min_genes {min_genes} filtering: {adata_combined.n_obs}")
        sc.pp.filter_genes(adata_combined, min_cells=min_cells)
        sc.pp.filter_cells(adata_combined, min_counts=min_counts)
        print(f"Cells remaining after min_counts {min_counts} filtering: {adata_combined.n_obs}")

        mito_genes = adata_combined.var_names.str.upper().str.startswith("MT-")
        adata_combined.obs["pct_counts_mt"] = (
            np.sum(adata_combined[:, mito_genes].X, axis=1).A1 /
            np.sum(adata_combined.X, axis=1).A1
        ) * 100

        adata_combined = adata_combined[adata_combined.obs["pct_counts_mt"] < mit_percent].copy()
        print(f"Cells remaining after gene-count filtering: {adata_combined.n_obs}")

    with tqdm(total=3, desc="Normalization steps") as pbar:
        sc.pp.normalize_total(adata_combined, target_sum=1e4); pbar.update(1)
        sc.pp.log1p(adata_combined); pbar.update(1)
        #adata_combined.X = adata_combined.X / np.log(2)
        pbar.update(1)

    if export_h5ad:
        adata_combined.write(output_dir+"/combined_qc_normalized.h5ad", compression="gzip")

    marker_genes = reference_df.index
    genes_present = reference_df.index.intersection(adata_combined.var_names)
    adata_filtered = adata_combined[:, genes_present].copy()

    missing_genes = set(reference_df.index) - set(genes_present)
    if missing_genes:
        print(f"Warning: {len(missing_genes)} marker genes not found in dataset and will be excluded out of {len(marker_genes)}.")
        print(f"(First 10) Missing genes: {sorted(list(missing_genes))[:10]}")

    if export_cptt:
        cptt_df = pd.DataFrame(
            adata_filtered.X,
            index=adata_filtered.obs_names,
            columns=adata_filtered.var_names
        ).T
        cptt_df.insert(0, "UID", cptt_df.index)
        cptt_df.to_csv(output_dir+"/CPTT_matrix.txt", sep="\t", index=False)

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
    similarities = 1 - cdist(query_norm.values, ref_norm.values, metric="cosine")
    best_matches = np.argmax(similarities, axis=1)
    assignments = ref_norm.index[best_matches]
    alignment_scores = similarities[np.arange(len(best_matches)), best_matches]

    match_df = pd.DataFrame({
        "CellBarcode": query_matrix.index,
        "AssignedCellState": assignments,
        "AlignmentScore": alignment_scores
    })

    ordered_match_df = pd.concat([
        match_df[match_df["AssignedCellState"] == pop].sort_values("AlignmentScore", ascending=False)
        for pop in reference_df.columns if pop in match_df["AssignedCellState"].values
    ], ignore_index=True)

    adata_combined = adata_combined[match_df.CellBarcode].copy()
    adata_combined.obs['AssignedCellState'] = match_df.set_index('CellBarcode').loc[adata_combined.obs_names]['AssignedCellState']

    ordered_match_df.to_csv(output_dir + "/cellHarmony_lite_assignments.txt", sep="\t", index=False)
    print(f"Cosine similarity computation completed in {time.time() - align_start_time:.2f} seconds.")
    print("Assignment summary (cells per reference state):")
    print(match_df["AssignedCellState"].value_counts().to_string())

    if generate_umap:
        os.chdir(output_dir)

        def downsample_cells_per_group(adata, groupby, max_cells=50):
            idx = []
            for group, count in adata.obs[groupby].value_counts().items():
                cells = adata.obs_names[adata.obs[groupby] == group]
                selected = np.random.choice(cells, min(len(cells), max_cells), replace=False)
                idx.extend(selected)
            return adata[idx].copy()

        adata_filtered = adata_combined[adata_combined.obs['AssignedCellState'].isin(
            adata_combined.obs['AssignedCellState'].value_counts()[lambda x: x >= 10].index
        )].copy()
        adata_filtered = downsample_cells_per_group(adata_filtered, 'AssignedCellState', max_cells=150)

        sc.tl.rank_genes_groups(adata_filtered, groupby='AssignedCellState', method='wilcoxon', use_raw=False)
        save_marker_genes(adata_filtered, 'AssignedCellState', os.path.join(output_dir, 'supervised_markers.txt'))

        deg_results = pd.DataFrame(adata_filtered.uns['rank_genes_groups']['names'])
        markers = list(set([gene for col in deg_results.columns for gene in deg_results[col][:10]]))
        excluded_prefixes = ('mt-', 'rp', 'xist')
        markers = [gene for gene in markers if not gene.lower().startswith(excluded_prefixes) and gene in adata_filtered.var_names]

        adata_markers = adata_filtered[:, markers].copy()
        sc.pp.pca(adata_markers, n_comps=50)
        sc.pp.neighbors(adata_markers)
        sc.tl.umap(adata_markers)

        # Reuse the selected markers on the full dataset (all cells)
        adata_all_cells = adata_combined[:, markers].copy()
        sc.pp.pca(adata_all_cells, n_comps=50)
        sc.pp.neighbors(adata_all_cells)
        sc.tl.umap(adata_all_cells)

        coords = adata_all_cells.obsm['X_umap']
        adata_combined.obs['UMAP-X'] = coords[:, 0]
        adata_combined.obs['UMAP-Y'] = coords[:, 1]
        adata_combined.obsm['X_umap'] = coords

        sc.pl.umap(adata_combined, color='AssignedCellState', 
            save=f"_UMAP.pdf", show=False, legend_loc='on data',
            legend_fontsize=3,legend_fontweight='normal')

        sc.pl.rank_genes_groups_heatmap(
            adata_filtered,
            #n_genes=5,
            show=False,
            save=f"_heatmap.pdf",
            standard_scale='var',
            dendrogram=False,
            swap_axes=True,
            var_group_rotation=90,
        )

    if unsupervised_cluster:
        os.chdir(output_dir)
        adata_unsup = adata_combined.copy()
        sc.pp.highly_variable_genes(adata_unsup, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata_unsup = adata_unsup[:, adata_unsup.var.highly_variable]
        sc.pp.scale(adata_unsup, max_value=10)
        sc.pp.pca(adata_unsup, n_comps=50)
        sc.pp.neighbors(adata_unsup)
        sc.tl.umap(adata_unsup)
        sc.tl.leiden(adata_unsup)

        sc.tl.rank_genes_groups(adata_unsup, groupby='leiden', method='wilcoxon', use_raw=False)
        save_marker_genes(adata_unsup, 'leiden', os.path.join(output_dir, 'unsupervised_markers.txt'))

        sc.pl.umap(adata_unsup, color='leiden', save="_unsupervised_umap.pdf", 
            show=False, legend_loc='on data', legend_fontsize=5,legend_fontweight='normal')

        sc.pl.rank_genes_groups_heatmap(
            adata_unsup,
            groupby='leiden',
            n_genes=5,
            show=False,
            save="_unsupervised_heatmap.pdf",
            standard_scale='var',
            dendrogram=False,
            swap_axes=True,
            var_group_rotation=0,
        )

    if save_adata:
        # If unsupervised clustering was run, merge leiden clusters into adata_combined before saving
        if unsupervised_cluster and 'leiden' in adata_unsup.obs.columns:
            adata_combined.obs['leiden'] = adata_unsup.obs['leiden']
        adata_combined.write(os.path.join(output_dir, "combined_with_umap_and_markers.h5ad"), compression="gzip")


    print(f"Analysis completed in {time.time() - start_time:.2f} seconds.")
    return ordered_match_df

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Compute cellHarmony from h5 files')
    parser.add_argument('--h5dir', type=str, default=None, help='the path to folder of h5 files')
    parser.add_argument('--refdir', type=str, default=None, help='cellHarmony reference states and genes')
    parser.add_argument('--outdir', type=str, default='output', help='output dir for the output file')
    parser.add_argument('--h5ad', type=str, default=None, help='existing h5ad')
    parser.add_argument('--cptt', action='store_true', help='export a dense tsv for ref gene normalized exp')
    parser.add_argument('--export_h5ad', action='store_true', help='export an h5ad with all counts and normalized exp')
    parser.add_argument('--min_genes', type=int, default=500, help='min_genes for scanpy QC')
    parser.add_argument('--min_cells', type=int, default=3, help='min_cells for scanpy QC')
    parser.add_argument('--min_counts', type=int, default=1000, help='min_counts for scanpy QC')
    parser.add_argument('--mit_percent', type=int, default=10, help='mit_percent for scanpy QC')
    parser.add_argument('--generate_umap', action='store_true', help='generate UMAP and marker analysis')
    parser.add_argument('--save_adata', action='store_true', help='save updated AnnData object')
    parser.add_argument('--unsupervised_cluster', action='store_true', help='perform unsupervised clustering analysis')
    args = parser.parse_args()

    h5_directory = args.h5dir
    cellharmony_ref = args.refdir
    output_dir = args.outdir
    h5ad_file = args.h5ad
    export_cptt = args.cptt
    export_h5ad = args.export_h5ad
    min_genes = args.min_genes
    min_cells = args.min_cells
    min_counts = args.min_counts
    mit_percent = args.mit_percent
    generate_umap = args.generate_umap
    save_adata = args.save_adata
    unsupervised_cluster = args.unsupervised_cluster

    h5_files = glob(os.path.join(h5_directory, "*.h5")) if '.h5' not in h5_directory else [h5_directory]

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
        generate_umap=generate_umap,
        save_adata=save_adata,
        unsupervised_cluster=unsupervised_cluster
    )
