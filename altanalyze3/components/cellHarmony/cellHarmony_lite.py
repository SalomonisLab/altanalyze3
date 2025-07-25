import os, time, sys, shutil, re
import tempfile
import pandas as pd
import numpy as np
import scipy.sparse as sp
from scipy.sparse import issparse
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
    unsupervised_cluster=False,
    append_obs_field=None,
    alignment_mode="classic"
):
    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)

    reference_df = pd.read_csv(cellharmony_ref, sep='\t', index_col=0)
    cell_populations = reference_df.columns.tolist()
    ref_name = os.path.basename(cellharmony_ref)[:-4]

    if h5ad_file is not None:
        adata_combined = sc.read_h5ad(h5ad_file)
        print(f"reimported adata shape: {adata_combined.shape} (cells x genes)")
    else:
        adata_list = []
        for path in tqdm(h5_files, desc="Loading input files"):
            if path.endswith(".h5"):
                sample_name = os.path.basename(path).replace(".h5", "")
                adata = sc.read_10x_h5(path)

            elif path.endswith(".h5ad"):
                adata = sc.read_h5ad(path)
                sample_name = os.path.basename(path).replace(".h5ad", "")
                group_name = sample_name.split('__')[0] if '__' in sample_name else '_'.join(sample_name.split('_')[:-1])
                if append_obs_field is not None:
                    if append_obs_field not in adata.obs.columns:
                        raise ValueError(f"Field '{append_obs_field}' not found in adata.obs columns for {sample_name}.")
                    appended_values = adata.obs[append_obs_field].astype(str).tolist()
                    adata.obs_names = [f"{bc}.{val}" for bc, val in zip(adata.obs_names, appended_values)]

            elif path.endswith((
                "_filtered_matrix.mtx.gz", "_counts.mtx.gz", "_matrix.mtx.gz",
                "filtered_matrix.mtx.gz", "matrix.mtx.gz", "matrix.mtx"
            )):
                # Determine suffix and prefix
                suffixes = [
                    "_filtered_matrix.mtx.gz", "_counts.mtx.gz", "_matrix.mtx.gz",
                    "filtered_matrix.mtx.gz", "matrix.mtx.gz", "matrix.mtx"
                ]
                suffix = next(s for s in suffixes if path.endswith(s))
                prefix = path[:-len(suffix)]

                # Determine file names
                if path.endswith(".gz"):
                    barcodes = prefix + ("_barcodes.tsv.gz" if any(path.endswith(f"_{s}") for s in ["filtered_matrix.mtx.gz", "counts.mtx.gz", "matrix.mtx.gz"]) else "barcodes.tsv.gz")
                    feature_options = [prefix + "_features.tsv.gz", prefix + "_genes.tsv.gz",
                                    prefix + "features.tsv.gz", prefix + "genes.tsv.gz"]
                else:
                    barcodes = prefix + ("_barcodes.tsv" if any(path.endswith(f"_{s}") for s in ["matrix.mtx"]) else "barcodes.tsv")
                    feature_options = [prefix + "_features.tsv", prefix + "_genes.tsv",
                                    prefix + "features.tsv", prefix + "genes.tsv"]

                # Select existing feature file
                features = next((f for f in feature_options if os.path.exists(f)), None)
                if not features:
                    raise FileNotFoundError(f"Missing features or genes file for prefix: {prefix}")
                if not os.path.exists(barcodes):
                    raise FileNotFoundError(f"Missing barcodes file for prefix: {prefix}")

                # Setup tmpdir and copy files
                tmpdir = tempfile.TemporaryDirectory()
                tmp_path = tmpdir.name

                # Normalize matrix file name inside tmpdir
                matrix_dest = os.path.join(tmp_path, "matrix.mtx.gz" if path.endswith(".gz") else "matrix.mtx")
                shutil.copy(path, matrix_dest)
                shutil.copy(barcodes, os.path.join(tmp_path, os.path.basename(barcodes)))
                shutil.copy(features, os.path.join(tmp_path, "features.tsv.gz" if features.endswith(".gz") else "features.tsv"))

                sample_name = os.path.basename(prefix)

                try:
                    adata = sc.read_10x_mtx(tmp_path, var_names='gene_symbols')
                except Exception:
                    from scipy.io import mmread
                    from anndata import AnnData 

                    matrix_file = "matrix.mtx.gz" if path.endswith(".gz") else "matrix.mtx"
                    matrix = mmread(os.path.join(tmp_path, matrix_file)).tocsr().T

                    barcodes_df = pd.read_csv(
                        os.path.join(tmp_path, os.path.basename(barcodes)), header=None
                    )
                    barcodes_list = barcodes_df[0].astype(str).tolist()

                    feature_file = "features.tsv.gz" if features.endswith(".gz") else "features.tsv"
                    genes_df = pd.read_csv(
                        os.path.join(tmp_path, feature_file), sep='\t', header=None
                    )

                    if genes_df.shape[1] != 2:
                        raise ValueError(f"Expected 2 columns in features.tsv for {sample_name}, got: {genes_df.shape[1]}")

                    gene_ids = genes_df[0].astype(str).tolist()
                    gene_symbols = genes_df[1].astype(str).tolist()

                    if matrix.shape[0] != len(barcodes_list):
                        raise ValueError(f"Mismatch between number of cells ({matrix.shape[0]}) and barcodes ({len(barcodes_list)})")
                    if matrix.shape[1] != len(gene_symbols):
                        raise ValueError(f"Mismatch between number of genes ({matrix.shape[1]}) and features ({len(gene_symbols)})")

                    adata = AnnData(X=matrix)
                    adata.obs_names = barcodes_list
                    adata.var_names = gene_symbols
                    adata.var["gene_ids"] = gene_ids


            elif os.path.isdir(path):
                adata = sc.read_10x_mtx(path, var_names='gene_symbols')
                sample_name = os.path.basename(os.path.normpath(path))

            else:
                raise ValueError(f"Unsupported input: {path}")

            group_name = sample_name.split('__')[0] if '__' in sample_name else '_'.join(sample_name.split('_')[:-1])
            adata.var_names_make_unique()
            if sample_name and len(sample_name) > 0:
                adata.obs_names = [f"{bc}.{sample_name}" for bc in adata.obs_names]
            else:
                adata.obs_names = list(adata.obs_names)  # Leave barcodes as-is

            adata.obs["sample"] = sample_name
            adata.obs["group"] = group_name
            adata.obs["Library"] = sample_name
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
        print(f"Cells remaining after mito-percent filtering: {adata_combined.n_obs}")

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
    
    print('Aligning cells to reference...',alignment_mode)
    align_start_time = time.time()
    ref_matrix = reference_df.loc[genes_present].T
    query_matrix = pd.DataFrame(
        adata_filtered.X.toarray() if sp.issparse(adata_filtered.X) else adata_filtered.X,
        index=adata_filtered.obs_names,
        columns=adata_filtered.var_names
    )

    if alignment_mode == "cosine":
        # Cosine method: L2 normalization + cosine similarity
        ref_norm = ref_matrix.div(np.linalg.norm(ref_matrix, axis=1), axis=0)
        query_norm = query_matrix.div(np.linalg.norm(query_matrix, axis=1), axis=0)
        similarities = 1 - cdist(query_norm.values, ref_norm.values, metric="cosine")
        alignment_scores = similarities[np.arange(len(similarities)), np.argmax(similarities, axis=1)]
        best_matches = np.argmax(similarities, axis=1)
        assignments = ref_matrix.index[best_matches]
        z_diff = None  # not computed

    elif alignment_mode == "classic":
        # Pearson correlation + z-score diff
        query_centered = query_matrix.sub(query_matrix.mean(axis=1), axis=0)
        ref_centered = ref_matrix.sub(ref_matrix.mean(axis=1), axis=0)
        query_std = query_centered.std(axis=1, ddof=0).replace(0, np.nan)
        ref_std = ref_centered.std(axis=1, ddof=0).replace(0, np.nan)
        norm_query = query_centered.div(query_std, axis=0)
        norm_ref = ref_centered.div(ref_std, axis=0)

        correlations = np.dot(norm_query.fillna(0).values, norm_ref.fillna(0).values.T) / query_matrix.shape[1]
        
        # Compute z-scores from correlations
        mean_corr = correlations.mean(axis=1)
        std_corr = correlations.std(axis=1)
        z_scores = (correlations - mean_corr[:, None]) / std_corr[:, None]
        
        # Assignments
        best_matches = np.argmax(z_scores, axis=1)
        second_matches = np.argsort(z_scores, axis=1)[:, -2]
        assignments = ref_matrix.index[best_matches]
        alignment_scores = correlations[np.arange(len(correlations)), best_matches]
        z_diff = z_scores[np.arange(len(z_scores)), best_matches] - z_scores[np.arange(len(z_scores)), second_matches]

    elif alignment_mode == "community":

        match_df = run_community_alignment(adata_filtered, ref_matrix)

        print("First few match_df CellBarcodes:", match_df['CellBarcode'].head().tolist())
        print("First few adata_combined.obs_names:", adata_combined.obs_names[:5].tolist())
        print("Do any barcodes match exactly?:", any(bc in match_df['CellBarcode'].values for bc in adata_combined.obs_names))

        assignments = match_df.set_index('CellBarcode').loc[adata_combined.obs_names, ref_name].values
        alignment_scores = match_df.set_index('CellBarcode').loc[adata_combined.obs_names, 'Similarity'].values
        z_diff = None
    else:
        raise ValueError(f"Invalid alignment_mode: {alignment_mode}")

    match_df = pd.DataFrame({
        "CellBarcode": query_matrix.index,
        ref_name: assignments,
        "AlignmentScore": alignment_scores
    })

    if alignment_mode == "community":
        ordered_match_df = match_df.sort_values(ref_name, ascending=False)
    else:
        ordered_match_df = pd.concat([
            match_df[match_df[ref_name] == pop].sort_values(ref_name, ascending=False)
            for pop in reference_df.columns if pop in match_df[ref_name].values
        ], ignore_index=True)


    adata_combined = adata_combined[match_df.CellBarcode].copy()
    adata_combined.obs[ref_name] = match_df.set_index('CellBarcode').loc[adata_combined.obs_names][ref_name]

    ordered_match_df.to_csv(output_dir + "/cellHarmony_lite_assignments.txt", sep="\t", index=False)
    print(f"Cosine similarity computation completed in {time.time() - align_start_time:.2f} seconds.")
    print("Assignment summary (cells per reference state):")
    print(match_df[ref_name].value_counts().to_string())

    if generate_umap:
        try:
            os.chdir(output_dir)

            def downsample_cells_per_group(adata, groupby, max_cells=50):
                idx = []
                for group, count in adata.obs[groupby].value_counts().items():
                    cells = adata.obs_names[adata.obs[groupby] == group]
                    selected = np.random.choice(cells, min(len(cells), max_cells), replace=False)
                    idx.extend(selected)
                return adata[idx].copy()

            adata_filtered = adata_combined[adata_combined.obs[ref_name].isin(
                adata_combined.obs[ref_name].value_counts()[lambda x: x >= 10].index
            )].copy()
            adata_filtered = downsample_cells_per_group(adata_filtered, ref_name, max_cells=150)

            sc.tl.rank_genes_groups(adata_filtered, groupby=ref_name, method='wilcoxon', use_raw=False)
            save_marker_genes(adata_filtered, ref_name, os.path.join(output_dir, 'supervised_markers.txt'))

            deg_results = pd.DataFrame(adata_filtered.uns['rank_genes_groups']['names'])
            markers = list(set([gene for col in deg_results.columns for gene in deg_results[col][:10]]))
            excluded_prefixes = ('mt-', 'rp', 'xist')
            markers = [gene for gene in markers if not gene.lower().startswith(excluded_prefixes) and gene in adata_filtered.var_names]

            adata_markers = adata_filtered[:, markers].copy()
            sc.pp.pca(adata_markers, n_comps=50)
            sc.pp.neighbors(adata_markers)
            sc.tl.umap(adata_markers)

            adata_all_cells = adata_combined[:, markers].copy()
            sc.pp.pca(adata_all_cells, n_comps=50)
            sc.pp.neighbors(adata_all_cells)
            sc.tl.umap(adata_all_cells)

            coords = adata_all_cells.obsm['X_umap']
            adata_combined.obs['UMAP-X'] = coords[:, 0]
            adata_combined.obs['UMAP-Y'] = coords[:, 1]
            adata_combined.obsm['X_umap'] = coords

            sc.pl.umap(adata_combined, color=ref_name, 
                save=f"_UMAP.pdf", show=False, legend_loc='on data',
                legend_fontsize=3,legend_fontweight='normal')

            sc.pl.rank_genes_groups_heatmap(
                adata_filtered,
                show=False,
                save=f"_heatmap.pdf",
                standard_scale='var',
                dendrogram=False,
                swap_axes=True,
                var_group_rotation=90,
            )
        except Exception as e:
            print(f"UMAP generation step skipped due to error: {e}")


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
            adata_combined.obs['scanpy-leiden'] = adata_unsup.obs['leiden']
        adata_combined.write(os.path.join(output_dir, "combined_with_umap_and_markers.h5ad"), compression="gzip")


    print(f"Analysis completed in {time.time() - start_time:.2f} seconds.")
    return ordered_match_df

def get_h5_and_mtx_files(folder_path):
    """
    Scans a directory for:
    - .h5 files
    - triplet-formatted 10x files using *_filtered_matrix.mtx.gz, *_counts.mtx.gz, or *_matrix.mtx.gz
    - Accepts either *_features.tsv.gz or *_genes.tsv.gz
    - .tar.gz or .tgz archives with filtered_feature_bc_matrix/

    Returns a list of .h5 files and representative paths to mtx inputs.
    """

    if '.' in folder_path[-5:]:
        return [folder_path]

    import tarfile
    import tempfile

    h5_files = sorted(glob(os.path.join(folder_path, "*.h5")))
 
    h5ad_files = sorted(glob(os.path.join(folder_path, "*.h5ad")))

    mtx_paths = []

    # Define suffix patterns to match
    suffixes = [
        "_filtered_matrix.mtx.gz", "_counts.mtx.gz", "_matrix.mtx.gz",
        "filtered_matrix.mtx.gz", "matrix.mtx.gz", "matrix.mtx"
    ]

    # Scan for all matching matrix files
    for suffix in suffixes:
        matrix_files = sorted(glob(os.path.join(folder_path, f"*{suffix}")))
        for matrix_path in matrix_files:
            prefix = matrix_path[:-len(suffix)]
            
            # Determine corresponding barcodes and features files
            if matrix_path.endswith(".gz"):
                barcode_file = prefix + "barcodes.tsv.gz"
                feature_options = [prefix + "features.tsv.gz", prefix + "genes.tsv.gz"]
            else:
                barcode_file = prefix + "barcodes.tsv"
                feature_options = [prefix + "features.tsv", prefix + "genes.tsv"]

            # Find existing feature file
            features = next((f for f in feature_options if os.path.exists(f)), None)

            # Check for barcode and feature file existence
            if os.path.exists(barcode_file) and features:
                mtx_paths.append(matrix_path)

    # Extract tar.gz/tgz archives with 10x-formatted content
    tarballs = sorted(glob(os.path.join(folder_path, "*.tar.gz")) + glob(os.path.join(folder_path, "*.tgz")))
    for archive in tarballs:
        with tarfile.open(archive, "r:gz") as tar:
            members = tar.getnames()
            for archive in tarballs:
                with tarfile.open(archive, "r:gz") as tar:
                    members = tar.getnames()
                    extract_dirs = set(os.path.dirname(m) for m in members if m.endswith((".mtx.gz", ".tsv.gz")))
                    for subdir in extract_dirs:
                        if (
                            any(m.endswith(f"{subdir}/matrix.mtx.gz") for m in members) and
                            any(m.endswith(f"{subdir}/barcodes.tsv.gz") for m in members) and
                            any(m.endswith(f"{subdir}/features.tsv.gz") for m in members) or
                            any(m.endswith(f"{subdir}/genes.tsv.gz") for m in members)
                        ):
                            tmpdir = tempfile.mkdtemp(prefix="mtx_")
                            tar.extractall(path=tmpdir)
                            extracted_path = os.path.join(tmpdir, subdir)
                            mtx_paths.append(extracted_path)
                            break  # assume only one valid subdir per archive
    return h5_files + mtx_paths + h5ad_files

########## Dedicated Functions for cellHarony-community re-implementation

def run_community_alignment(adata, ref_matrix, num_neighbors=10, resolution=1.0, louvain_level=0):
    """
    Community alignment using sparse KNN + hierarchical Louvain + correlation to reference.

    Parameters:
        adata: AnnData (filtered, normalized)
        ref_matrix: DataFrame (cell states x genes)
        num_neighbors: int, neighbors for KNN
        resolution: float, Louvain resolution
        louvain_level: int, hierarchical level for Louvain

    Returns:
        match_df: DataFrame with CellBarcode, AssignedCellState, Similarity
    """
    import networkx as nx
    import community as community_louvain
    import scipy.sparse as sp
    import numpy as np
    import pandas as pd

    genes = ref_matrix.columns.intersection(adata.var_names)
    if len(genes) == 0:
        raise ValueError("No overlap between reference genes and AnnData var_names.")

    query_X = adata[:, genes].X
    if not sp.issparse(query_X):
        query_X = sp.csr_matrix(query_X)

    """
    from sklearn.neighbors import NearestNeighbors
    knn = NearestNeighbors(n_neighbors=num_neighbors, metric="cosine").fit(query_X)
    knn_graph = knn.kneighbors_graph(query_X, mode="connectivity")
    G = nx.from_scipy_sparse_array(knn_graph)
    """

    from annoy import AnnoyIndex
    neighbor_dict = {}
    f = query_X.shape[1]
    index = AnnoyIndex(f, 'angular')
    for i in range(query_X.shape[0]):
        index.add_item(i, query_X[i].toarray().ravel())
    index.build(100)  # num_trees = 100
    for i in range(query_X.shape[0]):
        neighbor_dict[i] = index.get_nns_by_item(i, num_neighbors)
    G = nx.from_dict_of_lists(neighbor_dict)
    # Run Louvain using hierarchical dendrogram
    dendrogram = community_louvain.generate_dendrogram(G, resolution=resolution)
    level = max(0, min(louvain_level, len(dendrogram) - 1))
    partition = community_louvain.partition_at_level(dendrogram, level)

    cluster_df = pd.DataFrame({'cluster': list(partition.values())}, index=adata.obs_names)
    n_query_partitions = len(set(partition.values()))

    # Reference clustering using Annoy + Louvain
    ref_X = ref_matrix[genes].values
    ref_neighbor_dict = {}
    f_ref = ref_X.shape[1]
    ref_index = AnnoyIndex(f_ref, 'angular')
    for i in range(ref_X.shape[0]):
        ref_index.add_item(i, ref_X[i])
    ref_index.build(100)
    for i in range(ref_X.shape[0]):
        ref_neighbor_dict[i] = ref_index.get_nns_by_item(i, num_neighbors)
    G_ref = nx.from_dict_of_lists(ref_neighbor_dict)
    dendrogram_ref = community_louvain.generate_dendrogram(G_ref, resolution=resolution)
    level_ref = max(0, min(louvain_level, len(dendrogram_ref) - 1))
    ref_partition = community_louvain.partition_at_level(dendrogram_ref, level_ref)
    n_ref_partitions = len(set(ref_partition.values()))

    print(f"number of reference partitions {n_ref_partitions}, number of query partitions {n_query_partitions}")

    centroids = []
    cluster_ids = []

    for cluster_id in sorted(cluster_df['cluster'].unique()):
        cell_idx = cluster_df[cluster_df['cluster'] == cluster_id].index
        rows = [adata.obs_names.get_loc(bc) for bc in cell_idx]
        subset = query_X[rows, :]
        centroid = subset.mean(axis=0).A1
        centroids.append(centroid)
        cluster_ids.append(cluster_id)

    centroids = np.vstack(centroids)

    # Correlation to reference
    centroid_norm = (centroids - centroids.mean(axis=1, keepdims=True)) / centroids.std(axis=1, keepdims=True)
    ref_dense = ref_matrix[genes].values
    ref_norm = (ref_dense - ref_dense.mean(axis=1, keepdims=True)) / ref_dense.std(axis=1, keepdims=True)

    centroid_norm = np.nan_to_num(centroid_norm)
    ref_norm = np.nan_to_num(ref_norm)

    corr_matrix = np.dot(centroid_norm, ref_norm.T) / centroid_norm.shape[1]

    best_idx = np.argmax(corr_matrix, axis=1)
    best_labels = ref_matrix.index[best_idx]
    best_scores = corr_matrix[np.arange(len(best_idx)), best_idx]

    records = []
    for cid, label, score in zip(cluster_ids, best_labels, best_scores):
        barcodes = cluster_df[cluster_df['cluster'] == cid].index
        for bc in barcodes:
            records.append({
                'CellBarcode': bc,
                ref_name: label,
                'Similarity': score
            })

    return pd.DataFrame(records)

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
    parser.add_argument('--append_obs', type=str, default=None, help='Field in .obs to append to the cell barcode (e.g., donor_id)')
    parser.add_argument('--alignment_mode', type=str, default="cosine", help='Alignment mode: "cosine" or "classic"')
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
    alignment_mode = args.alignment_mode
    unsupervised_cluster = args.unsupervised_cluster
    append_obs_field = args.append_obs

    #h5_files = glob(os.path.join(h5_directory, "*.h5")) if '.h5' not in h5_directory else [h5_directory]
    h5_files = get_h5_and_mtx_files(h5_directory)

    if len(h5_files)==0:
        print ("No compatible h5, h5ad or .mtx files identified")
        sys.exit()
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
        unsupervised_cluster=unsupervised_cluster,
        append_obs_field=append_obs_field,
        alignment_mode=alignment_mode
    )
