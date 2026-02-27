import os, time, sys, shutil, re
import tempfile
from pathlib import Path
from datetime import datetime
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
import matplotlib.colors as mcolors
import warnings
warnings.filterwarnings("ignore", message="Variable names are not unique. To make them unique, call `.var_names_make_unique`.")
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*invalid value encountered in log2.*")
warnings.filterwarnings("ignore", category=FutureWarning, message=".*default backend for leiden will be igraph.*")
warnings.filterwarnings("ignore", category=UserWarning, message="Received a view of an AnnData. Making a copy.")

plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['figure.facecolor'] = 'white'


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
        return len(data)

    def flush(self):
        for stream in self.streams:
            stream.flush()

    @property
    def encoding(self):
        for stream in self.streams:
            if hasattr(stream, "encoding"):
                return stream.encoding
        return "utf-8"


def normalize_adata(adata, show_progress=False):
    """Log-normalize the AnnData object in-place."""
    if show_progress:
        with tqdm(total=2, desc="Normalization steps") as pbar:
            sc.pp.normalize_total(adata, target_sum=1e4); pbar.update(1)
            sc.pp.log1p(adata); pbar.update(1)
    else:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

def save_marker_genes(adata, groupby, output_file):
    deg = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    deg.to_csv(output_file, sep='\t', index=False)


def ensure_category_palette(adata, column):
    if column not in adata.obs:
        return
    obs_col = adata.obs[column].astype(str)
    adata.obs[column] = pd.Categorical(obs_col)
    categories = adata.obs[column].cat.categories
    num_cats = len(categories)
    base_palette = list(sc.pl.palettes.default_102)
    if num_cats <= len(base_palette):
        palette = base_palette[:num_cats]
    else:
        cmap = plt.get_cmap("gist_ncar", num_cats)
        palette = [mcolors.to_hex(cmap(i)) for i in range(num_cats)]
    adata.uns[f"{column}_colors"] = palette

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
    alignment_mode="classic",
    min_alignment_score=None,
    gene_translation_file=None,
    metacell_align=False,
    metacell_target_size=50,
    metacell_min_size=25,
    metacell_max_size=100,
    metacell_algorithm="kmeans",
    metacell_neighbors=30,
    metacell_hvg=3000,
    metacell_pcs=50,
    metacell_random_count=50,
    metacell_random_cells=5,
    metacell_random_replacement=False,
    metacell_random_state=0,
    ambient_correct_cutoff=None,
    concat_on_disk=False,
    concat_batch_size=50,
    verbose_import=False
):
    start_time = time.time()
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    def _rss_mb():
        import resource
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform == "darwin":
            return rss / (1024 * 1024)
        return rss / 1024.0

    def _log_step(message, elapsed=None):
        timing = f" time={elapsed:.2f}s" if elapsed is not None else ""
        print(f"[mem][timing] {message}{timing} rss={_rss_mb():.1f}MB")

    def _filter_cells_min_genes(adata, min_genes):
        if min_genes is None:
            return
        X = adata.X
        if sp.issparse(X):
            X = X.tocsr()
            n_genes = np.diff(X.indptr)
        else:
            n_genes = np.sum(X > 0, axis=1)
        keep = np.asarray(n_genes).ravel() >= min_genes
        adata._inplace_subset_obs(keep)

    def _filter_genes_min_cells(adata, min_cells):
        if min_cells is None:
            return
        X = adata.X
        if sp.issparse(X):
            X = X.tocsr()
            n_cells = np.bincount(X.indices, minlength=X.shape[1])
        else:
            n_cells = np.sum(X > 0, axis=0)
        keep = np.asarray(n_cells).ravel() >= min_cells
        adata._inplace_subset_var(keep)

    reference_df = pd.read_csv(cellharmony_ref, sep='\t', index_col=0)
    cell_populations = reference_df.columns.tolist()
    ref_name = os.path.basename(cellharmony_ref)[:-4]

    # Optional gene translation table (e.g., Ensembl → Symbol)
    translation_map = None
    if gene_translation_file is not None and os.path.exists(gene_translation_file):
        translation_map = load_gene_translation(gene_translation_file)


    if h5ad_file is not None:
        adata_combined = sc.read_h5ad(h5ad_file)
        apply_gene_translation(adata_combined, translation_map, os.path.basename(h5ad_file))
        adata_combined.var_names_make_unique()
        print(f"reimported adata shape: {adata_combined.shape} (cells x genes)")
    else:
        if concat_batch_size and not concat_on_disk:
            print("[warn] --concat_batch_size ignored unless --concat_on_disk is set.")
            concat_batch_size = 0
        use_on_disk = bool(concat_on_disk)
        adata_list = []
        load_start = time.time()
        _log_step("start loading inputs")

        def _load_adata(path, sample_name_override=None):
            if path.endswith(".h5"):
                sample_name = os.path.basename(path).replace(".h5", "")
                adata_local = sc.read_10x_h5(path)
            elif path.endswith(".h5ad"):
                adata_local = sc.read_h5ad(path)
                sample_name = os.path.basename(path).replace(".h5ad", "")
                group_name = sample_name.split('__')[0] if '__' in sample_name else '_'.join(sample_name.split('_')[:-1])
                if append_obs_field is not None:
                    if append_obs_field not in adata_local.obs.columns:
                        raise ValueError(f"Field '{append_obs_field}' not found in adata.obs columns for {sample_name}.")
                    appended_values = adata_local.obs[append_obs_field].astype(str).tolist()
                    adata_local.obs_names = [f"{bc}.{val}" for bc, val in zip(adata_local.obs_names, appended_values)]
            elif path.endswith((
                "_filtered_matrix.mtx.gz", "_counts.mtx.gz", "_matrix.mtx.gz",
                "filtered_matrix.mtx.gz", "matrix.mtx.gz", "matrix.mtx", "quants_mat.mtx"
            )):
                suffixes = [
                    "_filtered_matrix.mtx.gz", "_counts.mtx.gz", "_matrix.mtx.gz",
                    "filtered_matrix.mtx.gz", "matrix.mtx.gz", "matrix.mtx", "quants_mat.mtx"
                ]
                suffix = next(s for s in suffixes if path.endswith(s))
                prefix = path[:-len(suffix)]
                if path.endswith(".gz"):
                    barcodes = prefix + ("_barcodes.tsv.gz" if any(path.endswith(f"_{s}") for s in ["filtered_matrix.mtx.gz", "counts.mtx.gz", "matrix.mtx.gz"]) else "barcodes.tsv.gz")
                    feature_options = [prefix + "_features.tsv.gz", prefix + "_genes.tsv.gz",
                                    prefix + "features.tsv.gz", prefix + "genes.tsv.gz"]
                else:
                    barcodes = prefix + (
                        "_barcodes.tsv"
                        if any(path.endswith(f"_{s}") for s in ["matrix.mtx"])
                        else (
                            "quants_mat_rows.txt"
                            if path.endswith("quants_mat.mtx") 
                            else "barcodes.tsv"
                        )
                    )
                    feature_options = [
                        prefix + "_features.tsv", prefix + "_genes.tsv",
                        prefix + "features.tsv", prefix + "genes.tsv",
                        prefix + "quants_mat_cols.txt" 
                    ]
                features = next((f for f in feature_options if os.path.exists(f)), None)
                if not features:
                    raise FileNotFoundError(f"Missing features or genes file for prefix: {prefix}")
                if not os.path.exists(barcodes):
                    raise FileNotFoundError(f"Missing barcodes file for prefix: {prefix}")
                tmpdir = tempfile.TemporaryDirectory()
                tmp_path = tmpdir.name
                matrix_dest = os.path.join(tmp_path, "matrix.mtx.gz" if path.endswith(".gz") else "matrix.mtx")
                shutil.copy(path, matrix_dest)
                shutil.copy(barcodes, os.path.join(tmp_path, os.path.basename(barcodes)))
                shutil.copy(features, os.path.join(tmp_path, "features.tsv.gz" if features.endswith(".gz") else "features.tsv"))
                if os.path.basename(path) in ("matrix.mtx", "matrix.mtx.gz"):
                    base_dir = os.path.dirname(path)
                    sample_name = os.path.basename(os.path.normpath(base_dir))
                else:
                    sample_name = os.path.basename(prefix)
                try:
                    adata_local = sc.read_10x_mtx(tmp_path, var_names='gene_symbols')
                except Exception:
                    from scipy.io import mmread
                    from anndata import AnnData 
                    matrix_file = "matrix.mtx.gz" if path.endswith(".gz") else "matrix.mtx"
                    matrix = mmread(os.path.join(tmp_path, matrix_file)).tocsr()
                    if "quants_mat" not in os.path.basename(path):
                        matrix = matrix.T
                    barcodes_df = pd.read_csv(
                        os.path.join(tmp_path, os.path.basename(barcodes)), header=None
                    )
                    barcodes_list = barcodes_df[0].astype(str).tolist()
                    feature_file = "features.tsv.gz" if features.endswith(".gz") else "features.tsv"
                    genes_df = pd.read_csv(
                        os.path.join(tmp_path, feature_file), sep='\t', header=None
                    )
                    if genes_df.shape[1] != 2:
                        gene_ids = genes_df[0].astype(str).tolist()
                        gene_symbols = genes_df[0].astype(str).tolist()
                    else:
                        gene_ids = genes_df[0].astype(str).tolist()
                        gene_symbols = genes_df[1].astype(str).tolist()
                    if matrix.shape[0] != len(barcodes_list):
                        raise ValueError(f"Mismatch between number of cells ({matrix.shape[0]}) and barcodes ({len(barcodes_list)})")
                    if matrix.shape[1] != len(gene_symbols):
                        raise ValueError(f"Mismatch between number of genes ({matrix.shape[1]}) and features ({len(gene_symbols)})")
                    adata_local = AnnData(X=matrix)
                    adata_local.obs_names = barcodes_list
                    adata_local.var_names = gene_symbols
                    adata_local.var["gene_ids"] = gene_ids
            elif os.path.isdir(path):
                adata_local = sc.read_10x_mtx(path, var_names='gene_symbols')
                sample_name = os.path.basename(os.path.normpath(path))
            else:
                raise ValueError(f"Unsupported input: {path}")

            if sample_name_override:
                sample_name = str(sample_name_override)

            if not issparse(adata_local.X):
                adata_local.X = sp.csr_matrix(adata_local.X)
            else:
                adata_local.X = adata_local.X.tocsr()

            apply_gene_translation(adata_local, translation_map, sample_name)
            adata_local.var_names_make_unique()

            group_name = sample_name.split('__')[0] if '__' in sample_name else '_'.join(sample_name.split('_')[:-1])
            if sample_name and len(sample_name) > 0:
                adata_local.obs_names = [f"{bc}.{sample_name}" for bc in adata_local.obs_names]
            else:
                adata_local.obs_names = list(adata_local.obs_names)

            adata_local.obs["sample"] = sample_name
            adata_local.obs["group"] = group_name
            adata_local.obs["Library"] = sample_name
            return adata_local, sample_name

        if use_on_disk:
            combined_path = os.path.join(output_dir, "combined_raw.h5ad")
            batch_size = int(concat_batch_size) if concat_batch_size else 50
            if batch_size <= 0:
                batch_size = 50
            concat_tmpdir = tempfile.TemporaryDirectory(prefix="concat_batches_", dir=output_dir)
            print(
                "[INFO] Performing memory-efficient on-disk concatenation in batches "
                f"of {batch_size}; this can take a while."
            )

            def _concat_on_disk(inputs, output_path):
                if hasattr(ad, "experimental") and hasattr(ad.experimental, "concat_on_disk"):
                    ad.experimental.concat_on_disk(
                        inputs,
                        output_path,
                        max_loaded_elems=12000000000,
                        label="sample",
                        join="outer",
                        fill_value=0
                    )
                    return
                raise RuntimeError(
                    "On-disk concatenation requires anndata.experimental.concat_on_disk. "
                    "Please upgrade anndata or disable concat_on_disk."
                )

            progress = tqdm(total=len(h5_files), desc="Loading input files")
            batch_files = []
            batch_index = 0
            concat_start = time.time()
            for start in range(0, len(h5_files), batch_size):
                batch_start = time.time()
                adata_list = []
                batch = h5_files[start:start + batch_size]
                for idx, entry in enumerate(batch, start=start + 1):
                    sample_start = time.time()
                    if isinstance(entry, (tuple, list)) and len(entry) >= 2:
                        path, sample_override = entry[0], entry[1]
                    else:
                        path, sample_override = entry, None
                    adata, sample_name = _load_adata(path, sample_name_override=sample_override)
                    adata_list.append(adata)
                    if verbose_import:
                        _log_step(
                            f"loaded {sample_name} shape={adata.shape}",
                            elapsed=time.time() - sample_start
                        )
                    progress.update(1)
                batch_concat_start = time.time()
                batch_adata = ad.concat(adata_list, label="sample", join="outer", fill_value=0)
                batch_path = os.path.join(concat_tmpdir.name, f"batch_concat_{batch_index}.h5ad")
                write_start = time.time()
                batch_adata.write(batch_path, compression="gzip")
                if verbose_import:
                    _log_step(
                        f"wrote batch {batch_index} -> {batch_path}",
                        elapsed=time.time() - write_start
                    )
                del batch_adata
                del adata_list
                batch_files.append(batch_path)
                batch_index += 1
            progress.close()
            _log_step("completed batched in-memory concat", elapsed=time.time() - concat_start)
            if len(batch_files) == 1:
                shutil.copy(batch_files[0], combined_path)
            else:
                print(f"[INFO] Merging {len(batch_files)} batch files on disk...")
                _concat_on_disk(batch_files, combined_path)
            for path in batch_files:
                try:
                    os.remove(path)
                except OSError:
                    pass
            read_start = time.time()
            adata_combined = sc.read_h5ad(combined_path)
            _log_step("reloaded combined h5ad", elapsed=time.time() - read_start)
            concat_tmpdir.cleanup()
        else:
            for entry in tqdm(h5_files, desc="Loading input files"):
                sample_start = time.time()
                if isinstance(entry, (tuple, list)) and len(entry) >= 2:
                    path, sample_override = entry[0], entry[1]
                else:
                    path, sample_override = entry, None
                adata, sample_name = _load_adata(path, sample_name_override=sample_override)
                adata_list.append(adata)
                if verbose_import:
                    _log_step(
                        f"loaded {sample_name} shape={adata.shape}",
                        elapsed=time.time() - sample_start
                    )
            _log_step("completed input loading", elapsed=time.time() - load_start)

        if not use_on_disk:
            concat_start = time.time()
            _log_step("start in-memory concat")
            adata_combined = ad.concat(adata_list, label="sample", join="outer", fill_value=0)
            _log_step("completed in-memory concat", elapsed=time.time() - concat_start)
        print(f"adata shape: {adata_combined.shape} (cells x genes)")

        if ambient_correct_cutoff is not None:
            soupx_dir = os.path.join(output_dir, "soupx")
            os.makedirs(soupx_dir, exist_ok=True)
            try:
                from altanalyze3.components.ambient_rna import soupx_correct
            except ImportError as exc:
                raise RuntimeError(
                    "SoupX ambient correction requested (--ambient_correct_cutoff) but the soupx module "
                    "is unavailable. Install the 'soupx' package to enable this option."
                ) from exc

            print(f"[INFO] Running SoupX ambient correction (rho={ambient_correct_cutoff})...")
            corrected = soupx_correct.process_anndata(
                adata_combined,
                rho=float(ambient_correct_cutoff),
                library_col="Library",
                outdir=Path(soupx_dir),
                write_individual=False,
                merged_filename="soupx_corrected_merged.h5ad",
            )
            if corrected is None:
                raise RuntimeError("SoupX correction did not return a corrected AnnData object.")
            adata_combined = corrected
            adata_combined.var_names_make_unique()
            print(f"[INFO] SoupX correction complete. Corrected adata shape: {adata_combined.shape} (cells x genes)")

        if sp.issparse(adata_combined.X):
            adata_combined.X = adata_combined.X.tocsr()
        _filter_cells_min_genes(adata_combined, min_genes)
        print(f"Cells remaining after min_genes {min_genes} filtering: {adata_combined.n_obs}")
        _filter_genes_min_cells(adata_combined, min_cells)
        sc.pp.filter_cells(adata_combined, min_counts=min_counts)
        print(f"Cells remaining after min_counts {min_counts} filtering: {adata_combined.n_obs}")

        mito_genes = adata_combined.var_names.str.upper().str.startswith("MT-")
        adata_combined.obs["pct_counts_mt"] = (
            np.sum(adata_combined[:, mito_genes].X, axis=1).A1 /
            np.sum(adata_combined.X, axis=1).A1
        ) * 100

        adata_combined = adata_combined[adata_combined.obs["pct_counts_mt"] < mit_percent].copy()
        print(f"Cells remaining after mito-percent filtering: {adata_combined.n_obs}")

    original_cell_adata = None
    metacell_membership = None
    metacell_aligned_adata = None
    cell_to_metacell_map = None
    cell_metacell_lookup = None

    if metacell_align:
        from altanalyze3.components.metacells.main import MetacellParams, assemble_metacells

        original_cell_adata = adata_combined.copy()

        metacell_params = MetacellParams(
            target_size=metacell_target_size,
            min_size=metacell_min_size,
            max_size=metacell_max_size,
            preserve_small=False,
            aggregation="sum",
            graph_algorithm=metacell_algorithm,
            resolution=None,
            resolution_steps=5,
            resolution_tolerance=0.15,
            n_neighbors=max(5, metacell_neighbors),
            neighbor_method="auto",
            neighbor_metric="euclidean",
            n_top_genes=max(200, metacell_hvg),
            n_pcs=max(10, metacell_pcs),
            pca_svd_solver="randomized",
            hvg_layer=None,
            expression_layer=None,
            use_raw=False,
            random_state=metacell_random_state,
            random_metacell_count=metacell_random_count if metacell_algorithm == "random" else None,
            random_cells_per_metacell=metacell_random_cells,
            random_sampling_with_replacement=metacell_random_replacement,
        )

        mc_start = time.time()

        metacell_adata, membership_df = assemble_metacells(
            adata_combined,
            params=metacell_params,
            boundary_columns=[],
            sample_column=None,
            cell_type_column=None,
        )

        mc_end = time.time()
        mc_duration = mc_end - mc_start
        n_cells = original_cell_adata.n_obs
        n_metacells = metacell_adata.n_obs
        size_stats = metacell_adata.obs['metacell_size']

        print(f"[metacell] generated {n_metacells} metacells from {n_cells} cells in {mc_duration:.1f}s")
        print(f"[metacell] size median={size_stats.median():.1f} range=({int(size_stats.min())}, {int(size_stats.max())})")

        if "gene_symbols" not in metacell_adata.var.columns:
            metacell_adata.var["gene_symbols"] = metacell_adata.var_names.astype(str)

        print(f"[metacell] algorithm={metacell_algorithm} target={metacell_target_size} min={metacell_min_size} max={metacell_max_size} neighbors={metacell_neighbors} hvg={metacell_hvg} pcs={metacell_pcs}")
        metacell_adata.uns["metacell_membership"] = membership_df
        metacell_path = os.path.join(output_dir, "metacells.h5ad")
        metacell_adata.write(metacell_path, compression="gzip")

        adata_combined = metacell_adata.copy()
        metacell_membership = membership_df.copy()
        cell_to_metacell_map = membership_df[['cell_barcode', 'metacell_id']].set_index('cell_barcode')['metacell_id']

    # retain the original counts in the h5ad in the new counts slot
    adata_combined.layers["counts"] = adata_combined.X.copy()

    normalize_adata(adata_combined, show_progress=True)
    metacell_aligned_adata = adata_combined.copy() if metacell_align else None

    if export_h5ad:
        adata_combined.write(output_dir+"/combined_qc_normalized.h5ad", compression="gzip")

    marker_genes = reference_df.index
    adata_filtered, genes_present, missing_genes = subset_to_reference_genes(adata_combined, marker_genes)

    if missing_genes:
        print(f"Warning: {len(missing_genes)} marker genes not found in dataset and will be excluded out of {len(marker_genes)}.")
        print(f"(First 10) Missing genes: {sorted(missing_genes)[:10]}")

    if export_cptt and not metacell_align:
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
    query_index = adata_filtered.obs_names

    if alignment_mode == "cosine":
        # Cosine method: L2 normalization + cosine similarity
        query_X = adata_filtered.X
        if sp.issparse(query_X):
            query_X = query_X.tocsr()
            row_norms = np.sqrt(query_X.multiply(query_X).sum(axis=1)).A1
            row_norms[row_norms == 0] = 1.0
            query_norm = sp.diags(1.0 / row_norms).dot(query_X)
        else:
            row_norms = np.linalg.norm(query_X, axis=1)
            row_norms[row_norms == 0] = 1.0
            query_norm = query_X / row_norms[:, None]

        ref_values = ref_matrix.values
        ref_norms = np.linalg.norm(ref_values, axis=1)
        ref_norms[ref_norms == 0] = 1.0
        ref_norm = ref_values / ref_norms[:, None]

        similarities = query_norm.dot(ref_norm.T)
        similarities = np.asarray(similarities)
        alignment_scores = similarities[np.arange(len(similarities)), np.argmax(similarities, axis=1)]
        best_matches = np.argmax(similarities, axis=1)
        assignments = ref_matrix.index[best_matches]
        z_diff = None  # not computed

    elif alignment_mode == "classic":
        # Pearson correlation + z-score diff
        query_matrix = pd.DataFrame(
            adata_filtered.X.toarray() if sp.issparse(adata_filtered.X) else adata_filtered.X,
            index=adata_filtered.obs_names,
            columns=adata_filtered.var_names
        )
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
        "CellBarcode": query_index,
        ref_name: assignments,
        "AlignmentScore": alignment_scores
    })

    if min_alignment_score is not None:
        before = match_df.shape[0]
        match_df = match_df[match_df["AlignmentScore"] >= min_alignment_score].copy()
        after = match_df.shape[0]
        print(f"[INFO] Applied min_alignment_score={min_alignment_score}. "
            f"Excluded {before - after} cells, kept {after}.")

    if metacell_align:
        if metacell_membership is None or original_cell_adata is None:
            raise RuntimeError("Metacell membership data not available for propagation.")

        metacell_assignments = match_df.set_index("CellBarcode")[ref_name]
        metacell_scores = match_df.set_index("CellBarcode")["AlignmentScore"]

        membership_df = metacell_membership.copy()
        membership_df = membership_df.rename(columns={"cell_barcode": "CellBarcode"})
        membership_df[ref_name] = membership_df["metacell_id"].map(metacell_assignments)
        membership_df["AlignmentScore"] = membership_df["metacell_id"].map(metacell_scores)
        membership_df = membership_df.dropna(subset=[ref_name])

        available_cells = original_cell_adata.obs_names.intersection(membership_df["CellBarcode"])
        membership_df = membership_df[membership_df["CellBarcode"].isin(available_cells)]

        match_df = membership_df[["CellBarcode", ref_name, "AlignmentScore"]].copy()
        match_df = match_df.drop_duplicates(subset="CellBarcode")

        cell_metacell_lookup = membership_df.set_index("CellBarcode")["metacell_id"]
        cell_to_metacell_map = cell_metacell_lookup

        adata_combined = original_cell_adata[match_df["CellBarcode"]].copy()
        adata_combined.layers["counts"] = adata_combined.X.copy()
        normalize_adata(adata_combined, show_progress=False)

        marker_genes = reference_df.index
        adata_filtered, genes_present, _ = subset_to_reference_genes(adata_combined, marker_genes)
        if alignment_mode == "classic":
            query_matrix = pd.DataFrame(
                adata_filtered.X.toarray() if sp.issparse(adata_filtered.X) else adata_filtered.X,
                index=adata_filtered.obs_names,
                columns=adata_filtered.var_names
            )

    if alignment_mode == "community":
        ordered_match_df = match_df.sort_values(ref_name, ascending=False)
    else:
        ordered_match_df = pd.concat([
            match_df[match_df[ref_name] == pop].sort_values(ref_name, ascending=False)
            for pop in reference_df.columns if pop in match_df[ref_name].values
        ], ignore_index=True)

    adata_combined = adata_combined[match_df.CellBarcode].copy()
    adata_combined.obs[ref_name] = match_df.set_index('CellBarcode').loc[adata_combined.obs_names][ref_name]

    if export_cptt and metacell_align:
        cptt_df = pd.DataFrame(
            adata_filtered.X,
            index=adata_filtered.obs_names,
            columns=adata_filtered.var_names
        ).T
        cptt_df.insert(0, "UID", cptt_df.index)
        cptt_df.to_csv(output_dir+"/CPTT_matrix.txt", sep="\t", index=False)

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

            sc.pp.pca(adata_all_cells, n_comps=50) ### Segmentation fault: 11 - leaked semaphore - if library mismatch

            sc.pp.neighbors(adata_all_cells)
            sc.tl.umap(adata_all_cells)

            coords = adata_all_cells.obsm['X_umap']
            adata_combined.obs['UMAP-X'] = coords[:, 0]
            adata_combined.obs['UMAP-Y'] = coords[:, 1]
            adata_combined.obsm['X_umap'] = coords

            ensure_category_palette(adata_combined, ref_name)
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
        if metacell_align:
            adata_unsup_metacell = metacell_aligned_adata.copy()
            sc.pp.highly_variable_genes(adata_unsup_metacell, min_mean=0.0125, max_mean=3, min_disp=0.5)
            adata_unsup_metacell = adata_unsup_metacell[:, adata_unsup_metacell.var.highly_variable].copy()
            sc.pp.scale(adata_unsup_metacell, max_value=10)
            sc.pp.pca(adata_unsup_metacell, n_comps=50)
            sc.pp.neighbors(adata_unsup_metacell)
            sc.tl.umap(adata_unsup_metacell)
            sc.tl.leiden(adata_unsup_metacell, flavor="leidenalg")
            print('[metacell] unsupervised clustering performed on metacells')
            ensure_category_palette(adata_unsup_metacell, 'leiden')

            metacell_cluster_map = adata_unsup_metacell.obs['leiden'].astype(str).to_dict()
            pd.DataFrame({
                "metacell_id": adata_unsup_metacell.obs_names,
                "leiden": adata_unsup_metacell.obs['leiden'].astype(str).values
            }).to_csv(os.path.join(output_dir, "metacell_leiden_clusters.tsv"), sep="\t", index=False)

            if cell_to_metacell_map is not None:
                cell_clusters = [
                    metacell_cluster_map.get(cell_to_metacell_map.get(bc), "NA")
                    for bc in adata_combined.obs_names
                ]
            else:
                cell_clusters = ["NA"] * adata_combined.n_obs
            adata_combined.obs['scanpy-leiden'] = pd.Series(cell_clusters, index=adata_combined.obs_names).astype(str)

            adata_unsup_cells = adata_combined.copy()
            sc.pp.highly_variable_genes(adata_unsup_cells, min_mean=0.0125, max_mean=3, min_disp=0.5)
            adata_unsup_cells = adata_unsup_cells[:, adata_unsup_cells.var.highly_variable].copy()
            sc.pp.scale(adata_unsup_cells, max_value=10)
            sc.pp.pca(adata_unsup_cells, n_comps=50)
            sc.pp.neighbors(adata_unsup_cells)
            sc.tl.umap(adata_unsup_cells)
            adata_unsup_cells.obs['leiden'] = adata_combined.obs['scanpy-leiden'].reindex(adata_unsup_cells.obs_names).astype(str)

            adata_combined.obsm['X_umap_leiden'] = adata_unsup_cells.obsm['X_umap']
            adata_combined.obs['UMAP_leiden_X'] = adata_unsup_cells.obsm['X_umap'][:, 0]
            adata_combined.obs['UMAP_leiden_Y'] = adata_unsup_cells.obsm['X_umap'][:, 1]

            unique_clusters = adata_unsup_cells.obs['leiden'].dropna().unique()
            if len(unique_clusters) > 1:
                sc.tl.rank_genes_groups(adata_unsup_cells, groupby='leiden', method='wilcoxon', use_raw=False)
                save_marker_genes(adata_unsup_cells, 'leiden', os.path.join(output_dir, 'unsupervised_markers.txt'))
                ensure_category_palette(adata_unsup_cells, 'leiden')

                sc.pl.umap(
                    adata_unsup_cells,
                    color='leiden',
                    save="_unsupervised_umap.pdf",
                    show=False,
                    legend_loc='on data',
                    legend_fontsize=5,
                    legend_fontweight='normal',
                )

                sc.pl.rank_genes_groups_heatmap(
                    adata_unsup_cells,
                    groupby='leiden',
                    n_genes=5,
                    show=False,
                    save="_unsupervised_heatmap.pdf",
                    standard_scale='var',
                    dendrogram=False,
                    swap_axes=True,
                    var_group_rotation=0,
                )
            else:
                print("[INFO] Skipping unsupervised marker discovery (insufficient metacell clusters).")
            adata_unsup = adata_unsup_cells
        else:
            adata_unsup = adata_combined.copy()
            sc.pp.highly_variable_genes(adata_unsup, min_mean=0.0125, max_mean=3, min_disp=0.5)
            adata_unsup = adata_unsup[:, adata_unsup.var.highly_variable].copy()
            sc.pp.scale(adata_unsup, max_value=10)
            sc.pp.pca(adata_unsup, n_comps=50)
            sc.pp.neighbors(adata_unsup)
            sc.tl.umap(adata_unsup)
            sc.tl.leiden(adata_unsup, flavor="leidenalg", resolution=0.5)
            n_clusters = adata_unsup.obs['leiden'].nunique()
            n_cells = adata_unsup.n_obs
            print(f"[cell] Unsupervised clustering completed: {n_clusters} clusters from {n_cells} cells (resolution=0.5)")
            print('[cell] unsupervised clustering performed on individual cells')

            sc.tl.rank_genes_groups(adata_unsup, groupby='leiden', method='wilcoxon', use_raw=False)
            save_marker_genes(adata_unsup, 'leiden', os.path.join(output_dir, 'unsupervised_markers.txt'))

            ensure_category_palette(adata_unsup, 'leiden')
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
            adata_combined.obs['scanpy-leiden'] = adata_unsup.obs['leiden']

    if save_adata:
        # If unsupervised clustering was run, merge leiden clusters into adata_combined before saving
        if unsupervised_cluster and 'leiden' in adata_unsup.obs.columns:
            adata_combined.obs['scanpy-leiden'] = adata_unsup.obs['leiden']
            coords_unsup = adata_unsup.obsm['X_umap']
            adata_combined.obsm['X_umap_leiden'] = coords_unsup
            adata_combined.obs['UMAP_leiden_X'] = coords_unsup[:, 0]
            adata_combined.obs['UMAP_leiden_Y'] = coords_unsup[:, 1]

            # Export leiden results to a tsv
            leiden_export = pd.DataFrame({
                "CellBarcode": adata_unsup.obs_names,
                "Leiden": adata_unsup.obs['leiden'].values,
                "UMAP_X": coords_unsup[:, 0],
                "UMAP_Y": coords_unsup[:, 1]
            })
            leiden_export.to_csv(
                os.path.join(output_dir, "unsupervised_leiden_clusters.tsv"),
                sep="\t", index=False
            )

        adata_combined.uns["lineage_order"] = reference_df.columns.tolist()
        adata_combined.write(os.path.join(output_dir, "combined_with_umap_and_markers.h5ad"), compression="gzip")


    print(f"Analysis completed in {time.time() - start_time:.2f} seconds.")
    return ordered_match_df

def load_gene_translation(translation_path):
    """Load a two-column gene translation file (e.g., Ensembl → Symbol)."""
    df = pd.read_csv(translation_path, sep='\t', header=None, names=["source", "target"])
    df = df.dropna().drop_duplicates()
    translation_map = dict(zip(df["source"].astype(str), df["target"].astype(str)))
    print(f"[INFO] Loaded gene translation table: {len(translation_map)} mappings")
    return translation_map


def apply_gene_translation(adata, translation_map, sample_label=None):
    """Apply a gene translation map to an AnnData object."""
    sample_label = sample_label or "dataset"

    if "gene_symbols" not in adata.var.columns:
        adata.var["gene_symbols"] = adata.var_names.astype(str)

    if translation_map is None:
        return 0

    original_names = adata.var_names.astype(str)
    translated_names = []
    updated_count = 0

    for name in original_names:
        translated = translation_map.get(name, name)
        if translated != name:
            updated_count += 1
        translated_names.append(str(translated))

    adata.var["source_gene_id"] = original_names
    adata.var["gene_symbols"] = translated_names
    adata.var_names = pd.Index(translated_names)

    if updated_count:
        print(f"[INFO] Applied gene translation for {sample_label}: {updated_count} IDs updated")
    else:
        print(f"[INFO] Gene translation file provided but no IDs matched for {sample_label}")

    return updated_count


def subset_to_reference_genes(adata, reference_genes):
    """Subset AnnData to genes present in the reference, using gene symbols when available."""
    if isinstance(reference_genes, pd.Index):
        reference_genes = reference_genes.tolist()

    if "gene_symbols" in adata.var.columns:
        symbol_values = adata.var["gene_symbols"].astype(str).tolist()
    else:
        symbol_values = adata.var_names.astype(str).tolist()

    symbol_to_var = {}
    for var_name, symbol in zip(adata.var_names.tolist(), symbol_values):
        if symbol is None:
            continue
        symbol_str = str(symbol)
        if not symbol_str or symbol_str.lower() == "nan":
            continue
        if symbol_str not in symbol_to_var:
            symbol_to_var[symbol_str] = var_name

    matched_symbols = [gene for gene in reference_genes if gene in symbol_to_var]
    missing_symbols = [gene for gene in reference_genes if gene not in symbol_to_var]

    if not matched_symbols:
        raise ValueError("No overlapping genes between the dataset and the reference after translation.")

    selected_var_names = [symbol_to_var[gene] for gene in matched_symbols]
    filtered = adata[:, selected_var_names].copy()
    filtered.var_names = pd.Index(matched_symbols)
    filtered.var["gene_symbols"] = matched_symbols

    return filtered, matched_symbols, missing_symbols

def get_h5_and_mtx_files(folder_path):
    """
    Recursively scans a directory for:
      - .h5 files
      - .h5ad files
      - 10x-style triplets: matrix.mtx (+ barcodes.tsv + genes.tsv/features.tsv)
      - .tar.gz/.tgz archives containing 10x-formatted content
    Returns a list of file paths (h5, h5ad, or matrix.mtx) or
    (path, sample_name) tuples when a two-column mapping file is provided.
    """

    import tarfile
    import tempfile

    results = []
    
    # Direct case: if folder_path is itself a file
    if os.path.isfile(folder_path):
        if folder_path.endswith((".h5", ".h5ad", ".mtx", ".mtx.gz")):
            return [folder_path]
        parsed_rows = []
        try:
            with open(folder_path, "r", encoding="utf-8") as handle:
                for line_num, line in enumerate(handle, start=1):
                    raw = line.strip()
                    if not raw or raw.startswith("#"):
                        continue
                    parts = raw.split("\t") if "\t" in raw else raw.split()
                    if len(parts) < 2:
                        raise ValueError(
                            f"Expected at least 2 columns in {folder_path} on line {line_num}"
                        )
                    parsed_rows.append((line_num, parts))
        except UnicodeDecodeError:
            parsed_rows = []
        if parsed_rows:
            verified_paths = []
            missing_paths = []
            for line_num, parts in parsed_rows:
                sample_name = str(parts[0]).strip()
                if not sample_name:
                    raise ValueError(
                        f"Missing sample name in {folder_path} on line {line_num}"
                    )
                raw_path = parts[1]
                expanded_path = os.path.expandvars(os.path.expanduser(raw_path))
                if not os.path.exists(expanded_path):
                    missing_paths.append(raw_path)
                else:
                    verified_paths.append((expanded_path, sample_name))
            if missing_paths:
                raise FileNotFoundError(
                    "Missing paths in the second column: " + ", ".join(missing_paths)
                )
            for verified, _sample_name in verified_paths:
                print(f"[cellHarmony] Verified input path: {verified}")
            return verified_paths

    # Walk all subdirectories
    for root, dirs, files in os.walk(folder_path):
        # h5 / h5ad
        for f in files:
            if f.endswith(".h5") or f.endswith(".h5ad"):
                results.append(os.path.join(root, f))

        # Check for matrix triplet
        if "matrix.mtx" in files or "matrix.mtx.gz" in files:
            # determine filenames in this folder
            gz = "matrix.mtx.gz" in files
            barcode_file = os.path.join(root, "barcodes.tsv.gz" if gz else "barcodes.tsv")
            feature_options = [
                os.path.join(root, "features.tsv.gz" if gz else "features.tsv"),
                os.path.join(root, "genes.tsv.gz" if gz else "genes.tsv"),
            ]
            features = next((f for f in feature_options if os.path.exists(f)), None)
            if os.path.exists(barcode_file) and features:
                mtx_name = "matrix.mtx.gz" if gz else "matrix.mtx"
                results.append(os.path.join(root, mtx_name))

        # Detect Alevin-Fry output folders
        if "quants_mat.mtx" in files:
            if all(os.path.exists(os.path.join(root, f)) for f in [
                "quants_mat.mtx", "quants_mat_cols.txt", "quants_mat_rows.txt"
            ]):
                results.append(root+"quants_mat.mtx")

        # Check for tarballs
        for f in files:
            if f.endswith((".tar.gz", ".tgz")):
                archive = os.path.join(root, f)
                with tarfile.open(archive, "r:gz") as tar:
                    members = tar.getnames()
                    extract_dirs = set(os.path.dirname(m) for m in members if m.endswith((".mtx.gz", ".tsv.gz")))
                    for subdir in extract_dirs:
                        if (
                            any(m.endswith(f"{subdir}/matrix.mtx.gz") for m in members) and
                            any(m.endswith(f"{subdir}/barcodes.tsv.gz") for m in members) and
                            (any(m.endswith(f"{subdir}/features.tsv.gz") for m in members) or
                             any(m.endswith(f"{subdir}/genes.tsv.gz") for m in members))
                        ):
                            tmpdir = tempfile.mkdtemp(prefix="mtx_")
                            tar.extractall(path=tmpdir)
                            extracted_path = os.path.join(tmpdir, subdir, "matrix.mtx.gz")
                            results.append(extracted_path)
                            break

    return sorted(results)


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
    parser.add_argument('--min_genes', type=int, default=200, help='min_genes for scanpy QC')
    parser.add_argument('--min_cells', type=int, default=3, help='min_cells for scanpy QC')
    parser.add_argument('--min_counts', type=int, default=500, help='min_counts for scanpy QC')
    parser.add_argument('--mit_percent', type=int, default=10, help='mit_percent for scanpy QC')
    parser.add_argument('--generate_umap', action='store_true', help='generate UMAP and marker analysis')
    parser.add_argument('--save_adata', action='store_true', help='save updated AnnData object')
    parser.add_argument('--unsupervised_cluster', action='store_true', help='perform unsupervised clustering analysis')
    parser.add_argument('--append_obs', type=str, default=None, help='Field in .obs to append to the cell barcode (e.g., donor_id)')
    parser.add_argument('--alignment_mode', type=str, default="cosine", help='Alignment mode: "cosine" or "classic"')
    parser.add_argument('--align_cutoff', type=float, default=None, help='Exclude cells with AlignmentScore below this threshold (default: include all)')
    parser.add_argument("--gene_translation",type=str, default=None, help='Optional two-column TSV mapping file (e.g. Ensembl→Symbol) for gene ID translation')
    parser.add_argument('--metacell-align', action='store_true', help='Aggregate cells into metacells prior to alignment')
    parser.add_argument('--metacell-target-size', type=int, default=50, help='Target number of cells per metacell')
    parser.add_argument('--metacell-min-size', type=int, default=25, help='Minimum cells per metacell')
    parser.add_argument('--metacell-max-size', type=int, default=100, help='Maximum cells per metacell')
    parser.add_argument('--metacell-algorithm', type=str, default='kmeans', choices=['kmeans', 'leiden', 'louvain', 'random'], help='Clustering algorithm for metacell construction')
    parser.add_argument('--metacell-neighbors', type=int, default=30, help='Number of neighbors for metacell graph construction')
    parser.add_argument('--metacell-hvg', type=int, default=3000, help='Number of highly variable genes for metacell clustering')
    parser.add_argument('--metacell-pcs', type=int, default=50, help='Number of PCs for metacell clustering')
    parser.add_argument('--metacell-random-count', type=int, default=50, help='Metacells per group when using random aggregation')
    parser.add_argument('--metacell-random-cells', type=int, default=5, help='Cells per metacell for random aggregation')
    parser.add_argument('--metacell-random-replacement', action='store_true', help='Sample cells with replacement for random metacells')
    parser.add_argument('--metacell-random-state', type=int, default=0, help='Random state for metacell construction')
    parser.add_argument('--ambient_correct_cutoff', type=float, default=None, help='Apply SoupX ambient RNA correction with the specified contamination fraction (rho).')
    parser.add_argument('--concat_on_disk', action='store_true', help='Concatenate inputs on disk to reduce peak memory usage')
    parser.add_argument('--concat_batch_size', type=int, default=50, help='Batch size for on-disk concatenation (default: 50)')
    parser.add_argument('--verbose_import', action='store_true', help='Log per-sample import timings during input loading')

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
    align_cutoff = args.align_cutoff
    gene_translation = args.gene_translation
    metacell_align = args.metacell_align
    metacell_target = args.metacell_target_size
    metacell_min = args.metacell_min_size
    metacell_max = args.metacell_max_size
    metacell_algorithm = args.metacell_algorithm
    metacell_neighbors = args.metacell_neighbors
    metacell_hvg = args.metacell_hvg
    metacell_pcs = args.metacell_pcs
    metacell_random_count = args.metacell_random_count
    metacell_random_cells = args.metacell_random_cells
    metacell_random_replacement = args.metacell_random_replacement
    metacell_random_state = args.metacell_random_state
    ambient_correct_cutoff = args.ambient_correct_cutoff
    concat_on_disk = args.concat_on_disk
    concat_batch_size = args.concat_batch_size
    verbose_import = args.verbose_import

    if h5ad_file:
        h5_files = []
    else:
        if not h5_directory:
            print("No --h5dir or --h5ad provided.")
            sys.exit()
        h5_files = get_h5_and_mtx_files(h5_directory)

    if len(h5_files) == 0:
        if not h5ad_file:
            print("No compatible h5, h5ad or .mtx files identified")
            sys.exit()

    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    log_path = os.path.join(log_dir, f"cellHarmony-lite_{timestamp}.log")
    log_file = open(log_path, "w")
    log_file.write("# cellHarmony_lite run parameters\n")
    log_file.write(f"command = {' '.join(sys.argv)}\n")
    for key, value in sorted(vars(args).items()):
        log_file.write(f"{key} = {value}\n")
    log_file.write("\n")
    log_file.flush()

    stdout_orig, stderr_orig = sys.stdout, sys.stderr
    sys.stdout = Tee(stdout_orig, log_file)
    sys.stderr = Tee(stderr_orig, log_file)

    try:
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
            alignment_mode=alignment_mode,
            min_alignment_score=align_cutoff,
            gene_translation_file=gene_translation,
            metacell_align=metacell_align,
            metacell_target_size=metacell_target,
            metacell_min_size=metacell_min,
            metacell_max_size=metacell_max,
            metacell_algorithm=metacell_algorithm,
            metacell_neighbors=metacell_neighbors,
            metacell_hvg=metacell_hvg,
            metacell_pcs=metacell_pcs,
            metacell_random_count=metacell_random_count,
            metacell_random_cells=metacell_random_cells,
            metacell_random_replacement=metacell_random_replacement,
            metacell_random_state=metacell_random_state,
            ambient_correct_cutoff=ambient_correct_cutoff,
            concat_on_disk=concat_on_disk,
            concat_batch_size=concat_batch_size,
            verbose_import=verbose_import
        )
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        sys.stdout = stdout_orig
        sys.stderr = stderr_orig
        log_file.close()
        print(f"[INFO] cellHarmony_lite log written to {log_path}")
