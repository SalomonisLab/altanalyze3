import os
import time
import sys
import shutil
import re
import tempfile
from pathlib import Path
from datetime import datetime
from glob import glob

import pandas as pd
import numpy as np
import scipy.sparse as sp
from scipy.sparse import issparse
import scanpy as sc
import anndata as ad
from tqdm import tqdm
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
    adata.var_names = translated_names
    adata.var_names_make_unique()
    print(f"[INFO] Applied gene translation for {sample_label}: {updated_count} IDs updated")
    return updated_count


def get_h5_and_mtx_files(folder_path):
    """Collect h5/h5ad/mtx inputs from a directory."""
    if folder_path is None:
        return []
    h5_files = glob(os.path.join(folder_path, "*.h5"))
    h5ad_files = glob(os.path.join(folder_path, "*.h5ad"))
    mtx_files = glob(os.path.join(folder_path, "*matrix.mtx*"))
    return sorted(h5_files + h5ad_files + mtx_files)


def combine_and_cluster(
    h5_files,
    h5ad_file=None,
    output_dir="output",
    export_h5ad=False,
    min_genes=500,
    min_cells=5,
    min_counts=1000,
    mit_percent=10,
    generate_umap=True,
    append_obs_field=None,
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
):
    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)

    translation_map = None
    if gene_translation_file is not None and os.path.exists(gene_translation_file):
        translation_map = load_gene_translation(gene_translation_file)

    if h5ad_file is not None:
        adata_combined = sc.read_h5ad(h5ad_file)
        apply_gene_translation(adata_combined, translation_map, os.path.basename(h5ad_file))
        adata_combined.var_names_make_unique()
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
                    adata = sc.read_10x_mtx(tmp_path, var_names='gene_symbols')
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

                    adata = AnnData(X=matrix)
                    adata.obs_names = barcodes_list
                    adata.var_names = gene_symbols
                    adata.var["gene_ids"] = gene_ids

            elif os.path.isdir(path):
                adata = sc.read_10x_mtx(path, var_names='gene_symbols')
                sample_name = os.path.basename(os.path.normpath(path))
            else:
                raise ValueError(f"Unsupported input: {path}")

            apply_gene_translation(adata, translation_map, sample_name)
            adata.var_names_make_unique()

            group_name = sample_name.split('__')[0] if '__' in sample_name else '_'.join(sample_name.split('_')[:-1])
            if sample_name and len(sample_name) > 0:
                adata.obs_names = [f"{bc}.{sample_name}" for bc in adata.obs_names]
            else:
                adata.obs_names = list(adata.obs_names)

            adata.obs["sample"] = sample_name
            adata.obs["group"] = group_name
            adata.obs["Library"] = sample_name
            adata_list.append(adata)

        adata_combined = ad.concat(adata_list, label="sample", join="outer", fill_value=0)
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

    original_cell_adata = None
    metacell_membership = None
    metacell_aligned_adata = None
    cell_to_metacell_map = None

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

    adata_combined.layers["counts"] = adata_combined.X.copy()

    normalize_adata(adata_combined, show_progress=True)
    metacell_aligned_adata = adata_combined.copy() if metacell_align else None

    if export_h5ad:
        adata_combined.write(os.path.join(output_dir, "combined_qc_normalized.h5ad"), compression="gzip")

    adata_unsup = None
    os.chdir(output_dir)
    if metacell_align:
        adata_unsup_metacell = metacell_aligned_adata.copy()
        sc.pp.highly_variable_genes(adata_unsup_metacell, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata_unsup_metacell = adata_unsup_metacell[:, adata_unsup_metacell.var.highly_variable].copy()
        sc.pp.scale(adata_unsup_metacell, max_value=10)
        sc.pp.pca(adata_unsup_metacell, n_comps=50)
        sc.pp.neighbors(adata_unsup_metacell)
        if generate_umap:
            sc.tl.umap(adata_unsup_metacell)
        sc.tl.leiden(adata_unsup_metacell, flavor="leidenalg")
        print('[metacell] unsupervised clustering performed on metacells')
        ensure_category_palette(adata_unsup_metacell, 'leiden')

        pd.DataFrame({
            "metacell_id": adata_unsup_metacell.obs_names,
            "leiden": adata_unsup_metacell.obs['leiden'].astype(str).values
        }).to_csv(os.path.join(output_dir, "metacell_leiden_clusters.tsv"), sep="\t", index=False)

        if cell_to_metacell_map is not None:
            metacell_cluster_map = adata_unsup_metacell.obs['leiden'].astype(str).to_dict()
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
        if generate_umap:
            sc.tl.umap(adata_unsup_cells)
        adata_unsup_cells.obs['leiden'] = adata_combined.obs['scanpy-leiden'].reindex(adata_unsup_cells.obs_names).astype(str)

        if generate_umap and 'X_umap' in adata_unsup_cells.obsm:
            adata_combined.obsm['X_umap_leiden'] = adata_unsup_cells.obsm['X_umap']
            adata_combined.obs['UMAP_leiden_X'] = adata_unsup_cells.obsm['X_umap'][:, 0]
            adata_combined.obs['UMAP_leiden_Y'] = adata_unsup_cells.obsm['X_umap'][:, 1]

        unique_clusters = adata_unsup_cells.obs['leiden'].dropna().unique()
        if len(unique_clusters) > 1:
            sc.tl.rank_genes_groups(adata_unsup_cells, groupby='leiden', method='wilcoxon', use_raw=False)
            save_marker_genes(adata_unsup_cells, 'leiden', os.path.join(output_dir, 'unsupervised_markers.txt'))
            ensure_category_palette(adata_unsup_cells, 'leiden')

            if generate_umap and 'X_umap' in adata_unsup_cells.obsm:
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
        if generate_umap:
            sc.tl.umap(adata_unsup)
        sc.tl.leiden(adata_unsup, flavor="leidenalg", resolution=0.5)
        n_clusters = adata_unsup.obs['leiden'].nunique()
        n_cells = adata_unsup.n_obs
        print(f"[cell] Unsupervised clustering completed: {n_clusters} clusters from {n_cells} cells (resolution=0.5)")
        print('[cell] unsupervised clustering performed on individual cells')

        sc.tl.rank_genes_groups(adata_unsup, groupby='leiden', method='wilcoxon', use_raw=False)
        save_marker_genes(adata_unsup, 'leiden', os.path.join(output_dir, 'unsupervised_markers.txt'))

        ensure_category_palette(adata_unsup, 'leiden')
        if generate_umap and 'X_umap' in adata_unsup.obsm:
            sc.pl.umap(
                adata_unsup,
                color='leiden',
                save="_unsupervised_umap.pdf",
                show=False,
                legend_loc='on data',
                legend_fontsize=5,
                legend_fontweight='normal'
            )

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

    if adata_unsup is not None and 'leiden' in adata_unsup.obs.columns:
        adata_combined.obs['scanpy-leiden'] = adata_unsup.obs['leiden']
        if generate_umap and 'X_umap' in adata_unsup.obsm:
            coords_unsup = adata_unsup.obsm['X_umap']
            adata_combined.obsm['X_umap_leiden'] = coords_unsup
            adata_combined.obs['UMAP_leiden_X'] = coords_unsup[:, 0]
            adata_combined.obs['UMAP_leiden_Y'] = coords_unsup[:, 1]

        leiden_export = pd.DataFrame({
            "CellBarcode": adata_unsup.obs_names,
            "Leiden": adata_unsup.obs['leiden'].values,
        })
        if generate_umap and 'X_umap' in adata_unsup.obsm:
            leiden_export["UMAP_X"] = adata_unsup.obsm['X_umap'][:, 0]
            leiden_export["UMAP_Y"] = adata_unsup.obsm['X_umap'][:, 1]
        leiden_export.to_csv(
            os.path.join(output_dir, "unsupervised_leiden_clusters.tsv"),
            sep="\t", index=False
        )

    adata_combined.write(os.path.join(output_dir, "combined_with_leiden.h5ad"), compression="gzip")

    print(f"Analysis completed in {time.time() - start_time:.2f} seconds.")
    return adata_combined


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run unsupervised Leiden clustering on scRNA-seq inputs')
    parser.add_argument('--h5dir', type=str, default=None, help='path to folder of h5/h5ad/mtx files')
    parser.add_argument('--outdir', type=str, default='output', help='output directory')
    parser.add_argument('--h5ad', type=str, default=None, help='existing h5ad to reuse')
    parser.add_argument('--export_h5ad', action='store_true', help='export an h5ad with all counts and normalized exp')
    parser.add_argument('--min_genes', type=int, default=200, help='min_genes for scanpy QC')
    parser.add_argument('--min_cells', type=int, default=3, help='min_cells for scanpy QC')
    parser.add_argument('--min_counts', type=int, default=500, help='min_counts for scanpy QC')
    parser.add_argument('--mit_percent', type=int, default=10, help='mit_percent for scanpy QC')
    parser.add_argument('--generate_umap', action=argparse.BooleanOptionalAction, default=True,
                        help='generate UMAP and marker analysis (default: True)')
    parser.add_argument('--append_obs', type=str, default=None, help='Field in .obs to append to the cell barcode (e.g., donor_id)')
    parser.add_argument("--gene_translation", type=str, default=None, help='Optional two-column TSV mapping file (e.g. Ensembl→Symbol) for gene ID translation')
    parser.add_argument('--metacell-align', action='store_true', help='Aggregate cells into metacells prior to clustering')
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

    args = parser.parse_args()

    h5_files = get_h5_and_mtx_files(args.h5dir)

    if len(h5_files) == 0 and args.h5ad is None:
        print("No compatible h5, h5ad or .mtx files identified")
        sys.exit()

    log_dir = os.path.join(args.outdir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    log_path = os.path.join(log_dir, f"leiden_cluster_{timestamp}.log")
    log_file = open(log_path, "w")
    log_file.write("# leiden_cluster run parameters\n")
    log_file.write(f"command = {' '.join(sys.argv)}\n")
    for key, value in sorted(vars(args).items()):
        log_file.write(f"{key} = {value}\n")
    log_file.write("\n")
    log_file.flush()

    stdout_orig, stderr_orig = sys.stdout, sys.stderr
    sys.stdout = Tee(stdout_orig, log_file)
    sys.stderr = Tee(stderr_orig, log_file)

    try:
        combine_and_cluster(
            h5_files=h5_files,
            h5ad_file=args.h5ad,
            output_dir=args.outdir,
            export_h5ad=args.export_h5ad,
            min_genes=args.min_genes,
            min_cells=args.min_cells,
            min_counts=args.min_counts,
            mit_percent=args.mit_percent,
            generate_umap=args.generate_umap,
            append_obs_field=args.append_obs,
            gene_translation_file=args.gene_translation,
            metacell_align=args.metacell_align,
            metacell_target_size=args.metacell_target_size,
            metacell_min_size=args.metacell_min_size,
            metacell_max_size=args.metacell_max_size,
            metacell_algorithm=args.metacell_algorithm,
            metacell_neighbors=args.metacell_neighbors,
            metacell_hvg=args.metacell_hvg,
            metacell_pcs=args.metacell_pcs,
            metacell_random_count=args.metacell_random_count,
            metacell_random_cells=args.metacell_random_cells,
            metacell_random_replacement=args.metacell_random_replacement,
            metacell_random_state=args.metacell_random_state,
            ambient_correct_cutoff=args.ambient_correct_cutoff,
        )
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        sys.stdout = stdout_orig
        sys.stderr = stderr_orig
        log_file.close()
        print(f"[INFO] leiden_cluster log written to {log_path}")
