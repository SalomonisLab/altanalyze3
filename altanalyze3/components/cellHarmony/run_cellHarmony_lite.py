#!/usr/bin/env python3
"""Compact local wrapper around cellHarmony-lite -- the same alignment cellHarmony-web runs.

Given a gene-level h5ad (cells x genes) and a centroid reference (the human bone-marrow
``Hs-MarrowAtlas-L3M.txt`` states file that the web uses by default), this assigns each cell
barcode to a reference population and writes a 2-column ``barcode<tab>cluster`` TSV -- the exact
format consumed downstream by ``isoform_matrix.return_cluster_order`` /
``isoform_matrix.import_barcode_clusters``.

It mirrors the web's ``flask/pipeline.py`` call to ``cellHarmony_lite.combine_and_align_h5`` with the
web's default QC/alignment options. The one addition is ``gene_translation_file``: BAM-derived gene
ids are Ensembl IDs while the BM reference is keyed by gene symbol, so a 2-column Ensembl->Symbol
table (e.g. ``Hs_Ensembl-annotations.txt``) is passed through to translate.
"""

from __future__ import annotations

import os

import pandas as pd

from . import cellHarmony_lite

# Default human bone-marrow centroid reference shipped with cellHarmony-web (hs_bm_reference).
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_BM_REFERENCE = os.path.join(
    _HERE, "flask", "references", "Human", "BoneMarrow", "Zhang-2024", "Hs-MarrowAtlas-L3M.txt"
)


def _normalize_translation_file(translation_file, output_dir):
    """cellHarmony_lite.load_gene_translation expects EXACTLY two columns (source<tab>target).

    Gene-annotation tables like Hs_Ensembl-annotations.txt carry extra columns (description, ...);
    read with only two ``names`` they collapse into a MultiIndex and yield 0 mappings. Write a clean
    2-column (Ensembl-id, symbol) copy and return its path. Returns the original path unchanged when
    it is already a usable 2-column table.
    """
    if not translation_file or not os.path.exists(translation_file):
        return translation_file
    df = pd.read_csv(translation_file, sep="\t", header=None, dtype=str)
    if df.shape[1] <= 2:
        return translation_file
    two = df.iloc[:, :2].dropna()
    two = two[two.iloc[:, 1].astype(str).str.len() > 0]
    out = os.path.join(output_dir, "_ensembl_to_symbol.txt")
    two.to_csv(out, sep="\t", header=False, index=False)
    return out


def run_cellharmony_lite(
    gene_h5ad_path,
    output_dir,
    cellharmony_ref=DEFAULT_BM_REFERENCE,
    gene_translation_file=None,
    barcode_cluster_out=None,
    sample_name=None,
    # web defaults (flask/pipeline.py)
    min_genes=200,
    min_cells=3,
    min_counts=500,
    mit_percent=10,
    alignment_mode="cosine",
    min_alignment_score=0.1,
    log=print,
):
    """Align one gene-level h5ad to the centroid reference; return the barcode->cluster TSV path.

    Parameters
    ----------
    gene_h5ad_path : str
        Gene-level h5ad (cells x genes) -- e.g. the output of aggregate.gene_aggregate.
    output_dir : str
        cellHarmony-lite working/output dir (writes ``cellHarmony_lite_assignments.txt`` here).
    cellharmony_ref : str
        Centroid/states reference TSV (default: human bone marrow Hs-MarrowAtlas-L3M.txt).
    gene_translation_file : str, optional
        2-column Ensembl->Symbol table; required when the h5ad var_names are Ensembl ids.
    barcode_cluster_out : str, optional
        Where to write the 2-column barcode<tab>cluster TSV. Defaults to
        ``<output_dir>/<library>_barcode_clusters.txt``.

    Returns
    -------
    str
        Path to the 2-column barcode->cluster TSV (input to return_cluster_order).
    """
    os.makedirs(output_dir, exist_ok=True)
    log(f"[cellHarmony] aligning {os.path.basename(gene_h5ad_path)} to "
        f"{os.path.basename(cellharmony_ref)}")

    # Ensure the Ensembl->Symbol table is a clean 2-column file (see _normalize_translation_file).
    gene_translation_file = _normalize_translation_file(gene_translation_file, output_dir)

    cellHarmony_lite.combine_and_align_h5(
        h5_files=[],
        h5ad_file=gene_h5ad_path,
        cellharmony_ref=cellharmony_ref,
        output_dir=output_dir,
        export_cptt=False,
        export_h5ad=False,
        min_genes=min_genes,
        min_cells=min_cells,
        min_counts=min_counts,
        mit_percent=mit_percent,
        generate_umap=False,
        save_adata=False,
        unsupervised_cluster=False,
        alignment_mode=alignment_mode,
        min_alignment_score=min_alignment_score,
        gene_translation_file=gene_translation_file,
        metacell_align=False,
        return_adata=False,
    )

    assignments_path = os.path.join(output_dir, "cellHarmony_lite_assignments.txt")
    if not os.path.exists(assignments_path):
        raise FileNotFoundError(f"cellHarmony-lite assignments missing: {assignments_path}")

    # reduce CellBarcode,<ref_name>,AlignmentScore  ->  "barcode.sample"<tab>cluster
    # The downstream isoform_matrix.import_barcode_clusters splits the first column on '.' to recover
    # (barcode, sample_name), so the barcode must be namespaced as "<barcode>.<sample_name>".
    df = pd.read_csv(assignments_path, sep="\t")
    barcode_col = df.columns[0]
    cluster_col = df.columns[1]  # the reference-population column
    base = os.path.basename(gene_h5ad_path).split(".h5ad")[0].replace("-gene", "")
    sample = sample_name or base
    out2 = pd.DataFrame({
        "barcode_cluster": df[barcode_col].astype(str) + "." + sample,
        "cluster": df[cluster_col].astype(str),
    })

    if barcode_cluster_out is None:
        barcode_cluster_out = os.path.join(output_dir, f"{base}_barcode_clusters.txt")
    out2.to_csv(barcode_cluster_out, sep="\t", header=False, index=False)
    log(f"[cellHarmony] wrote {len(out2):,} barcode->cluster rows -> {barcode_cluster_out}")
    return barcode_cluster_out
