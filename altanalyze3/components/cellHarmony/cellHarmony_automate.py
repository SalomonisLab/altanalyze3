#!/usr/bin/env python3
"""Automate per-sample cell-type annotation for a long-read workflow.

For each sample in a ``sample_dict`` (as returned by
``long_read.isoform_automate.import_metadata``), this:
  1. aggregates the BAM-derived molecule-level h5ad to gene-level counts
     (``aggregate.gene_aggregate.aggregate_to_gene_level`` -> ``<library>-gene.h5ad``), then
  2. assigns each cell barcode to a reference population with cellHarmony-lite
     (``run_cellHarmony_lite.run_cellharmony_lite``, default human bone-marrow centroid reference).

It returns ``barcode_cluster_dirs`` -- the list of per-sample ``barcode.sample<tab>cluster`` TSVs
consumed by ``isoform_matrix.return_cluster_order`` / ``import_barcode_clusters`` -- so the master
workflow script reduces to a single call.
"""

from __future__ import annotations

import os

from ..aggregate.gene_aggregate import aggregate_to_gene_level
from .run_cellHarmony_lite import run_cellharmony_lite, DEFAULT_BM_REFERENCE


def annotate_samples(sample_dict, cellharmony_ref=DEFAULT_BM_REFERENCE,
                     gene_translation_file=None, output_dir=None, log=print, logger=None):
    """Aggregate to gene level and run cellHarmony-lite for every sample.

    Parameters
    ----------
    sample_dict : dict
        ``{uid: [ {library, matrix, ...}, ... ]}`` from ``import_metadata``. ``matrix`` is the
        BAM-derived molecule-level h5ad path; ``library`` names the sample.
    cellharmony_ref : str
        Centroid/states reference TSV (default: human bone marrow Hs-MarrowAtlas-L3M.txt).
    gene_translation_file : str, optional
        2-column Ensembl->Symbol table (BAM gene ids are Ensembl; the BM reference uses symbols).
    output_dir : str, optional
        Root for cellHarmony outputs; each sample writes to ``<output_dir>/<library>/``. Defaults to
        a ``cellHarmony`` directory under the current working directory.
    logger : WorkflowLogger, optional
        When given, records a per-sample time/memory/multiprocessing entry for the gene-aggregation
        and cellHarmony sub-steps. Both sub-steps are single-process. ``logger.log`` is used for
        progress messages so they land in the workflow log.

    Returns
    -------
    list[str]
        ``barcode_cluster_dirs`` -- one ``barcode.sample<tab>cluster`` TSV per sample.
    """
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), "cellHarmony")
    _log = logger.log if logger is not None else log

    barcode_cluster_dirs = []
    for uid, samples in sample_dict.items():
        for s in samples:
            library = s["library"]
            _log(f"[cellHarmony-automate] {uid} / {library}")
            if logger is not None:
                with logger.step("2_gene_aggregation", sample=library, multiprocessing=False):
                    gene_h5ad = aggregate_to_gene_level(s["matrix"], log=_log)
                with logger.step("3_cellHarmony_lite", sample=library, multiprocessing=False,
                                 note=os.path.basename(cellharmony_ref)):
                    bc_clusters = run_cellharmony_lite(
                        gene_h5ad, output_dir=os.path.join(output_dir, library),
                        cellharmony_ref=cellharmony_ref,
                        gene_translation_file=gene_translation_file,
                        sample_name=library, log=_log,
                    )
            else:
                gene_h5ad = aggregate_to_gene_level(s["matrix"], log=_log)
                bc_clusters = run_cellharmony_lite(
                    gene_h5ad, output_dir=os.path.join(output_dir, library),
                    cellharmony_ref=cellharmony_ref,
                    gene_translation_file=gene_translation_file,
                    sample_name=library, log=_log,
                )
            barcode_cluster_dirs.append(bc_clusters)
    return barcode_cluster_dirs
