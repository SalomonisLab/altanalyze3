#!/usr/bin/env python3
"""Post-pruning annotation for scTriangulate: MarkerFinder, cell-state naming, HOPACH order.

``ScTriangulate.lazy_run`` stops at ``obs['pruned']``. Those labels carry the source annotation
and the source cluster id, for example ``ICGS3_RNA@C3``, and nothing about the biology. Every
downstream reader then has to repeat the same three steps by hand. This module runs them, and
``cli.py`` runs it by default.

Step 4  MarkerFinder on the pruned clusters.
        Calls altanalyze3.components.visualization.marker_heatmap_h5ad
        .generate_marker_heatmap_from_adata, the same function the marker-heatmap CLI calls, with
        the same arguments. It writes the marker table, the redundant-marker table, the fold
        matrix, the centroid matrix and the heatmap.

Step 5  Cell-state naming by gene-set enrichment.
        Calls altanalyze3.components.clustering.ICGS.biomarker_enrichment, the same function ICGS3
        runs, on the step-4 markers against the GO-Elite BioMarkers cell-state sets, then
        ICGS.clean_biomarker_prediction_labels for the label cleanup. Each cluster keeps its
        source id, so a name reads ``T cell (ICGS3_RNA@C11)``.

Step 6  HOPACH on the centroids, then MarkerFinder again.
        Calls altanalyze3.components.clustering.hopach.hopach on the step-4 centroid matrix, with
        the clusters as objects. The order goes into ``uns['lineage_order']``, which
        marker_heatmap_h5ad._resolve_cluster_order reads, so the second heatmap draws its blocks
        in HOPACH order. When the object also carries antibody features, the module writes a
        matching antibody-only heatmap over the same cells, states and order.

Nothing here reimplements a method. Every step calls the validated altanalyze3 function.
"""

from __future__ import annotations

import json
import os
import re
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

__all__ = ["annotate_pruned_clusters", "DEFAULT_ADT_PREFIX"]

DEFAULT_ADT_PREFIX = "AB_"

#: Labels a curated annotation uses to mean "no call". A cluster whose dominant --annotate-lead
#: label is one of these takes the enriched cell-state term instead, because "unassigned" names
#: nothing. QC calls such as "Doublet" and "Dying" are NOT here: they do name the cluster.
UNINFORMATIVE_LEAD_LABELS = frozenset({
    "", "na", "n/a", "nan", "none", "null", "unassigned", "unassign", "unknown",
    "unlabeled", "unlabelled", "unclassified", "uncharacterized", "undetermined", "other",
})

_SPECIES_TO_ICGS = {
    "human": "Hs", "hs": "Hs", "homo sapiens": "Hs", "Hs": "Hs",
    "mouse": "Mm", "mm": "Mm", "mus musculus": "Mm", "Mm": "Mm",
}


def _icgs_species(species: str) -> str:
    key = str(species).strip()
    resolved = _SPECIES_TO_ICGS.get(key) or _SPECIES_TO_ICGS.get(key.lower())
    if resolved is None:
        raise ValueError(
            f"Unsupported species '{species}'. GO-Elite BioMarkers ship for human and mouse only."
        )
    return resolved


def _sanitize(label: str) -> str:
    """Collapse whitespace and drop the characters that break a TSV or a figure label."""
    text = re.sub(r"\s+", " ", str(label)).strip()
    return text.replace("\t", " ").replace("\n", " ")


def _predictions_from_enrichment(
    markers: pd.DataFrame,
    background: Sequence[str],
    clusters: Sequence[str],
    species: str,
    outdir: str,
    biomarker_file: Optional[str],
    max_fdr: float,
    min_overlap: int,
    log,
) -> Dict[str, str]:
    """Cluster -> cell-state label, from the GO-Elite BioMarkers enrichment ICGS3 uses.

    ``biomarker_enrichment`` returns the single best term per cluster whether or not that term is
    significant. Naming a cluster from a weak term invents a cell type: a mouse thymus run
    produced "Retina Rheaume et al." at FDR 0.027 and "Marrow Oligodendrocyte" at FDR 0.049.
    Two gates drop those, and the caller then names the cluster from its own strongest marker.

    Both defaults come from a measurement, not from taste. Across the 79 clusters of three
    scTriangulate runs on one mouse thymus dataset, the best-term FDR ran from 2.1e-71 to 5.3e-02.
    Every term that FDR 0.05 accepted and FDR 1e-5 rejected -- 18 of them -- overlapped the term
    set by only 2 to 8 of about 36 to 60 markers. Every term a reader recognised overlapped by 10
    to 36. FDR 1e-5 keeps 60 of 79 clusters (76%); 0.05 keeps 78 of 79 (99%), including the two
    named above. ``min_overlap`` guards the same failure independently, because an overlap of 2
    genes names nothing however small its FDR.

    Returns an empty mapping when no BioMarkers file exists for the species, which lets the caller
    fall back instead of failing the whole run.
    """
    from ..clustering.ICGS import (
        ICGS3Config,
        biomarker_enrichment,
        clean_biomarker_prediction_labels,
        _default_biomarker_file,
    )

    resolved = biomarker_file or _default_biomarker_file(_icgs_species(species))
    if not resolved or not os.path.exists(resolved):
        log(f"[annotate] no GO-Elite BioMarkers file for species '{species}'; "
            f"falling back to the top marker gene for every name")
        return {}
    log(f"[annotate] gene-set enrichment against {resolved}")

    query = markers.loc[~markers["_is_adt"], ["Gene", "cluster"]].rename(
        columns={"Gene": "marker", "cluster": "top_cluster"}
    )
    if query.empty:
        log("[annotate] no non-antibody markers to enrich; falling back to the top marker feature")
        return {}

    config = ICGS3Config(input_paths=["<in-memory>"], output_dir=outdir,
                         species=_icgs_species(species), biomarker_file=resolved)
    labels = biomarker_enrichment(query, list(background), config, outdir)
    if labels is None or labels.empty:
        log("[annotate] enrichment returned no term; falling back to the top marker feature")
        return {}

    # clean_biomarker_prediction_labels appends "_c<cluster>"; the cluster id is carried in the
    # composed name already, so strip the suffix it just added rather than reimplement its rules.
    labels = clean_biomarker_prediction_labels(labels)
    out: Dict[str, str] = {}
    rejected: List[str] = []
    for _, row in labels.iterrows():
        cluster = str(row["cluster"])
        fdr = float(row["fdr"]) if "fdr" in labels.columns and pd.notna(row.get("fdr")) else 1.0
        overlap = int(row["overlap"]) if "overlap" in labels.columns and pd.notna(row.get("overlap")) else 0
        if fdr > max_fdr or overlap < min_overlap:
            why = []
            if fdr > max_fdr:
                why.append(f"FDR {fdr:.2g}")
            if overlap < min_overlap:
                why.append(f"overlap {overlap}")
            rejected.append(f"{cluster} ({', '.join(why)})")
            continue
        prediction = str(row["cell_type_prediction"])
        suffix = f"_c{cluster}"
        if prediction.endswith(suffix):
            prediction = prediction[: -len(suffix)]
        prediction = _sanitize(prediction)
        if prediction:
            out[cluster] = prediction
    if rejected:
        log(f"[annotate] {len(rejected)} of {len(clusters)} clusters have no term at FDR <= "
            f"{max_fdr:g}, so they take a marker-derived name: "
            f"{', '.join(rejected[:6])}{' ...' if len(rejected) > 6 else ''}")
    missing = [c for c in clusters if c not in out and not any(r.startswith(c + " ") for r in rejected)]
    if missing:
        log(f"[annotate] {len(missing)} of {len(clusters)} clusters got no enriched term at all: "
            f"{', '.join(map(str, missing[:6]))}{' ...' if len(missing) > 6 else ''}")
    return out


def _benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Return BH-adjusted p-values for a 1-D array. Order and length are preserved."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * n / (np.arange(n) + 1.0)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(ranked, 1.0)
    return out


def lead_enrichment_table(
    labels_in: pd.Series,
    obs: pd.DataFrame,
    lead_columns: Sequence[str],
    *,
    clusters: Sequence[str],
) -> pd.DataFrame:
    """Score every (cluster, reference label) barcode overlap with a hypergeometric test.

    One test per cluster and reference label. ``a`` is the number of barcodes the cluster and
    the reference label share. The null draws ``len(cluster)`` barcodes from the ``N`` annotated
    barcodes without replacement, so the right-tail probability of seeing ``a`` or more is
    ``hypergeom.sf(a - 1, N, n_label, n_cluster)``. ``fold_enrichment`` is the observed share of
    the cluster over the reference label's share of all barcodes. Benjamini-Hochberg runs once
    over every test of every reference column, so the two references compete on one scale.
    """
    from scipy.stats import hypergeom

    rows = []
    for column in lead_columns:
        values = obs[column].astype(str)
        annotated = ~values.str.strip().str.lower().isin(UNINFORMATIVE_LEAD_LABELS)
        total = int(annotated.sum())
        if total == 0:
            continue
        label_counts = values[annotated].value_counts()
        for cluster in clusters:
            in_cluster = (labels_in == cluster) & annotated
            n_cluster = int(in_cluster.sum())
            if n_cluster == 0:
                continue
            observed = values[in_cluster].value_counts()
            for label, overlap in observed.items():
                n_label = int(label_counts[label])
                overlap = int(overlap)
                pvalue = float(hypergeom.sf(overlap - 1, total, n_label, n_cluster))
                expected = n_cluster * n_label / total
                rows.append({
                    "cluster": cluster,
                    "reference": column,
                    "reference_label": label,
                    "overlap": overlap,
                    "n_cluster": n_cluster,
                    "n_reference_label": n_label,
                    "n_annotated": total,
                    "expected_overlap": expected,
                    "fold_enrichment": (overlap / expected) if expected > 0 else np.inf,
                    "fraction_of_cluster": overlap / n_cluster,
                    "pvalue": pvalue,
                })
    table = pd.DataFrame(rows)
    if table.empty:
        return table
    table["FDR"] = _benjamini_hochberg(table["pvalue"].to_numpy())
    return table.sort_values(
        ["cluster", "FDR", "fold_enrichment", "overlap"],
        ascending=[True, True, False, False], kind="mergesort",
    ).reset_index(drop=True)


def _predictions_from_lead_enrichment(
    table: pd.DataFrame,
    clusters: Sequence[str],
    *,
    max_fdr: float,
    min_overlap: int,
    log,
) -> Dict[str, Dict[str, object]]:
    """Pick, per cluster, the single best-enriched reference label across every reference."""
    picks: Dict[str, Dict[str, object]] = {}
    if table.empty:
        return picks
    rejected = 0
    for cluster in clusters:
        sub = table[table["cluster"] == cluster]
        if sub.empty:
            continue
        best = sub.iloc[0]
        if float(best["FDR"]) > max_fdr or int(best["overlap"]) < min_overlap:
            rejected += 1
            continue
        picks[str(cluster)] = {
            "label": _sanitize(str(best["reference_label"])),
            "reference": str(best["reference"]),
            "FDR": float(best["FDR"]),
            "fold_enrichment": float(best["fold_enrichment"]),
            "overlap": int(best["overlap"]),
            "fraction_of_cluster": float(best["fraction_of_cluster"]),
        }
    if rejected:
        log(f"[annotate] {rejected} of {len(clusters)} clusters reach no reference label at "
            f"FDR <= {max_fdr:g} and overlap >= {min_overlap}; those take the enriched term instead")
    by_reference = pd.Series([v["reference"] for v in picks.values()], dtype=object).value_counts().to_dict()
    if by_reference:
        log("[annotate] winning reference per cluster: "
            + ", ".join(f"{k}={v}" for k, v in sorted(by_reference.items())))
    return picks


def _resolve_marker_validate_scaling(adata, adt_prefix: str, layer, use_raw: bool, log) -> bool:
    """Decide whether MarkerFinder's whole-matrix depth-normalization guard can run.

    MarkerFinder correlates each feature against a 0/1 cluster indicator, a statistic that per-cell
    sequencing depth confounds, so ``detect_input_scaling`` refuses any matrix whose cell totals
    are not constant after inverting a log. A CITE-seq object fails that test for a reason the
    guard does not cover: the antibody block carries no sequencing depth and is normalized on its
    own scale, so appending it moves each cell total by the cell's antibody load. Measured on this
    thymus object, 195 antibody columns beside 32,346 RNA columns lift the coefficient of variation
    of the cell totals from about 0 to 0.017, seventeen times the 1e-3 tolerance.

    So the guard is switched off only when the RNA block ALONE passes it. When the RNA block fails,
    the guard stays on and MarkerFinder refuses the matrix, which is the correct outcome.
    """
    from ..cellHarmony.markerFinder import detect_input_scaling, _format_scaling_report

    names = np.asarray([str(v) for v in adata.var_names])
    is_adt = np.char.startswith(names, adt_prefix)
    if not is_adt.any():
        return True

    if use_raw and adata.raw is not None:
        matrix = adata.raw.X
    elif layer:
        matrix = adata.layers[layer]
    else:
        matrix = adata.X
    rna_report = detect_input_scaling(matrix[:, ~is_adt])
    full_report = detect_input_scaling(matrix)
    if rna_report["status"] == "ok" and full_report["status"] != "ok":
        log(f"[annotate] MarkerFinder scaling guard: the {int(is_adt.sum())} '{adt_prefix}' "
            f"features put the whole matrix at {_format_scaling_report(full_report)}, but the "
            f"{int((~is_adt).sum())} RNA features alone are {rna_report['transform']} "
            f"(cv {rna_report['cv_cell_total']:.2g}). The antibody block carries no sequencing "
            f"depth, so the whole-matrix check is switched off for this run.")
        return False
    if rna_report["status"] != "ok":
        log(f"[annotate] MarkerFinder scaling guard stays ON: the RNA features are "
            f"{rna_report['transform']} ({_format_scaling_report(rna_report)})")
    return True


def _fallback_label(cluster: str, markers: pd.DataFrame, adt_prefix: str) -> str:
    """Name a cluster from its strongest marker when enrichment gives nothing."""
    grp = markers[markers["cluster"] == cluster]
    if grp.empty:
        return "Unresolved"
    top = grp.sort_values("Fold", ascending=False)["Gene"].astype(str).iloc[0]
    if top.startswith(adt_prefix):
        top = top[len(adt_prefix):]
        top = re.sub(r"_TotalA$", "", top)
    return _sanitize(top) + "-high"


def annotate_pruned_clusters(
    adata,
    outdir: str,
    *,
    cluster_key: str = "pruned",
    species: str = "human",
    top_n: int = 60,
    cells_per_cluster: int = 100,
    seed: int = 0,
    layer: Optional[str] = None,
    use_raw: bool = False,
    covariate_columns: Optional[Sequence[str]] = None,
    lead_annotation: Optional[str] = None,
    lead_mode: str = "dominant",
    lead_max_fdr: float = 0.05,
    lead_min_overlap: int = 10,
    biomarker_file: Optional[str] = None,
    enrichment_max_fdr: float = 1e-5,
    enrichment_min_overlap: int = 5,
    adt_prefix: str = DEFAULT_ADT_PREFIX,
    hopach_distance: str = "cor",
    hopach_kmax: int = 9,
    hopach_kmin: int = 2,
    hopach_mincluster: int = 5,
    hopach_random_state: int = 0,
    name_key: str = "cluster_name",
    hopach_key: str = "hopach_cluster",
    log=print,
) -> Dict[str, object]:
    """Run steps 4 to 6 on a pruned scTriangulate object. Modifies ``adata`` in place.

    Parameters
    ----------
    adata : AnnData
        A pruned scTriangulate object. ``obs[cluster_key]`` must hold the pruned labels.
    outdir : str
        Where the MarkerFinder, GO-Elite and HOPACH folders go.
    lead_annotation : str, optional
        An obs column holding a trusted annotation, for example a published cell type. When given,
        the dominant label of that column names the cluster and the enriched term becomes the
        fallback. Use it when a curated annotation competes in ``--query``. A dominant label in
        ``UNINFORMATIVE_LEAD_LABELS``, such as "unassigned", never names a cluster.
    enrichment_max_fdr : float
        Reject an enriched cell-state term at a worse FDR than this. Default 1e-5, measured; see
        ``_predictions_from_enrichment``.
    enrichment_min_overlap : int
        Reject a term that overlaps the cluster markers by fewer genes than this. Default 5.

    Returns
    -------
    dict
        Paths written, the name mapping, the HOPACH order and the HOPACH diagnostics.
    """
    from ..clustering.hopach import hopach
    from ..visualization.marker_heatmap_h5ad import generate_marker_heatmap_from_adata

    outdir = os.path.abspath(outdir)
    marker_dir = os.path.join(outdir, "MarkerFinder", "pruned")
    ordered_dir = os.path.join(outdir, "MarkerFinder", "hopach_ordered")
    adt_dir = os.path.join(outdir, "MarkerFinder", "hopach_ordered_ADT")
    hopach_dir = os.path.join(outdir, "hopach")
    for path in (marker_dir, ordered_dir, hopach_dir):
        os.makedirs(path, exist_ok=True)

    if cluster_key not in adata.obs:
        raise KeyError(
            f"annotate_pruned_clusters needs obs['{cluster_key}']. "
            f"Available: {', '.join(map(str, adata.obs.columns))}"
        )
    labels_in = adata.obs[cluster_key].astype(str)
    clusters_in = sorted(labels_in.unique())
    n_cells = int(adata.n_obs)
    log(f"[annotate] step 4: MarkerFinder on {len(clusters_in)} '{cluster_key}' clusters "
        f"across {n_cells} cells, top {top_n} per cluster, {cells_per_cluster} cells shown")

    covariates = [c for c in (covariate_columns or []) if c in adata.obs.columns]
    validate_scaling = _resolve_marker_validate_scaling(adata, adt_prefix, layer, use_raw, log)
    step4 = generate_marker_heatmap_from_adata(
        adata,
        cluster_key=cluster_key,
        out=os.path.join(marker_dir, "sctri_pruned_marker_heatmap.pdf"),
        top_n=top_n,
        marker_method="markerfinder",
        use_raw=use_raw,
        layer=layer,
        cells_per_cluster=cells_per_cluster,
        seed=seed,
        species=species,
        covariate_columns=covariates,
        validate_scaling=validate_scaling,
    )

    markers = pd.read_csv(step4["markers_tsv"], sep="\t")
    if markers.empty:
        raise ValueError(f"MarkerFinder wrote no marker to {step4['markers_tsv']}")
    markers["_is_adt"] = markers["Gene"].astype(str).str.startswith(adt_prefix)
    per_cluster = markers.groupby("cluster").size()
    short = per_cluster[per_cluster < top_n]
    if len(short):
        log(f"[annotate] {len(short)} of {markers['cluster'].nunique()} clusters own fewer than "
            f"{top_n} unique markers: "
            + ", ".join(f"{k}={v}" for k, v in short.sort_values().items()))

    background = [str(v) for v in adata.var_names if not str(v).startswith(adt_prefix)]
    log(f"[annotate] step 5: cell-state naming; enrichment background = {len(background)} features "
        f"({int(adata.n_vars) - len(background)} antibody features excluded)")
    predictions = _predictions_from_enrichment(
        markers, background, clusters_in, species, outdir, biomarker_file,
        enrichment_max_fdr, enrichment_min_overlap, log
    )

    lead: Dict[str, str] = {}
    lead_detail: Dict[str, Dict[str, object]] = {}
    lead_columns: List[str] = []
    lead_table_path: Optional[str] = None
    if lead_annotation:
        lead_columns = [c.strip() for c in str(lead_annotation).split(",") if c.strip()]
        missing = [c for c in lead_columns if c not in adata.obs.columns]
        if missing:
            raise KeyError(
                f"lead_annotation {missing} are not obs columns. "
                f"Available: {', '.join(map(str, adata.obs.columns))}"
            )
        if lead_mode not in ("dominant", "enrichment"):
            raise ValueError(f"lead_mode must be 'dominant' or 'enrichment'; got '{lead_mode}'")
        if lead_mode == "dominant":
            if len(lead_columns) != 1:
                raise ValueError(
                    "lead_mode 'dominant' takes exactly one obs column; "
                    f"got {len(lead_columns)}. Use lead_mode 'enrichment' to let several compete."
                )
            column = lead_columns[0]
            log(f"[annotate] naming led by the dominant obs['{column}'] label; "
                "the enriched term is the fallback")
            skipped = 0
            for cluster in clusters_in:
                sub = adata.obs.loc[labels_in == cluster, column].astype(str)
                if not len(sub):
                    continue
                label = _sanitize(sub.value_counts().index[0])
                if label.strip().lower() in UNINFORMATIVE_LEAD_LABELS:
                    skipped += 1
                    continue
                lead[cluster] = label
                lead_detail[cluster] = {"label": label, "reference": column, "mode": "dominant"}
            if skipped:
                log(f"[annotate] {skipped} of {len(clusters_in)} clusters have an uninformative "
                    f"dominant '{column}' label; those take the enriched term instead")
        else:
            log(f"[annotate] naming led by barcode-overlap enrichment over "
                f"{', '.join(lead_columns)}; hypergeometric, BH over every test, "
                f"accept at FDR <= {lead_max_fdr:g} and overlap >= {lead_min_overlap}")
            lead_table = lead_enrichment_table(
                labels_in, adata.obs, lead_columns, clusters=clusters_in)
            lead_table_path = os.path.join(hopach_dir, "lead_reference_enrichment.tsv")
            lead_table.to_csv(lead_table_path, sep="\t", index=False, float_format="%.6g")
            log(f"[annotate] {len(lead_table)} cluster-by-reference-label tests -> {lead_table_path}")
            picks = _predictions_from_lead_enrichment(
                lead_table, clusters_in,
                max_fdr=lead_max_fdr, min_overlap=lead_min_overlap, log=log)
            for cluster, detail in picks.items():
                lead[cluster] = str(detail["label"])
                lead_detail[cluster] = dict(detail, mode="enrichment")

    names: Dict[str, str] = {}
    source: Dict[str, str] = {}
    for cluster in clusters_in:
        label = lead.get(cluster) or predictions.get(cluster)
        source[cluster] = ("lead" if lead.get(cluster)
                           else ("enrichment" if predictions.get(cluster) else "top-marker"))
        if not label:
            label = _fallback_label(cluster, markers, adt_prefix)
        names[cluster] = f"{label} ({cluster})"
    if len(set(names.values())) != len(names):
        # Two clusters can enrich for the same term. The source id keeps every name unique, so a
        # collision here means the id was lost, which would silently merge two clusters.
        raise ValueError("composed cluster names are not unique; the source id was lost")

    log(f"[annotate] step 6: HOPACH on the {len(clusters_in)} centroids, d='{hopach_distance}'")
    centroids = pd.read_csv(step4["centroids_tsv"], sep="\t", index_col=0)
    if centroids.isna().any().any():
        raise ValueError(f"{int(centroids.isna().sum().sum())} NaN in {step4['centroids_tsv']}")
    cols = [str(c) for c in centroids.columns]
    unknown = [c for c in cols if c not in names]
    if unknown:
        raise ValueError(f"centroid columns absent from obs['{cluster_key}']: {unknown}")
    result = hopach(
        centroids.T.to_numpy(dtype=float),
        d=hopach_distance, kmax=hopach_kmax, kmin=hopach_kmin,
        mincluster=hopach_mincluster, collapse=True, ordering="medoid",
        random_state=hopach_random_state,
    )
    order_idx = np.asarray(result.order)
    if sorted(order_idx.tolist()) != list(range(len(cols))):
        raise ValueError("HOPACH order is not a permutation of the centroid columns")
    hop_labels = np.asarray(result.clust.labels)
    ordered = [cols[i] for i in order_idx]
    remap: Dict[int, str] = {}
    for i in order_idx:
        if hop_labels[i] not in remap:
            remap[hop_labels[i]] = f"H{len(remap) + 1}"
    groups = {c: remap[hop_labels[cols.index(c)]] for c in cols}
    log(f"[annotate] HOPACH: k={result.clust.k}, MSS={result.clust.mss:.4f}, "
        f"sizes={result.clust.sizes.tolist()}")

    lineage = [names[c] for c in ordered]
    table = pd.DataFrame({
        "hopach_order": np.arange(1, len(ordered) + 1),
        "sctriangulate_cluster": ordered,
        "cluster_name": [names[c].rsplit(" (", 1)[0] for c in ordered],
        "named_label": lineage,
        "hopach_cluster": [groups[c] for c in ordered],
        "name_source": [source[c] for c in ordered],
        "n_cells": [int((labels_in == c).sum()) for c in ordered],
    })
    if lead_detail:
        table["lead_reference"] = [lead_detail.get(c, {}).get("reference", "") for c in ordered]
        table["lead_FDR"] = [lead_detail.get(c, {}).get("FDR", np.nan) for c in ordered]
        table["lead_fold_enrichment"] = [
            lead_detail.get(c, {}).get("fold_enrichment", np.nan) for c in ordered]
        table["lead_overlap"] = [lead_detail.get(c, {}).get("overlap", np.nan) for c in ordered]
        table["lead_fraction_of_cluster"] = [
            lead_detail.get(c, {}).get("fraction_of_cluster", np.nan) for c in ordered]
    if int(table["n_cells"].sum()) != n_cells:
        raise ValueError(f"cluster sizes sum to {int(table['n_cells'].sum())}, expected {n_cells}")
    table_path = os.path.join(hopach_dir, "hopach_centroid_clusters.tsv")
    table.to_csv(table_path, sep="\t", index=False)
    ordered_centroids_path = os.path.join(hopach_dir, "hopach_ordered_centroids.tsv")
    centroids[ordered].rename(columns=names).to_csv(
        ordered_centroids_path, sep="\t", float_format="%.4g"
    )

    adata.obs[name_key] = pd.Categorical(labels_in.map(names), categories=lineage, ordered=True)
    adata.obs[hopach_key] = pd.Categorical(
        labels_in.map(groups), categories=[f"H{i + 1}" for i in range(result.clust.k)], ordered=True
    )
    if adata.obs[name_key].isna().any():
        raise ValueError(f"{int(adata.obs[name_key].isna().sum())} cells received no name")
    adata.uns["lineage_order"] = lineage
    adata.uns["sctriangulate_annotation"] = json.dumps({
        "cluster_key": cluster_key, "species": species, "top_n": top_n,
        "cells_per_cluster": cells_per_cluster, "lead_annotation": lead_annotation,
        "lead_mode": lead_mode, "lead_columns": lead_columns,
        "lead_max_fdr": lead_max_fdr, "lead_min_overlap": lead_min_overlap,
        "lead_table": lead_table_path,
        "enrichment_max_fdr": enrichment_max_fdr,
        "enrichment_min_overlap": enrichment_min_overlap,
        "hopach": {"distance": hopach_distance, "kmax": hopach_kmax, "kmin": hopach_kmin,
                   "mincluster": hopach_mincluster, "k": int(result.clust.k),
                   "mss": float(result.clust.mss), "sizes": result.clust.sizes.tolist()},
        "name_source_counts": pd.Series(source).value_counts().to_dict(),
    })

    log(f"[annotate] step 6: MarkerFinder again, blocks in HOPACH order")
    ordered_covariates = list(dict.fromkeys(list(covariates) + [hopach_key]))
    step6 = generate_marker_heatmap_from_adata(
        adata,
        cluster_key=name_key,
        out=os.path.join(ordered_dir, "sctri_hopach_marker_heatmap.pdf"),
        top_n=top_n,
        marker_method="markerfinder",
        use_raw=use_raw,
        layer=layer,
        cells_per_cluster=cells_per_cluster,
        seed=seed,
        species=species,
        covariate_columns=ordered_covariates,
        validate_scaling=validate_scaling,
    )

    adt_outputs = None
    adt_mask = np.asarray([str(v).startswith(adt_prefix) for v in adata.var_names])
    if adt_mask.any() and (~adt_mask).any():
        os.makedirs(adt_dir, exist_ok=True)
        log(f"[annotate] antibody companion heatmap: {int(adt_mask.sum())} '{adt_prefix}' features, "
            f"same cells, same states, same order")
        adt = adata[:, adt_mask].copy()
        adt.var_names = pd.Index([str(v)[len(adt_prefix):] for v in adt.var_names])
        adt.var["gene_symbols"] = list(adt.var_names)
        adt.uns["lineage_order"] = lineage
        adt_outputs = generate_marker_heatmap_from_adata(
            adt,
            cluster_key=name_key,
            out=os.path.join(adt_dir, "sctri_hopach_ADT_marker_heatmap.pdf"),
            top_n=top_n,
            marker_method="markerfinder",
            use_raw=False,
            layer=layer,
            cells_per_cluster=cells_per_cluster,
            seed=seed,
            species=species,
            covariate_columns=ordered_covariates,
            validate_scaling=validate_scaling,
        )

    log(f"[annotate] wrote {table_path}")
    return {
        "cluster_table": table_path,
        "ordered_centroids": ordered_centroids_path,
        "marker_outputs_pruned": step4,
        "marker_outputs_hopach": step6,
        "marker_outputs_adt": adt_outputs,
        "names": names,
        "lineage_order": lineage,
        "hopach_groups": groups,
        "hopach_k": int(result.clust.k),
        "hopach_mss": float(result.clust.mss),
        "name_source": source,
    }
