"""Build the scALABLE side-artifacts a bundle needs, using scALABLE's own backends.

A precomputed bundle carries the matrix, the embedding, the markers and the DEG
tables. Three scALABLE views need artifacts a bundle does not carry:

  MarkerHeatmap     a genes x cells matrix npz
                    (format: marker_heatmap_h5ad.py:67 `_write_heatmap_cache`)
  MarkerNetwork     one NetPerspective interaction TSV per cell state
                    (generator: marker_heatmap_h5ad.py:757 `_export_marker_networks`)
  Cell communication  a fastComm scores TSV
                    (generator: altanalyze3.components.fastComm.cli `score`)

This module produces the marker networks by calling the same generator the cellHarmony
pipeline calls, and reads finished runs for the other two. It implements no scoring, no
layout and no statistic of its own.

Two further views are fed from artifacts the cellHarmony DIFFERENTIAL already wrote,
and must never be recomputed here:

  Differential GO terms      <contrast>/GeneSetEnrichment/GOElite_<tag>.tsv
                             (writer: flask/pipeline.py:2507)
  Differential networks      <contrast>/interaction-plots/<population>_interactions.tsv
                             (writer: flask/pipeline.py:2431-2456)

Pass `--goelite-from` and `--diff-networks-from` to ingest those files. The older
`--goelite` flag recomputes enrichment from the bundle DEG table and does NOT reproduce
the pipeline; it is kept only for a bundle with no differential run behind it.

Every viewer input is REQUIRED, and each one has a named opt-out
--------------------------------------------------------------
Before 2026-08-12 the marker heatmap and the cell communication view were optional
side effects of two optional flags: no `--marker-gct` meant no `heatmap_cache` key,
no `--fastcomm-dir` meant no `fastcomm_analysis` key, and the run still exited 0. The
COPD `assets_v3` and `assets_v4` builds omitted both flags, so the viewer served an
Explore tab with no marker heatmap and no cell communication and reported no error.

`check_completeness` now closes that hole. Each required input is either present or
explicitly waived by its own flag, and a run that satisfies neither exits non-zero:

  marker heatmap      auto-discovered beside --markers-tsv     waive: --no-marker-heatmap
  marker networks     built from the redundant marker table    waive: --skip-networks
  cell communication  auto-discovered, or --fastcomm-dir       waive: --no-fastcomm
  differential nets   --diff-networks-from                     waive: --skip-networks
  GO terms            --goelite-from                           waive: --no-goelite

Run:
  PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 \
  /opt/homebrew/opt/python@3.11/bin/python3.11 \
    -m altanalyze3.components.visualization.scalable_viewer.prepare_assets \
      --bundle-dir <dir> --prefix <prefix> --out <asset dir> \
      --markers-tsv <marker table> \
      [--marker-heatmap-npz <MarkerFinder *_fold_matrix.npz>] [--marker-gct <GCT>] \
      [--goelite-from <differential root>] [--diff-networks-from <differential root>] \
      [--fastcomm-dir <fastComm run dir>]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from . import bundle as B
from . import bundle_meta


#: Markers per cell state that `generate_marker_heatmap_from_adata` hands the network
#: exporter. marker_heatmap_h5ad.py:962-967 builds `redundant_selected` with a literal
#: 250, and marker_heatmap_h5ad.py:1176 passes that frame to `_export_marker_networks`.
#: The `network_top_n` argument only caps the frame, so it must not sit below 250.
REDUNDANT_MARKERS_PER_STATE = 250


def build_marker_networks(markers_tsv: str, out_dir: str, *, species: str = "human",
                          network_top_n: int = REDUNDANT_MARKERS_PER_STATE,
                          network_jobs: int = 1) -> List[Dict[str, str]]:
    """Per-cell-state marker networks, via marker_heatmap_h5ad's own exporter.

    The exporter is the function the cellHarmony pipeline reaches
    (flask/pipeline.py:1492 -> generate_marker_heatmap_from_adata(export_networks=True)
    -> marker_heatmap_h5ad.py:757). It is called directly here because the top-level
    entry point would re-run markerFinder on the source h5ad, and the bundle already
    holds that marker table.

    Feed this the REDUNDANT marker table, the one the pipeline feeds it. The pipeline
    writes two tables from one markerFinder run:

      <base>_markers.tsv                    the unique table, ~25 genes per state, the
                                            heatmap rows. NOT the network input.
      <base>_redundant_markers.tsv          250 genes per state, built at
                                            marker_heatmap_h5ad.py:962 and handed to the
                                            exporter at marker_heatmap_h5ad.py:1176.

    NetPerspective keeps a gene only when the gene interacts with another supplied gene.
    A 25-gene state rarely holds one such pair, so the unique table yields almost no
    networks. Reproducing the pipeline means passing the redundant table.
    """
    from altanalyze3.components.visualization import marker_heatmap_h5ad as MH

    stats = pd.read_csv(markers_tsv, sep="\t")
    required = {"Gene", "Fold", "cluster"}
    missing = required.difference(stats.columns)
    if missing:
        raise ValueError(f"{markers_tsv} is missing marker columns {sorted(missing)}; "
                         f"found {list(stats.columns)}")
    os.makedirs(out_dir, exist_ok=True)
    networks = MH._export_marker_networks(
        stats, out_dir, network_top_n=network_top_n, network_jobs=network_jobs, species=species
    )
    return [n for n in networks if n and os.path.isfile(str(n.get("tsv", "")))]


#: File-name suffixes `generate_marker_heatmap_from_adata` gives its heatmap artifacts
#: when the caller does not override them (marker_heatmap_h5ad.py:857-875), keyed by the
#: `marker_analysis` key each one fills (flask/pipeline.py:614-619).
MARKERFINDER_HEATMAP_SUFFIXES = {
    "heatmap_cache": "_fold_matrix.npz",       # the npz app.py:1294 loads
    "heatmap_tsv": "_fold_matrix.tsv",         # row-scaled fold matrix, optional
    "expression_tsv": "_exp_matrix.tsv",       # column-scaled expression, optional
}


def find_markerfinder_heatmap_cache(markers_tsv: str) -> Optional[str]:
    """Locate the MarkerFinder heatmap npz that sits beside the marker table.

    `generate_marker_heatmap_from_adata` writes the marker heatmap cache itself, at
    marker_heatmap_h5ad.py:1113-1128, through `_write_heatmap_cache`. That npz holds the
    four arrays `_get_marker_heatmap_cache_entry` reads (webapp/app.py:1294-1304):
    `matrix`, `row_ids`, `col_ids`, `col_barcodes`. It is the SAME artifact
    flask/pipeline.py:1527 registers as `marker_analysis["heatmap_cache"]`, so ingesting
    it reproduces the pipeline exactly.

    A Morpheus GCT is therefore never needed. `--marker-gct` remains for a project whose
    only surviving artifact is a GCT.

    The npz keeps its default name even when the caller renames the heatmap TSVs, because
    marker_heatmap_h5ad.py:868 only falls back on `--heatmap-cache`. The COPD run renamed
    both TSVs and still wrote `COPD_ms_marker_heatmap_fold_matrix.npz`. This globs on the
    suffix rather than on a reconstructed base name for that reason.
    """
    directory = Path(markers_tsv).parent
    hits = sorted(directory.glob(f"*{MARKERFINDER_HEATMAP_SUFFIXES['heatmap_cache']}"))
    if not hits:
        return None
    if len(hits) == 1:
        return str(hits[0])
    # More than one MarkerFinder run wrote into this directory. Prefer the npz whose base
    # name prefixes the marker table's, and never pick one silently at random.
    base = Path(markers_tsv).name
    for candidate in hits:
        stem = candidate.name[: -len(MARKERFINDER_HEATMAP_SUFFIXES["heatmap_cache"])]
        if base.startswith(stem.split("_marker")[0]):
            return str(candidate)
    raise ValueError(
        f"{directory} holds {len(hits)} MarkerFinder heatmap npz files and none matches "
        f"the marker table {markers_tsv}: {[str(h) for h in hits]}; pass --marker-heatmap-npz")


def ingest_marker_heatmap_cache(npz_path: str, out_npz: str, *, markers_tsv: Optional[str] = None,
                                clusters_tsv: Optional[str] = None) -> Dict[str, object]:
    """Copy the MarkerFinder heatmap npz into the asset directory and check it.

    Copies verbatim, then asserts the four arrays exist and agree in shape, so a
    truncated or foreign npz fails here instead of failing in the browser. Reports the
    marker-table overlap and the bundle-barcode overlap, both with their denominators.
    """
    import shutil

    import numpy as np

    source = Path(npz_path)
    if not source.is_file():
        raise FileNotFoundError(f"marker heatmap npz not found: {npz_path}")
    os.makedirs(os.path.dirname(out_npz), exist_ok=True)
    shutil.copy2(source, out_npz)

    with np.load(out_npz, allow_pickle=False) as npz:
        missing = [k for k in ("matrix", "row_ids", "col_ids", "col_barcodes") if k not in npz.files]
        if missing:
            raise AssertionError(
                f"{npz_path} is not a MarkerFinder heatmap cache: arrays {missing} are absent, "
                f"found {list(npz.files)}. Expected the npz `_write_heatmap_cache` writes "
                f"(marker_heatmap_h5ad.py:67).")
        matrix = np.asarray(npz["matrix"], dtype=np.float32)
        row_ids = np.asarray(npz["row_ids"], dtype=str)
        col_ids = np.asarray(npz["col_ids"], dtype=str)
        col_barcodes = np.asarray(npz["col_barcodes"], dtype=str)
    if matrix.shape != (row_ids.size, col_ids.size):
        raise AssertionError(
            f"{npz_path} matrix is {matrix.shape} but carries {row_ids.size} row ids and "
            f"{col_ids.size} column ids")
    if col_barcodes.size != col_ids.size:
        raise AssertionError(
            f"{npz_path} carries {col_barcodes.size} column barcodes for {col_ids.size} columns")

    # `row_ids` are "<cluster>:<gene>" (marker_heatmap_h5ad.py:69).
    genes = {value.split(":", 1)[1] if ":" in value else value for value in row_ids}
    states = {value.split(":", 1)[0] for value in row_ids if ":" in value}
    info: Dict[str, object] = {
        "path": out_npz, "source": str(source), "provenance": "MarkerFinder (ingested)",
        "n_rows": int(matrix.shape[0]), "n_cols": int(matrix.shape[1]),
        "n_genes": len(genes), "n_states": len(states),
    }

    if markers_tsv and os.path.isfile(markers_tsv):
        table_genes = set(pd.read_csv(markers_tsv, sep="\t", usecols=["Gene"])["Gene"].astype(str))
        shared = genes & table_genes
        info["marker_table_genes"] = len(table_genes)
        info["marker_table_overlap"] = len(shared)
        print(f"[assets]   heatmap genes shared with {markers_tsv}: {len(shared)} of "
              f"{len(genes)} heatmap genes, {len(table_genes)} marker-table genes", flush=True)
        if not shared:
            raise AssertionError(
                f"{npz_path} shares no gene with the marker table {markers_tsv}; the two come "
                f"from different runs")

    if clusters_tsv and os.path.isfile(clusters_tsv):
        bundle_barcodes = set(
            pd.read_csv(clusters_tsv, sep="\t", usecols=[0]).iloc[:, 0].astype(str))
        matched = len(set(col_barcodes) & bundle_barcodes)
        info["bundle_barcodes"] = len(bundle_barcodes)
        info["bundle_barcode_overlap"] = matched
        print(f"[assets]   heatmap columns matching bundle barcodes: {matched} of "
              f"{col_barcodes.size} ({len(bundle_barcodes)} barcodes in the bundle)", flush=True)
        if matched == 0:
            # Not fatal: `/marker/heatmap.tsv` with no display filter never touches the
            # barcodes (app.py:4848-4870). Only the filtered path intersects them
            # (app.py:3658), and it returns 404 when the intersection is empty.
            print(f"[WARN] no heatmap column barcode matches a bundle barcode. The marker "
                  f"heatmap renders, but the two display filters return 404 on it. This is "
                  f"expected when MarkerFinder ran on pseudobulk columns rather than on the "
                  f"bundle's own cells.", flush=True)
    return info


def find_redundant_markers_tsv(markers_tsv: str) -> Optional[str]:
    """Locate the 250-per-state marker table beside the unique one.

    `generate_marker_heatmap_from_adata` writes both files into one directory and names
    the second `<base>_redundant_markers.tsv` (marker_heatmap_h5ad.py:875).
    """
    directory = Path(markers_tsv).parent
    hits = sorted(directory.glob("*_redundant_markers.tsv"))
    return str(hits[0]) if hits else None


def write_marker_fold_lookup(unique_tsv: str, redundant_tsv: Optional[str], out_path: str) -> Dict[str, object]:
    """One fold-change table covering every gene either marker table names.

    `_build_marker_network_payload` (app.py:2217) colours each network node by looking the
    gene up in `marker_analysis["markers_tsv"]`, filtered to that cell state. Networks now
    come from the redundant table, and the two tables do not contain each other: the COPD
    run holds 121 of 964 unique-table (state, gene) pairs that the redundant table omits.
    Either table alone therefore leaves nodes uncoloured at fold change 0.

    This writes the union, keyed on (cluster, Gene), and keeps the redundant row when both
    tables carry the same pair.
    """
    unique = pd.read_csv(unique_tsv, sep="\t")
    frames = [unique]
    if redundant_tsv and os.path.isfile(redundant_tsv):
        frames.insert(0, pd.read_csv(redundant_tsv, sep="\t"))
    merged = pd.concat(frames, ignore_index=True)
    for column in ("cluster", "Gene"):
        if column not in merged.columns:
            raise ValueError(f"marker table is missing the {column!r} column; found {list(merged.columns)}")
        merged[column] = merged[column].astype(str)
    merged = merged.drop_duplicates(subset=["cluster", "Gene"], keep="first")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    merged.to_csv(out_path, sep="\t", index=False)

    # Conservation check: the union must name every pair either input names.
    written = set(zip(merged["cluster"], merged["Gene"]))
    for label, path in (("unique", unique_tsv), ("redundant", redundant_tsv)):
        if not path or not os.path.isfile(path):
            continue
        source = pd.read_csv(path, sep="\t")
        pairs = set(zip(source["cluster"].astype(str), source["Gene"].astype(str)))
        missing = pairs.difference(written)
        if missing:
            raise AssertionError(
                f"fold lookup lost {len(missing)} of {len(pairs)} (cluster, Gene) pairs "
                f"from the {label} table {path}")
    return {"path": out_path, "n_rows": int(len(merged)),
            "n_pairs": len(written), "n_genes": int(merged["Gene"].nunique()),
            "unique_tsv": unique_tsv, "redundant_tsv": redundant_tsv or ""}


def build_differential_networks(deg_tsv: str, out_dir: str, *, species: str = "Hs") -> List[Dict[str, str]]:
    """Per-population differential interaction networks.

    This is the same NetPerspective call the cellHarmony differential makes at
    flask/pipeline.py:2431-2456: same gene table columns, same `pval_column="fdr"`,
    same `max_genes=None`. Only the DEG frame is supplied from the bundle instead of
    from a live differential run.
    """
    from altanalyze3.components.visualization import NetPerspective

    detailed = pd.read_csv(deg_tsv, sep="\t")
    if "population" not in detailed.columns:
        raise ValueError(f"{deg_tsv} has no 'population' column")
    interactions = NetPerspective.load_interaction_data(species=species)
    os.makedirs(out_dir, exist_ok=True)
    networks: List[Dict[str, str]] = []
    used: set = set()
    for index, (population, pop_df) in enumerate(detailed.groupby("population"), start=1):
        stats = pop_df.copy().reset_index(drop=True)
        stats["gene"] = stats["gene"].astype(str)
        stats["log2fc"] = pd.to_numeric(stats.get("log2fc"), errors="coerce")
        keep = [c for c in ("gene", "log2fc", "fdr", "pval") if c in stats.columns]
        selected = stats.loc[:, keep].dropna(subset=["gene", "log2fc"]).drop_duplicates(subset=["gene"])
        if selected.empty or selected["gene"].nunique() < 2:
            continue
        network_id = NetPerspective.safe_component(str(population), fallback=f"network_{index}")
        if network_id in used:
            network_id = f"{network_id}_{index}"
        used.add(network_id)
        prefix = os.path.join(out_dir, network_id)
        try:
            pdf_path, png_path, tsv_path = NetPerspective.generate_network_for_genes(
                selected, interactions, prefix,
                gene_column="gene", fold_change_column="log2fc",
                pval_column="fdr" if "fdr" in selected.columns else ("pval" if "pval" in selected.columns else None),
                max_genes=None,
            )
        except NetPerspective.NetworkGenerationError as exc:
            print(f"[WARN] no differential interaction edges for {population}: {exc}", flush=True)
            continue
        networks.append({"id": network_id, "population": str(population),
                         "pdf": pdf_path, "png": png_path, "tsv": tsv_path})
    return networks


def build_goelite(deg_tsv: str, background_genes_tsv: str, out_dir: str, *,
                  comparison_tag: str, species: str = "Hs") -> Optional[str]:
    """GO-Elite enrichment per cell state, through cellHarmony's own runner.

    Calls `cellHarmony_differential.run_goelite_for_clusters` (the function
    flask/pipeline.py:2297 calls) with a `de_store` built from the bundle's DEG table.
    The runner writes `GOElite_<tag>.tsv`, which is the artifact
    `_get_differential_go_table` (app.py:1072) reads.
    """
    import numpy as np
    from altanalyze3.components.cellHarmony import cellHarmony_differential as CD
    from altanalyze3.components.visualization import NetPerspective

    detailed = pd.read_csv(deg_tsv, sep="\t")
    assigned = pd.DataFrame({
        "gene": detailed["gene"].astype(str),
        "group": "local",
        "population_or_pattern": detailed["population"].astype(str),
        "key_log2fc": pd.to_numeric(detailed["log2fc"], errors="coerce"),
    }).dropna(subset=["key_log2fc"])
    de_store = {"assigned_groups": assigned}

    background = pd.read_csv(background_genes_tsv, sep="\t")
    col = "symbol" if "symbol" in background.columns else background.columns[-1]
    background_genes = background[col].astype(str)

    os.makedirs(out_dir, exist_ok=True)
    payload = CD.run_goelite_for_clusters(
        de_store=de_store,
        background_genes=background_genes,
        outdir=out_dir,
        comparison_tag=comparison_tag,
        species=species,
    )
    if not payload:
        return None
    path = os.path.join(out_dir, f"GOElite_{NetPerspective.safe_component(comparison_tag)}.tsv")
    return path if os.path.isfile(path) else None


# =================================================================================
# Pipeline artifact ingestion
#
# GO-Elite tables and differential interaction networks are STANDARD ARTIFACTS of the
# cellHarmony differential run. The run writes them at flask/pipeline.py:2431-2456
# (networks) and flask/pipeline.py:2507 (GO-Elite), into
#
#   <contrast dir>/GeneSetEnrichment/GOElite_<safe(tag)>.tsv
#   <contrast dir>/interaction-plots/<safe(population)>_interactions.tsv
#
# The viewer must ship THOSE files. `build_goelite` and `build_differential_networks`
# above recompute an approximation from the bundle's DEG_detailed table, and the two do
# not agree: the pipeline enriches the ASSIGNED-GROUP gene lists
# (DEGs/DEG_assigned_groups_<tag>.tsv, which carry the `local`, `co-regulated` and
# `overall` groups), while the recompute enriches every DEG_detailed row with a finite
# log2fc. On the COPD run that inflated the AT1__down query from 57 genes to 211 and
# dropped all four `coreg_*` populations. Prefer these ingest functions.
# =================================================================================

#: Subdirectory names the cellHarmony differential writes under each contrast directory.
PIPELINE_GOELITE_SUBDIR = "GeneSetEnrichment"
PIPELINE_NETWORK_SUBDIR = "interaction-plots"


def _count_data_rows(path: str) -> int:
    """Non-empty lines minus the header, matching `awk 'END{print NR}'` minus one.

    Counts by iteration so a file without a trailing newline is not undercounted.
    """
    total = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.strip():
                total += 1
    return max(total - 1, 0)


def find_pipeline_contrast_dir(root: str, comparison: str) -> Optional[str]:
    """Locate the cellHarmony contrast directory for `comparison` under `root`.

    The differential writes `<root>/<run name>/<comparison tag>/`, where the leaf
    directory is named for the comparison tag the bundle also carries (for example
    `.../differential_ms/copd_vs_noncopd/COPD_vs_non-COPD/`). Accepts `root` pointing at
    the run collection, at one run, or at the contrast directory itself.
    """
    base = Path(root)
    if not base.is_dir():
        raise FileNotFoundError(f"pipeline artifact root not found: {root}")
    if base.name == comparison:
        return str(base)
    direct = base / comparison
    if direct.is_dir():
        return str(direct)
    hits = sorted(p for p in base.glob(f"*/{comparison}") if p.is_dir())
    if len(hits) > 1:
        raise ValueError(
            f"{root} holds {len(hits)} directories named {comparison!r}: "
            f"{[str(h) for h in hits]}; pass the contrast directory itself")
    return str(hits[0]) if hits else None


def ingest_goelite(contrast_dir: str, out_dir: str, *, comparison_tag: str) -> Dict[str, object]:
    """Copy the differential run's own GO-Elite table into the asset directory.

    Reads `<contrast_dir>/GeneSetEnrichment/GOElite_<safe(tag)>.tsv`, the exact file
    flask/pipeline.py:2507 registers as the `goelite_tsv` artifact and the file
    `_get_differential_go_table` (webapp/app.py:1072) reads. Copies it verbatim; this
    function computes no enrichment.
    """
    import shutil

    from altanalyze3.components.visualization import NetPerspective

    name = f"GOElite_{NetPerspective.safe_component(comparison_tag)}.tsv"
    source = Path(contrast_dir) / PIPELINE_GOELITE_SUBDIR / name
    if not source.is_file():
        raise FileNotFoundError(
            f"the differential run under {contrast_dir} has no {PIPELINE_GOELITE_SUBDIR}/{name}; "
            f"re-run the cellHarmony differential with GO terms enabled, or drop --goelite-from")
    os.makedirs(out_dir, exist_ok=True)
    target = Path(out_dir) / name
    shutil.copy2(source, target)

    source_rows = _count_data_rows(str(source))
    target_rows = _count_data_rows(str(target))
    if source_rows != target_rows:
        raise AssertionError(
            f"GO-Elite copy lost rows: {source} has {source_rows}, {target} has {target_rows}")
    frame = pd.read_csv(target, sep="\t")
    populations = sorted(str(v) for v in frame["population"].dropna().unique()) \
        if "population" in frame.columns else []
    n_selected = int((frame["selected"].astype(str) == "True").sum()) \
        if "selected" in frame.columns else -1
    return {"path": str(target), "source": str(source), "n_terms": source_rows,
            "n_selected": n_selected, "n_populations": len(populations),
            "populations": populations}


def ingest_differential_networks(contrast_dir: str, deg_tsv: str, out_dir: str) -> List[Dict[str, str]]:
    """Copy the differential run's own interaction networks into the asset directory.

    Rebuilds the population -> filename map the pipeline used, rather than guessing a
    population back out of a filename: flask/pipeline.py:2431-2443 groups DEG_detailed by
    `population`, names each network `NetPerspective.safe_component(population)` and
    de-duplicates a repeated id as `<id>_<index>`. This walks the same groups in the same
    order and copies the `<id>.pdf`, `<id>.png` and `<id>_interactions.tsv` the run wrote.
    A population with no file is reported, never invented.
    """
    import shutil

    from altanalyze3.components.visualization import NetPerspective

    net_dir = Path(contrast_dir) / PIPELINE_NETWORK_SUBDIR
    if not net_dir.is_dir():
        raise FileNotFoundError(f"no {PIPELINE_NETWORK_SUBDIR}/ under {contrast_dir}")
    detailed = pd.read_csv(deg_tsv, sep="\t")
    if "population" not in detailed.columns:
        raise ValueError(f"{deg_tsv} has no 'population' column")
    os.makedirs(out_dir, exist_ok=True)

    networks: List[Dict[str, str]] = []
    used: set = set()
    absent: List[str] = []
    for index, (population, _pop_df) in enumerate(detailed.groupby("population"), start=1):
        network_id = NetPerspective.safe_component(str(population), fallback=f"network_{index}")
        if network_id in used:
            network_id = f"{network_id}_{index}"
        used.add(network_id)
        source_tsv = net_dir / f"{network_id}_interactions.tsv"
        if not source_tsv.is_file():
            absent.append(str(population))
            continue
        entry: Dict[str, str] = {"id": network_id, "population": str(population)}
        target_tsv = Path(out_dir) / source_tsv.name
        shutil.copy2(source_tsv, target_tsv)
        source_rows = _count_data_rows(str(source_tsv))
        target_rows = _count_data_rows(str(target_tsv))
        if source_rows != target_rows:
            raise AssertionError(
                f"network copy lost rows: {source_tsv} has {source_rows}, "
                f"{target_tsv} has {target_rows}")
        entry["tsv"] = str(target_tsv)
        entry["source_tsv"] = str(source_tsv)
        entry["n_edges"] = source_rows
        for kind in ("pdf", "png"):
            source_plot = net_dir / f"{network_id}.{kind}"
            if source_plot.is_file():
                target_plot = Path(out_dir) / source_plot.name
                shutil.copy2(source_plot, target_plot)
                entry[kind] = str(target_plot)
        networks.append(entry)
    if absent:
        # The pipeline skips a population with under two genes or no interaction edge
        # (flask/pipeline.py:2438, 2459). Reported, never silently dropped.
        print(f"[assets]   {len(absent)} of {len(used)} populations have no network file "
              f"in {net_dir}: {absent}", flush=True)
    return networks


# =================================================================================
# Cell communication: find the finished fastComm run
# =================================================================================

#: Directory name the fastComm runs of a scalable_viewer deployment live under. A run
#: for dataset `<id>` is `<root>/fastComm/<id>/`, holding `fastcomm_scores.tsv`.
FASTCOMM_SUBDIR = "fastComm"


def find_fastcomm_run(dataset_id: str, *, bundle_dir: str, out_dir: str) -> Optional[str]:
    """Locate a finished fastComm run for this dataset, or return None.

    A fastComm run is not written beside the bundle, so this searches the conventional
    layout of a scalable_viewer deployment: the bundle root and the asset root are
    siblings of `fastComm/`, and each holds one directory per dataset.

      <scalable_viewer>/bundles_v3/<id>/        the bundle
      <scalable_viewer>/assets_v5/<id>/         the assets
      <scalable_viewer>/fastComm/<id>/          the fastComm run

    Every candidate it tried is printed, so a miss is diagnosable rather than silent.
    """
    bundle = Path(bundle_dir).resolve()
    out = Path(out_dir).resolve()
    candidates = [
        bundle.parent.parent / FASTCOMM_SUBDIR / dataset_id,
        out.parent.parent / FASTCOMM_SUBDIR / dataset_id,
        bundle.parent.parent / FASTCOMM_SUBDIR,
        out.parent.parent / FASTCOMM_SUBDIR,
    ]
    seen: List[str] = []
    for candidate in candidates:
        text = str(candidate)
        if text in seen:
            continue
        seen.append(text)
        if (candidate / "fastcomm_scores.tsv").is_file():
            print(f"[assets] fastComm run discovered at {candidate}", flush=True)
            return text
    print(f"[assets] no fastComm run found; tried {seen}", flush=True)
    return None


def verify_fastcomm_matches_bundle(run_dir: str, bundle_meta_json: Dict, states: List[str],
                                   cluster_key: str,
                                   populations: Optional[List[str]] = None) -> Dict[str, object]:
    """Fail when the fastComm run was not scored on the cells this bundle serves.

    A fastComm run is a separate job, so nothing links it to a bundle except the paths an
    operator types. The COPD deployment hit exactly that: the run under
    `scalable_viewer/fastComm/COPD-metacells/` scored 161,432 cells of
    `COPD_metacells.deid.h5ad`, while `bundles_v3` serves 123,076 cells of
    `COPD_metacells.metasample.h5ad`. Both carry 39 cell states, so no count in the viewer
    would have looked wrong, and the communication scores would have described other cells.

    `fastComm.cli score` writes `summary.json` with `n_cells`, `n_states` and `state_key`
    (the fields this reads). A run with no `summary.json` cannot be checked and is
    rejected, because an unverifiable provenance is not a provenance.
    """
    summary_path = Path(run_dir) / "summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(
            f"{run_dir} carries no summary.json, so its provenance cannot be checked against "
            f"the bundle. Re-run `fastComm.cli score` with --summary-json, or pass "
            f"--allow-stale-fastcomm to ship it unchecked.")
    with open(summary_path, "r") as fh:
        summary = json.load(fh)

    sv = bundle_meta_json.get("scalable_viewer", {})
    bundle_cells = int(sv.get("n_cells") or bundle_meta_json.get("n_cells") or 0)
    run_cells = int(summary.get("n_cells") or 0)
    run_states = int(summary.get("n_states") or 0)
    run_state_key = str(summary.get("state_key") or "")

    problems: List[str] = []
    if run_state_key and run_state_key != cluster_key:
        problems.append(f"state key: bundle {cluster_key!r}, fastComm {run_state_key!r}")
    if bundle_cells and run_cells and bundle_cells != run_cells:
        problems.append(f"cells: bundle {bundle_cells}, fastComm {run_cells}")

    # A fastComm state set SMALLER than the bundle's is normal: `--min-cells` drops a
    # state with too few cells to score. A state fastComm names and the bundle does not
    # is not normal - it means the two ran on different data.
    bundle_states = {str(value) for value in (states or [])}
    scored_states = {str(value) for value in (populations or [])}
    unknown = sorted(scored_states - bundle_states) if bundle_states and scored_states else []
    dropped = sorted(bundle_states - scored_states) if bundle_states and scored_states else []
    if unknown:
        problems.append(
            f"{len(unknown)} cell state(s) scored by fastComm are absent from the bundle: {unknown}")
    if problems:
        raise AssertionError(
            f"the fastComm run at {run_dir} does not match this bundle - "
            + "; ".join(problems)
            + f". Re-score the bundle's own h5ad ({bundle_meta_json.get('source_h5ad')}) with "
              f"`fastComm.cli score`, or pass --allow-stale-fastcomm to ship it on purpose.")
    if dropped:
        print(f"[assets]   {len(dropped)} of {len(bundle_states)} bundle cell states carry no "
              f"fastComm score: {dropped}. `fastComm.cli score --min-cells` drops a state with "
              f"too few cells.", flush=True)
    return {"summary_json": str(summary_path), "n_cells": run_cells, "n_states": run_states,
            "state_key": run_state_key, "bundle_n_cells": bundle_cells,
            "bundle_n_states": len(bundle_states), "n_states_scored": len(scored_states),
            "states_without_scores": dropped, "verified": True}


# =================================================================================
# Completeness: every viewer input is present, or explicitly waived
# =================================================================================

def check_completeness(assets: Dict[str, object], *, waived: Dict[str, str]) -> List[str]:
    """Report every viewer input, and name the ones that are absent and not waived.

    Returns the list of failures. `main` raises on a non-empty list. The report prints
    whether or not it passes, so a build log always states what the viewer will serve.

    `waived` maps an input key to the flag the operator passed to skip it. An input in
    `waived` is reported as WAIVED and never fails the run.
    """
    differential = assets.get("differential") or {}
    n_contrasts = len(differential)
    n_diff_nets = sum(len(entry.get("networks") or []) for entry in differential.values())
    n_goelite = sum(1 for entry in differential.values() if entry.get("goelite_tsv"))
    fastcomm = assets.get("fastcomm_analysis") or {}
    heatmap_info = assets.get("heatmap_cache_info") or {}

    checks = [
        ("marker_heatmap", "Explore / MarkerHeatmap", bool(assets.get("heatmap_cache")),
         f"{heatmap_info.get('n_rows', 0)} rows x {heatmap_info.get('n_cols', 0)} columns"),
        ("marker_networks", "Explore / MarkerNetwork", bool(assets.get("networks")),
         f"{len(assets.get('networks') or [])} cell states"),
        ("fastcomm", "Explore / Cell communication", bool(fastcomm.get("enabled")),
         f"{fastcomm.get('n_rows', 0)} scored edges over "
         f"{len(fastcomm.get('populations') or [])} cell states"),
        ("diff_networks", "Differential / Network", n_diff_nets > 0,
         f"{n_diff_nets} networks over {n_contrasts} contrasts"),
        ("goelite", "Differential / GO Terms", n_goelite > 0,
         f"{n_goelite} of {n_contrasts} contrasts"),
    ]
    failures: List[str] = []
    print("[assets] ---- completeness ----", flush=True)
    for key, label, present, detail in checks:
        if key in waived:
            print(f"[assets]   WAIVED  {label:32s} ({waived[key]})", flush=True)
            continue
        if present:
            print(f"[assets]   OK      {label:32s} {detail}", flush=True)
            continue
        print(f"[assets]   MISSING {label:32s} -", flush=True)
        failures.append(label)
    return failures


def top_marker_per_state(markers_tsv: str, order_tsv: Optional[str] = None) -> Dict[str, object]:
    """The default dot-plot gene set: the single best marker of every cell state.

    "Best" is the marker table's own top row per cluster, ranked by FDR then Fold -
    the same ranking `_export_marker_networks` applies (marker_heatmap_h5ad.py:775).
    Cell states follow the canonical order file when one is given.
    """
    stats = pd.read_csv(markers_tsv, sep="\t")
    stats["cluster"] = stats["cluster"].astype(str)
    stats["Gene"] = stats["Gene"].astype(str)
    stats["Fold"] = pd.to_numeric(stats.get("Fold"), errors="coerce")
    fdr_col = "FDR p-value" if "FDR p-value" in stats.columns else None
    if fdr_col:
        stats[fdr_col] = pd.to_numeric(stats[fdr_col], errors="coerce")
        ranked = stats.sort_values(["cluster", fdr_col, "Fold"], ascending=[True, True, False])
    else:
        ranked = stats.sort_values(["cluster", "Fold"], ascending=[True, False])
    best = ranked.dropna(subset=["Gene", "Fold"]).drop_duplicates(subset=["cluster"], keep="first")

    order: List[str] = []
    if order_tsv and os.path.isfile(order_tsv):
        frame = pd.read_csv(order_tsv, sep="\t")
        col = "cell_state" if "cell_state" in frame.columns else frame.columns[1]
        order = [str(v) for v in frame[col].tolist()]
    by_state = dict(zip(best["cluster"], best["Gene"]))
    ordered_states = [s for s in order if s in by_state] + \
                     [s for s in by_state if s not in set(order)]
    genes: List[str] = []
    pairs: List[Dict[str, str]] = []
    for state in ordered_states:
        gene = str(by_state[state])
        pairs.append({"state": state, "gene": gene})
        if gene not in genes:
            genes.append(gene)
    return {"pairs": pairs, "genes": genes, "n_states": len(pairs),
            "n_unique_genes": len(genes), "markers_tsv": markers_tsv,
            "order_tsv": order_tsv or ""}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundle-dir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--out", required=True, help="asset directory (never the bundle directory)")
    ap.add_argument("--marker-heatmap-npz", default=None,
                    help="MarkerFinder heatmap cache (*_fold_matrix.npz). Defaults to the one "
                         "beside --markers-tsv, which generate_marker_heatmap_from_adata wrote")
    ap.add_argument("--no-marker-heatmap", action="store_true",
                    help="build without the Explore marker heatmap, on purpose")
    ap.add_argument("--marker-gct", default=None,
                    help="per-cell marker GCT to convert to the npz. Only for a project whose "
                         "MarkerFinder *_fold_matrix.npz is gone; prefer --marker-heatmap-npz")
    ap.add_argument("--markers-tsv", default=None, help="marker table with Gene/Fold/cluster")
    ap.add_argument("--network-markers-tsv", default=None,
                    help="redundant marker table (250 genes per state) used for the marker "
                         "networks; defaults to *_redundant_markers.tsv beside --markers-tsv")
    ap.add_argument("--order-tsv", default=None, help="canonical cell-state order TSV")
    ap.add_argument("--study-id", default=None,
                    help="LungMAP study id this bundle belongs to, for example "
                         "lmdata:LMEX0000009416. The Study tab reads it. Without it the "
                         "Study tab reports that no id is configured, and never shows "
                         "another study's record")
    ap.add_argument("--fastcomm-dir", default=None,
                    help="finished fastComm run directory. Defaults to the conventional "
                         "<scalable_viewer>/fastComm/<dataset id>/ beside the bundle root")
    ap.add_argument("--no-fastcomm", action="store_true",
                    help="build without the Explore cell communication view, on purpose")
    ap.add_argument("--allow-stale-fastcomm", action="store_true",
                    help="ship a fastComm run whose cell count, cell states or state key "
                         "disagree with the bundle. Off by default: a run scored on other "
                         "cells describes other cells")
    ap.add_argument("--fastcomm-sample-key", default="",
                    help="obs column used as the fastComm per-sample split key")
    ap.add_argument("--species", default="human")
    ap.add_argument("--skip-networks", action="store_true")
    ap.add_argument("--goelite", action="store_true",
                    help="RECOMPUTE GO-Elite per contrast from the bundle DEG table "
                         "(downloads go-basic.obo and goa.gaf.gz). This does NOT reproduce the "
                         "cellHarmony differential's own enrichment - it enriches every "
                         "DEG_detailed row instead of the assigned-group gene lists, and it "
                         "omits the coreg_* populations. Prefer --goelite-from.")
    ap.add_argument("--goelite-from", default=None,
                    help="ingest the cellHarmony differential's OWN GeneSetEnrichment/ tables. "
                         "Point at the directory holding the runs (for example "
                         ".../analysis/differential_ms); each contrast is found as "
                         "<root>/*/<comparison tag>/GeneSetEnrichment/GOElite_<tag>.tsv")
    ap.add_argument("--diff-networks-from", default=None,
                    help="ingest the cellHarmony differential's OWN interaction-plots/ networks "
                         "instead of recomputing them; same root layout as --goelite-from")
    ap.add_argument("--no-goelite", action="store_true",
                    help="build without the Differential GO Terms view, on purpose")
    a = ap.parse_args(argv)

    if a.goelite and a.goelite_from:
        ap.error("--goelite recomputes enrichment and --goelite-from ingests the pipeline's own; "
                 "pass only one")
    if a.no_marker_heatmap and (a.marker_heatmap_npz or a.marker_gct):
        ap.error("--no-marker-heatmap waives the marker heatmap; drop it or drop the heatmap input")
    if a.no_fastcomm and a.fastcomm_dir:
        ap.error("--no-fastcomm waives cell communication; drop it or drop --fastcomm-dir")

    paths = B.BundlePaths(a.bundle_dir, a.prefix)
    meta = B.read_metadata(paths.metadata)
    cluster_key = meta.get("cluster_key", "cell_state")
    markers_tsv = a.markers_tsv or meta.get("scalable_viewer", {}).get("markers", {}).get("source")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    assets: Dict[str, object] = {"bundle_dir": paths.bundle_dir, "prefix": paths.prefix}

    # The Study tab reads this. `study_ids_for_dataset` (scalable_app.py) prefers the
    # manifest over the bundle metadata, because rebuilding a manifest costs minutes and
    # rebuilding a bundle costs hours.
    study_id = a.study_id or str(meta.get("scalable_viewer", {}).get("study_id") or "").strip()
    if study_id:
        assets["study_id"] = study_id
        print(f"[assets] LungMAP study id: {study_id}", flush=True)
    else:
        print(f"[WARN] no --study-id and none in the bundle metadata. The Study tab will "
              f"report that no id is configured, and will show no study record.", flush=True)

    # ---- Explore / MarkerHeatmap ------------------------------------------------
    # Preference order, most faithful first:
    #   1. --marker-heatmap-npz            the operator names the MarkerFinder npz
    #   2. auto-discovery beside the marker table   the npz MarkerFinder already wrote
    #   3. --marker-gct                    reshape a Morpheus GCT, for a project with no npz
    # Option 2 is what the cellHarmony pipeline registers (flask/pipeline.py:1527), so a
    # standard MarkerFinder run needs no flag at all.
    heatmap_npz = a.marker_heatmap_npz
    if not heatmap_npz and not a.marker_gct and not a.no_marker_heatmap and markers_tsv:
        heatmap_npz = find_markerfinder_heatmap_cache(markers_tsv)
        if heatmap_npz:
            print(f"[assets] MarkerFinder heatmap npz discovered beside {markers_tsv}: "
                  f"{heatmap_npz}", flush=True)
    if heatmap_npz:
        npz = str(out / f"{a.prefix}_marker_heatmap_cache.npz")
        info = ingest_marker_heatmap_cache(
            heatmap_npz, npz, markers_tsv=markers_tsv, clusters_tsv=paths.clusters)
        assets["heatmap_cache"] = npz
        assets["heatmap_cache_info"] = info
        print(f"[assets] marker heatmap npz {npz}: {info['n_rows']} rows x {info['n_cols']} "
              f"columns, {info['n_genes']} genes over {info['n_states']} cell states, "
              f"ingested from {info['source']}", flush=True)
    elif a.marker_gct:
        npz = str(out / f"{a.prefix}_marker_heatmap_cache.npz")
        info = bundle_meta.gct_to_marker_cache(a.marker_gct, npz)
        info["provenance"] = "Morpheus GCT (reshaped)"
        info["source"] = a.marker_gct
        assets["heatmap_cache"] = npz
        assets["heatmap_cache_info"] = info
        print(f"[assets] marker heatmap npz {npz}: {info['n_rows']} genes x {info['n_cols']} cells, "
              f"reshaped from the GCT {a.marker_gct}", flush=True)

    if markers_tsv and os.path.isfile(markers_tsv):
        assets["markers_tsv"] = markers_tsv
        assets["dotplot_default"] = top_marker_per_state(markers_tsv, a.order_tsv)
        d = assets["dotplot_default"]
        print(f"[assets] dot-plot default: {d['n_states']} states -> {d['n_unique_genes']} unique genes",
              flush=True)
        # Networks come from the redundant table, the frame the pipeline passes the
        # exporter (marker_heatmap_h5ad.py:1176). Node colours come from the union of both
        # tables, because neither one contains the other.
        network_markers_tsv = a.network_markers_tsv or find_redundant_markers_tsv(markers_tsv)
        if network_markers_tsv and not os.path.isfile(network_markers_tsv):
            raise FileNotFoundError(f"--network-markers-tsv not found: {network_markers_tsv}")
        assets["network_markers_tsv"] = network_markers_tsv or ""
        lookup = write_marker_fold_lookup(
            markers_tsv, network_markers_tsv, str(out / f"{a.prefix}_marker_fold_lookup.tsv"))
        assets["marker_fold_lookup_tsv"] = lookup["path"]
        assets["marker_fold_lookup_info"] = lookup
        print(f"[assets] fold lookup {lookup['path']}: {lookup['n_rows']} rows, "
              f"{lookup['n_genes']} genes", flush=True)
        if not a.skip_networks:
            net_dir = str(out / f"{a.prefix}_marker_networks")
            if not network_markers_tsv:
                print(f"[WARN] no *_redundant_markers.tsv beside {markers_tsv}; marker networks "
                      f"fall back to the unique table and will be sparse", flush=True)
            source_tsv = network_markers_tsv or markers_tsv
            per_state = pd.read_csv(source_tsv, sep="\t").groupby("cluster").size()
            print(f"[assets] marker network input {source_tsv}: {len(per_state)} states, "
                  f"{int(per_state.min())}-{int(per_state.max())} genes per state", flush=True)
            networks = build_marker_networks(source_tsv, net_dir, species=a.species)
            assets["networks"] = networks
            print(f"[assets] marker networks: {len(networks)} of {len(per_state)} states "
                  f"written under {net_dir}", flush=True)

    # Per-contrast differential artifacts: interaction networks and GO-Elite terms.
    per_state = [c for c in meta.get("scalable_viewer", {}).get("deg", {}).get("comparisons", [])
                 if c.get("kind") == "per_cell_state"]
    differential: Dict[str, Dict[str, object]] = {}
    for comparison in per_state:
        entry: Dict[str, object] = {}
        tag = comparison["comparison"]

        # Where the cellHarmony differential put this contrast's own artifacts.
        goelite_contrast_dir = (find_pipeline_contrast_dir(a.goelite_from, tag)
                                if a.goelite_from else None)
        network_contrast_dir = (find_pipeline_contrast_dir(a.diff_networks_from, tag)
                                if a.diff_networks_from else None)

        if not a.skip_networks:
            net_dir = str(out / f"{a.prefix}_diff_networks" / tag)
            if a.diff_networks_from:
                if not network_contrast_dir:
                    raise FileNotFoundError(
                        f"--diff-networks-from {a.diff_networks_from} holds no contrast directory "
                        f"named {tag!r}")
                nets = ingest_differential_networks(network_contrast_dir, comparison["path"], net_dir)
                edges = sum(int(n.get("n_edges", 0)) for n in nets)
                entry["networks"] = nets
                entry["networks_source"] = str(Path(network_contrast_dir) / PIPELINE_NETWORK_SUBDIR)
                entry["networks_provenance"] = "cellHarmony differential (ingested)"
                print(f"[assets] differential networks {tag}: {len(nets)} networks, {edges} edges "
                      f"ingested from {entry['networks_source']}", flush=True)
            else:
                nets = build_differential_networks(comparison["path"], net_dir, species="Hs")
                entry["networks"] = nets
                entry["networks_provenance"] = "recomputed by prepare_assets"
                print(f"[assets] differential networks {tag}: {len(nets)} (RECOMPUTED)", flush=True)

        if a.goelite_from:
            if not goelite_contrast_dir:
                raise FileNotFoundError(
                    f"--goelite-from {a.goelite_from} holds no contrast directory named {tag!r}")
            go_dir = str(out / f"{a.prefix}_goelite" / tag)
            info = ingest_goelite(goelite_contrast_dir, go_dir, comparison_tag=tag)
            entry["goelite_tsv"] = info["path"]
            entry["goelite_info"] = info
            entry["goelite_provenance"] = "cellHarmony differential (ingested)"
            print(f"[assets] GO-Elite {tag}: {info['n_terms']} terms "
                  f"({info['n_selected']} selected) across {info['n_populations']} populations, "
                  f"ingested from {info['source']}", flush=True)
        if a.goelite:
            go_dir = str(out / f"{a.prefix}_goelite" / comparison["comparison"])
            try:
                go_tsv = build_goelite(comparison["path"], paths.genes, go_dir,
                                       comparison_tag=comparison["comparison"], species="Hs")
            except Exception as exc:                       # recorded, never silent
                go_tsv = None
                print(f"[WARN] GO-Elite failed for {comparison['comparison']}: "
                      f"{type(exc).__name__}: {exc}", flush=True)
            if go_tsv:
                entry["goelite_tsv"] = go_tsv
                entry["goelite_provenance"] = "recomputed by prepare_assets"
            print(f"[assets] GO-Elite {comparison['comparison']}: {go_tsv} (RECOMPUTED)", flush=True)
        if entry:
            differential[comparison["id"]] = entry
    if differential:
        assets["differential"] = differential

    # ---- Explore / Cell communication -------------------------------------------
    dataset_id = str(meta.get("scalable_viewer", {}).get("id") or Path(a.bundle_dir).name)
    fastcomm_dir = a.fastcomm_dir
    if not fastcomm_dir and not a.no_fastcomm:
        fastcomm_dir = find_fastcomm_run(dataset_id, bundle_dir=a.bundle_dir, out_dir=str(out))
    if fastcomm_dir:
        assets["fastcomm_dir"] = fastcomm_dir
        analysis = bundle_meta.fastcomm_analysis_from_dir(fastcomm_dir, cluster_key,
                                                          sample_key=a.fastcomm_sample_key)
        if not analysis.get("enabled"):
            raise FileNotFoundError(
                f"{fastcomm_dir} is not a finished fastComm run: {analysis.get('message')}")
        states = list(meta.get("scalable_viewer", {}).get("states") or [])
        if a.allow_stale_fastcomm:
            assets["fastcomm_provenance"] = {"verified": False, "waived": "--allow-stale-fastcomm"}
            print(f"[WARN] fastComm provenance NOT checked (--allow-stale-fastcomm). The scores "
                  f"may describe cells this bundle does not serve.", flush=True)
        else:
            provenance = verify_fastcomm_matches_bundle(
                fastcomm_dir, meta, states, cluster_key,
                populations=analysis.get("populations"))
            assets["fastcomm_provenance"] = provenance
            print(f"[assets] fastComm provenance verified: same {cluster_key} key, "
                  f"{provenance['n_cells']} cells as the bundle, "
                  f"{provenance['n_states_scored']} of {provenance['bundle_n_states']} "
                  f"cell states scored", flush=True)
        assets["fastcomm_analysis"] = analysis
        print(f"[assets] fastComm: enabled={analysis.get('enabled')} "
              f"rows={analysis.get('n_rows')} states={len(analysis.get('populations', []))} "
              f"from {fastcomm_dir}", flush=True)

    manifest = out / f"{a.prefix}_assets.json"
    with open(manifest, "w") as fh:
        json.dump(assets, fh, indent=1)
    print(f"[assets] manifest {manifest}", flush=True)

    # ---- the gate that stops a view being dropped without anyone noticing --------
    waived: Dict[str, str] = {}
    if a.no_marker_heatmap:
        waived["marker_heatmap"] = "--no-marker-heatmap"
    if a.no_fastcomm:
        waived["fastcomm"] = "--no-fastcomm"
    if a.no_goelite:
        waived["goelite"] = "--no-goelite"
    if a.skip_networks:
        waived["marker_networks"] = "--skip-networks"
        waived["diff_networks"] = "--skip-networks"
    failures = check_completeness(assets, waived=waived)
    if failures:
        print(f"[assets] FAILED: {len(failures)} viewer input(s) absent and not waived: "
              f"{failures}. The manifest above was written, so nothing is lost, but the viewer "
              f"would serve those views empty. Supply the input, or pass its opt-out flag on "
              f"purpose.", flush=True)
        return 2
    print("[assets] every viewer input is present or explicitly waived", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
