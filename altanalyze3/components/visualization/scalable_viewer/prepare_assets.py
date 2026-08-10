"""Build the scALABLE side-artifacts a bundle needs, using scALABLE's own backends.

A precomputed bundle carries the matrix, the embedding, the markers and the DEG
tables. Three scALABLE views need artifacts a bundle does not carry:

  MarkerHeatmap     a genes x cells matrix npz
                    (format: marker_heatmap_h5ad.py:67 `_write_heatmap_cache`)
  MarkerNetwork     one NetPerspective interaction TSV per cell state
                    (generator: marker_heatmap_h5ad.py:757 `_export_marker_networks`)
  Cell communication  a fastComm scores TSV
                    (generator: altanalyze3.components.fastComm.cli `score`)

This module produces the first two by calling the same generators the cellHarmony
pipeline calls, and reads a finished fastComm run for the third. It implements no
scoring, no layout and no statistic of its own.

Two further views are fed from artifacts the cellHarmony DIFFERENTIAL already wrote,
and must never be recomputed here:

  Differential GO terms      <contrast>/GeneSetEnrichment/GOElite_<tag>.tsv
                             (writer: flask/pipeline.py:2507)
  Differential networks      <contrast>/interaction-plots/<population>_interactions.tsv
                             (writer: flask/pipeline.py:2431-2456)

Pass `--goelite-from` and `--diff-networks-from` to ingest those files. The older
`--goelite` flag recomputes enrichment from the bundle DEG table and does NOT reproduce
the pipeline; it is kept only for a bundle with no differential run behind it.

Run:
  PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3 \
  /opt/homebrew/opt/python@3.11/bin/python3.11 \
    -m altanalyze3.components.visualization.scalable_viewer.prepare_assets \
      --bundle-dir <dir> --prefix <prefix> --out <asset dir> \
      [--marker-gct <per-cell GCT>] [--markers-tsv <marker table>] \
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
    ap.add_argument("--marker-gct", default=None, help="per-cell marker GCT to convert to the npz")
    ap.add_argument("--markers-tsv", default=None, help="marker table with Gene/Fold/cluster")
    ap.add_argument("--network-markers-tsv", default=None,
                    help="redundant marker table (250 genes per state) used for the marker "
                         "networks; defaults to *_redundant_markers.tsv beside --markers-tsv")
    ap.add_argument("--order-tsv", default=None, help="canonical cell-state order TSV")
    ap.add_argument("--fastcomm-dir", default=None, help="finished fastComm run directory")
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
    a = ap.parse_args(argv)

    if a.goelite and a.goelite_from:
        ap.error("--goelite recomputes enrichment and --goelite-from ingests the pipeline's own; "
                 "pass only one")

    paths = B.BundlePaths(a.bundle_dir, a.prefix)
    meta = B.read_metadata(paths.metadata)
    cluster_key = meta.get("cluster_key", "cell_state")
    markers_tsv = a.markers_tsv or meta.get("scalable_viewer", {}).get("markers", {}).get("source")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    assets: Dict[str, object] = {"bundle_dir": paths.bundle_dir, "prefix": paths.prefix}

    if a.marker_gct:
        npz = str(out / f"{a.prefix}_marker_heatmap_cache.npz")
        info = bundle_meta.gct_to_marker_cache(a.marker_gct, npz)
        assets["heatmap_cache"] = npz
        assets["heatmap_cache_info"] = info
        print(f"[assets] marker heatmap npz {npz}: {info['n_rows']} genes x {info['n_cols']} cells",
              flush=True)

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

    if a.fastcomm_dir:
        assets["fastcomm_dir"] = a.fastcomm_dir
        analysis = bundle_meta.fastcomm_analysis_from_dir(a.fastcomm_dir, cluster_key,
                                                         sample_key=a.fastcomm_sample_key)
        assets["fastcomm_analysis"] = analysis
        print(f"[assets] fastComm: enabled={analysis.get('enabled')} "
              f"rows={analysis.get('n_rows')} states={len(analysis.get('populations', []))}", flush=True)

    manifest = out / f"{a.prefix}_assets.json"
    with open(manifest, "w") as fh:
        json.dump(assets, fh, indent=1)
    print(f"[assets] manifest {manifest}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
