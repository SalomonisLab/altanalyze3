"""Utilities for generating simplified gene regulatory network visualizations.

This module centralizes logic that was previously embedded inside ad-hoc
figure generation scripts.  It accepts differential expression results,
filters them against a curated interaction table, and exports network plots
plus tabular summaries describing the inferred interactions.

Example usage
-------------

>>> from altanalyze3.components.visualization import NetPerspective as npv
>>> interactions = npv.load_interaction_data()
>>> outputs = npv.generate_network_for_genes(deg_df, interactions, "output/network")

The call above produces a PDF, PNG, and TSV inside ``output``.  The function
raises :class:`NetworkGenerationError` when a plot cannot be generated because
no interactions remain after filtering; callers may catch the exception to
silence non-critical failures.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd

__all__ = [
    "INTERACTION_TABLE",
    "NetworkGenerationError",
    "sanitize_filename",
    "safe_component",
    "load_interaction_data",
    "generate_network_for_genes",
]


INTERACTION_TABLE = os.path.join(
    os.path.dirname(__file__),
    "interactions",
    "Hs_Ensembl-TF-BioGRID-Pathway.txt",
)

# Allow default loading of both human and mouse interaction resources when available.
INTERACTION_TABLES = [
    INTERACTION_TABLE,
    os.path.join(
        os.path.dirname(__file__),
        "interactions",
        "Mm_Ensembl-TF-BioGRID-Pathway.txt",
    ),
]


class NetworkGenerationError(RuntimeError):
    """Raised when a network cannot be generated for the provided genes."""


def sanitize_filename(path: str, max_length: int = 200) -> str:
    """Return a filesystem-safe version of ``path``.

    Only the filename component is sanitized; directory segments are preserved.
    """

    directory, filename = os.path.split(path)
    filename = re.sub(r"[^A-Za-z0-9._\-]", "_", filename)
    filename = filename.strip()
    if not filename:
        filename = "network"
    if len(filename) > max_length:
        filename = filename[:max_length]
    return os.path.join(directory, filename)


def safe_component(value: str, fallback: str = "value") -> str:
    """Sanitize an arbitrary string for inclusion in a filename component."""

    if value is None:
        return fallback
    cleaned = re.sub(r"[^A-Za-z0-9._\-]", "_", str(value)).strip("._-")
    return cleaned or fallback


def _load_single_interaction_table(path: str) -> pd.DataFrame:
    """Load a single interaction table from ``path`` with validation."""

    df = pd.read_csv(path, sep="\t")
    required = {"Symbol1", "InteractionType", "Symbol2"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Interaction file '{path}' is missing required columns: {sorted(missing)}")

    df["Symbol1"] = df["Symbol1"].astype(str)
    df["Symbol2"] = df["Symbol2"].astype(str)
    df["InteractionType"] = df["InteractionType"].astype(str)
    if "Source" not in df.columns:
        df["Source"] = ""

    return df


def load_interaction_data(interaction_path: Optional[str] = None) -> pd.DataFrame:
    """
    Load curated interaction tables used for GRN visualisations.

    When ``interaction_path`` is ``None`` the loader looks for both human (Hs_)
    and mouse (Mm_) resources in the ``interactions`` directory and concatenates
    any that are present.  Callers may pass either a single file path or an
    iterable of paths to override this behaviour explicitly.
    """

    if interaction_path is None:
        paths = INTERACTION_TABLES
    elif isinstance(interaction_path, (list, tuple, set)):
        paths = list(interaction_path)
    else:
        paths = [interaction_path]

    loaded_tables = []
    for path in paths:
        if not os.path.exists(path):
            continue
        loaded_tables.append(_load_single_interaction_table(path))

    if not loaded_tables:
        raise FileNotFoundError(
            "No interaction tables were found. Checked paths: "
            + ", ".join(paths)
        )

    df = pd.concat(loaded_tables, axis=0, ignore_index=True)
    return df


@dataclass
class _PreparedGenes:
    genes: pd.DataFrame
    fold_change_map: Dict[str, float]


def _prepare_gene_table(
    gene_table: pd.DataFrame,
    gene_column: str,
    fold_change_column: str,
    pval_column: Optional[str],
    max_genes: Optional[int],
) -> _PreparedGenes:
    """Normalise the gene statistics table used for plotting."""

    if gene_column not in gene_table.columns:
        gene_table = gene_table.reset_index().rename(columns={gene_table.index.name or "index": gene_column})

    working = gene_table.copy()
    working = working.dropna(subset=[gene_column]).copy()
    working[gene_column] = working[gene_column].astype(str)

    # Ensure fold change column exists and numeric
    if fold_change_column not in working.columns:
        raise ValueError(f"Gene table is missing required column '{fold_change_column}'.")
    working[fold_change_column] = pd.to_numeric(working[fold_change_column], errors="coerce").fillna(0.0)

    ascending = []
    sort_cols: List[str] = []
    if pval_column and pval_column in working.columns:
        working[pval_column] = pd.to_numeric(working[pval_column], errors="coerce").fillna(1.0)
        sort_cols.append(pval_column)
        ascending.append(True)

    sort_cols.append(fold_change_column)
    ascending.append(False)  # Sort by absolute fold when used below

    # Sort by p-value (if provided) then by absolute log fold change
    working = working.sort_values(
        by=sort_cols,
        ascending=ascending,
        key=lambda s: np.abs(s) if s.name == fold_change_column else s,
    )

    if max_genes and len(working) > max_genes:
        working = working.head(max_genes)

    fc_map = working.set_index(gene_column)[fold_change_column].to_dict()
    working = working.drop_duplicates(subset=[gene_column])

    return _PreparedGenes(genes=working, fold_change_map=fc_map)


def _filter_interactions(
    interaction_data: pd.DataFrame,
    genes: Iterable[str],
    *,
    exclude_biogrid: bool,
) -> pd.DataFrame:
    """Return the subset of interactions connecting the provided genes."""

    genes = set(str(g) for g in genes)
    if not genes:
        return pd.DataFrame(columns=interaction_data.columns)

    subset = interaction_data[
        interaction_data["Symbol1"].isin(genes) & interaction_data["Symbol2"].isin(genes)
    ].copy()

    subset = subset[subset["Symbol1"] != subset["Symbol2"]]

    if exclude_biogrid and "Source" in subset.columns:
        mask = ~subset["Source"].str.contains("BioGRID", na=False)
        subset = subset[mask]

    return subset


def generate_network_for_genes(
    gene_table: pd.DataFrame,
    interaction_data: pd.DataFrame,
    output_prefix: str,
    *,
    gene_column: str = "gene",
    fold_change_column: str = "log2fc",
    pval_column: Optional[str] = None,
    edge_color_map: Optional[Dict[str, str]] = None,
    exclude_biogrid: bool = True,
    max_genes: Optional[int] = 75,
) -> List[str]:
    """Generate PDF/PNG network plots and a TSV summary for the given genes.

    Returns a list of file paths written to disk.  Raises
    :class:`NetworkGenerationError` when no network could be generated.
    """

    prepared = _prepare_gene_table(
        gene_table,
        gene_column=gene_column,
        fold_change_column=fold_change_column,
        pval_column=pval_column,
        max_genes=max_genes,
    )

    filtered = _filter_interactions(
        interaction_data,
        prepared.genes[gene_column],
        exclude_biogrid=exclude_biogrid,
    )

    if filtered.empty or filtered["Symbol1"].nunique() + filtered["Symbol2"].nunique() < 2:
        raise NetworkGenerationError("No qualifying interactions among the supplied genes.")

    try:
        from igraph import Graph, plot as igraph_plot
    except ImportError as exc:  # pragma: no cover - depends on optional igraph
        raise ImportError(
            "python-igraph is required to create gene regulatory network plots."
        ) from exc

    edge_color_map = edge_color_map or {
        "transcriptional_target": "red",
        "Arrow": "grey",
        "Tbar": "blue",
    }

    nodes = sorted(set(filtered["Symbol1"]) | set(filtered["Symbol2"]))
    g = Graph(directed=True)
    g.add_vertices(nodes)
    g.vs["name"] = nodes

    edges = list(zip(filtered["Symbol1"], filtered["Symbol2"]))
    g.add_edges(edges)
    g.es["color"] = [edge_color_map.get(t, "gray") for t in filtered["InteractionType"]]

    fc_map = prepared.fold_change_map
    g.vs["color"] = [
        "lightcoral" if fc_map.get(name, 0) > 0 else
        "deepskyblue" if fc_map.get(name, 0) < 0 else
        "gray"
        for name in g.vs["name"]
    ]
    g.vs["label"] = g.vs["name"]

    layout = g.layout("fr")

    pdf_path = sanitize_filename(f"{output_prefix}.pdf")
    png_path = sanitize_filename(f"{output_prefix}.png")
    os.makedirs(os.path.dirname(pdf_path) or ".", exist_ok=True)

    igraph_plot(
        g,
        target=pdf_path,
        layout=layout,
        vertex_size=20,
        bbox=(1400, 1400),
        margin=50,
        vertex_label_size=11,
        edge_arrow_size=0.5,
        vertex_frame_width=0,
    )

    igraph_plot(
        g,
        target=png_path,
        layout=layout,
        vertex_size=20,
        bbox=(1400, 1400),
        margin=50,
        vertex_label_size=11,
        edge_arrow_size=0.5,
        vertex_frame_width=0,
    )

    # Create a concise summary table describing the downstream targets by sign.
    summary = filtered.copy()
    summary["Symbol2_LogFC"] = summary["Symbol2"].map(fc_map).fillna(0.0)
    summary["Direction"] = summary["Symbol2_LogFC"].apply(lambda x: "up" if x > 0 else ("down" if x < 0 else "neutral"))
    summary = (
        summary.groupby(["Symbol1", "InteractionType", "Direction"])["Symbol2"]
        .apply(lambda vals: "|".join(sorted(set(vals))))
        .reset_index()
    )

    summary_path = sanitize_filename(f"{output_prefix}_interactions.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)

    return [pdf_path, png_path, summary_path]
