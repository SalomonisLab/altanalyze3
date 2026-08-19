"""rna2grn.grn_io — parse the GRN connection-score matrix and match its
sample x cell-state columns to RNA pseudobulks.

The GRN matrix (tf_to_gene_connection_scores_*.txt) has columns:
    TF_CON, Gene, <sampleA_cellstate1>, <sampleA_cellstate2>, ...
Each data column is one (sample, cell-state) pseudobulk; each row is a
(TF -> target Gene) edge; the value is a connection score.

GRN column names join the sample and the cell-state with ``_`` and replace
spaces inside cell-state names with ``_`` as well, e.g.
    "1009_AfInv16_29M_Intermediate_Mono-1"
maps to sample "1009_AfInv16_29M" + cell-state "Intermediate Mono-1".

The RNA pseudobulk h5ad keys pseudobulks as ``"<Sample>|<CellState>"`` with the
original (space-containing) cell-state names, so the matching strategy is:

  1. Strip the longest cell-state suffix (from the known cell-state vocabulary,
     with spaces -> underscores) off each GRN column -> (sample_token, cellstate).
  2. Normalise/override the sample token to an RNA ``Sample`` value.
  3. Look up the RNA pseudobulk ``"<rna_sample>|<cellstate>"``.

This module is intentionally dependency-light (numpy/pandas only).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# Control-aggregate GRN samples that have no single RNA sample (they are built
# from an average of control-sample pseudobulks of the same cell state).
CONTROL_AGG_SAMPLES = ("Multiome_WT", "RC2_TEA")

# Manual GRN-sample -> RNA-sample overrides where automatic normalisation fails
# (the GRN file abbreviated or re-delimited the sample id).
SAMPLE_OVERRIDES: Dict[str, str] = {
    "287_AfN_P53F": "287_AfInv_P53F_N1c",
    "3114_EInv16_45F": "3114_EInv_45F",
    "1702_Et922_66F": "1702-Et922-66F",
    "98_31_AfInv16_34F": "98-31_AfInv16_34F",
    "5381_v2_AfInv16_P49M": "5381v2_AfInv16_P49M",
    # CCHMC AML CITE samples (replicate suffixes _1/_2 collapse to one merged
    # RNA pseudobulk; replicate GRN columns are averaged downstream).
    "AML_7_1": "AML-7_CCHMC",
    "AML_7_2": "AML-7_CCHMC",
    "AML_12_1": "AML-12_CCHMC",
    "AML_12_2": "AML-12_CCHMC",
    "AML_14": "AML-14_CCHMC",
    # "AML_13" intentionally absent: the only RNA "AML-13" is a different
    # Colorado cohort patient; the CCHMC AML-13 RNA is not in the atlas -> drop.
}

# GRN samples to drop outright (no defensible RNA partner).
DROP_SAMPLES = ("AML_13",)


def _norm_sample(s: str) -> str:
    return s.lower().replace("-", "_").replace("__", "_").strip("_")


@dataclass
class GrnMatrix:
    edges: pd.DataFrame            # columns: TF, Gene
    edge_ids: List[str]           # "TF|Gene"
    columns: List[str]            # original GRN data-column names
    values: np.ndarray            # (n_edges, n_columns) float32

    @property
    def n_edges(self) -> int:
        return self.values.shape[0]


def load_grn_matrix(path: str) -> GrnMatrix:
    """Load the GRN connection-score matrix file."""
    df = pd.read_csv(path, sep="\t")
    tf = df.iloc[:, 0].astype(str).values
    gene = df.iloc[:, 1].astype(str).values
    cols = df.columns[2:].tolist()
    vals = df[cols].to_numpy(dtype=np.float32)
    edges = pd.DataFrame({"TF": tf, "Gene": gene})
    edge_ids = [f"{t}|{g}" for t, g in zip(tf, gene)]
    return GrnMatrix(edges=edges, edge_ids=edge_ids, columns=cols, values=vals)


def parse_grn_columns(
    columns: Sequence[str], cellstate_vocab: Sequence[str]
) -> Dict[str, Tuple[Optional[str], Optional[str]]]:
    """Split each GRN column into (sample_token, cell_state).

    cellstate_vocab uses original (space-containing) names. Matching is done by
    longest cell-state suffix after replacing spaces with underscores.
    """
    cs_us = sorted(
        ((c.replace(" ", "_"), c) for c in cellstate_vocab), key=lambda x: -len(x[0])
    )
    out: Dict[str, Tuple[Optional[str], Optional[str]]] = {}
    for col in columns:
        match = (None, None)
        for us, orig in cs_us:
            if col == us:
                match = ("", orig)
                break
            if col.endswith("_" + us):
                match = (col[: len(col) - len(us) - 1], orig)
                break
        out[col] = match
    return out


def resolve_sample(sample_token: str, rna_samples: Sequence[str]) -> Optional[str]:
    """Map a GRN sample token to an RNA Sample value (or None)."""
    if sample_token in DROP_SAMPLES:
        return None
    if sample_token in SAMPLE_OVERRIDES:
        return SAMPLE_OVERRIDES[sample_token]
    if sample_token in rna_samples:
        return sample_token
    norm_map = {_norm_sample(s): s for s in rna_samples}
    if _norm_sample(sample_token) in norm_map:
        return norm_map[_norm_sample(sample_token)]
    return None


@dataclass
class MatchResult:
    # one row per matched GRN column
    table: pd.DataFrame                     # grn_column, sample_token, cell_state,
                                            # rna_sample, rna_pseudobulk, group
    unmatched: List[str] = field(default_factory=list)
    control_columns: List[str] = field(default_factory=list)


def classify_group(sample_token: str) -> str:
    if sample_token == "Multiome_WT":
        return "Multiome_control"
    if sample_token == "RC2_TEA":
        return "TEA_control"
    if sample_token.startswith("AML_"):
        return "AML_CCHMC"
    return "patient"


def match_columns(
    grn_columns: Sequence[str],
    cellstate_vocab: Sequence[str],
    rna_samples: Sequence[str],
    rna_pseudobulks: Sequence[str],
) -> MatchResult:
    """Build the GRN-column -> RNA-pseudobulk match table."""
    pbset = set(rna_pseudobulks)
    parsed = parse_grn_columns(grn_columns, cellstate_vocab)
    rows = []
    unmatched: List[str] = []
    control_cols: List[str] = []
    for col, (stok, cs) in parsed.items():
        if cs is None:
            unmatched.append(col)
            continue
        group = classify_group(stok)
        if stok in CONTROL_AGG_SAMPLES:
            control_cols.append(col)
            rows.append(
                dict(grn_column=col, sample_token=stok, cell_state=cs,
                     rna_sample=None, rna_pseudobulk=None, group=group)
            )
            continue
        rs = resolve_sample(stok, rna_samples)
        pb = f"{rs}|{cs}" if rs is not None else None
        if pb is not None and pb in pbset:
            rows.append(
                dict(grn_column=col, sample_token=stok, cell_state=cs,
                     rna_sample=rs, rna_pseudobulk=pb, group=group)
            )
        else:
            unmatched.append(col)
    table = pd.DataFrame(rows)
    return MatchResult(table=table, unmatched=unmatched, control_columns=control_cols)


# ---------------------------------------------------------------------------
# Lung (LungMAP IPF/control) reference
# ---------------------------------------------------------------------------
# The lung GRN file uses a different column grammar from the leukemia file:
#     "<Group>.<CellState>"   e.g. "Control.AT1", "IPF.KRT5-_KRT17"
# and the paired pseudobulk matrix keys columns as
#     "<CellState>|<Group>"   e.g. "AT1|Control", "KRT5-/KRT17+|IPF"
# The GRN writer mangles the cell-state name: it drops "+", then replaces " "
# and "/" with "_". ``mangle_cellstate`` reproduces that mangling exactly, so the
# match is an equality test on the mangled key, not a fuzzy suffix search.

LUNG_GROUPS = ("Control", "IPF")


def mangle_cellstate(cell_state: str) -> str:
    """Reproduce the lung GRN writer's cell-state name mangling.

    ``"Secretory - SCGB1A1+/MUC5B+"`` -> ``"Secretory_-_SCGB1A1_MUC5B"``
    """
    return str(cell_state).replace("+", "").replace(" ", "_").replace("/", "_")


def lung_pseudobulk_key(cell_state: str, group: str) -> str:
    """The GRN column name a ``"<CellState>|<Group>"`` pseudobulk maps to."""
    return f"{group}.{mangle_cellstate(cell_state)}"


def match_lung_columns(
    grn_columns: Sequence[str],
    rna_pseudobulks: Sequence[str],
    groups: Sequence[str] = LUNG_GROUPS,
) -> MatchResult:
    """Match lung GRN columns to ``"<CellState>|<Group>"`` pseudobulk columns.

    Returns a MatchResult whose table carries one row per matched GRN column with
    ``rna_sample`` = the group (Control / IPF) and ``cell_state`` = the original,
    un-mangled cell-state name from the RNA matrix.
    """
    gset = set(str(g) for g in groups)
    key_to_pb: Dict[str, List[str]] = {}
    pb_meta: Dict[str, Tuple[str, str]] = {}
    for pb in rna_pseudobulks:
        pb = str(pb)
        if "|" not in pb:
            continue
        cs, grp = pb.rsplit("|", 1)
        if grp not in gset:
            continue
        key_to_pb.setdefault(lung_pseudobulk_key(cs, grp), []).append(pb)
        pb_meta[pb] = (cs, grp)

    collisions = {k: v for k, v in key_to_pb.items() if len(v) > 1}
    if collisions:
        raise ValueError(
            "lung cell-state mangling is not injective; colliding GRN keys: "
            f"{ {k: v for k, v in list(collisions.items())[:5]} }"
        )

    rows: List[dict] = []
    unmatched: List[str] = []
    for col in grn_columns:
        col = str(col)
        hits = key_to_pb.get(col)
        if not hits:
            unmatched.append(col)
            continue
        pb = hits[0]
        cs, grp = pb_meta[pb]
        rows.append(
            dict(grn_column=col, sample_token=grp, cell_state=cs,
                 rna_sample=grp, rna_pseudobulk=pb, group=grp)
        )
    return MatchResult(table=pd.DataFrame(rows), unmatched=unmatched,
                       control_columns=[])
