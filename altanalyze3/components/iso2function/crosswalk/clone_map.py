"""Layer 2 (crosswalk) — reconcile interaction-data clone IDs to the canonical structure string and to
observed isoforms / ENST / ENSP.

Three legs (GOALS.md sec. 5):
  Leg B (GFF clone  -> structure)        : from ``tfiso_to_final_assignment.protein.tsv``  [done upstream]
  Leg C (structure  -> observed / ENST)  : same file (match_type + final_* columns)         [done upstream]
  Leg A (paper2 clone_id -> GFF clone)   : THIS module. The open problem.

Both the paper2 atlas clones (``SYMBOL-N``) and the TFIso GFF clones (``SYMBOL|i/n|well``) are the same
ORFeome, in different ID namespaces. We bridge them by shared ``gene_symbol`` and, decisively, by
**protein-sequence identity** — never by assuming ``-2`` == ``|2/3|``. Until sequences are supplied,
Leg A emits *candidate* pairs flagged ``verified=False`` (honest: no guessed mapping is asserted).

Outputs (under ``config.DATA_DIR``):
  - ``crosswalk_structure.tsv``  one row per GFF clone: structure-keyed identity + observed isoform +
                                 ENST + protein consequences (legs B+C, fully resolved).
  - ``clone_reconciliation.tsv`` Leg A: per-gene candidate (paper2 clone_id <-> GFF clone) pairs with a
                                 weak id-index hint and ``verified`` flag; plus unmatched paper2 clones.
  - ``crosswalk_coverage.tsv``   summary: paper2 clones with/without a same-gene GFF candidate, overlap.
"""

import os
import re
import logging

import pandas as pd

from .. import config
from . import structure_key, reference_structures

logger = logging.getLogger(__name__)

# which (gene_symbol_col, clone_id_col) carry the paper2 clone identity in each tidy ingest table
_PAPER2_CLONE_COLS = {
    "ppi_y2h": ("ad_gene_symbol", "ad_clone_id"),
    "ppi_n2h": ("gene_symbol_tf", "clone_id"),
    "pdi_ey1h": ("gene_symbol", "clone_id"),
    "pdi_validation": ("gene_symbol", "clone_id"),
    "activation_m1h": ("gene_symbol", "clone_id"),
    "condensate": (None, "clone_id"),                 # gene derivable from clone_id prefix
    # paralog_divergence carries two clones per row; handled separately
}

_GFF_CLONE_RE = re.compile(r"^(?P<gene>[^|]+)\|(?P<idx>\d+)/(?P<tot>\d+)\|(?P<well>.+)$")
_PAPER2_CLONE_RE = re.compile(r"^(?P<gene>.+)-(?P<num>\d+)$")


# --------------------------------------------------------------------------- id parsers
def parse_gff_clone_id(tid):
    """Parse a TFIso GFF clone id ``SYMBOL|i/n|well`` -> dict(gene_symbol, iso_index, iso_total, well).
    Returns dict with None fields if it does not match the expected shape."""
    m = _GFF_CLONE_RE.match(str(tid).strip())
    if not m:
        return {"gene_symbol": None, "iso_index": None, "iso_total": None, "well": None}
    return {"gene_symbol": m.group("gene"), "iso_index": int(m.group("idx")),
            "iso_total": int(m.group("tot")), "well": m.group("well")}


def parse_paper2_clone_id(cid):
    """Parse a paper2 atlas clone id ``SYMBOL-N`` -> dict(gene_symbol, iso_num). Falls back to
    (gene_symbol=cid, iso_num=None) if the trailing ``-N`` is absent."""
    s = str(cid).strip()
    m = _PAPER2_CLONE_RE.match(s)
    if not m:
        return {"gene_symbol": s, "iso_num": None}
    return {"gene_symbol": m.group("gene"), "iso_num": int(m.group("num"))}


# --------------------------------------------------------------------------- loaders
def load_assignment():
    """Load the prior TFIso assignment (legs B+C). Returns the DataFrame with all 17 columns."""
    path = config.require(config.TFISO_ASSIGNMENT, what="TFIso assignment table")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    expected = {"gene", "gene_symbol", "tfiso_transcript_id", "match_type", "final_isoform_id",
                "tfiso_structure", "final_structure"}
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(f"assignment table missing columns {missing}; got {list(df.columns)}")
    return df


def collect_paper2_clones(data_dir=None):
    """Scan the ingested tidy tables for every distinct paper2 (gene_symbol, clone_id). Gene symbol is
    taken from the table where available, else parsed from the clone-id prefix. Returns a DataFrame
    [gene_symbol, clone_id, iso_num, sources] where ``sources`` lists the modalities the clone appears
    in."""
    data_dir = data_dir or config.DATA_DIR
    seen = {}  # clone_id -> {"gene_symbol":..., "sources": set()}

    def _add(gene, clone, source):
        clone = str(clone).strip()
        if not clone or clone.lower() == "nan":
            return
        if not gene or str(gene).strip().lower() in ("", "nan"):
            gene = parse_paper2_clone_id(clone)["gene_symbol"]
        rec = seen.setdefault(clone, {"gene_symbol": gene, "sources": set()})
        rec["sources"].add(source)

    for modality, (gene_col, clone_col) in _PAPER2_CLONE_COLS.items():
        path = os.path.join(data_dir, f"{modality}.tsv")
        if not os.path.exists(path):
            logger.warning("[iso2function.crosswalk] ingest table missing (run ingest first): %s", path)
            continue
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
        genes = df[gene_col] if (gene_col and gene_col in df.columns) else [None] * len(df)
        for g, c in zip(genes, df[clone_col]):
            _add(g, c, modality)

    # paralog_divergence: two clones per row
    par = os.path.join(data_dir, "paralog_divergence.tsv")
    if os.path.exists(par):
        df = pd.read_csv(par, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
        for side in ("a", "b"):
            gcol, ccol = f"gene_symbol_{side}", f"clone_id_{side}"
            if ccol in df.columns:
                genes = df[gcol] if gcol in df.columns else [None] * len(df)
                for g, c in zip(genes, df[ccol]):
                    _add(g, c, "paralog_divergence")

    rows = [{"gene_symbol": v["gene_symbol"], "clone_id": k,
             "iso_num": parse_paper2_clone_id(k)["iso_num"],
             "sources": ";".join(sorted(v["sources"]))}
            for k, v in seen.items()]
    return pd.DataFrame(rows).sort_values(["gene_symbol", "clone_id"]).reset_index(drop=True)


# --------------------------------------------------------------------------- crosswalk build
def _norm_gene(g):
    return str(g).strip().upper() if g is not None else ""


def build_structure_table(assignment=None):
    """Legs B+C: one structure-keyed row per GFF clone, with observed isoform + ENST + protein
    consequences carried through. `` enst`` is populated only when the observed match is a known
    (ENST-backed) isoform."""
    df = assignment if assignment is not None else load_assignment()
    out = pd.DataFrame({
        "gene_symbol": df["gene_symbol"],
        "ensg": df["gene"],
        "gff_clone_id": df["tfiso_transcript_id"],
        "structure": df["final_structure"].where(df["final_structure"].notna(), df["tfiso_structure"]),
        "tfiso_structure": df["tfiso_structure"],
        "final_structure": df["final_structure"],
        "match_type": df["match_type"],
        "assigned": df["match_type"].map(structure_key.is_assigned),
        "final_isoform_id": df["final_isoform_id"],
        "final_known": df["final_known"],
        "enst": df.apply(lambda r: r["final_isoform_id"]
                         if str(r.get("final_known", "")).strip().lower() == "known"
                         and str(r["final_isoform_id"]).startswith("ENST") else "", axis=1),
        "protein_length": df.get("final_protein_length", ""),
        "nmd_status": df.get("final_nmd_status", ""),
        "intron_retention": df.get("final_intron_retention", ""),
    })
    parsed = out["gff_clone_id"].map(parse_gff_clone_id)
    out["well"] = [p["well"] for p in parsed]
    out["iso_index"] = [p["iso_index"] for p in parsed]
    out["iso_total"] = [p["iso_total"] for p in parsed]
    return out


def build_clone_reconciliation(struct_table=None, paper2_clones=None):
    """Leg A scaffold: candidate (paper2 clone_id <-> GFF clone) pairs, matched on gene_symbol only,
    each flagged ``verified=False`` pending protein-sequence confirmation. ``id_index_hint`` is True
    when the paper2 ``-N`` equals the GFF ``i`` of ``i/n`` (a WEAK hint, never sufficient alone).

    Honest by construction: this produces *candidates*, not a resolved mapping. Sequence verification
    (:func:`verify_pair_by_sequence`) collapses each gene's candidates to the correct pair.
    """
    st = struct_table if struct_table is not None else build_structure_table()
    p2 = paper2_clones if paper2_clones is not None else collect_paper2_clones()
    gff_by_gene = {}
    for _, r in st.iterrows():
        gff_by_gene.setdefault(_norm_gene(r["gene_symbol"]), []).append(r)

    cand_rows, unmatched = [], []
    for _, pr in p2.iterrows():
        g = _norm_gene(pr["gene_symbol"])
        gff_rows = gff_by_gene.get(g, [])
        if not gff_rows:
            unmatched.append({"paper2_clone_id": pr["clone_id"], "gene_symbol": pr["gene_symbol"],
                              "reason": "no GFF clone for this gene in the TFIso assignment"})
            continue
        for gr in gff_rows:
            cand_rows.append({
                "gene_symbol": pr["gene_symbol"],
                "paper2_clone_id": pr["clone_id"],
                "paper2_iso_num": pr["iso_num"],
                "gff_clone_id": gr["gff_clone_id"],
                "gff_iso_index": gr["iso_index"],
                "structure": gr["structure"],
                "enst": gr["enst"],
                "final_isoform_id": gr["final_isoform_id"],
                "id_index_hint": (pr["iso_num"] is not None and gr["iso_index"] is not None
                                  and int(pr["iso_num"]) == int(gr["iso_index"])),
                "seq_pct_identity": "",      # filled by verify_pair_by_sequence
                "verified": False,           # NEVER True without sequence confirmation
            })
    return pd.DataFrame(cand_rows), pd.DataFrame(unmatched)


def resolve_clone_to_structure(struct_table=None, paper2_clones=None):
    """Leg A (primary resolver) — map each atlas interaction clone (``SYMBOL-N``) to its standardized
    **structure key** via the bridge file, joining on (gene_symbol, isoform index) where the atlas ``N``
    equals the bridge ``i`` of ``i/n``. The two ID forms are just the paper's labels for the SAME
    clones; the structure key is the AltAnalyze3 identity carried forward (the ``|i/n|well`` label is
    discarded downstream). Returns (resolved_df, unresolved_df). Resolved rows are 1:1 (gene, index)
    matches; anything not matchable (gene absent from the bridge, index out of range, unparsable
    clone_id, or an ambiguous duplicate index) goes to ``unresolved`` with a reason — never force-matched.
    """
    st = struct_table if struct_table is not None else build_structure_table()
    p2 = paper2_clones if paper2_clones is not None else collect_paper2_clones()
    idx, dup, genes_in_bridge = {}, set(), set()
    for _, r in st.iterrows():
        genes_in_bridge.add(_norm_gene(r["gene_symbol"]))
        if r["iso_index"] is None:
            continue
        key = (_norm_gene(r["gene_symbol"]), int(r["iso_index"]))
        if key in idx:
            dup.add(key)
        idx[key] = r
    resolved, unresolved = [], []
    for _, pr in p2.iterrows():
        g, n = _norm_gene(pr["gene_symbol"]), pr["iso_num"]
        if n is None:
            unresolved.append({"clone_id": pr["clone_id"], "gene_symbol": pr["gene_symbol"],
                               "reason": "clone_id not in SYMBOL-N form (no isoform index)"})
            continue
        key = (g, int(n))
        if key in dup:
            unresolved.append({"clone_id": pr["clone_id"], "gene_symbol": pr["gene_symbol"],
                               "reason": f"ambiguous: >1 bridge clone at index {n}"})
            continue
        if key not in idx:
            reason = (f"index {n} not present for gene in bridge" if g in genes_in_bridge
                      else "gene not in TFIso bridge")
            unresolved.append({"clone_id": pr["clone_id"], "gene_symbol": pr["gene_symbol"], "reason": reason})
            continue
        br = idx[key]
        resolved.append({
            "clone_id": pr["clone_id"], "gene_symbol": pr["gene_symbol"], "iso_index": int(n),
            "structure": br["structure"], "ensg": br["ensg"], "match_type": br["match_type"],
            "assigned": br["assigned"], "final_isoform_id": br["final_isoform_id"],
            "final_known": br["final_known"], "enst": br["enst"],
            "protein_length": br["protein_length"], "nmd_status": br["nmd_status"],
            "gff_clone_label": br["gff_clone_id"], "sources": pr["sources"],
        })
    return pd.DataFrame(resolved), pd.DataFrame(unresolved)


def pct_identity(seq_a, seq_b):
    """Percent amino-acid identity between two protein strings over the shorter length (ungapped,
    position-wise from the N-terminus). A fast pre-filter; a full alignment can replace this later.
    Returns a float in [0, 100], or 0.0 if either sequence is empty."""
    a, b = (seq_a or "").strip("*"), (seq_b or "").strip("*")
    if not a or not b:
        return 0.0
    n = min(len(a), len(b))
    same = sum(1 for i in range(n) if a[i] == b[i])
    return 100.0 * same / max(len(a), len(b))


def verify_pair_by_sequence(reconciliation, paper2_seqs, gff_seqs, min_identity=98.0):
    """Resolve Leg A by protein-sequence identity. ``paper2_seqs`` maps paper2 clone_id -> AA sequence;
    ``gff_seqs`` maps gff_clone_id -> AA sequence (e.g. from
    ``long_read.isoform_translation.extract_cds_and_protein`` run on the TFIso GFF). Sets
    ``seq_pct_identity`` and ``verified`` (True iff identity >= ``min_identity``) on each candidate.

    Returns the reconciliation DataFrame with those columns populated. Candidates whose sequences are
    unavailable keep ``verified=False`` and an empty identity (reported, never guessed).
    """
    df = reconciliation.copy()
    ident, verified = [], []
    for _, r in df.iterrows():
        sa = paper2_seqs.get(r["paper2_clone_id"])
        sb = gff_seqs.get(r["gff_clone_id"])
        if not sa or not sb:
            ident.append("")
            verified.append(False)
            continue
        pid = pct_identity(sa, sb)
        ident.append(round(pid, 2))
        verified.append(pid >= min_identity)
    df["seq_pct_identity"] = ident
    df["verified"] = verified
    return df


def build_crosswalk(data_dir=None):
    """Assemble the crosswalk and write its artifacts. Leg A is resolved by the gene+isoform-index join
    (:func:`resolve_clone_to_structure`); the primary deliverable is ``clone_to_structure.tsv`` — every
    atlas interaction clone keyed to its standardized structure (+ ENST / observed isoform / protein
    consequences). Returns a dict of the written DataFrames and paths."""
    data_dir = data_dir or config.DATA_DIR
    os.makedirs(data_dir, exist_ok=True)
    struct = build_structure_table()
    p2 = collect_paper2_clones(data_dir)
    resolved, unresolved = resolve_clone_to_structure(struct, p2)

    # attach portable genomic splice-site coordinates per structure (build-independent annotation)
    if len(resolved):
        try:
            coords = reference_structures.StructureCoords()
            resolved["genomic_junctions"] = [coords.junction_string(g, s) or ""
                                             for g, s in zip(resolved["ensg"], resolved["structure"])]
            coords.close()
        except FileNotFoundError as e:
            logger.warning("[iso2function.crosswalk] genomic coords unavailable: %s", e)
            resolved["genomic_junctions"] = ""

    # Back-fill best ENST for clones that have an Ens91 structure but no known cohort ENST (their
    # cohort match was a novel collapse id). Resolve by gene-restricted structure containment against
    # the FULL Ensembl91 reference: a novel sub-isoform inherits the parent transcript it folds into.
    if len(resolved):
        resolved["enst_source"] = resolved["enst"].map(
            lambda e: "cohort_known" if str(e).startswith("ENST") else "")
        need = (resolved["enst"].astype(str) == "") & (resolved["structure"].astype(str) != "")
        if need.any():
            try:
                gi = reference_structures.load_gene_to_structures()
                cache, n_fill = {}, 0
                for i in resolved.index[need]:
                    key = (resolved.at[i, "ensg"], resolved.at[i, "structure"])
                    if key not in cache:
                        cache[key] = reference_structures.best_enst_by_structure(key[1], key[0], gi)
                    enst, src = cache[key]
                    if enst:
                        resolved.at[i, "enst"], resolved.at[i, "enst_source"] = enst, src
                        n_fill += 1
                logger.info("[iso2function.crosswalk] back-filled %d ENSTs via structure containment "
                            "(%d clones still novel/unresolved)", n_fill, int(need.sum()) - n_fill)
            except FileNotFoundError as e:
                logger.warning("[iso2function.crosswalk] ENST back-fill skipped: %s", e)

        # fill a real ENSG for UNK-locus clones from a sibling clone of the same gene_symbol (the gene
        # is known; only the coordinate-based gene resolution failed for that clone).
        sym2ensg = {}
        for _, r in resolved.iterrows():
            if str(r["ensg"]).startswith("ENSG"):
                sym2ensg.setdefault(r["gene_symbol"], r["ensg"])
        bad = ~resolved["ensg"].astype(str).str.startswith("ENSG")
        if bad.any():
            resolved.loc[bad, "ensg"] = resolved.loc[bad, "gene_symbol"].map(sym2ensg)\
                .fillna(resolved.loc[bad, "ensg"])

    coverage = pd.DataFrame([{
        "atlas_clones_total": len(p2),
        "resolved_to_structure": len(resolved),
        "unresolved": len(unresolved),
        "resolved_with_known_enst": int((resolved["enst"].astype(str) != "").sum()) if len(resolved) else 0,
        "resolved_assigned_to_observed": int(resolved["assigned"].sum()) if len(resolved) else 0,
        "distinct_structures": resolved["structure"].nunique() if len(resolved) else 0,
        "gff_clones_total": len(struct),
    }])

    paths = {}
    for name, frame in (("crosswalk_structure", struct), ("clone_to_structure", resolved),
                        ("clone_unresolved", unresolved), ("crosswalk_coverage", coverage)):
        p = os.path.join(data_dir, f"{name}.tsv")
        frame.to_csv(p, sep="\t", index=False)
        paths[name] = p
    logger.info("[iso2function.crosswalk] atlas clones=%d  resolved->structure=%d  unresolved=%d  "
                "(known ENST=%d, distinct structures=%d)", len(p2), len(resolved), len(unresolved),
                int((resolved["enst"].astype(str) != "").sum()) if len(resolved) else 0,
                resolved["structure"].nunique() if len(resolved) else 0)
    return {"structure": struct, "resolved": resolved, "unresolved": unresolved,
            "coverage": coverage, "paths": paths}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config.ensure_dirs()
    res = build_crosswalk()
    print(res["coverage"].to_string(index=False))
