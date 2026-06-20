"""UniProt isoform -> ENST(by sequence) -> Ensembl91 structure + genomic coords + domain-contact PPIs.

Implements the requested UniProt layer, reusing the Daedalus-processed UniProt resources (no download):
  - ``uniprot_isoform_proteins.tsv``  : human isoform protein SEQUENCES (canonical + alternative)
  - ``uniprot_isoform_features.tsv``   : per-isoform domain / interface feature counts & spans
                                         (role_ppi_interface_*, role_dna_interface_*, domain_count, ...)
  - ``uniprot_entries.tsv``            : per-accession curated interaction_partners + ensembl xrefs

Pipeline:
  1. Crosswalk each UniProt isoform to ENST by PROTEIN-SEQUENCE identity vs GENCODE (UniProt's own
     ENST xref is unreliable at isoform level), then to the Ensembl91 structure + genomic splice
     coords via the reference tables -> ``uniprot_isoform_map.tsv``.
  2. Attach domain-contact features + the gene's UniProt-curated PPI partners, and infer per-isoform
     PPI impact from PPI-INTERFACE RETENTION vs the gene's canonical isoform -> ``uniprot_isoform_function.tsv``.

Evidence labelling (high-confidence discipline): the PPI partner LIST is UniProt-curated (measured,
gene-level); the per-isoform GAIN/LOSS call is ``domain_inferred`` (from interface retention) and is
kept in its own column, never conflated with the measured atlas/paper1 PPIs.
"""

import os
import re
import logging

import pandas as pd

from . import config
from .crosswalk import sequence_map, reference_structures

logger = logging.getLogger(__name__)

_SPLIT = re.compile(r"[;,|]")


def load_human_isoforms():
    df = pd.read_csv(config.require(config.UNIPROT_ISOFORM_PROTEINS, what="uniprot isoform proteins"),
                     sep="\t", dtype=str).fillna("")
    return df[df["species"] == "human"].copy()


def load_isoform_features():
    df = pd.read_csv(config.require(config.UNIPROT_ISOFORM_FEATURES, what="uniprot isoform features"),
                     sep="\t", dtype=str).fillna("")
    return df[df["species"] == "human"].copy() if "species" in df.columns else df


def load_interaction_partners():
    """{primary_accession -> [partner tokens]} from uniprot_entries.interaction_partners (curated)."""
    df = pd.read_csv(config.require(config.UNIPROT_ENTRIES, what="uniprot entries"),
                     sep="\t", dtype=str).fillna("")
    out = {}
    for _, r in df.iterrows():
        raw = str(r.get("interaction_partners", "")).strip()
        if not raw or raw.lower() == "nan":
            continue
        parts = [p.strip() for p in _SPLIT.split(raw) if p.strip() and p.strip().lower() != "nan"]
        if parts:
            out[r["primary_accession"]] = sorted(set(parts))
    logger.info("[iso2function.uniprot] entries with curated interaction partners: %d", len(out))
    return out


def build_uniprot_crosswalk(data_dir=None, min_identity=99.0):
    """Sequence-match human UniProt isoforms to ENST -> Ensembl91 structure + genomic coords. Writes
    ``uniprot_isoform_map.tsv``."""
    data_dir = data_dir or config.DATA_DIR
    os.makedirs(data_dir, exist_ok=True)
    iso = load_human_isoforms()
    query = {r["isoform_id"]: r["protein_sequence"] for _, r in iso.iterrows() if r["protein_sequence"]}
    gene_of_query = {r["isoform_id"]: str(r["gene_name"]).upper() for _, r in iso.iterrows()}
    meta = {r["isoform_id"]: {"gene_symbol": r["gene_name"], "accession": r["primary_accession"],
                              "is_canonical": r["is_canonical"]} for _, r in iso.iterrows()}

    gencode = sequence_map.load_gencode_protein_index()
    gene_of_target = {prot: m["gene"].upper() for prot, m in gencode.items()}
    res = sequence_map.crosswalk_sequences(query, gencode, min_identity=min_identity,
                                           gene_of_query=gene_of_query, gene_of_target=gene_of_target)

    e2s = reference_structures.load_enst_to_structure()
    coords = reference_structures.StructureCoords()
    rows = []
    for r in res:
        iid = r["query_id"]
        ensts = r.get("ensts") or ([r.get("enst")] if r.get("enst") else [])
        enst, struct, ensg = reference_structures.resolve_structure([str(e).split(".")[0] for e in ensts], e2s)
        rows.append({
            "isoform_id": iid, **meta.get(iid, {}),
            "matched": r["matched"], "identity": r["identity"],
            "match_mode": ("exact" if r["identity"] == 100.0 else
                           (f"clone>={min_identity:g}" if r["matched"] else "unmatched")),
            "enst": enst, "ensp": r.get("ensp", ""), "ensg": ensg,
            "structure": struct, "has_structure": bool(struct),
            "genomic_junctions": (coords.junction_string(ensg, struct) or "") if struct else "",
        })
    coords.close()
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(data_dir, "uniprot_isoform_map.tsv"), sep="\t", index=False)
    logger.info("[iso2function.uniprot] crosswalk: %d isoforms, %d ENST-mapped, %d -> structure",
                len(out), int(out["matched"].sum()), int(out["has_structure"].sum()))
    return out


def build_uniprot_function(data_dir=None):
    """Join crosswalk + isoform domain features + gene-level UniProt-curated PPIs, and infer per-isoform
    PPI impact from PPI-interface retention vs the gene's canonical isoform. Writes
    ``uniprot_isoform_function.tsv``."""
    data_dir = data_dir or config.DATA_DIR
    mp_path = os.path.join(data_dir, "uniprot_isoform_map.tsv")
    if not os.path.exists(mp_path):
        build_uniprot_crosswalk(data_dir)
    mp = pd.read_csv(mp_path, sep="\t", dtype=str).fillna("")
    # drop_duplicates: a few isoform_ids appear twice in the features table; without this, .loc returns
    # a Series (multi-row) and writes newline-containing cells. keep first -> clean scalar per isoform.
    feats = load_isoform_features().drop_duplicates("isoform_id").set_index("isoform_id")
    iso = load_human_isoforms().drop_duplicates("isoform_id").set_index("isoform_id")
    partners = load_interaction_partners()

    def _f(iid, col):
        return feats.loc[iid][col] if (iid in feats.index and col in feats.columns) else ""

    # canonical PPI-interface span per gene (for the retention comparison)
    canon_span = {}
    for iid, r in iso.iterrows():
        if str(r.get("is_canonical")) == "1":
            try:
                canon_span[str(r["gene_name"]).upper()] = float(_f(iid, "role_ppi_interface_span") or 0)
            except ValueError:
                pass

    rows = []
    for _, r in mp.iterrows():
        iid = r["isoform_id"]; gene = str(r.get("gene_symbol", "")).upper()
        acc = r.get("accession", "")
        ppi_span = _f(iid, "role_ppi_interface_span")
        try:
            ppi_span_f = float(ppi_span or 0)
        except ValueError:
            ppi_span_f = 0.0
        c_span = canon_span.get(gene)
        is_canon = str(r.get("is_canonical")) == "1"
        if is_canon or c_span is None:
            ppi_call = "reference" if is_canon else "no_canonical_ref"
        elif ppi_span_f < c_span:
            ppi_call = "predicted_ppi_loss"        # lost interface vs canonical -> partners at risk
        elif ppi_span_f > c_span:
            ppi_call = "predicted_ppi_gain"
        else:
            ppi_call = "predicted_retained"
        plist = partners.get(acc, [])
        rows.append({
            "isoform_id": iid, "gene_symbol": r.get("gene_symbol", ""), "accession": acc,
            "is_canonical": r.get("is_canonical", ""), "enst": r.get("enst", ""),
            "structure": r.get("structure", ""), "genomic_junctions": r.get("genomic_junctions", ""),
            "domain_count": _f(iid, "domain_count"),
            "ppi_interface_count": _f(iid, "role_ppi_interface_count"),
            "ppi_interface_span": ppi_span,
            "dna_interface_count": _f(iid, "role_dna_interface_count"),
            "zinc_finger_count": _f(iid, "zinc_finger_feature_count"),
            "dimerization_count": _f(iid, "role_dimerization_count"),
            "evidence": "domain_inferred",
            "ppi_impact_call": ppi_call,
            "n_uniprot_partners": len(plist),
            "uniprot_partners": ";".join(plist),
        })
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(data_dir, "uniprot_isoform_function.tsv"), sep="\t", index=False)
    logger.info("[iso2function.uniprot] uniprot_isoform_function.tsv: %d isoforms (%d w/ structure, "
                "%d predicted_ppi_loss, %d with curated partners)", len(out),
                int((out["structure"].astype(str) != "").sum()),
                int((out["ppi_impact_call"] == "predicted_ppi_loss").sum()),
                int((out["n_uniprot_partners"] > 0).sum()))
    return out


def build_all(data_dir=None):
    build_uniprot_crosswalk(data_dir)
    return build_uniprot_function(data_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config.ensure_dirs()
    build_all()
