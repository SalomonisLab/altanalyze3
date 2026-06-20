"""Resolve an ENST (and Ens91 structure) for EVERY Vidal2025 atlas clone, iteratively, via a cascade
that prefers transcripts with an Ensembl91 structure:

  1. CloneList authoritative ENST     -- SuppTable_CloneList.ensembl_transcript_ids (paper-assigned)
  2. existing resolved ENST           -- whatever clone_to_structure already had (cohort/protein/containment)
  3. CloneList CDS sequence match     -- translate cds_seq -> GENCODE protein -> ENST (novel clones)
  4. MANE Select (gene)               -- the gene's standard "main select" transcript as the reference ENST

At each step, when several ENSTs are candidates, the one present in the Ens91 reference
(ENST_reference_structures) is preferred so it carries an Ens91 structure. Updates clone_to_structure.tsv
in place (enst, enst_source, and structure for clones that had none). Idempotent.
"""

import os
import gzip
import logging
import re

import pandas as pd

from .. import config
from . import sequence_map, reference_structures

logger = logging.getLogger(__name__)
_SPLIT = re.compile(r"[;,|/ ]+")


def _load_clonelist():
    df = pd.read_csv(config.require(config.SUPP_CLONELIST, what="SuppTable_CloneList"),
                     sep="\t", dtype=str, keep_default_na=False, na_values=[""]).fillna("")
    out = {}
    for _, r in df.iterrows():
        ensts = [e.split(".")[0] for e in _SPLIT.split(str(r.get("ensembl_transcript_ids", "")))
                 if e.startswith("ENST")]
        out[r["clone_id"]] = {"ensts": ensts, "cds": r.get("cds_seq", ""),
                              "gene": r.get("gene_symbol", "")}
    logger.info("[iso2function.clonelist] CloneList: %d clones (%d with >=1 ENST)", len(out),
                sum(1 for v in out.values() if v["ensts"]))
    return out


def _load_mane():
    """{base_ENSG -> base_ENST} and {symbol -> base_ENST} for MANE Select transcripts."""
    path = config.require(config.MANE_SUMMARY, what="MANE summary")
    by_ensg, by_sym = {}, {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        h = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= h.get("MANE_status", 0) or f[h["MANE_status"]] != "MANE Select":
                continue
            enst = f[h["Ensembl_nuc"]].split(".")[0]
            ensg = f[h["Ensembl_Gene"]].split(".")[0]
            sym = f[h["symbol"]]
            if enst:
                by_ensg.setdefault(ensg, enst)
                by_sym.setdefault(sym, enst)
    logger.info("[iso2function.clonelist] MANE Select: %d genes (ENSG), %d symbols", len(by_ensg), len(by_sym))
    return by_ensg, by_sym


def _prefer_ens91(ensts, e2s):
    """Pick the first ENST that has an Ens91 structure; else the first ENST; else ''."""
    ensts = [e for e in ensts if e]
    for e in ensts:
        if e in e2s:
            return e
    return ensts[0] if ensts else ""


def resolve(data_dir=None):
    data_dir = data_dir or config.DATA_DIR
    path = os.path.join(data_dir, "clone_to_structure.tsv")
    cts = pd.read_csv(config.require(path, what="clone_to_structure"), sep="\t", dtype=str,
                      keep_default_na=False, na_values=[""]).fillna("")
    clist = _load_clonelist()
    mane_ensg, mane_sym = _load_mane()
    gencode = sequence_map.load_gencode_protein_index()      # protein -> {ensts,...}
    e2s = reference_structures.load_enst_to_structure()      # base ENST -> {structure, gene} (Ens91)

    src_counts = {}
    for i, r in cts.iterrows():
        cid, gene, ensg = r["clone_id"], r.get("gene_symbol", ""), r.get("ensg", "")
        cl = clist.get(cid, {})
        enst, source = "", ""
        # 1. CloneList authoritative ENST (prefer one with an Ens91 structure)
        if cl.get("ensts"):
            enst, source = _prefer_ens91(cl["ensts"], e2s), "clonelist"
        # 2. existing resolved ENST
        if not enst and str(r.get("enst", "")).startswith("ENST"):
            enst, source = r["enst"], (r.get("enst_source") or "existing")
        # 3. CloneList CDS -> protein -> GENCODE -> ENST (novel clones)
        if not enst and cl.get("cds"):
            prot = sequence_map.translate_orf(cl["cds"])
            meta = gencode.get(sequence_map._norm_protein(prot)) if prot else None
            if meta:
                enst, source = _prefer_ens91(meta.get("ensts", []), e2s), "clonelist_cds_seq"
        # 4. MANE Select for the gene (the standard "main select isoform" reference)
        if not enst:
            m = mane_ensg.get(str(ensg).split(".")[0]) or mane_sym.get(gene)
            if m:
                enst, source = m, "mane_select_ref"
        if enst:
            cts.at[i, "enst"] = enst
            cts.at[i, "enst_source"] = source
            # ensure an Ens91 structure: keep clone's own if present, else take the ENST's Ens91 structure
            if not str(r.get("structure", "")).strip():
                cts.at[i, "structure"] = e2s.get(enst, {}).get("structure", "")
            if not str(r.get("ensg", "")).startswith("ENSG"):
                cts.at[i, "ensg"] = e2s.get(enst, {}).get("gene", "") or ensg
            src_counts[source] = src_counts.get(source, 0) + 1

    cts.to_csv(path, sep="\t", index=False)
    have = int(cts["enst"].astype(str).str.startswith("ENST").sum())
    logger.info("[iso2function.clonelist] clone_to_structure: %d/%d clones have an ENST (sources this run: %s)",
                have, len(cts), src_counts)
    miss = cts[~cts["enst"].astype(str).str.startswith("ENST")]
    if len(miss):
        logger.warning("[iso2function.clonelist] %d clones STILL without ENST: %s", len(miss),
                       list(miss["clone_id"])[:10])
    return cts


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    resolve()
