"""Add the most accurate ENST to Vidal2025 atlas clones by PROTEIN-SEQUENCE matching.

The atlas clones were keyed to structures via the bridge; ENSTs came from cohort-known exact matches or
gene-restricted structure containment (parent transcript). This step upgrades the non-cohort-known
clones to a SEQUENCE-VERIFIED ENST: take each clone's actual protein (precomputed in the MDS-AML
`protein_sequences.fasta`, keyed by its final_isoform_id), match it to GENCODE proteins (exact, then
gene-restricted clone-variant >=99% gap-aware), and resolve to the Ensembl91 reference ENST (trying all
of the protein's ENSTs). Only ACCURATE tiers are applied (protein_exact / protein_clone) — no
approximate nearest match; clones with neither keep their structure-based ENST.

Updates ``clone_to_structure.tsv`` in place (``enst``, ``enst_source``). Heavy (loads the 500 MB protein
FASTA once); run on demand, not part of the default crosswalk.
"""

import os
import logging

import pandas as pd

from .. import config
from . import sequence_map, reference_structures

logger = logging.getLogger(__name__)

ACCURATE = {"cohort_known", "protein_exact", "protein_clone"}  # don't override these with weaker tiers


def load_protein_fasta(path=None):
    """{final_isoform_id (first header token) -> normalized protein} from protein_sequences.fasta."""
    path = config.require(path or config.MDS_AML_PROTEIN_FASTA, what="MDS-AML protein FASTA")
    idx, hdr, seq = {}, None, []

    def _flush(h, s):
        if not h:
            return
        key = h[1:].split()[0]
        p = sequence_map._norm_protein("".join(s))
        if key and p:
            idx.setdefault(key, p)

    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                _flush(hdr, seq); hdr, seq = line.strip(), []
            else:
                seq.append(line.strip())
        _flush(hdr, seq)
    logger.info("[iso2function.atlas_protein] protein FASTA: %d sequences", len(idx))
    return idx


def _protein_for(final_id, pfa):
    """Look up a clone's protein by its final_isoform_id, tolerating header-token variants."""
    f = str(final_id)
    for k in (f, f.split()[0] if f else "", f.split(".")[0]):
        if k and k in pfa:
            return pfa[k]
    return ""


def backfill_atlas_enst_by_protein(data_dir=None, min_identity=99.0):
    """Upgrade non-cohort-known atlas clones to sequence-verified ENSTs. Returns the updated DataFrame
    and rewrites ``clone_to_structure.tsv``."""
    data_dir = data_dir or config.DATA_DIR
    path = os.path.join(data_dir, "clone_to_structure.tsv")
    cts = pd.read_csv(config.require(path, what="clone_to_structure"), sep="\t", dtype=str,
                      keep_default_na=False, na_values=[""]).fillna("")
    if "enst_source" not in cts.columns:
        cts["enst_source"] = cts["enst"].map(lambda e: "cohort_known" if str(e).startswith("ENST") else "")

    pfa = load_protein_fasta()
    gencode = sequence_map.load_gencode_protein_index()       # protein -> {ensts, gene}
    e2s = reference_structures.load_enst_to_structure()       # base ENST -> {structure, gene}
    gene_prot = {}
    for prot, m in gencode.items():
        gene_prot.setdefault(m["gene"].upper(), []).append(prot)

    n_exact = n_clone = 0
    for i, r in cts.iterrows():
        if str(r.get("enst_source")) in ACCURATE:             # already sequence/cohort accurate
            continue
        prot = _protein_for(r.get("final_isoform_id", ""), pfa)
        if not prot:
            continue
        meta, tier = None, None
        if prot in gencode:
            meta, tier = gencode[prot], "protein_exact"
        else:                                                  # gene-restricted clone-variant >=99%
            best, bid = None, 0.0
            for cand in gene_prot.get(str(r.get("gene_symbol", "")).upper(), []):
                pid = sequence_map.pct_identity(prot, cand)
                if pid > bid:
                    best, bid = cand, pid
            if best and bid >= min_identity:
                meta, tier = gencode.get(best), "protein_clone"
        if not meta:
            continue
        enst, struct, ensg = reference_structures.resolve_structure(meta["ensts"], e2s)
        if enst:
            cts.at[i, "enst"] = enst
            cts.at[i, "enst_source"] = tier
            if not str(r.get("structure", "")).strip() and struct:
                cts.at[i, "structure"] = struct
            if not str(r.get("ensg", "")).startswith("ENSG") and ensg:
                cts.at[i, "ensg"] = ensg
            n_exact += int(tier == "protein_exact"); n_clone += int(tier == "protein_clone")

    cts.to_csv(path, sep="\t", index=False)
    src = cts["enst_source"].value_counts().to_dict()
    logger.info("[iso2function.atlas_protein] upgraded %d clones (protein_exact=%d, protein_clone=%d); "
                "enst_source now: %s", n_exact + n_clone, n_exact, n_clone, src)
    logger.info("[iso2function.atlas_protein] atlas clones with an ENST: %d / %d",
                int((cts["enst"].astype(str).str.startswith("ENST")).sum()), len(cts))
    return cts


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    backfill_atlas_enst_by_protein()
