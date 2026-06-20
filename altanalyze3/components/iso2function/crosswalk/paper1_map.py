"""Crosswalk paper1 (Yang 2016) isoforms to the canonical key by PROTEIN SEQUENCE.

paper1 isoforms are labeled ``SYMBOL_N`` with no Ensembl id, but ship full ORF nucleotide sequences.
We translate each ORF and match it to the GENCODE protein set (exact, then high-identity same-gene) to
assign a sequence-verified ENST/ENSP, then attach the structure key where we already have it
(``crosswalk_structure.tsv``, the TFIso bridge's known ENSTs). For paper1's non-TF / non-bridge genes
the ENST is the sequence-verified canonical key (a structure string can be derived later by running
gff_process on the reference). No ID guessing — every assignment is sequence-backed.

Output ``paper1_isoform_map.tsv``: paper1 Isoform_ID -> gene, ENST/ENSP, structure (if known),
percent identity, match mode (exact / fuzzy>=thr / unmatched).
"""

import os
import logging

import pandas as pd

from .. import config
from . import sequence_map, reference_structures

logger = logging.getLogger(__name__)


def build_paper1_crosswalk(data_dir=None, min_identity=99.0):
    """Crosswalk every paper1 ORF isoform to ENST by sequence. Returns the map DataFrame and writes
    ``paper1_isoform_map.tsv`` + appends a coverage row. ``min_identity`` (<100) enables a same-gene
    fuzzy fallback for isoforms whose protein drifted from GENCODE (build/annotation differences)."""
    data_dir = data_dir or config.DATA_DIR
    orf_path = os.path.join(data_dir, "paper1_orfs.tsv")
    if not os.path.exists(orf_path):
        raise FileNotFoundError("paper1_orfs.tsv not found — run `iso2function ingest` first.")
    orf = pd.read_csv(orf_path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    seq_col = next(c for c in orf.columns if "Reading_Frame" in c or "Sequence" in c)

    query_proteins, gene_of_query, meta = {}, {}, {}
    for _, r in orf.iterrows():
        iid = r["Isoform_ID"]
        prot = sequence_map.translate_orf(r.get(seq_col))
        if not prot:
            continue
        query_proteins[iid] = prot
        gene_of_query[iid] = str(r.get("Gene_Symbol", "")).upper()
        meta[iid] = {"gene_symbol": r.get("Gene_Symbol", ""), "category": r.get("Isoform_Category", ""),
                     "novelty": r.get("Isoform_Novelty", "")}

    gencode = sequence_map.load_gencode_protein_index()
    gene_of_target = {prot: m["gene"].upper() for prot, m in gencode.items()}
    res = sequence_map.crosswalk_sequences(query_proteins, gencode, min_identity=min_identity,
                                           gene_of_query=gene_of_query, gene_of_target=gene_of_target)

    # comprehensive ENST -> Ensembl91 structure (244k) + structure -> genomic splice coords (reused)
    e2s = reference_structures.load_enst_to_structure()
    coords = reference_structures.StructureCoords()
    rows = []
    for r in res:
        iid = r["query_id"]
        # try ALL ENSTs of the matched protein so the build-gap ENST that has a structure is found
        ensts = r.get("ensts") or ([r.get("enst")] if r.get("enst") else [])
        enst, struct, ensg = reference_structures.resolve_structure([str(e).split(".")[0] for e in ensts], e2s)
        junc = coords.junction_string(ensg, struct) if struct else ""
        rows.append({
            "isoform_id": iid, **meta.get(iid, {}),
            "matched": r["matched"], "identity": r["identity"],
            "match_mode": ("exact" if r["identity"] == 100.0 else
                           (f"clone>={min_identity:g}" if r["matched"] else "unmatched")),
            "enst": enst, "ensp": r.get("ensp", ""), "ensg": ensg, "matched_gene": r.get("gene", ""),
            "structure": struct, "has_structure": bool(struct),
            "genomic_junctions": junc or "",
        })
    coords.close()
    out = pd.DataFrame(rows)
    path = os.path.join(data_dir, "paper1_isoform_map.tsv")
    out.to_csv(path, sep="\t", index=False)
    n_exact = int((out["match_mode"] == "exact").sum())
    n_clone = int(out["match_mode"].str.startswith("clone").sum())
    n_struct = int(out["has_structure"].sum())
    n_coords = int((out["genomic_junctions"].astype(str) != "").sum())
    logger.info("[iso2function.crosswalk.paper1] %d isoforms: %d exact + %d clone-variant = %d "
                "ENST-mapped; %d -> Ensembl91 structure; %d with genomic splice coords; %d unmatched "
                "(novel alt isoforms with no reference match)",
                len(out), n_exact, n_clone, n_exact + n_clone, n_struct, n_coords,
                int((~out["matched"]).sum()))
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config.ensure_dirs()
    df = build_paper1_crosswalk()
    print(df["match_mode"].value_counts().to_string())
