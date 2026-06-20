"""Resolve structureless isoforms against the MDS-AML full-length catalog, and build a gene-indexed
isoform sequence/structure store for fast viewer retrieval.

Part 1 (resolve): every isoform in isoform_interactions.txt that still lacks an Ens91 structure has its
protein matched (exact) against the cohort's full-length isoform proteins (protein_sequences.fasta,
known+novel) -> final_isoform_id -> Ens91 structure (transcript_associations.txt). This recovers
structures for isoforms the Ensembl-based mapping missed but that were observed full-length in the cohort.

Part 2 (index): write ``isoform_seq_index.sqlite`` -- table isoforms(gene, source, source_isoform_id,
final_isoform_id, best_ENST, structure, protein) with an index on gene and on source_isoform_id, so the
viewer retrieves a gene's isoforms (sequence + structure) with `WHERE gene=?` instantly.

Protein per isoform: atlas (Vidal2025) via final_isoform_id -> protein_sequences.fasta; Yang (Vidal2016)
by translating its ORF (paper1_orfs); UniProt from the Daedalus uniprot_isoform_proteins table.
"""

import os
import sqlite3
import logging

import pandas as pd

from . import config
from .crosswalk import sequence_map

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- protein sources
def _yang_proteins(data_dir):
    """{Yang Isoform_ID -> protein} by translating paper1 ORFs."""
    p = os.path.join(data_dir, "paper1_orfs.tsv")
    if not os.path.exists(p):
        return {}
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, na_values=[""]).fillna("")
    col = next((c for c in df.columns if "Reading_Frame" in c or "Sequence" in c), None)
    return {r["Isoform_ID"]: sequence_map.translate_orf(r.get(col, "")) for _, r in df.iterrows()} if col else {}


def _uniprot_proteins():
    """{UniProt id -> protein} keyed by both isoform_id and bare accession (canonical)."""
    p = config.resolve_shared(config.UNIPROT_ISOFORM_PROTEINS)
    if not os.path.exists(p):
        return {}
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, na_values=[""]).fillna("")
    df = df[df["species"] == "human"] if "species" in df.columns else df
    out = {}
    for _, r in df.iterrows():
        seq = r.get("protein_sequence", "")
        if seq:
            out[r.get("isoform_id", "")] = seq
            out.setdefault(r.get("primary_accession", ""), seq)
    return out


def _atlas_final_ids(data_dir):
    """{atlas clone_id -> final_isoform_id} from clone_to_structure."""
    p = os.path.join(data_dir, "clone_to_structure.tsv")
    if not os.path.exists(p):
        return {}
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, na_values=[""]).fillna("")
    return {r["clone_id"]: r.get("final_isoform_id", "") for _, r in df.iterrows()}


# --------------------------------------------------------------------------- catalog loaders (heavy)
def load_cohort_proteins(needed_finals=None, path=None):
    """Single pass over the cohort full-length protein FASTA (~2.1M). Returns
    (protein->final exact map, final->raw-protein for the ``needed_finals`` set only) — the second
    keeps memory bounded (we only retain the proteins of the isoforms we will actually index)."""
    path = config.require(path or config.MDS_AML_PROTEIN_FASTA, what="cohort protein FASTA")
    needed = set(needed_finals or ())
    p2f, f2p, hdr, seq = {}, {}, None, []

    def _flush(h, s):
        if not h:
            return
        fid = h[1:].split()[0]
        raw = "".join(s)
        norm = sequence_map._norm_protein(raw)
        if norm:
            p2f.setdefault(norm, fid)
            if fid in needed:
                f2p.setdefault(fid, raw.strip().rstrip("*"))

    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                _flush(hdr, seq); hdr, seq = line.strip(), []
            else:
                seq.append(line.strip())
        _flush(hdr, seq)
    logger.info("[iso2function.seq_index] cohort proteins: %d protein->final, %d final->protein (needed)",
                len(p2f), len(f2p))
    return p2f, f2p


def load_final_to_structure(path=None):
    """{final_isoform_id -> (gene_ENSG, structure)} from transcript_associations.txt (~2.14M)."""
    path = config.require(path or config.TRANSCRIPT_ASSOCIATIONS, what="transcript_associations")
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                gene, _strand, structure, fid = parts[0], parts[1], parts[2], parts[3]
                out.setdefault(fid, (gene, structure))
    logger.info("[iso2function.seq_index] final->structure index: %d", len(out))
    return out


# --------------------------------------------------------------------------- main
def build(data_dir=None):
    data_dir = data_dir or config.DATA_DIR
    inter = pd.read_csv(config.require(os.path.join(data_dir, "isoform_interactions.txt"),
                                       what="isoform_interactions.txt"),
                        sep="\t", dtype=str, keep_default_na=False, na_values=[""]).fillna("")
    iso = inter.drop_duplicates("source_isoform_id").copy()

    # protein per interaction isoform
    yang = _yang_proteins(data_dir)
    up = _uniprot_proteins()
    atlas_fid = _atlas_final_ids(data_dir)
    prot2final, final2prot = load_cohort_proteins(needed_finals=set(atlas_fid.values()))
    final2struct = load_final_to_structure()

    def protein_of(row):
        s, sid = row["source"], row["source_isoform_id"]
        if s.startswith("Vidal2016"):
            return yang.get(sid, "")
        if s == "UniProt":
            return up.get(sid, "") or up.get(str(sid).split("-")[0], "")
        if s.startswith("Vidal2025"):
            return final2prot.get(atlas_fid.get(sid, ""), "")
        return ""

    rows, n_resolved = [], 0
    for _, r in iso.iterrows():
        prot = protein_of(r)
        structure, fid, best_enst = r.get("isoform_structure", ""), "", r.get("best_ENST", "")
        # Part 1: resolve structureless by cohort protein match
        if not structure and prot:
            f = prot2final.get(sequence_map._norm_protein(prot))
            if f:
                gene, structure = final2struct.get(f, ("", ""))
                fid = f
                if structure:
                    n_resolved += 1
        # carry the cohort protein for the index (use the matched full-length protein if resolved)
        rows.append({"gene": r.get("Symbol", ""), "ensg": r.get("ENSG", ""),
                     "source": r.get("source", ""), "source_isoform_id": r.get("source_isoform_id", ""),
                     "final_isoform_id": fid, "best_ENST": best_enst, "structure": structure,
                     "protein": prot})
    out = pd.DataFrame(rows)
    res_path = os.path.join(data_dir, "structureless_resolved.tsv")
    out[out["final_isoform_id"] != ""].to_csv(res_path, sep="\t", index=False)
    logger.info("[iso2function.seq_index] resolved %d structureless isoforms to a cohort full-length "
                "structure -> %s", n_resolved, res_path)

    # Part 2: gene-indexed SQLite
    db = os.path.join(data_dir, "isoform_seq_index.sqlite")
    if os.path.exists(db):
        os.remove(db)
    conn = sqlite3.connect(db)
    conn.execute("""CREATE TABLE isoforms (gene TEXT, ensg TEXT, source TEXT, source_isoform_id TEXT,
                    final_isoform_id TEXT, best_ENST TEXT, structure TEXT, protein TEXT)""")
    conn.executemany("INSERT INTO isoforms VALUES (?,?,?,?,?,?,?,?)",
                     [(r["gene"], r["ensg"], r["source"], r["source_isoform_id"], r["final_isoform_id"],
                       r["best_ENST"], r["structure"], r["protein"]) for r in rows])
    conn.execute("CREATE INDEX idx_gene ON isoforms(gene)")
    conn.execute("CREATE INDEX idx_sid ON isoforms(source_isoform_id)")
    conn.commit()
    ngenes = conn.execute("SELECT COUNT(DISTINCT gene) FROM isoforms").fetchone()[0]
    nseq = conn.execute("SELECT COUNT(*) FROM isoforms WHERE protein!=''").fetchone()[0]
    conn.close()
    logger.info("[iso2function.seq_index] %s: %d isoforms across %d genes, %d with protein (indexed by gene)",
                db, len(rows), ngenes, nseq)
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config.ensure_dirs()
    build()
