"""Protein-sequence -> structure-key crosswalk (shared by paper1 and UniProt).

Some isoform sources (paper1 ORFs, UniProt isoforms) carry a protein/ORF SEQUENCE but no Ensembl
transcript id. This module maps such a sequence to our canonical structure key by protein-sequence
identity against a TARGET INDEX of isoforms whose structure IS known (e.g. translated Ensembl/GENCODE
transcripts, or the TFIso bridge clones). Sequence identity is the arbiter — never an ID guess.

The matcher is source-agnostic: callers supply
  - a query map   {query_id -> protein_sequence}
  - a target index {protein_sequence -> {"transcript_id":..., "structure":..., "enst":...}}
and get back a per-query assignment (exact match, or best high-identity match above a threshold, or
unresolved). The target index is built separately (``build_target_index_from_fasta`` /
``build_target_index_from_translation``) from whatever protein resource is available.
"""

import os
import logging
from difflib import SequenceMatcher

from Bio.Seq import Seq

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- translation
def translate_orf(nt_seq):
    """Translate an ORF NUCLEOTIDE sequence (ATG..stop) to protein, trimming a trailing stop. Returns
    '' for empty/invalid input. Length not a multiple of 3 is truncated to the last full codon (paper1
    ORFs are clean ORFs, but this is defensive)."""
    if not nt_seq:
        return ""
    s = str(nt_seq).strip().upper().replace("U", "T")
    s = "".join(c for c in s if c in "ACGTN")
    s = s[: len(s) - (len(s) % 3)]
    if not s:
        return ""
    return str(Seq(s).translate(to_stop=True))


def _norm_protein(p):
    return (p or "").strip().rstrip("*").upper()


# --------------------------------------------------------------------------- identity
def pct_identity(a, b):
    """Gap-aware percent similarity between two proteins (difflib ratio * 100). Unlike an ungapped
    N-terminal comparison, this tolerates the indels that distinguish isoforms / hORFeome clone variants
    (a 1-aa N-terminal shift no longer collapses identity to ~0). Returns 0 if either is empty. Used
    only to recognize near-identical (clone-variant) proteins; exact matches use a dict lookup."""
    a, b = _norm_protein(a), _norm_protein(b)
    if not a or not b:
        return 0.0
    return 100.0 * SequenceMatcher(None, a, b).ratio()


# --------------------------------------------------------------------------- target index builders
def build_target_index_from_fasta(fasta_path, id_to_meta=None):
    """Build {normalized_protein -> meta} from a protein FASTA. ``id_to_meta`` optionally maps a FASTA
    record id (e.g. ENST/ENSP) to a metadata dict (e.g. {'structure':..., 'enst':...}); records with no
    metadata still index with {'transcript_id': record_id}. Requires BioPython SeqIO."""
    from Bio import SeqIO
    index = {}
    id_to_meta = id_to_meta or {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        key = _norm_protein(str(rec.seq))
        if not key:
            continue
        meta = dict(id_to_meta.get(rec.id, {}))
        meta.setdefault("transcript_id", rec.id)
        index.setdefault(key, meta)
    logger.info("[iso2function.sequence_map] target index: %d distinct proteins from %s",
                len(index), os.path.basename(str(fasta_path)))
    return index


def build_target_index_from_records(records):
    """Build {normalized_protein -> meta} from an iterable of (protein_seq, meta_dict)."""
    index = {}
    for seq, meta in records:
        key = _norm_protein(seq)
        if key:
            index.setdefault(key, dict(meta))
    return index


def load_gencode_protein_index(fasta_path=None):
    """Build {normalized_protein -> {'enst','ensts','ensp','gene'}} from the GENCODE pc_translations
    FASTA (header: >ENSP|ENST|ENSG|OTT|OTT|name-2xx|SYMBOL|len). The SAME protein sequence is encoded by
    MANY transcripts, so ``ensts`` accumulates ALL base ENSTs for that protein (critical: lets the
    structure lookup find whichever ENST exists in the Ensembl91 reference table, instead of discarding
    it). ``enst`` is the first seen. Versions stripped to the stable base id."""
    import gzip
    from .. import config
    fasta_path = fasta_path or config.require(config.GENCODE_PROTEIN_FASTA, what="GENCODE protein FASTA")
    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    index, hdr, seq = {}, None, []

    def _flush(h, s):
        if not h:
            return
        prot = _norm_protein("".join(s))
        if not prot:
            return
        parts = h[1:].split("|")
        ensp = parts[0].split(".")[0] if parts else ""
        enst = next((p.split(".")[0] for p in parts if p.startswith("ENST")), "")
        gene = parts[6] if len(parts) > 6 else ""
        rec = index.get(prot)
        if rec is None:
            index[prot] = {"enst": enst, "ensts": [enst] if enst else [], "ensp": ensp, "gene": gene}
        elif enst and enst not in rec["ensts"]:
            rec["ensts"].append(enst)

    with opener(fasta_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                _flush(hdr, seq); hdr, seq = line.strip(), []
            else:
                seq.append(line.strip())
        _flush(hdr, seq)
    logger.info("[iso2function.sequence_map] GENCODE protein index: %d distinct proteins", len(index))
    return index


# --------------------------------------------------------------------------- matcher
def match_one(query_protein, target_index, min_identity=100.0, allowed_targets=None):
    """Match a single query protein to ``target_index``. Exact (normalized) match first; if none and
    ``min_identity`` < 100, scan ``allowed_targets`` (a subset of proteins to compare against, e.g. the
    same gene's proteins — REQUIRED for the fuzzy path to stay tractable) for the best identity >=
    ``min_identity``. Returns (meta, identity) or (None, 0.0)."""
    q = _norm_protein(query_protein)
    if not q:
        return None, 0.0
    if q in target_index:
        return target_index[q], 100.0
    if min_identity >= 100.0 or not allowed_targets:
        return None, 0.0
    best, best_id = None, 0.0
    for tp in allowed_targets:
        pid = pct_identity(q, tp)
        if pid > best_id:
            best, best_id = target_index.get(_norm_protein(tp)), pid
    return (best, best_id) if best_id >= min_identity else (None, 0.0)


def crosswalk_sequences(query_proteins, target_index, min_identity=100.0, gene_of_query=None,
                        gene_of_target=None):
    """Map a dict {query_id -> protein} to target metadata by sequence. With ``min_identity`` < 100 a
    fuzzy fallback runs per query restricted to same-gene targets (needs ``gene_of_query`` /
    ``gene_of_target`` maps to build the allowed set), so it never does an all-vs-all scan. Returns a
    list of dicts: {query_id, matched (bool), identity, ...target meta}."""
    same_gene = None
    if min_identity < 100.0 and gene_of_query and gene_of_target:
        same_gene = {}
        for prot, g in gene_of_target.items():
            same_gene.setdefault(g, []).append(prot)
    out = []
    for qid, prot in query_proteins.items():
        allowed = None
        if same_gene is not None and gene_of_query and qid in gene_of_query:
            allowed = same_gene.get(gene_of_query[qid])
        meta, pid = match_one(prot, target_index, min_identity=min_identity, allowed_targets=allowed)
        row = {"query_id": qid, "matched": meta is not None, "identity": round(pid, 2)}
        if meta:
            row.update(meta)
        out.append(row)
    n_ok = sum(1 for r in out if r["matched"])
    logger.info("[iso2function.sequence_map] crosswalk: %d/%d query proteins matched (min_identity=%s)",
                n_ok, len(out), min_identity)
    return out
