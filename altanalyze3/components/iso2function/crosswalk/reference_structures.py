"""Map any ENST to its Ensembl91 structure string, and any structure to its genomic splice-site
coordinates — the two portable representations every isoform annotation carries.

Reuses the MDS-AML collapse outputs (no recompute):
  - ``ENST_reference_structures.tsv`` (gene, structure, ENST; 244k transcripts) -> ENST -> structure
  - ``structure_coords.sqlite`` (gene, structure, chrom, strand, blob) -> structure -> genomic exons,
    decoded with ``bam.coord_store.unpack_exons`` -> splice-site (donor/acceptor) coordinates.

Rationale (user requirement): the Ensembl91 structure string is the canonical key for the CURRENT
annotation build; the genomic splice-site coordinates are an alternative, build-independent annotation
so a user who rebuilds an isoform reference against a DIFFERENT database build can re-derive structures
from the raw junctions. Both are attached to every crosswalked isoform.
"""

import os
import re
import sqlite3
import logging

import pandas as pd

from .. import config
from ...bam import coord_store

logger = logging.getLogger(__name__)


def _norm_struct(s):
    """Strip novel-site coordinate suffixes (E8.3_52052622 -> E8.3) so a novel isoform's structure can
    be compared against the clean reference exon vocabulary."""
    return re.sub(r"_\d+", "", str(s or ""))


def load_gene_to_structures(path=None):
    """{ENSG -> [(normalized_structure, base_ENST)]} from ENST_reference_structures — for resolving an
    isoform's best ENST by structure containment within its own gene."""
    path = config.require(path or config.ENST_REFERENCE_STRUCTURES, what="ENST reference structures")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""]).fillna("")
    idx = {}
    for _, r in df.iterrows():
        enst = str(r.get("ENST", "")).split(".")[0]
        if enst.startswith("ENST") and r.get("structure"):
            idx.setdefault(r.get("gene", ""), []).append((_norm_struct(r["structure"]), enst))
    logger.info("[iso2function.reference_structures] gene->structures index: %d genes", len(idx))
    return idx


def best_enst_by_structure(structure, ensg, gene_index):
    """Resolve an isoform's best reference ENST by structure: returns (enst, relation) where the
    isoform structure (within gene ``ensg``) is 'exact' to, or 'contained_in', a reference transcript;
    else ('', ''). 'contained_in' yields the parent transcript the isoform folds into — a meaningful
    best-ENST for a novel sub-isoform. Requires structure_key (imported lazily to avoid a cycle)."""
    from . import structure_key
    cs = _norm_struct(structure)
    if not cs or not ensg:
        return "", ""
    parent = ""
    for ref_struct, enst in gene_index.get(ensg, ()):
        mt = structure_key.match_structures(cs, ref_struct)
        if mt == "exact":
            return enst, "ref_exact"
        if mt == "contained_in" and not parent:
            parent = enst
    return (parent, "ref_contained_in") if parent else ("", "")


def load_enst_to_structure(path=None):
    """Return {base_ENST -> {'structure':..., 'gene':ENSG}} from ENST_reference_structures.tsv
    (version-stripped ENST keys)."""
    path = config.require(path or config.ENST_REFERENCE_STRUCTURES, what="ENST reference structures")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    out = {}
    for _, r in df.iterrows():
        enst = str(r.get("ENST", "")).split(".")[0]
        if enst.startswith("ENST") and r.get("structure"):
            out.setdefault(enst, {"structure": r["structure"], "gene": r.get("gene", "")})
    logger.info("[iso2function.reference_structures] ENST->structure: %d transcripts", len(out))
    return out


def resolve_structure(ensts, e2s):
    """Given candidate base ENSTs (all transcripts encoding one protein) and the ENST->{structure,gene}
    map, return (chosen_enst, structure, ensg) for the FIRST ENST that has an Ensembl91 structure; if
    none do, return (first_enst, '', ''). This is the fix for the build gap: one of a protein's many
    ENSTs is usually present in the reference even when the first-seen (newest-build) one is not."""
    ensts = [e for e in (ensts or []) if e]
    for e in ensts:
        if e in e2s:
            return e, e2s[e]["structure"], e2s[e].get("gene", "")
    return (ensts[0] if ensts else ""), "", ""


def load_exon_reference(path=None):
    """{(ENSG, exon_token) -> (chrom, strand, start, stop)} from the Ensembl91 exon reference
    (Hs_Ensembl_exon.txt). Maps any Ens91 E/I structure token to its genomic coordinates — the
    universal source of per-isoform exon coordinates (works for any structure, not just cohort-observed
    ones)."""
    path = config.require(path or config.ENSEMBL91_EXON, what="Ensembl91 exon reference")
    ref = {}
    with open(path) as fh:
        next(fh, None)  # header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 6:
                try:
                    ref[(f[0], f[1])] = (f[2], f[3], int(f[4]), int(f[5]))
                except ValueError:
                    continue
    logger.info("[iso2function.reference_structures] Ens91 exon reference: %d (gene,exon) entries", len(ref))
    return ref


def exon_coords_for_structure(ensg, structure, exon_ref):
    """Map an Ens91 structure to 'chrom:strand:(s1,e1),(s2,e2),...' genomic exon tuples via the exon
    reference. Novel-site coordinate suffixes (E8.3_52052622) are stripped to the base token for the
    region lookup. Tokens absent from the reference are skipped; '' if none resolve. Blocks are sorted
    by genomic start."""
    if not (ensg and structure):
        return ""
    chrom = strand = None
    blocks = []
    for tok in str(structure).split("|"):
        base = _norm_struct(tok)
        hit = exon_ref.get((ensg, base))
        if hit:
            chrom, strand = hit[0], hit[1]
            blocks.append((hit[2], hit[3]))
    if not blocks:
        return ""
    blocks = sorted(set(blocks))
    return f"{chrom}:{strand}:" + ",".join(f"({s},{e})" for s, e in blocks)


class StructureCoords:
    """Lazy reader over ``structure_coords.sqlite`` mapping (gene, structure) -> genomic coordinates."""

    def __init__(self, db_path=None):
        self.db_path = config.require(db_path or config.STRUCTURE_COORDS_SQLITE,
                                      what="structure coords sqlite")
        self._conn = sqlite3.connect(f"file:{self.db_path}?mode=ro", uri=True)

    def exons(self, gene, structure):
        """Return (chrom, strand, [(start,end),...]) for a (gene, structure), or None if absent.
        Exons are genomic, sorted by start."""
        row = self._conn.execute(
            "SELECT chrom, strand, blob FROM structures WHERE gene=? AND structure=?",
            (gene, structure)).fetchone()
        if not row:
            return None
        chrom, strand, blob = row
        exons = sorted(coord_store.unpack_exons(blob), key=lambda x: x[0])
        return chrom, strand, exons

    def splice_sites(self, gene, structure):
        """Return the list of intron (donor, acceptor) genomic coordinate pairs between consecutive
        exons (genomic order). Empty for single-exon. Donor = 5' intron base, acceptor = 3' intron base,
        oriented by strand. Returns None if the structure is not in the store."""
        got = self.exons(gene, structure)
        if got is None:
            return None
        chrom, strand, exons = got
        if len(exons) < 2:
            return {"chrom": chrom, "strand": strand, "junctions": []}
        junctions = []
        for (s1, e1), (s2, e2) in zip(exons[:-1], exons[1:]):
            # intron spans e1..s2 (1-based inclusive exon ends); orient donor/acceptor by strand
            donor, acceptor = (e1, s2) if strand == "+" else (s2, e1)
            junctions.append((donor, acceptor))
        return {"chrom": chrom, "strand": strand, "junctions": junctions}

    def junction_string(self, gene, structure):
        """Compact, portable splice-site string: 'chrom:strand:donor1-acceptor1,donor2-acceptor2,...'.
        '' if single-exon; None if structure absent from the store."""
        ss = self.splice_sites(gene, structure)
        if ss is None:
            return None
        if not ss["junctions"]:
            return f"{ss['chrom']}:{ss['strand']}:"
        joined = ",".join(f"{d}-{a}" for d, a in ss["junctions"])
        return f"{ss['chrom']}:{ss['strand']}:{joined}"

    def close(self):
        try:
            self._conn.close()
        except Exception:
            pass
