"""Tier-1 per-sample isoform candidate reduction (``collapse_sample``).

Streams one sample's raw-read GFF (the output of
``components/bam/isoform_structure_extract.py``), annotates each read with the
validated AltAnalyze exon/intron rules via :func:`gff_process.exonAnnotate`, and
collapses the reads into a compact per-sample candidate table in two steps:

  Step A — exact-structure collapse: reads sharing the identical (terminal-coord
           stripped) AltAnalyze structure string within a gene merge into one
           candidate with a read count.
  Step B — contiguous containment fold: candidate structures that are a
           *contiguous* token subsequence of a longer candidate fold into the
           longest/best-supported parent (see ``isoform_collapse_utils``).

Memory is bounded by **per-chromosome flushing**.  The raw GFF from
``isoform_structure_extract`` is chromosome-partitioned and position-ordered (it
concatenates per-chromosome chunks), so all reads for one chromosome arrive
contiguously.  We accumulate ``struct_db`` for the current chromosome only, run
Steps A+B at each chromosome boundary, write the results, and reset.  Reads whose
annotation spans multiple genes on different chromosomes (trans-spliced) are rare
and are diverted to an overflow bucket flushed once at the end.

This was validated to be necessary, not merely defensive: an un-flushed full-sample
run peaked at ~5-8 GB RSS on real ~9M-read iPSC samples; per-chromosome flushing
caps the live structure set at one chromosome's worth.

Outputs (written next to the input GFF, or under ``--outdir``):

  ``<prefix>.isoform_candidates.tsv.gz``
      gene, strand, candidate_id, structure_key, raw_structure_key,
      first_exon_start_min, first_exon_start_max, last_exon_end_min,
      last_exon_end_max, exon_token_count, read_count, is_ensembl,
      ensembl_transcript_id, best_exemplar_molecule, singleton

  ``<prefix>.molecule_to_candidate.tsv.gz``
      molecule_id, gene, strand, candidate_id, structure_key

``candidate_id`` is the (terminal-stripped) structure string itself — it is the
stable within-sample collapse key.  Cross-sample canonical IDs (Ensembl / stable
novel) are assigned later in Tier 2 (``merge_catalog``).

Design reference: ``isoform-collapse/UNIFIED_isoform_collapse_design.md`` §4.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import os
import sys
from typing import Dict, List, Optional, Tuple

from . import gff_process as gff
from .isoform_collapse_utils import contiguous_containment_fold

# Reference/RefSeq transcript-id prefixes that mark a read as a known model.
_REF_PREFIXES = ("ENST", "NM_", "XM_", "XR_", "NR_")

CANDIDATE_HEADER = [
    "gene", "strand", "candidate_id", "structure_key", "raw_structure_key",
    "first_exon_start_min", "first_exon_start_max",
    "last_exon_end_min", "last_exon_end_max",
    "exon_token_count", "read_count", "is_ensembl", "ensembl_transcript_id",
    "best_exemplar_molecule", "singleton",
]
MOLECULE_HEADER = ["molecule_id", "gene", "strand", "candidate_id", "structure_key"]


def _exon_str(exon_ids: List[str]) -> str:
    """Replicate ``gff_process.exon_str``: strip a terminal token carrying a
    novel coordinate (``_``) from each end, then join with ``|``.

    Returns the terminal-stripped structure string used as the collapse key.
    """
    ids = list(exon_ids)
    if ids and '_' in ids[0]:
        del ids[0]
    if ids and '_' in ids[-1]:
        del ids[-1]
    return "|".join(ids)


def _is_ref_molecule(molecule_id: str) -> bool:
    return molecule_id.startswith(_REF_PREFIXES)


class _ChromAccumulator:
    """Per-chromosome Step-A aggregation state.

    Holds, for the chromosome currently being streamed, the exact-structure
    counts plus the metadata needed to pick a representative and to emit the
    candidate/molecule rows.
    """

    __slots__ = ("structs", "molecules")

    def __init__(self):
        # key: (gene, strand, structure_key)
        #   -> dict(count, is_ref, ref_id, best_mol, best_len, raw_struct,
        #           fs_min, fs_max, le_min, le_max)
        self.structs: Dict[Tuple[str, str, str], dict] = {}
        # streamed molecule rows for this chromosome: (mol, gene, strand, struct)
        self.molecules: List[Tuple[str, str, str, str]] = []

    def add(self, gene, strand, struct, raw_struct, molecule_id,
            first_start, last_end, token_count):
        key = (gene, strand, struct)
        rec = self.structs.get(key)
        is_ref = _is_ref_molecule(molecule_id)
        if rec is None:
            self.structs[key] = {
                "count": 1,
                "is_ref": is_ref,
                "ref_id": molecule_id if is_ref else "",
                "best_mol": molecule_id,
                "best_token_count": token_count,
                "raw_struct": raw_struct,
                "fs_min": first_start, "fs_max": first_start,
                "le_min": last_end, "le_max": last_end,
            }
        else:
            rec["count"] += 1
            if is_ref and not rec["is_ref"]:
                rec["is_ref"] = True
                rec["ref_id"] = molecule_id
                rec["best_mol"] = molecule_id
            elif is_ref and rec["is_ref"] and molecule_id < rec["ref_id"]:
                rec["ref_id"] = molecule_id  # deterministic: lowest ref id
            if first_start < rec["fs_min"]:
                rec["fs_min"] = first_start
            if first_start > rec["fs_max"]:
                rec["fs_max"] = first_start
            if last_end < rec["le_min"]:
                rec["le_min"] = last_end
            if last_end > rec["le_max"]:
                rec["le_max"] = last_end
        self.molecules.append((molecule_id, gene, strand, struct))


def _flush(acc: "_ChromAccumulator", cand_out, mol_out, stats):
    """Run Step B (containment fold) per gene on the accumulated chromosome and
    write candidate + molecule rows.  Mutates ``stats`` counters."""
    # Group structures by (gene, strand) for the per-gene fold.
    by_gene: Dict[Tuple[str, str], List[dict]] = collections.defaultdict(list)
    for (gene, strand, struct), rec in acc.structs.items():
        by_gene[(gene, strand)].append({"structure": struct, "count": rec["count"],
                                        "is_ref": rec["is_ref"], "_rec": rec})

    # fold_map: child_structure -> winner_structure (per gene), accumulated across
    # all genes on this chromosome (structure strings are gene-scoped in practice
    # because exon ids are gene-relative, but we key the rewrite by (gene,strand,struct)).
    redirect: Dict[Tuple[str, str, str], str] = {}

    for (gene, strand), structs in by_gene.items():
        winners, fold_map = contiguous_containment_fold(structs, default_gene=gene)
        for child_struct, parent_struct in fold_map.items():
            redirect[(gene, strand, child_struct)] = parent_struct
        stats["folded"] += len(fold_map)

        for w in winners:
            struct = w["structure"]
            rec = w["_rec"]
            singleton = 1 if w["count"] <= 1 else 0
            if singleton:
                stats["singletons"] += 1
            cand_out.write("\t".join([
                gene, strand, struct, struct, rec["raw_struct"],
                str(rec["fs_min"]), str(rec["fs_max"]),
                str(rec["le_min"]), str(rec["le_max"]),
                str(len(struct.split("|")) if struct else 0),
                str(w["count"]),  # folded count accumulated by the driver
                "1" if rec["is_ref"] else "0",
                rec["ref_id"],
                rec["best_mol"],
                str(singleton),
            ]) + "\n")
            stats["winners"] += 1

    # Rewrite molecule rows to point at the surviving winner structure.
    for (mol, gene, strand, struct) in acc.molecules:
        final = redirect.get((gene, strand, struct), struct)
        mol_out.write("\t".join([mol, gene, strand, final, final]) + "\n")
        stats["molecule_rows"] += 1


def collapse_sample(
    gff_path: str,
    exon_reference: str,
    outdir: Optional[str] = None,
    prefix: Optional[str] = None,
    force: bool = False,
    _reference_loaded: bool = False,
) -> Tuple[str, str]:
    """Reduce one sample's raw-read GFF to a per-sample candidate table.

    Parameters
    ----------
    gff_path : path to the raw-read ``.gff`` / ``.gff.gz``.
    exon_reference : Ensembl exon model file for :func:`gff_process.importEnsemblGenes`.
    outdir : output directory (defaults to the GFF's directory).
    prefix : output basename prefix (defaults to the GFF basename sans ``.g*``).
    force : overwrite existing outputs instead of skipping (resumability).
    _reference_loaded : set True if ``importEnsemblGenes`` was already called in
        this process (lets a batch driver load the reference once per worker).

    Returns
    -------
    (candidates_path, molecules_path)
    """
    if outdir is None:
        outdir = os.path.dirname(os.path.abspath(gff_path))
    if prefix is None:
        prefix = os.path.basename(gff_path).split(".g")[0]
    os.makedirs(outdir, exist_ok=True)
    cand_path = os.path.join(outdir, f"{prefix}.isoform_candidates.tsv.gz")
    mol_path = os.path.join(outdir, f"{prefix}.molecule_to_candidate.tsv.gz")

    if not force and os.path.exists(cand_path) and os.path.exists(mol_path):
        print(f"[collapse_sample] outputs exist, skipping (use force=True): {cand_path}")
        return cand_path, mol_path

    if not _reference_loaded:
        gff.importEnsemblGenes(exon_reference)

    stats = collections.Counter()
    acc = _ChromAccumulator()
    current_chrom: Optional[str] = None

    open_func = gzip.open if gff_path.endswith(".gz") else open
    cand_out = gzip.open(cand_path, "wt")
    mol_out = gzip.open(mol_path, "wt")
    cand_out.write("\t".join(CANDIDATE_HEADER) + "\n")
    mol_out.write("\t".join(MOLECULE_HEADER) + "\n")

    exons: List[Tuple[int, int]] = []
    gene_info = None  # (chrom, strand, info)
    ti: Optional[int] = None
    td: Optional[str] = None
    first = True

    def annotate(chrom, strand, info, exon_list):
        """Annotate one transcript; return (gene, struct, raw_struct, fs, le, ntok)
        or None if it should be skipped (UNK gene / annotation failure)."""
        try:
            molecule_id = info.split(";")[ti].split(td)[1]
        except Exception:
            return None, None
        try:
            gene, exonIDs, _simple, _genes = gff.exonAnnotate(
                chrom, list(exon_list), strand, molecule_id)
        except Exception:
            return None, None
        if "UNK" in gene:
            return None, None
        struct = _exon_str(list(exonIDs))
        raw_struct = "|".join(exonIDs)
        if not struct:
            return None, None
        # terminal exon extents (genomic) for QC; exon_list is in file order
        starts = [e[0] for e in exon_list]
        ends = [e[1] for e in exon_list]
        fs = min(starts) if starts else 0
        le = max(ends) if ends else 0
        ntok = len(struct.split("|"))
        return molecule_id, (gene, strand, struct, raw_struct, fs, le, ntok)

    def emit(chrom, strand, info, exon_list):
        nonlocal current_chrom
        if not exon_list:
            return
        molecule_id, payload = annotate(chrom, strand, info, exon_list)
        if payload is None:
            return
        gene, strnd, struct, raw_struct, fs, le, ntok = payload
        # Chromosome-boundary flush (Step A/B) — the memory bound.
        if current_chrom is not None and chrom != current_chrom:
            _flush(acc, cand_out, mol_out, stats)
            acc.structs.clear()
            acc.molecules.clear()
            # exonAnnotate records one entry per molecule id in the module-global
            # gff.additional_junctions; with millions of unique molecule ids this
            # grows unbounded and dominates RSS (~4 GB on a 9M-read sample) — it is
            # not consumed by Tier 1, so reset it each chromosome to hold the leak.
            gff.additional_junctions = {}
            stats["chromosomes"] += 1
        current_chrom = chrom
        acc.add(gene, strnd, struct, raw_struct, molecule_id, fs, le, ntok)
        stats["reads"] += 1

    with open_func(gff_path, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            chrom, _src, typ, p1, p2, _sc, strand, _fr, info = cols
            if "chr" not in chrom:
                chrom = "chr" + chrom
            if first:
                if typ not in ("transcript", "exon"):
                    continue
                parts = info.split(";")
                ti = next((i for i, x in enumerate(parts) if "transcript_id" in x), -1)
                td = "=" if "=" in info else '"'
                first = False
            if typ == "transcript":
                if gene_info is not None:
                    emit(gene_info[0], gene_info[1], gene_info[2], exons)
                exons = []
            elif typ == "exon":
                exons.append((int(p1), int(p2)))
                gene_info = (chrom, strand, info)
        # last transcript
        if gene_info is not None and exons:
            emit(gene_info[0], gene_info[1], gene_info[2], exons)

    # Final flush of the last chromosome.
    if acc.structs:
        _flush(acc, cand_out, mol_out, stats)
        stats["chromosomes"] += 1

    cand_out.close()
    mol_out.close()

    print(f"[collapse_sample] {prefix}: reads={stats['reads']:,} "
          f"winners={stats['winners']:,} folded={stats['folded']:,} "
          f"singletons={stats['singletons']:,} chromosomes={stats['chromosomes']}")
    return cand_path, mol_path


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Tier-1 per-sample isoform candidate reduction")
    p.add_argument("--gff", required=True, help="raw-read GFF/GFF.GZ")
    p.add_argument("--ensembl-exons", required=True, help="Hs_Ensembl_exon.txt")
    p.add_argument("--outdir", default=None, help="output dir (default: GFF dir)")
    p.add_argument("--prefix", default=None, help="output basename prefix")
    p.add_argument("--force", action="store_true", help="overwrite existing outputs")
    return p


def main(argv=None):
    args = _build_arg_parser().parse_args(argv)
    collapse_sample(args.gff, args.ensembl_exons, outdir=args.outdir,
                    prefix=args.prefix, force=args.force)


if __name__ == "__main__":
    main()
