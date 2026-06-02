"""Obtain Ensembl/GENCODE reference transcript structures as AltAnalyze exon-structure strings.

EVIDENCE-BASED DESIGN NOTE (do not re-introduce a bespoke annotator here):
  Reference (ENST) structures MUST be produced by the SAME ``gff_process`` run that processes the
  sample read GFFs -- i.e. the GENCODE/Ensembl reference GFF is submitted as just another input
  GFF to ``gff_process.consolidateLongReadGFFs``. Annotating the reference GFF *in isolation*
  assigns sub-exon version tokens (E5.3 vs E5.4, U1.1_<coord> prefixes, ...) from a DIFFERENT
  exon-token namespace than the combined sample run, so an ENST structure derived alone will NOT
  exactly match the same ENST as it appears alongside the reads -- silently breaking ENST-priority.
  This was verified directly for ENST00000262407 (ITGA2B): the isolated annotation and the combined
  ``transcript_associations.txt`` produced DIFFERENT structure strings for the identical transcript.

  Therefore this module does NOT annotate. It consumes the ``transcript_associations.txt`` that
  ``gff_process`` writes (columns: gene, strand, structure, transcript_id, source), keeps the
  structure string verbatim (same token namespace as the reads), filters to reference transcripts
  (ENST*/NM_*/NR_*), and strips the transcript-id version suffix (".6") for display/quantification.

Output: in-memory dict ``{gene: {structure: ENST}}`` and a cached TSV
  ``gene <tab> structure <tab> ENST`` (version stripped).
"""

from __future__ import annotations

import os

from .. import gff_process as gff


#: transcript-id prefixes that denote a curated reference model (vs a novel read molecule).
REFERENCE_ID_PREFIXES = ("ENST", "NM_", "NR_", "XM_", "XR_")


def _is_reference_id(tid):
    return any(tid.startswith(p) for p in REFERENCE_ID_PREFIXES)


def _strip_version(tid):
    """ENST00000262407.6 -> ENST00000262407 (display id only; structure is kept verbatim)."""
    return tid.split('.')[0]


def _load_from_transcript_associations(ta_path):
    """Parse a gff_process transcript_associations.txt into {gene: {structure: ENST}}.

    Columns: gene, strand, structure, transcript_id, source. Only reference transcripts are kept.
    The structure is taken VERBATIM (same exon-token namespace as the reads in the same run); only
    the transcript-id version is stripped.
    """
    out = {}
    with open(ta_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            gene, _strand, struct, tid = parts[0], parts[1], parts[2], parts[3]
            if not struct or 'UNK' in gene or not _is_reference_id(tid):
                continue
            out.setdefault(gene, {})[struct] = _strip_version(tid)
    return out


def annotate_reference(ref_gff, ref_exons, cache_path=None, force=False, log=print,
                       transcript_associations=None):
    """Return {gene: {structure: ENST}} for the reference transcripts.

    The structures come from ``gff_process`` (NOT re-derived here), so they share the read
    isoforms' exon-token namespace and ENST-priority matching is exact.

    Resolution order:
      1. cache_path TSV (gene, structure, ENST) if present and not ``force``;
      2. an explicit ``transcript_associations`` path (a gff_process output);
      3. the gff-output/transcript_associations.txt next to ``ref_gff`` if already produced;
      4. otherwise run ``gff_process.consolidateLongReadGFFs([ref_gff], ref_exons, mode='Ensembl')``
         -- submitting the reference GFF as just another input GFF -- and read its output.
    """
    # 1. our own cached {gene, structure, ENST}
    if cache_path and not force and os.path.exists(cache_path):
        out = {}
        with open(cache_path) as f:
            next(f, None)
            for line in f:
                g, struct, enst = line.rstrip("\n").split("\t")
                out.setdefault(g, {})[struct] = enst
        log(f"[reference] loaded cached ENST structures: {sum(len(v) for v in out.values()):,} "
            f"transcripts, {len(out):,} genes")
        return out

    # 2/3. an existing gff_process transcript_associations.txt
    ta_path = transcript_associations
    if not ta_path and ref_gff:
        candidate = os.path.join(os.path.dirname(ref_gff), 'gff-output', 'transcript_associations.txt')
        if os.path.exists(candidate):
            ta_path = candidate

    # 4. produce it via the SAME gff_process workflow the samples use (reference GFF as an input GFF)
    if (not ta_path or not os.path.exists(ta_path)) and ref_gff:
        log(f"[reference] no gff_process output found; running consolidateLongReadGFFs on {ref_gff}")
        ta_path = gff.consolidateLongReadGFFs([ref_gff], ref_exons, mode='Ensembl')

    if not ta_path or not os.path.exists(ta_path):
        raise FileNotFoundError(
            "Could not obtain reference structures: no cache, no transcript_associations.txt, "
            "and no ref_gff to run gff_process on.")

    out = _load_from_transcript_associations(ta_path)

    if cache_path:
        os.makedirs(os.path.dirname(cache_path) or '.', exist_ok=True)
        with open(cache_path, 'w') as o:
            o.write("gene\tstructure\tENST\n")
            for g, m in out.items():
                for struct, enst in m.items():
                    o.write(f"{g}\t{struct}\t{enst}\n")
        log(f"[reference] loaded {sum(len(v) for v in out.values()):,} reference structures "
            f"({len(out):,} genes) from {os.path.basename(ta_path)} -> cached {cache_path}")
    return out
