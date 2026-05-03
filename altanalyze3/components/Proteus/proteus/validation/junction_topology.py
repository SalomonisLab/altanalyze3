"""
junction_topology.py: Validate that isoform-distinguishing junction peptides
fall within transmembrane or surface-accessible regions.

Conceptual workflow
-------------------
AltAnalyze3 (long-read) produces novel transcript GFF files. For each
reference–alternative isoform pair:

  1. The **junction peptide** is the amino acid sequence spanning the novel
     exon-exon junction that is *unique to the alternative isoform* — i.e.,
     absent from the reference protein and all other GENCODE isoforms.

  2. We map that junction peptide onto the full-length alternative protein
     sequence to find its start/end residue positions.

  3. We predict the transmembrane topology of the full-length protein using
     SequenceTopologyEncoder (Kyte-Doolittle, same as Daedalus).

  4. We test whether the junction peptide overlaps a TM helix, a predicted
     N-terminal signal region, or an extracellular loop (gaps between TM helices
     or beyond the last TM helix on the extracellular side).

  5. The result is a per-pair **topology_overlap** flag and a **surface_evidence**
     string describing how the unique peptide relates to the membrane.

This is the correct framing for surfaceome validation:
  - A junction peptide *within* a TM helix → novel topology / topology change
  - A junction peptide *in an extracellular loop* → potentially surface-exposed
    and accessible to cell surface proteomics / antibodies
  - A junction peptide *in a cytoplasmic loop* → not surface-accessible
  - A junction peptide causing NMD → protein absent from surface

Usage
-----
from proteus.validation.junction_topology import JunctionTopologyValidator

validator = JunctionTopologyValidator()

results = validator.validate(
    pairs=[
        {
            "pair_id": "ENST00000001.1_vs_ALT1",
            "reference_protein_id": "ENSP00000001.1",
            "alternative_protein_seq": "MKTII...",
            "junction_peptide": "GGSLKV",   # the unique spanning peptide
            "has_alt_protein": True,         # False if NMD
        },
        ...
    ]
)
validator.to_tsv(results, "junction_topology_validation.tsv")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


# Topology category codes
TOPO_UNKNOWN         = "unknown"
TOPO_SOLUBLE         = "soluble"
TOPO_TM_HELIX        = "transmembrane_helix"
TOPO_EXTRACELLULAR   = "extracellular"
TOPO_CYTOPLASMIC     = "cytoplasmic"
TOPO_SIGNAL_PEPTIDE  = "signal_peptide"
TOPO_NMD             = "nmd_absent"


@dataclass
class JunctionTopologyResult:
    """
    Result of topology mapping for a single isoform pair.

    Fields
    ------
    pair_id : str
    junction_peptide : str
        The unique junction-spanning peptide sequence.
    junction_start_aa : int
        0-based start residue in the alternative protein sequence.
    junction_end_aa : int
        0-based end residue (exclusive).
    junction_found : bool
        Whether the junction peptide was located in the alternative protein.
    alt_protein_tm_segments : list[(int,int)]
        Predicted TM helix segments (1-based, start, end).
    alt_seq_signal_peptide : bool
    topology_of_junction : str
        One of: transmembrane_helix, extracellular, cytoplasmic,
        signal_peptide, soluble, nmd_absent, unknown.
    is_surface_accessible : bool
        True if topology is extracellular or transmembrane_helix.
    confidence : float
        0-1 score based on how strongly the topology is predicted.
    evidence_summary : str
        Human-readable description.
    """

    pair_id: str = ""
    junction_peptide: str = ""
    junction_start_aa: int = -1
    junction_end_aa: int = -1
    junction_found: bool = False
    alt_protein_tm_segments: List[Tuple[int, int]] = field(default_factory=list)
    alt_seq_signal_peptide: bool = False
    topology_of_junction: str = TOPO_UNKNOWN
    is_surface_accessible: bool = False
    confidence: float = 0.0
    evidence_summary: str = ""


class JunctionTopologyValidator:
    """
    Maps AltAnalyze3 junction-spanning peptides onto predicted protein topology
    to determine if the junction is surface-accessible.

    The key insight: for surfaceome validation, what matters is not just whether
    the protein overall is a membrane protein, but whether the *unique peptide
    region that distinguishes the isoform* (the junction-spanning sequence) is
    in a transmembrane or extracellular position — making it a meaningful
    structural change from the reference isoform.
    """

    def __init__(self, tm_threshold: float = 1.6, signal_score_min: float = 2.0) -> None:
        self.tm_threshold = tm_threshold
        self.signal_score_min = signal_score_min

    def _assign_topology(
        self,
        start_aa: int,
        end_aa: int,
        tm_segments: List[Tuple[int, int]],
        has_signal: bool,
        seq_len: int,
    ) -> Tuple[str, float]:
        """
        Given a peptide range [start_aa, end_aa) (0-based) and predicted TM
        segments (1-based), determine the topology category of the peptide.

        Logic:
        - If in signal peptide region (aa 1-35 and signal predicted): signal_peptide
        - If ≥50% of residues overlap a TM helix: transmembrane_helix
        - If between TM helices or beyond last TM on the N-terminal side (for type I)
          or C-terminal side (for type II): extracellular
        - If in cytoplasmic loops: cytoplasmic
        - If no TM helices: soluble
        """
        if not tm_segments:
            # No TM helices predicted
            if has_signal and start_aa < 35:
                return TOPO_SIGNAL_PEPTIDE, 0.7
            return TOPO_SOLUBLE, 0.8

        n_residues = end_aa - start_aa
        if n_residues <= 0:
            return TOPO_UNKNOWN, 0.0

        # Convert junction to 1-based for comparison with tm_segments
        j_start_1 = start_aa + 1
        j_end_1   = end_aa          # inclusive end in 1-based

        # Count residues overlapping TM helices
        tm_overlap = 0
        for tm_s, tm_e in tm_segments:
            overlap_start = max(j_start_1, tm_s)
            overlap_end   = min(j_end_1,   tm_e)
            if overlap_end >= overlap_start:
                tm_overlap += overlap_end - overlap_start + 1

        tm_fraction = tm_overlap / n_residues

        if tm_fraction >= 0.5:
            return TOPO_TM_HELIX, min(0.5 + tm_fraction, 1.0)

        # Check signal peptide
        if has_signal and j_start_1 <= 35:
            return TOPO_SIGNAL_PEPTIDE, 0.7

        # Classify as extracellular vs cytoplasmic using a simple alternating
        # topology model. For single-pass type I proteins, the region before the
        # first TM is extracellular; for multi-pass, loops alternate starting
        # with extracellular on the N-terminal side.
        n_tm = len(tm_segments)
        peptide_mid = (j_start_1 + j_end_1) / 2

        # Find which inter-TM loop the peptide falls in
        loop_idx = 0  # 0 = before first TM, 1 = between TM1-TM2, etc.
        for tm_s, tm_e in tm_segments:
            if peptide_mid > tm_e:
                loop_idx += 1
            else:
                break

        # Extracellular vs cytoplasmic assignment:
        # For type I single-pass: N-terminus (loop 0) extracellular, C (loop 1) cytoplasmic
        # For multi-pass: even loops extracellular, odd loops cytoplasmic (canonical)
        is_extracellular = (loop_idx % 2 == 0)
        confidence = 0.55 if n_tm >= 2 else 0.45  # lower confidence for single-pass

        if is_extracellular:
            return TOPO_EXTRACELLULAR, confidence
        else:
            return TOPO_CYTOPLASMIC, confidence

    def validate_one(
        self,
        pair_id: str,
        alternative_protein_seq: str,
        junction_peptide: str,
        has_alt_protein: bool = True,
    ) -> JunctionTopologyResult:
        """
        Validate topology of a single junction peptide.

        Parameters
        ----------
        pair_id : str
        alternative_protein_seq : str
            Full translated sequence of the alternative isoform.
        junction_peptide : str
            The unique amino acid sequence spanning the novel junction.
            Should be 6-30 aa, long enough to be unique.
        has_alt_protein : bool
            False if the alternative isoform is predicted to undergo NMD.
        """
        from proteus.encoders.topology import predict_tm_segments, predict_signal_peptide

        result = JunctionTopologyResult(pair_id=pair_id, junction_peptide=junction_peptide)

        if not has_alt_protein:
            result.topology_of_junction = TOPO_NMD
            result.is_surface_accessible = False
            result.confidence = 1.0
            result.evidence_summary = "NMD predicted — protein absent from surface"
            return result

        seq = alternative_protein_seq.upper().strip()
        pep = junction_peptide.upper().strip()

        if not seq or not pep:
            result.topology_of_junction = TOPO_UNKNOWN
            result.evidence_summary = "Missing sequence or junction peptide"
            return result

        # Locate junction peptide in alternative protein
        pos = seq.find(pep)
        if pos == -1:
            result.junction_found = False
            result.topology_of_junction = TOPO_UNKNOWN
            result.evidence_summary = (
                f"Junction peptide '{pep[:10]}...' not found in alternative protein "
                f"(len={len(seq)}). Check frame or translation."
            )
            return result

        result.junction_found = True
        result.junction_start_aa = pos
        result.junction_end_aa   = pos + len(pep)

        # Predict TM topology of the full alternative protein
        tm_segs = predict_tm_segments(seq, threshold=self.tm_threshold)
        sig_candidate, sig_score = predict_signal_peptide(seq)
        has_signal = sig_score >= self.signal_score_min

        result.alt_protein_tm_segments = tm_segs
        result.alt_seq_signal_peptide  = bool(has_signal)

        # Assign topology
        topo, conf = self._assign_topology(
            pos, pos + len(pep), tm_segs, has_signal, len(seq)
        )
        result.topology_of_junction = topo
        result.confidence = float(conf)
        result.is_surface_accessible = topo in (TOPO_TM_HELIX, TOPO_EXTRACELLULAR, TOPO_SIGNAL_PEPTIDE)

        # Build evidence summary
        tm_str = f"{len(tm_segs)} TM helices" if tm_segs else "no TM helices"
        result.evidence_summary = (
            f"Junction aa {pos+1}-{pos+len(pep)} | "
            f"{tm_str} | "
            f"signal={int(has_signal)} | "
            f"topology={topo} | "
            f"surface_accessible={result.is_surface_accessible}"
        )

        return result

    def validate(self, pairs: List[Dict[str, Any]]) -> List[JunctionTopologyResult]:
        """
        Validate a list of isoform pairs.

        Each dict must have keys:
          pair_id, alternative_protein_seq, junction_peptide
          Optional: has_alt_protein (default True)

        Returns results sorted by is_surface_accessible desc, confidence desc.
        """
        results = []
        for pair in pairs:
            r = self.validate_one(
                pair_id=str(pair.get("pair_id", "")),
                alternative_protein_seq=str(pair.get("alternative_protein_seq", "")),
                junction_peptide=str(pair.get("junction_peptide", "")),
                has_alt_protein=bool(pair.get("has_alt_protein", True)),
            )
            results.append(r)

        results.sort(key=lambda r: (r.is_surface_accessible, r.confidence), reverse=True)
        return results

    def to_dataframe(self, results: List[JunctionTopologyResult]):
        """Convert results to a pandas DataFrame."""
        import pandas as pd
        rows = []
        for r in results:
            rows.append({
                "pair_id": r.pair_id,
                "junction_peptide": r.junction_peptide,
                "junction_start_aa": r.junction_start_aa,
                "junction_end_aa": r.junction_end_aa,
                "junction_found": int(r.junction_found),
                "n_tm_helices": len(r.alt_protein_tm_segments),
                "tm_segments": ";".join(f"{s}-{e}" for s, e in r.alt_protein_tm_segments),
                "alt_seq_signal_peptide": int(r.alt_seq_signal_peptide),
                "topology_of_junction": r.topology_of_junction,
                "is_surface_accessible": int(r.is_surface_accessible),
                "confidence": round(r.confidence, 3),
                "evidence_summary": r.evidence_summary,
            })
        return pd.DataFrame(rows)

    def to_tsv(self, results: List[JunctionTopologyResult], output_path) -> None:
        """Write results to TSV file."""
        df = self.to_dataframe(results)
        df.to_csv(output_path, sep="\t", index=False)
        n_surface = sum(1 for r in results if r.is_surface_accessible)
        print(f"[JunctionTopologyValidator] {len(results)} pairs evaluated, "
              f"{n_surface} with surface-accessible junctions.")
        print(f"  Written to: {output_path}")
