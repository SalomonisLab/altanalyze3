"""
topology.py: Sequence-based transmembrane topology and signal peptide prediction.

Uses the same Kyte-Doolittle hydropathy algorithm as Daedalus
build_gencode_protein_sequence_features.py so it can be applied to novel
AltAnalyze3 isoforms not present in the GENCODE FASTA.

CPU-only, no external dependencies beyond numpy.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

import numpy as np


# Kyte-Doolittle hydropathy scale
_HYDRO: Dict[str, float] = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4,
    "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5,
    "P": -1.6, "Q": -3.5, "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2,
    "W": -0.9, "Y": -1.3,
}
_SMALL_AA = set("AGSCTV")
_GLYCO_RE = re.compile(r"N[^P][ST]")


def _prefix_hydropathy(seq: str) -> List[float]:
    prefix = [0.0]
    total = 0.0
    for aa in seq:
        total += _HYDRO.get(aa, 0.0)
        prefix.append(total)
    return prefix


def _window_hydropathy(prefix: List[float], start: int, end: int) -> float:
    end = min(end, len(prefix) - 1)
    length = end - start
    if length <= 0:
        return 0.0
    return (prefix[end] - prefix[start]) / length


def _merge_segments(segments: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not segments:
        return []
    segments.sort()
    merged = [segments[0]]
    for start, end in segments[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end + 1:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def predict_tm_segments(
    seq: str,
    window: int = 19,
    threshold: float = 1.6,
) -> List[Tuple[int, int]]:
    """
    Predict transmembrane helix segments using a sliding Kyte-Doolittle window.

    Returns list of (start, end) 1-based residue positions.
    """
    prefix = _prefix_hydropathy(seq)
    segments: List[Tuple[int, int]] = []
    for i in range(0, max(len(seq) - window + 1, 0)):
        frag = seq[i: i + window]
        hydro = _window_hydropathy(prefix, i, i + window)
        hydrophobic_fraction = sum(aa in "AILMFWVYC" for aa in frag) / len(frag)
        if hydro >= threshold and hydrophobic_fraction >= 0.68 and "P" not in frag[:15]:
            segments.append((i + 1, i + window))
    return _merge_segments(segments)


def predict_signal_peptide(seq: str) -> Tuple[int, float]:
    """
    Predict signal peptide from N-terminal hydropathy and cleavage motifs.

    Returns (is_signal: int, score: float).
    Score components: positive N-terminal charge, core hydrophobicity, cleavage motif.
    """
    nterm = seq[:35]
    if len(nterm) < 15:
        return 0, 0.0
    positive_n = sum(aa in "KR" for aa in nterm[:5])
    prefix = _prefix_hydropathy(nterm)
    best_hydro = max(
        (_window_hydropathy(prefix, i, i + 8) for i in range(2, max(len(nterm) - 7, 3))),
        default=0.0,
    )
    cleavage_bonus = 0.0
    for cut in range(15, min(len(nterm) - 1, 30)):
        tri = nterm[max(cut - 3, 0): cut]
        if len(tri) == 3 and tri[0] in _SMALL_AA and tri[2] in _SMALL_AA:
            cleavage_bonus = 1.0
            break
    score = 0.0
    if positive_n >= 1:
        score += 1.0
    if best_hydro >= 1.8:
        score += 1.0
    score += cleavage_bonus
    return int(score >= 2.0), score


def predict_topology_type(tm_segments: List[Tuple[int, int]], seq_len: int) -> int:
    """
    Classify predicted membrane topology type.

    Returns:
        0 = soluble (0 TM)
        1 = single-pass type I (N-out, C-in, TM near N-term)
        2 = single-pass type II (N-in, C-out, TM near N-term)
        3 = single-pass type III (TM near C-term)
        4 = multi-pass (≥2 TM)
        5 = GPI-anchored proxy (very C-terminal single TM)
    """
    n = len(tm_segments)
    if n == 0:
        return 0
    if n >= 2:
        return 4
    start, end = tm_segments[0]
    midpoint = (start + end) / 2
    rel = midpoint / max(seq_len, 1)
    if rel < 0.35:
        return 1  # Type I or II
    if rel > 0.80:
        return 5  # GPI proxy / C-terminal anchor
    return 3  # Type III


class SequenceTopologyEncoder:
    """
    Computes a 9-dimensional topology feature vector from a protein sequence.

    Mirrors the Daedalus build_gencode_protein_sequence_features.py algorithm
    so it can be applied to any sequence at inference time, including novel
    AltAnalyze3 isoforms not present in the pre-computed GENCODE tables.

    Feature layout [9]:
    [0] seq_tm_count (log1p of predicted TM helices)
    [1] seq_signal_candidate (binary)
    [2] seq_signal_score (clamped 0-1)
    [3] seq_tm_span_fraction (total TM span / protein length)
    [4] seq_max_hydropathy_19 (clamped 0-1, from window=19 scan)
    [5] seq_n_terminal_hydropathy (mean hydropathy first 8 aa, clamped 0-1)
    [6] seq_glyco_motif_count (log1p)
    [7] seq_cysteine_fraction (cysteine count / length)
    [8] seq_topology_type (one of 0-5, divided by 5 for normalisation)
    """

    FEATURE_DIM: int = 9

    def compute_features(self, seq: str) -> np.ndarray:
        """
        Compute topology features from a protein sequence string.

        Parameters
        ----------
        seq : str
            Protein sequence (uppercase single-letter codes).

        Returns
        -------
        np.ndarray
            Shape [9], dtype float32.
        """
        seq = seq.upper().strip()
        features = np.zeros(self.FEATURE_DIM, dtype=np.float32)

        if not seq:
            return features

        n = len(seq)
        prefix = _prefix_hydropathy(seq)

        # TM prediction
        tm_segs = predict_tm_segments(seq)
        tm_count = len(tm_segs)
        tm_span = sum(e - s for s, e in tm_segs)

        # Signal prediction
        sig_candidate, sig_score = predict_signal_peptide(seq)

        # Max hydropathy in any 19-aa window
        max_hydro_19 = max(
            (_window_hydropathy(prefix, i, i + 19) for i in range(max(n - 18, 1))),
            default=0.0,
        )

        # N-terminal hydropathy (first 8 aa)
        n_term_hydro = _window_hydropathy(prefix, 0, min(8, n))

        # Topology type
        topo_type = predict_topology_type(tm_segs, n)

        features[0] = float(np.log1p(tm_count))
        features[1] = float(sig_candidate)
        features[2] = float(min(sig_score / 3.0, 1.0))           # max score is 3
        features[3] = float(tm_span / n) if n > 0 else 0.0
        features[4] = float(np.clip((max_hydro_19 + 4.5) / 9.0, 0, 1))  # normalise -4.5..4.5 → 0..1
        features[5] = float(np.clip((n_term_hydro + 4.5) / 9.0, 0, 1))
        features[6] = float(np.log1p(len(_GLYCO_RE.findall(seq))))
        features[7] = float(seq.count("C") / n)
        features[8] = float(topo_type / 5.0)

        return features
