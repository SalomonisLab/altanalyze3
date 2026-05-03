"""
Sequence encoding utilities for Proteus.

Provides 6-track RNA encoding and protein sequence processing utilities.
"""

from __future__ import annotations

from typing import List, Optional

import numpy as np


# Nucleotide to one-hot index mapping
_NT_TO_IDX = {
    "A": 0, "a": 0,
    "C": 1, "c": 1,
    "G": 2, "g": 2,
    "U": 3, "u": 3,
    "T": 3, "t": 3,  # T treated as U (RNA)
}


def encode_rna_6track(
    sequence: str,
    splice_sites: Optional[List[int]] = None,
    cds_start: Optional[int] = None,
    cds_end: Optional[int] = None,
) -> np.ndarray:
    """
    Encode an RNA/DNA sequence using the 6-track scheme expected by Orthrus.

    Track layout:
    - Track 0: A indicator (1.0 if position is A, else 0.0)
    - Track 1: C indicator
    - Track 2: G indicator
    - Track 3: U/T indicator
    - Track 4: splice junction indicator (1.0 at splice donor/acceptor positions)
    - Track 5: CDS codon start indicator (1.0 at first nt of each codon in CDS)

    Ambiguous nucleotides (N, -, etc.) are encoded as all-zeros in tracks 0-3.

    Parameters
    ----------
    sequence : str
        RNA or DNA nucleotide sequence. Case-insensitive. T is treated as U.
    splice_sites : list[int] | None
        0-based positions of splice donors and acceptors. Track 4 = 1.0 here.
        If None, track 4 is all zeros.
    cds_start : int | None
        0-based start position of the CDS (first nucleotide of start codon).
    cds_end : int | None
        0-based end position of the CDS (exclusive; position after stop codon).

    Returns
    -------
    np.ndarray
        Shape [6, L], dtype float32.
    """
    L = len(sequence)
    encoding = np.zeros((6, L), dtype=np.float32)

    # Tracks 0-3: one-hot nucleotide encoding
    for i, nt in enumerate(sequence):
        idx = _NT_TO_IDX.get(nt)
        if idx is not None:
            encoding[idx, i] = 1.0
        # N or other ambiguous → all zeros (already initialized to 0)

    # Track 4: splice junction indicators
    if splice_sites:
        for pos in splice_sites:
            if 0 <= pos < L:
                encoding[4, pos] = 1.0

    # Track 5: codon start positions within CDS
    if cds_start is not None and cds_end is not None:
        cds_start = max(0, cds_start)
        cds_end = min(L, cds_end)
        # Mark the first nucleotide of each codon
        for codon_start in range(cds_start, cds_end - 2, 3):
            encoding[5, codon_start] = 1.0

    return encoding


def center_crop_protein(sequence: str, max_length: int = 1024) -> str:
    """
    Center-crop a protein sequence to at most max_length amino acids.

    The crop is centered on the sequence: equal amounts are removed from
    the N-terminus and C-terminus (with the extra residue kept on the
    N-terminal side for even-length differences).

    Parameters
    ----------
    sequence : str
        Amino acid sequence (1-letter codes).
    max_length : int
        Maximum output length. Default 1024 (ESM2 context window).

    Returns
    -------
    str
        Cropped sequence of length <= max_length.

    Examples
    --------
    >>> center_crop_protein("ACDEFGHIKLMNPQRSTVWY" * 100, max_length=1024)
    # returns middle 1024 residues
    """
    L = len(sequence)
    if L <= max_length:
        return sequence

    # Center crop: keep middle max_length residues
    offset = (L - max_length) // 2
    return sequence[offset : offset + max_length]
