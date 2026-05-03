"""
DisorderEncoder: Intrinsic disorder feature extractor.

Attempts to use STARLING first, falls back to metapredict.
Computes a 50-dimensional feature vector summarizing per-residue
disorder scores for downstream use in the Proteus delta model.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import List, Optional

import numpy as np


class DisorderEncoder:
    """
    Computes intrinsic disorder feature vectors from protein sequences.

    Attempts to load STARLING (preferred); falls back to metapredict if
    STARLING is unavailable. Returns zero features if both are unavailable.

    Feature vector [50-dim]:
    - 10 binned disorder fractions (0-0.1, 0.1-0.2, ..., 0.9-1.0)
    - fraction of N-terminal 20% with score > 0.5
    - fraction of C-terminal 20% with score > 0.5
    - max contiguous IDR length (normalized by total length)
    - number of IDR regions (normalized by length)
    - mean disorder score
    - std disorder score
    - fraction with score > 0.5
    - fraction with score > 0.8
    - 10 positional decile means
    - 10 positional decile fractions (> 0.5)
    - 6 reserved zeros (PTM-site disorder, populated downstream)
    Total: 10 + 2 + 2 + 2 + 2 + 10 + 10 + 6 = 44 → padded to 50
    """

    FEATURE_DIM: int = 50

    def __init__(self, backend: str = "auto") -> None:
        """
        Parameters
        ----------
        backend : str
            "auto" tries STARLING then metapredict. "starling" or "metapredict"
            forces a specific backend.
        """
        self.backend = backend
        self._predictor = None
        self._backend_name: Optional[str] = None
        self._load_backend()

    def _load_backend(self) -> None:
        """Load disorder prediction backend."""
        if self.backend in ("auto", "starling"):
            try:
                import starling  # noqa: F401
                self._backend_name = "starling"
                return
            except ImportError:
                if self.backend == "starling":
                    warnings.warn(
                        "STARLING not available. Install with: pip install starling",
                        UserWarning,
                        stacklevel=2,
                    )
                    return

        if self.backend in ("auto", "metapredict"):
            try:
                import metapredict as meta  # noqa: F401
                self._backend_name = "metapredict"
                return
            except ImportError:
                pass

        warnings.warn(
            "Neither STARLING nor metapredict is available. "
            "DisorderEncoder will return zero features. "
            "Install metapredict: pip install metapredict",
            UserWarning,
            stacklevel=2,
        )
        self._backend_name = None

    def predict_disorder(self, sequence: str) -> np.ndarray:
        """
        Predict per-residue disorder scores for a protein sequence.

        Parameters
        ----------
        sequence : str
            Amino acid sequence (1-letter codes).

        Returns
        -------
        np.ndarray
            Shape [L], values in [0, 1]. Higher = more disordered.
            Returns zeros if no backend available.
        """
        if not sequence:
            return np.array([], dtype=np.float32)

        if self._backend_name == "starling":
            try:
                import starling
                scores = starling.predict(sequence)
                return np.array(scores, dtype=np.float32)
            except Exception as exc:
                warnings.warn(f"STARLING prediction failed: {exc}. Using zeros.", stacklevel=2)
                return np.zeros(len(sequence), dtype=np.float32)

        elif self._backend_name == "metapredict":
            try:
                import metapredict as meta
                scores = meta.predict(sequence)
                return np.array(scores, dtype=np.float32)
            except Exception as exc:
                warnings.warn(f"metapredict prediction failed: {exc}. Using zeros.", stacklevel=2)
                return np.zeros(len(sequence), dtype=np.float32)

        return np.zeros(len(sequence), dtype=np.float32)

    def compute_features(
        self,
        sequence: str,
        per_residue_scores: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Compute the 50-dimensional disorder feature vector.

        Parameters
        ----------
        sequence : str
            Amino acid sequence.
        per_residue_scores : np.ndarray | None
            Pre-computed per-residue scores. If None, computed internally.

        Returns
        -------
        np.ndarray
            Shape [50].
        """
        if not sequence:
            return np.zeros(self.FEATURE_DIM, dtype=np.float32)

        if per_residue_scores is None:
            per_residue_scores = self.predict_disorder(sequence)

        scores = per_residue_scores.astype(np.float32)
        L = len(scores)

        if L == 0:
            return np.zeros(self.FEATURE_DIM, dtype=np.float32)

        features: List[float] = []

        # --- 10 binned disorder fractions ---
        bin_edges = np.linspace(0.0, 1.0, 11)
        for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
            frac = float(np.mean((scores >= lo) & (scores < hi)))
            features.append(frac)

        # --- N-terminal and C-terminal disorder fractions ---
        n_term_size = max(1, int(0.2 * L))
        c_term_size = max(1, int(0.2 * L))
        nterm_frac = float(np.mean(scores[:n_term_size] > 0.5))
        cterm_frac = float(np.mean(scores[-c_term_size:] > 0.5))
        features.append(nterm_frac)
        features.append(cterm_frac)

        # --- Max contiguous IDR length (normalized) ---
        binary = (scores > 0.5).astype(int)
        max_idr_len = 0
        n_idr_regions = 0
        current_len = 0
        in_idr = False
        for val in binary:
            if val == 1:
                current_len += 1
                if not in_idr:
                    n_idr_regions += 1
                    in_idr = True
                max_idr_len = max(max_idr_len, current_len)
            else:
                current_len = 0
                in_idr = False

        features.append(float(max_idr_len) / L)
        features.append(float(n_idr_regions) / L)

        # --- Global statistics ---
        features.append(float(np.mean(scores)))            # mean
        features.append(float(np.std(scores)))             # std
        features.append(float(np.mean(scores > 0.5)))      # fraction > 0.5
        features.append(float(np.mean(scores > 0.8)))      # fraction > 0.8

        # --- 10 positional decile means ---
        decile_size = max(1, L // 10)
        for d in range(10):
            start = d * decile_size
            end = (d + 1) * decile_size if d < 9 else L
            chunk = scores[start:end]
            features.append(float(np.mean(chunk)) if len(chunk) > 0 else 0.0)

        # --- 10 positional decile fractions (> 0.5) ---
        for d in range(10):
            start = d * decile_size
            end = (d + 1) * decile_size if d < 9 else L
            chunk = scores[start:end]
            features.append(float(np.mean(chunk > 0.5)) if len(chunk) > 0 else 0.0)

        # --- 6 reserved zeros (PTM-site disorder) ---
        features.extend([0.0] * 6)

        feat_array = np.array(features, dtype=np.float32)

        # Ensure exactly FEATURE_DIM elements
        if len(feat_array) < self.FEATURE_DIM:
            feat_array = np.concatenate([
                feat_array,
                np.zeros(self.FEATURE_DIM - len(feat_array), dtype=np.float32)
            ])
        elif len(feat_array) > self.FEATURE_DIM:
            feat_array = feat_array[: self.FEATURE_DIM]

        return feat_array

    def extract_and_cache(
        self,
        protein_ids: List[str],
        sequences: List[Optional[str]],
        output_dir: Path,
    ) -> None:
        """
        Compute and cache disorder features for a list of proteins.

        Parameters
        ----------
        protein_ids : list[str]
            Protein identifiers.
        sequences : list[str | None]
            Amino acid sequences. None or "" yields zero features.
        output_dir : Path
            Directory to save {protein_id}.npy files.
        """
        from tqdm import tqdm

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        to_process = []
        cached_count = 0
        for pid, seq in zip(protein_ids, sequences):
            out_path = output_dir / f"{pid}.npy"
            if out_path.exists():
                cached_count += 1
            else:
                to_process.append((pid, seq or ""))

        print(
            f"[DisorderEncoder] {cached_count} already cached, "
            f"{len(to_process)} to compute."
        )

        for pid, seq in tqdm(to_process, desc="Computing disorder features"):
            out_path = output_dir / f"{pid}.npy"
            try:
                feats = self.compute_features(seq)
                np.save(str(out_path), feats)
            except Exception as exc:
                warnings.warn(f"Failed to compute disorder features for {pid!r}: {exc}", stacklevel=2)
                np.save(str(out_path), np.zeros(self.FEATURE_DIM, dtype=np.float32))

    @staticmethod
    def load_cached(protein_id: str, cache_dir: Path) -> Optional[np.ndarray]:
        """
        Load cached disorder features.

        Parameters
        ----------
        protein_id : str
            Protein identifier.
        cache_dir : Path
            Directory containing {protein_id}.npy files.

        Returns
        -------
        np.ndarray | None
            Shape [50], or None if not found.
        """
        path = Path(cache_dir) / f"{protein_id}.npy"
        if not path.exists():
            return None
        return np.load(str(path)).astype(np.float32)
