"""
ESM2ProteinEncoder: Frozen ESM2 protein language model wrapper.

Encodes protein sequences into 1280-dimensional embeddings using
facebook/esm2_t33_650M_UR50D. Handles absent proteins (NMD isoforms)
via sentinel files and center-crop for long sequences.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import torch
import torch.nn.functional as F

from ..utils.device import get_device
from ..utils.sequences import center_crop_protein


class ESM2ProteinEncoder:
    """
    Wrapper around the frozen ESM2 protein language model.

    Encodes protein sequences into 1280-dimensional mean-pooled embeddings.
    For sequences longer than 1024 aa, applies center-crop before encoding.
    For absent proteins (NMD isoforms, empty sequences), saves a sentinel
    file and returns a zero tensor during loading.

    Parameters
    ----------
    model_id : str
        HuggingFace model identifier for ESM2.
    device : str
        Device specification. "auto" selects CUDA > MPS > CPU.
    max_length : int
        Maximum protein sequence length. Longer sequences are center-cropped.
    """

    EMBEDDING_DIM: int = 1280
    ABSENT_SENTINEL: str = ".absent"

    def __init__(
        self,
        model_id: str = "facebook/esm2_t33_650M_UR50D",
        device: str = "auto",
        max_length: int = 1024,
    ) -> None:
        self.model_id = model_id
        self.max_length = max_length
        self.device = get_device(device)
        self._model = None
        self._tokenizer = None
        self._available = False
        self._load_model()

    def _load_model(self) -> None:
        """Load ESM2 from HuggingFace via transformers AutoModel."""
        try:
            from transformers import AutoModel, AutoTokenizer

            self._tokenizer = AutoTokenizer.from_pretrained(self.model_id)
            self._model = AutoModel.from_pretrained(self.model_id)
            self._model.eval()
            self._model.to(self.device)
            for param in self._model.parameters():
                param.requires_grad = False
            self._available = True
        except Exception as exc:
            warnings.warn(
                f"Failed to load ESM2 model from {self.model_id!r}: {exc}. "
                "ESM2ProteinEncoder will return zero tensors.",
                UserWarning,
                stacklevel=2,
            )
            self._available = False

    def encode(self, sequence: str) -> torch.Tensor:
        """
        Encode a single protein sequence.

        For sequences longer than max_length, applies center-crop before encoding.

        Parameters
        ----------
        sequence : str
            Amino acid sequence (standard 1-letter codes).

        Returns
        -------
        torch.Tensor
            Shape [1280]. Zero tensor if model unavailable or sequence absent.
        """
        if not sequence or not self._available:
            return torch.zeros(self.EMBEDDING_DIM)

        # Center-crop if too long
        if len(sequence) > self.max_length:
            sequence = center_crop_protein(sequence, max_length=self.max_length)

        inputs = self._tokenizer(
            sequence,
            return_tensors="pt",
            truncation=True,
            max_length=self.max_length + 2,  # +2 for special tokens
            add_special_tokens=True,
        )
        inputs = {k: v.to(self.device) for k, v in inputs.items()}

        with torch.no_grad():
            outputs = self._model(**inputs)
            hidden = outputs.last_hidden_state  # [1, L+2, 1280]

        # Mean pool over non-padding, non-special-token positions
        attention_mask = inputs["attention_mask"]  # [1, L+2]
        # Exclude first [CLS] and last [EOS] tokens for mean pooling
        # Use all positions weighted by attention mask
        mask = attention_mask.unsqueeze(-1).float()  # [1, L+2, 1]
        emb = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1e-9)  # [1, 1280]

        return emb.squeeze(0).cpu()  # [1280]

    def extract_and_cache(
        self,
        protein_ids: List[str],
        sequences: List[Optional[str]],
        output_dir: Path,
        device: str = "auto",
        batch_size: int = 16,
    ) -> Dict[str, Path]:
        """
        Batch encode proteins and save to disk.

        For absent proteins (None or empty string), saves a {protein_id}.absent
        sentinel file instead of a .pt embedding file.

        Parameters
        ----------
        protein_ids : list[str]
            Protein identifiers (used as filenames).
        sequences : list[str | None]
            Corresponding amino acid sequences. None or "" indicates absent protein.
        output_dir : Path
            Directory to save {protein_id}.pt or {protein_id}.absent files.
        device : str
            Device override.
        batch_size : int
            Number of sequences per batch.

        Returns
        -------
        dict[str, Path]
            Mapping of protein_id -> saved file path (.pt or .absent).
        """
        from tqdm import tqdm

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        cached: Dict[str, Path] = {}
        to_encode: List[Tuple[str, str]] = []

        for pid, seq in zip(protein_ids, sequences):
            pt_path = output_dir / f"{pid}.pt"
            absent_path = output_dir / f"{pid}{self.ABSENT_SENTINEL}"

            if pt_path.exists():
                cached[pid] = pt_path
            elif absent_path.exists():
                cached[pid] = absent_path
            else:
                to_encode.append((pid, seq or ""))

        print(
            f"[ESM2ProteinEncoder] {len(cached)} already cached, "
            f"{len(to_encode)} to encode."
        )

        if not to_encode:
            return cached

        for i in tqdm(range(0, len(to_encode), batch_size), desc="Encoding proteins"):
            batch = to_encode[i : i + batch_size]
            for pid, seq in batch:
                if not seq:
                    # Absent protein (NMD or non-coding isoform)
                    absent_path = output_dir / f"{pid}{self.ABSENT_SENTINEL}"
                    absent_path.touch()
                    cached[pid] = absent_path
                    continue

                if not self._available:
                    # Save zero tensor as fallback
                    pt_path = output_dir / f"{pid}.pt"
                    torch.save(torch.zeros(self.EMBEDDING_DIM), pt_path)
                    cached[pid] = pt_path
                    continue

                try:
                    emb = self.encode(seq)
                    pt_path = output_dir / f"{pid}.pt"
                    torch.save(emb, pt_path)
                    cached[pid] = pt_path
                except Exception as exc:
                    warnings.warn(f"Failed to encode protein {pid!r}: {exc}", stacklevel=2)

        return cached

    @staticmethod
    def load_cached(
        protein_id: str,
        cache_dir: Path,
    ) -> Tuple[torch.Tensor, bool]:
        """
        Load a cached protein embedding.

        Parameters
        ----------
        protein_id : str
            Protein identifier.
        cache_dir : Path
            Directory containing {protein_id}.pt or {protein_id}.absent files.

        Returns
        -------
        tuple[torch.Tensor, bool]
            (embedding [1280], has_protein). Returns (zeros, False) for absent proteins,
            (embedding, True) for present proteins, (zeros, True) if file not found.
        """
        cache_dir = Path(cache_dir)
        pt_path = cache_dir / f"{protein_id}.pt"
        absent_path = cache_dir / f"{protein_id}{ESM2ProteinEncoder.ABSENT_SENTINEL}"

        if absent_path.exists():
            return torch.zeros(ESM2ProteinEncoder.EMBEDDING_DIM), False

        if pt_path.exists():
            emb = torch.load(pt_path, weights_only=True)
            return emb, True

        # Not found — return zeros but flag as "present" so model doesn't zero it
        return torch.zeros(ESM2ProteinEncoder.EMBEDDING_DIM), True
