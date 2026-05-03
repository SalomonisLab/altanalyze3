"""
OrthrusRNAEncoder: Frozen Orthrus foundation model wrapper for RNA transcript encoding.

Orthrus uses a 6-track Mamba SSM architecture trained on paired RNA transcripts.
This module handles loading, encoding, caching, and graceful fallback.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import torch

from ..utils.device import get_device
from ..utils.sequences import encode_rna_6track


class OrthrusRNAEncoder:
    """
    Wrapper around the frozen Orthrus RNA foundation model.

    Encodes RNA transcript sequences into 512-dimensional embeddings using
    6-track input encoding (4-hot nucleotide + splice junction + CDS track).

    The model is loaded from HuggingFace and kept frozen (eval mode, no grad).
    If mamba_ssm is unavailable, returns zero tensors with a warning.

    Parameters
    ----------
    model_id : str
        HuggingFace model identifier for Orthrus.
    device : str
        Device specification. "auto" selects CUDA > MPS > CPU.
    max_length : int
        Maximum sequence length in nucleotides. Longer sequences are truncated.
    """

    EMBEDDING_DIM: int = 512
    N_TRACKS: int = 6

    def __init__(
        self,
        model_id: str = "quietflamingo/orthrus-large-6-track",
        device: str = "auto",
        max_length: int = 12288,
    ) -> None:
        self.model_id = model_id
        self.max_length = max_length
        self.device = get_device(device)
        self._model = None
        self._available = False
        self._load_model()

    def _load_model(self) -> None:
        """Attempt to load Orthrus from HuggingFace; fall back gracefully on failure."""
        try:
            import mamba_ssm  # noqa: F401 — required for Orthrus
            from transformers import AutoModel

            self._model = AutoModel.from_pretrained(
                self.model_id,
                trust_remote_code=True,
            )
            self._model.eval()
            self._model.to(self.device)
            for param in self._model.parameters():
                param.requires_grad = False
            self._available = True
        except ImportError:
            warnings.warn(
                "mamba_ssm is not installed. OrthrusRNAEncoder will return zero tensors. "
                "Install mamba_ssm (requires CUDA) for full functionality: "
                "pip install mamba-ssm",
                UserWarning,
                stacklevel=2,
            )
            self._available = False
        except Exception as exc:
            warnings.warn(
                f"Failed to load Orthrus model from {self.model_id!r}: {exc}. "
                "OrthrusRNAEncoder will return zero tensors.",
                UserWarning,
                stacklevel=2,
            )
            self._available = False

    def encode(
        self,
        sequence: str,
        splice_sites: Optional[List[int]] = None,
        cds_start: Optional[int] = None,
        cds_end: Optional[int] = None,
    ) -> torch.Tensor:
        """
        Encode a single RNA transcript sequence.

        Parameters
        ----------
        sequence : str
            RNA/DNA nucleotide sequence (A/C/G/U/T, case-insensitive).
        splice_sites : list[int] | None
            0-based nucleotide positions of splice donors and acceptors.
            Track 4 will be set to 1.0 at these positions.
        cds_start : int | None
            0-based start position of the CDS. Used for track 5 (codon positions).
        cds_end : int | None
            0-based end position of the CDS (exclusive).

        Returns
        -------
        torch.Tensor
            Shape [512]. Zero tensor if model unavailable or sequence empty.
        """
        if not sequence or not self._available:
            if not sequence:
                warnings.warn("Empty sequence passed to OrthrusRNAEncoder.encode()", stacklevel=2)
            return torch.zeros(self.EMBEDDING_DIM)

        sequence = sequence.upper().replace("T", "U")

        if len(sequence) > self.max_length:
            warnings.warn(
                f"Sequence length {len(sequence)} exceeds max_length {self.max_length}. "
                "Truncating from 3' end.",
                UserWarning,
                stacklevel=2,
            )
            sequence = sequence[: self.max_length]

        # Build 6-track encoding: [6, L]
        track = encode_rna_6track(
            sequence,
            splice_sites=splice_sites,
            cds_start=cds_start,
            cds_end=cds_end,
        )  # np.ndarray [6, L]

        # Convert to tensor [1, 6, L] for batch dimension
        x = torch.from_numpy(track).float().unsqueeze(0).to(self.device)  # [1, 6, L]

        with torch.no_grad():
            # Orthrus expects [B, C, L] input and returns per-position hidden states
            output = self._model(x)  # dict or tensor
            if isinstance(output, dict):
                hidden = output.get("last_hidden_state", output.get("hidden_states"))
            elif hasattr(output, "last_hidden_state"):
                hidden = output.last_hidden_state
            else:
                hidden = output

            # Mean pool over sequence positions → [512]
            emb = hidden.squeeze(0).mean(dim=0)  # [512]

        return emb.cpu()

    def extract_and_cache(
        self,
        transcript_ids: List[str],
        sequences: List[str],
        output_dir: Path,
        device: str = "auto",
        batch_size: int = 32,
        splice_sites_map: Optional[Dict[str, List[int]]] = None,
        cds_map: Optional[Dict[str, tuple]] = None,
    ) -> Dict[str, Path]:
        """
        Batch encode transcripts and save to disk, skipping already-cached IDs.

        Parameters
        ----------
        transcript_ids : list[str]
            Transcript identifiers (used as filenames).
        sequences : list[str]
            Corresponding nucleotide sequences.
        output_dir : Path
            Directory to save {transcript_id}.pt files.
        device : str
            Device override for this call.
        batch_size : int
            Number of sequences to process per batch.
        splice_sites_map : dict | None
            Optional mapping transcript_id -> list of splice site positions.
        cds_map : dict | None
            Optional mapping transcript_id -> (cds_start, cds_end).

        Returns
        -------
        dict[str, Path]
            Mapping of transcript_id -> saved .pt file path.
        """
        from tqdm import tqdm

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        cached: Dict[str, Path] = {}
        to_encode: List[tuple] = []

        for tid, seq in zip(transcript_ids, sequences):
            out_path = output_dir / f"{tid}.pt"
            if out_path.exists():
                cached[tid] = out_path
            else:
                to_encode.append((tid, seq))

        print(
            f"[OrthrusRNAEncoder] {len(cached)} already cached, "
            f"{len(to_encode)} to encode."
        )

        if not to_encode or not self._available:
            if not self._available and to_encode:
                warnings.warn(
                    "Orthrus model unavailable. Saving zero embeddings for all uncached transcripts.",
                    UserWarning,
                    stacklevel=2,
                )
                for tid, _ in to_encode:
                    out_path = output_dir / f"{tid}.pt"
                    torch.save(torch.zeros(self.EMBEDDING_DIM), out_path)
                    cached[tid] = out_path
            return cached

        for i in tqdm(range(0, len(to_encode), batch_size), desc="Encoding RNA"):
            batch = to_encode[i : i + batch_size]
            for tid, seq in batch:
                out_path = output_dir / f"{tid}.pt"
                splice_sites = (splice_sites_map or {}).get(tid)
                cds_coords = (cds_map or {}).get(tid)
                cds_start = cds_coords[0] if cds_coords else None
                cds_end = cds_coords[1] if cds_coords else None

                try:
                    emb = self.encode(seq, splice_sites=splice_sites, cds_start=cds_start, cds_end=cds_end)
                    torch.save(emb, out_path)
                    cached[tid] = out_path
                except Exception as exc:
                    warnings.warn(f"Failed to encode transcript {tid!r}: {exc}", stacklevel=2)

        return cached

    @staticmethod
    def load_cached(transcript_id: str, cache_dir: Path) -> Optional[torch.Tensor]:
        """
        Load a cached embedding.

        Parameters
        ----------
        transcript_id : str
            Transcript identifier.
        cache_dir : Path
            Directory containing {transcript_id}.pt files.

        Returns
        -------
        torch.Tensor | None
            [512] embedding, or None if not found.
        """
        path = Path(cache_dir) / f"{transcript_id}.pt"
        if not path.exists():
            return None
        return torch.load(path, weights_only=True)
