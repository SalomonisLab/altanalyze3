"""
ProteusDataset: Loads isoform pairs and their cached embeddings.

Reads from a TSV file and lazily loads pre-cached embeddings from disk.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import torch
from torch.utils.data import Dataset


def _impute_alt_structural(
    ref_str: "torch.Tensor",
    ref_dis: "torch.Tensor",
    alt_dis: "torch.Tensor",
) -> "torch.Tensor":
    """
    Impute structural features for a novel alternative isoform that has no
    cached ENSP entry in the Daedalus tables.

    Strategy: scale annotation-derived features from the reference by an
    estimated domain-retention fraction derived from the disorder profile
    length ratio (a proxy for protein length ratio when actual lengths
    aren't stored in the feature vectors).

    Features that are NOT scaled (binary gene-level or sequence-derived):
      [5]  has_signal_peptide
      [6]  is_kinase
      [7]  is_transcription_factor
      [8]  is_membrane_protein
      [19-27] sequence-predicted topology block
      [45-46] HPA binary flags

    Parameters
    ----------
    ref_str : Tensor [96]
    ref_dis : Tensor [50]  — reference disorder features (contains length proxy)
    alt_dis : Tensor [50]  — alternative disorder features

    Returns
    -------
    Tensor [96]  — imputed alt structural features
    """
    import torch

    # Estimate length ratio from disorder vector length proxy
    # DisorderEncoder feature [0] = mean disorder (0-1, not length-dependent)
    # Feature [2] = n_idr_count (correlates with length)
    # Best proxy: ratio of sums of disorder feature magnitudes
    ref_magnitude = ref_dis.norm().item()
    alt_magnitude = alt_dis.norm().item()
    if ref_magnitude > 1e-6 and alt_magnitude > 1e-6:
        length_ratio = float(torch.clamp(
            torch.tensor(alt_magnitude / ref_magnitude), 0.1, 1.5
        ))
    else:
        length_ratio = 1.0  # no info, assume full retention

    alt_str = ref_str.clone()

    # Indices to NOT scale (binary gene-level or sequence-derived)
    BINARY_OR_SEQ_INDICES = {5, 6, 7, 8, 19, 20, 21, 22, 23, 24, 25, 26, 27, 45, 46}

    scale = torch.ones(ref_str.shape[0])
    for i in range(ref_str.shape[0]):
        if i not in BINARY_OR_SEQ_INDICES:
            scale[i] = length_ratio

    alt_str = ref_str * scale

    # Zero out sequence-predicted topology block [19-27] — it must come from
    # the actual alt sequence, not be copied from the reference
    alt_str[19:28] = 0.0

    return alt_str


@dataclass
class ProteusRecord:
    """
    Typed container for a single isoform pair record.

    Attributes
    ----------
    reference_transcript_id : str
    alt_transcript_id : str
    reference_protein_id : str
    alt_protein_id : str
    gene_id : str
    gene_name : str
    split : str
        "train", "val", or "test".
    label_preservation : int
        1 = preserved, 0 = not preserved, -1 = unlabeled.
    has_alt_protein : bool
        False if alternative isoform is NMD or non-coding.
    is_membrane_pair : bool
    is_surface_pair : bool
    is_kinase_pair : bool
    is_tf_pair : bool
    is_localization_pair : bool
        True if protein has known subcellular localization (HPA/UniProt).
    is_signaling_pair : bool
        True if protein is a receptor tyrosine kinase, GPCR, or other
        transmembrane signaling receptor.
    is_ppi_pair : bool
        True if protein has known PPI binding interfaces (UniProt or BioGRID).
    deviation_class : int
        0-7 deviation class, or -1 if unknown.
    label_localization : int
        1 = localization maintained, 0 = changed, -1 = unknown.
    label_tm_insertion : int
        1 = TM insertion competent, 0 = disrupted, -1 = unknown.
    label_signaling : int
        1 = signaling competent, 0 = disrupted, -1 = unknown.
    label_ppi : int
        1 = ≥1 PPI interface retained, 0 = all interfaces lost, -1 = unknown.
    """

    reference_transcript_id: str = ""
    alt_transcript_id: str = ""
    reference_protein_id: str = ""
    alt_protein_id: str = ""
    gene_id: str = ""
    gene_name: str = ""
    split: str = "train"
    label_preservation: int = -1
    has_alt_protein: bool = True
    # Original task masks
    is_membrane_pair: bool = False
    is_surface_pair: bool = False
    is_kinase_pair: bool = False
    is_tf_pair: bool = False
    # New task masks (v1.1)
    is_localization_pair: bool = False
    is_signaling_pair: bool = False
    is_ppi_pair: bool = False
    deviation_class: int = -1
    # New labels (v1.1)
    label_localization: int = -1
    label_tm_insertion: int = -1
    label_signaling: int = -1
    label_ppi: int = -1
    extra: Dict[str, Any] = field(default_factory=dict)


class ProteusDataset(Dataset):
    """
    Dataset of isoform pairs for Proteus training/evaluation.

    Loads pre-cached embeddings lazily from disk. Missing embeddings are
    replaced with zero tensors of the correct shape.

    Parameters
    ----------
    tsv_path : Path
        Path to TSV file (proteus_train.tsv, proteus_val.tsv, etc.).
    cache_dir : Path
        Root directory containing rna_embeddings/, protein_embeddings/,
        disorder_features/, structural_features/ subdirectories.
    split : str | None
        If provided, filter rows to only this split (e.g., "train").
    """

    RNA_DIM = 512
    PROTEIN_DIM = 1280
    DISORDER_DIM = 50
    STRUCTURAL_DIM = 96

    def __init__(
        self,
        tsv_path: Path,
        cache_dir: Path,
        split: Optional[str] = None,
    ) -> None:
        import pandas as pd

        self.cache_dir = Path(cache_dir)
        self.rna_cache = self.cache_dir / "rna_embeddings"
        self.protein_cache = self.cache_dir / "protein_embeddings"
        self.disorder_cache = self.cache_dir / "disorder_features"
        self.structural_cache = self.cache_dir / "structural_features"

        df = pd.read_csv(tsv_path, sep="\t", low_memory=False)

        if split is not None and "split" in df.columns:
            df = df[df["split"] == split].reset_index(drop=True)

        self.records: List[ProteusRecord] = []
        for _, row in df.iterrows():
            rec = self._parse_row(row)
            self.records.append(rec)

        print(
            f"[ProteusDataset] Loaded {len(self.records)} records "
            f"from {tsv_path}" + (f" (split={split})" if split else "")
        )

    @staticmethod
    def _safe_str(val, default: str = "") -> str:
        if val is None or (isinstance(val, float) and np.isnan(val)):
            return default
        return str(val).strip()

    @staticmethod
    def _safe_int(val, default: int = -1) -> int:
        try:
            if val is None or (isinstance(val, float) and np.isnan(val)):
                return default
            return int(val)
        except (ValueError, TypeError):
            return default

    @staticmethod
    def _safe_bool(val, default: bool = False) -> bool:
        if val is None or (isinstance(val, float) and np.isnan(val)):
            return default
        if isinstance(val, bool):
            return val
        if isinstance(val, (int, float)):
            return bool(val)
        return str(val).lower().strip() in ("1", "true", "yes", "y")

    def _parse_row(self, row) -> ProteusRecord:
        """Parse a DataFrame row into a ProteusRecord."""
        # Determine label
        weak_label = self._safe_str(row.get("weak_label", ""))
        if weak_label in ("1", "preserved", "reference_preferred"):
            label = 1
        elif weak_label in ("0", "not_preserved", "diverged"):
            label = 0
        else:
            try:
                label = int(float(weak_label))
            except (ValueError, TypeError):
                label = -1

        # NMD/absent protein detection
        has_alt_protein = True
        nmd_col = row.get("nmd_triggered", row.get("is_nmd", None))
        if self._safe_bool(nmd_col):
            has_alt_protein = False
        alt_prot_id = self._safe_str(
            row.get("alternative_protein_id", row.get("alt_protein_id", ""))
        )
        if not alt_prot_id or alt_prot_id.lower() in ("nan", "none", ""):
            has_alt_protein = False

        return ProteusRecord(
            reference_transcript_id=self._safe_str(row.get("reference_transcript_id", "")),
            alt_transcript_id=self._safe_str(
                row.get("alternative_transcript_id", row.get("alt_transcript_id", ""))
            ),
            reference_protein_id=self._safe_str(row.get("reference_protein_id", "")),
            alt_protein_id=alt_prot_id,
            gene_id=self._safe_str(row.get("gene_id", "")),
            gene_name=self._safe_str(row.get("gene_name", "")),
            split=self._safe_str(row.get("split", "train")),
            label_preservation=label,
            has_alt_protein=has_alt_protein,
            # Original task masks
            is_membrane_pair=self._safe_bool(row.get("is_membrane_pair", False)),
            is_surface_pair=self._safe_bool(row.get("is_surface_pair", False)),
            is_kinase_pair=self._safe_bool(row.get("is_kinase_pair", False)),
            is_tf_pair=self._safe_bool(row.get("is_tf_pair", False)),
            # New task masks (v1.1)
            is_localization_pair=self._safe_bool(row.get("is_localization_pair", False)),
            is_signaling_pair=self._safe_bool(row.get("is_signaling_pair", False)),
            is_ppi_pair=self._safe_bool(row.get("is_ppi_pair", False)),
            deviation_class=self._safe_int(row.get("deviation_class", -1)),
            # New labels (v1.1)
            label_localization=self._safe_int(row.get("label_localization", -1)),
            label_tm_insertion=self._safe_int(row.get("label_tm_insertion", -1)),
            label_signaling=self._safe_int(row.get("label_signaling", -1)),
            label_ppi=self._safe_int(row.get("label_ppi", -1)),
        )

    def _load_rna(self, transcript_id: str) -> torch.Tensor:
        """Load RNA embedding or return zeros."""
        if not transcript_id:
            return torch.zeros(self.RNA_DIM)
        path = self.rna_cache / f"{transcript_id}.pt"
        if path.exists():
            return torch.load(path, weights_only=True)
        return torch.zeros(self.RNA_DIM)

    def _load_protein(self, protein_id: str) -> tuple[torch.Tensor, bool]:
        """Load protein embedding. Returns (embedding, has_protein)."""
        if not protein_id:
            return torch.zeros(self.PROTEIN_DIM), False
        pt_path = self.protein_cache / f"{protein_id}.pt"
        absent_path = self.protein_cache / f"{protein_id}.absent"
        if absent_path.exists():
            return torch.zeros(self.PROTEIN_DIM), False
        if pt_path.exists():
            return torch.load(pt_path, weights_only=True), True
        return torch.zeros(self.PROTEIN_DIM), True  # missing but not marked absent

    def _load_disorder(self, protein_id: str) -> torch.Tensor:
        """Load disorder features or return zeros."""
        if not protein_id:
            return torch.zeros(self.DISORDER_DIM)
        path = self.disorder_cache / f"{protein_id}.npy"
        if path.exists():
            arr = np.load(str(path)).astype(np.float32)
            return torch.from_numpy(arr)
        return torch.zeros(self.DISORDER_DIM)

    def _load_structural(self, protein_id: str) -> torch.Tensor:
        """Load structural features or return zeros."""
        if not protein_id:
            return torch.zeros(self.STRUCTURAL_DIM)
        path = self.structural_cache / f"{protein_id}.pt"
        if path.exists():
            t = torch.load(path, weights_only=True)
            if isinstance(t, torch.Tensor):
                return t.float()
            return torch.from_numpy(t).float()
        # Also check .npy format (some extraction scripts write numpy)
        npy_path = self.structural_cache / f"{protein_id}.npy"
        if npy_path.exists():
            return torch.from_numpy(np.load(str(npy_path)).astype(np.float32))
        return torch.zeros(self.STRUCTURAL_DIM)

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, idx: int) -> Dict[str, Any]:
        """
        Load all embeddings for a single isoform pair.

        Returns
        -------
        dict
            Keys: ref_rna_emb, alt_rna_emb, ref_protein_emb, alt_protein_emb,
            ref_disorder_raw, alt_disorder_raw, ref_structural_raw, alt_structural_raw,
            has_alt_protein, label_preservation, deviation_class,
            task_masks (dict), gene_id, reference_transcript_id, alt_transcript_id
        """
        rec = self.records[idx]

        ref_rna = self._load_rna(rec.reference_transcript_id)
        alt_rna = self._load_rna(rec.alt_transcript_id)

        ref_prot, ref_has_prot = self._load_protein(rec.reference_protein_id)
        alt_prot, alt_has_prot = self._load_protein(rec.alt_protein_id)

        # Override has_alt_protein with flag from metadata
        has_alt_protein = rec.has_alt_protein and alt_has_prot

        ref_dis = self._load_disorder(rec.reference_protein_id)
        alt_dis = self._load_disorder(rec.alt_protein_id)

        ref_str = self._load_structural(rec.reference_protein_id)
        alt_str = self._load_structural(rec.alt_protein_id)

        # Structural feature imputation for novel alternative isoforms
        # ---------------------------------------------------------------
        # If the alt protein has no cached structural features (ENSP not in
        # Daedalus tables — e.g., a novel AltAnalyze3 isoform), we impute
        # the annotation-derived features [0-18, 28-95] from the reference
        # scaled by the length ratio (alt_len / ref_len), while keeping the
        # sequence-predicted topology block [19-27] at zero (it should be
        # re-computed from the alt sequence using extract_sequence_topology_features.py).
        #
        # Rationale: a novel isoform that retains 70% of the reference length
        # likely retains ~70% of its PPI partners, domain count, etc. This is
        # much more informative than an all-zero vector, which makes the loss
        # gate fire indiscriminately on every annotation feature regardless of
        # whether any specific domain was actually disrupted.
        #
        # Features excluded from scaling (they are binary or sequence-derived):
        #   [5]  has_signal_peptide — binary, not length-proportional
        #   [6]  is_kinase — binary, gene-level
        #   [7]  is_transcription_factor — binary
        #   [8]  is_membrane_protein — binary
        #   [19-27] sequence-predicted topology — should come from alt sequence
        #   [45-46] HPA extracellular / membrane — binary, gene-level
        if alt_str.sum() == 0 and ref_str.sum() > 0 and has_alt_protein:
            alt_str = _impute_alt_structural(ref_str, ref_dis, alt_dis)

        # Task masks: 1.0 if this task applies to this sample, 0.0 if not
        task_masks: Dict[str, torch.Tensor] = {
            "global_preservation":     torch.tensor(1.0),
            "topology_preserved":      torch.tensor(float(rec.is_membrane_pair)),
            "surface_retained":        torch.tensor(float(rec.is_surface_pair)),
            "kinase_competent":        torch.tensor(float(rec.is_kinase_pair)),
            "tf_competent":            torch.tensor(float(rec.is_tf_pair)),
            "disorder_preserved":      torch.tensor(1.0),
            # New v1.1 masks
            "localization_preserved":  torch.tensor(float(rec.is_localization_pair)),
            "tm_insertion_competent":  torch.tensor(float(rec.is_membrane_pair)),
            "signaling_competent":     torch.tensor(float(rec.is_signaling_pair)),
            "ppi_interface_preserved": torch.tensor(float(rec.is_ppi_pair)),
        }

        return {
            "ref_rna_emb": ref_rna,
            "alt_rna_emb": alt_rna,
            "ref_protein_emb": ref_prot,
            "alt_protein_emb": alt_prot,
            "ref_disorder_raw": ref_dis,
            "alt_disorder_raw": alt_dis,
            "ref_structural_raw": ref_str,
            "alt_structural_raw": alt_str,
            "has_alt_protein": torch.tensor(has_alt_protein, dtype=torch.bool),
            "label_preservation":  torch.tensor(rec.label_preservation, dtype=torch.long),
            "deviation_class":     torch.tensor(rec.deviation_class, dtype=torch.long),
            # New v1.1 labels
            "label_localization":  torch.tensor(rec.label_localization, dtype=torch.long),
            "label_tm_insertion":  torch.tensor(rec.label_tm_insertion, dtype=torch.long),
            "label_signaling":     torch.tensor(rec.label_signaling, dtype=torch.long),
            "label_ppi":           torch.tensor(rec.label_ppi, dtype=torch.long),
            "task_masks": task_masks,
            "gene_id": rec.gene_id,
            "gene_name": rec.gene_name,
            "reference_transcript_id": rec.reference_transcript_id,
            "alt_transcript_id": rec.alt_transcript_id,
        }
