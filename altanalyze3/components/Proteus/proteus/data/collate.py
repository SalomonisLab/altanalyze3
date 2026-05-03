"""
Custom collate function for Proteus DataLoader.

Handles stacking of embedding tensors, task masks, and metadata strings.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import torch
from torch import Tensor


def proteus_collate_fn(batch: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Collate a list of ProteusDataset items into a batched dict.

    All embedding tensors are stacked into [B, D] tensors.
    None embeddings are replaced with zeros.
    String metadata fields (gene_id, transcript IDs) are kept as lists.
    Task masks are stacked into [B] tensors.

    Parameters
    ----------
    batch : list[dict]
        List of dicts from ProteusDataset.__getitem__().

    Returns
    -------
    dict
        Batched tensors and metadata.
    """
    if not batch:
        return {}

    # Tensor fields to stack
    tensor_fields = [
        "ref_rna_emb",
        "alt_rna_emb",
        "ref_protein_emb",
        "alt_protein_emb",
        "ref_disorder_raw",
        "alt_disorder_raw",
        "ref_structural_raw",
        "alt_structural_raw",
        "has_alt_protein",
        "label_preservation",
        "deviation_class",
    ]

    # String/metadata fields
    list_fields = [
        "gene_id",
        "reference_transcript_id",
        "alt_transcript_id",
    ]

    # Expected dims for zero fallback
    dim_defaults = {
        "ref_rna_emb": 512,
        "alt_rna_emb": 512,
        "ref_protein_emb": 1280,
        "alt_protein_emb": 1280,
        "ref_disorder_raw": 50,
        "alt_disorder_raw": 50,
        "ref_structural_raw": 48,
        "alt_structural_raw": 48,
        "has_alt_protein": None,      # scalar
        "label_preservation": None,   # scalar
        "deviation_class": None,      # scalar
    }

    result: Dict[str, Any] = {}

    # Stack tensor fields
    for field in tensor_fields:
        tensors = []
        for item in batch:
            val = item.get(field)
            if val is None:
                # Create zero fallback
                dim = dim_defaults.get(field)
                if dim is not None:
                    val = torch.zeros(dim)
                else:
                    val = torch.tensor(0)
            tensors.append(val if isinstance(val, Tensor) else torch.tensor(val))
        result[field] = torch.stack(tensors, dim=0)

    # Collect string fields as lists
    for field in list_fields:
        result[field] = [item.get(field, "") for item in batch]

    # Stack task masks: each item has a dict of task_name -> scalar tensor
    task_mask_keys = set()
    for item in batch:
        tm = item.get("task_masks", {})
        task_mask_keys.update(tm.keys())

    if task_mask_keys:
        stacked_masks: Dict[str, Tensor] = {}
        for task in task_mask_keys:
            mask_list = []
            for item in batch:
                tm = item.get("task_masks", {})
                mask_list.append(tm.get(task, torch.tensor(0.0)))
            stacked_masks[task] = torch.stack(mask_list, dim=0)  # [B]
        result["task_masks"] = stacked_masks

    return result
