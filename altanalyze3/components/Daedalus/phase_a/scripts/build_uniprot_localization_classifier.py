#!/usr/bin/env python3
"""DeepLoc-style subcellular-localization classifier trained on UniProt
entries' subcellular_location annotations.

Workflow:
    1. Mine all human/mouse UniProt entries with subcellular_locations
       comments. Map each entry's primary localization to one of 10 classes.
    2. Get ESM2 t30 mean-pooled embedding (640-dim) for each entry's canonical
       sequence. We reuse the existing isoform-embeddings parquet for
       canonical isoforms (those covering primary_accessions); we extract
       additional embeddings on the fly for accessions not yet seen.
    3. Train a small MLP (640 -> 256 -> 10) with class-weighted cross-entropy
       on a stratified hold-out.
    4. Apply the trained head to all 5,882 catalog isoforms (using their
       reconstructed-protein ESM2 t30 embeddings) and emit per-isoform
       10-class probabilities as features.

Outputs:
    phase_a/data/interim/uniprot_localization_classifier.pt
    phase_a/data/interim/uniprot_isoform_localization_probs.tsv
"""
from __future__ import annotations

import argparse
import gzip
import json
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, TensorDataset

PHASE_A = Path(__file__).resolve().parents[1] / "data"
RAW_JSONL = PHASE_A / "raw" / "uniprot_reviewed_human_mouse.jsonl.gz"
CANONICAL_EMB_PARQUET = PHASE_A / "interim" / "uniprot_canonical_esm2_t30_embeddings.parquet"
ISOFORM_EMB_PARQUET = PHASE_A / "interim" / "uniprot_isoform_esm2_t30_embeddings.parquet"
OUT_HEAD = PHASE_A / "interim" / "uniprot_localization_classifier.pt"
OUT_PROBS = PHASE_A / "interim" / "uniprot_isoform_localization_probs.tsv"
ESM_MODEL_NAME = "facebook/esm2_t30_150M_UR50D"

LOC_CLASSES = [
    "cell_membrane",
    "cytoplasm",
    "nucleus",
    "endoplasmic_reticulum",
    "golgi",
    "mitochondrion",
    "lysosome",
    "secreted",
    "other_membrane",  # ER/Golgi/lysosome MEMBRANE — TM proteins not at PM
    "other",
]


def _classify_loc(loc_str: str) -> str:
    s = (loc_str or "").lower()
    if not s:
        return "other"
    # Order matters: PM-likes first, then organellar membranes, then soluble
    if "cell membrane" in s or "plasma membrane" in s or "cell surface" in s or "apical" in s or "basolateral" in s:
        return "cell_membrane"
    if any(t in s for t in ("endoplasmic reticulum membrane", "er membrane")):
        return "other_membrane"
    if "golgi" in s and "membrane" in s:
        return "other_membrane"
    if any(t in s for t in ("lysosome membrane", "endosome membrane", "vesicle membrane",
                             "mitochondrion membrane", "mitochondrial membrane",
                             "lipid droplet", "peroxisome membrane")):
        return "other_membrane"
    if "endoplasmic reticulum" in s:
        return "endoplasmic_reticulum"
    if "golgi" in s:
        return "golgi"
    if "mitochondri" in s:
        return "mitochondrion"
    if "lysosom" in s or "endosom" in s:
        return "lysosome"
    if "nucleus" in s or "nuclear" in s or "chromosome" in s:
        return "nucleus"
    if "secreted" in s or "extracellular" in s:
        return "secreted"
    if "cytoplasm" in s or "cytosol" in s or "cytoskeleton" in s:
        return "cytoplasm"
    return "other"


def _load_training_labels() -> dict[str, str]:
    """Map primary_accession -> localization class using the FIRST listed
    location string (the "primary" location)."""
    labels: dict[str, str] = {}
    with gzip.open(RAW_JSONL, "rt") as f:
        for line in f:
            d = json.loads(line)
            acc = d.get("primaryAccession")
            primary_loc = None
            for c in d.get("comments") or []:
                if c.get("commentType") != "SUBCELLULAR LOCATION":
                    continue
                for loc in c.get("subcellularLocations") or []:
                    v = (loc.get("location") or {}).get("value")
                    if v:
                        primary_loc = v
                        break
                if primary_loc:
                    break
            if primary_loc:
                labels[acc] = _classify_loc(primary_loc)
    return labels


class _LocHead(nn.Module):
    def __init__(self, in_dim: int = 640, hid: int = 256, n_classes: int = len(LOC_CLASSES)):
        super().__init__()
        self.body = nn.Sequential(
            nn.Linear(in_dim, hid),
            nn.LayerNorm(hid),
            nn.GELU(),
            nn.Dropout(0.2),
            nn.Linear(hid, hid // 2),
            nn.LayerNorm(hid // 2),
            nn.GELU(),
            nn.Dropout(0.2),
            nn.Linear(hid // 2, n_classes),
        )

    def forward(self, x):
        return self.body(x)


def _build_canonical_embedding_table(labels: dict[str, str]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Use the freshly-extracted 37k canonical UniProt embeddings keyed by
    primary_accession. Each entry already carries `primary_loc`."""
    df = pd.read_parquet(CANONICAL_EMB_PARQUET)
    esm_cols = [c for c in df.columns if c.startswith("esm2_")]

    accs = []; embs = []; ys = []
    for _, r in df.iterrows():
        acc = str(r["primary_accession"])
        # Prefer the parquet's own primary_loc, fall back to mined labels
        loc = str(r.get("primary_loc") or "") or labels.get(acc, "")
        cls = _classify_loc(loc) if loc else "other"
        if cls not in LOC_CLASSES:
            continue
        accs.append(acc)
        embs.append(np.array([r[c] for c in esm_cols], dtype=np.float32))
        ys.append(LOC_CLASSES.index(cls))
    X = np.stack(embs)
    y = np.array(ys, dtype=np.int64)
    print(f"[info] training set: {len(accs)} accessions  "
          f"(class dist: {dict(Counter(LOC_CLASSES[i] for i in y).most_common())})")
    return X, y, accs


def _train(X: np.ndarray, y: np.ndarray, epochs: int = 30, lr: float = 1e-3, batch: int = 64) -> _LocHead:
    rng = np.random.default_rng(0)
    idx_train, idx_val = train_test_split(np.arange(len(y)), test_size=0.15, stratify=y, random_state=0)
    X_train, y_train = X[idx_train], y[idx_train]
    X_val, y_val = X[idx_val], y[idx_val]

    cnt = Counter(y_train.tolist())
    weights = np.array([1.0 / max(cnt.get(i, 1), 1) for i in range(len(LOC_CLASSES))], dtype=np.float32)
    weights = weights / weights.sum() * len(LOC_CLASSES)
    cw = torch.tensor(weights, dtype=torch.float32)
    print(f"[info] class weights: {dict(zip(LOC_CLASSES, weights.round(3)))}")

    model = _LocHead(in_dim=X.shape[1])
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    crit = nn.CrossEntropyLoss(weight=cw)

    ds = TensorDataset(torch.tensor(X_train), torch.tensor(y_train))
    loader = DataLoader(ds, batch_size=batch, shuffle=True)

    best_val_acc = 0.0
    best_state = None
    stale = 0
    Xv = torch.tensor(X_val); yv = torch.tensor(y_val)
    for ep in range(epochs):
        model.train()
        for xb, yb in loader:
            opt.zero_grad(set_to_none=True)
            logits = model(xb)
            loss = crit(logits, yb)
            loss.backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            val_logits = model(Xv)
            val_acc = (val_logits.argmax(1) == yv).float().mean().item()
        if val_acc > best_val_acc:
            best_val_acc = val_acc
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
            if stale >= 5:
                break
        if ep % 3 == 0:
            print(f"  epoch {ep}: val_acc={val_acc:.3f}  best={best_val_acc:.3f}", flush=True)
    if best_state is not None:
        model.load_state_dict(best_state)
    print(f"[info] best val_acc: {best_val_acc:.3f}")
    return model


def _apply_to_all_isoforms(model: _LocHead) -> pd.DataFrame:
    df = pd.read_parquet(ISOFORM_EMB_PARQUET)
    esm_cols = [c for c in df.columns if c.startswith("esm2_")]
    X = df[esm_cols].to_numpy(dtype=np.float32)
    model.eval()
    with torch.no_grad():
        logits = model(torch.tensor(X))
        probs = torch.softmax(logits, dim=1).numpy()
    out = pd.DataFrame(probs, columns=[f"loc_prob_{c}" for c in LOC_CLASSES])
    out["species"] = df["species"].values
    out["gene_name"] = df["gene_name"].values
    out["primary_accession"] = df["primary_accession"].values
    out["isoform_id"] = df["isoform_id"].values
    return out[["species","gene_name","primary_accession","isoform_id"] + [f"loc_prob_{c}" for c in LOC_CLASSES]]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--epochs", type=int, default=30)
    args = ap.parse_args()

    print(f"[info] {time.strftime('%H:%M:%S')} loading labels", flush=True)
    labels = _load_training_labels()
    print(f"[info] {len(labels)} accessions with subcellular_locations", flush=True)

    print(f"[info] {time.strftime('%H:%M:%S')} building training set", flush=True)
    X, y, _ = _build_canonical_embedding_table(labels)

    print(f"[info] {time.strftime('%H:%M:%S')} training", flush=True)
    model = _train(X, y, epochs=args.epochs)
    torch.save(model.state_dict(), OUT_HEAD)
    print(f"[ok] saved trained head -> {OUT_HEAD}")

    print(f"[info] {time.strftime('%H:%M:%S')} applying to all isoforms", flush=True)
    out = _apply_to_all_isoforms(model)
    out.to_csv(OUT_PROBS, sep="\t", index=False)
    print(f"[ok] wrote {OUT_PROBS}  rows={len(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
