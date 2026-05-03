#!/usr/bin/env python3
"""ESM2 t30 mean-pooled embeddings for every reviewed UniProt entry's
canonical sequence in human + mouse. Used as the training set for the
DeepLoc-style localization classifier.

Inputs:
    phase_a/data/raw/uniprot_reviewed_human_mouse.jsonl.gz

Output:
    phase_a/data/interim/uniprot_canonical_esm2_t30_embeddings.parquet
"""
from __future__ import annotations

import argparse
import gzip
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from transformers import AutoModel, AutoTokenizer

PHASE_A = Path(__file__).resolve().parents[1] / "data"
RAW_JSONL = PHASE_A / "raw" / "uniprot_reviewed_human_mouse.jsonl.gz"
OUT_PATH = PHASE_A / "interim" / "uniprot_canonical_esm2_t30_embeddings.parquet"

MODEL_NAME = "facebook/esm2_t30_150M_UR50D"
MAX_LEN = 1022


def _load_canonical_sequences():
    rows = []
    with gzip.open(RAW_JSONL, "rt") as f:
        for line in f:
            d = json.loads(line)
            acc = d.get("primaryAccession")
            organism = (d.get("organism") or {}).get("scientificName", "")
            seq = (d.get("sequence") or {}).get("value", "")
            gene = ""
            for g in d.get("genes") or []:
                gn = (g.get("geneName") or {}).get("value")
                if gn:
                    gene = gn
                    break
            primary_loc = ""
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
            if seq:
                rows.append({
                    "primary_accession": acc, "organism": organism, "gene_name": gene,
                    "primary_loc": primary_loc, "sequence": seq, "length": len(seq),
                })
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=OUT_PATH)
    ap.add_argument("--batch-size", type=int, default=4)
    args = ap.parse_args()

    df = _load_canonical_sequences()
    print(f"[info] {len(df)} canonical sequences loaded", flush=True)

    tok = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = AutoModel.from_pretrained(MODEL_NAME)
    model.eval()
    print(f"[info] model loaded, hidden={model.config.hidden_size}", flush=True)

    H = model.config.hidden_size
    embeddings = np.zeros((len(df), H), dtype=np.float32)
    sequences = df["sequence"].tolist()
    order = sorted(range(len(sequences)), key=lambda i: len(sequences[i]))

    t_start = time.time()
    for batch_start in range(0, len(order), args.batch_size):
        idxs = order[batch_start: batch_start + args.batch_size]
        seqs = [sequences[i][:MAX_LEN] for i in idxs]
        with torch.no_grad():
            inputs = tok(seqs, return_tensors="pt", padding=True, truncation=True, max_length=MAX_LEN + 2)
            out = model(**inputs)
            hidden = out.last_hidden_state
            mask = inputs["attention_mask"].unsqueeze(-1).float()
            mask[:, 0] = 0
            for b in range(len(idxs)):
                last = int(inputs["attention_mask"][b].sum().item()) - 1
                if last >= 0:
                    mask[b, last] = 0
            num = (hidden * mask).sum(dim=1)
            den = mask.sum(dim=1).clamp_min(1)
            emb = (num / den).cpu().numpy()
        for j, i in enumerate(idxs):
            embeddings[i] = emb[j]
        done = batch_start + len(idxs)
        if (batch_start // args.batch_size) % 100 == 0 or done == len(order):
            rate = done / max(time.time() - t_start, 1e-6)
            remain = (len(order) - done) / max(rate, 1e-6)
            print(f"[progress] {done}/{len(order)}  rate={rate:.1f}/s  ETA={remain/60:.1f} min", flush=True)

    cols = [f"esm2_{i:03d}" for i in range(H)]
    out_df = pd.DataFrame(embeddings, columns=cols)
    for c in ["primary_accession", "organism", "gene_name", "primary_loc", "length"]:
        out_df[c] = df[c].values
    out_df = out_df[["primary_accession","organism","gene_name","primary_loc","length"] + cols]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(args.out, index=False)
    print(f"[ok] wrote {args.out}  rows={len(out_df)}  dim={H}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
