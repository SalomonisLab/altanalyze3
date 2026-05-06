#!/usr/bin/env python3
"""Mean-pooled ESM2 embeddings (480-dim) for every catalog isoform's
reconstructed protein sequence. The embeddings encode general
protein-language-model knowledge of structural propensity, surface vs
internal residue patterns, and folding-relevant motifs that aren't
captured by the count-statistic feature pipeline.

Inputs:
    data_ingest/data/interim/uniprot_isoform_proteins.tsv

Output:
    data_ingest/data/interim/uniprot_isoform_esm2_t12_embeddings.parquet
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from transformers import AutoModel, AutoTokenizer

INTERIM_DIR = Path(__file__).resolve().parents[1] / "data" / "interim"
PROT_PATH = INTERIM_DIR / "uniprot_isoform_proteins.tsv"
OUT_PATH = INTERIM_DIR / "uniprot_isoform_esm2_t12_embeddings.parquet"

MODEL_NAME = "facebook/esm2_t12_35M_UR50D"
EMB_DIM = 480
MAX_LEN = 1022  # ESM2 max with [CLS] + [SEP] = 1024


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-name", default=MODEL_NAME)
    ap.add_argument("--out", type=Path, default=OUT_PATH)
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--max-len", type=int, default=MAX_LEN)
    args = ap.parse_args()

    proteins = pd.read_csv(PROT_PATH, sep="\t")
    print(f"[info] {len(proteins)} isoforms", flush=True)

    print(f"[info] loading {args.model_name} ...", flush=True)
    t0 = time.time()
    tok = AutoTokenizer.from_pretrained(args.model_name)
    model = AutoModel.from_pretrained(args.model_name)
    model.eval()
    device = torch.device("cpu")
    model.to(device)
    print(f"[info] model loaded in {time.time()-t0:.1f}s, hidden_dim={model.config.hidden_size}", flush=True)

    embeddings = np.zeros((len(proteins), model.config.hidden_size), dtype=np.float32)
    seq_lens = np.zeros(len(proteins), dtype=np.int32)

    # Sort indices by sequence length to amortise padding cost within batch
    sequences = proteins["protein_sequence"].astype(str).tolist()
    order = sorted(range(len(sequences)), key=lambda i: len(sequences[i]))

    t_start = time.time()
    for batch_start in range(0, len(order), args.batch_size):
        batch_idx = order[batch_start: batch_start + args.batch_size]
        batch_seqs = [sequences[i][: args.max_len] for i in batch_idx]
        with torch.no_grad():
            inputs = tok(batch_seqs, return_tensors="pt", padding=True, truncation=True, max_length=args.max_len + 2)
            inputs = {k: v.to(device) for k, v in inputs.items()}
            out = model(**inputs)
            hidden = out.last_hidden_state  # (B, L, H)
            attn_mask = inputs["attention_mask"].unsqueeze(-1).float()
            # Skip [CLS] (pos 0) and [SEP] (last). Use mask to avoid padding.
            # Strategy: mean-pool over residues with attention=1 excluding the
            # special tokens. Simpler: use attention mask with first/last set to 0.
            mask = attn_mask.clone()
            mask[:, 0] = 0  # CLS
            for b, idx in enumerate(batch_idx):
                # Find last non-pad position; that's [SEP]
                last = int(inputs["attention_mask"][b].sum().item()) - 1
                if last >= 0:
                    mask[b, last] = 0
            num = (hidden * mask).sum(dim=1)
            den = mask.sum(dim=1).clamp_min(1)
            mean_emb = (num / den).cpu().numpy()
        for j, idx in enumerate(batch_idx):
            embeddings[idx] = mean_emb[j]
            seq_lens[idx] = len(sequences[idx])

        done = batch_start + len(batch_idx)
        if (batch_start // args.batch_size) % 25 == 0 or done == len(order):
            rate = done / max(time.time() - t_start, 1e-6)
            remain = (len(order) - done) / max(rate, 1e-6)
            print(f"[progress] {done}/{len(order)}  rate={rate:.1f}/s  ETA={remain/60:.1f} min", flush=True)

    # Write
    cols = [f"esm2_{i:03d}" for i in range(model.config.hidden_size)]
    df = pd.DataFrame(embeddings, columns=cols)
    df["species"] = proteins["species"].values
    df["gene_name"] = proteins["gene_name"].values
    df["primary_accession"] = proteins["primary_accession"].values
    df["isoform_id"] = proteins["isoform_id"].values
    df["protein_length"] = seq_lens
    df = df[["species", "gene_name", "primary_accession", "isoform_id", "protein_length"] + cols]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    print(f"[ok] wrote {args.out}  rows={len(df)}  cols={len(df.columns)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
