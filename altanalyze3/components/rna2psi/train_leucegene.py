"""Train the Leucegene AML rna2psi bundle.

Targets: AltAnalyze MultiPath-PSI events that (a) appear in a variant/covariate PSI.* file
with |dPSI| > 0.2 and (b) are present in the 75p PSI EventAnnotation matrix.
Features: protein-coding Ensembl genes from the Leucegene Kallisto **counts** matrix.

Usage:
    python -m components.rna2psi.train_leucegene --out-dir <dir> [--max-events N]
"""
from __future__ import annotations

import argparse
import glob
import gzip
import os
import pickle
import time

import numpy as np
import pandas as pd

from ._train import train_rna2psi

BASE = "/Users/saljh8/Desktop/dataAnalysis/SalomonisLab/Leucegene/hg38-Leucegene"
PSI75 = f"{BASE}/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-75p.txt"
COUNTS = f"{BASE}/counts.Leucegene.txt"
VARIANT_DIRS = [f"{BASE}/MajorCovariates/Events-dPSI_0.1_adjp",
                f"{BASE}/Variants/Events-dPSI_0.1_adjp"]
BIOTYPES = "/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart72/ensembl/Hs/Hs_Ensembl_transcript-biotypes.txt"
SYMBOLS = "/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart72/ensembl/Hs/Ensembl-Symbol.txt"


def candidate_uids(dpsi_cut: float = 0.2) -> set:
    uids = set()
    for d in VARIANT_DIRS:
        for f in sorted(glob.glob(os.path.join(d, "PSI.*.txt"))):
            df = pd.read_csv(f, sep="\t", usecols=lambda c: c in ("UID", "dPSI"))
            if "dPSI" not in df.columns:
                continue
            keep = df["UID"][df["dPSI"].abs() > dpsi_cut]
            uids.update(keep.dropna().astype(str))
    return uids


def protein_coding_genes() -> set:
    bt = pd.read_csv(BIOTYPES, sep="\t")
    col = bt.columns
    return set(bt.loc[bt[col[2]] == "protein_coding", col[0]].astype(str))


def symbol_map() -> dict:
    sm = pd.read_csv(SYMBOLS, sep="\t")
    return dict(zip(sm.iloc[:, 0].astype(str), sm.iloc[:, 1].astype(str)))


def _strip_sample(name: str) -> str:
    import re
    return re.sub(r"_[0-9]+$", "", re.sub(r"\.bed$", "", str(name)))   # SRS..._1.bed -> SRS...


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--max-events", type=int, default=0, help="0 = all; else random subset (smoke test)")
    ap.add_argument("--n-candidates", type=int, default=15)
    ap.add_argument("--max-features", type=int, default=5)
    ap.add_argument("--cv", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    t0 = time.time()

    print("[load] candidate variant events (|dPSI|>0.2) ...", flush=True)
    cands = candidate_uids()
    print(f"       {len(cands)} unique candidate UIDs")

    print("[load] PSI 75p matrix ...", flush=True)
    psi = pd.read_csv(PSI75, sep="\t", index_col=0)          # NA + empty -> NaN by default
    psi = psi[~psi.index.duplicated(keep="first")]
    psi = psi.loc[psi.index.intersection(cands)]
    psi.columns = [_strip_sample(c) for c in psi.columns]
    print(f"       {psi.shape[0]} target events x {psi.shape[1]} samples in matrix")

    print("[load] counts matrix ...", flush=True)
    counts = pd.read_csv(COUNTS, sep="\t", index_col=0)
    counts = counts[~counts.index.duplicated(keep="first")]

    shared = [s for s in psi.columns if s in set(counts.columns)]
    psi = psi[shared]; counts = counts[shared]
    print(f"       {len(shared)} shared AML samples")

    pc = protein_coding_genes()
    feature_genes = [g for g in counts.index if g in pc]
    print(f"       {len(feature_genes)} protein-coding feature genes present in counts")

    if args.max_events and args.max_events < psi.shape[0]:
        rng = np.random.default_rng(args.seed)
        pick = rng.choice(psi.shape[0], size=args.max_events, replace=False)
        psi = psi.iloc[np.sort(pick)]
        print(f"[smoke] restricted to {psi.shape[0]} random events")

    bundle, perf = train_rna2psi(
        counts, psi, feature_genes=feature_genes, symbol_map=symbol_map(),
        n_candidates=args.n_candidates, max_features=args.max_features,
        cv=args.cv, random_state=args.seed)

    tag = f"_smoke{args.max_events}" if args.max_events else ""
    bundle_path = os.path.join(args.out_dir, f"rna2psi_leucegene_bundle{tag}.pkl.gz")
    perf_path = os.path.join(args.out_dir, f"rna2psi_per_event_performance{tag}.csv")
    with gzip.open(bundle_path, "wb") as fh:
        pickle.dump(bundle, fh, protocol=4)
    perf.to_csv(perf_path, index=False)

    v = perf[perf["n_features"] > 0]
    print("\n==== SUMMARY ====")
    print(f"events (targets)              : {len(perf)}")
    print(f"events with a fitted model    : {len(v)}")
    print(f"median held-out CV Spearman   : {np.nanmedian(v['cv_spearman']):.3f}")
    print(f"imputable  (Sp>0.3)           : {int((v['cv_spearman'] > 0.3).sum())}")
    print(f"strong     (Sp>0.5)           : {int((v['cv_spearman'] > 0.5).sum())}")
    print(f"estimator mix                 : {v['estimator'].value_counts().to_dict()}")
    print(f"bundle -> {bundle_path}")
    print(f"perf   -> {perf_path}")
    print(f"elapsed: {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
