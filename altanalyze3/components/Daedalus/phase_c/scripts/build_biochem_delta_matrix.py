#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from peptides import Peptide


ROOT = Path(__file__).resolve().parents[3]
PHASE_A = ROOT / "Daedalus" / "phase_a" / "data"
PHASE_B = ROOT / "Daedalus" / "phase_b" / "data" / "processed"
PHASE_C = ROOT / "Daedalus" / "phase_c" / "data" / "processed"

PAIR_PATH = PHASE_A / "interim" / "priority_pair_subsets.tsv"
BASE_MATRIX = PHASE_B / "reference_delta_matrix.parquet"
UNIPROT_FEATURES = PHASE_A / "interim" / "uniprot_features.tsv"
UNIPROT_MAP = PHASE_A / "interim" / "uniprot_gencode_map.tsv"
GTF_PROTEINS = {
    "human": PHASE_A / "raw" / "gencode" / "gencode.v49.pc_translations.fa.gz",
    "mouse": PHASE_A / "raw" / "gencode" / "gencode.vM38.pc_translations.fa.gz",
}

OUT_TSV = PHASE_C / "biochem_delta_matrix.tsv"
OUT_PARQUET = PHASE_C / "biochem_delta_matrix.parquet"
OUT_SCHEMA = PHASE_C / "biochem_delta_schema.json"

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
ROLE_PATTERNS = {
    "signal": lambda ftype, desc: ftype == "Signal",
    "tm": lambda ftype, desc: ftype in {"Transmembrane", "Intramembrane"},
    "extracellular": lambda ftype, desc: ftype == "Topological domain" and any(x in desc for x in ["extracellular", "luminal", "lumenal", "external"]),
    "cytoplasmic": lambda ftype, desc: ftype == "Topological domain" and any(x in desc for x in ["cytoplasmic", "cytosolic", "intracellular"]),
    "glycosylation": lambda ftype, desc: ftype == "Glycosylation",
    "disulfide": lambda ftype, desc: ftype == "Disulfide bond",
    "active_site": lambda ftype, desc: ftype == "Active site",
    "phospho": lambda ftype, desc: ftype == "Modified residue" and "phospho" in desc,
    "ppi_interface": lambda ftype, desc: ftype in {"Binding site", "Site", "Region", "Domain", "Motif", "Coiled coil"} and any(
        x in desc for x in ["binding", "interact", "dimer", "partner", "sh2", "sh3", "pdz", "ww domain"]
    ),
    "dimerization": lambda ftype, desc: ftype in {"Region", "Domain", "Motif", "Coiled coil"} and any(
        x in desc for x in ["dimer", "oligomer", "coiled-coil", "leucine zipper", "homodimer", "heterodimer"]
    ),
    "dna_interface": lambda ftype, desc: ftype in {"DNA binding", "Zinc finger"} or any(
        x in desc for x in ["dna-binding", "dna binding", "homeobox", "bzip", "helix-loop-helix", "forkhead", "ets", "hmg", "helix-turn-helix"]
    ),
    "nls_region": lambda ftype, desc: any(x in desc for x in ["nuclear localization", "nls"]),
    "nes_region": lambda ftype, desc: any(x in desc for x in ["nuclear export", "nes"]),
    "activation_region": lambda ftype, desc: any(x in desc for x in ["activation", "transactivation"]),
    "repression_region": lambda ftype, desc: any(x in desc for x in ["repression", "corepressor"]),
    "kinase_region": lambda ftype, desc: any(x in desc for x in ["kinase", "catalytic", "atp", "substrate"]),
    "trafficking_region": lambda ftype, desc: any(x in desc for x in ["trafficking", "sorting", "internalization", "basolateral"]),
}
ROLE_NAMES = list(ROLE_PATTERNS)
MOTIF_PATTERNS = {
    "nls_basic": re.compile(r"[KR]{3,}"),
    "leucine_zipper": re.compile(r"L.{6}L.{6}L"),
    "dileucine": re.compile(r"[DEQRNSTAPG].?LL"),
    "yxxphi": re.compile(r"Y..[AILMFWV]"),
    "gxxxg": re.compile(r"G...G"),
    "kdel_tail": re.compile(r"KDEL$"),
    "kkxx_tail": re.compile(r"KK..$"),
}


def _clean_seq(seq: str) -> str:
    seq = (seq or "").upper()
    return "".join(ch for ch in seq if ch in set(AA_ORDER))


def _safe_float(value: float | int | None) -> float:
    if value is None:
        return 0.0
    if isinstance(value, (float, np.floating)):
        if math.isnan(value) or math.isinf(value):
            return 0.0
    return float(value)


def _parse_int(value: object) -> int:
    try:
        if value is None:
            return 0
        if isinstance(value, (float, np.floating)) and (math.isnan(value) or math.isinf(value)):
            return 0
        return int(float(value))
    except Exception:
        return 0


def _parse_fasta_map() -> dict[str, str]:
    seqs: dict[str, str] = {}
    for species, path in GTF_PROTEINS.items():
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                parts = record.id.split("|")
                if len(parts) < 2:
                    continue
                transcript_id = parts[1]
                seqs[transcript_id] = _clean_seq(str(record.seq))
    return seqs


def _load_accession_map() -> dict[str, list[str]]:
    df = pd.read_csv(UNIPROT_MAP, sep="\t")
    tx_to_acc: dict[str, set[str]] = defaultdict(set)
    for row in df.itertuples(index=False):
        tx = getattr(row, "gencode_transcript_id", "")
        acc = getattr(row, "primary_accession", "")
        if tx and acc:
            tx_to_acc[tx].add(acc)
    return {tx: sorted(accs) for tx, accs in tx_to_acc.items()}


def _load_uniprot_feature_index() -> dict[str, list[tuple[str, str, int, int]]]:
    df = pd.read_csv(UNIPROT_FEATURES, sep="\t")
    index: dict[str, list[tuple[str, str, int, int]]] = defaultdict(list)
    for row in df.itertuples(index=False):
        begin = _parse_int(getattr(row, "begin", 0))
        end = _parse_int(getattr(row, "end", 0))
        if begin <= 0 or end <= 0 or end < begin:
            continue
        index[getattr(row, "primary_accession")].append(
            (
                getattr(row, "feature_type", ""),
                str(getattr(row, "description", "") or "").lower(),
                begin,
                end,
            )
        )
    return index


def _segment_change(ref_seq: str, alt_seq: str) -> dict[str, float | int | str]:
    ref_len = len(ref_seq)
    alt_len = len(alt_seq)
    prefix = 0
    max_prefix = min(ref_len, alt_len)
    while prefix < max_prefix and ref_seq[prefix] == alt_seq[prefix]:
        prefix += 1
    suffix = 0
    while (
        suffix < (ref_len - prefix)
        and suffix < (alt_len - prefix)
        and ref_seq[ref_len - 1 - suffix] == alt_seq[alt_len - 1 - suffix]
    ):
        suffix += 1
    ref_changed = ref_seq[prefix : ref_len - suffix if suffix else ref_len]
    alt_changed = alt_seq[prefix : alt_len - suffix if suffix else alt_len]
    if prefix == 0 and suffix > 0:
        change_class = "n_term"
    elif prefix > 0 and suffix == 0:
        change_class = "c_term"
    elif prefix > 0 and suffix > 0:
        change_class = "internal"
    elif ref_changed == alt_changed:
        change_class = "identical"
    else:
        change_class = "complex"
    return {
        "common_prefix_len": prefix,
        "common_suffix_len": suffix,
        "ref_changed_len": len(ref_changed),
        "alt_changed_len": len(alt_changed),
        "ref_changed_frac": len(ref_changed) / ref_len if ref_len else 0.0,
        "alt_changed_frac": len(alt_changed) / alt_len if alt_len else 0.0,
        "change_class": change_class,
        "is_n_term_change": int(change_class == "n_term"),
        "is_c_term_change": int(change_class == "c_term"),
        "is_internal_change": int(change_class == "internal"),
        "is_complex_change": int(change_class == "complex"),
        "is_truncation": int(alt_len < ref_len),
        "is_extension": int(alt_len > ref_len),
        "ref_change_start": prefix + 1 if ref_changed else 0,
        "ref_change_end": ref_len - suffix if ref_changed else 0,
        "ref_changed_segment": ref_changed,
        "alt_changed_segment": alt_changed,
    }


def _pep_descriptor(prefix: str, seq: str) -> dict[str, float]:
    seq = _clean_seq(seq)
    if not seq:
        fields = {
            f"{prefix}_len": 0.0,
            f"{prefix}_aliphatic_index": 0.0,
            f"{prefix}_boman": 0.0,
            f"{prefix}_charge": 0.0,
            f"{prefix}_hydrophobicity": 0.0,
            f"{prefix}_instability_index": 0.0,
            f"{prefix}_isoelectric_point": 0.0,
        }
        for i in range(5):
            fields[f"{prefix}_atchley_{i+1}"] = 0.0
        for i in range(10):
            fields[f"{prefix}_kidera_{i+1}"] = 0.0
        for i in range(8):
            fields[f"{prefix}_vhse_{i+1}"] = 0.0
        return fields
    pep = Peptide(seq)
    fields = {
        f"{prefix}_len": float(len(seq)),
        f"{prefix}_aliphatic_index": _safe_float(pep.aliphatic_index()),
        f"{prefix}_boman": _safe_float(pep.boman()),
        f"{prefix}_charge": _safe_float(pep.charge()),
        f"{prefix}_hydrophobicity": _safe_float(pep.hydrophobicity()),
        f"{prefix}_instability_index": _safe_float(pep.instability_index()),
        f"{prefix}_isoelectric_point": _safe_float(pep.isoelectric_point()),
    }
    for i, value in enumerate(pep.atchley_factors(), start=1):
        fields[f"{prefix}_atchley_{i}"] = _safe_float(value)
    for i, value in enumerate(pep.kidera_factors(), start=1):
        fields[f"{prefix}_kidera_{i}"] = _safe_float(value)
    for i, value in enumerate(pep.vhse_scales(), start=1):
        fields[f"{prefix}_vhse_{i}"] = _safe_float(value)
    return fields


def _motif_counts(prefix: str, seq: str) -> dict[str, float]:
    seq = _clean_seq(seq)
    out: dict[str, float] = {}
    for name, pattern in MOTIF_PATTERNS.items():
        out[f"{prefix}_motif_{name}_count"] = float(len(pattern.findall(seq))) if seq else 0.0
    return out


def _tail_descriptors(prefix: str, seq: str, n: int = 30) -> dict[str, float]:
    seq = _clean_seq(seq)
    n_seq = seq[:n]
    c_seq = seq[-n:] if len(seq) > n else seq
    out = {}
    out.update(_pep_descriptor(f"{prefix}_nterm", n_seq))
    out.update(_pep_descriptor(f"{prefix}_cterm", c_seq))
    return out


def _overlap(a1: int, a2: int, b1: int, b2: int) -> int:
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def _feature_role(ftype: str, desc: str) -> list[str]:
    roles: list[str] = []
    for name, fn in ROLE_PATTERNS.items():
        if fn(ftype, desc):
            roles.append(name)
    return roles


def _typed_overlap(accessions: list[str], feature_index: dict[str, list[tuple[str, str, int, int]]], start: int, end: int, changed_len: int) -> dict[str, float]:
    out = {f"overlap_{role}_count": 0.0 for role in ROLE_NAMES}
    out.update({f"overlap_{role}_span": 0.0 for role in ROLE_NAMES})
    out.update({f"overlap_{role}_frac": 0.0 for role in ROLE_NAMES})
    if not accessions or start <= 0 or end <= 0 or changed_len <= 0:
        return out
    for acc in accessions:
        for ftype, desc, begin, finish in feature_index.get(acc, []):
            ov = _overlap(start, end, begin, finish)
            if ov <= 0:
                continue
            for role in _feature_role(ftype, desc):
                out[f"overlap_{role}_count"] += 1.0
                out[f"overlap_{role}_span"] += float(ov)
    for role in ROLE_NAMES:
        out[f"overlap_{role}_frac"] = out[f"overlap_{role}_span"] / float(changed_len)
    return out


def _delta_features(ref: dict[str, float], alt: dict[str, float], out_prefix: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for key, ref_val in ref.items():
        base = key.removeprefix("ref_").removeprefix("alt_")
        alt_key = key.replace("ref_", "alt_") if key.startswith("ref_") else f"alt_{base}"
        if alt_key in alt:
            out[f"{out_prefix}_{base}_delta"] = _safe_float(alt[alt_key]) - _safe_float(ref_val)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Build Phase C biochemical delta matrix.")
    parser.add_argument("--max-pairs", type=int, default=0, help="Optional cap on eligible pairs for a faster benchmark-aligned run.")
    parser.add_argument("--sample-seed", type=int, default=0, help="Random seed used when --max-pairs is set.")
    args = parser.parse_args()

    PHASE_C.mkdir(parents=True, exist_ok=True)
    base = pd.read_parquet(BASE_MATRIX)
    task_cols = [
        "task_global_seed",
        "task_membrane_seed",
        "task_surface_seed",
        "task_kinase_seed",
        "task_tf_seed",
    ]
    eligible = base[
        base["label_preservation_seed"].notna() & (base[task_cols].sum(axis=1) > 0)
    ].copy()
    if args.max_pairs and len(eligible) > args.max_pairs:
        eligible = eligible.sample(n=args.max_pairs, random_state=args.sample_seed).reset_index(drop=True)
    pair_df = pd.read_csv(PAIR_PATH, sep="\t", low_memory=False)
    seq_map = _parse_fasta_map()
    tx_to_acc = _load_accession_map()
    feature_index = _load_uniprot_feature_index()

    pair_df = pair_df[
        [
            "species",
            "gene_id",
            "gene_name",
            "reference_transcript_id",
            "alternative_transcript_id",
            "family_class",
            "reference_source",
        ]
    ].copy()
    merge_cols = [
        "species",
        "gene_id",
        "gene_name",
        "reference_transcript_id",
        "alternative_transcript_id",
        "family_class",
        "reference_source",
    ]
    pair_df = pair_df.merge(eligible[merge_cols], on=merge_cols, how="inner")
    print(f"[status] eligible_pairs={len(eligible)} merged_pairs={len(pair_df)}", flush=True)
    pair_df["reference_seq"] = pair_df["reference_transcript_id"].map(seq_map).fillna("")
    pair_df["alternative_seq"] = pair_df["alternative_transcript_id"].map(seq_map).fillna("")
    pair_df["reference_accessions"] = pair_df["reference_transcript_id"].map(tx_to_acc).apply(lambda x: x if isinstance(x, list) else [])

    transcript_cache: dict[str, dict[str, float]] = {}
    for transcript_id in pd.unique(
        pd.concat([pair_df["reference_transcript_id"], pair_df["alternative_transcript_id"]], ignore_index=True)
    ):
        seq = seq_map.get(transcript_id, "")
        prefix = f"tx::{transcript_id}"
        cached = {}
        cached.update(_pep_descriptor("full", seq))
        cached.update(_tail_descriptors("seq", seq))
        cached.update(_motif_counts("full", seq))
        transcript_cache[transcript_id] = cached
    print(f"[status] transcript_cache={len(transcript_cache)}", flush=True)

    rows: list[dict[str, float | str]] = []
    for idx, row in enumerate(pair_df.itertuples(index=False), start=1):
        ref_seq = getattr(row, "reference_seq", "") or ""
        alt_seq = getattr(row, "alternative_seq", "") or ""
        seg = _segment_change(ref_seq, alt_seq)
        ref_full = {
            f"ref_{k}": v for k, v in transcript_cache.get(getattr(row, "reference_transcript_id"), {}).items()
        }
        alt_full = {
            f"alt_{k}": v for k, v in transcript_cache.get(getattr(row, "alternative_transcript_id"), {}).items()
        }
        ref_changed = _pep_descriptor("ref_changed", seg["ref_changed_segment"])
        alt_changed = _pep_descriptor("alt_changed", seg["alt_changed_segment"])
        change_motifs = _motif_counts("alt_changed", seg["alt_changed_segment"])
        overlap = _typed_overlap(
            getattr(row, "reference_accessions", []),
            feature_index,
            int(seg["ref_change_start"]),
            int(seg["ref_change_end"]),
            int(seg["ref_changed_len"]),
        )
        feature_row: dict[str, float | str] = {
            "species": getattr(row, "species"),
            "gene_id": getattr(row, "gene_id"),
            "gene_name": getattr(row, "gene_name"),
            "reference_transcript_id": getattr(row, "reference_transcript_id"),
            "alternative_transcript_id": getattr(row, "alternative_transcript_id"),
            "family_class": getattr(row, "family_class"),
            "reference_source": getattr(row, "reference_source"),
            "change_class": str(seg["change_class"]),
            "common_prefix_len": float(seg["common_prefix_len"]),
            "common_suffix_len": float(seg["common_suffix_len"]),
            "ref_changed_len": float(seg["ref_changed_len"]),
            "alt_changed_len": float(seg["alt_changed_len"]),
            "ref_changed_frac": float(seg["ref_changed_frac"]),
            "alt_changed_frac": float(seg["alt_changed_frac"]),
            "is_n_term_change": float(seg["is_n_term_change"]),
            "is_c_term_change": float(seg["is_c_term_change"]),
            "is_internal_change": float(seg["is_internal_change"]),
            "is_complex_change": float(seg["is_complex_change"]),
            "is_truncation": float(seg["is_truncation"]),
            "is_extension": float(seg["is_extension"]),
            "ref_change_start": float(seg["ref_change_start"]),
            "ref_change_end": float(seg["ref_change_end"]),
        }
        feature_row.update(ref_full)
        feature_row.update(alt_full)
        feature_row.update(ref_changed)
        feature_row.update(alt_changed)
        feature_row.update(change_motifs)
        feature_row.update(_delta_features(ref_full, alt_full, "full"))
        feature_row.update(_delta_features(ref_changed, alt_changed, "changed"))
        feature_row.update(overlap)
        rows.append(feature_row)
        if idx % 5000 == 0:
            print(f"[status] encoded_pairs={idx}", flush=True)

    feature_df = pd.DataFrame(rows)
    out = eligible.merge(feature_df, on=merge_cols, how="inner", suffixes=("", "_phasec"))
    out.to_csv(OUT_TSV, sep="\t", index=False)
    out.to_parquet(OUT_PARQUET, index=False)

    excluded = {
        "species",
        "gene_id",
        "gene_name",
        "split",
        "family_class",
        "reference_source",
        "reference_transcript_id",
        "alternative_transcript_id",
        "label_preservation_seed",
        "task_global_seed",
        "task_membrane_seed",
        "task_surface_seed",
        "task_kinase_seed",
        "task_tf_seed",
        "change_class",
    }
    numeric_columns = [c for c in out.columns if c not in excluded]
    schema = {
        "id_columns": ["gene_id", "gene_name", "reference_transcript_id", "alternative_transcript_id"],
        "categorical_columns": ["family_class", "reference_source", "change_class"],
        "feature_columns": numeric_columns,
        "label_column": "label_preservation_seed",
        "task_columns": [
            "task_global_seed",
            "task_membrane_seed",
            "task_surface_seed",
            "task_kinase_seed",
            "task_tf_seed",
        ],
        "description": "Phase C segment-aware biochemical feature matrix built on top of the Phase B reference-delta matrix.",
    }
    with OUT_SCHEMA.open("w", encoding="utf-8") as handle:
        json.dump(schema, handle, indent=2)
    print(f"[ok] Wrote {OUT_TSV} rows={len(out)}")
    print(f"[ok] Wrote {OUT_PARQUET}")
    print(f"[ok] Wrote {OUT_SCHEMA}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
