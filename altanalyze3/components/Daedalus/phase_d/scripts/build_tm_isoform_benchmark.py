#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


PHASE_D = Path(__file__).resolve().parents[1]
PROCESSED = PHASE_D / "data" / "processed"
PHASE_A = Path(__file__).resolve().parents[2] / "phase_a" / "data" / "interim"

BASE_BENCH_PATH = PROCESSED / "phase_d_benchmark_instances.parquet"
BASE_SCHEMA_PATH = PROCESSED / "phase_d_benchmark_schema.json"
CATALOG_PATH = PHASE_A / "tm_negatives" / "uniprot_isoform_catalog.tsv"
CATALOG_FULL_PATH = PHASE_A / "tm_negatives" / "uniprot_isoform_catalog_full.tsv"
PAIR_PATH = PHASE_A / "isoform_pair_candidates.with_splits.tsv"
MAP_PATH = PHASE_A / "uniprot_gencode_map.tsv"
TRANSCRIPT_CATALOG_PATH = PHASE_A / "transcript_supervision_catalog.tsv"
ISOFORM_RESOLVED_PATH = PHASE_A / "tm_negatives" / "uniprot_isoform_to_ensp_resolved.tsv"

OUT_TSV = PROCESSED / "tm_isoform_benchmark_instances.tsv"
OUT_PARQUET = PROCESSED / "tm_isoform_benchmark_instances.parquet"
OUT_SCHEMA = PROCESSED / "tm_isoform_benchmark_schema.json"
OUT_SUMMARY = PROCESSED / "tm_isoform_benchmark_summary.tsv"

TASK_NAME = "cell_surface_localization"
NEG_GROUP_TM_RETAINED = "tm_retained_non_surface_true_negative"
NEG_GROUP_TM_LOST = "tm_lost_no_tm_control"
POS_GROUP = "uniprot_tm_positive"

# Physicochemical descriptor prefixes that ablation showed contribute
# negligible discriminative value; dropped from the feature matrix.
_PHYSCHEM_PREFIXES = (
    "ref_full_atchley", "ref_full_kidera", "ref_full_vhse",
    "ref_full_aliphatic", "ref_full_boman", "ref_full_charge",
    "ref_full_hydrophobicity", "ref_full_instability", "ref_full_isoelectric",
    "ref_seq_nterm_atchley", "ref_seq_nterm_kidera", "ref_seq_nterm_vhse",
    "ref_seq_nterm_aliphatic", "ref_seq_nterm_boman", "ref_seq_nterm_charge",
    "ref_seq_nterm_hydrophobicity", "ref_seq_nterm_instability", "ref_seq_nterm_isoelectric",
    "ref_seq_cterm_atchley", "ref_seq_cterm_kidera", "ref_seq_cterm_vhse",
    "ref_seq_cterm_aliphatic", "ref_seq_cterm_boman", "ref_seq_cterm_charge",
    "ref_seq_cterm_hydrophobicity", "ref_seq_cterm_instability", "ref_seq_cterm_isoelectric",
    "alt_full_atchley", "alt_full_kidera", "alt_full_vhse",
    "alt_full_aliphatic", "alt_full_boman", "alt_full_charge",
    "alt_full_hydrophobicity", "alt_full_instability", "alt_full_isoelectric",
    "alt_seq_nterm_atchley", "alt_seq_nterm_kidera", "alt_seq_nterm_vhse",
    "alt_seq_nterm_aliphatic", "alt_seq_nterm_boman", "alt_seq_nterm_charge",
    "alt_seq_nterm_hydrophobicity", "alt_seq_nterm_instability", "alt_seq_nterm_isoelectric",
    "alt_seq_cterm_atchley", "alt_seq_cterm_kidera", "alt_seq_cterm_vhse",
    "alt_seq_cterm_aliphatic", "alt_seq_cterm_boman", "alt_seq_cterm_charge",
    "alt_seq_cterm_hydrophobicity", "alt_seq_cterm_instability", "alt_seq_cterm_isoelectric",
    "full_full_atchley", "full_full_kidera", "full_full_vhse",
    "full_full_aliphatic", "full_full_boman", "full_full_charge",
    "full_full_hydrophobicity", "full_full_instability", "full_full_isoelectric",
    "full_seq_nterm_atchley", "full_seq_nterm_kidera", "full_seq_nterm_vhse",
    "full_seq_nterm_aliphatic", "full_seq_nterm_boman", "full_seq_nterm_charge",
    "full_seq_nterm_hydrophobicity", "full_seq_nterm_instability", "full_seq_nterm_isoelectric",
    "full_seq_cterm_atchley", "full_seq_cterm_kidera", "full_seq_cterm_vhse",
    "full_seq_cterm_aliphatic", "full_seq_cterm_boman", "full_seq_cterm_charge",
    "full_seq_cterm_hydrophobicity", "full_seq_cterm_instability", "full_seq_cterm_isoelectric",
    "changed_changed_atchley", "changed_changed_kidera", "changed_changed_vhse",
    "changed_changed_aliphatic", "changed_changed_boman", "changed_changed_charge",
    "changed_changed_hydrophobicity", "changed_changed_instability", "changed_changed_isoelectric",
    "ref_changed_atchley", "ref_changed_kidera", "ref_changed_vhse",
    "ref_changed_aliphatic", "ref_changed_boman", "ref_changed_charge",
    "ref_changed_hydrophobicity", "ref_changed_instability", "ref_changed_isoelectric",
    "alt_changed_atchley", "alt_changed_kidera", "alt_changed_vhse",
    "alt_changed_aliphatic", "alt_changed_boman", "alt_changed_charge",
    "alt_changed_hydrophobicity", "alt_changed_instability", "alt_changed_isoelectric",
)


def _norm_species(value: object) -> str:
    value = str(value or "").strip().lower()
    if value in {"homo sapiens", "human", "hs"}:
        return "human"
    if value in {"mus musculus", "mouse", "mm"}:
        return "mouse"
    return value


def _pair_key(frame: pd.DataFrame) -> pd.Series:
    return frame.apply(
        lambda row: (
            str(row["species"]),
            str(row["gene_name"]),
            str(row["reference_transcript_id"]),
            str(row["alternative_transcript_id"]),
        ),
        axis=1,
    )


def _parse_shared_accessions(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    return {token for token in str(value).split(";") if token}


def _load_full_catalog() -> pd.DataFrame:
    """Full per-isoform catalog (positives + negatives in TM+CM genes)."""
    cat = pd.read_csv(CATALOG_FULL_PATH, sep="\t")
    cat["species"] = cat["organism"].map(_norm_species)
    cat["primary_accession"] = cat["isoform_id"].astype(str).map(lambda x: x.split("-")[0])
    return cat


def _acc_to_transcript_map() -> dict[tuple[str, str], list[str]]:
    mapping = pd.read_csv(MAP_PATH, sep="\t")
    mapping["species"] = mapping["gencode_species"].map(_norm_species)
    return (
        mapping.dropna(subset=["primary_accession", "gencode_transcript_id"])
        .groupby(["primary_accession", "species"])["gencode_transcript_id"]
        .apply(lambda s: sorted({str(x) for x in s if str(x)}))
        .to_dict()
    )


def _load_isoform_resolution() -> dict[tuple[str, str, str], dict[str, str | int | float]]:
    """Load the alignment-resolved (species, gene, isoform_id) -> ENSP/ENST map.

    Returns a dict keyed by (species, gene_name, isoform_id) with values
    {"ensp_id": ..., "enst_id": ..., "identity_pct": ...}.
    """
    if not ISOFORM_RESOLVED_PATH.exists():
        return {}
    df = pd.read_csv(ISOFORM_RESOLVED_PATH, sep="\t")
    out: dict[tuple[str, str, str], dict] = {}
    for r in df.itertuples(index=False):
        key = (str(r.species), str(r.gene_name), str(r.isoform_id))
        out[key] = {
            "ensp_id": str(r.ensp_id),
            "enst_id": str(r.enst_id),
            "identity_pct": float(r.identity_pct),
            "assignment_method": str(r.assignment_method),
            "isoform_protein_length": int(r.isoform_protein_length),
        }
    return out


def _transcript_length_lookup() -> dict[tuple[str, str], int]:
    """(species, gencode_transcript_id) -> protein_length.

    Used to resolve a UniProt isoform suffix (e.g. O14788-2) to the gencode
    transcript whose translated protein length matches the isoform's
    catalog `protein_length`.
    """
    ts = pd.read_csv(TRANSCRIPT_CATALOG_PATH, sep="\t", low_memory=False)
    ts["species"] = ts["species"].map(_norm_species)
    ts["protein_length"] = pd.to_numeric(ts["protein_length"], errors="coerce")
    ts = ts.dropna(subset=["protein_length"])
    return {
        (str(r.species), str(r.transcript_id)): int(r.protein_length)
        for r in ts.itertuples(index=False)
    }


def _resolve_isoform_to_transcript(
    iso_row,
    candidate_txs: list[str],
    tx_lengths: dict[tuple[str, str], int],
) -> str | None:
    """Pick the gencode transcript whose protein length matches the catalog
    isoform's `protein_length`. Falls back to the closest length within 5%
    of the target. Returns None when no candidate exists."""
    if not candidate_txs:
        return None
    target = int(iso_row.protein_length) if pd.notna(iso_row.protein_length) else None
    species = str(iso_row.species)
    if target is None:
        return candidate_txs[0]
    exact = [tx for tx in candidate_txs if tx_lengths.get((species, tx)) == target]
    if exact:
        return sorted(exact)[0]
    scored = [
        (abs(tx_lengths.get((species, tx), 0) - target), tx)
        for tx in candidate_txs
        if (species, tx) in tx_lengths
    ]
    if not scored:
        return None
    scored.sort()
    best_diff, best_tx = scored[0]
    if target > 0 and best_diff / target <= 0.05:
        return best_tx
    return None


def _map_isoforms_to_pairs(
    isoforms: pd.DataFrame,
    pairs: pd.DataFrame,
    acc_to_txs: dict[tuple[str, str], list[str]],
    tx_lengths: dict[tuple[str, str], int],
    role: str,
) -> pd.DataFrame:
    """Map catalog isoforms to alternative-transcript pair rows.

    Each catalog isoform suffix (e.g. O14788-2 = 244 aa) is resolved to the
    single gencode transcript whose protein length matches; that transcript
    becomes the pair's alternative. Pairs are kept only when the pair's
    shared_uniprot_accessions also include the isoform's primary accession
    (so the reference and alternative both map to the same UniProt entry).
    """
    prefix = role  # "positive" | "negative"
    rows: list[dict[str, object]] = []
    for iso in isoforms.itertuples(index=False):
        candidate_txs = list(acc_to_txs.get((str(iso.primary_accession), str(iso.species)), []))
        if not candidate_txs:
            continue
        resolved_tx = _resolve_isoform_to_transcript(iso, candidate_txs, tx_lengths)
        if resolved_tx is None:
            continue
        sub = pairs[
            (pairs["species"] == str(iso.species))
            & (pairs["gene_name"].astype(str) == str(iso.gene_name))
            & (pairs["alternative_transcript_id"].astype(str) == resolved_tx)
        ]
        if sub.empty:
            continue
        sub = sub[
            sub["shared_uniprot_accessions"].map(
                lambda value: str(iso.primary_accession) in _parse_shared_accessions(value)
            )
        ]
        if sub.empty:
            continue
        for pair in sub.itertuples(index=False):
            rows.append({
                "species": str(iso.species),
                "gene_name": str(iso.gene_name),
                "reference_transcript_id": str(pair.reference_transcript_id),
                "alternative_transcript_id": str(pair.alternative_transcript_id),
                f"{prefix}_isoform_id": str(iso.isoform_id),
                f"{prefix}_isoform_locations": str(iso.isoform_locations or ""),
                f"{prefix}_tm_count": int(iso.tm_count),
                f"{prefix}_signal_peptide": str(iso.signal_peptide or ""),
                f"{prefix}_alt_seq_changes": str(iso.alt_seq_changes or ""),
                f"{prefix}_protein_length": int(iso.protein_length) if pd.notna(iso.protein_length) else 0,
                "primary_accession": str(iso.primary_accession),
            })
    return pd.DataFrame(rows)


def _aggregate_pairs(mapped: pd.DataFrame, role: str) -> pd.DataFrame:
    if mapped.empty:
        return mapped
    mapped = mapped.copy()
    mapped["pair_key"] = _pair_key(mapped)
    prefix = role
    return (
        mapped.groupby("pair_key", as_index=False)
        .agg(
            species=("species", "first"),
            gene_name=("gene_name", "first"),
            reference_transcript_id=("reference_transcript_id", "first"),
            alternative_transcript_id=("alternative_transcript_id", "first"),
            primary_accessions=("primary_accession", lambda s: ";".join(sorted(set(s)))),
            **{
                f"{prefix}_isoform_ids": (f"{prefix}_isoform_id", lambda s: ";".join(sorted(set(s)))),
                f"{prefix}_isoform_locations": (f"{prefix}_isoform_locations", lambda s: ";".join(sorted({x for x in s if x}))),
                f"{prefix}_tm_counts": (f"{prefix}_tm_count", lambda s: ";".join(str(x) for x in sorted(set(int(v) for v in s)))),
                f"{prefix}_signal_peptides": (f"{prefix}_signal_peptide", lambda s: ";".join(sorted({x for x in s if x}))),
                f"{prefix}_alt_seq_changes": (f"{prefix}_alt_seq_changes", lambda s: " || ".join(sorted({x for x in s if x}))),
                f"{prefix}_protein_lengths": (f"{prefix}_protein_length", lambda s: ";".join(str(x) for x in sorted(set(int(v) for v in s)))),
                "supporting_catalog_isoform_count": (f"{prefix}_isoform_id", pd.Series.nunique),
            },
        )
    )


def _build_alignment_pairs(
    catalog: pd.DataFrame,
    resolution: dict[tuple[str, str, str], dict],
    base_task: pd.DataFrame,
    role: str,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Build (ref_ENST, alt_ENST) pairs using the alignment-resolved
    isoform_id -> ENST map. The alternative ENST comes from the catalog
    isoform's own resolution; the reference ENST comes from the catalog's
    `canonical_isoform_id` resolution. A pair is emitted only when the
    phase_d feature matrix has a row at that exact (ref_ENST, alt_ENST)."""
    base_keys = set(base_task["pair_key"].tolist()) if "pair_key" in base_task.columns else set()
    rows: list[dict[str, object]] = []
    n_no_alt_resolution = 0
    n_no_ref_resolution = 0
    n_same_tx = 0
    n_no_features = 0

    for iso in catalog.itertuples(index=False):
        species = str(iso.species)
        gene = str(iso.gene_name)
        iso_id = str(iso.isoform_id)
        canon_iso_id = str(iso.canonical_isoform_id)

        alt = resolution.get((species, gene, iso_id))
        if alt is None:
            n_no_alt_resolution += 1
            continue
        ref = resolution.get((species, gene, canon_iso_id))
        if ref is None:
            n_no_ref_resolution += 1
            continue
        alt_tx = alt["enst_id"]
        ref_tx = ref["enst_id"]
        if alt_tx == ref_tx:
            n_same_tx += 1
            continue
        pair_key = (species, gene, ref_tx, alt_tx)
        if pair_key not in base_keys:
            n_no_features += 1
            continue
        rows.append({
            "species": species,
            "gene_name": gene,
            "reference_transcript_id": ref_tx,
            "alternative_transcript_id": alt_tx,
            f"{role}_isoform_id": iso_id,
            f"{role}_isoform_locations": str(iso.isoform_locations or ""),
            f"{role}_tm_count": int(iso.tm_count) if pd.notna(iso.tm_count) else 0,
            f"{role}_signal_peptide": str(iso.signal_peptide or ""),
            f"{role}_alt_seq_changes": str(iso.alt_seq_changes or ""),
            f"{role}_protein_length": int(iso.protein_length) if pd.notna(iso.protein_length) else 0,
            f"{role}_identity_pct": float(alt["identity_pct"]),
            "primary_accession": str(iso.primary_accession),
        })

    mapped = pd.DataFrame(rows)
    stats = {
        "catalog_rows": int(len(catalog)),
        "mapped_rows": int(len(mapped)),
        "no_alt_resolution": int(n_no_alt_resolution),
        "no_ref_resolution": int(n_no_ref_resolution),
        "ref_eq_alt": int(n_same_tx),
        "no_feature_row": int(n_no_features),
    }
    if mapped.empty:
        return mapped, stats
    agg = _aggregate_pairs(mapped, role=role)
    joined = base_task[["pair_key"]].drop_duplicates().merge(agg, on="pair_key", how="inner")
    stats["mapped_unique_pairs"] = int(len(agg))
    stats["feature_joined_unique_pairs"] = int(len(joined))
    return joined, stats


def _build_catalog_tm_negative_pairs(base_task: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    """Isoform-level true negatives: TM retained, Cell Membrane lost (per
    UniProt isoform catalog), with ENST resolved by sequence alignment."""
    catalog = _load_full_catalog()
    catalog = catalog[
        (catalog["canonical_has_cell_membrane"] == 1)
        & (catalog["canonical_has_transmembrane"] == 1)
        & (catalog["is_canonical"] == 0)
        & (catalog["has_cell_membrane"] == 0)
        & (catalog["has_transmembrane"] == 1)
    ].copy()
    catalog["primary_accession"] = catalog["isoform_id"].astype(str).str.split("-").str[0]
    resolution = _load_isoform_resolution()
    return _build_alignment_pairs(catalog, resolution, base_task, role="negative")


def _build_catalog_tm_positive_pairs(
    base_task: pd.DataFrame, excluded_keys: set,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Isoform-level positives: every catalog isoform with TM + Cell Membrane,
    paired against the catalog's canonical_isoform_id (the gene's reference
    surface protein). ENST resolution comes from the alignment table."""
    catalog = _load_full_catalog()
    catalog = catalog[
        (catalog["canonical_has_cell_membrane"] == 1)
        & (catalog["canonical_has_transmembrane"] == 1)
        & (catalog["has_cell_membrane"] == 1)
        & (catalog["has_transmembrane"] == 1)
    ].copy()
    catalog["primary_accession"] = catalog["isoform_id"].astype(str).str.split("-").str[0]
    resolution = _load_isoform_resolution()
    joined, stats = _build_alignment_pairs(catalog, resolution, base_task, role="positive")
    if not joined.empty and excluded_keys:
        joined = joined[~joined["pair_key"].isin(excluded_keys)].copy()
        stats["after_excluded_keys"] = int(len(joined))
    return joined, stats


def _sample_controls(
    base_task: pd.DataFrame,
    excluded_keys: set[tuple[str, str, str, str]],
    split_caps: dict[str, int] | None,
    random_seed: int,
) -> pd.DataFrame:
    candidates = base_task[
        (base_task["reference_tm_feature_count"] > 0)
        & (base_task["reference_is_membrane"] > 0)
        & (base_task["alternative_predicted_tm_count"] == 0)
        & (base_task["overlap_tm_count"] == 0)
        & (base_task["alt_has_protein"] > 0)
    ].copy()
    candidates = candidates[~candidates["pair_key"].isin(excluded_keys)].copy()

    if split_caps is None:
        return candidates.copy().reset_index(drop=True)

    sampled: list[pd.DataFrame] = []
    for split, cap in split_caps.items():
        if cap <= 0:
            continue
        sub = candidates[candidates["split"] == split].copy()
        if sub.empty:
            continue
        n_take = min(int(cap), len(sub))
        sampled.append(sub.sample(n=n_take, random_state=random_seed))
    if not sampled:
        return candidates.iloc[0:0].copy()
    return pd.concat(sampled, axis=0, ignore_index=True)


def _cap_positives(
    positives: pd.DataFrame,
    neg_ratio: float,
    n_negatives_by_split: dict[str, int],
    random_seed: int,
) -> pd.DataFrame:
    """Optionally sub-sample positives per split at a fixed positives:negatives ratio."""
    if float(neg_ratio) <= 0 or not n_negatives_by_split:
        return positives.reset_index(drop=True)
    sampled: list[pd.DataFrame] = []
    for split, n_neg in n_negatives_by_split.items():
        cap = int(np.ceil(float(neg_ratio) * n_neg))
        sub = positives[positives["split"] == split].copy()
        if sub.empty:
            continue
        n_take = min(cap, len(sub))
        sampled.append(sub.sample(n=n_take, random_state=random_seed))
    if not sampled:
        return positives.iloc[0:0].copy()
    return pd.concat(sampled, axis=0, ignore_index=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build TM isoform benchmark with strict UniProt positives and curated true negatives."
    )
    parser.add_argument("--task", default=TASK_NAME)
    parser.add_argument(
        "--control-multiplier",
        type=float,
        default=0.0,
        help="No-TM controls per retained-TM negative per split. <=0 keeps all.",
    )
    parser.add_argument(
        "--neg-ratio",
        type=float,
        default=0.0,
        help="Positives per retained-TM negative per split (e.g. 1=1:1, 2=2:1, 4=4:1). <=0 keeps all strict positives.",
    )
    parser.add_argument("--random-seed", type=int, default=0)
    args = parser.parse_args()

    with BASE_SCHEMA_PATH.open() as handle:
        base_schema = json.load(handle)

    # Prune ablation-irrelevant physicochemical features
    all_feature_cols = base_schema["feature_columns"]
    feature_columns = [
        c for c in all_feature_cols
        if not any(c.startswith(p) for p in _PHYSCHEM_PREFIXES)
    ]
    print(f"[info] Features: {len(all_feature_cols)} total -> {len(feature_columns)} after physchem pruning")

    base = pd.read_parquet(BASE_BENCH_PATH)
    base_task = base[base["task"] == args.task].copy()
    base_task["pair_key"] = _pair_key(base_task)
    if base_task["pair_key"].duplicated().any():
        dupes = int(base_task["pair_key"].duplicated().sum())
        raise ValueError(f"{args.task} has {dupes} duplicated pair rows in {BASE_BENCH_PATH}")

    _NEG_META_COLS = (
        "negative_isoform_ids", "negative_isoform_locations", "negative_signal_peptides",
        "negative_alt_seq_changes", "negative_tm_counts",
    )
    _POS_META_COLS = (
        "positive_isoform_ids", "positive_isoform_locations", "positive_signal_peptides",
        "positive_alt_seq_changes", "positive_tm_counts",
    )

    def _fill_missing_meta(frame: pd.DataFrame, cols: tuple[str, ...], fill: str = "") -> None:
        for col in cols:
            if col not in frame.columns:
                frame[col] = fill

    # --- True negatives: per-isoform TM retained, Cell Membrane lost ---
    tm_retained, neg_mapping_stats = _build_catalog_tm_negative_pairs(base_task)
    tm_retained_keys = set(tm_retained["pair_key"].tolist()) if not tm_retained.empty else set()
    tm_retained_rows = base_task.merge(
        tm_retained.drop(columns=["species", "gene_name", "reference_transcript_id", "alternative_transcript_id"]),
        on="pair_key", how="inner",
    )
    tm_retained_rows["benchmark_group"] = NEG_GROUP_TM_RETAINED
    tm_retained_rows["benchmark_label_source"] = "uniprot_isoform_catalog"
    tm_retained_rows["is_true_negative"] = 1
    tm_retained_rows["label_binary"] = int(0)
    tm_retained_rows["label_state"] = "negative"
    _fill_missing_meta(tm_retained_rows, _POS_META_COLS)

    # --- Isoform-level positives: per-isoform TM retained, Cell Membrane retained ---
    positives, pos_mapping_stats = _build_catalog_tm_positive_pairs(base_task, tm_retained_keys)
    positives = base_task.merge(
        positives.drop(columns=["species", "gene_name", "reference_transcript_id", "alternative_transcript_id"]),
        on="pair_key", how="inner",
    )
    positives["benchmark_group"] = POS_GROUP
    positives["benchmark_label_source"] = "uniprot_isoform_catalog"
    positives["is_true_negative"] = 0
    positives["label_binary"] = 1
    positives["label_state"] = "positive"
    _fill_missing_meta(positives, _NEG_META_COLS)

    # Cap positives per split if a ratio is requested
    n_retained_by_split = tm_retained_rows["split"].value_counts().to_dict()
    positives = _cap_positives(positives, float(args.neg_ratio), n_retained_by_split, args.random_seed)

    # --- TM-loss structural controls (no TM at all on the alternative) ---
    tm_retained_counts = tm_retained_rows["split"].value_counts().to_dict()
    control_caps = None
    if float(args.control_multiplier) > 0:
        control_caps = {
            split: int(np.ceil(float(args.control_multiplier) * count))
            for split, count in tm_retained_counts.items()
        }
    pos_keys = set(positives["pair_key"].tolist()) if "pair_key" in positives.columns else set()
    tm_lost_rows = _sample_controls(base_task, tm_retained_keys | pos_keys, control_caps, args.random_seed)
    tm_lost_rows["benchmark_group"] = NEG_GROUP_TM_LOST
    tm_lost_rows["benchmark_label_source"] = "tm_loss_control"
    tm_lost_rows["is_true_negative"] = 0
    tm_lost_rows["label_binary"] = int(0)
    tm_lost_rows["label_state"] = "negative"
    for col in _NEG_META_COLS + _POS_META_COLS + ("primary_accessions",):
        tm_lost_rows[col] = ""
    tm_lost_rows["supporting_catalog_isoform_count"] = 0

    # Ensure primary_accessions exists on all blocks
    for frame in (tm_retained_rows, positives):
        if "primary_accessions" not in frame.columns:
            frame["primary_accessions"] = ""

    mapping_stats = {"negatives": neg_mapping_stats, "positives": pos_mapping_stats}

    out = pd.concat([positives, tm_retained_rows, tm_lost_rows], axis=0, ignore_index=True)
    out = out.drop(columns=["pair_key"]).sort_values(
        ["split", "benchmark_group", "species", "gene_name",
         "reference_transcript_id", "alternative_transcript_id"]
    ).reset_index(drop=True)
    out.to_parquet(OUT_PARQUET, index=False)
    out.to_csv(OUT_TSV, sep="\t", index=False)

    summary = (
        out.groupby(["split", "benchmark_group", "label_state"], dropna=False)
        .size()
        .reset_index(name="n_rows")
        .sort_values(["split", "benchmark_group", "label_state"])
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    schema = {
        "benchmark_name": "tm_isoform_benchmark_v5",
        "task": args.task,
        "positive_definition": "isoform_catalog_tm_and_cell_membrane_length_resolved",
        "negative_definition": "isoform_catalog_tm_retained_cell_membrane_lost_length_resolved",
        "negative_reference_fallback": "catalog_canonical_isoform_id_length_resolved_when_pair_candidates_missing",
        "isoform_to_transcript_resolution": "protein_length_match_within_5pct",
        "source_phase_d_benchmark": str(BASE_BENCH_PATH),
        "source_catalog": str(CATALOG_PATH),
        "source_catalog_full": str(CATALOG_FULL_PATH),
        "source_pair_candidates": str(PAIR_PATH),
        "source_uniprot_map": str(MAP_PATH),
        "feature_columns": feature_columns,
        "categorical_columns": ["family_class", "reference_source", "change_class"],
        "split_column": "split",
        "label_column": "label_binary",
        "label_state_column": "label_state",
        "benchmark_group_column": "benchmark_group",
        "sampling": {
            "control_multiplier": float(args.control_multiplier),
            "neg_ratio": float(args.neg_ratio),
            "random_seed": int(args.random_seed),
        },
        "mapping_stats": mapping_stats,
        "split_counts": {
            split: {
                POS_GROUP: int((out["split"].eq(split) & out["benchmark_group"].eq(POS_GROUP)).sum()),
                NEG_GROUP_TM_RETAINED: int(
                    (out["split"].eq(split) & out["benchmark_group"].eq(NEG_GROUP_TM_RETAINED)).sum()
                ),
                NEG_GROUP_TM_LOST: int(
                    (out["split"].eq(split) & out["benchmark_group"].eq(NEG_GROUP_TM_LOST)).sum()
                ),
            }
            for split in sorted(out["split"].astype(str).unique().tolist())
        },
    }
    OUT_SCHEMA.write_text(json.dumps(schema, indent=2) + "\n", encoding="utf-8")

    print(f"[ok] Wrote {OUT_PARQUET} rows={len(out)}")
    print(f"[ok] Wrote {OUT_TSV}")
    print(f"[ok] Wrote {OUT_SUMMARY}")
    print(f"[ok] Wrote {OUT_SCHEMA}")
    for split, split_df in out.groupby("split"):
        pos_n = int((split_df["label_binary"] == 1).sum())
        neg_ret = int((split_df["benchmark_group"] == NEG_GROUP_TM_RETAINED).sum())
        neg_lost = int((split_df["benchmark_group"] == NEG_GROUP_TM_LOST).sum())
        print(f"  split={split}  positives={pos_n}  tm_retained_neg={neg_ret}  tm_lost_ctrl={neg_lost}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
