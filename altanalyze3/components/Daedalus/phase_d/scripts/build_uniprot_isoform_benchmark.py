#!/usr/bin/env python3
"""UniProt-only isoform benchmark.

For each catalog non-canonical isoform we emit one row pairing it with its
gene's canonical isoform. Features are sequence- and UniProt-derived only;
no gencode/Ensembl involvement. Labels:

    label_binary = 1  positive — alternative isoform retains TM AND
                                  Cell Membrane localization
    label_binary = 0  negative — alternative isoform retains TM but
                                  loses Cell Membrane localization
                                  (the canonical's locations are TM+CM)

Splits are gene-grouped: every isoform of a gene is in the same split.
Where the existing phase_d benchmark covers the gene we reuse its split
assignment for comparability; otherwise we hash the gene name into one of
{train, val, test} with an 80/10/10 split.

Output schema: 154-feature parquet matching the existing kfold runner's
expectations as closely as possible.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
DAEDALUS = ROOT / "Daedalus"
PHASE_A = DAEDALUS / "phase_a" / "data" / "interim"
PHASE_D = DAEDALUS / "phase_d" / "data" / "processed"

ISOFORM_FEATURES = PHASE_A / "uniprot_isoform_features.tsv"
ISOFORM_PROTEINS = PHASE_A / "uniprot_isoform_proteins.tsv"
TRAFFICKING_FEATURES = PHASE_A / "uniprot_isoform_trafficking_features.tsv"
ESM_EMBEDDINGS = PHASE_A / "uniprot_isoform_esm2_t33_embeddings.parquet"
LOCALIZATION_PROBS = PHASE_A / "uniprot_isoform_localization_probs.tsv"
CANONICAL_PM_EVIDENCE = PHASE_A / "uniprot_canonical_pm_evidence.tsv"
EXISTING_BENCH = PHASE_D / "phase_d_benchmark_instances.parquet"

OUT_PARQUET = PHASE_D / "uniprot_isoform_benchmark_instances.parquet"
OUT_TSV = PHASE_D / "uniprot_isoform_benchmark_instances.tsv"
OUT_SCHEMA = PHASE_D / "uniprot_isoform_benchmark_schema.json"
OUT_SUMMARY = PHASE_D / "uniprot_isoform_benchmark_summary.tsv"

POS_GROUP = "uniprot_tm_positive"
NEG_GROUP = "tm_retained_non_surface_true_negative"

# Columns in the per-isoform feature table that should NOT participate as
# either reference or alternative features (they're identifiers / labels).
META_COLS = {
    "species", "gene_name", "primary_accession", "isoform_id",
    "is_canonical", "has_cell_membrane", "has_transmembrane",
    "isoform_locations",
}


def _hash_split(species: str, gene: str, train_frac: float = 0.8, val_frac: float = 0.1) -> str:
    """Deterministic gene-grouped split assignment by md5 hash bucket."""
    digest = hashlib.md5(f"{species}|{gene}".encode()).hexdigest()
    bucket = int(digest[:8], 16) / 0xFFFFFFFF
    if bucket < train_frac:
        return "train"
    if bucket < train_frac + val_frac:
        return "val"
    return "test"


def _existing_gene_splits() -> dict[tuple[str, str], str]:
    """Pull (species, gene) -> split from the existing phase_d benchmark for
    consistency with prior kfold/stacking runs."""
    if not EXISTING_BENCH.exists():
        return {}
    b = pd.read_parquet(EXISTING_BENCH)
    b = b[b["task"] == "cell_surface_localization"]
    out: dict[tuple[str, str], str] = {}
    for r in b[["species", "gene_name", "split"]].drop_duplicates().itertuples(index=False):
        out[(str(r.species).lower(), str(r.gene_name))] = str(r.split)
    return out


_SOLUBLE_TOKENS = (
    "secreted", "cytoplasm", "cytosol", "nucleus", "nuclear",
    "endomembrane system",
)
_OTHER_MEMBRANE_TOKENS = (
    # Any "<organelle> membrane" that is NOT plasma membrane. These isoforms
    # are TM proteins on intracellular membranes; they're indistinguishable
    # from PM proteins by sequence/feature signal and add label noise.
    "endoplasmic reticulum membrane", "er membrane",
    "golgi", "lysosome membrane", "endosome membrane",
    "mitochondrion", "mitochondrial",
    "lipid droplet", "peroxisome membrane",
    "vesicle membrane", "vesicular membrane",
    "myelin membrane", "ruffle membrane",
    "podosome", "phagosome",
    "early endosome membrane", "late endosome membrane",
)
_PM_LIKE_TOKENS = (
    "cell membrane", "plasma membrane", "cell surface", "cell junction",
    "apical", "basolateral",
)


def _classify_negative_locations(loc_str: str) -> str:
    """Classify a non-canonical isoform's `isoform_locations` string into:
        soluble  — clearly non-membrane (secreted, cytoplasm, nucleus, ER lumen)
        organellar_membrane — TM protein but on a non-PM organelle membrane
        pm_like  — annotated at PM-adjacent locations (suspicious negative)
        missing  — no annotation
    """
    if loc_str is None:
        return "missing"
    s = str(loc_str).lower().strip()
    if not s or s == "nan":
        return "missing"
    if any(tok in s for tok in _PM_LIKE_TOKENS):
        return "pm_like"
    if any(tok in s for tok in _OTHER_MEMBRANE_TOKENS):
        return "organellar_membrane"
    if "membrane" in s:
        # generic "membrane" without organelle qualifier — still ambiguous
        return "organellar_membrane"
    if any(tok in s for tok in _SOLUBLE_TOKENS):
        return "soluble"
    # Anything else (extracellular matrix variants, paranodal junction etc.)
    # — default to soluble; rare and not PM.
    return "soluble"


def _label_for_isoform(row, *, exclude_organellar: bool = True,
                       exclude_missing: bool = True,
                       exclude_pm_like: bool = True) -> tuple[str, int]:
    """Return (benchmark_group, label_binary) for a non-canonical isoform.

    Positives: alt isoform retains TM AND has Cell Membrane localization
    (per ``has_cell_membrane`` from the catalog).
    Negatives: subdivided by ``isoform_locations`` into soluble (kept),
    organellar_membrane (excluded by default — these are TM proteins on
    other organelles and add label noise), pm_like (excluded — likely
    misannotated as negative), and missing (excluded by default).
    """
    if row["has_transmembrane"] != 1:
        return ("excluded_tm_lost", -1)
    if row["has_cell_membrane"] == 1:
        return (POS_GROUP, 1)
    raw_loc = row.get("isoform_locations", None)
    if pd.isna(raw_loc):
        raw_loc = None
    cls = _classify_negative_locations(raw_loc)
    if cls == "soluble":
        return (NEG_GROUP, 0)
    if cls == "organellar_membrane" and not exclude_organellar:
        return (NEG_GROUP, 0)
    if cls == "missing" and not exclude_missing:
        return (NEG_GROUP, 0)
    if cls == "pm_like" and not exclude_pm_like:
        return (NEG_GROUP, 0)
    return (f"excluded_{cls}", -1)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep-organellar-negatives", action="store_true",
                    help="Include alt isoforms on non-PM organellar membranes as negatives "
                         "(default: excluded — they're indistinguishable from PM proteins by feature signal).")
    ap.add_argument("--keep-missing-negatives", action="store_true",
                    help="Include alt isoforms with no `isoform_locations` annotation as negatives "
                         "(default: excluded — annotation is required).")
    ap.add_argument("--keep-pm-like-negatives", action="store_true",
                    help="Include alt isoforms tagged Cell-surface/Cell-junction/etc as negatives "
                         "(default: excluded — almost certainly mislabelled).")
    ap.add_argument("--strict-canonical-pm", action="store_true",
                    help="Drop genes whose canonical reference's UniProt subcellular_location "
                         "lists both PM and another major compartment (Cytoplasm / Nucleus / "
                         "Mitochondrion / Lysosome). These ambiguous canonicals add label "
                         "noise because the 'reference' itself is multi-compartment.")
    args = ap.parse_args()

    feats = pd.read_csv(ISOFORM_FEATURES, sep="\t")
    print(f"[info] loaded {len(feats)} isoform feature rows ({len(feats.columns)} cols)")

    # Merge trafficking biology features (topology integrity, ER export /
    # retention motifs, folding QC) keyed by (species, gene, isoform_id).
    if TRAFFICKING_FEATURES.exists():
        traff = pd.read_csv(TRAFFICKING_FEATURES, sep="\t")
        traff_meta = {"species", "gene_name", "primary_accession", "isoform_id", "is_canonical"}
        traff_cols = [c for c in traff.columns if c not in traff_meta]
        print(f"[info] merging {len(traff_cols)} trafficking biology features")
        feats = feats.merge(
            traff[["species", "gene_name", "isoform_id"] + traff_cols],
            on=["species", "gene_name", "isoform_id"], how="left",
        )
        for c in traff_cols:
            feats[c] = feats[c].fillna(0)
        print(f"[info] post-merge feature rows: {len(feats)} ({len(feats.columns)} cols)")
    else:
        print(f"[warn] trafficking features file not found at {TRAFFICKING_FEATURES}; skipping")

    # Identify the canonical row per (species, gene)
    canonicals = feats[feats["is_canonical"] == 1].copy()
    canonicals = canonicals.drop_duplicates(subset=["species", "gene_name"])
    canon_idx: dict[tuple[str, str], dict] = {
        (str(r["species"]).lower(), str(r["gene_name"])): r.to_dict()
        for _, r in canonicals.iterrows()
    }
    print(f"[info] canonical rows indexed: {len(canon_idx)}")

    # Optional stricter filter: drop genes whose canonical UniProt entry has
    # ambiguous subcellular_location (PM + another compartment). This removes
    # the "reference is multi-compartment" source of label noise.
    excluded_ambiguous_canonicals = 0
    if args.strict_canonical_pm and CANONICAL_PM_EVIDENCE.exists():
        ev = pd.read_csv(CANONICAL_PM_EVIDENCE, sep="\t")
        # Clean reference = has PM AND no other major compartment
        clean_acc = set(ev[(ev["has_pm"] == True) & (ev["has_other_compartment"] == False)]["primary_accession"].astype(str))
        new_canon_idx = {}
        for k, v in canon_idx.items():
            acc = str(v.get("primary_accession", ""))
            if acc in clean_acc:
                new_canon_idx[k] = v
            else:
                excluded_ambiguous_canonicals += 1
        canon_idx = new_canon_idx
        print(f"[info] strict_canonical_pm: kept {len(canon_idx)} clean-PM canonicals "
              f"(excluded {excluded_ambiguous_canonicals} ambiguous)")

    feature_cols = [c for c in feats.columns if c not in META_COLS]
    print(f"[info] feature columns: {len(feature_cols)}")

    existing_splits = _existing_gene_splits()
    print(f"[info] gene-split assignments inherited from existing benchmark: {len(existing_splits)}")

    # Load ESM2 protein-language-model embeddings (480-dim mean-pooled).
    # Indexed by (species, gene_name, isoform_id); used as alternative + delta
    # features against the canonical isoform's embedding.
    esm_alt_idx: dict[tuple[str, str, str], np.ndarray] = {}
    esm_canon_by_gene: dict[tuple[str, str], np.ndarray] = {}
    esm_cols: list[str] = []
    if ESM_EMBEDDINGS.exists():
        esm_df = pd.read_parquet(ESM_EMBEDDINGS)
        esm_cols = [c for c in esm_df.columns if c.startswith("esm2_")]
        print(f"[info] loaded {len(esm_df)} ESM2 embeddings ({len(esm_cols)} dims)")
        # Build per-isoform map and canonical-per-gene map
        for r in esm_df.itertuples(index=False):
            arr = np.array([getattr(r, c) for c in esm_cols], dtype=np.float32)
            sp = str(r.species).lower(); gn = str(r.gene_name); iso = str(r.isoform_id)
            esm_alt_idx[(sp, gn, iso)] = arr
        # Canonical per gene = the isoform_id that is_canonical==1 in feats
        for _, row in feats[feats["is_canonical"] == 1].drop_duplicates(["species", "gene_name"]).iterrows():
            sp = str(row["species"]).lower(); gn = str(row["gene_name"])
            iso = str(row["isoform_id"])
            arr = esm_alt_idx.get((sp, gn, iso))
            if arr is not None:
                esm_canon_by_gene[(sp, gn)] = arr
        print(f"[info] canonical ESM2 embeddings indexed for {len(esm_canon_by_gene)} genes")
    else:
        print(f"[warn] ESM2 embeddings file not found at {ESM_EMBEDDINGS}; skipping")

    # Load DeepLoc-style 10-class subcellular-localization probabilities
    loc_alt_idx: dict[tuple[str, str, str], np.ndarray] = {}
    loc_canon_by_gene: dict[tuple[str, str], np.ndarray] = {}
    loc_cols: list[str] = []
    if LOCALIZATION_PROBS.exists():
        loc_df = pd.read_csv(LOCALIZATION_PROBS, sep="\t")
        loc_cols = [c for c in loc_df.columns if c.startswith("loc_prob_")]
        print(f"[info] loaded {len(loc_df)} DeepLoc-style localization-probability rows ({len(loc_cols)} classes)")
        for r in loc_df.itertuples(index=False):
            arr = np.array([getattr(r, c) for c in loc_cols], dtype=np.float32)
            sp = str(r.species).lower(); gn = str(r.gene_name); iso = str(r.isoform_id)
            loc_alt_idx[(sp, gn, iso)] = arr
        for _, row in feats[feats["is_canonical"] == 1].drop_duplicates(["species", "gene_name"]).iterrows():
            sp = str(row["species"]).lower(); gn = str(row["gene_name"])
            iso = str(row["isoform_id"])
            arr = loc_alt_idx.get((sp, gn, iso))
            if arr is not None:
                loc_canon_by_gene[(sp, gn)] = arr
        print(f"[info] canonical DeepLoc probs indexed for {len(loc_canon_by_gene)} genes")
    else:
        print(f"[warn] localization probs file not found at {LOCALIZATION_PROBS}; skipping")

    rows: list[dict] = []
    skipped_no_canonical = 0
    excluded_counts: dict[str, int] = {}
    for _, row in feats[feats["is_canonical"] == 0].iterrows():
        key = (str(row["species"]).lower(), str(row["gene_name"]))
        canon = canon_idx.get(key)
        if canon is None:
            skipped_no_canonical += 1
            continue
        group, label = _label_for_isoform(
            row,
            exclude_organellar=not args.keep_organellar_negatives,
            exclude_missing=not args.keep_missing_negatives,
            exclude_pm_like=not args.keep_pm_like_negatives,
        )
        if label == -1:
            excluded_counts[group] = excluded_counts.get(group, 0) + 1
            continue
        # Categorical stubs retained for compatibility with the existing
        # kfold runner's encoder; in UniProt mode they carry no signal.
        change_class = "unknown"
        out_row = {
            "species": str(row["species"]).lower(),
            "gene_name": str(row["gene_name"]),
            "primary_accession": str(row["primary_accession"]),
            "reference_isoform_id": str(canon["isoform_id"]),
            "alternative_isoform_id": str(row["isoform_id"]),
            "benchmark_group": group,
            "label_binary": int(label),
            "label_state": "positive" if label == 1 else "negative",
            "is_true_negative": int(label == 0),
            "alt_isoform_locations": str(row["isoform_locations"] or ""),
            "ref_isoform_locations": str(canon.get("isoform_locations") or ""),
            "family_class": "unknown",
            "reference_source": "uniprot_canonical",
            "change_class": change_class,
        }
        # Reference-keyed features (count statistics + trafficking features)
        for c in feature_cols:
            ref_v = canon.get(c)
            alt_v = row[c]
            try:
                ref_f = float(ref_v) if ref_v is not None else 0.0
                alt_f = float(alt_v) if alt_v is not None else 0.0
            except (TypeError, ValueError):
                ref_f, alt_f = 0.0, 0.0
            out_row[f"reference_{c}"] = ref_f
            out_row[f"alternative_{c}"] = alt_f
            out_row[f"delta_{c}"] = alt_f - ref_f

        # ESM2 protein-language-model embeddings (alt + delta only;
        # reference is constant per gene and would inflate the matrix).
        if esm_cols:
            alt_emb = esm_alt_idx.get((str(row["species"]).lower(),
                                        str(row["gene_name"]),
                                        str(row["isoform_id"])))
            ref_emb = esm_canon_by_gene.get(key)
            if alt_emb is not None and ref_emb is not None:
                delta = alt_emb - ref_emb
                for i, c in enumerate(esm_cols):
                    out_row[f"alternative_{c}"] = float(alt_emb[i])
                    out_row[f"delta_{c}"] = float(delta[i])
            else:
                for c in esm_cols:
                    out_row[f"alternative_{c}"] = 0.0
                    out_row[f"delta_{c}"] = 0.0

        # DeepLoc-style subcellular-localization probabilities
        # (10 classes: alt + delta-from-canonical)
        if loc_cols:
            alt_loc = loc_alt_idx.get((str(row["species"]).lower(),
                                        str(row["gene_name"]),
                                        str(row["isoform_id"])))
            ref_loc = loc_canon_by_gene.get(key)
            if alt_loc is not None and ref_loc is not None:
                delta = alt_loc - ref_loc
                for i, c in enumerate(loc_cols):
                    out_row[f"alternative_{c}"] = float(alt_loc[i])
                    out_row[f"delta_{c}"] = float(delta[i])
            else:
                for c in loc_cols:
                    out_row[f"alternative_{c}"] = 0.0
                    out_row[f"delta_{c}"] = 0.0

        out_row["split"] = existing_splits.get(key, _hash_split(*key))
        rows.append(out_row)

    out = pd.DataFrame(rows)
    print(f"[info] benchmark rows: {len(out)}  (skipped no_canonical={skipped_no_canonical})")
    if excluded_counts:
        for k, v in sorted(excluded_counts.items()):
            print(f"  excluded[{k}]: {v}")

    # Sort and write outputs
    out = out.sort_values(
        ["split", "benchmark_group", "species", "gene_name",
         "reference_isoform_id", "alternative_isoform_id"]
    ).reset_index(drop=True)

    OUT_PARQUET.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUT_PARQUET, index=False)
    out.to_csv(OUT_TSV, sep="\t", index=False)

    summary = (
        out.groupby(["split", "benchmark_group", "label_state"], dropna=False)
        .size().reset_index(name="n_rows")
        .sort_values(["split", "benchmark_group", "label_state"])
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    feature_columns_full = []
    for c in feature_cols:
        for prefix in ("reference_", "alternative_", "delta_"):
            feature_columns_full.append(f"{prefix}{c}")
    # ESM2 features (alt + delta only)
    for c in esm_cols:
        feature_columns_full.append(f"alternative_{c}")
        feature_columns_full.append(f"delta_{c}")
    # DeepLoc-style localization probabilities (alt + delta only)
    for c in loc_cols:
        feature_columns_full.append(f"alternative_{c}")
        feature_columns_full.append(f"delta_{c}")

    schema = {
        "benchmark_name": "uniprot_isoform_benchmark_v1",
        "feature_source": "uniprot_only",
        "task": "cell_surface_localization",
        "positive_definition": "isoform_catalog_tm_and_cell_membrane",
        "negative_definition": "isoform_catalog_tm_retained_cell_membrane_lost",
        "feature_columns": feature_columns_full,
        "categorical_columns": ["family_class", "reference_source", "change_class"],
        "split_column": "split",
        "label_column": "label_binary",
        "label_state_column": "label_state",
        "benchmark_group_column": "benchmark_group",
        "row_id_columns": [
            "species", "gene_name", "primary_accession",
            "reference_isoform_id", "alternative_isoform_id",
        ],
        "split_counts": {
            split: {
                POS_GROUP: int(((out["split"] == split) & (out["benchmark_group"] == POS_GROUP)).sum()),
                NEG_GROUP: int(((out["split"] == split) & (out["benchmark_group"] == NEG_GROUP)).sum()),
            }
            for split in sorted(out["split"].unique().tolist())
        },
    }
    OUT_SCHEMA.write_text(json.dumps(schema, indent=2) + "\n", encoding="utf-8")
    print(f"[ok] wrote {OUT_PARQUET}")
    print(f"[ok] wrote {OUT_TSV}")
    print(f"[ok] wrote {OUT_SUMMARY}")
    print(f"[ok] wrote {OUT_SCHEMA}")
    for split, sdf in out.groupby("split"):
        pos_n = int((sdf["benchmark_group"] == POS_GROUP).sum())
        neg_n = int((sdf["benchmark_group"] == NEG_GROUP).sum())
        print(f"  split={split}  positives={pos_n}  negatives={neg_n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
