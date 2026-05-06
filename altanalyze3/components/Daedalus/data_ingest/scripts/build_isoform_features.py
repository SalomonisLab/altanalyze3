#!/usr/bin/env python3
"""Per-isoform feature matrix derived entirely from UniProt + the
reconstructed protein sequence (no gencode/Ensembl dependencies).

For every isoform in `uniprot_isoform_proteins.tsv` we emit:

  Sequence-derived (computed from the reconstructed protein):
      protein_length
      predicted_tm_count, predicted_tm_total_span, predicted_tm_max_hydropathy_19
      predicted_signal_candidate, predicted_signal_score
      n_terminal_hydropathy_max8
      glyco_motif_count
      cysteine_count, cysteine_fraction
      motif_<name>_count for the 7 trafficking motifs

  UniProt-feature counts intersected with the isoform's surviving residues:
      tm_feature_count, signal_feature_count, signal_max_end
      extracellular_topology_aa, cytoplasmic_topology_aa
      glycosylation_count, disulfide_count, modified_residue_count,
      phospho_feature_count, ppi_binding_feature_count,
      dna_binding_feature_count, dna_region_feature_count,
      zinc_finger_feature_count, domain_count, motif_count

  Entry-level UniProt flags:
      is_membrane_protein, has_signal_peptide, is_kinase,
      is_transcription_factor, interaction_partner_count

Inputs:
    data_ingest/data/interim/uniprot_isoform_proteins.tsv
    data_ingest/data/interim/uniprot_features.tsv
    data_ingest/data/interim/uniprot_entries.tsv
    data_ingest/data/interim/uniprot_functional_priors.tsv (for is_kinase / is_tf)

Output:
    data_ingest/data/interim/uniprot_isoform_features.tsv
"""
from __future__ import annotations

import argparse
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

INTERIM_DIR = Path(__file__).resolve().parents[1] / "data" / "interim"
PROT_PATH = INTERIM_DIR / "uniprot_isoform_proteins.tsv"
ENTRY_PATH = INTERIM_DIR / "uniprot_entries.tsv"
FEATURE_PATH = INTERIM_DIR / "uniprot_features.tsv"
FUNCTIONAL_PATH = INTERIM_DIR / "uniprot_functional_priors.tsv"
OUT_PATH = INTERIM_DIR / "uniprot_isoform_features.tsv"

HYDRO = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4, "H": -3.2, "I": 4.5,
    "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5, "S": -0.8,
    "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}
SMALL_AA = set("AGSCTV")
HYDROPHOBIC_AA = set("AILMFWVYC")
GLYCO_RE = re.compile(r"N[^P][ST]")

MOTIF_PATTERNS = {
    "nls_basic": re.compile(r"[KR]{3,}"),
    "leucine_zipper": re.compile(r"L.{6}L.{6}L"),
    "dileucine": re.compile(r"[DEQRNSTAPG].?LL"),
    "yxxphi": re.compile(r"Y..[AILMFWV]"),
    "gxxxg": re.compile(r"G...G"),
    "kdel_tail": re.compile(r"KDEL$"),
    "kkxx_tail": re.compile(r"KK..$"),
}

# Map UniProt feature_type + description fragment → role
def _feature_role(ftype: str, description: str) -> str | None:
    desc = (description or "").lower()
    if ftype == "Signal":
        return "signal"
    if ftype in {"Transmembrane", "Intramembrane"}:
        return "tm"
    if ftype == "Topological domain":
        if any(x in desc for x in ("extracellular", "luminal", "lumenal", "external")):
            return "extracellular"
        if any(x in desc for x in ("cytoplasmic", "cytosolic", "intracellular")):
            return "cytoplasmic"
    if ftype == "Glycosylation":
        return "glycosylation"
    if ftype == "Disulfide bond":
        return "disulfide"
    if ftype == "Active site":
        return "active_site"
    if ftype == "Modified residue" and "phospho" in desc:
        return "phospho"
    if ftype in {"Binding site", "Site", "Region", "Domain", "Motif", "Coiled coil"}:
        if any(x in desc for x in ("binding", "interact", "dimer", "partner", "sh2", "sh3", "pdz", "ww domain")):
            return "ppi_interface"
    if ftype == "DNA binding":
        return "dna_interface"
    if ftype == "Zinc finger":
        return "dna_interface"
    return None


def _prefix_hydropathy(seq: str) -> list[float]:
    prefix = [0.0]
    total = 0.0
    for aa in seq:
        total += HYDRO.get(aa, 0.0)
        prefix.append(total)
    return prefix


def _window_hydropathy(prefix: list[float], start: int, end: int) -> float:
    end = min(end, len(prefix) - 1)
    length = end - start
    if length <= 0:
        return 0.0
    return (prefix[end] - prefix[start]) / length


def _merge_segments(segments: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not segments:
        return []
    segments.sort()
    merged = [segments[0]]
    for start, end in segments[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end + 1:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def _predict_tm_segments(seq: str, window: int = 19, threshold: float = 1.6) -> list[tuple[int, int]]:
    prefix = _prefix_hydropathy(seq)
    segments: list[tuple[int, int]] = []
    for i in range(0, max(len(seq) - window + 1, 0)):
        frag = seq[i: i + window]
        hydro = _window_hydropathy(prefix, i, i + window)
        hydrophobic_fraction = sum(aa in HYDROPHOBIC_AA for aa in frag) / len(frag)
        if hydro >= threshold and hydrophobic_fraction >= 0.68 and "P" not in frag[:15]:
            segments.append((i + 1, i + window))
    return _merge_segments(segments)


def _predict_signal(seq: str) -> tuple[int, float]:
    nterm = seq[:35]
    if len(nterm) < 15:
        return 0, 0.0
    positive_n = sum(aa in "KR" for aa in nterm[:5])
    prefix = _prefix_hydropathy(nterm)
    best_hydro = max(
        (_window_hydropathy(prefix, i, i + 8) for i in range(2, max(len(nterm) - 7, 3))),
        default=0.0,
    )
    cleavage_bonus = 0.0
    for cut in range(15, min(len(nterm) - 1, 30)):
        tri = nterm[max(cut - 3, 0): cut]
        if len(tri) == 3 and tri[0] in SMALL_AA and tri[2] in SMALL_AA:
            cleavage_bonus = 1.0
            break
    score = 0.0
    if positive_n >= 1:
        score += 1.0
    if best_hydro >= 1.8:
        score += 1.0
    score += cleavage_bonus
    return int(score >= 2.0), score


def _sequence_features(seq: str) -> dict:
    seq = seq.upper()
    if not seq:
        return {}
    tm_segments = _predict_tm_segments(seq)
    sig_call, sig_score = _predict_signal(seq)
    seq_prefix = _prefix_hydropathy(seq)
    nterm = seq[:35]
    nterm_prefix = _prefix_hydropathy(nterm)
    out = {
        "protein_length": len(seq),
        "predicted_tm_count": len(tm_segments),
        "predicted_tm_total_span": sum(end - start + 1 for start, end in tm_segments),
        "predicted_tm_max_hydropathy_19": max(
            (_window_hydropathy(seq_prefix, i, i + 19) for i in range(0, max(len(seq) - 18, 1))),
            default=0.0,
        ),
        "predicted_signal_candidate": sig_call,
        "predicted_signal_score": sig_score,
        "n_terminal_hydropathy_max8": max(
            (_window_hydropathy(nterm_prefix, i, i + 8) for i in range(0, max(len(nterm) - 7, 1))),
            default=0.0,
        ),
        "glyco_motif_count": len(GLYCO_RE.findall(seq)),
        "cysteine_count": seq.count("C"),
        "cysteine_fraction": seq.count("C") / len(seq),
    }
    for name, pat in MOTIF_PATTERNS.items():
        out[f"motif_{name}_count"] = len(pat.findall(seq))
    return out


def _changed_segments_from_alt_seq(alt_seq_changes: str) -> list[tuple[int, int]]:
    """Parse the catalog's `alt_seq_changes` string (pipe-separated VSP entries
    like '1-73: in isoform 2') into list of (begin, end) coordinates on the
    canonical sequence that are deleted/altered in this isoform. Used to
    invert UniProt feature counts: a feature is kept only if it does NOT
    intersect any changed segment."""
    if not alt_seq_changes:
        return []
    segs: list[tuple[int, int]] = []
    for piece in str(alt_seq_changes).split("|"):
        piece = piece.strip()
        if not piece:
            continue
        # "1-73: in isoform 2" -> begin=1 end=73
        head = piece.split(":", 1)[0].strip()
        if "-" in head:
            try:
                a, b = head.split("-", 1)
                segs.append((int(a), int(b)))
            except ValueError:
                continue
    return segs


def _feature_intersects_changed(begin: int, end: int, changed: list[tuple[int, int]]) -> bool:
    for cs, ce in changed:
        if not (end < cs or begin > ce):
            return True
    return False


def _load_uniprot_features() -> pd.DataFrame:
    return pd.read_csv(FEATURE_PATH, sep="\t", low_memory=False)


def _load_entry_flags() -> dict[str, dict]:
    """Per-accession entry-level flags: is_membrane_protein, has_signal_peptide,
    interaction_partner_count, plus is_kinase/is_tf from functional_priors."""
    entries = pd.read_csv(ENTRY_PATH, sep="\t")
    flags: dict[str, dict] = {}
    for r in entries.itertuples(index=False):
        acc = str(r.primary_accession)
        partners = {x for x in str(getattr(r, "interaction_partners", "") or "").split(";") if x}
        sub_locs = str(getattr(r, "subcellular_locations", "") or "").lower()
        flags[acc] = {
            "is_membrane_protein": int("membrane" in sub_locs),
            "has_signal_peptide": 0,  # filled below from features
            "interaction_partner_count": len(partners),
            "is_kinase": 0,
            "is_transcription_factor": 0,
        }
    if FUNCTIONAL_PATH.exists():
        func = pd.read_csv(FUNCTIONAL_PATH, sep="\t")
        for r in func.itertuples(index=False):
            acc = str(r.primary_accession)
            if acc in flags:
                flags[acc]["is_kinase"] = int(getattr(r, "is_kinase", 0) or 0)
                flags[acc]["is_transcription_factor"] = int(
                    getattr(r, "is_transcription_factor", 0) or 0
                )
    return flags


def _isoform_uniprot_features(
    iso_id: str,
    accession: str,
    iso_length: int,
    changed_segments: list[tuple[int, int]],
    features_by_acc: dict[str, list[tuple[str, str, int, int]]],
) -> dict:
    """Aggregate UniProt feature counts/spans for a single isoform.

    A canonical-coordinate feature (begin, end) is treated as 'present' on
    this isoform iff its interval does NOT intersect any changed segment.
    For canonical isoforms, changed_segments=[] so all features are kept.
    """
    bucket: Counter = Counter()
    span_bucket: Counter = Counter()
    role_counts: Counter = Counter()
    role_spans: Counter = Counter()
    for ftype, desc, begin, end in features_by_acc.get(accession, []):
        if begin <= 0 or end < begin:
            continue
        if _feature_intersects_changed(begin, end, changed_segments):
            continue
        span = end - begin + 1
        if ftype == "Signal":
            bucket["signal_feature_count"] += 1
            bucket["signal_max_end"] = max(bucket["signal_max_end"], end)
        elif ftype in {"Transmembrane", "Intramembrane"}:
            bucket["tm_feature_count"] += 1
            bucket["tm_total_span"] += span
        elif ftype == "Topological domain":
            d = desc.lower()
            if any(x in d for x in ("extracellular", "luminal", "lumenal", "external")):
                bucket["extracellular_topology_aa"] += span
                bucket["extracellular_topology_count"] += 1
            elif any(x in d for x in ("cytoplasmic", "cytosolic", "intracellular")):
                bucket["cytoplasmic_topology_aa"] += span
                bucket["cytoplasmic_topology_count"] += 1
        elif ftype == "Glycosylation":
            bucket["glycosylation_count"] += 1
        elif ftype == "Disulfide bond":
            bucket["disulfide_count"] += 1
        elif ftype == "Modified residue":
            bucket["modified_residue_count"] += 1
            if "phospho" in desc.lower():
                bucket["phospho_feature_count"] += 1
        elif ftype in {"Binding site", "Site", "Region", "Domain", "Motif"}:
            d = desc.lower()
            if any(x in d for x in ("binding", "interact", "dimer", "partner", "sh2", "sh3", "pdz", "ww domain")):
                bucket["ppi_binding_feature_count"] += 1
            if any(x in d for x in ("dna-binding", "dna binding", "homeobox", "helix-loop-helix", "bzip")):
                bucket["dna_region_feature_count"] += 1
            if ftype == "Domain":
                bucket["domain_count"] += 1
            if ftype == "Motif":
                bucket["motif_count"] += 1
        elif ftype == "DNA binding":
            bucket["dna_binding_feature_count"] += 1
        elif ftype == "Zinc finger":
            bucket["zinc_finger_feature_count"] += 1

        role = _feature_role(ftype, desc)
        if role:
            role_counts[role] += 1
            role_spans[role] += span

    out = {}
    for k in (
        "signal_feature_count", "signal_max_end",
        "tm_feature_count", "tm_total_span",
        "extracellular_topology_aa", "extracellular_topology_count",
        "cytoplasmic_topology_aa", "cytoplasmic_topology_count",
        "glycosylation_count", "disulfide_count", "modified_residue_count",
        "phospho_feature_count", "ppi_binding_feature_count",
        "dna_binding_feature_count", "dna_region_feature_count",
        "zinc_finger_feature_count", "domain_count", "motif_count",
    ):
        out[k] = int(bucket.get(k, 0))
    # Surviving role counts/spans (used downstream in overlap features —
    # for the alt isoform these represent UniProt features that survive
    # the alt's splice changes; for canonical they're the entry totals)
    for role in (
        "signal", "tm", "extracellular", "cytoplasmic", "glycosylation",
        "disulfide", "active_site", "phospho", "ppi_interface", "dimerization",
        "dna_interface", "nls_region", "nes_region", "activation_region",
        "repression_region", "kinase_region", "trafficking_region",
    ):
        out[f"role_{role}_count"] = int(role_counts.get(role, 0))
        out[f"role_{role}_span"] = int(role_spans.get(role, 0))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=OUT_PATH)
    args = ap.parse_args()

    proteins = pd.read_csv(PROT_PATH, sep="\t")
    print(f"[info] {len(proteins)} isoforms to featurize")

    feat = _load_uniprot_features()
    feat["primary_accession"] = feat["primary_accession"].astype(str)
    features_by_acc: dict[str, list[tuple[str, str, int, int]]] = defaultdict(list)
    for r in feat.itertuples(index=False):
        try:
            begin = int(r.begin) if pd.notna(r.begin) else -1
            end = int(r.end) if pd.notna(r.end) else -1
        except (TypeError, ValueError):
            continue
        if begin <= 0 or end < begin:
            continue
        features_by_acc[str(r.primary_accession)].append(
            (str(r.feature_type or ""), str(r.description or ""), begin, end)
        )
    print(f"[info] features cached for {len(features_by_acc)} accessions")

    entry_flags = _load_entry_flags()
    # Fill in has_signal_peptide from features (canonical-level: any Signal feature)
    for acc, feats in features_by_acc.items():
        if any(f[0] == "Signal" for f in feats):
            if acc in entry_flags:
                entry_flags[acc]["has_signal_peptide"] = 1

    rows: list[dict] = []
    for r in proteins.itertuples(index=False):
        seq = str(r.protein_sequence or "")
        if not seq:
            continue
        seq_feats = _sequence_features(seq)
        changed = _changed_segments_from_alt_seq(str(r.alt_seq_changes or ""))
        uniprot_feats = _isoform_uniprot_features(
            str(r.isoform_id), str(r.primary_accession), len(seq), changed, features_by_acc,
        )
        flags = entry_flags.get(str(r.primary_accession), {
            "is_membrane_protein": 0,
            "has_signal_peptide": 0,
            "interaction_partner_count": 0,
            "is_kinase": 0,
            "is_transcription_factor": 0,
        })
        rows.append({
            "species": str(r.species),
            "gene_name": str(r.gene_name),
            "primary_accession": str(r.primary_accession),
            "isoform_id": str(r.isoform_id),
            "is_canonical": int(r.is_canonical),
            "has_cell_membrane": int(r.has_cell_membrane),
            "has_transmembrane": int(r.has_transmembrane),
            "isoform_locations": str(r.isoform_locations or ""),
            **seq_feats,
            **uniprot_feats,
            **flags,
        })

    out_df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"[ok] wrote {args.out}  rows={len(out_df)}  cols={len(out_df.columns)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
