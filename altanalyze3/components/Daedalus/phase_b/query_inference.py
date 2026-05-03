from __future__ import annotations

import csv
import gzip
import json
import math
import re
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import torch

from Daedalus.phase_b.models import ReferenceDeltaNet


ROOT = Path(__file__).resolve().parents[1]
PHASE_A_INTERIM = ROOT / "phase_a" / "data" / "interim"
PHASE_A_RAW = ROOT / "phase_a" / "data" / "raw"
PHASE_B_CHECKPOINTS = ROOT / "phase_b" / "checkpoints"
DEFAULT_CHECKPOINT = PHASE_B_CHECKPOINTS / "reference_delta_multitask.pt"

HYDRO = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}
SMALL_AA = set("AGSCTV")
GLYCO_RE = re.compile(r"N[^P][ST]")


def _base_id(identifier: str) -> str:
    return identifier.split(".", 1)[0] if identifier else ""


def _parse_int(value: Any) -> int:
    if value in (None, ""):
        return 0
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def _parse_float(value: Any) -> float:
    if value in (None, ""):
        return 0.0
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _safe_div(numer: float, denom: float, default: float = 0.0) -> float:
    if abs(denom) < 1e-12:
        return default
    return numer / denom


def _clamp01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))


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
        frag = seq[i : i + window]
        hydro = _window_hydropathy(prefix, i, i + window)
        hydrophobic_fraction = sum(aa in "AILMFWVYC" for aa in frag) / len(frag)
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
        tri = nterm[max(cut - 3, 0) : cut]
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


def _reference_source(is_mane_select: bool, appris_tag: str, protein_length: int, n_uniprot_links: int) -> str:
    if is_mane_select:
        return "MANE_SELECT"
    if appris_tag == "PRINCIPAL:1":
        return "APPRIS_PRINCIPAL_1"
    if appris_tag.startswith("PRINCIPAL"):
        return "APPRIS_PRINCIPAL"
    if n_uniprot_links > 0 and protein_length > 0:
        return "UNIPROT_SUPPORTED_LONGEST"
    if protein_length > 0:
        return "LONGEST_PROTEIN"
    return "LONGEST_TRANSCRIPT"


def _parse_appris_rank(tag: str) -> tuple[int, int]:
    if not tag:
        return (9, 9)
    if tag.startswith("PRINCIPAL:"):
        try:
            return (0, int(tag.split(":", 1)[1]))
        except ValueError:
            return (0, 9)
    if tag.startswith("ALTERNATIVE:"):
        try:
            return (1, int(tag.split(":", 1)[1]))
        except ValueError:
            return (1, 9)
    return (9, 9)


def _iter_tsv(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:  # type: ignore[arg-type]
        reader = csv.DictReader(handle, delimiter="\t")
        yield from reader


def _load_transcript_catalog() -> list[dict[str, str]]:
    return list(_iter_tsv(PHASE_A_INTERIM / "transcript_supervision_catalog.tsv"))


def _load_gene_catalog() -> dict[tuple[str, str], dict[str, str]]:
    rows: dict[tuple[str, str], dict[str, str]] = {}
    for row in _iter_tsv(PHASE_A_INTERIM / "gene_supervision_catalog.tsv"):
        rows[(row["species"], _base_id(row["gencode_gene_id"]))] = row
    return rows


def _load_gencode_reference_lengths() -> dict[tuple[str, str], dict[str, str]]:
    rows: dict[tuple[str, str], dict[str, str]] = {}
    for row in _iter_tsv(PHASE_A_INTERIM / "gencode_transcript_reference.tsv"):
        rows[(row["species"], row["transcript_id"])] = row
    return rows


def _lookup_gene_row(
    species: str,
    gene_id: str | None,
    gene_name: str | None,
    gene_catalog: dict[tuple[str, str], dict[str, str]],
) -> dict[str, str]:
    if gene_id:
        row = gene_catalog.get((species, _base_id(gene_id)))
        if row:
            return row
    if gene_name:
        for (row_species, _), row in gene_catalog.items():
            if row_species == species and row.get("gene_name", "") == gene_name:
                return row
    raise ValueError(f"Could not resolve gene for species={species!r}, gene_id={gene_id!r}, gene_name={gene_name!r}")


def select_reference_transcript(
    species: str,
    gene_id: str | None = None,
    gene_name: str | None = None,
    reference_transcript_id: str | None = None,
) -> dict[str, Any]:
    transcript_rows = _load_transcript_catalog()
    gene_catalog = _load_gene_catalog()
    gene_row = _lookup_gene_row(species, gene_id, gene_name, gene_catalog)
    gene_base = _base_id(gene_row["gencode_gene_id"])
    candidates = [
        row
        for row in transcript_rows
        if row["species"] == species and _base_id(row["gene_id"]) == gene_base
    ]
    if not candidates:
        raise ValueError(f"No transcript candidates found for gene={gene_row['gene_name']} ({gene_base})")
    if reference_transcript_id:
        selected = next(
            (row for row in candidates if _base_id(row["transcript_id"]) == _base_id(reference_transcript_id)),
            None,
        )
        if selected is None:
            raise ValueError(f"Reference transcript {reference_transcript_id} not found for gene {gene_row['gene_name']}")
    else:
        candidates.sort(
            key=lambda row: (
                -_parse_int(row.get("is_mane_select", "0")),
                _parse_appris_rank(row.get("appris_tag", "")),
                -_parse_int(row.get("n_uniprot_links", "0")),
                -_parse_int(row.get("protein_length", "0")),
                row["transcript_id"],
            )
        )
        selected = candidates[0]
    protein_length = _parse_int(selected.get("protein_length", "0"))
    source = _reference_source(
        _parse_int(selected.get("is_mane_select", "0")) == 1,
        selected.get("appris_tag", ""),
        protein_length,
        _parse_int(selected.get("n_uniprot_links", "0")),
    )
    return {
        "gene": gene_row,
        "reference": selected,
        "reference_source": source,
        "selection_reason": {
            "is_mane_select": _parse_int(selected.get("is_mane_select", "0")),
            "appris_tag": selected.get("appris_tag", ""),
            "n_uniprot_links": _parse_int(selected.get("n_uniprot_links", "0")),
            "protein_length": protein_length,
        },
    }


def derive_query_sequence_features(protein_seq: str) -> dict[str, Any]:
    seq = protein_seq.strip().upper()
    tm_segments = _predict_tm_segments(seq)
    signal_candidate, signal_score = _predict_signal(seq)
    seq_prefix = _prefix_hydropathy(seq)
    nterm_prefix = _prefix_hydropathy(seq[:35])
    return {
        "protein_length": len(seq),
        "predicted_tm_count": len(tm_segments),
        "predicted_tm_total_span": sum(end - start + 1 for start, end in tm_segments),
        "predicted_tm_max_hydropathy_19": max(
            (_window_hydropathy(seq_prefix, i, i + 19) for i in range(0, max(len(seq) - 18, 1))),
            default=0.0,
        ),
        "predicted_signal_candidate": signal_candidate,
        "predicted_signal_score": signal_score,
        "n_terminal_hydropathy_max8": max(
            (_window_hydropathy(nterm_prefix, i, i + 8) for i in range(0, max(len(seq[:35]) - 7, 1))),
            default=0.0,
        ),
        "glyco_motif_count": len(GLYCO_RE.findall(seq)),
        "cysteine_count": seq.count("C"),
        "cysteine_fraction": _safe_div(seq.count("C"), len(seq), 0.0),
        "tm_segments": tm_segments,
    }


def _family_class(reference: dict[str, str], alt_features: dict[str, Any]) -> str:
    if _parse_int(reference.get("is_kinase", "0")) == 1:
        return "kinase"
    if _parse_int(reference.get("is_transcription_factor", "0")) == 1:
        return "transcription_factor"
    if _parse_int(reference.get("is_membrane_protein", "0")) == 1 or alt_features["predicted_tm_count"] > 0 or alt_features["predicted_signal_candidate"] == 1:
        return "membrane"
    return "other"


def _load_biogrid_partners(gene_name: str, max_partners: int = 25) -> list[dict[str, Any]]:
    zip_path = PHASE_A_RAW / "BIOGRID-ALL-LATEST.tab3.zip"
    if not zip_path.exists() or not gene_name:
        return []
    partner_counts: dict[str, dict[str, int]] = defaultdict(lambda: {"all": 0, "physical": 0})
    with zipfile.ZipFile(zip_path) as zf:
        member = next((name for name in zf.namelist() if name.endswith(".tab3.txt")), None)
        if member is None:
            return []
        with zf.open(member) as handle:
            text = (line.decode("utf-8", errors="replace") for line in handle)
            reader = csv.DictReader(text, delimiter="\t")
            for row in reader:
                a = row.get("Official Symbol Interactor A", "")
                b = row.get("Official Symbol Interactor B", "")
                if gene_name not in {a, b}:
                    continue
                partner = b if a == gene_name else a
                if not partner:
                    continue
                partner_counts[partner]["all"] += 1
                if row.get("Experimental System Type", "").lower() == "physical":
                    partner_counts[partner]["physical"] += 1
    ranked = sorted(
        (
            {"partner": partner, "all_interactions": counts["all"], "physical_interactions": counts["physical"]}
            for partner, counts in partner_counts.items()
        ),
        key=lambda row: (-row["physical_interactions"], -row["all_interactions"], row["partner"]),
    )
    return ranked[:max_partners]


def _build_numeric_row(
    reference_context: dict[str, Any],
    alt_features: dict[str, Any],
    alt_protein_length: int,
    alt_transcript_length: int,
) -> dict[str, float]:
    reference = reference_context["reference"]
    gene = reference_context["gene"]
    ref_protein_length = _parse_int(reference.get("protein_length", "0"))
    ref_transcript_length = _parse_int(reference_context["reference_lengths"].get("transcript_length", "0"))
    reference_predicted_tm_count = _parse_int(reference.get("predicted_tm_count", "0"))
    reference_predicted_tm_total_span = _parse_int(reference.get("predicted_tm_total_span", "0"))
    reference_predicted_signal_candidate = _parse_int(reference.get("predicted_signal_candidate", "0"))
    reference_predicted_signal_score = _parse_float(reference.get("predicted_signal_score", "0"))
    return {
        "reference_protein_length": float(ref_protein_length),
        "alternative_protein_length": float(alt_protein_length),
        "reference_transcript_length": float(ref_transcript_length),
        "alternative_transcript_length": float(alt_transcript_length),
        "reference_uniprot_links": float(_parse_int(reference.get("n_uniprot_links", "0"))),
        "reference_is_membrane": float(_parse_int(reference.get("is_membrane_protein", "0"))),
        "reference_has_signal_peptide": float(_parse_int(reference.get("has_signal_peptide", "0"))),
        "reference_is_kinase": float(_parse_int(reference.get("is_kinase", "0"))),
        "reference_is_tf": float(_parse_int(reference.get("is_transcription_factor", "0"))),
        "gene_has_appris_principal_1": float(_parse_int(gene.get("has_appris_principal_1", "0"))),
        "gene_clinvar_pathogenic_splice_count": float(_parse_int(gene.get("clinvar_pathogenic_splice_count", "0"))),
        "gene_hpa_has_extracellular_annotation": float(_parse_int(gene.get("hpa_has_extracellular_annotation", "0"))),
        "gene_biogrid_partner_count": float(_parse_int(gene.get("biogrid_partner_count", "0"))),
        "gene_biogrid_physical_partner_count": float(_parse_int(gene.get("biogrid_physical_partner_count", "0"))),
        "protein_length_delta": float(alt_protein_length - ref_protein_length),
        "protein_length_ratio": float(_safe_div(alt_protein_length, ref_protein_length, 0.0)),
        "transcript_length_delta": float(alt_transcript_length - ref_transcript_length),
        "ref_alt_same_protein_length": float(int(alt_protein_length == ref_protein_length and alt_protein_length > 0)),
        "alt_has_protein": float(int(alt_protein_length > 0)),
        "reference_tm_feature_count": float(_parse_int(reference.get("tm_feature_count", "0"))),
        "reference_signal_feature_count": float(_parse_int(reference.get("signal_feature_count", "0"))),
        "reference_extracellular_topology_aa": float(_parse_int(reference.get("extracellular_topology_aa", "0"))),
        "reference_cytoplasmic_topology_aa": float(_parse_int(reference.get("cytoplasmic_topology_aa", "0"))),
        "reference_glycosylation_count": float(_parse_int(reference.get("glycosylation_count", "0"))),
        "reference_disulfide_count": float(_parse_int(reference.get("disulfide_count", "0"))),
        "reference_phospho_feature_count": float(_parse_int(reference.get("phospho_feature_count", "0"))),
        "reference_ppi_binding_feature_count": float(_parse_int(reference.get("ppi_binding_feature_count", "0"))),
        "reference_dna_binding_feature_count": float(_parse_int(reference.get("dna_binding_feature_count", "0"))),
        "reference_dna_region_feature_count": float(_parse_int(reference.get("dna_region_feature_count", "0"))),
        "reference_zinc_finger_feature_count": float(_parse_int(reference.get("zinc_finger_feature_count", "0"))),
        "reference_domain_count": float(_parse_int(reference.get("domain_count", "0"))),
        "reference_motif_count": float(_parse_int(reference.get("motif_count", "0"))),
        "alternative_predicted_tm_count": float(alt_features["predicted_tm_count"]),
        "reference_predicted_tm_count": float(reference_predicted_tm_count),
        "predicted_tm_count_delta": float(alt_features["predicted_tm_count"] - reference_predicted_tm_count),
        "reference_predicted_tm_total_span": float(reference_predicted_tm_total_span),
        "alternative_predicted_tm_total_span": float(alt_features["predicted_tm_total_span"]),
        "predicted_tm_total_span_delta": float(alt_features["predicted_tm_total_span"] - reference_predicted_tm_total_span),
        "reference_predicted_signal_candidate": float(reference_predicted_signal_candidate),
        "alternative_predicted_signal_candidate": float(alt_features["predicted_signal_candidate"]),
        "reference_predicted_signal_score": float(reference_predicted_signal_score),
        "alternative_predicted_signal_score": float(alt_features["predicted_signal_score"]),
        "predicted_signal_score_delta": float(alt_features["predicted_signal_score"] - reference_predicted_signal_score),
        "alternative_glyco_motif_count": float(alt_features["glyco_motif_count"]),
        "alternative_cysteine_count": float(alt_features["cysteine_count"]),
    }


def _encode_numeric(
    numeric_row: dict[str, float],
    numeric_columns: list[str],
    mean: np.ndarray,
    std: np.ndarray,
) -> np.ndarray:
    raw = np.array([[float(numeric_row.get(column, 0.0)) for column in numeric_columns]], dtype=np.float32)
    return (raw - mean.astype(np.float32)) / std.astype(np.float32)


def _task_applicability(reference: dict[str, str], alt_features: dict[str, Any], gene: dict[str, str]) -> dict[str, bool]:
    membrane_like = (
        _parse_int(reference.get("is_membrane_protein", "0")) == 1
        or alt_features["predicted_tm_count"] > 0
        or alt_features["predicted_signal_candidate"] == 1
    )
    return {
        "global": True,
        "membrane": membrane_like,
        "surface": membrane_like,
        "kinase": _parse_int(reference.get("is_kinase", "0")) == 1,
        "transcription_factor": _parse_int(reference.get("is_transcription_factor", "0")) == 1,
    }


def _compute_task_probabilities(
    checkpoint: dict[str, Any],
    numeric_row: dict[str, float],
    species: str,
    family_class: str,
    reference_source: str,
) -> dict[str, float]:
    metadata = checkpoint["metadata"]
    numeric_columns = metadata["numeric_columns"]
    mean = np.asarray(metadata["mean"], dtype=np.float32)
    std = np.asarray(metadata["std"], dtype=np.float32)
    model = ReferenceDeltaNet(
        n_numeric=len(numeric_columns),
        n_tasks=len(metadata["task_map"]),
        n_species=len(metadata["species_map"]),
        n_families=len(metadata["family_map"]),
        n_reference_sources=len(metadata["reference_source_map"]),
    )
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    species_id = metadata["species_map"].get(species, 0)
    family_id = metadata["family_map"].get(family_class, 0)
    refsrc_id = metadata["reference_source_map"].get(reference_source, 0)
    numeric = _encode_numeric(numeric_row, numeric_columns, mean, std)
    results: dict[str, float] = {}
    with torch.no_grad():
        for task_name, task_id in sorted(metadata["task_map"].items(), key=lambda item: item[1]):
            logits = model(
                torch.tensor(numeric, dtype=torch.float32),
                torch.tensor([task_id], dtype=torch.long),
                torch.tensor([species_id], dtype=torch.long),
                torch.tensor([family_id], dtype=torch.long),
                torch.tensor([refsrc_id], dtype=torch.long),
            )
            results[str(task_name)] = float(torch.sigmoid(logits).item())
    return results


def _compute_channel_evidence(
    reference_context: dict[str, Any],
    alt_features: dict[str, Any],
    numeric_row: dict[str, float],
    task_probs: dict[str, float],
    include_biogrid_partners: bool = False,
) -> dict[str, Any]:
    reference = reference_context["reference"]
    gene = reference_context["gene"]
    global_p = task_probs.get("global", 0.5)
    membrane_p = task_probs.get("membrane", global_p)
    surface_p = task_probs.get("surface", membrane_p)
    kinase_p = task_probs.get("kinase", global_p)
    tf_p = task_probs.get("transcription_factor", global_p)

    ppi_reference_weight = math.log1p(
        _parse_int(reference.get("ppi_binding_feature_count", "0"))
        + _parse_int(gene.get("biogrid_partner_count", "0"))
    ) / 5.0
    truncation = _clamp01(1.0 - numeric_row["protein_length_ratio"])
    topology_shift = _clamp01(
        abs(numeric_row["predicted_tm_count_delta"]) / 4.0
        + abs(numeric_row["predicted_signal_score_delta"]) / 3.0
    )
    ppi_impact = _clamp01(0.45 * (1.0 - global_p) + 0.30 * topology_shift + 0.15 * truncation + 0.10 * ppi_reference_weight)

    dna_ref_weight = math.log1p(
        _parse_int(reference.get("dna_binding_feature_count", "0"))
        + _parse_int(reference.get("dna_region_feature_count", "0"))
        + _parse_int(reference.get("zinc_finger_feature_count", "0"))
    ) / 5.0
    pdi_impact = _clamp01(0.55 * (1.0 - tf_p) + 0.25 * truncation + 0.20 * dna_ref_weight)
    localization_shift = _clamp01(0.45 * (1.0 - membrane_p) + 0.35 * topology_shift + 0.20 * (1.0 - surface_p))
    kinase_signal = _clamp01(
        0.65 * kinase_p
        + 0.20 * (1.0 - truncation)
        + 0.15 * _clamp01(math.log1p(_parse_int(reference.get("phospho_feature_count", "0"))) / 4.0)
    )
    dna_binding = _clamp01(
        0.70 * tf_p
        + 0.15 * (1.0 - truncation)
        + 0.15 * _clamp01(math.log1p(_parse_int(reference.get("dna_binding_feature_count", "0"))) / 3.0)
    )

    ref_signal = max(
        _parse_int(reference.get("has_signal_peptide", "0")),
        int(_parse_int(reference.get("signal_feature_count", "0")) > 0),
        _parse_int(reference.get("predicted_signal_candidate", "0")),
    )
    signal_retained = _clamp01(
        0.50 * float(alt_features["predicted_signal_candidate"])
        + 0.20 * (1.0 if ref_signal else 0.5)
        + 0.20 * _clamp01(alt_features["predicted_signal_score"] / 3.0)
        + 0.10 * membrane_p
    )
    tm_insertion = _clamp01(
        0.45 * membrane_p
        + 0.20 * _clamp01(alt_features["predicted_tm_count"] / 3.0)
        + 0.20 * _clamp01(alt_features["predicted_tm_total_span"] / 60.0)
        + 0.15 * _clamp01(alt_features["predicted_tm_max_hydropathy_19"] / 2.5)
    )
    tm_fold = _clamp01(
        0.45 * tm_insertion
        + 0.15 * signal_retained
        + 0.15 * _clamp01(alt_features["glyco_motif_count"] / 3.0)
        + 0.15 * _clamp01(alt_features["cysteine_fraction"] / 0.05)
        + 0.10 * surface_p
    )

    return {
        "ppi_impact_risk_score": ppi_impact,
        "pdi_impact_risk_score": pdi_impact,
        "dpi_impact_risk_score": pdi_impact,
        "localization_shift_risk_score": localization_shift,
        "tf_dna_binding_retained_score": dna_binding if _parse_int(reference.get("is_transcription_factor", "0")) == 1 else None,
        "kinase_signaling_retained_score": kinase_signal if _parse_int(reference.get("is_kinase", "0")) == 1 else None,
        "signal_peptide_retained_score": signal_retained,
        "tm_insertion_support_score": tm_insertion,
        "tm_fold_support_score": tm_fold,
        "biogrid_partner_lookup_performed": bool(include_biogrid_partners),
        "known_biogrid_partners": _load_biogrid_partners(gene.get("gene_name", "")) if include_biogrid_partners else [],
    }


def load_query_protein_sequence(path: Path | None, inline_sequence: str | None) -> str:
    if inline_sequence:
        return inline_sequence.strip().replace(" ", "").replace("\n", "").upper()
    if path is None:
        raise ValueError("Either --alt-protein-seq or --alt-protein-fasta is required.")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:  # type: ignore[arg-type]
        seq: list[str] = []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    break
                continue
            seq.append(line)
    result = "".join(seq).upper()
    if not result:
        raise ValueError(f"No protein sequence found in {path}")
    return result


def predict_query_isoform(
    species: str,
    alt_protein_seq: str,
    alt_transcript_length: int | None,
    gene_id: str | None = None,
    gene_name: str | None = None,
    reference_transcript_id: str | None = None,
    alternative_transcript_id: str | None = None,
    checkpoint_path: Path = DEFAULT_CHECKPOINT,
    include_biogrid_partners: bool = False,
) -> dict[str, Any]:
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    reference_context = select_reference_transcript(
        species=species,
        gene_id=gene_id,
        gene_name=gene_name,
        reference_transcript_id=reference_transcript_id,
    )
    gencode_lengths = _load_gencode_reference_lengths()
    reference_context["reference_lengths"] = gencode_lengths.get(
        (species, reference_context["reference"]["transcript_id"]),
        {"transcript_length": "0"},
    )
    alt_features = derive_query_sequence_features(alt_protein_seq)
    alt_transcript_length_value = alt_transcript_length or max(3 * alt_features["protein_length"], 0)
    family_class = _family_class(reference_context["reference"], alt_features)
    numeric_row = _build_numeric_row(reference_context, alt_features, alt_features["protein_length"], alt_transcript_length_value)
    task_probs = _compute_task_probabilities(
        checkpoint=checkpoint,
        numeric_row=numeric_row,
        species=species,
        family_class=family_class,
        reference_source=reference_context["reference_source"],
    )
    applicability = _task_applicability(reference_context["reference"], alt_features, reference_context["gene"])
    evidence = _compute_channel_evidence(
        reference_context,
        alt_features,
        numeric_row,
        task_probs,
        include_biogrid_partners=include_biogrid_partners,
    )
    return {
        "species": species,
        "gene_id": reference_context["gene"]["gencode_gene_id"],
        "gene_name": reference_context["gene"]["gene_name"],
        "reference_selection": {
            "reference_transcript_id": reference_context["reference"]["transcript_id"],
            "reference_transcript_name": reference_context["reference"].get("transcript_name", ""),
            "reference_source": reference_context["reference_source"],
            "selection_reason": reference_context["selection_reason"],
        },
        "query_isoform": {
            "alternative_transcript_id": alternative_transcript_id or "query_alt",
            "protein_length": alt_features["protein_length"],
            "transcript_length": alt_transcript_length_value,
            "predicted_tm_segments": [{"start_aa": s, "end_aa": e} for s, e in alt_features["tm_segments"]],
        },
        "family_class": family_class,
        "task_applicability": applicability,
        "task_probabilities": task_probs,
        "channel_evidence": evidence,
        "feature_snapshot": {
            "reference_is_membrane": _parse_int(reference_context["reference"].get("is_membrane_protein", "0")),
            "reference_is_kinase": _parse_int(reference_context["reference"].get("is_kinase", "0")),
            "reference_is_transcription_factor": _parse_int(reference_context["reference"].get("is_transcription_factor", "0")),
            "reference_ppi_binding_feature_count": _parse_int(reference_context["reference"].get("ppi_binding_feature_count", "0")),
            "reference_dna_binding_feature_count": _parse_int(reference_context["reference"].get("dna_binding_feature_count", "0")),
            "reference_phospho_feature_count": _parse_int(reference_context["reference"].get("phospho_feature_count", "0")),
            "alternative_predicted_signal_candidate": alt_features["predicted_signal_candidate"],
            "alternative_predicted_signal_score": alt_features["predicted_signal_score"],
            "alternative_predicted_tm_count": alt_features["predicted_tm_count"],
            "alternative_predicted_tm_total_span": alt_features["predicted_tm_total_span"],
            "alternative_glyco_motif_count": alt_features["glyco_motif_count"],
            "alternative_cysteine_fraction": alt_features["cysteine_fraction"],
        },
    }


def format_prediction_json(result: dict[str, Any]) -> str:
    return json.dumps(result, indent=2, sort_keys=False)
