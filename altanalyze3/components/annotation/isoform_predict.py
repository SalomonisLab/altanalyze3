"""Annotate splicing-event tables with higher-confidence isoform associations.

This module combines:
1) long-read significance output containing a `Feature` column with competing junctions
2) ISV summarize TSV output with top cluster structures per gene/event
3) transcript_associations and translated isoform outputs from gff-output

The result is the input event table plus added columns describing the best
supported isoform for each side of the event, associated ORF/protein sequences,
and the minimal protein differences between the two selected isoforms.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd
from Bio import SeqIO

if __package__ in {None, ""}:  # Support direct `python path/to/script.py ...`
    _THIS_FILE = os.path.abspath(__file__)
    _REPO_ROOT = os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(_THIS_FILE)))
    )
    if _REPO_ROOT not in sys.path:
        sys.path.insert(0, _REPO_ROOT)
    from altanalyze3.components.annotation.junction_isoform import (  # type: ignore
        load_protein_summary,
        _is_potential_nmd_status,
    )
else:
    from .junction_isoform import load_protein_summary, _is_potential_nmd_status


_TOP_CLUSTER_FIELD_RE = re.compile(
    r"^top_cluster_(\d+)_(label|structure|transcript_id|structure_n|tx_n)$"
)


@dataclass
class ClusterSummary:
    rank: int
    label: str
    structure: str
    transcript_id: str
    structure_n: str
    tx_n: str


@dataclass
class IsvSummaryRow:
    gene_label: str
    gene_symbol: str
    supplementary_annotation: str
    cell_types: str
    conditions: str
    clusters: List[ClusterSummary]
    raw: Dict[str, str]


@dataclass
class TranscriptAssociation:
    gene_id: str
    strand: str
    structure: str
    isoform_id: str
    sample: str


@dataclass
class IsoformMatch:
    junction: str
    junction_key: str
    summary_gene: str
    cluster_rank: Optional[int]
    cluster_label: str
    cluster_structure: str
    matched_structure: str
    isoform_id: str
    protein_length: int
    nmd_status: str
    sequence_source: str
    orf_sequence: str
    protein_sequence: str


def _split_gene_label(label: str) -> Tuple[str, str]:
    value = str(label or "").strip()
    if "_" not in value:
        return value, ""
    base, suffix = value.rsplit("_", 1)
    if suffix and re.match(r"^[EIU]\d+(?:\.\d+)?[A-Za-z0-9.-]*$", suffix):
        return base, suffix
    return value, ""


def _clean_string(value: object) -> str:
    if value is None:
        return ""
    if pd.isna(value):
        return ""
    return str(value).strip()


def _strip_version(identifier: str) -> str:
    value = _clean_string(identifier)
    if not value:
        return ""
    return value.split(".", 1)[0]


def _normalize_isoform_id(identifier: str) -> str:
    value = _clean_string(identifier)
    if not value:
        return ""
    primary = value.split("|", 1)[0].strip()
    primary = primary.replace("(NMD)", "").replace("(NMD-Potential)", "").strip()
    return primary


def _isoform_lookup_keys(identifier: str) -> List[str]:
    value = _normalize_isoform_id(identifier)
    if not value:
        return []
    keys = [value, _strip_version(value)]
    if "_" in value:
        tail = value.split("_")[-1]
        keys.extend([tail, _strip_version(tail)])
    return [k for i, k in enumerate(keys) if k and k not in keys[:i]]


def _normalize_structure(structure: str) -> str:
    return "|".join(tok.strip() for tok in _clean_string(structure).split("|") if tok.strip())


def _structure_tokens(structure: str) -> List[str]:
    return [tok.strip() for tok in _clean_string(structure).split("|") if tok.strip()]


def _normalize_junction_component(token: str) -> str:
    value = _clean_string(token)
    if ":" in value:
        value = value.split(":")[-1]
    return value


def _normalize_junction_key(junction: str) -> str:
    value = _clean_string(junction)
    if not value:
        return ""
    if ":" in value:
        value = value.split(":", 1)[1]
    parts = value.split("-")
    if len(parts) != 2:
        return value
    return "-".join(_normalize_junction_component(part) for part in parts)


def _junctions_from_structure(structure: str) -> List[str]:
    tokens = _structure_tokens(structure)
    return [
        _normalize_junction_key(f"{tokens[i]}-{tokens[i + 1]}")
        for i in range(len(tokens) - 1)
    ]


def _token_set(structure: str) -> set:
    return set(_structure_tokens(structure))


def _jaccard(set_a: set, set_b: set) -> float:
    if not set_a and not set_b:
        return 1.0
    if not set_a or not set_b:
        return 0.0
    union = set_a | set_b
    if not union:
        return 0.0
    return len(set_a & set_b) / float(len(union))


def _query_coverage(query_set: set, reference_set: set) -> float:
    if not query_set:
        return 0.0
    return len(query_set & reference_set) / float(len(query_set))


def _protein_length_for_isoform(
    gene_id: str,
    isoform_id: str,
    protein_summary: Dict[Tuple[str, str], Tuple[int, str, str]],
) -> Tuple[int, str]:
    for key in _isoform_lookup_keys(isoform_id):
        length, nmd_status, _longest = protein_summary.get((gene_id, key), (0, "", ""))
        try:
            return int(length), _clean_string(nmd_status)
        except Exception:
            continue
    return 0, ""


def _is_ensembl_transcript(identifier: str) -> bool:
    return _normalize_isoform_id(identifier).startswith("ENST")


def _extract_record_lookup_keys(record) -> List[str]:
    description = _clean_string(record.description)
    candidates: List[str] = []
    candidates.extend(
        re.findall(r"(ENST\d+(?:\.\d+)?)", description)
    )
    candidates.extend(
        re.findall(r"(ENSP\d+(?:\.\d+)?)", description)
    )
    for tag in ("transcript_id:", "transcript:", "isoform:", "protein_id:"):
        if tag in description:
            fragment = description.split(tag, 1)[1].split()[0].strip(";")
            candidates.append(fragment)
    primary = description.split(None, 1)[0].strip(";")
    if primary:
        candidates.append(primary)
        if "_" in primary:
            candidates.append(primary.split("_")[-1])
    keys: List[str] = []
    for candidate in candidates:
        for key in _isoform_lookup_keys(candidate):
            if key not in keys:
                keys.append(key)
    return keys


def load_sequence_fasta(file_path: str) -> Dict[str, str]:
    lookup: Dict[str, str] = {}
    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq).strip()
            if not sequence:
                continue
            for key in _extract_record_lookup_keys(record):
                lookup.setdefault(key, sequence)
    return lookup


def load_ensembl_protein_tsv(file_path: Optional[str]) -> Dict[str, str]:
    if not file_path:
        return {}
    df = pd.read_csv(file_path, sep="\t", dtype=str).fillna("")
    tx_col = next(
        (col for col in df.columns if col in {"ensembl_transcript_id", "transcript_id", "ENST"}),
        None,
    )
    seq_col = next(
        (col for col in df.columns if col in {"peptide", "protein_sequence", "sequence"}),
        None,
    )
    if tx_col is None or seq_col is None:
        raise ValueError(
            "Ensembl protein TSV must contain transcript and sequence columns, e.g. "
            "'ensembl_transcript_id' and 'peptide'."
        )
    lookup: Dict[str, str] = {}
    for _, row in df.iterrows():
        tx_id = _strip_version(row.get(tx_col, ""))
        sequence = _clean_string(row.get(seq_col, ""))
        if tx_id and sequence:
            lookup[tx_id] = sequence
    return lookup


def load_transcript_associations_by_gene(file_path: str) -> Dict[str, List[TranscriptAssociation]]:
    by_gene: Dict[str, List[TranscriptAssociation]] = defaultdict(list)
    with open(file_path, "r") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            gene_id, strand, structure, isoform_id, sample = parts[:5]
            assoc = TranscriptAssociation(
                gene_id=_clean_string(gene_id),
                strand=_clean_string(strand),
                structure=_normalize_structure(structure),
                isoform_id=_normalize_isoform_id(isoform_id),
                sample=_clean_string(sample),
            )
            if assoc.gene_id and assoc.structure and assoc.isoform_id:
                by_gene[assoc.gene_id].append(assoc)
    return by_gene


def load_isv_summary(file_path: str) -> List[IsvSummaryRow]:
    df = pd.read_csv(file_path, sep="\t", dtype=str).fillna("")
    rows: List[IsvSummaryRow] = []
    for _, row in df.iterrows():
        gene_label = _clean_string(row.get("gene") or row.get("gene_symbol"))
        gene_symbol = _clean_string(row.get("gene_symbol"))
        if not gene_symbol:
            gene_symbol, _suffix = _split_gene_label(gene_label)
        supplementary_annotation = _clean_string(row.get("supplementary_annotation"))
        if not supplementary_annotation:
            _base, supplementary_annotation = _split_gene_label(gene_label)
        clusters_by_rank: Dict[int, Dict[str, str]] = defaultdict(dict)
        for col in df.columns:
            match = _TOP_CLUSTER_FIELD_RE.match(col)
            if not match:
                continue
            rank = int(match.group(1))
            field = match.group(2)
            clusters_by_rank[rank][field] = _clean_string(row.get(col))
        clusters = [
            ClusterSummary(
                rank=rank,
                label=fields.get("label", ""),
                structure=_normalize_structure(fields.get("structure", "")),
                transcript_id=_normalize_isoform_id(fields.get("transcript_id", "")),
                structure_n=fields.get("structure_n", ""),
                tx_n=fields.get("tx_n", ""),
            )
            for rank, fields in sorted(clusters_by_rank.items())
            if _normalize_structure(fields.get("structure", ""))
        ]
        rows.append(
            IsvSummaryRow(
                gene_label=gene_label,
                gene_symbol=gene_symbol,
                supplementary_annotation=supplementary_annotation,
                cell_types=_clean_string(row.get("cell_types")),
                conditions=_clean_string(row.get("conditions")),
                clusters=clusters,
                raw={col: _clean_string(row.get(col)) for col in df.columns},
            )
        )
    return rows


def _feature_junctions(feature: str) -> Tuple[str, str]:
    parts = _clean_string(feature).split("|")
    if len(parts) != 2:
        raise ValueError(f"Expected Feature to contain two competing junctions separated by '|': {feature}")
    return parts[0].strip(), parts[1].strip()


def _gene_id_from_junction(junction: str) -> str:
    return _clean_string(junction).split(":", 1)[0]


def _event_gene_labels(row: pd.Series) -> List[str]:
    labels: List[str] = []
    for key in ("gene", "Gene_Symbol", "gene_symbol", "Symbol"):
        value = _clean_string(row.get(key))
        if value and value not in labels:
            labels.append(value)
    return labels


def _summary_row_overlap_score(summary: IsvSummaryRow, junction_keys: Sequence[str]) -> Tuple[int, int]:
    score = 0
    best_rank = 999999
    for cluster in summary.clusters:
        cluster_junctions = set(_junctions_from_structure(cluster.structure))
        local_hits = sum(1 for junction in junction_keys if junction in cluster_junctions)
        if local_hits:
            score += local_hits
            best_rank = min(best_rank, cluster.rank)
    return score, -best_rank


def choose_isv_summary_row(event_row: pd.Series, isv_rows: Sequence[IsvSummaryRow]) -> Optional[IsvSummaryRow]:
    feature = _clean_string(event_row.get("Feature"))
    if not feature:
        return None
    gene_labels = _event_gene_labels(event_row)
    base_labels = []
    for label in gene_labels:
        base, _suffix = _split_gene_label(label)
        if base and base not in base_labels:
            base_labels.append(base)

    candidates: List[IsvSummaryRow] = []
    for summary in isv_rows:
        if summary.gene_label in gene_labels or summary.gene_symbol in gene_labels:
            candidates.append(summary)
    if not candidates:
        for summary in isv_rows:
            if summary.gene_symbol in base_labels or summary.gene_label in base_labels:
                candidates.append(summary)
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]

    junction_keys = [_normalize_junction_key(part) for part in _feature_junctions(feature)]
    ranked = sorted(
        candidates,
        key=lambda summary: (
            -_summary_row_overlap_score(summary, junction_keys)[0],
            -_summary_row_overlap_score(summary, junction_keys)[1],
            summary.gene_label,
        ),
    )
    return ranked[0]


def _structure_match_rank(
    assoc: TranscriptAssociation,
    cluster: ClusterSummary,
    junction_key: str,
) -> Tuple:
    assoc_structure = _normalize_structure(assoc.structure)
    cluster_structure = _normalize_structure(cluster.structure)
    assoc_junctions = set(_junctions_from_structure(assoc_structure))
    cluster_junctions = set(_junctions_from_structure(cluster_structure))
    assoc_tokens = _token_set(assoc_structure)
    cluster_tokens = _token_set(cluster_structure)
    exact_structure = int(assoc_structure != cluster_structure)
    transcript_match = int(
        _normalize_isoform_id(assoc.isoform_id) != _normalize_isoform_id(cluster.transcript_id)
    )
    junction_hit = int(junction_key not in assoc_junctions)
    junction_query_coverage = _query_coverage(cluster_junctions, assoc_junctions)
    token_query_coverage = _query_coverage(cluster_tokens, assoc_tokens)
    return (
        junction_hit,
        cluster.rank,
        -junction_query_coverage,
        -token_query_coverage,
        exact_structure,
        transcript_match,
        -len(assoc_tokens),
        assoc.isoform_id,
    )


def _lookup_sequence(sequence_db: Dict[str, str], isoform_id: str) -> str:
    for key in _isoform_lookup_keys(isoform_id):
        if key in sequence_db:
            return sequence_db[key]
    return ""


def _choose_junction_candidate_pool(
    candidates: Sequence[TranscriptAssociation],
) -> Tuple[List[TranscriptAssociation], str]:
    ensembl = [assoc for assoc in candidates if _is_ensembl_transcript(assoc.isoform_id)]
    if ensembl:
        return ensembl, "ensembl_junction_match"
    return list(candidates), "non_ensembl_junction_match"


def _candidate_selection_priority(
    assoc: TranscriptAssociation,
    gene_id: str,
    protein_summary: Dict[Tuple[str, str], Tuple[int, str, str]],
    orf_sequences: Dict[str, str],
    protein_sequences: Dict[str, str],
) -> Tuple:
    protein_seq = _lookup_sequence(protein_sequences, assoc.isoform_id)
    orf_seq = _lookup_sequence(orf_sequences, assoc.isoform_id)
    length, nmd_status = _resolved_protein_length(
        gene_id,
        assoc.isoform_id,
        protein_summary,
        protein_seq,
        orf_seq,
    )
    has_coding = int(not (length > 0))
    is_nmd = int(_is_potential_nmd_status(nmd_status))
    return (has_coding, is_nmd, -length, assoc.isoform_id)


def _resolved_protein_length(
    gene_id: str,
    isoform_id: str,
    protein_summary: Dict[Tuple[str, str], Tuple[int, str, str]],
    protein_seq: str = "",
    orf_seq: str = "",
) -> Tuple[int, str]:
    length, nmd_status = _protein_length_for_isoform(gene_id, isoform_id, protein_summary)
    if length > 0:
        return length, nmd_status
    if protein_seq:
        return len(protein_seq), nmd_status
    if orf_seq:
        return max(0, len(orf_seq) // 3), nmd_status
    return 0, nmd_status


def _best_cluster_for_junction(
    assoc: TranscriptAssociation,
    relevant_clusters: Sequence[ClusterSummary],
    junction_key: str,
) -> Tuple[Tuple, Optional[ClusterSummary]]:
    if not relevant_clusters:
        return (
            (
                0,
                999999,
                0.0,
                0.0,
                1,
                1,
                -len(_token_set(assoc.structure)),
                assoc.isoform_id,
            ),
            None,
        )
    ranked = [
        (_structure_match_rank(assoc, cluster, junction_key), cluster)
        for cluster in relevant_clusters
    ]
    return min(ranked, key=lambda item: item[0])


def _choose_best_transcript_match(
    gene_id: str,
    junction: str,
    summary: Optional[IsvSummaryRow],
    transcript_db: Dict[str, List[TranscriptAssociation]],
    protein_summary: Dict[Tuple[str, str], Tuple[int, str, str]],
    orf_sequences: Dict[str, str],
    protein_sequences: Dict[str, str],
    ensembl_protein_sequences: Dict[str, str],
) -> Optional[IsoformMatch]:
    associations = transcript_db.get(gene_id, [])
    if not associations:
        return None

    junction_key = _normalize_junction_key(junction)
    candidates = [
        assoc for assoc in associations
        if junction_key in set(_junctions_from_structure(assoc.structure))
    ]
    if not candidates:
        return None

    candidates, junction_match_source = _choose_junction_candidate_pool(candidates)
    coding_candidates = [
        assoc for assoc in candidates
        if _candidate_selection_priority(
            assoc, gene_id, protein_summary, orf_sequences, protein_sequences
        )[0] == 0
    ]
    if coding_candidates:
        candidates = coding_candidates

    matching_clusters = []
    if summary is not None:
        matching_clusters = [
            cluster for cluster in summary.clusters
            if junction_key in set(_junctions_from_structure(cluster.structure))
        ]
        if not matching_clusters:
            matching_clusters = list(summary.clusters)

    scored_candidates: List[Tuple[Tuple, TranscriptAssociation, Optional[ClusterSummary]]] = []
    for assoc in candidates:
        best_rank, best_cluster = _best_cluster_for_junction(assoc, matching_clusters, junction_key)
        scored_candidates.append((best_rank, assoc, best_cluster))

    if not scored_candidates:
        return None

    best_structure_rank = min(rank for rank, _assoc, _cluster in scored_candidates)
    top_candidates = [
        (assoc, cluster) for rank, assoc, cluster in scored_candidates if rank == best_structure_rank
    ]

    best_assoc, best_cluster = sorted(
        top_candidates,
        key=lambda item: _candidate_selection_priority(
            item[0], gene_id, protein_summary, orf_sequences, protein_sequences
        ),
    )[0]
    protein_sequence = _lookup_sequence(protein_sequences, best_assoc.isoform_id)
    orf_sequence = _lookup_sequence(orf_sequences, best_assoc.isoform_id)
    protein_length, nmd_status = _resolved_protein_length(
        gene_id,
        best_assoc.isoform_id,
        protein_summary,
        protein_sequence,
        orf_sequence,
    )
    sequence_source = (
        f"{junction_match_source}|protein_sequences.fasta"
        if protein_sequence else junction_match_source
    )
    if _is_ensembl_transcript(best_assoc.isoform_id):
        ensembl_seq = _lookup_sequence(ensembl_protein_sequences, best_assoc.isoform_id)
        if ensembl_seq:
            protein_sequence = ensembl_seq
            sequence_source = f"{junction_match_source}|ensembl_protein_tsv"
            protein_length = len(protein_sequence)
    return IsoformMatch(
        junction=junction,
        junction_key=junction_key,
        summary_gene=summary.gene_label if summary is not None else "",
        cluster_rank=best_cluster.rank if best_cluster is not None else None,
        cluster_label=best_cluster.label if best_cluster is not None else "",
        cluster_structure=best_cluster.structure if best_cluster is not None else "",
        matched_structure=best_assoc.structure,
        isoform_id=best_assoc.isoform_id,
        protein_length=protein_length,
        nmd_status=nmd_status,
        sequence_source=sequence_source,
        orf_sequence=orf_sequence,
        protein_sequence=protein_sequence,
    )


def _sliding_window(sequence: str, window_size: int) -> Dict[str, int]:
    return {
        sequence[i:i + window_size]: i
        for i in range(max(0, len(sequence) - window_size + 1))
    }


def _refine_unique_sequences(seq1: str, seq2: str, unique_matches: Dict[str, int]) -> List[Tuple[str, str]]:
    refined_sequences = []
    sorted_matches = sorted(unique_matches.items(), key=lambda x: x[1])
    merged_sequence = None
    for peptide, position in sorted_matches:
        if merged_sequence is None:
            merged_sequence = [peptide, position, position + len(peptide)]
        elif position <= merged_sequence[2]:
            overlap_start = max(merged_sequence[2] - position, 0)
            merged_sequence[0] += peptide[overlap_start:]
            merged_sequence[2] = position + len(peptide)
        else:
            refined_sequences.append(merged_sequence)
            merged_sequence = [peptide, position, position + len(peptide)]
    if merged_sequence:
        refined_sequences.append(merged_sequence)

    final_sequences = []
    for merged_peptide, start, end in refined_sequences:
        while start < len(seq1) and start < len(seq2) and seq1[start] == seq2[start]:
            start += 1
        while end > start and end <= len(seq1) and end <= len(seq2) and seq1[end - 1] == seq2[end - 1]:
            end -= 1
        refined_peptide = seq1[start:end]
        flanking_region = seq1[max(0, start - 10):min(len(seq1), end + 10)]
        if refined_peptide:
            final_sequences.append((refined_peptide, flanking_region))
    return final_sequences


def _find_unique_sequences_with_refinement(seq1: str, seq2: str, window_size: int = 5) -> List[Tuple[str, str]]:
    if not seq1 or not seq2 or len(seq1) < window_size:
        return []
    peptides1 = _sliding_window(seq1, window_size)
    peptides2 = _sliding_window(seq2, window_size)
    unique_to_seq1 = {
        peptide: position
        for peptide, position in peptides1.items()
        if peptide not in peptides2
    }
    return _refine_unique_sequences(seq1, seq2, unique_to_seq1)


def _minimal_difference(seq1: str, seq2: str) -> Tuple[str, str]:
    unique = _find_unique_sequences_with_refinement(seq1, seq2)
    if not unique:
        return "", ""
    best = min(unique, key=lambda item: (len(item[0]), item[0]))
    return best


def annotate_high_confidence_isoforms(
    event_tsv: str,
    isv_summary_tsv: str,
    transcript_associations_tsv: str,
    protein_summary_tsv: str,
    protein_sequences_fasta: str,
    orf_sequences_fasta: str,
    output_tsv: Optional[str] = None,
    ensembl_protein_tsv: Optional[str] = None,
) -> pd.DataFrame:
    print(f"[isoform_predict] Loading event table: {event_tsv}")
    events_df = pd.read_csv(event_tsv, sep="\t", dtype=str).fillna("")
    print(f"[isoform_predict] Loaded {len(events_df)} event rows")

    print(f"[isoform_predict] Loading ISV summary: {isv_summary_tsv}")
    isv_rows = load_isv_summary(isv_summary_tsv)
    print(f"[isoform_predict] Loaded {len(isv_rows)} ISV summary rows")

    print(f"[isoform_predict] Loading transcript associations: {transcript_associations_tsv}")
    transcript_db = load_transcript_associations_by_gene(transcript_associations_tsv)
    print(f"[isoform_predict] Loaded transcript associations for {len(transcript_db)} genes")

    print(f"[isoform_predict] Loading protein summary: {protein_summary_tsv}")
    protein_summary = load_protein_summary(protein_summary_tsv)
    print(f"[isoform_predict] Loaded {len(protein_summary)} protein summary entries")

    print(f"[isoform_predict] Loading protein sequences: {protein_sequences_fasta}")
    protein_sequences = load_sequence_fasta(protein_sequences_fasta)
    print(f"[isoform_predict] Loaded {len(protein_sequences)} protein sequence identifiers")

    print(f"[isoform_predict] Loading ORF sequences: {orf_sequences_fasta}")
    orf_sequences = load_sequence_fasta(orf_sequences_fasta)
    print(f"[isoform_predict] Loaded {len(orf_sequences)} ORF sequence identifiers")

    if ensembl_protein_tsv:
        print(f"[isoform_predict] Loading Ensembl protein TSV: {ensembl_protein_tsv}")
    ensembl_protein_sequences = load_ensembl_protein_tsv(ensembl_protein_tsv)
    if ensembl_protein_tsv:
        print(f"[isoform_predict] Loaded {len(ensembl_protein_sequences)} Ensembl protein sequences")

    annotated_rows: List[Dict[str, str]] = []
    print("[isoform_predict] Annotating events")
    progress_every = 25 if len(events_df) <= 250 else 100
    for idx, (_, row) in enumerate(events_df.iterrows(), start=1):
        feature = _clean_string(row.get("Feature"))
        row_out = row.to_dict()
        if not feature:
            annotated_rows.append(row_out)
            if idx % progress_every == 0 or idx == len(events_df):
                print(f"[isoform_predict] Processed {idx}/{len(events_df)} rows")
            continue

        junction1, junction2 = _feature_junctions(feature)
        summary_row = choose_isv_summary_row(row, isv_rows)
        gene_id_1 = _gene_id_from_junction(junction1)
        gene_id_2 = _gene_id_from_junction(junction2)
        match1 = _choose_best_transcript_match(
            gene_id_1,
            junction1,
            summary_row,
            transcript_db,
            protein_summary,
            orf_sequences,
            protein_sequences,
            ensembl_protein_sequences,
        )
        match2 = _choose_best_transcript_match(
            gene_id_2,
            junction2,
            summary_row,
            transcript_db,
            protein_summary,
            orf_sequences,
            protein_sequences,
            ensembl_protein_sequences,
        )

        protein1 = match1.protein_sequence if match1 else ""
        protein2 = match2.protein_sequence if match2 else ""
        diff1, flank1 = _minimal_difference(protein1, protein2)
        diff2, flank2 = _minimal_difference(protein2, protein1)

        def _populate(prefix: str, match: Optional[IsoformMatch]) -> None:
            row_out[f"{prefix}_junction"] = match.junction if match else ""
            row_out[f"{prefix}_junction_key"] = match.junction_key if match else ""
            row_out[f"{prefix}_summary_gene"] = match.summary_gene if match else ""
            row_out[f"{prefix}_cluster_rank"] = "" if match is None or match.cluster_rank is None else str(match.cluster_rank)
            row_out[f"{prefix}_cluster_label"] = match.cluster_label if match else ""
            row_out[f"{prefix}_cluster_structure"] = match.cluster_structure if match else ""
            row_out[f"{prefix}_matched_structure"] = match.matched_structure if match else ""
            row_out[f"{prefix}_isoform_id"] = match.isoform_id if match else ""
            row_out[f"{prefix}_protein_length"] = "" if match is None else str(match.protein_length)
            row_out[f"{prefix}_nmd_status"] = match.nmd_status if match else ""
            row_out[f"{prefix}_sequence_source"] = match.sequence_source if match else ""
            row_out[f"{prefix}_orf_sequence"] = match.orf_sequence if match else ""
            row_out[f"{prefix}_protein_sequence"] = match.protein_sequence if match else ""

        _populate("high_confidence_isoform_1", match1)
        _populate("high_confidence_isoform_2", match2)
        row_out["minimal_protein_difference_1"] = diff1
        row_out["minimal_protein_difference_1_flank"] = flank1
        row_out["minimal_protein_difference_2"] = diff2
        row_out["minimal_protein_difference_2_flank"] = flank2
        row_out["isv_summary_gene"] = summary_row.gene_label if summary_row else ""
        row_out["isv_summary_cell_types"] = summary_row.cell_types if summary_row else ""
        row_out["isv_summary_conditions"] = summary_row.conditions if summary_row else ""
        annotated_rows.append(row_out)
        if idx % progress_every == 0 or idx == len(events_df):
            print(f"[isoform_predict] Processed {idx}/{len(events_df)} rows")

    result_df = pd.DataFrame(annotated_rows)
    if output_tsv is None:
        base, ext = os.path.splitext(event_tsv)
        output_tsv = f"{base}-high_confidence{ext or '.tsv'}"
    print(f"[isoform_predict] Writing output: {output_tsv}")
    result_df.to_csv(output_tsv, sep="\t", index=False)
    print(f"[isoform_predict] Completed: wrote {len(result_df)} rows")
    return result_df


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Annotate a splicing-event TSV with higher-confidence isoform "
            "associations using ISV summary output and gff-output translation files."
        )
    )
    parser.add_argument("event_tsv", help="Input splicing-event TSV with a Feature column.")
    parser.add_argument("isv_summary_tsv", help="ISV summarize TSV output.")
    parser.add_argument("transcript_associations_tsv", help="gff-output/transcript_associations.txt")
    parser.add_argument("protein_summary_tsv", help="protein_summary.txt")
    parser.add_argument("protein_sequences_fasta", help="protein_sequences.fasta")
    parser.add_argument("orf_sequences_fasta", help="orf_sequences.fasta")
    parser.add_argument(
        "--ensembl-protein-tsv",
        default=None,
        help=(
            "Optional Ensembl transcript-to-protein sequence TSV containing "
            "columns such as ensembl_transcript_id and peptide."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output TSV path. Defaults to <event_tsv>-high_confidence.tsv",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    annotate_high_confidence_isoforms(
        event_tsv=args.event_tsv,
        isv_summary_tsv=args.isv_summary_tsv,
        transcript_associations_tsv=args.transcript_associations_tsv,
        protein_summary_tsv=args.protein_summary_tsv,
        protein_sequences_fasta=args.protein_sequences_fasta,
        orf_sequences_fasta=args.orf_sequences_fasta,
        output_tsv=args.output,
        ensembl_protein_tsv=args.ensembl_protein_tsv,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
