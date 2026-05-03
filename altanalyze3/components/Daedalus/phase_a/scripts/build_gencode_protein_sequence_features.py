#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path
import os


RAW = Path(__file__).resolve().parents[1] / "data" / "raw" / "gencode"
OUT = Path(__file__).resolve().parents[1] / "data" / "interim" / "gencode_protein_sequence_features.tsv"

RESOURCES = [
    ("human", RAW / "gencode.v49.pc_translations.fa.gz"),
    ("mouse", RAW / "gencode.vM38.pc_translations.fa.gz"),
]

HYDRO = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4, "H": -3.2, "I": 4.5,
    "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5, "S": -0.8,
    "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}
SMALL_AA = set("AGSCTV")
GLYCO_RE = re.compile(r"N[^P][ST]")


def _hydropathy(window: str) -> float:
    if not window:
        return 0.0
    return sum(HYDRO.get(aa, 0.0) for aa in window) / len(window)


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
    best_hydro = max((_window_hydropathy(prefix, i, i + 8) for i in range(2, max(len(nterm) - 7, 3))), default=0.0)
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


def _parse_fasta(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header: str | None = None
        seq: list[str] = []
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header is not None:
            yield header, "".join(seq)


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    temp_out = OUT.with_suffix(".tmp.tsv")
    fieldnames = [
        "species",
        "protein_id",
        "transcript_id",
        "gene_id",
        "protein_length",
        "predicted_tm_count",
        "predicted_tm_total_span",
        "predicted_tm_max_hydropathy_19",
        "predicted_signal_candidate",
        "predicted_signal_score",
        "n_terminal_hydropathy_max8",
        "glyco_motif_count",
        "cysteine_count",
        "cysteine_fraction",
    ]
    with temp_out.open("w", encoding="utf-8", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        rows = 0
        for species, path in RESOURCES:
            for header, seq in _parse_fasta(path):
                parts = header.split("|")
                if len(parts) < 3:
                    continue
                protein_id, transcript_id, gene_id = parts[:3]
                tm_segments = _predict_tm_segments(seq)
                signal_candidate, signal_score = _predict_signal(seq)
                seq_prefix = _prefix_hydropathy(seq)
                nterm = seq[:35]
                nterm_prefix = _prefix_hydropathy(nterm)
                writer.writerow(
                    {
                        "species": species,
                        "protein_id": protein_id,
                        "transcript_id": transcript_id,
                        "gene_id": gene_id,
                        "protein_length": len(seq),
                        "predicted_tm_count": len(tm_segments),
                        "predicted_tm_total_span": sum(end - start + 1 for start, end in tm_segments),
                        "predicted_tm_max_hydropathy_19": max((_window_hydropathy(seq_prefix, i, i + 19) for i in range(0, max(len(seq) - 18, 1))), default=0.0),
                        "predicted_signal_candidate": signal_candidate,
                        "predicted_signal_score": f"{signal_score:.3f}",
                        "n_terminal_hydropathy_max8": f"{max((_window_hydropathy(nterm_prefix, i, i + 8) for i in range(0, max(len(nterm) - 7, 1))), default=0.0):.3f}",
                        "glyco_motif_count": len(GLYCO_RE.findall(seq)),
                        "cysteine_count": seq.count("C"),
                        "cysteine_fraction": f"{(seq.count('C') / len(seq)):.6f}" if seq else "0.000000",
                    }
                )
                rows += 1
    os.replace(temp_out, OUT)
    print(f"[ok] Wrote {OUT} rows={rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
