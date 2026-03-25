"""Cross-annotate splicing events against other event result files."""

from __future__ import annotations

import argparse
import csv
import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd


@dataclass
class EventFileInfo:
    path: str
    format_name: str
    feature_column: str
    effect_column: str
    comparison_label: str
    folder_label: str


@dataclass
class EventOccurrence:
    comparison_label: str
    folder_label: str
    path: str
    format_name: str
    feature: str
    left_junction: str
    right_junction: str
    effect_sign: int


def _clean_string(value: object) -> str:
    if value is None:
        return ""
    if pd.isna(value):
        return ""
    return str(value).strip()


def _detect_effect_column(header: Sequence[str]) -> str:
    cleaned = [_clean_string(col) for col in header]
    if "dPSI" in cleaned:
        return "dPSI"
    diffmeans = [col for col in cleaned if col.startswith("DiffMeans_")]
    if diffmeans:
        return diffmeans[0]
    delta = [col for col in cleaned if col.lower() in {"dpsi", "delta_psi", "deltapsi", "diffmeans"}]
    if delta:
        return delta[0]
    return ""


def _detect_event_file_format(path: str) -> Optional[Tuple[str, str, str]]:
    try:
        with open(path, "r", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
    except Exception:
        return None
    header_set = {_clean_string(col) for col in header}
    effect_column = _detect_effect_column(header)
    if "Feature" in header_set:
        return "altanalyze3", "Feature", effect_column
    if "UID" in header_set:
        return "altanalyze2", "UID", effect_column
    return None


def _infer_comparison_label(path: str, format_name: str) -> str:
    base = os.path.splitext(os.path.basename(path))[0]
    if format_name == "altanalyze2":
        if base.startswith("PSI."):
            base = base[4:]
        if "_vs_" in base:
            base = base.split("_vs_", 1)[0]
        return base
    if format_name == "altanalyze3":
        for suffix in ("_significant-results", "-significant-results"):
            if base.endswith(suffix):
                return base[: -len(suffix)]
        return base
    return base


def _infer_folder_label(path: str, root_abs: str) -> str:
    parent_dir = os.path.basename(os.path.dirname(path))
    if parent_dir:
        return parent_dir
    fallback = os.path.basename(root_abs.rstrip(os.sep))
    return fallback or root_abs


def _column_safe_label(label: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", _clean_string(label)).strip("_")
    return cleaned or "unknown"


def _iter_event_files(folder_roots: Sequence[str]) -> List[EventFileInfo]:
    infos: List[EventFileInfo] = []
    for root in folder_roots:
        root_abs = os.path.abspath(root)
        for dirpath, _dirnames, filenames in os.walk(root_abs):
            for filename in sorted(filenames):
                if not filename.endswith((".tsv", ".txt")):
                    continue
                path = os.path.join(dirpath, filename)
                detected = _detect_event_file_format(path)
                if not detected:
                    continue
                format_name, feature_column, effect_column = detected
                infos.append(
                    EventFileInfo(
                        path=os.path.abspath(path),
                        format_name=format_name,
                        feature_column=feature_column,
                        effect_column=effect_column,
                        comparison_label=_infer_comparison_label(path, format_name),
                        folder_label=_infer_folder_label(path, root_abs),
                    )
                )
    return infos


def _feature_from_row(row: pd.Series, feature_column: str) -> str:
    return _clean_string(row.get(feature_column))


def _extract_junction_feature_text(feature: str) -> str:
    text = _clean_string(feature)
    if not text:
        return ""
    if "|" not in text:
        return text
    ensg_match = re.search(r"(ENS[A-Z0-9]*G\d+:[^\t]+?\|[^\t]+)$", text)
    if ensg_match:
        return ensg_match.group(1).strip()
    if ": " in text:
        _prefix, remainder = text.split(": ", 1)
        if "|" in remainder:
            return remainder.strip()
    if ":" in text and not text.startswith("ENS"):
        _prefix, remainder = text.split(":", 1)
        if remainder.strip().startswith("ENS") and "|" in remainder:
            return remainder.strip()
    return text


def _split_feature(feature: str) -> Tuple[str, str]:
    normalized = _extract_junction_feature_text(feature)
    parts = normalized.split("|")
    if len(parts) != 2:
        raise ValueError(f"Expected exactly two competing junctions in feature: {feature}")
    return parts[0].strip(), parts[1].strip()


def _canonical_feature_key(feature: str) -> str:
    left, right = _split_feature(feature)
    ordered = sorted([left, right])
    return "|".join(ordered)


def _reversed_feature_key(feature: str) -> str:
    left, right = _split_feature(feature)
    return f"{right}|{left}"


def _effect_sign_from_row(row: pd.Series, effect_column: str) -> int:
    if not effect_column:
        return 0
    raw = _clean_string(row.get(effect_column))
    if not raw:
        return 0
    try:
        value = float(raw)
    except ValueError:
        return 0
    if value > 0:
        return 1
    if value < 0:
        return -1
    return 0


def _junction_match_relations(
    query_left: str,
    query_right: str,
    ref_left: str,
    ref_right: str,
) -> List[Tuple[str, bool]]:
    relations: List[Tuple[str, bool]] = []
    if query_left == ref_left:
        relations.append(("q1_to_r1", True))
    if query_left == ref_right:
        relations.append(("q1_to_r2", False))
    if query_right == ref_left:
        relations.append(("q2_to_r1", False))
    if query_right == ref_right:
        relations.append(("q2_to_r2", True))
    return relations


def _classify_relation_status(query_sign: int, ref_sign: int, same_position: bool) -> str:
    if query_sign == 0 or ref_sign == 0:
        return "unknown"
    if same_position:
        return "match" if query_sign == ref_sign else "opposite"
    return "match" if query_sign != ref_sign else "opposite"


def _relation_priority(relation_name: str) -> int:
    priorities = {
        "q1_to_r1": 4,
        "q2_to_r2": 4,
        "q1_to_r2": 3,
        "q2_to_r1": 3,
    }
    return priorities.get(relation_name, 0)


def _status_priority(status: str) -> int:
    if status == "match":
        return 3
    if status == "opposite":
        return 2
    return 1


def _index_comparison_events(
    comparison_files: Sequence[EventFileInfo],
    input_path: str,
) -> Dict[str, List[EventOccurrence]]:
    junction_index: Dict[str, List[EventOccurrence]] = {}
    input_abs = os.path.abspath(input_path)
    for info in comparison_files:
        if info.path == input_abs:
            continue
        print(f"[event_cross_comparison] Indexing {info.comparison_label}: {info.path}")
        df = pd.read_csv(info.path, sep="\t", dtype=str).fillna("")
        for _, row in df.iterrows():
            feature = _feature_from_row(row, info.feature_column)
            if not feature:
                continue
            try:
                left, right = _split_feature(feature)
            except ValueError:
                continue
            occurrence = EventOccurrence(
                comparison_label=info.comparison_label,
                folder_label=info.folder_label,
                path=info.path,
                format_name=info.format_name,
                feature=_extract_junction_feature_text(feature),
                left_junction=left,
                right_junction=right,
                effect_sign=_effect_sign_from_row(row, info.effect_column),
            )
            for junction in {left, right}:
                junction_index.setdefault(junction, []).append(occurrence)
    return junction_index


def annotate_event_cross_comparison(
    input_tsv: str,
    comparison_roots: Sequence[str],
    output_tsv: Optional[str] = None,
) -> pd.DataFrame:
    input_detected = _detect_event_file_format(input_tsv)
    if input_detected is None:
        raise ValueError(
            "Input TSV format not recognized. Expected AltAnalyze3 with 'Feature' or AltAnalyze2 with 'UID'."
        )
    input_format_name, input_feature_column, input_effect_column = input_detected
    print(f"[event_cross_comparison] Input format: {input_format_name} ({input_feature_column})")
    print(f"[event_cross_comparison] Loading input: {input_tsv}")
    input_df = pd.read_csv(input_tsv, sep="\t", dtype=str).fillna("")
    print(f"[event_cross_comparison] Loaded {len(input_df)} input rows")

    comparison_files = _iter_event_files(comparison_roots)
    print(f"[event_cross_comparison] Found {len(comparison_files)} comparison files")
    junction_index = _index_comparison_events(comparison_files, input_tsv)
    print(f"[event_cross_comparison] Indexed {len(junction_index)} unique junctions")
    folder_labels = sorted({info.folder_label for info in comparison_files})
    folder_columns = {
        label: (
            f"cross_comparison_matches_{_column_safe_label(label)}",
            f"cross_comparison_match_count_{_column_safe_label(label)}",
        )
        for label in folder_labels
    }

    annotated_rows: List[Dict[str, str]] = []
    progress_every = 25 if len(input_df) <= 250 else 100
    for idx, (_, row) in enumerate(input_df.iterrows(), start=1):
        row_out = row.to_dict()
        feature = _feature_from_row(row, input_feature_column)
        if not feature:
            annotated_rows.append(row_out)
            continue
        try:
            canonical_key = _canonical_feature_key(feature)
            left, right = _split_feature(feature)
        except ValueError:
            row_out["cross_comparison_feature_key"] = ""
            for matches_column, count_column in folder_columns.values():
                row_out[matches_column] = ""
                row_out[count_column] = "0"
            annotated_rows.append(row_out)
            continue

        query_sign = _effect_sign_from_row(row, input_effect_column)
        candidate_occurrences: List[EventOccurrence] = []
        seen_occurrences: set[Tuple[str, str, str]] = set()
        for junction in {left, right}:
            for occ in junction_index.get(junction, []):
                key = (occ.path, occ.feature, occ.comparison_label)
                if key in seen_occurrences:
                    continue
                seen_occurrences.add(key)
                candidate_occurrences.append(occ)

        best_by_folder_dataset: Dict[Tuple[str, str], Tuple[str, str, str]] = {}
        for occ in candidate_occurrences:
            relations = _junction_match_relations(left, right, occ.left_junction, occ.right_junction)
            if not relations:
                continue
            best_relation_name = ""
            best_status = ""
            best_priority = (-1, -1)
            for relation_name, same_position in relations:
                status = _classify_relation_status(query_sign, occ.effect_sign, same_position)
                priority = (_status_priority(status), _relation_priority(relation_name))
                if priority > best_priority:
                    best_priority = priority
                    best_relation_name = relation_name
                    best_status = status
            detail_label = f"{occ.comparison_label}({best_status})"
            dataset_key = (occ.folder_label, occ.comparison_label)
            replacement = (best_status, best_relation_name, detail_label)
            existing = best_by_folder_dataset.get(dataset_key)
            if existing is None:
                best_by_folder_dataset[dataset_key] = replacement
            else:
                existing_priority = (_status_priority(existing[0]), _relation_priority(existing[1]))
                new_priority = (_status_priority(best_status), _relation_priority(best_relation_name))
                if new_priority > existing_priority:
                    best_by_folder_dataset[dataset_key] = replacement

        row_out["cross_comparison_feature_key"] = canonical_key
        for folder_label, (matches_column, count_column) in folder_columns.items():
            folder_matches = [
                value[2]
                for (matched_folder, _comparison_label), value in sorted(best_by_folder_dataset.items())
                if matched_folder == folder_label
            ]
            row_out[matches_column] = "|".join(folder_matches)
            row_out[count_column] = str(len(folder_matches))
        annotated_rows.append(row_out)

        if idx % progress_every == 0 or idx == len(input_df):
            print(f"[event_cross_comparison] Processed {idx}/{len(input_df)} rows")

    result_df = pd.DataFrame(annotated_rows)
    if output_tsv is None:
        base, ext = os.path.splitext(input_tsv)
        output_tsv = f"{base}-cross_comparison{ext or '.tsv'}"
    print(f"[event_cross_comparison] Writing output: {output_tsv}")
    result_df.to_csv(output_tsv, sep="\t", index=False)
    print(f"[event_cross_comparison] Completed: wrote {len(result_df)} rows")
    return result_df


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Annotate a splicing-event TSV with cross-dataset/condition event matches "
            "from one or more folders of AltAnalyze2/AltAnalyze3 result files."
        )
    )
    parser.add_argument("input_tsv", help="Input event TSV (AltAnalyze2, AltAnalyze3, or isoform_predict output).")
    parser.add_argument(
        "comparison_roots",
        nargs="+",
        help="One or more folder roots to scan recursively for event result files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output TSV path. Defaults to <input_tsv>-cross_comparison.tsv",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    annotate_event_cross_comparison(
        input_tsv=args.input_tsv,
        comparison_roots=args.comparison_roots,
        output_tsv=args.output,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
