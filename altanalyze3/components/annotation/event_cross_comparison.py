"""Cross-annotate splicing events against other event result files."""

from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd


@dataclass
class EventFileInfo:
    path: str
    format_name: str
    feature_column: str
    comparison_label: str
    root_label: str


@dataclass
class EventOccurrence:
    comparison_label: str
    root_label: str
    path: str
    format_name: str
    feature: str


def _clean_string(value: object) -> str:
    if value is None:
        return ""
    if pd.isna(value):
        return ""
    return str(value).strip()


def _detect_event_file_format(path: str) -> Optional[Tuple[str, str]]:
    try:
        with open(path, "r", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
    except Exception:
        return None
    header_set = {_clean_string(col) for col in header}
    if "Feature" in header_set:
        return "altanalyze3", "Feature"
    if "UID" in header_set:
        return "altanalyze2", "UID"
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


def _iter_event_files(folder_roots: Sequence[str]) -> List[EventFileInfo]:
    infos: List[EventFileInfo] = []
    for root in folder_roots:
        root_abs = os.path.abspath(root)
        root_label = os.path.basename(root_abs.rstrip(os.sep)) or root_abs
        for dirpath, _dirnames, filenames in os.walk(root_abs):
            for filename in sorted(filenames):
                if not filename.endswith((".tsv", ".txt")):
                    continue
                path = os.path.join(dirpath, filename)
                detected = _detect_event_file_format(path)
                if not detected:
                    continue
                format_name, feature_column = detected
                infos.append(
                    EventFileInfo(
                        path=os.path.abspath(path),
                        format_name=format_name,
                        feature_column=feature_column,
                        comparison_label=_infer_comparison_label(path, format_name),
                        root_label=root_label,
                    )
                )
    return infos


def _feature_from_row(row: pd.Series, feature_column: str) -> str:
    return _clean_string(row.get(feature_column))


def _split_feature(feature: str) -> Tuple[str, str]:
    parts = _clean_string(feature).split("|")
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


def _index_comparison_events(
    comparison_files: Sequence[EventFileInfo],
    input_path: str,
) -> Dict[str, List[EventOccurrence]]:
    event_index: Dict[str, List[EventOccurrence]] = {}
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
                canonical_key = _canonical_feature_key(feature)
            except ValueError:
                continue
            occurrences = event_index.setdefault(canonical_key, [])
            occurrences.append(
                EventOccurrence(
                    comparison_label=info.comparison_label,
                    root_label=info.root_label,
                    path=info.path,
                    format_name=info.format_name,
                    feature=feature,
                )
            )
    return event_index


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
    input_format_name, input_feature_column = input_detected
    print(f"[event_cross_comparison] Input format: {input_format_name} ({input_feature_column})")
    print(f"[event_cross_comparison] Loading input: {input_tsv}")
    input_df = pd.read_csv(input_tsv, sep="\t", dtype=str).fillna("")
    print(f"[event_cross_comparison] Loaded {len(input_df)} input rows")

    comparison_files = _iter_event_files(comparison_roots)
    print(f"[event_cross_comparison] Found {len(comparison_files)} comparison files")
    event_index = _index_comparison_events(comparison_files, input_tsv)
    print(f"[event_cross_comparison] Indexed {len(event_index)} unique event keys")

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
            reversed_feature = f"{right}|{left}"
        except ValueError:
            row_out["cross_comparison_present_elsewhere"] = "FALSE"
            row_out["cross_comparison_match_count"] = "0"
            row_out["cross_comparison_matches"] = ""
            row_out["cross_comparison_root_labels"] = ""
            row_out["cross_comparison_same_orientation_matches"] = ""
            row_out["cross_comparison_reversed_orientation_matches"] = ""
            row_out["cross_comparison_feature_key"] = ""
            annotated_rows.append(row_out)
            continue

        occurrences = event_index.get(canonical_key, [])
        all_matches: List[str] = []
        same_orientation: List[str] = []
        reversed_orientation: List[str] = []
        root_labels: List[str] = []
        for occ in occurrences:
            detail = f"{occ.root_label}:{occ.comparison_label}"
            if detail not in all_matches:
                all_matches.append(detail)
            if occ.root_label not in root_labels:
                root_labels.append(occ.root_label)
            if occ.feature == feature:
                same_orientation.append(detail)
            elif occ.feature == reversed_feature:
                reversed_orientation.append(detail)

        row_out["cross_comparison_present_elsewhere"] = "TRUE" if all_matches else "FALSE"
        row_out["cross_comparison_match_count"] = str(len(all_matches))
        row_out["cross_comparison_matches"] = "|".join(all_matches)
        row_out["cross_comparison_root_labels"] = "|".join(root_labels)
        row_out["cross_comparison_same_orientation_matches"] = "|".join(sorted(set(same_orientation)))
        row_out["cross_comparison_reversed_orientation_matches"] = "|".join(sorted(set(reversed_orientation)))
        row_out["cross_comparison_feature_key"] = canonical_key
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
