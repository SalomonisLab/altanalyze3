"""Extract original-scale LungMAP targets without altering negative log values.

Uses cached XLSX values and formulas; requires no Excel reader dependency.
The workbook's own fold formulas establish its intended log2 interpretation.
This does not recover raw intensities or validate RNA-to-lipid predictions.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def read_sheet(path):
    with zipfile.ZipFile(path) as archive:
        strings = ["".join(node.itertext()) for node in
                   ET.fromstring(archive.read("xl/sharedStrings.xml"))]
        root = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    rows, formulas = {}, {}
    for row in root.findall(".//m:row", NS):
        values = {}
        for cell in row:
            address = cell.get("r")
            column = re.sub(r"\d", "", address)
            formula = cell.find("m:f", NS)
            if formula is not None and formula.text:
                formulas[address] = formula.text
            value = cell.find("m:v", NS)
            if value is None:
                continue
            if cell.get("t") == "s":
                text = strings[int(value.text)]
                values[column] = None if text in {"NA", "NaN", "nan", ""} else text
            elif cell.get("t") == "e":
                values[column] = None
            else:
                values[column] = float(value.text)
        rows[int(row.get("r"))] = values
    return rows, formulas


def extract(path, output_dir):
    rows, formulas = read_sheet(path)
    # Select real sample columns by identifiers, excluding averages and folds.
    columns = [c for c, v in rows[3].items()
               if isinstance(v, str) and re.fullmatch(r"Sample_\d+_\w+", v)]
    if len(columns) != 90:
        raise ValueError("Expected two blocks of 45 labeled samples.")
    left, right = columns[:45], columns[45:]
    if [rows[3][c] for c in left] != [rows[3][c] for c in right]:
        raise ValueError("Sample blocks disagree.")
    if formulas.get("BA4") != "(B4-$AY4)/$AZ4":
        raise ValueError("Unexpected standardization formula.")
    if formulas.get("AO4") != "IF(AN4<0,-1*2^(ABS(AN4)),2^AN4)":
        raise ValueError("Expected workbook log2-difference-to-fold formula.")
    lipid_rows = [(i, r) for i, r in rows.items() if i >= 4 and r.get("A")]
    lipids = [r["A"] for _, r in lipid_rows]
    if len(set(lipids)) != len(lipids):
        raise ValueError("Duplicate feature IDs require explicit resolution.")
    errors, observed, missing, scales = [], [], [], []
    for i, row in lipid_rows:
        scales.append([row["A"], row.get("AY"), row.get("AZ")])
        for a, b in zip(left, right):
            x, z = row.get(a), row.get(b)
            if x is None:
                missing.append([row["A"], rows[3][a], a + str(i)])
                continue
            if not math.isfinite(x):
                raise ValueError("Nonfinite measurement.")
            observed.append(x)
            if z is not None and row.get("AY") is not None and row.get("AZ") is not None:
                errors.append(abs(x - (z * row["AZ"] + row["AY"])))
    if not errors or max(errors) > 1e-10:
        raise ValueError("Standardized block fails inversion check.")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    def write(name, header, data):
        with (output_dir / name).open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(data)

    write("targets_native_log2.csv", ["sample"] + lipids,
          [[rows[3][c]] + [r.get(c) for _, r in lipid_rows] for c in left])
    complete_rows = [(i, r) for i, r in lipid_rows if all(r.get(c) is not None for c in left)]
    write("targets_native_log2_complete.csv", ["sample"] + [r["A"] for _, r in complete_rows],
          [[rows[3][c]] + [r[c] for _, r in complete_rows] for c in left])
    write("stored_standardization.csv", ["lipid", "stored_mean", "stored_sd"], scales)
    write("missing_values.csv", ["lipid", "sample", "source_cell"], missing)
    write("sample_metadata.csv", ["sample", "source_donor", "cell_type", "age_days", "source_column"],
          [[rows[3][c], "_".join(rows[3][c].split("_")[:2]), rows[2][c], rows[1][c], c] for c in left])
    report = {
        "source": str(Path(path).resolve()),
        "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        "samples": len(left), "donors": len({rows[3][c].split("_")[1] for c in left}),
        "lipid_features": len(lipids), "observed_values": len(observed),
        "complete_lipid_features": len(complete_rows),
        "missing_values": len(missing), "negative_values_preserved": sum(x < 0 for x in observed),
        "native_range": [min(observed), max(observed)],
        "standardization_comparisons": len(errors), "max_inverse_error": max(errors),
        "scale_evidence": {k: formulas.get(k) for k in ["AN4", "AO4", "AY4", "AZ4", "BA4"]},
        "interpretation": "Relative log2 measurements, supported by workbook fold formulas; absolute intensities and detailed original normalization are unavailable.",
        "important": [
            "AY/AZ summarize B:AX, which includes calculated columns. Do not recompute these from only the 45 samples to invert the right block.",
            "No clipping, offset, logging, normalization or missing-value imputation was applied.",
            "Keep positive/negative ion-mode suffixes and ambiguous feature annotations; do not merge by stripping them.",
            "Sample_N to DNNN RNA matching is a candidate identifier mapping that requires provenance confirmation.",
            "This workbook has 9 donors; the shipped model used 10. This extraction does not reproduce its training targets.",
        ],
    }
    (output_dir / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    print(json.dumps(extract(args.workbook, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
