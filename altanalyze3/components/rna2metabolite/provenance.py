"""Audit AML metabolite identities against the original study workbook.

Keep model IDs stable. A source synonym or molecular formula is evidence to
investigate, not permission to rename an unidentified mass-spectrometry feature.
This command reads local source files and writes only the requested audit folder.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile

SOURCE_DOI = "10.1038/s43018-026-01175-6"
SOURCE_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs43018-026-01175-6/MediaObjects/"
    "43018_2026_1175_MOESM3_ESM.xlsx"
)
NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def workbook_rows(path, sheet):
    """Yield (Excel row number, {column letter: literal value}) without pandas."""
    with zipfile.ZipFile(path) as archive:
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            strings = ["".join(e.itertext()) for e in ET.fromstring(
                archive.read("xl/sharedStrings.xml")).findall("s:si", NS)]
        relations = {e.attrib["Id"]: e.attrib["Target"] for e in ET.fromstring(
            archive.read("xl/_rels/workbook.xml.rels"))}
        sheets = ET.fromstring(archive.read("xl/workbook.xml"))
        element = next(e for e in sheets.findall("s:sheets/s:sheet", NS)
                       if e.attrib["name"] == sheet)
        target = relations[element.attrib[
            "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]]
        target = target.lstrip("/") if target.startswith("/") else "xl/" + target
        for row in ET.fromstring(archive.read(target)).findall("s:sheetData/s:row", NS):
            values = {}
            for cell in row:
                node = cell.find("s:v", NS)
                value = node.text or "" if node is not None else ""
                if cell.attrib.get("t") == "s":
                    value = strings[int(value)]
                elif cell.attrib.get("t") == "inlineStr":
                    node = cell.find("s:is", NS)
                    value = "".join(node.itertext()) if node is not None else ""
                values[re.sub(r"\d+", "", cell.attrib["r"])] = value
            yield int(row.attrib["r"]), values


def unknown(name):
    return name.lower().startswith("unknown")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit(workbook, selection_file, pdc_table, output):
    """Join by the exact original name AND the assay selected during training."""
    source = {}
    for sheet, assay in [("Table21", "HILIC"), ("Table22", "RP")]:
        rows = iter(workbook_rows(workbook, sheet))
        _, headers = next(rows)
        for number, cells in rows:
            name = cells.get("A", "").strip()
            if not name:
                continue
            # Exclude patient measurements; retain literal annotation headers,
            # including the apparent HILIC H:K header inconsistency.
            metadata = {header: cells.get(col, "") for col, header in headers.items()
                        if not re.search(r"C3[LN][._-]\d{5}", header)}
            key = (name, assay)
            if key in source:
                raise ValueError(f"Ambiguous source identity: {key}")
            source[key] = (sheet, number, metadata)
    with open(pdc_table, newline="") as handle:
        pdc_rows = list(csv.reader(handle))
    pdc_ids = {row[0] for row in pdc_rows[1:] if row}
    ica_ids = {cells.get("A", "") for n, cells in workbook_rows(workbook, "Table24") if n > 1}
    with open(selection_file, newline="") as handle:
        selection = list(csv.DictReader(handle, delimiter="\t"))
    records = []
    for selected in selection:
        name = selected[next(iter(selected))]
        assay = selected["protocol"]
        sheet, number, raw = source[(name, assay)]
        pdc_id = assay + ":" + name
        if pdc_id not in pdc_ids:
            raise ValueError(f"Source target missing from deposited PDC table: {pdc_id}")
        candidate = raw.get("Synonyms.or.Isomers", "")
        notes = []
        if assay == "HILIC":
            notes.append("Table21 H:K headers appear inconsistent with their values; retained verbatim, not reinterpreted")
        if unknown(name) and candidate:
            notes.append("Candidate/isomer text only; not a confirmed identity")
        if name == "Unknown0 1621":
            notes.append("Sarcosine candidate conflicts with source C3 H9 N O2 and m/z 92.07052; requires author clarification")
        records.append({
            "model_feature_id": name, "source_name": raw["Metabolite"],
            "identification_status": "unidentified_in_source" if unknown(name) else "named_in_source_confidence_not_standardized",
            "assay": assay, "source_sheet": sheet, "source_excel_row": number,
            "source_cell": f"{sheet}!A{number}", "mz": raw.get("mz", ""),
            "retention_time_min": raw.get("Rt.min", ""), "source_tag": raw.get("Tags", ""),
            "source_synonyms_or_isomers": candidate,
            "source_kegg_id": raw.get("KEGG.ID..If.applicable.", ""),
            "source_formula": raw.get("Formula", ""),
            "pdc_feature_id": pdc_id, "in_paper_ica_table24": pdc_id in ica_ids,
            "training_detection_fraction": selected["detection"],
            "training_median_signal": selected["median_signal"],
            "source_doi": SOURCE_DOI, "source_url": SOURCE_URL,
            "notes": "; ".join(notes),
            "raw_source_annotations_json": json.dumps(raw, ensure_ascii=False, sort_keys=True),
        })
    if len({r["model_feature_id"] for r in records}) != len(records):
        raise ValueError("Duplicate selected target names")
    unidentified = [r for r in records if r["identification_status"] == "unidentified_in_source"]
    report = {
        "source_doi": SOURCE_DOI, "source_url": SOURCE_URL,
        "inputs": {key: {"path": str(Path(path).resolve()), "sha256": sha256(path)}
                   for key, path in [("workbook", workbook), ("training_selection", selection_file), ("pdc_table", pdc_table)]},
        "source_measurements": len(source), "pdc_measurements": len(pdc_ids),
        "pdc_samples": len(pdc_rows[0]) - 1,
        "source_pdc_feature_ids_identical": {assay + ":" + name for name, assay in source} == pdc_ids,
        "model_targets": len(records), "named_targets": len(records) - len(unidentified),
        "unidentified_targets": len(unidentified),
        "unidentified_with_candidate_text": sum(bool(r["source_synonyms_or_isomers"]) for r in unidentified),
        "unidentified_with_formula": sum(bool(r["source_formula"]) for r in unidentified),
        "unidentified_with_kegg": sum(bool(r["source_kegg_id"]) for r in unidentified),
        "named_with_kegg": sum(bool(r["source_kegg_id"]) for r in records if r not in unidentified),
        "paper_ica_measurements": len(ica_ids),
        "paper_ica_unknown_measurements": sum(unknown(x.split(":", 1)[-1]) for x in ica_ids),
        "selected_measurements_in_paper_ica": sum(r["in_paper_ica_table24"] for r in records),
        "confirmed_unknown_renamings": 0,
    }
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for filename, subset in [("aml_metabolite_provenance.tsv", records),
                             ("unidentified_features_for_authors.tsv", unidentified)]:
        with (output / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(subset)
    (output / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", required=True)
    parser.add_argument("--selection", required=True)
    parser.add_argument("--pdc-table", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.workbook, args.selection, args.pdc_table, args.output), indent=2))


if __name__ == "__main__":
    main()
