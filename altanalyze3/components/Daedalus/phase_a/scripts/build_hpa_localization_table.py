#!/usr/bin/env python3
from __future__ import annotations

import csv
import io
import zipfile
from pathlib import Path


INPUT_FILE = Path(__file__).resolve().parents[1] / "data" / "raw" / "subcellular_location.tsv.zip"
OUTPUT_FILE = Path(__file__).resolve().parents[1] / "data" / "interim" / "hpa_localization.tsv"


def _split_locations(value: str) -> list[str]:
    if not value:
        return []
    return [part.strip() for part in value.split(";") if part.strip()]


def main() -> int:
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(INPUT_FILE) as archive:
        member = archive.namelist()[0]
        with archive.open(member) as handle:
            text_stream = io.TextIOWrapper(handle, encoding="utf-8", newline="")
            reader = csv.DictReader(text_stream, delimiter="\t")
            fieldnames = [
                "ensembl_gene_id",
                "gene_name",
                "reliability",
                "main_locations",
                "additional_locations",
                "extracellular_locations",
                "all_locations",
                "is_membrane_or_surface_annotated",
                "has_extracellular_annotation"
            ]
            with OUTPUT_FILE.open("w", encoding="utf-8", newline="") as out_handle:
                writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                for row in reader:
                    main_locations = _split_locations(row.get("Main location", ""))
                    additional_locations = _split_locations(row.get("Additional location", ""))
                    extracellular_locations = _split_locations(row.get("Extracellular location", ""))
                    all_locations = sorted(set(main_locations + additional_locations + extracellular_locations))
                    is_membrane_or_surface = any(
                        location in {
                            "Plasma membrane",
                            "Cell Junctions",
                            "Vesicles",
                            "Golgi apparatus",
                            "Endoplasmic reticulum"
                        }
                        for location in all_locations
                    )
                    writer.writerow(
                        {
                            "ensembl_gene_id": row["Gene"],
                            "gene_name": row["Gene name"],
                            "reliability": row["Reliability"],
                            "main_locations": ";".join(main_locations),
                            "additional_locations": ";".join(additional_locations),
                            "extracellular_locations": ";".join(extracellular_locations),
                            "all_locations": ";".join(all_locations),
                            "is_membrane_or_surface_annotated": int(is_membrane_or_surface),
                            "has_extracellular_annotation": int(bool(extracellular_locations))
                        }
                    )
    print(f"[ok] Wrote {OUTPUT_FILE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
