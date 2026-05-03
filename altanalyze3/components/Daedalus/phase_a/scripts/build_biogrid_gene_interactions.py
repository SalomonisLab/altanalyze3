#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path
from zipfile import ZipFile


RAW = Path(__file__).resolve().parents[1] / "data" / "raw" / "BIOGRID-ALL-LATEST.tab3.zip"
OUT = Path(__file__).resolve().parents[1] / "data" / "interim" / "biogrid_gene_interactions.tsv"

TAXON_TO_SPECIES = {"9606": "human", "10090": "mouse"}


def _open_tab3_member(zf: ZipFile) -> str:
    for name in zf.namelist():
        lowered = name.lower()
        if lowered.endswith(".tab3.txt") or lowered.endswith(".tab3.tsv"):
            return name
    raise FileNotFoundError("No tab3 table found in BioGRID archive.")


def main() -> int:
    partner_sets: dict[tuple[str, str], set[str]] = defaultdict(set)
    physical_partner_sets: dict[tuple[str, str], set[str]] = defaultdict(set)
    interaction_counts: Counter = Counter()
    physical_counts: Counter = Counter()

    with ZipFile(RAW) as zf:
        member = _open_tab3_member(zf)
        with zf.open(member) as handle:
            reader = csv.DictReader((line.decode("utf-8", errors="replace") for line in handle), delimiter="\t")
            for row in reader:
                tax_a = str(row.get("Organism ID Interactor A", row.get("Organism Interactor A", "")))
                tax_b = str(row.get("Organism ID Interactor B", row.get("Organism Interactor B", "")))
                if tax_a != tax_b:
                    continue
                species = TAXON_TO_SPECIES.get(tax_a)
                if not species:
                    continue
                gene_a = (row.get("Official Symbol Interactor A") or "").strip()
                gene_b = (row.get("Official Symbol Interactor B") or "").strip()
                if not gene_a or not gene_b or gene_a == "-" or gene_b == "-":
                    continue
                key_a = (species, gene_a)
                key_b = (species, gene_b)
                partner_sets[key_a].add(gene_b)
                partner_sets[key_b].add(gene_a)
                interaction_counts[key_a] += 1
                interaction_counts[key_b] += 1
                if (row.get("Experimental System Type") or "").strip().lower() == "physical":
                    physical_partner_sets[key_a].add(gene_b)
                    physical_partner_sets[key_b].add(gene_a)
                    physical_counts[key_a] += 1
                    physical_counts[key_b] += 1

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "species",
            "gene_name",
            "biogrid_partner_count",
            "biogrid_interaction_count",
            "biogrid_physical_partner_count",
            "biogrid_physical_interaction_count",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        keys = sorted(set(partner_sets) | set(physical_partner_sets) | set(interaction_counts) | set(physical_counts))
        for species, gene_name in keys:
            key = (species, gene_name)
            writer.writerow(
                {
                    "species": species,
                    "gene_name": gene_name,
                    "biogrid_partner_count": len(partner_sets.get(key, set())),
                    "biogrid_interaction_count": interaction_counts.get(key, 0),
                    "biogrid_physical_partner_count": len(physical_partner_sets.get(key, set())),
                    "biogrid_physical_interaction_count": physical_counts.get(key, 0),
                }
            )
    print(f"[ok] Wrote {OUT} rows={len(keys)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
