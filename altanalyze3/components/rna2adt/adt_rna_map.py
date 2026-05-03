"""Curated mapping from CITE-seq ADT names to canonical RNA gene partners.

Each ADT name (the literal column in the bone marrow atlas, prefixed ``Hu.``)
maps to one or more HGNC official symbols whose protein product the antibody
binds. For multi-subunit proteins (e.g. CD3 complex, CD8 dimer, CD45 isoforms)
all relevant chains are listed; the targeted regressor will weight them
through fitting.

The map is the source of truth — the writer in ``write_map_to_tsv`` produces
``configs/adt_rna_map.tsv`` for inspection / hand-edit, and the loader
``load_curated_adt_rna_map`` reads it back. Genes not present in the input
RNA matrix are reported by ``audit_against_genes`` so callers can either
fall back or surface them to the user for manual correction.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


CURATED_MAP: Dict[str, Tuple[str, ...]] = {
    # --- Lineage / core surface ---
    "Hu.CD2":          ("CD2",),
    "Hu.CD3":          ("CD3D", "CD3E", "CD3G", "CD247"),
    "Hu.CD4":          ("CD4",),
    "Hu.CD4_RPA.T4":   ("CD4",),
    "Hu.CD5":          ("CD5",),
    "Hu.CD7":          ("CD7",),
    "Hu.CD8":          ("CD8A", "CD8B"),
    "Hu.CD9":          ("CD9",),
    "Hu.CD10":         ("MME",),
    "Hu.CD11a":        ("ITGAL",),
    "Hu.CD11b":        ("ITGAM",),
    "Hu.CD11c":        ("ITGAX",),
    "Hu.CD13":         ("ANPEP",),
    "Hu.CD14_M5E2":    ("CD14",),
    "Hu.CD15_W6D3":    ("FUT4",),
    "Hu.CD16":         ("FCGR3A", "FCGR3B"),
    "Hu.CD18":         ("ITGB2",),
    "Hu.CD19":         ("CD19",),
    "Hu.CD22":         ("CD22",),  # not in panel but kept for safety
    "Hu.CD24":         ("CD24",),
    "Hu.CD25":         ("IL2RA",),
    "Hu.CD26":         ("DPP4",),
    "Hu.CD27":         ("CD27",),
    "Hu.CD28":         ("CD28",),
    "Hu.CD29":         ("ITGB1",),
    "Hu.CD32":         ("FCGR2A", "FCGR2B", "FCGR2C"),
    "Hu.CD33":         ("CD33",),
    "Hu.CD34":         ("CD34",),
    "Hu.CD35":         ("CR1",),
    "Hu.CD36":         ("CD36",),
    "Hu.CD37":         ("CD37",),
    "Hu.CD38_HIT2":    ("CD38",),
    "Hu.CD41":         ("ITGA2B",),
    "Hu.CD42b":        ("GP1BA",),
    "Hu.CD43":         ("SPN",),
    "Hu.CD45_2D1":     ("PTPRC",),
    "Hu.CD45RA":       ("PTPRC",),  # isoforms — same gene, alternative splicing
    "Hu.CD45RB":       ("PTPRC",),
    "Hu.CD45RO":       ("PTPRC",),
    "Hu.CD47":         ("CD47",),
    "Hu.CD49b":        ("ITGA2",),
    "Hu.CD52":         ("CD52",),
    "Hu.CD54":         ("ICAM1",),
    "Hu.CD55":         ("CD55",),
    "Hu.CD56":         ("NCAM1",),
    "Hu.CD57":         ("B3GAT1",),
    "Hu.CD58":         ("CD58",),
    "Hu.CD59":         ("CD59",),
    "Hu.CD61":         ("ITGB3",),
    "Hu.CD62L":        ("SELL",),
    "Hu.CD62P":        ("SELP",),
    "Hu.CD63":         ("CD63",),
    "Hu.CD64":         ("FCGR1A", "FCGR1B"),
    "Hu.CD69":         ("CD69",),
    "Hu.CD71":         ("TFRC",),
    "Hu.CD72":         ("CD72",),
    "Hu.CD73":         ("NT5E",),
    "Hu.CD81":         ("CD81",),
    "Hu.CD82":         ("CD82",),
    "Hu.CD83":         ("CD83",),
    "Hu.CD84":         ("CD84",),
    "Hu.CD85g":        ("LILRB4",),
    "Hu.CD90":         ("THY1",),
    "Hu.CD93":         ("CD93",),
    "Hu.CD98":         ("SLC3A2",),

    # --- Cytokine / growth-factor receptors ---
    "Hu.CD101":        ("CD101",),
    "Hu.CD102":        ("ICAM2",),
    "Hu.CD103":        ("ITGAE",),
    "Hu.CD105_43A3":   ("ENG",),
    "Hu.CD106":        ("VCAM1",),
    "Hu.CD109":        ("CD109",),
    "Hu.CD110":        ("MPL",),
    "Hu.CD112":        ("NECTIN2",),
    "Hu.CD115":        ("CSF1R",),
    "Hu.CD116":        ("CSF2RA",),
    "Hu.CD117":        ("KIT",),
    "Hu.CD119":        ("IFNGR1",),
    "Hu.CD123":        ("IL3RA",),
    "Hu.CD127":        ("IL7R",),
    "Hu.CD133_S16016B":("PROM1",),
    "Hu.CD135":        ("FLT3",),
    "Hu.CD138_DL.101": ("SDC1",),
    "Hu.CD140b":       ("PDGFRB",),
    "Hu.CD141":        ("THBD",),
    "Hu.CD150":        ("SLAMF1",),
    "Hu.CD151":        ("CD151",),
    "Hu.CD154":        ("CD40LG",),
    "Hu.CD155":        ("PVR",),
    "Hu.CD158e1":      ("KIR3DL1",),
    "Hu.CD158f":       ("KIR2DL5A", "KIR2DL5B"),
    "Hu.CD162":        ("SELPLG",),
    "Hu.CD163":        ("CD163",),
    "Hu.CD164":        ("CD164",),
    "Hu.CD172a":       ("SIRPA",),
    "Hu.CD177":        ("CD177",),
    "Hu.CD183":        ("CXCR3",),
    "Hu.CD185":        ("CXCR5",),
    "Hu.CD186":        ("CXCR6",),
    "Hu.CD192":        ("CCR2",),
    "Hu.CD1a":         ("CD1A",),
    "Hu.CD1d":         ("CD1D",),
    "Hu.CD200":        ("CD200",),
    "Hu.CD201":        ("PROCR",),
    "Hu.CD202b":       ("TEK",),
    "Hu.CD205":        ("LY75",),
    "Hu.CD226_TX25":   ("CD226",),
    "Hu.CD235a":       ("GYPA",),
    "Hu.CD271":        ("NGFR",),
    "Hu.CD274":        ("CD274",),
    "Hu.CD279":        ("PDCD1",),
    "Hu.CD304":        ("NRP1",),
    "Hu.CD305_LAIR1":  ("LAIR1",),
    "Hu.CD309":        ("KDR",),
    "Hu.CD325":        ("CDH2",),
    "Hu.CD326":        ("EPCAM",),
    "Hu.CD335":        ("NCR1",),
    "Hu.CD354":        ("TREM1",),
    "Hu.CD366":        ("HAVCR2",),

    # --- Non-CD-numbered ---
    "Hu.C5L2":         ("C5AR2",),
    "Hu.CLEC1B":       ("CLEC1B",),
    "Hu.Cadherin.11":  ("CDH11",),
    "Hu.FR.b":         ("FOLR2",),
    "Hu.FceRIa":       ("FCER1A",),
    "Hu.GARP":         ("LRRC32",),
    "Hu.GPR56":        ("ADGRG1",),
    "Hu.Galectin.9":   ("LGALS9",),
    "Hu.HLA.ABC":      ("HLA-A", "HLA-B", "HLA-C"),
    "Hu.HLA.DR.DP.DQ": ("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1"),
    "Hu.KLRG1":        ("KLRG1",),
    "Hu.TIM.4":        ("TIMD4",),
    "Hu.TSPAN33":      ("TSPAN33",),

    # --- Isotype controls (no biologically meaningful RNA partner) ---
    "Hu.IgG.Fc":       (),  # isotype control
}


DEFAULT_MAP_TSV = Path(__file__).parent / "configs" / "adt_rna_map.tsv"


def write_map_to_tsv(path: Path = DEFAULT_MAP_TSV) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for adt, genes in sorted(CURATED_MAP.items()):
        rows.append({
            "adt_name": adt,
            "rna_genes": ",".join(genes),
            "notes": "isotype" if not genes else "",
        })
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["adt_name", "rna_genes", "notes"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def load_curated_adt_rna_map(path: Path | None = None) -> Dict[str, List[str]]:
    """Load the ADT->RNA mapping, preferring the on-disk TSV (so users can
    override) and falling back to the in-code CURATED_MAP."""
    if path is None:
        path = DEFAULT_MAP_TSV
    if not path.exists():
        return {adt: list(genes) for adt, genes in CURATED_MAP.items()}
    out: Dict[str, List[str]] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            adt = (row.get("adt_name") or "").strip()
            if not adt:
                continue
            raw = (row.get("rna_genes") or "").strip()
            genes = [g.strip() for g in raw.split(",") if g.strip()]
            out[adt] = genes
    return out


def audit_against_genes(
    adt_names: Iterable[str],
    rna_genes: Sequence[str],
    *,
    map_obj: Dict[str, List[str]] | None = None,
) -> Dict[str, Dict[str, object]]:
    """Resolve every ADT against an RNA panel, splitting matched/missing.

    Returns a dict keyed by ADT name with fields:
        matched   : list of mapped genes that ARE in rna_genes
        missing   : list of mapped genes that are NOT in rna_genes
        unmapped  : True if the ADT is not in the curated map at all
        isotype   : True if the curated entry is empty (intentional null)
    """
    if map_obj is None:
        map_obj = load_curated_adt_rna_map()
    rna_set = set(rna_genes)
    audit: Dict[str, Dict[str, object]] = {}
    for adt in adt_names:
        entry = map_obj.get(adt)
        if entry is None:
            audit[adt] = {"matched": [], "missing": [], "unmapped": True, "isotype": False}
            continue
        if not entry:
            audit[adt] = {"matched": [], "missing": [], "unmapped": False, "isotype": True}
            continue
        matched = [g for g in entry if g in rna_set]
        missing = [g for g in entry if g not in rna_set]
        audit[adt] = {"matched": matched, "missing": missing, "unmapped": False, "isotype": False}
    return audit


def print_audit_report(audit: Dict[str, Dict[str, object]]) -> None:
    unmapped = sorted(adt for adt, info in audit.items() if info["unmapped"])
    isotypes = sorted(adt for adt, info in audit.items() if info["isotype"])
    no_match = sorted(adt for adt, info in audit.items()
                      if not info["unmapped"] and not info["isotype"] and not info["matched"])
    partial = sorted(
        (adt, info["matched"], info["missing"]) for adt, info in audit.items()
        if info["matched"] and info["missing"]
    )
    fully_matched = sorted(adt for adt, info in audit.items()
                           if info["matched"] and not info["missing"]
                           and not info["unmapped"] and not info["isotype"])

    def _hdr(label: str, n: int) -> str:
        return f"\n=== {label} ({n}) ==="

    print(_hdr("Fully matched", len(fully_matched)))
    for adt in fully_matched:
        print(f"  {adt}\t{','.join(audit[adt]['matched'])}")
    if partial:
        print(_hdr("Partially matched (some genes not in atlas)", len(partial)))
        for adt, matched, missing in partial:
            print(f"  {adt}\tmatched={','.join(matched)}\tmissing={','.join(missing)}")
    if no_match:
        print(_hdr("Mapped but NO genes found in atlas — needs manual symbol", len(no_match)))
        for adt in no_match:
            print(f"  {adt}\t(curated genes: {','.join(load_curated_adt_rna_map().get(adt, []))})")
    if unmapped:
        print(_hdr("Not in curated map — please supply RNA partners", len(unmapped)))
        for adt in unmapped:
            print(f"  {adt}")
    if isotypes:
        print(_hdr("Isotype controls (intentionally no RNA partner)", len(isotypes)))
        for adt in isotypes:
            print(f"  {adt}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-tsv", action="store_true", help="Write configs/adt_rna_map.tsv")
    parser.add_argument("--audit-h5ad", type=Path, help="Audit the curated map against an .h5ad's var_names")
    args = parser.parse_args()
    if args.write_tsv:
        path = write_map_to_tsv()
        print(f"Wrote {path}")
    if args.audit_h5ad:
        import anndata as ad
        a = ad.read_h5ad(args.audit_h5ad, backed="r")
        var_names = [str(v) for v in a.var_names]
        adt_names = [v for v in var_names if v.startswith("Hu.")]
        rna_genes = [v for v in var_names if not v.startswith("Hu.")]
        report = audit_against_genes(adt_names, rna_genes)
        print_audit_report(report)
