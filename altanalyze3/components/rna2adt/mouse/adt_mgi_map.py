"""Curated mapping from murine CITE-seq ADT names to canonical MGI gene partners.

Mirrors ``components.rna2adt.adt_rna_map`` but for the mouse panel. The
input column names follow the ``ADT-<marker>[_<clone>]`` convention seen in
the training atlas
(``adata_combined_60k_rna_adt.h5ad``); for each ADT we list its mouse
encoding gene(s) using MGI official symbols (mouse convention =
sentence-case, e.g. ``Cd19``, ``Itgam``, ``Ly6a``).

Multi-subunit complexes list every chain (CD3 = ``Cd3d, Cd3e, Cd3g, Cd247``;
CD8 = ``Cd8a, Cd8b1``; CD16/32 = ``Fcgr3, Fcgr4``; FcγR2 = ``Fcgr2b``).
Cross-reactive antibodies that bind both human and mouse epitopes are
treated as mouse-only here since the training atlas is murine.

Isotype controls (Mouse_IgG, Rat_IgG, Armenian_Hamster_IgG) have no
biological RNA partner and are intentionally excluded from the model's
output set during whitelist generation.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


_ADT_PREFIX = "ADT-"


def strip_prefix(name: str) -> str:
    """Drop the ``ADT-`` prefix and the species-cross-reactivity wrapper.

    e.g. 'ADT-anti_mouse_human_CD11b' -> 'CD11b',
         'ADT-CD8a' -> 'CD8a',
         'ADT-CD183_CXCR3' -> 'CD183_CXCR3'.
    Notably we DO NOT strip the trailing ``_<clone>`` because some clones
    carry biologically meaningful information (e.g. CD45_1 vs CD45_2 are
    different alleles, not different antibodies of the same protein).
    """
    n = name
    if n.startswith(_ADT_PREFIX):
        n = n[len(_ADT_PREFIX):]
    n = re.sub(r'^anti_(?:mouse|human|rat)(?:_(?:mouse|human|rat))?_', '', n)
    return n


# --- Curated map: ADT raw name (with ADT- prefix, as in atlas var) -> mouse gene(s) ---
CURATED_MAP: Dict[str, Tuple[str, ...]] = {
    # T cell core
    "ADT-CD2":                            ("Cd2",),
    "ADT-CD3":                            ("Cd3d", "Cd3e", "Cd3g", "Cd247"),
    "ADT-CD4":                            ("Cd4",),
    "ADT-CD5":                            ("Cd5",),
    "ADT-CD8a":                           ("Cd8a", "Cd8b1"),
    "ADT-CD25":                           ("Il2ra",),
    "ADT-CD27":                           ("Cd27",),  # via 'anti_mouse_rat_human_CD27' (cleaned to human_CD27)
    "ADT-anti_mouse_rat_human_CD27":      ("Cd27",),
    "ADT-CD28":                           ("Cd28",),
    "ADT-CD45":                           ("Ptprc",),
    "ADT-CD45_1":                         ("Ptprc",),  # Ly5.1 allele
    "ADT-CD45_2":                         ("Ptprc",),  # Ly5.2 allele
    "ADT-CD69":                           ("Cd69",),
    "ADT-CD90_2":                         ("Thy1",),
    "ADT-CD127":                          ("Il7r",),
    "ADT-CD159a_NKG2AB6":                 ("Klrc1",),  # NKG2A; antibody clone is NKG2A-specific despite NKG2AB6 label
    "ADT-CD223_LAG_3":                    ("Lag3",),
    "ADT-CD279_PD_1":                     ("Pdcd1",),
    "ADT-CD274_B7_H1_PD_L1":              ("Cd274",),
    "ADT-CD357_GITR":                     ("Tnfrsf18",),
    "ADT-CD95_Fas":                       ("Fas",),
    "ADT-Ly108":                          ("Slamf6",),
    "ADT-NK_1_1":                         ("Klrb1c",),
    "ADT-Ly_49A":                         ("Klra1",),
    "ADT-CD183_CXCR3":                    ("Cxcr3",),
    "ADT-CD182_CXCR2":                    ("Cxcr2",),
    "ADT-CCR3_CD193":                     ("Ccr3",),
    "ADT-CX3CR1":                         ("Cx3cr1",),
    "ADT-CD23":                           ("Fcer2a",),
    "ADT-CD22":                           ("Cd22",),

    # B cell core
    "ADT-CD19":                           ("Cd19",),
    "ADT-CD20":                           ("Ms4a1",),
    "ADT-IgD":                            ("Ighd",),
    "ADT-IgM":                            ("Ighm",),
    "ADT-CD79b_Ig":                       ("Cd79b",),
    "ADT-CD21_CD35_CR2_CR1":              ("Cr2",),  # mouse Cr2 encodes both CD21 and CD35 by alt splicing
    "ADT-anti_mouse_human_CD45R_B220":    ("Ptprc",),  # B220 = CD45R isoform of Ptprc
    "ADT-CD93_AA4_1_early_B_lineage":     ("Cd93",),

    # Monocyte/myeloid/granulocyte
    "ADT-CD14":                           ("Cd14",),
    "ADT-anti_mouse_human_CD11b":         ("Itgam",),
    "ADT-CD11a":                          ("Itgal",),
    "ADT-CD11c":                          ("Itgax",),
    "ADT-CD16_32":                        ("Fcgr3", "Fcgr2b"),  # CD16=Fcgr3, CD32=Fcgr2b; Fcgr4 is "CD16-2" not bound by 2.4G2 clone
    "ADT-CD64":                           ("Fcgr1",),  # not in panel
    "ADT-CD163":                          ("Cd163",),
    "ADT-CD172a_SIRP":                    ("Sirpa",),
    "ADT-Ly_6C":                          ("Ly6c1", "Ly6c2"),
    "ADT-Ly_6G":                          ("Ly6g",),
    "ADT-F4_80":                          ("Adgre1",),
    "ADT-CD169_Siglec_1":                 ("Siglec1",),
    "ADT-CD170_Siglec_F":                 ("Siglecf",),
    "ADT-Siglec_H":                       ("Siglech",),
    "ADT-CD317_BST2_PDCA_1":              ("Bst2",),
    "ADT-CD371_CLEC12A":                  ("Clec12a",),
    "ADT-CD301b":                         ("Mgl2",),  # CD301b is Mgl2; CD301a is Mgl1
    "ADT-CD107a_LAMP_1":                  ("Lamp1",),
    "ADT-CD86":                           ("Cd86",),
    "ADT-CD40":                           ("Cd40",),
    "ADT-I_A_I_E":                        ("H2-Aa", "H2-Ab1", "H2-Eb1"),  # MHC class II
    "ADT-CD205_DEC_205":                  ("Ly75",),
    "ADT-CD200_OX2":                      ("Cd200",),
    "ADT-CD200R_OX2R":                    ("Cd200r1",),
    "ADT-CD270_HVEM":                     ("Tnfrsf14",),

    # HSC / progenitor / erythroid / megakaryocyte
    "ADT-CD117_c_kit":                    ("Kit",),
    "ADT-CD135":                          ("Flt3",),
    "ADT-Ly_6A_E_Sca_1":                  ("Ly6a", "Ly6e"),
    "ADT-CD150_SLAM":                     ("Slamf1",),
    "ADT-CD34":                           ("Cd34",),  # not in panel but reserved
    "ADT-CD48":                           ("Cd48",),
    "ADT-CD201_EPCR":                     ("Procr",),
    "ADT-CD140a":                         ("Pdgfra",),
    "ADT-CD41":                           ("Itga2b",),
    "ADT-CD9":                            ("Cd9",),
    "ADT-CD71":                           ("Tfrc",),
    "ADT-TER_119_Erythroid_Cells":        ("Gypa",),  # Ter-119 epitope on glycophorin-A complex; Ly76 retired in current annotation
    "ADT-CD43":                           ("Spn",),
    "ADT-CD63":                           ("Cd63",),

    # Endothelial / stromal
    "ADT-CD102":                          ("Icam2",),
    "ADT-CD105":                          ("Eng",),
    "ADT-CD106":                          ("Vcam1",),
    "ADT-CD31":                           ("Pecam1",),
    "ADT-CD51":                           ("Itgav",),
    "ADT-CD54":                           ("Icam1",),
    "ADT-CD55_DAF":                       ("Cd55",),
    "ADT-ESAM":                           ("Esam",),
    "ADT-CD300LG_Nepmucin":               ("Cd300lg",),
    "ADT-CD155_PVR":                      ("Pvr",),

    # Integrins (other)
    "ADT-CD49a":                          ("Itga1",),
    "ADT-CD49b":                          ("Itga2",),
    "ADT-CD49d":                          ("Itga4",),
    "ADT-anti_human_mouse_CD49f":         ("Itga6",),
    "ADT-anti_human_mouse_integrin_7":    ("Itgb7",),
    "ADT-anti_mouse_rat_CD29":            ("Itgb1",),
    "ADT-anti_mouse_rat_CD61":            ("Itgb3",),
    "ADT-anti_mouse_rat_CD81":            ("Cd81",),

    # Adhesion / migration
    "ADT-CD62L":                          ("Sell",),
    "ADT-anti_mouse_human_CD44":          ("Cd44",),

    # Other / immunometabolism / inhibitory
    "ADT-CD24":                           ("Cd24a",),
    "ADT-CD26_DPP_4":                     ("Dpp4",),
    "ADT-CD36":                           ("Cd36",),
    "ADT-CD38":                           ("Cd38",),
    "ADT-CD73":                           ("Nt5e",),
    "ADT-CD85k_gp49_Receptor":            ("Lilr4b", "Lilrb4a"),  # gp49 = Lilrb4 family
    "ADT-CD1d_CD1_1_Ly_38":               ("Cd1d1",),  # antibody is Cd1d1-specific; Cd1d2 minor
    "ADT-PIR_A_B":                        ("Pira2", "Pirb"),  # paired Ig-like receptor A/B family
    "ADT-IL_33R__IL1RL1_ST2":             ("Il1rl1",),
    "ADT-Tim_4":                          ("Timd4",),

    # Isotype controls — intentionally empty (no RNA partner)
    "ADT-Mouse_IgG1___isotype_Ctrl":      (),
    "ADT-Mouse_IgG2a___isotype_Ctrl":     (),
    "ADT-Mouse_IgG2b___isotype_Ctrl":     (),
    "ADT-Rat_IgG1___Isotype_Ctrl":        (),
    "ADT-Rat_IgG1___isotype_Ctrl":        (),
    "ADT-Rat_IgG2a___Isotype_Ctrl":       (),
    "ADT-Rat_IgG2b___Isotype_Ctrl":       (),
    "ADT-Rat_IgG2c___Isotype_Ctrl":       (),
    "ADT-Armenian_Hamster_IgG_Isotype_Ctrl": (),
}


DEFAULT_MAP_TSV = Path(__file__).parent / "configs" / "adt_mgi_map.tsv"


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


def load_curated_adt_mgi_map(path: Path | None = None) -> Dict[str, List[str]]:
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
    if map_obj is None:
        map_obj = load_curated_adt_mgi_map()
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
            print(f"  {adt}\t(curated genes: {','.join(load_curated_adt_mgi_map().get(adt, []))})")
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
    parser.add_argument("--write-tsv", action="store_true", help="Write configs/adt_mgi_map.tsv")
    parser.add_argument("--audit-h5ad", type=Path, help="Audit the curated map against an .h5ad's var_names")
    args = parser.parse_args()
    if args.write_tsv:
        path = write_map_to_tsv()
        print(f"Wrote {path}")
    if args.audit_h5ad:
        import anndata as ad
        a = ad.read_h5ad(args.audit_h5ad, backed="r")
        var_names = [str(v) for v in a.var_names]
        adt_names = [v for v in var_names if v.startswith("ADT-")]
        rna_genes = [v for v in var_names if not v.startswith("ADT-")]
        report = audit_against_genes(adt_names, rna_genes)
        print_audit_report(report)
