"""Curated ADT -> HGNC gene map for the human lung (COVID TotalVI) CITE-seq panel.

The panel ships 56 antibodies whose column names carry a BioLegend catalogue
suffix (``CD235ab_A0196``). ``clean_adt_name`` strips that suffix and any
punctuation, and ``atlas_var_name`` adds the ``Hu.`` prefix the human bone
marrow atlas uses, so the lung bundle emits the same name style as
``rna2adt_bm_bundle.pkl``.

Every entry whose marker also exists in the bone marrow panel reuses the
partner gene(s) already curated in ``components.rna2adt.adt_rna_map`` so the
two human models agree. New entries name the HGNC symbol of the gene whose
protein the antibody binds; multi-subunit complexes list every chain.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


# Panel as released, in file-column order. Source of truth:
# /Users/saljh8/Dropbox/Transfer/COVID-TotalVI/denoised_ADT_HTC_covid_all_samples.txt
PANEL_RAW: Tuple[str, ...] = (
    "CD235ab_A0196", "CD45_A0391", "CD31_A0124", "CD144_A0400", "CD324_A0135",
    "CD326_A0123", "CD271_A0816", "CD49f_A0070", "CD66c & e_A0188", "CD55_A0383",
    "CD140a_A0128", "CD140b_A0129", "CD73_A0577", "CD90_A0060", "CD45RA_A0063",
    "CD19_A0050", "CD20_A0100", "CD138_A0831", "CD27_A0154", "CD38_A0410",
    "IgD_A0384", "IgM_A0136", "CD3_A0034", "CD8_A0046", "CD4_A0072",
    "CD127_A0390", "CD69_A0146", "CD103_A0145", "CD25_A0085", "CCR4_A0071",
    "CCR6_A0143", "CCR7_A0148", "CD161_A0149", "CXCR3_A0140", "CXCR5_A0144",
    "CD335_A0101", "CD117_A0061", "FceRIa_A0352", "CD15_A0392", "HLA-DR_A0159",
    "CD1c_A0160", "CD141_A0163", "CD303_A0370", "CD123_A064", "CD206_A0205",
    "CD169_A0206", "CD14_A0051", "CD274 (PD-L1)_A0007", "CD279 (PD-1)_A0088",
    "CD49a_A0575", "CD41_A0353", "Podoplanin_A0127", "VEGFR3_A0865",
    "CD186_A0804", "CCR5_A0141", "CD16_A0083",
)

_CATALOGUE_SUFFIX = re.compile(r"_A\d+$")
_PUNCTUATION = re.compile(r"[^A-Za-z0-9-]+")

ATLAS_PREFIX = "Hu."


def clean_adt_name(raw: str) -> str:
    """``"CD274 (PD-L1)_A0007"`` -> ``"CD274_PD-L1"``. Deterministic."""
    value = _CATALOGUE_SUFFIX.sub("", str(raw).strip())
    value = value.replace("(", " ").replace(")", " ")
    value = _PUNCTUATION.sub("_", value)
    return value.strip("_")


def atlas_var_name(raw: str) -> str:
    """Panel column -> the var name the training atlas and the bundle use."""
    return f"{ATLAS_PREFIX}{clean_adt_name(raw)}"


def strip_prefix(name: str) -> str:
    value = str(name).strip()
    return value[len(ATLAS_PREFIX):] if value.startswith(ATLAS_PREFIX) else value


# Keys are the atlas var names produced by ``atlas_var_name``.
# "BM" marks a partner set copied verbatim from components.rna2adt.adt_rna_map.
CURATED_MAP: Dict[str, Tuple[str, ...]] = {
    # --- Erythroid / platelet ---
    "Hu.CD235ab":   ("GYPA", "GYPB"),          # BM Hu.CD235a = GYPA; the ab clone binds A and B
    "Hu.CD41":      ("ITGA2B",),               # BM
    # --- Pan-leukocyte ---
    "Hu.CD45":      ("PTPRC",),                # BM Hu.CD45_2D1
    "Hu.CD45RA":    ("PTPRC",),                # BM; isoform of the same gene
    # --- Endothelium / lymphatics ---
    "Hu.CD31":      ("PECAM1",),
    "Hu.CD144":     ("CDH5",),
    "Hu.VEGFR3":    ("FLT4",),
    "Hu.Podoplanin": ("PDPN",),
    # --- Epithelium ---
    "Hu.CD324":     ("CDH1",),
    "Hu.CD326":     ("EPCAM",),                # BM
    "Hu.CD49f":     ("ITGA6",),
    "Hu.CD66c_e":   ("CEACAM6", "CEACAM5"),    # CD66c = CEACAM6, CD66e = CEACAM5
    # --- Mesenchyme ---
    "Hu.CD271":     ("NGFR",),                 # BM
    "Hu.CD140a":    ("PDGFRA",),
    "Hu.CD140b":    ("PDGFRB",),               # BM
    "Hu.CD73":      ("NT5E",),                 # BM
    "Hu.CD90":      ("THY1",),                 # BM
    "Hu.CD55":      ("CD55",),                 # BM
    # --- B / plasma ---
    "Hu.CD19":      ("CD19",),                 # BM
    "Hu.CD20":      ("MS4A1",),
    "Hu.CD138":     ("SDC1",),                 # BM Hu.CD138_DL.101
    "Hu.CD27":      ("CD27",),                 # BM
    "Hu.CD38":      ("CD38",),                 # BM Hu.CD38_HIT2
    "Hu.IgD":       ("IGHD",),
    "Hu.IgM":       ("IGHM",),
    # --- T / NK ---
    "Hu.CD3":       ("CD3D", "CD3E", "CD3G", "CD247"),   # BM
    "Hu.CD8":       ("CD8A", "CD8B"),          # BM
    "Hu.CD4":       ("CD4",),                  # BM
    "Hu.CD127":     ("IL7R",),                 # BM
    "Hu.CD69":      ("CD69",),                 # BM
    "Hu.CD103":     ("ITGAE",),                # BM
    "Hu.CD25":      ("IL2RA",),                # BM
    "Hu.CD161":     ("KLRB1",),
    "Hu.CD335":     ("NCR1",),                 # BM
    "Hu.CD279_PD-1": ("PDCD1",),               # BM Hu.CD279
    "Hu.CD274_PD-L1": ("CD274",),              # BM Hu.CD274
    "Hu.CD49a":     ("ITGA1",),
    # --- Chemokine receptors ---
    "Hu.CCR4":      ("CCR4",),
    "Hu.CCR5":      ("CCR5",),
    "Hu.CCR6":      ("CCR6",),
    "Hu.CCR7":      ("CCR7",),
    "Hu.CXCR3":     ("CXCR3",),                # BM Hu.CD183
    "Hu.CXCR5":     ("CXCR5",),                # BM Hu.CD185
    "Hu.CD186":     ("CXCR6",),                # BM
    # --- Myeloid / DC / mast ---
    "Hu.CD117":     ("KIT",),                  # BM
    "Hu.FceRIa":    ("FCER1A",),               # BM
    "Hu.CD15":      ("FUT4",),                 # BM Hu.CD15_W6D3
    "Hu.HLA-DR":    ("HLA-DRA", "HLA-DRB1", "HLA-DRB5"),  # BM Hu.HLA.DR.DP.DQ, DR chains only
    "Hu.CD1c":      ("CD1C",),
    "Hu.CD141":     ("THBD",),                 # BM
    "Hu.CD303":     ("CLEC4C",),
    "Hu.CD123":     ("IL3RA",),                # BM
    "Hu.CD206":     ("MRC1",),
    "Hu.CD169":     ("SIGLEC1",),
    "Hu.CD14":      ("CD14",),                 # BM Hu.CD14_M5E2
    "Hu.CD16":      ("FCGR3A", "FCGR3B"),      # BM
}


def load_curated_adt_rna_map() -> Dict[str, Tuple[str, ...]]:
    return dict(CURATED_MAP)


def panel_atlas_names() -> List[str]:
    return [atlas_var_name(name) for name in PANEL_RAW]


def write_map_to_tsv(path: Path | None = None) -> Path:
    path = path or (Path(__file__).parent / "configs" / "adt_rna_map.tsv")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["adt_raw", "adt_atlas", "adt_clean", "rna_genes"])
        for raw in PANEL_RAW:
            atlas = atlas_var_name(raw)
            writer.writerow([raw, atlas, strip_prefix(atlas),
                             ",".join(CURATED_MAP.get(atlas, ()))])
    return path


def audit_against_genes(
    adt_names: Sequence[str],
    rna_genes: Iterable[str],
    map_obj: Dict[str, Tuple[str, ...]] | None = None,
) -> Dict[str, Dict[str, object]]:
    curated = map_obj if map_obj is not None else load_curated_adt_rna_map()
    rna_set = set(str(gene) for gene in rna_genes)
    audit: Dict[str, Dict[str, object]] = {}
    for adt in adt_names:
        entry = curated.get(str(adt))
        if entry is None:
            audit[str(adt)] = {"matched": [], "missing": [], "unmapped": True}
            continue
        audit[str(adt)] = {
            "matched": [g for g in entry if g in rna_set],
            "missing": [g for g in entry if g not in rna_set],
            "unmapped": False,
        }
    return audit


def print_audit_report(audit: Dict[str, Dict[str, object]]) -> None:
    unmapped = sorted(a for a, i in audit.items() if i["unmapped"])
    no_match = sorted(a for a, i in audit.items() if not i["unmapped"] and not i["matched"])
    partial = sorted((a, i["matched"], i["missing"]) for a, i in audit.items()
                     if i["matched"] and i["missing"])
    full = sorted(a for a, i in audit.items()
                  if i["matched"] and not i["missing"] and not i["unmapped"])
    print(f"\n=== Fully matched ({len(full)}) ===")
    for adt in full:
        print(f"  {adt}\t{','.join(audit[adt]['matched'])}")
    if partial:
        print(f"\n=== Partially matched ({len(partial)}) ===")
        for adt, matched, missing in partial:
            print(f"  {adt}\tmatched={','.join(matched)}\tmissing={','.join(missing)}")
    if no_match:
        print(f"\n=== Mapped but NO gene in atlas ({len(no_match)}) ===")
        for adt in no_match:
            print(f"  {adt}\t(curated: {','.join(CURATED_MAP.get(adt, ()))})")
    if unmapped:
        print(f"\n=== Not in curated map ({len(unmapped)}) ===")
        for adt in unmapped:
            print(f"  {adt}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-tsv", action="store_true")
    parser.add_argument("--audit-h5ad", type=Path,
                        help="RNA .h5ad whose var_names are checked for every curated partner")
    args = parser.parse_args()
    if args.write_tsv:
        print(f"Wrote {write_map_to_tsv()}")
    if args.audit_h5ad:
        import h5py
        with h5py.File(args.audit_h5ad, "r") as handle:
            genes = [v.decode() if isinstance(v, bytes) else str(v)
                     for v in handle["var/_index"][:]]
        print_audit_report(audit_against_genes(panel_atlas_names(), genes))
