"""Structure-keyed functional annotation — the product layer.

Joins the resolved clone->structure crosswalk to every association (PPI partners, PDI DNA targets, M1H
activation call, condensate behavior) and to the authors' per-switch functional categories (Data_S8),
producing two structure-keyed tables that are the deliverable for AltAnalyze / ISV integration:

  - ``isoform_function.tsv`` : one row per resolved isoform (keyed on the structure string) with its
    gene, ENST/observed-isoform, NMD status, M1H call, and its detected PPI partners + PDI DNA targets
    (counts + lists), plus how many partners/targets were assayed (tested) so coverage is explicit.
  - ``switch_function.tsv``  : one row per reference->alternative isoform switch (both ends mapped to
    their structures) carrying the authors' PDI/PPI/M1H/localization categories and DBD loss.

Only binary (called) associations are used: Y2H for PPI, eY1H for PDI (the authors' Boolean calls);
M1H uses the verbatim >=1/<=-1 rule. N2H/luciferase scores are excluded (no cutoff). The clone label is
dropped from the identity — the structure string is the key, exactly as intended.
"""

import os
import logging

import pandas as pd

from . import config

logger = logging.getLogger(__name__)


def _load(data_dir, name):
    path = os.path.join(data_dir, f"{name}.tsv")
    if not os.path.exists(path):
        logger.warning("[iso2function.associate] missing %s (run ingest/crosswalk first)", path)
        return None
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])


def _detected_mask(df, col="detected"):
    return df[col].astype(str).str.strip().str.lower().isin(config.TRUE_TOKENS)


def build_isoform_function_table(data_dir=None):
    """Build the structure-keyed per-isoform functional annotation. Returns the DataFrame and writes
    ``isoform_function.tsv``. Requires ``clone_to_structure.tsv`` (run crosswalk first)."""
    data_dir = data_dir or config.DATA_DIR
    cts = _load(data_dir, "clone_to_structure")
    if cts is None or not len(cts):
        raise FileNotFoundError("clone_to_structure.tsv not found — run `iso2function crosswalk` first.")

    # ---- detected PPI partners per clone (Y2H, binary) ----
    y2h = _load(data_dir, "ppi_y2h")
    ppi_det, ppi_tested = {}, {}
    if y2h is not None and len(y2h):
        det = y2h[_detected_mask(y2h)]
        ppi_det = det.groupby("ad_clone_id")["db_gene_symbol"].apply(
            lambda s: sorted(set(x for x in s if x))).to_dict()
        ppi_tested = y2h.groupby("ad_clone_id")["db_gene_symbol"].nunique().to_dict()

    # ---- detected PDI DNA targets per clone (eY1H, binary) ----
    ey1h = _load(data_dir, "pdi_ey1h")
    pdi_det, pdi_tested = {}, {}
    if ey1h is not None and len(ey1h):
        det = ey1h[_detected_mask(ey1h)]
        pdi_det = det.groupby("clone_id")["bait_id"].apply(
            lambda s: sorted(set(x for x in s if x))).to_dict()
        pdi_tested = ey1h.groupby("clone_id")["bait_id"].nunique().to_dict()

    # ---- per-clone M1H call + condensate reference/alternative ----
    m1h = _load(data_dir, "activation_m1h")
    m1h_call = dict(zip(m1h["clone_id"], m1h.get("m1h_call", ""))) if m1h is not None else {}
    m1h_mean = dict(zip(m1h["clone_id"], m1h.get("M1H_mean", ""))) if m1h is not None else {}
    cond = _load(data_dir, "condensate")
    cond_refalt = (dict(zip(cond["clone_id"], cond.get("cloned_reference_or_alternative", "")))
                   if cond is not None else {})

    rows = []
    for _, r in cts.iterrows():
        cid = r["clone_id"]
        ppi = ppi_det.get(cid, [])
        pdi = pdi_det.get(cid, [])
        rows.append({
            "structure": r["structure"],
            "gene_symbol": r["gene_symbol"],
            "clone_id": cid,                              # provenance only; structure is the key
            "enst": r.get("enst", ""),
            "final_isoform_id": r.get("final_isoform_id", ""),
            "match_type": r.get("match_type", ""),
            "nmd_status": r.get("nmd_status", ""),
            "m1h_call": m1h_call.get(cid, ""),
            "m1h_mean": m1h_mean.get(cid, ""),
            "condensate_ref_alt": cond_refalt.get(cid, ""),
            "n_ppi_detected": len(ppi),
            "n_ppi_tested": int(ppi_tested.get(cid, 0)),
            "ppi_partners": ";".join(ppi),
            "n_pdi_detected": len(pdi),
            "n_pdi_tested": int(pdi_tested.get(cid, 0)),
            "pdi_targets": ";".join(pdi),
        })
    out = pd.DataFrame(rows)
    path = os.path.join(data_dir, "isoform_function.tsv")
    out.to_csv(path, sep="\t", index=False)
    logger.info("[iso2function.associate] isoform_function.tsv: %d isoforms (%d with >=1 PPI, "
                "%d with >=1 PDI, %d M1H-called)", len(out),
                int((out["n_ppi_detected"] > 0).sum()), int((out["n_pdi_detected"] > 0).sum()),
                int((out["m1h_call"].isin(["activator", "repressor"])).sum()))
    return out


def build_switch_table(data_dir=None):
    """Build the structure-keyed reference->alternative switch table from the authors' Data_S8
    categories, mapping both isoforms to their structures. Returns the DataFrame and writes
    ``switch_function.tsv``."""
    data_dir = data_dir or config.DATA_DIR
    fc = _load(data_dir, "functional_categories")
    cts = _load(data_dir, "clone_to_structure")
    if fc is None or cts is None:
        raise FileNotFoundError("need functional_categories.tsv + clone_to_structure.tsv (run ingest+crosswalk).")
    s = cts.drop_duplicates("clone_id").set_index("clone_id")
    struct = s["structure"].to_dict()
    enst = s["enst"].to_dict() if "enst" in s.columns else {}
    fid = s["final_isoform_id"].to_dict() if "final_isoform_id" in s.columns else {}

    rows = []
    for _, r in fc.iterrows():
        ref, alt = r.get("reference_isoform"), r.get("alternative_isoform")
        rows.append({
            "gene_symbol": r.get("gene_symbol"),
            "ref_clone": ref, "ref_structure": struct.get(ref, ""), "ref_enst": enst.get(ref, ""),
            "ref_final_isoform_id": fid.get(ref, ""),
            "alt_clone": alt, "alt_structure": struct.get(alt, ""), "alt_enst": enst.get(alt, ""),
            "alt_final_isoform_id": fid.get(alt, ""),
            "both_mapped": bool(struct.get(ref) and struct.get(alt)),
            "DBD_pct_lost_in_alt": r.get("DBD_pct_lost_in_alt", ""),
            "PDI_category": r.get("PDI_category", ""),
            "PPI_category": r.get("PPI_category", ""),
            "M1H_activation_category": r.get("M1H_activation_category", ""),
            "localization_category": r.get("localization_category", ""),
            "alt_iso_classification": r.get("alt_iso_classification", ""),
        })
    out = pd.DataFrame(rows)
    path = os.path.join(data_dir, "switch_function.tsv")
    out.to_csv(path, sep="\t", index=False)
    logger.info("[iso2function.associate] switch_function.tsv: %d switches (%d with both isoforms "
                "mapped to a structure)", len(out), int(out["both_mapped"].sum()))
    return out


def build_paper1_function(data_dir=None):
    """Structure/ENST-keyed paper1 per-isoform PPI annotation (measured, includes non-TF). Joins the
    paper1 PPI calls to the paper1->ENST crosswalk. Returns the DataFrame and writes
    ``paper1_isoform_function.tsv``. Requires ``paper1_ppi.tsv`` (ingest) + ``paper1_isoform_map.tsv``
    (crosswalk)."""
    data_dir = data_dir or config.DATA_DIR
    ppi = _load(data_dir, "paper1_ppi")
    m = _load(data_dir, "paper1_isoform_map")
    if ppi is None or m is None:
        raise FileNotFoundError("need paper1_ppi.tsv (ingest) + paper1_isoform_map.tsv (crosswalk).")
    mp = m.drop_duplicates("isoform_id").set_index("isoform_id").fillna("")
    det = ppi[ppi["detected"].astype(str).str.lower() == "true"]
    partners = det.groupby("Isoform_ID")["Interactor_Symbol"].apply(
        lambda s: sorted(set(x for x in s if x))).to_dict()
    tested = ppi.groupby("Isoform_ID")["Interactor_Symbol"].nunique().to_dict()
    rows = []
    for iid in ppi["Isoform_ID"].dropna().unique():
        meta = mp.loc[iid].to_dict() if iid in mp.index else {}
        plist = partners.get(iid, [])
        rows.append({
            "isoform_id": iid, "gene_symbol": meta.get("gene_symbol", iid.split("_")[0]),
            "enst": meta.get("enst", ""), "structure": meta.get("structure", ""),
            "identity": meta.get("identity", ""), "match_mode": meta.get("match_mode", "unmatched"),
            "source": "paper1", "evidence": "measured",
            "n_ppi_detected": len(plist), "n_ppi_tested": int(tested.get(iid, 0)),
            "ppi_partners": ";".join(plist),
        })
    out = pd.DataFrame(rows)
    path = os.path.join(data_dir, "paper1_isoform_function.tsv")
    out.to_csv(path, sep="\t", index=False)
    enst_mapped = int((out["enst"].astype(str).str.strip().str.startswith("ENST")).sum())
    logger.info("[iso2function.associate] paper1_isoform_function.tsv: %d isoforms (%d ENST-mapped, "
                "%d novel/unmapped, %d with >=1 detected PPI)", len(out), enst_mapped,
                len(out) - enst_mapped, int((out["n_ppi_detected"] > 0).sum()))
    return out


def build_all(data_dir=None):
    """Build the structure-keyed annotation tables (paper2 atlas + paper1). Returns
    (isoform_function, switch_function, paper1_isoform_function)."""
    iso = build_isoform_function_table(data_dir)
    sw = build_switch_table(data_dir)
    try:
        p1 = build_paper1_function(data_dir)
    except FileNotFoundError as e:
        logger.warning("[iso2function.associate] skipping paper1 product: %s", e)
        p1 = None
    return iso, sw, p1


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config.ensure_dirs()
    iso, sw = build_all()
    print(iso.head().to_string(index=False))
