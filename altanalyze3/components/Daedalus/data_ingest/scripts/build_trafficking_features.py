#!/usr/bin/env python3
"""Biology-aware features for cell-surface trafficking prediction.

For each catalog isoform we emit features that capture the structural
arithmetic of trafficking failure modes — signal peptide presence and
integrity, transmembrane topology coherence, ER export / retention motif
balance in the cytoplasmic tail, glycosylation/disulfide retention in the
extracellular domain, and PDZ/GPI/furin/polybasic signal status. The
intent is to give the autoencoder a continuous structural manifold to
compress, rather than count statistics dominated by axis-aligned splits.

Inputs:
    data_ingest/data/interim/uniprot_isoform_proteins.tsv
    data_ingest/data/interim/uniprot_features.tsv

Output:
    data_ingest/data/interim/uniprot_isoform_trafficking_features.tsv
"""
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

INTERIM_DIR = Path(__file__).resolve().parents[1] / "data" / "interim"
PROT_PATH = INTERIM_DIR / "uniprot_isoform_proteins.tsv"
FEAT_PATH = INTERIM_DIR / "uniprot_features.tsv"
OUT_PATH = INTERIM_DIR / "uniprot_isoform_trafficking_features.tsv"

# ---------- Trafficking motif library (all matched on the protein sequence) ----------
# ER export motifs (forward trafficking from ER -> Golgi -> PM)
RE_DI_ACIDIC = re.compile(r"[DE].{0,2}[DE]")          # generic di-acidic export
RE_FCYENEV = re.compile(r"FCYENEV")                   # KIR-style ER export
RE_FXYXXL = re.compile(r"F.Y..[ILV]")                 # FxYxxΦ ER export

# ER retention / retrieval signals (block PM delivery)
RE_KDEL_TAIL = re.compile(r"KDEL$")                   # soluble ER retention
RE_KKXX_TAIL = re.compile(r"KK..$")                   # type-I ER retrieval
RE_RXR_RETENTION = re.compile(r"R.R")                 # Kir/RxR retention motif

# Endolysosomal sorting (degradation pathway)
RE_YXXPHI = re.compile(r"Y..[AILMFWV]")               # YxxΦ tyrosine-based
RE_DI_LEUCINE = re.compile(r"[DE].{2,3}LL")           # [DE]xxxLL or [DE]xxLL

# PDZ-binding C-terminal classes
RE_PDZ_CLASS_I = re.compile(r"[ST].[ILV]$")           # class I: S/T-X-Φ
RE_PDZ_CLASS_II = re.compile(r"[YFW].[ILV]$")         # class II: Φ-X-Φ

# Polybasic / GPI / furin
RE_POLYBASIC_TAIL = re.compile(r"[KR]{4,}$")
RE_FURIN = re.compile(r"[KR].{0,2}[KR]R")             # PC convertase recognition
RE_GPI_OMEGA = re.compile(r"[GASCND].{0,15}$")        # small omega residue + hydrophobic stretch downstream
RE_HYDROPHOBIC_RUN = re.compile(r"[AILMFWVY]{15,}")   # contiguous hydrophobic stretch (rough GPI / TM)

CYTO_TAIL_LEN = 50                                    # last N residues = "cytoplasmic tail"
NEAR_TM_WINDOW = 10                                   # within N residues of last TM


def _parse_alt_seg(text: str) -> list[tuple[int, int]]:
    """Parse the catalog's pipe-separated alt_seq_changes string into
    (begin, end) coordinates on the canonical sequence that are
    deleted/replaced in the alt isoform."""
    if pd.isna(text) or not text:
        return []
    out: list[tuple[int, int]] = []
    for piece in str(text).split("|"):
        piece = piece.strip()
        if not piece:
            continue
        head = piece.split(":", 1)[0].strip()
        if "-" in head:
            try:
                a, b = head.split("-", 1)
                out.append((int(a), int(b)))
            except ValueError:
                continue
    return out


def _segments_overlap(a: tuple[int, int], b_list: list[tuple[int, int]]) -> bool:
    a0, a1 = a
    return any(not (b1 < a0 or b0 > a1) for b0, b1 in b_list)


def _segment_intersection_len(a: tuple[int, int], b: tuple[int, int]) -> int:
    a0, a1 = a; b0, b1 = b
    lo = max(a0, b0); hi = min(a1, b1)
    return max(0, hi - lo + 1)


def _signal_peptide_range(features: list[tuple[str, str, int, int]]) -> tuple[int, int] | None:
    for ftype, _desc, b, e in features:
        if ftype == "Signal":
            return (b, e)
    return None


def _tm_segments(features: list[tuple[str, str, int, int]]) -> list[tuple[int, int]]:
    return [(b, e) for ftype, _desc, b, e in features if ftype in {"Transmembrane", "Intramembrane"}]


def _topology_segments(features: list[tuple[str, str, int, int]], kind: str) -> list[tuple[int, int]]:
    # kind in {"extracellular", "cytoplasmic"}
    if kind == "extracellular":
        keys = ("extracellular", "luminal", "lumenal", "external")
    else:
        keys = ("cytoplasmic", "cytosolic", "intracellular")
    out = []
    for ftype, desc, b, e in features:
        if ftype != "Topological domain":
            continue
        d = (desc or "").lower()
        if any(k in d for k in keys):
            out.append((b, e))
    return out


def _glyco_sites(features: list[tuple[str, str, int, int]]) -> list[int]:
    return [b for ftype, _desc, b, e in features if ftype == "Glycosylation"]


def _disulfide_pairs(features: list[tuple[str, str, int, int]]) -> list[tuple[int, int]]:
    # In UniProt features.tsv each disulfide is encoded as one row with begin/end
    # = the two cysteine positions.
    return [(b, e) for ftype, _desc, b, e in features if ftype == "Disulfide bond"]


def _propeptide_ranges(features: list[tuple[str, str, int, int]]) -> list[tuple[int, int]]:
    return [(b, e) for ftype, _desc, b, e in features if ftype == "Propeptide"]


def _chain_ranges(features: list[tuple[str, str, int, int]]) -> list[tuple[int, int]]:
    return [(b, e) for ftype, _desc, b, e in features if ftype == "Chain"]


def _last_tm_end_in_alt(canonical_tms: list[tuple[int, int]],
                         deleted: list[tuple[int, int]],
                         alt_len: int) -> int | None:
    """Return the alt-coord position of the last TM that survives. If no TM
    survives, return None."""
    if not canonical_tms:
        return None
    surviving = [tm for tm in canonical_tms if not _segments_overlap(tm, deleted)]
    if not surviving:
        return None
    canonical_last = max(surviving, key=lambda x: x[1])
    # Map canonical end -> alt-coord end by subtracting deleted residues before it
    deleted_before = sum(_segment_intersection_len((1, canonical_last[1]), d) for d in deleted)
    alt_pos = canonical_last[1] - deleted_before
    return min(alt_pos, alt_len)


def _features_for_isoform(
    iso_row: pd.Series,
    features: list[tuple[str, str, int, int]],
    canonical_seq: str,
) -> dict:
    """Compute trafficking features for a single isoform (canonical or alt).

    For canonical, ``deleted = []`` and we evaluate the feature panel directly
    on the canonical sequence + features. For alt isoforms we compute the
    same panel but masked by the splice-deleted segments and on the alt's
    own reconstructed protein sequence.
    """
    is_canonical = int(iso_row["is_canonical"]) == 1
    deleted = [] if is_canonical else _parse_alt_seg(iso_row.get("alt_seq_changes", ""))
    seq = str(iso_row.get("protein_sequence", "") or "")
    seq_len = len(seq)
    if not seq_len:
        return {}

    sp_range = _signal_peptide_range(features)
    tms = _tm_segments(features)
    extr = _topology_segments(features, "extracellular")
    cyto = _topology_segments(features, "cytoplasmic")
    glyco_sites = _glyco_sites(features)
    disulf_pairs = _disulfide_pairs(features)
    propeps = _propeptide_ranges(features)
    chains = _chain_ranges(features)

    # ---- signal peptide integrity ----
    if sp_range is not None:
        sp_intact = int(not _segments_overlap(sp_range, deleted))
        sp_partial_aa = sum(_segment_intersection_len(sp_range, d) for d in deleted)
        sp_total_len = sp_range[1] - sp_range[0] + 1
        sp_retained_frac = max(0.0, (sp_total_len - sp_partial_aa) / max(sp_total_len, 1))
    else:
        sp_intact = 0
        sp_partial_aa = 0
        sp_retained_frac = 0.0

    # ---- TM topology integrity ----
    tm_lost = sum(1 for tm in tms if _segments_overlap(tm, deleted)
                  and all(_segment_intersection_len(tm, d) >= (tm[1]-tm[0]+1)*0.5 for d in deleted if _segment_intersection_len(tm, d) > 0))
    # simpler: count TMs >50% overlapped by any single deletion
    tm_lost = 0; tm_partial = 0
    for tm in tms:
        max_overlap = max((_segment_intersection_len(tm, d) for d in deleted), default=0)
        tm_len = tm[1] - tm[0] + 1
        if max_overlap >= tm_len * 0.5:
            tm_lost += 1
        elif max_overlap > 0:
            tm_partial += 1
    tm_remaining = len(tms) - tm_lost
    canonical_tm_count = len(tms)
    topology_parity_violated = int(canonical_tm_count > 0 and (tm_remaining % 2) != (canonical_tm_count % 2))

    extr_lost_aa = sum(_segment_intersection_len(seg, d) for seg in extr for d in deleted)
    cyto_lost_aa = sum(_segment_intersection_len(seg, d) for seg in cyto for d in deleted)
    extr_total_aa = sum(e - b + 1 for b, e in extr) or 1
    cyto_total_aa = sum(e - b + 1 for b, e in cyto) or 1
    extr_retained_frac = max(0.0, (extr_total_aa - extr_lost_aa) / extr_total_aa)
    cyto_retained_frac = max(0.0, (cyto_total_aa - cyto_lost_aa) / cyto_total_aa)

    # ---- glycosylation in extracellular domain ----
    glyco_in_extr = sum(1 for site in glyco_sites if any(b <= site <= e for b, e in extr))
    glyco_in_extr_lost = sum(1 for site in glyco_sites if any(b <= site <= e for b, e in extr)
                              and any(d0 <= site <= d1 for d0, d1 in deleted))
    glyco_extr_retained_frac = (
        (glyco_in_extr - glyco_in_extr_lost) / glyco_in_extr if glyco_in_extr else 1.0
    )

    # ---- disulfide partner orphaning ----
    disulf_total = len(disulf_pairs)
    disulf_orphaned = 0
    for c1, c2 in disulf_pairs:
        c1_lost = any(d0 <= c1 <= d1 for d0, d1 in deleted)
        c2_lost = any(d0 <= c2 <= d1 for d0, d1 in deleted)
        if c1_lost ^ c2_lost:
            disulf_orphaned += 1
    disulf_orphaned_frac = (disulf_orphaned / disulf_total) if disulf_total else 0.0

    # ---- propeptide / chain integrity ----
    propep_intact = int(bool(propeps) and not any(_segments_overlap(p, deleted) for p in propeps))
    chain_disrupted = int(any(_segments_overlap(c, deleted) for c in chains)) if chains else 0

    # ---- furin cleavage site ----
    furin_canonical = bool(RE_FURIN.search(canonical_seq)) if canonical_seq else False
    furin_in_seq = bool(RE_FURIN.search(seq))
    furin_intact = int(not furin_canonical or furin_in_seq)

    # ---- cytoplasmic tail / sorting motifs ----
    # Define cytoplasmic tail as the C-terminal CYTO_TAIL_LEN residues
    tail = seq[-CYTO_TAIL_LEN:]
    # Near-TM window: residues immediately after the last surviving TM
    last_tm_alt_end = _last_tm_end_in_alt(tms, deleted, seq_len) if not is_canonical else (max(tm[1] for tm in tms) if tms else None)
    if last_tm_alt_end is not None and last_tm_alt_end < seq_len:
        near_tm = seq[last_tm_alt_end: last_tm_alt_end + NEAR_TM_WINDOW]
    else:
        near_tm = ""

    # Counts on alt sequence
    di_acidic_tail = len(RE_DI_ACIDIC.findall(tail))
    fcyenev_seq = int(bool(RE_FCYENEV.search(seq)))
    fxyxxl_seq = int(bool(RE_FXYXXL.search(seq)))
    kdel_tail = int(bool(RE_KDEL_TAIL.search(seq)))
    kkxx_tail = int(bool(RE_KKXX_TAIL.search(seq)))
    rxr_tail = int(bool(RE_RXR_RETENTION.search(tail)))
    yxxphi_near_tm = int(bool(RE_YXXPHI.search(near_tm))) if near_tm else 0
    dileucine_tail = int(bool(RE_DI_LEUCINE.search(tail)))
    pdz_class_i = int(bool(RE_PDZ_CLASS_I.search(seq)))
    pdz_class_ii = int(bool(RE_PDZ_CLASS_II.search(seq)))
    polybasic_tail = int(bool(RE_POLYBASIC_TAIL.search(tail)))
    hydrophobic_run_cterm = int(bool(RE_HYDROPHOBIC_RUN.search(seq[-30:])))

    # Canonical-ref motif presence (for "lost" / "gained" booleans on alts)
    if is_canonical:
        ref_kdel = kdel_tail
        ref_kkxx = kkxx_tail
        ref_rxr = rxr_tail
        ref_yxxphi_neartm = yxxphi_near_tm
        ref_dileucine = dileucine_tail
        ref_pdzI = pdz_class_i
        ref_pdzII = pdz_class_ii
        ref_polybasic = polybasic_tail
        ref_hydro_cterm = hydrophobic_run_cterm
        ref_fcyenev = fcyenev_seq
        ref_fxyxxl = fxyxxl_seq
    else:
        ref_tail = canonical_seq[-CYTO_TAIL_LEN:] if canonical_seq else ""
        ref_kdel = int(bool(RE_KDEL_TAIL.search(canonical_seq))) if canonical_seq else 0
        ref_kkxx = int(bool(RE_KKXX_TAIL.search(canonical_seq))) if canonical_seq else 0
        ref_rxr = int(bool(RE_RXR_RETENTION.search(ref_tail))) if ref_tail else 0
        ref_dileucine = int(bool(RE_DI_LEUCINE.search(ref_tail))) if ref_tail else 0
        ref_pdzI = int(bool(RE_PDZ_CLASS_I.search(canonical_seq))) if canonical_seq else 0
        ref_pdzII = int(bool(RE_PDZ_CLASS_II.search(canonical_seq))) if canonical_seq else 0
        ref_polybasic = int(bool(RE_POLYBASIC_TAIL.search(ref_tail))) if ref_tail else 0
        ref_hydro_cterm = int(bool(RE_HYDROPHOBIC_RUN.search(canonical_seq[-30:]))) if canonical_seq else 0
        ref_fcyenev = int(bool(RE_FCYENEV.search(canonical_seq))) if canonical_seq else 0
        ref_fxyxxl = int(bool(RE_FXYXXL.search(canonical_seq))) if canonical_seq else 0
        # Last TM near-window in canonical
        if tms:
            cano_last_tm_end = max(tm[1] for tm in tms)
            cano_near_tm = canonical_seq[cano_last_tm_end: cano_last_tm_end + NEAR_TM_WINDOW]
            ref_yxxphi_neartm = int(bool(RE_YXXPHI.search(cano_near_tm))) if cano_near_tm else 0
        else:
            ref_yxxphi_neartm = 0

    out = {
        # Signal peptide integrity
        "sp_intact": sp_intact,
        "sp_retained_frac": round(sp_retained_frac, 4),
        # TM topology
        "tm_segments_lost": tm_lost,
        "tm_segments_partial": tm_partial,
        "tm_segments_retained": tm_remaining,
        "topology_parity_violated": topology_parity_violated,
        "extracellular_retained_frac": round(extr_retained_frac, 4),
        "cytoplasmic_retained_frac": round(cyto_retained_frac, 4),
        "extracellular_lost_aa": extr_lost_aa,
        "cytoplasmic_lost_aa": cyto_lost_aa,
        # Folding QC
        "glyco_extracellular_retained_frac": round(glyco_extr_retained_frac, 4),
        "disulfide_orphaned_count": disulf_orphaned,
        "disulfide_orphaned_frac": round(disulf_orphaned_frac, 4),
        "propeptide_intact": propep_intact,
        "chain_disrupted": chain_disrupted,
        "furin_intact": furin_intact,
        # Trafficking motifs (alt sequence)
        "di_acidic_tail_count": di_acidic_tail,
        "fcyenev_present": fcyenev_seq,
        "fxyxxl_present": fxyxxl_seq,
        "kdel_tail_present": kdel_tail,
        "kkxx_tail_present": kkxx_tail,
        "rxr_retention_tail_present": rxr_tail,
        "yxxphi_near_tm_present": yxxphi_near_tm,
        "dileucine_tail_present": dileucine_tail,
        "pdz_class_I_present": pdz_class_i,
        "pdz_class_II_present": pdz_class_ii,
        "polybasic_tail_present": polybasic_tail,
        "hydrophobic_run_cterm_present": hydrophobic_run_cterm,
        # Lost / gained relative to canonical (0 for canonical itself)
        "kdel_lost_vs_canonical": int(ref_kdel and not kdel_tail),
        "kdel_gained_vs_canonical": int(kdel_tail and not ref_kdel),
        "kkxx_lost_vs_canonical": int(ref_kkxx and not kkxx_tail),
        "kkxx_gained_vs_canonical": int(kkxx_tail and not ref_kkxx),
        "rxr_lost_vs_canonical": int(ref_rxr and not rxr_tail),
        "rxr_gained_vs_canonical": int(rxr_tail and not ref_rxr),
        "yxxphi_neartm_lost_vs_canonical": int(ref_yxxphi_neartm and not yxxphi_near_tm),
        "yxxphi_neartm_gained_vs_canonical": int(yxxphi_near_tm and not ref_yxxphi_neartm),
        "dileucine_lost_vs_canonical": int(ref_dileucine and not dileucine_tail),
        "dileucine_gained_vs_canonical": int(dileucine_tail and not ref_dileucine),
        "pdzI_lost_vs_canonical": int(ref_pdzI and not pdz_class_i),
        "pdzI_gained_vs_canonical": int(pdz_class_i and not ref_pdzI),
        "pdzII_lost_vs_canonical": int(ref_pdzII and not pdz_class_ii),
        "pdzII_gained_vs_canonical": int(pdz_class_ii and not ref_pdzII),
        "polybasic_lost_vs_canonical": int(ref_polybasic and not polybasic_tail),
        "polybasic_gained_vs_canonical": int(polybasic_tail and not ref_polybasic),
        "hydro_cterm_lost_vs_canonical": int(ref_hydro_cterm and not hydrophobic_run_cterm),
        "hydro_cterm_gained_vs_canonical": int(hydrophobic_run_cterm and not ref_hydro_cterm),
        "fcyenev_lost_vs_canonical": int(ref_fcyenev and not fcyenev_seq),
        "fxyxxl_lost_vs_canonical": int(ref_fxyxxl and not fxyxxl_seq),
    }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=OUT_PATH)
    args = ap.parse_args()

    proteins = pd.read_csv(PROT_PATH, sep="\t")
    print(f"[info] {len(proteins)} isoforms in proteins TSV")

    # Index canonical sequences by (species, gene_name) for "vs canonical" calculations
    canonical_seqs: dict[tuple[str, str], str] = {}
    for r in proteins[proteins["is_canonical"] == 1].itertuples(index=False):
        canonical_seqs[(str(r.species), str(r.gene_name))] = str(r.protein_sequence or "")

    # Cache canonical-coordinate features per accession
    feat_df = pd.read_csv(FEAT_PATH, sep="\t", low_memory=False)
    features_by_acc: dict[str, list[tuple[str, str, int, int]]] = defaultdict(list)
    for r in feat_df.itertuples(index=False):
        try:
            b = int(r.begin) if pd.notna(r.begin) else -1
            e = int(r.end) if pd.notna(r.end) else -1
        except (TypeError, ValueError):
            continue
        if b <= 0 or e < b:
            continue
        features_by_acc[str(r.primary_accession)].append(
            (str(r.feature_type or ""), str(r.description or ""), b, e)
        )
    print(f"[info] feature cache built for {len(features_by_acc)} accessions")

    rows: list[dict] = []
    for r in proteins.itertuples(index=False):
        feats = features_by_acc.get(str(r.primary_accession), [])
        canonical = canonical_seqs.get((str(r.species), str(r.gene_name)), "")
        out = _features_for_isoform(
            pd.Series({
                "is_canonical": r.is_canonical,
                "alt_seq_changes": r.alt_seq_changes,
                "protein_sequence": r.protein_sequence,
            }),
            feats, canonical,
        )
        if not out:
            continue
        rows.append({
            "species": str(r.species),
            "gene_name": str(r.gene_name),
            "primary_accession": str(r.primary_accession),
            "isoform_id": str(r.isoform_id),
            "is_canonical": int(r.is_canonical),
            **out,
        })

    out_df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    n_feat = len(out_df.columns) - 5
    print(f"[ok] wrote {args.out}  rows={len(out_df)}  trafficking_features={n_feat}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
