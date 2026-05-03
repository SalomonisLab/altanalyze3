"""
reference_selector.py: Automated APPRIS/MANE SELECT reference isoform selection.

When a user provides a single novel isoform from AltAnalyze3 long-read output,
Proteus needs to pair it with the correct reference (canonical) isoform for
that gene.  This module resolves that reference automatically.

Priority order:
  1. MANE SELECT (NCBI + Ensembl consensus canonical)
  2. APPRIS Principal 1
  3. APPRIS Principal 2 (if no P1)
  4. Longest CDS isoform (fallback)

Data sources (Daedalus Phase A interim):
  - gencode_protein_reference.tsv   (protein_id, transcript_id, gene_id, is_mane_select,
                                     appris_principal_score, protein_length)
  - gene_supervision_catalog.tsv    (gene_id, gene_name, mane_select_transcript_id,
                                     appris_principal_transcript_id)

Usage
-----
from proteus.data.reference_selector import ReferenceSelector

selector = ReferenceSelector()
selector.load("/path/to/daedalus/data/interim")

# Single transcript
ref = selector.get_reference("ENST00000504297.5")
# → {"gene_id": "ENSG00000146648.17", "ref_transcript_id": "ENST00000275493.7",
#    "ref_protein_id": "ENSP00000275493.3", "selection_method": "mane_select"}

# Batch (AltAnalyze3 junction TSV)
pairs = selector.build_pairs(novel_transcript_ids=["ENST...", "ENST..."])
# → list of dicts suitable for ProteusDataset
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional


class ReferenceSelector:
    """
    Selects the canonical reference isoform for any transcript in a gene.

    Stores three lookup tables built from Daedalus interim:
      _gene_of_transcript : transcript_id → gene_id
      _reference_of_gene  : gene_id → reference record (transcript_id, protein_id, method)
      _protein_of_transcript : transcript_id → protein_id

    Parameters
    ----------
    priority : list[str]
        Ordered list of selection methods to try.
        Default: ["mane_select", "appris_principal_1", "appris_principal_2", "longest_cds"]
    """

    def __init__(
        self,
        priority: Optional[List[str]] = None,
    ) -> None:
        self.priority = priority or [
            "mane_select",
            "appris_principal_1",
            "appris_principal_2",
            "longest_cds",
        ]
        self._gene_of_transcript: Dict[str, str] = {}        # ENST → ENSG
        self._protein_of_transcript: Dict[str, str] = {}     # ENST → ENSP
        self._reference_of_gene: Dict[str, Dict] = {}        # ENSG → ref record
        self._gene_name: Dict[str, str] = {}                  # ENSG → gene name
        self._loaded = False

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    def load(self, daedalus_interim_dir: Path) -> None:
        """
        Build lookup tables from Daedalus Phase A interim directory.

        Parameters
        ----------
        daedalus_interim_dir : Path
            Path to Daedalus data/interim/.
        """
        import pandas as pd

        d = Path(daedalus_interim_dir)

        # ---- Primary source: gencode_protein_reference.tsv ----
        gpr_path = d / "gencode_protein_reference.tsv"
        if not gpr_path.exists():
            warnings.warn(f"[ReferenceSelector] gencode_protein_reference.tsv not found in {d}. "
                          "Reference selection will only use gene_supervision_catalog.")
            gpr_df = None
        else:
            gpr_df = pd.read_csv(gpr_path, sep="\t", low_memory=False)

        # ---- Secondary source: gene_supervision_catalog.tsv ----
        sup_path = d / "gene_supervision_catalog.tsv"
        sup_df = None
        if sup_path.exists():
            sup_df = pd.read_csv(sup_path, sep="\t", low_memory=False)

        # Build transcript → gene + protein lookups from gencode_protein_reference
        if gpr_df is not None:
            for _, row in gpr_df.iterrows():
                tid = str(row.get("transcript_id", "") or "").strip()
                gid = str(row.get("gene_id", "") or "").strip()
                pid = str(row.get("protein_id", "") or "").strip()
                if tid and gid:
                    self._gene_of_transcript[tid] = gid
                if tid and pid:
                    self._protein_of_transcript[tid] = pid

            # Build gene → reference record
            # Prefer rows flagged is_mane_select=1, then by appris_principal_score
            for gid, grp in gpr_df.groupby("gene_id"):
                gid = str(gid)
                ref_record = self._select_from_group(grp)
                if ref_record:
                    self._reference_of_gene[gid] = ref_record

        # Override/supplement with gene_supervision_catalog
        if sup_df is not None:
            for _, row in sup_df.iterrows():
                gid = str(row.get("gene_id", "") or "").strip()
                gname = str(row.get("gene_name", "") or "").strip()
                if gid and gname:
                    self._gene_name[gid] = gname

                # If catalog has explicit MANE/APPRIS columns, use them
                mane_tid = str(row.get("mane_select_transcript_id", "") or "").strip()
                appris_tid = str(row.get("appris_principal_transcript_id", "") or "").strip()

                if gid and mane_tid and mane_tid not in ("", "nan", "None"):
                    pid = self._protein_of_transcript.get(mane_tid, "")
                    self._reference_of_gene[gid] = {
                        "gene_id": gid,
                        "ref_transcript_id": mane_tid,
                        "ref_protein_id": pid,
                        "selection_method": "mane_select",
                        "gene_name": gname,
                    }
                elif gid and appris_tid and appris_tid not in ("", "nan", "None"):
                    if gid not in self._reference_of_gene:
                        pid = self._protein_of_transcript.get(appris_tid, "")
                        self._reference_of_gene[gid] = {
                            "gene_id": gid,
                            "ref_transcript_id": appris_tid,
                            "ref_protein_id": pid,
                            "selection_method": "appris_principal_1",
                            "gene_name": gname,
                        }

        self._loaded = True
        print(
            f"[ReferenceSelector] Loaded: "
            f"{len(self._gene_of_transcript):,} transcripts mapped to genes, "
            f"{len(self._reference_of_gene):,} genes with reference isoforms."
        )

    def _select_from_group(self, grp) -> Optional[Dict]:
        """
        Given a DataFrame group for one gene, select the reference isoform
        according to the priority order.
        """
        # Try MANE SELECT first
        if "is_mane_select" in grp.columns:
            mane = grp[grp["is_mane_select"].astype(str).isin(("1", "True", "true"))]
            if not mane.empty:
                row = mane.iloc[0]
                return self._row_to_record(row, "mane_select")

        # Try APPRIS principal_1
        if "appris_principal_score" in grp.columns:
            try:
                scored = grp.copy()
                scored["_score"] = scored["appris_principal_score"].apply(
                    lambda x: self._appris_priority(str(x))
                )
                best = scored.sort_values("_score").iloc[0]
                score_val = best["_score"]
                if score_val <= 2:
                    method = "appris_principal_1" if score_val == 1 else "appris_principal_2"
                    return self._row_to_record(best, method)
            except Exception:
                pass

        # Fallback: longest CDS (largest protein_length)
        if "protein_length" in grp.columns:
            try:
                best = grp.sort_values("protein_length", ascending=False).iloc[0]
                return self._row_to_record(best, "longest_cds")
            except Exception:
                pass

        if len(grp) > 0:
            return self._row_to_record(grp.iloc[0], "first_available")

        return None

    @staticmethod
    def _appris_priority(label: str) -> int:
        """Map APPRIS label to sort key (lower = more canonical)."""
        mapping = {
            "principal1": 1, "principal_1": 1, "p1": 1,
            "principal2": 2, "principal_2": 2, "p2": 2,
            "principal3": 3, "principal_3": 3, "p3": 3,
            "principal4": 4, "principal_4": 4, "p4": 4,
            "principal5": 5, "principal_5": 5, "p5": 5,
            "alternative1": 6, "alternative2": 7,
        }
        return mapping.get(label.lower().replace(" ", "").replace("-", ""), 99)

    @staticmethod
    def _row_to_record(row, method: str) -> Dict:
        return {
            "gene_id": str(row.get("gene_id", "") or ""),
            "ref_transcript_id": str(row.get("transcript_id", "") or ""),
            "ref_protein_id": str(row.get("protein_id", "") or ""),
            "selection_method": method,
            "gene_name": str(row.get("gene_name", "") or ""),
        }

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_reference(self, transcript_id: str) -> Optional[Dict[str, str]]:
        """
        Return the reference isoform record for the gene containing transcript_id.

        Parameters
        ----------
        transcript_id : str
            Any transcript ID (reference OR alternative). The gene is looked up
            and the canonical reference for that gene is returned.

        Returns
        -------
        dict | None
            {"gene_id", "ref_transcript_id", "ref_protein_id",
             "selection_method", "gene_name"}
            Returns None if the gene cannot be found.
        """
        if not self._loaded:
            raise RuntimeError("Call .load(daedalus_interim_dir) first.")

        # Try direct lookup
        gene_id = self._gene_of_transcript.get(transcript_id, "")

        # Try base ID (strip version suffix)
        if not gene_id:
            base = transcript_id.split(".")[0]
            for tid, gid in self._gene_of_transcript.items():
                if tid.split(".")[0] == base:
                    gene_id = gid
                    break

        if not gene_id:
            return None

        return self._reference_of_gene.get(gene_id)

    def get_reference_by_gene(self, gene_id: str) -> Optional[Dict[str, str]]:
        """Return reference isoform for a gene ID directly."""
        if not self._loaded:
            raise RuntimeError("Call .load(daedalus_interim_dir) first.")
        return self._reference_of_gene.get(gene_id)

    def build_pairs(
        self,
        novel_transcript_ids: List[str],
        novel_protein_ids: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """
        Build reference-alternative pairs for a list of novel transcript IDs.

        Parameters
        ----------
        novel_transcript_ids : list[str]
            Novel alternative transcript IDs from AltAnalyze3 output.
        novel_protein_ids : list[str] | None
            Corresponding protein IDs. If None, left empty.

        Returns
        -------
        list[dict]
            Each dict has: pair_id, reference_transcript_id, reference_protein_id,
            alternative_transcript_id, alternative_protein_id, gene_id,
            gene_name, selection_method.
            Pairs where the reference could not be found are omitted
            (a warning is printed).
        """
        if not self._loaded:
            raise RuntimeError("Call .load(daedalus_interim_dir) first.")

        pairs = []
        n_failed = 0
        pids = novel_protein_ids or [""] * len(novel_transcript_ids)

        for alt_tid, alt_pid in zip(novel_transcript_ids, pids):
            ref = self.get_reference(alt_tid)
            if ref is None:
                n_failed += 1
                continue

            # Skip self-pairing (novel transcript IS the reference)
            if ref["ref_transcript_id"] == alt_tid:
                continue

            pairs.append({
                "pair_id": f"{ref['ref_transcript_id']}_vs_{alt_tid}",
                "reference_transcript_id": ref["ref_transcript_id"],
                "reference_protein_id": ref["ref_protein_id"],
                "alternative_transcript_id": alt_tid,
                "alternative_protein_id": alt_pid,
                "gene_id": ref["gene_id"],
                "gene_name": ref["gene_name"],
                "selection_method": ref["selection_method"],
            })

        if n_failed > 0:
            warnings.warn(
                f"[ReferenceSelector] {n_failed} novel transcripts could not be "
                f"mapped to a gene/reference — excluded from pairs."
            )

        method_counts: Dict[str, int] = {}
        for p in pairs:
            m = p["selection_method"]
            method_counts[m] = method_counts.get(m, 0) + 1

        print(f"[ReferenceSelector] Built {len(pairs):,} pairs. "
              f"Method breakdown: {method_counts}")
        return pairs

    def to_dataframe(self, pairs: List[Dict[str, Any]]):
        """Convert list of pair dicts to pandas DataFrame."""
        import pandas as pd
        return pd.DataFrame(pairs)
