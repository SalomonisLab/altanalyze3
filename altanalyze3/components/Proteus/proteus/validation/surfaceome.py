"""
surfaceome.py: Validation of novel isoform surface/membrane topology predictions.

Provides SurfaceomeValidator which cross-checks Proteus model predictions
against orthogonal lines of evidence that a protein is correctly folded
in the plasma membrane:

  1. Sequence-predicted TM helices  (SequenceTopologyEncoder)
  2. UniProt topology annotations   (extracellular_topology_count, tm_total_span)
  3. HPA extracellular annotation   (has_extracellular_annotation)
  4. BioGRID physical interactions  (membrane complex partners)
  5. Model surface_retained score   (from ProteusModel)

Each isoform receives a composite confidence score and a set of evidence flags.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np


@dataclass
class SurfaceEvidence:
    """
    Evidence record for a single isoform's surface/membrane status.

    Fields
    ------
    protein_id : str
        Ensembl protein ID of the alternative isoform.
    gene_name : str
    model_surface_score : float
        Raw logit from ProteusModel surface_retained head (higher = more likely surface).
    model_surface_prob : float
        Sigmoid-transformed probability.
    seq_tm_count : int
        Sequence-predicted TM helices for the alternative isoform.
    seq_signal_peptide : bool
    uniprot_tm_total_span : int
        UniProt-annotated TM span for the reference protein.
    uniprot_extracellular_count : int
        UniProt extracellular topology segment count.
    hpa_extracellular : bool
        HPA extracellular annotation present for this gene.
    biogrid_physical_partners : int
        Number of physical BioGRID interaction partners.
    confidence_score : float
        Composite score 0-1 (weighted sum of evidence lines).
    is_validated : bool
        True if confidence_score >= validation_threshold.
    evidence_flags : list[str]
        Human-readable list of positive evidence lines.
    """

    protein_id: str
    gene_name: str = ""
    model_surface_score: float = 0.0
    model_surface_prob: float = 0.0
    seq_tm_count: int = 0
    seq_signal_peptide: bool = False
    uniprot_tm_total_span: int = 0
    uniprot_extracellular_count: int = 0
    hpa_extracellular: bool = False
    biogrid_physical_partners: int = 0
    confidence_score: float = 0.0
    is_validated: bool = False
    evidence_flags: List[str] = field(default_factory=list)


class SurfaceomeValidator:
    """
    Validates novel isoform surface/membrane predictions using multi-source evidence.

    Usage
    -----
    validator = SurfaceomeValidator()
    validator.load(daedalus_interim_dir)

    # After running ProteusModel inference on novel isoforms:
    results = validator.validate(
        predictions=model_outputs,          # dict: protein_id → surface_logit
        sequences=novel_sequences,          # dict: protein_id → aa sequence
    )

    # Write to TSV
    validator.to_tsv(results, "surfaceome_validation.tsv")
    """

    # Weights for composite confidence score
    _WEIGHTS = {
        "model_prob": 0.35,         # strongest signal
        "seq_tm": 0.20,             # sequence-predicted TM
        "uniprot_tm": 0.15,         # UniProt TM annotation of reference
        "uniprot_extracell": 0.10,  # UniProt extracellular topology
        "hpa": 0.12,                # HPA membrane annotation
        "biogrid": 0.08,            # BioGRID membrane complex partners
    }

    def __init__(
        self,
        validation_threshold: float = 0.55,
    ) -> None:
        self.validation_threshold = validation_threshold
        self._region_priors: Dict[str, Dict] = {}   # ENSP → region_priors row
        self._hpa: Dict[str, Dict] = {}             # ENSP → HPA row
        self._biogrid: Dict[str, Dict] = {}         # ENSP → biogrid row
        self._gene_names: Dict[str, str] = {}       # ENSP → gene_name
        self._loaded = False

    def load(self, daedalus_interim_dir: Path) -> None:
        """
        Load reference tables from Daedalus Phase A interim directory.
        """
        import pandas as pd
        d = Path(daedalus_interim_dir)

        # Load gencode_protein_reference for ENSP → gene_name bridge
        gp_path = d / "gencode_protein_reference.tsv"
        if gp_path.exists():
            gp = pd.read_csv(gp_path, sep="\t", low_memory=False,
                             usecols=lambda c: c in ("gene_id", "protein_id"))
        else:
            gp = None

        # gene_supervision_catalog for gene_id → gene_name
        gene_id_to_name: Dict[str, str] = {}
        sup_path = d / "gene_supervision_catalog.tsv"
        if sup_path.exists():
            try:
                sup = pd.read_csv(sup_path, sep="\t", low_memory=False,
                                  usecols=lambda c: c in ("gene_id", "gene_name"))
                gene_id_to_name = dict(zip(sup["gene_id"].astype(str), sup["gene_name"].astype(str)))
            except Exception:
                pass

        # uniprot_region_priors → ENSP mapping
        rp_path = d / "uniprot_region_priors.tsv"
        map_path = d / "uniprot_gencode_map.tsv"
        if rp_path.exists() and map_path.exists():
            rp = pd.read_csv(rp_path, sep="\t", low_memory=False)
            gmap = pd.read_csv(map_path, sep="\t", low_memory=False,
                               usecols=lambda c: c in ("primary_accession",
                                                        "gencode_protein_id",
                                                        "gene_name"))
            acc_dict = {str(r["primary_accession"]): r.to_dict() for _, r in rp.iterrows()}
            for _, mrow in gmap.iterrows():
                acc = str(mrow.get("primary_accession", ""))
                pid = str(mrow.get("gencode_protein_id", ""))
                gn  = str(mrow.get("gene_name", ""))
                if acc in acc_dict and pid and pid != "nan":
                    self._region_priors[pid] = acc_dict[acc]
                    if gn and gn != "nan":
                        self._gene_names[pid] = gn

        # HPA → ENSP bridge
        hpa_path = d / "hpa_localization.tsv"
        if hpa_path.exists() and gp is not None:
            hpa = pd.read_csv(hpa_path, sep="\t", low_memory=False)
            hpa_by_gid: Dict[str, Dict] = {
                str(r["ensembl_gene_id"]): r.to_dict() for _, r in hpa.iterrows()
            }
            for _, row in gp.iterrows():
                pid = str(row.get("protein_id", ""))
                gid_versioned = str(row.get("gene_id", ""))
                gid_base = gid_versioned.split(".")[0]  # HPA uses unversioned IDs
                hpa_row = hpa_by_gid.get(gid_versioned) or hpa_by_gid.get(gid_base)
                if pid and hpa_row:
                    self._hpa[pid] = hpa_row
                if pid and gid_versioned and gid_versioned in gene_id_to_name and pid not in self._gene_names:
                    self._gene_names[pid] = gene_id_to_name[gid_versioned]

        # BioGRID → ENSP bridge
        bg_path = d / "biogrid_gene_interactions.tsv"
        if bg_path.exists():
            bg = pd.read_csv(bg_path, sep="\t", low_memory=False)
            bg_by_name = {str(r["gene_name"]): r.to_dict() for _, r in bg.iterrows()}
            for pid, gn in self._gene_names.items():
                if gn in bg_by_name:
                    self._biogrid[pid] = bg_by_name[gn]

        self._loaded = True
        print(f"[SurfaceomeValidator] Loaded: "
              f"region_priors={len(self._region_priors):,}, "
              f"HPA={len(self._hpa):,}, "
              f"BioGRID={len(self._biogrid):,} ENSP entries.")

    def _sigmoid(self, x: float) -> float:
        return float(1.0 / (1.0 + np.exp(-x)))

    def validate_one(
        self,
        protein_id: str,
        model_surface_logit: float,
        sequence: Optional[str] = None,
    ) -> SurfaceEvidence:
        """
        Produce a SurfaceEvidence record for a single protein.

        Parameters
        ----------
        protein_id : str
            Alternative isoform protein ID.
        model_surface_logit : float
            Raw logit from ProteusModel surface_retained head.
        sequence : str | None
            Amino acid sequence of the alternative isoform (for TM prediction).
        """
        from proteus.encoders.topology import SequenceTopologyEncoder, predict_tm_segments, predict_signal_peptide

        ev = SurfaceEvidence(protein_id=protein_id)
        ev.gene_name = self._gene_names.get(protein_id, "")
        ev.model_surface_score = model_surface_logit
        ev.model_surface_prob = self._sigmoid(model_surface_logit)

        # Sequence-based TM prediction for the novel isoform
        if sequence:
            tm_segs = predict_tm_segments(sequence)
            ev.seq_tm_count = len(tm_segs)
            ev.seq_signal_peptide = bool(predict_signal_peptide(sequence)[0])

        # UniProt reference topology
        rp = self._region_priors.get(protein_id, {})
        ev.uniprot_tm_total_span    = int(rp.get("tm_total_span", 0) or 0)
        ev.uniprot_extracellular_count = int(rp.get("extracellular_topology_count", 0) or 0)

        # HPA
        hpa = self._hpa.get(protein_id, {})
        ev.hpa_extracellular = bool(hpa.get("has_extracellular_annotation", False))

        # BioGRID
        bg = self._biogrid.get(protein_id, {})
        ev.biogrid_physical_partners = int(bg.get("biogrid_physical_partner_count", 0) or 0)

        # Composite confidence score
        weights = self._WEIGHTS
        score = weights["model_prob"] * ev.model_surface_prob
        score += weights["seq_tm"]        * min(ev.seq_tm_count / 4.0, 1.0)
        score += weights["uniprot_tm"]    * min(ev.uniprot_tm_total_span / 120.0, 1.0)
        score += weights["uniprot_extracell"] * min(ev.uniprot_extracellular_count / 3.0, 1.0)
        score += weights["hpa"]           * float(ev.hpa_extracellular)
        score += weights["biogrid"]       * min(ev.biogrid_physical_partners / 10.0, 1.0)
        ev.confidence_score = float(score)
        ev.is_validated = ev.confidence_score >= self.validation_threshold

        # Evidence flags
        flags: List[str] = []
        if ev.model_surface_prob >= 0.5:
            flags.append(f"model_surface_prob={ev.model_surface_prob:.3f}")
        if ev.seq_tm_count > 0:
            flags.append(f"seq_tm_helices={ev.seq_tm_count}")
        if ev.seq_signal_peptide:
            flags.append("seq_signal_peptide")
        if ev.uniprot_tm_total_span > 0:
            flags.append(f"uniprot_tm_span={ev.uniprot_tm_total_span}aa")
        if ev.uniprot_extracellular_count > 0:
            flags.append(f"uniprot_extracellular_domains={ev.uniprot_extracellular_count}")
        if ev.hpa_extracellular:
            flags.append("HPA_extracellular")
        if ev.biogrid_physical_partners > 0:
            flags.append(f"biogrid_physical_partners={ev.biogrid_physical_partners}")
        ev.evidence_flags = flags

        return ev

    def validate(
        self,
        predictions: Dict[str, float],
        sequences: Optional[Dict[str, str]] = None,
    ) -> List[SurfaceEvidence]:
        """
        Validate a batch of model surface predictions.

        Parameters
        ----------
        predictions : dict
            {protein_id: model_surface_logit}
        sequences : dict | None
            {protein_id: aa_sequence} for novel TM prediction.

        Returns
        -------
        list[SurfaceEvidence]
            Sorted by confidence_score descending.
        """
        if not self._loaded:
            raise RuntimeError("Call .load(daedalus_interim_dir) first.")

        results = []
        seqs = sequences or {}
        for pid, logit in predictions.items():
            ev = self.validate_one(pid, logit, seqs.get(pid))
            results.append(ev)

        results.sort(key=lambda r: r.confidence_score, reverse=True)
        return results

    def to_dataframe(self, results: List[SurfaceEvidence]):
        """Convert validation results to a pandas DataFrame."""
        import pandas as pd
        rows = []
        for r in results:
            rows.append({
                "protein_id": r.protein_id,
                "gene_name": r.gene_name,
                "model_surface_prob": round(r.model_surface_prob, 4),
                "seq_tm_count": r.seq_tm_count,
                "seq_signal_peptide": int(r.seq_signal_peptide),
                "uniprot_tm_total_span": r.uniprot_tm_total_span,
                "uniprot_extracellular_count": r.uniprot_extracellular_count,
                "hpa_extracellular": int(r.hpa_extracellular),
                "biogrid_physical_partners": r.biogrid_physical_partners,
                "confidence_score": round(r.confidence_score, 4),
                "is_validated": int(r.is_validated),
                "evidence_flags": "|".join(r.evidence_flags),
            })
        return pd.DataFrame(rows)

    def to_tsv(self, results: List[SurfaceEvidence], output_path: Path) -> None:
        """Write validation results to a TSV file."""
        df = self.to_dataframe(results)
        df.to_csv(output_path, sep="\t", index=False)
        n_validated = df["is_validated"].sum()
        print(f"[SurfaceomeValidator] {len(df)} isoforms evaluated, "
              f"{n_validated} validated (confidence ≥ {self.validation_threshold}).")
        print(f"  Written to: {output_path}")
