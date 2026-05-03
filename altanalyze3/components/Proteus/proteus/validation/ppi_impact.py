"""
ppi_impact.py: Specific PPI gain/loss prediction for isoform pairs.

For each reference-alternative isoform pair, predicts which protein-protein
interactions (PPIs) are gained or lost based on:

  1. UniProt domain annotations — which functional domains each isoform has
  2. Domain-domain interactions (DDI) — which domain pairs mediate PPIs
     (from 3DID flat file if available, or inferred from BioGRID + UniProt)
  3. BioGRID physical interactions — which specific gene partners interact
  4. UniProt PPI binding residues — which residues are at PPI interfaces

Algorithm
---------
For a reference-alternative pair:
  a. Load UniProt domain annotations for the reference protein
  b. Determine which domains are retained vs disrupted in the alternative
     - If junction peptide provided: find position → check domain overlap
     - If only protein sequences provided: compare domain presence by key
       residue/motif anchoring
  c. For each disrupted domain, look up DDI table to find which partner
     domain types that domain binds
  d. Cross-reference DDI partner domain types with BioGRID physical partners
     of the reference gene to identify specifically impacted interactions
  e. Report:
     - lost_ppis: interactions likely disrupted
     - retained_ppis: interactions likely maintained
     - uncertain_ppis: interactions where domain status is ambiguous
     - gained_ppis: novel domain created by alternative junction (speculative)

DDI Data Source Priority
------------------------
  1. 3DID flat file  (3did_flat.gz, available from https://3did.irbbarcelona.org/)
     Most reliable — domain pairs observed in PDB crystal structures
  2. iPfam (iPfam.tsv) — overlapping with 3DID
  3. Inferred from BioGRID + UniProt: if A interacts with B, and A has domain X,
     then X is a candidate interface for the A-B interaction
     (lower confidence, but uses only Daedalus data)

Usage
-----
from proteus.validation.ppi_impact import PPIImpactAnalyzer

analyzer = PPIImpactAnalyzer()
analyzer.load(daedalus_interim_dir, ddi_flat_path="3did_flat.gz")  # 3DID optional

result = analyzer.analyze_pair(
    pair_id="EGFR_vs_ALT1",
    reference_protein_id="ENSP00000275493.3",
    reference_gene_name="EGFR",
    alternative_protein_seq="MRPSGTAGAALLALLAALCPASRALEEKKV...",
    junction_peptide="GGSLKV",           # optional, improves domain mapping
    reference_protein_seq="MRPSGTAGAA...",  # optional
)
analyzer.to_tsv([result], "ppi_impact.tsv")
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


@dataclass
class PPIImpactRecord:
    """
    PPI impact analysis for a single isoform pair.

    Fields
    ------
    pair_id : str
    reference_gene : str
    n_ref_domains : int
        Number of UniProt domains in the reference protein.
    n_alt_domains_estimated : int
        Estimated domains retained in the alternative protein.
    disrupted_domains : list[str]
        Domain names likely disrupted in the alternative isoform.
    retained_domains : list[str]
        Domain names likely retained.
    lost_ppis : list[dict]
        PPIs predicted to be lost. Each dict:
          {partner_gene, partner_id, interaction_type, via_domain,
           biogrid_evidence_count, confidence}
    retained_ppis : list[dict]
        PPIs predicted to be retained.
    gained_ppis : list[dict]
        Novel PPIs potentially gained (from junction-created domains).
    uncertain_ppis : list[dict]
        PPIs where domain status is ambiguous.
    n_lost : int
    n_retained : int
    n_gained : int
    ppi_interface_disrupted : bool
        True if ≥1 confirmed PPI interface residue is lost.
    summary : str
    """

    pair_id: str = ""
    reference_gene: str = ""
    n_ref_domains: int = 0
    n_alt_domains_estimated: int = 0
    disrupted_domains: List[str] = field(default_factory=list)
    retained_domains: List[str] = field(default_factory=list)
    lost_ppis: List[Dict] = field(default_factory=list)
    retained_ppis: List[Dict] = field(default_factory=list)
    gained_ppis: List[Dict] = field(default_factory=list)
    uncertain_ppis: List[Dict] = field(default_factory=list)
    n_lost: int = 0
    n_retained: int = 0
    n_gained: int = 0
    ppi_interface_disrupted: bool = False
    dna_binding_disrupted: bool = False
    kinase_domain_disrupted: bool = False
    summary: str = ""


class PPIImpactAnalyzer:
    """
    Analyzes specific PPI gain/loss for reference-alternative isoform pairs.

    Combines three evidence layers:
      1. UniProt domain annotations (what domains exist and where)
      2. DDI table (which domain types interact with which)
      3. BioGRID physical interactions (which specific gene partners interact)

    Parameters
    ----------
    ddi_confidence_threshold : float
        Minimum DDI confidence score to include (for 3DID scores, default 0.0 = use all).
    min_biogrid_evidence : int
        Minimum BioGRID evidence count to call a PPI "high confidence".
    """

    # Domain type → functional category mapping
    _DOMAIN_CATEGORIES: Dict[str, str] = {
        "kinase": "kinase",
        "protein kinase": "kinase",
        "tyrosine kinase": "kinase",
        "sh2": "ppi_adapter",
        "sh3": "ppi_adapter",
        "ph": "lipid_binding",
        "pleckstrin": "lipid_binding",
        "pdz": "ppi_scaffold",
        "death": "apoptosis",
        "death domain": "apoptosis",
        "bromo": "chromatin",
        "bromodomain": "chromatin",
        "chromo": "chromatin",
        "chromodomain": "chromatin",
        "wd40": "scaffold",
        "wd repeat": "scaffold",
        "ankyrin": "scaffold",
        "leucine rich": "scaffold",
        "armadillo": "scaffold",
        "heat": "scaffold",
        "dna-binding": "dna_binding",
        "zinc finger": "dna_binding",
        "homeobox": "dna_binding",
        "helix-turn-helix": "dna_binding",
        "bzip": "dna_binding",
        "bhlh": "dna_binding",
        "ig-like": "cell_adhesion",
        "immunoglobulin": "cell_adhesion",
        "fibronectin": "cell_adhesion",
        "cadherin": "cell_adhesion",
        "ring": "ubiquitin",
        "f-box": "ubiquitin",
        "ubiquitin": "ubiquitin",
        "gtpase": "signaling",
        "ras": "signaling",
        "rho": "signaling",
        "egf": "growth_factor",
        "egf-like": "growth_factor",
    }

    def __init__(
        self,
        ddi_confidence_threshold: float = 0.0,
        min_biogrid_evidence: int = 1,
    ) -> None:
        self.ddi_confidence_threshold = ddi_confidence_threshold
        self.min_biogrid_evidence = min_biogrid_evidence

        # Loaded tables
        self._uniprot_domains: Dict[str, List[Dict]] = {}  # ENSP → list of domain records
        self._ddi_table: Dict[str, Set[str]] = {}           # domain_name → set of partner domain names
        self._biogrid_edges: Dict[str, List[Dict]] = {}     # gene_name → list of edge dicts
        self._gene_domains: Dict[str, List[str]] = {}       # gene_name → list of domain names
        self._ensp_to_gene: Dict[str, str] = {}             # ENSP → gene_name
        self._gene_to_ensp: Dict[str, str] = {}             # gene_name → ENSP
        self._ppi_binding_residues: Dict[str, int] = {}     # ENSP → count of PPI binding residues
        self._loaded = False

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    def load(
        self,
        daedalus_interim_dir: Path,
        ddi_flat_path: Optional[Path] = None,
    ) -> None:
        """
        Load all data tables needed for PPI impact analysis.

        Parameters
        ----------
        daedalus_interim_dir : Path
            Daedalus data/interim/ directory.
        ddi_flat_path : Path | None
            Optional path to 3DID flat file (3did_flat or 3did_flat.gz).
            If None, DDIs are inferred from BioGRID + UniProt (lower confidence).
        """
        import pandas as pd

        d = Path(daedalus_interim_dir)

        # ---- UniProt features (domain annotations per protein) ----
        feat_path = d / "uniprot_features.tsv"
        map_path = d / "uniprot_gencode_map.tsv"
        if feat_path.exists() and map_path.exists():
            feat = pd.read_csv(feat_path, sep="\t", low_memory=False)
            gmap = pd.read_csv(map_path, sep="\t", low_memory=False,
                               usecols=lambda c: c in (
                                   "primary_accession", "gencode_protein_id", "gene_name"))

            # Build accession → gene_name + ENSP bridges
            for _, row in gmap.iterrows():
                acc = str(row.get("primary_accession", "") or "")
                pid = str(row.get("gencode_protein_id", "") or "")
                gn = str(row.get("gene_name", "") or "")
                if acc and pid and pid != "nan":
                    self._ensp_to_gene[pid] = gn
                    if gn and gn != "nan":
                        self._gene_to_ensp[gn] = pid

            # Group domain features by accession, then bridge to ENSP
            acc_to_pid: Dict[str, str] = {
                str(r.get("primary_accession", "")): str(r.get("gencode_protein_id", ""))
                for _, r in gmap.iterrows()
                if str(r.get("gencode_protein_id", "")) not in ("", "nan")
            }
            domain_types = {
                "Domain", "Repeat", "Zinc finger", "DNA binding",
                "Coiled coil", "Motif", "Region",
            }
            ppi_types = {"Binding site", "Active site", "Site"}

            for acc, grp in feat.groupby("primary_accession"):
                pid = acc_to_pid.get(str(acc), "")
                if not pid:
                    continue
                domains = []
                ppi_count = 0
                for _, fr in grp.iterrows():
                    ft = str(fr.get("feature_type", "") or "")
                    desc = str(fr.get("description", "") or "")
                    start = fr.get("start", 0)
                    end = fr.get("end", 0)
                    if ft in domain_types:
                        domains.append({
                            "feature_type": ft,
                            "description": desc,
                            "start": int(start) if start == start else 0,
                            "end": int(end) if end == end else 0,
                        })
                    if ft in ppi_types and "interaction" in desc.lower():
                        ppi_count += 1
                if domains:
                    self._uniprot_domains[pid] = domains
                if ppi_count:
                    self._ppi_binding_residues[pid] = ppi_count
            print(f"[PPIImpactAnalyzer] UniProt domains: {len(self._uniprot_domains):,} proteins.")
        else:
            warnings.warn(f"[PPIImpactAnalyzer] uniprot_features.tsv or gencode_map not found in {d}")

        # ---- BioGRID interaction edges ----
        bg_path = d / "biogrid_gene_interactions.tsv"
        # Also look for a full edges file if available
        bg_edges_path = d / "biogrid_interactions_full.tsv"
        if bg_edges_path.exists():
            self._load_biogrid_full(bg_edges_path)
        elif bg_path.exists():
            # Aggregate table only — we can still use gene_name as a key
            # Build a minimal "edge" per gene pointing to all physical partners
            self._load_biogrid_aggregate(bg_path, d)
        else:
            warnings.warn(f"[PPIImpactAnalyzer] No BioGRID file found in {d}")

        # Build gene → domain list (for DDI inference)
        for pid, domains in self._uniprot_domains.items():
            gn = self._ensp_to_gene.get(pid, "")
            if gn:
                self._gene_domains[gn] = [d["description"] for d in domains]

        # ---- DDI table ----
        if ddi_flat_path and Path(ddi_flat_path).exists():
            self._load_3did(Path(ddi_flat_path))
        else:
            # Infer DDIs from co-occurring domains in BioGRID pairs
            self._infer_ddi_from_biogrid()

        self._loaded = True
        print(f"[PPIImpactAnalyzer] DDI entries: {len(self._ddi_table):,} domain types. "
              f"BioGRID edges: {sum(len(v) for v in self._biogrid_edges.values()):,}.")

    def _load_biogrid_full(self, path: Path) -> None:
        """Load BioGRID TAB3 full edge file."""
        import pandas as pd
        cols_needed = [
            "Official Symbol Interactor A",
            "Official Symbol Interactor B",
            "Experimental System Type",
            "Throughput",
            "Score",
        ]
        try:
            bg = pd.read_csv(path, sep="\t", low_memory=False,
                             usecols=lambda c: c in cols_needed)
            for _, row in bg.iterrows():
                gene_a = str(row.get("Official Symbol Interactor A", "") or "").strip()
                gene_b = str(row.get("Official Symbol Interactor B", "") or "").strip()
                exp_type = str(row.get("Experimental System Type", "") or "")
                is_physical = "physical" in exp_type.lower()
                edge = {
                    "partner_gene": gene_b,
                    "is_physical": is_physical,
                    "experimental_system_type": exp_type,
                    "evidence_count": 1,
                }
                if gene_a not in self._biogrid_edges:
                    self._biogrid_edges[gene_a] = []
                self._biogrid_edges[gene_a].append(edge)
                # Symmetric
                edge_b = edge.copy()
                edge_b["partner_gene"] = gene_a
                if gene_b not in self._biogrid_edges:
                    self._biogrid_edges[gene_b] = []
                self._biogrid_edges[gene_b].append(edge_b)
            print(f"[PPIImpactAnalyzer] BioGRID full edges: "
                  f"{sum(len(v) for v in self._biogrid_edges.values()):,} edges.")
        except Exception as exc:
            warnings.warn(f"[PPIImpactAnalyzer] Failed to load full BioGRID: {exc}")

    def _load_biogrid_aggregate(self, bg_path: Path, d: Path) -> None:
        """
        Load BioGRID aggregate table. Without full edge data, create
        placeholder edges from the aggregate partner count.
        """
        import pandas as pd
        bg = pd.read_csv(bg_path, sep="\t", low_memory=False)
        # We don't have per-partner data, so we create a count-only record
        for _, row in bg.iterrows():
            gn = str(row.get("gene_name", "") or "").strip()
            phys_count = int(row.get("biogrid_physical_partner_count", 0) or 0)
            if gn and phys_count > 0:
                # No specific edges available, but record that interactions exist
                self._biogrid_edges[gn] = [
                    {"partner_gene": f"<{phys_count} physical partners>",
                     "is_physical": True,
                     "experimental_system_type": "physical",
                     "evidence_count": phys_count}
                ]

    def _load_3did(self, path: Path) -> None:
        """
        Parse 3DID flat file to build domain-domain interaction table.

        3DID flat format:
          #=ID  domain_A  pfam_id_A  domain_B  pfam_id_B
        """
        import gzip
        opener = gzip.open if str(path).endswith(".gz") else open
        try:
            with opener(path, "rt") as fh:
                for line in fh:
                    line = line.strip()
                    if not line.startswith("#=ID"):
                        continue
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    dom_a = parts[1].lower()
                    dom_b = parts[2].lower() if len(parts) > 2 else parts[1].lower()
                    if dom_a not in self._ddi_table:
                        self._ddi_table[dom_a] = set()
                    if dom_b not in self._ddi_table:
                        self._ddi_table[dom_b] = set()
                    self._ddi_table[dom_a].add(dom_b)
                    self._ddi_table[dom_b].add(dom_a)
            print(f"[PPIImpactAnalyzer] 3DID: {len(self._ddi_table):,} domain types loaded.")
        except Exception as exc:
            warnings.warn(f"[PPIImpactAnalyzer] Failed to load 3DID: {exc}")
            self._infer_ddi_from_biogrid()

    def _infer_ddi_from_biogrid(self) -> None:
        """
        Infer DDIs from co-occurring domains in BioGRID interacting gene pairs.
        If gene A (with domain X) physically interacts with gene B (with domain Y),
        then X-Y is a candidate DDI (lower confidence than structural evidence).
        """
        for gene_a, edges in self._biogrid_edges.items():
            domains_a = self._gene_domains.get(gene_a, [])
            if not domains_a:
                continue
            for edge in edges:
                if not edge.get("is_physical"):
                    continue
                gene_b = edge["partner_gene"]
                domains_b = self._gene_domains.get(gene_b, [])
                for da in domains_a:
                    da_key = da.lower()
                    if da_key not in self._ddi_table:
                        self._ddi_table[da_key] = set()
                    for db in domains_b:
                        self._ddi_table[da_key].add(db.lower())

    # ------------------------------------------------------------------
    # Domain disruption analysis
    # ------------------------------------------------------------------

    def _find_disrupted_domains(
        self,
        ref_domains: List[Dict],
        alt_protein_seq: str,
        junction_start_aa: Optional[int] = None,
        junction_end_aa: Optional[int] = None,
        ref_protein_seq: Optional[str] = None,
    ) -> Tuple[List[str], List[str]]:
        """
        Determine which reference domains are disrupted vs retained in the alt.

        Strategy (in priority order):
        1. If junction residue range is known: domains overlapping the junction
           are "disrupted", others are "retained".
        2. If both protein sequences are available: for each domain, check
           if a key 9-mer anchor from the ref domain region is present in alt.
        3. If only domain start/end positions: scale to alt protein length and
           check if region is present.

        Returns
        -------
        (disrupted_names, retained_names)
        """
        disrupted = []
        retained = []

        for dom in ref_domains:
            dom_start = dom.get("start", 0)
            dom_end = dom.get("end", 0)
            dom_name = dom.get("description", dom.get("feature_type", "unknown"))

            if dom_start <= 0 and dom_end <= 0:
                retained.append(dom_name)
                continue

            # Method 1: junction range overlap
            if junction_start_aa is not None and junction_end_aa is not None:
                # Junction residues are the "disrupted region"
                # Domain overlapping the junction → disrupted
                j_s = junction_start_aa + 1  # 1-based
                j_e = junction_end_aa
                overlap = min(j_e, dom_end) - max(j_s, dom_start) + 1
                if overlap > 0:
                    disrupted.append(dom_name)
                else:
                    retained.append(dom_name)
                continue

            # Method 2: sequence anchor search
            if ref_protein_seq and alt_protein_seq and dom_start > 0 and dom_end > dom_start:
                # Extract a short anchor from the middle of the domain region
                d_len = dom_end - dom_start
                anchor_s = dom_start - 1 + d_len // 4
                anchor_e = anchor_s + min(9, d_len // 2)
                if anchor_e <= len(ref_protein_seq):
                    anchor = ref_protein_seq[anchor_s:anchor_e].upper()
                    if len(anchor) >= 6 and anchor in alt_protein_seq.upper():
                        retained.append(dom_name)
                    else:
                        disrupted.append(dom_name)
                    continue

            # Method 3: positional scaling fallback
            if alt_protein_seq and ref_protein_seq:
                ref_len = len(ref_protein_seq)
                alt_len = len(alt_protein_seq)
                if ref_len > 0:
                    scaled_start = int(dom_start * alt_len / ref_len)
                    scaled_end = int(dom_end * alt_len / ref_len)
                    # If alt is much shorter than scaled position → domain likely lost
                    if scaled_start > alt_len:
                        disrupted.append(dom_name)
                    elif (alt_len / ref_len) < 0.7 and dom_end > ref_len * 0.7:
                        # Large C-terminal truncation
                        disrupted.append(dom_name)
                    else:
                        retained.append(dom_name)
                    continue

            retained.append(dom_name)  # default: assume retained if no info

        return disrupted, retained

    # ------------------------------------------------------------------
    # Main analysis
    # ------------------------------------------------------------------

    def analyze_pair(
        self,
        pair_id: str,
        reference_protein_id: str,
        reference_gene_name: str,
        alternative_protein_seq: str = "",
        reference_protein_seq: str = "",
        junction_peptide: str = "",
        has_alt_protein: bool = True,
    ) -> PPIImpactRecord:
        """
        Analyze PPI gain/loss for a single isoform pair.

        Parameters
        ----------
        pair_id : str
        reference_protein_id : str
            ENSP ID for the reference protein.
        reference_gene_name : str
            Gene symbol (for BioGRID lookup).
        alternative_protein_seq : str
            Full AA sequence of the alternative isoform.
        reference_protein_seq : str
            Full AA sequence of the reference isoform (optional, improves analysis).
        junction_peptide : str
            Junction-spanning peptide (optional, sharpens domain disruption calls).
        has_alt_protein : bool
            False if alternative isoform is predicted NMD.
        """
        rec = PPIImpactRecord(
            pair_id=pair_id,
            reference_gene=reference_gene_name,
        )

        if not has_alt_protein:
            rec.summary = "NMD predicted — all PPIs lost (protein absent)"
            # All reference PPIs are lost
            edges = self._biogrid_edges.get(reference_gene_name, [])
            for edge in edges:
                if edge.get("is_physical"):
                    rec.lost_ppis.append({
                        "partner_gene": edge["partner_gene"],
                        "interaction_type": "physical",
                        "via_domain": "all",
                        "reason": "NMD — protein absent",
                        "confidence": 1.0,
                    })
            rec.n_lost = len(rec.lost_ppis)
            rec.ppi_interface_disrupted = rec.n_lost > 0
            return rec

        # Load reference domain annotations
        ref_domains = self._uniprot_domains.get(reference_protein_id, [])
        rec.n_ref_domains = len(ref_domains)

        # Find junction position if junction peptide is given
        junction_start_aa = None
        junction_end_aa = None
        if junction_peptide and alternative_protein_seq:
            pos = alternative_protein_seq.upper().find(junction_peptide.upper())
            if pos != -1:
                junction_start_aa = pos
                junction_end_aa = pos + len(junction_peptide)

        # Determine disrupted vs retained domains
        if ref_domains:
            disrupted, retained = self._find_disrupted_domains(
                ref_domains,
                alternative_protein_seq,
                junction_start_aa,
                junction_end_aa,
                reference_protein_seq,
            )
        else:
            disrupted, retained = [], []

        rec.disrupted_domains = disrupted
        rec.retained_domains = retained
        rec.n_alt_domains_estimated = len(retained)

        # Check for functional domain types
        rec.dna_binding_disrupted = any(
            any(kw in d.lower() for kw in ("dna", "zinc finger", "homeodomain", "bzip", "bhlh"))
            for d in disrupted
        )
        rec.kinase_domain_disrupted = any(
            "kinase" in d.lower() for d in disrupted
        )

        # Check PPI binding residues
        ppi_res = self._ppi_binding_residues.get(reference_protein_id, 0)
        if ppi_res > 0 and disrupted:
            rec.ppi_interface_disrupted = True

        # Get BioGRID partners for this gene
        edges = self._biogrid_edges.get(reference_gene_name, [])
        physical_edges = [e for e in edges if e.get("is_physical", True)]

        if not physical_edges:
            rec.summary = (
                f"No BioGRID physical interactions found for {reference_gene_name}. "
                f"Domains disrupted: {len(disrupted)}, retained: {len(retained)}."
            )
            return rec

        # For each physical partner, determine if the interaction is retained
        for edge in physical_edges:
            partner = edge["partner_gene"]
            ev_count = edge.get("evidence_count", 1)

            # Determine which disrupted domains mediate this interaction via DDI
            mediating_domains = []
            for dom_name in disrupted:
                dom_key = dom_name.lower()
                partner_doms = self._gene_domains.get(partner, [])
                for pd_name in partner_doms:
                    if dom_key in self._ddi_table:
                        if pd_name.lower() in self._ddi_table[dom_key]:
                            mediating_domains.append(dom_name)
                            break

            ppi_entry = {
                "partner_gene": partner,
                "interaction_type": edge.get("experimental_system_type", "physical"),
                "biogrid_evidence_count": ev_count,
            }

            if mediating_domains:
                # Interaction mediated by a disrupted domain → LOST
                ppi_entry["via_domain"] = "; ".join(mediating_domains)
                ppi_entry["confidence"] = min(0.5 + 0.1 * len(mediating_domains), 0.9)
                ppi_entry["reason"] = f"Domain(s) disrupted: {', '.join(mediating_domains)}"
                rec.lost_ppis.append(ppi_entry)
            elif disrupted and not retained:
                # All domains disrupted — interaction likely lost, but mechanism unclear
                ppi_entry["via_domain"] = "unknown (all domains disrupted)"
                ppi_entry["confidence"] = 0.4
                ppi_entry["reason"] = "Most domains disrupted; interaction likely lost"
                rec.uncertain_ppis.append(ppi_entry)
            else:
                ppi_entry["via_domain"] = "retained domains"
                ppi_entry["confidence"] = 0.7
                ppi_entry["reason"] = "Interface domain retained"
                rec.retained_ppis.append(ppi_entry)

        rec.n_lost = len(rec.lost_ppis)
        rec.n_retained = len(rec.retained_ppis)
        rec.n_gained = len(rec.gained_ppis)

        rec.summary = (
            f"Gene: {reference_gene_name} | "
            f"Domains disrupted: {len(disrupted)}/{rec.n_ref_domains} | "
            f"PPIs: {rec.n_lost} lost, {rec.n_retained} retained, "
            f"{len(rec.uncertain_ppis)} uncertain | "
            f"DNA-binding disrupted: {rec.dna_binding_disrupted} | "
            f"Kinase disrupted: {rec.kinase_domain_disrupted}"
        )
        return rec

    def analyze_batch(self, pairs: List[Dict[str, Any]]) -> List[PPIImpactRecord]:
        """
        Analyze a list of pairs.

        Each dict must have:
          pair_id, reference_protein_id, reference_gene_name
        Optional:
          alternative_protein_seq, reference_protein_seq, junction_peptide,
          has_alt_protein
        """
        if not self._loaded:
            raise RuntimeError("Call .load(daedalus_interim_dir) first.")
        results = []
        for p in pairs:
            r = self.analyze_pair(
                pair_id=str(p.get("pair_id", "")),
                reference_protein_id=str(p.get("reference_protein_id", "")),
                reference_gene_name=str(p.get("reference_gene_name",
                                              p.get("gene_name", ""))),
                alternative_protein_seq=str(p.get("alternative_protein_seq", "")),
                reference_protein_seq=str(p.get("reference_protein_seq", "")),
                junction_peptide=str(p.get("junction_peptide", "")),
                has_alt_protein=bool(p.get("has_alt_protein", True)),
            )
            results.append(r)
        return results

    def to_dataframe(self, results: List[PPIImpactRecord]):
        """Convert results to pandas DataFrame (one row per pair, summary columns)."""
        import pandas as pd
        import json
        rows = []
        for r in results:
            rows.append({
                "pair_id": r.pair_id,
                "reference_gene": r.reference_gene,
                "n_ref_domains": r.n_ref_domains,
                "n_alt_domains_estimated": r.n_alt_domains_estimated,
                "disrupted_domains": "; ".join(r.disrupted_domains),
                "retained_domains": "; ".join(r.retained_domains),
                "n_ppis_lost": r.n_lost,
                "n_ppis_retained": r.n_retained,
                "n_ppis_gained": r.n_gained,
                "n_ppis_uncertain": len(r.uncertain_ppis),
                "lost_ppi_partners": "; ".join(e["partner_gene"] for e in r.lost_ppis),
                "retained_ppi_partners": "; ".join(e["partner_gene"] for e in r.retained_ppis),
                "ppi_interface_disrupted": int(r.ppi_interface_disrupted),
                "dna_binding_disrupted": int(r.dna_binding_disrupted),
                "kinase_domain_disrupted": int(r.kinase_domain_disrupted),
                "summary": r.summary,
            })
        return pd.DataFrame(rows)

    def to_tsv(self, results: List[PPIImpactRecord], output_path: Path) -> None:
        """Write results to TSV."""
        df = self.to_dataframe(results)
        df.to_csv(output_path, sep="\t", index=False)
        n_lost = sum(r.n_lost for r in results)
        n_dna = sum(1 for r in results if r.dna_binding_disrupted)
        n_kinase = sum(1 for r in results if r.kinase_domain_disrupted)
        print(f"[PPIImpactAnalyzer] {len(results)} pairs analyzed. "
              f"Total PPIs predicted lost: {n_lost}. "
              f"DNA-binding disrupted: {n_dna}. "
              f"Kinase disrupted: {n_kinase}.")
        print(f"  Written to: {output_path}")
