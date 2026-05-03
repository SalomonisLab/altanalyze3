"""
SyntheticExonSkipGenerator: Generates synthetic exon-skip pairs for SSL pretraining.

Reads GENCODE GTF and UniProt annotation tables to enumerate all possible
single-exon-skip variants and compute rule-based functional change labels.
"""

from __future__ import annotations

import gzip
import re
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def _parse_attribute(attr_str: str, key: str) -> Optional[str]:
    """Parse a GTF attribute string and return the value for a given key."""
    pattern = rf'{key}\s+"([^"]+)"'
    m = re.search(pattern, attr_str)
    return m.group(1) if m else None


def _parse_gtf(gtf_path: Path, feature_type: str = "exon") -> pd.DataFrame:
    """
    Parse a GENCODE GTF file (plain or gzipped) and return a DataFrame
    of the specified feature type.

    Parameters
    ----------
    gtf_path : Path
        Path to GTF file (.gtf or .gtf.gz).
    feature_type : str
        Feature type to extract (default: "exon").

    Returns
    -------
    pd.DataFrame
        Columns: chrom, start, end, strand, transcript_id, gene_id,
        transcript_type, exon_number, transcript_name
    """
    rows = []
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(str(gtf_path), "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != feature_type:
                continue
            attrs = parts[8]
            rows.append({
                "chrom": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "transcript_id": _parse_attribute(attrs, "transcript_id") or "",
                "gene_id": _parse_attribute(attrs, "gene_id") or "",
                "transcript_type": _parse_attribute(attrs, "transcript_type") or "",
                "exon_number": _parse_attribute(attrs, "exon_number") or "",
                "transcript_name": _parse_attribute(attrs, "transcript_name") or "",
            })
    return pd.DataFrame(rows)


def _parse_cds(gtf_path: Path) -> pd.DataFrame:
    """Parse CDS features from GTF."""
    return _parse_gtf(gtf_path, feature_type="CDS")


def _compute_nmd(
    exons: List[Tuple[int, int]],
    skip_idx: int,
    cds_start_genomic: int,
    strand: str,
) -> bool:
    """
    Determine if an exon skip triggers NMD via the canonical 50/55-nt rule.

    NMD is triggered if the exon skip creates a PTC > 55 nt upstream of the
    last exon-exon junction (EJC). We approximate this using frame tracking
    on the remaining exon sequence lengths.

    Parameters
    ----------
    exons : list[tuple[int, int]]
        List of (start, end) genomic coordinates (0-based, half-open) in
        transcript order (5' to 3', already sorted for strand).
    skip_idx : int
        Index of the skipped exon in the sorted exon list.
    cds_start_genomic : int
        Genomic coordinate of the CDS start (approximate).
    strand : str
        "+" or "-".

    Returns
    -------
    bool
        True if skip likely triggers NMD.
    """
    # Build list of exons without the skipped one
    remaining = [e for i, e in enumerate(exons) if i != skip_idx]

    if len(remaining) < 2:
        return False

    # Compute transcript lengths of remaining exons
    exon_lengths = [e[1] - e[0] for e in remaining]
    total_len = sum(exon_lengths)

    # The skipped exon length
    skipped_len = exons[skip_idx][1] - exons[skip_idx][0]

    # If skip causes frameshift (skipped length not divisible by 3),
    # a PTC will be introduced early in the remaining sequence.
    frame_shift = skipped_len % 3

    if frame_shift == 0:
        # In-frame skip — no frameshift, no PTC from this skip alone
        return False

    # With a frameshift, the PTC is created at the first stop codon in the
    # new frame. Conservatively estimate PTC is introduced early.
    # NMD rule: PTC > 55 nt upstream of last EJC.
    # Last EJC is 20-24 nt upstream of the last exon-exon junction.
    # If there are >= 2 remaining exons and the last exon is < total_len - 55,
    # then NMD is likely.

    if len(remaining) < 2:
        return False

    last_exon_len = exon_lengths[-1]
    dist_from_last_ejc = total_len - sum(exon_lengths[:-1])  # = last exon length

    # Frameshift + PTC expected early → NMD if last EJC is > 55 nt downstream of PTC
    # Conservative: assume PTC is in the first third of post-skip transcript
    # and check if last_exon_len >= 55
    return last_exon_len >= 55


def _process_transcript(
    tid: str,
    gene_id: str,
    strand: str,
    exon_coords: List[Tuple[int, int]],
    protein_id: str,
    feat_row: Optional[Any],
) -> List[Dict[str, Any]]:
    """
    Generate all single-exon-skip variants for one transcript.

    Parameters
    ----------
    tid : str
        Transcript ID.
    gene_id : str
        Gene ID.
    strand : str
        "+"/"-".
    exon_coords : list[tuple]
        Sorted exon (start, end) coordinates in transcript order.
    protein_id : str
        Corresponding UniProt/Ensembl protein ID.
    feat_row : row | None
        UniProt feature row for this protein.

    Returns
    -------
    list[dict]
        One dict per possible single-exon-skip.
    """
    results = []

    if len(exon_coords) < 3:
        return results  # Need >= 3 exons to skip an internal one

    # Find CDS start (approximate: first exon)
    cds_start_approx = exon_coords[0][0]

    # Get protein features
    n_tm_helices = 0.0
    has_domain = False
    has_active_site = False

    if feat_row is not None:
        try:
            n_tm_helices = float(feat_row.get("n_tm_helices", 0) or 0)
            n_domains = float(feat_row.get("n_domains", 0) or 0)
            n_active = float(feat_row.get("n_active_sites", 0) or 0)
            has_domain = n_domains > 0
            has_active_site = n_active > 0
        except (TypeError, ValueError):
            pass

    total_cds_len = sum(e[1] - e[0] for e in exon_coords)

    # Skip internal exons only (not first or last)
    for skip_idx in range(1, len(exon_coords) - 1):
        skipped_start, skipped_end = exon_coords[skip_idx]
        skipped_len = skipped_end - skipped_start

        # 1. domain_disrupted: heuristic using exon position + domain presence
        #    Use: if protein has domain(s) AND skipped exon is in the "functional middle"
        #    (exons 20-80% of total CDS), flag as domain-disrupted.
        cumulative_before = sum(
            e[1] - e[0] for e in exon_coords[:skip_idx]
        )
        rel_pos = cumulative_before / max(total_cds_len, 1)
        domain_disrupted = (has_domain or has_active_site) and (0.1 < rel_pos < 0.9)

        # 2. nmd_triggered: frame tracking
        nmd_triggered = _compute_nmd(
            exon_coords, skip_idx, cds_start_approx, strand
        )

        # 3. tm_lost: heuristic — protein has TM helices AND skipped exon is in middle third
        tm_lost = (n_tm_helices > 0) and (0.33 < rel_pos < 0.67)

        # 4. disorder_delta: 0.0 (populated if disorder cache is available)
        disorder_delta = 0.0

        # Frame preservation
        frame_preserved = (skipped_len % 3 == 0)

        alt_tid = f"{tid}__skip_exon{skip_idx}"

        results.append({
            "reference_transcript_id": tid,
            "alternative_transcript_id": alt_tid,
            "reference_protein_id": protein_id,
            "gene_id": gene_id,
            "skip_exon_index": skip_idx,
            "skip_exon_start_nt": skipped_start,
            "skip_exon_end_nt": skipped_end,
            "skipped_exon_len": skipped_len,
            "domain_disrupted": int(domain_disrupted),
            "nmd_triggered": int(nmd_triggered),
            "tm_lost": int(tm_lost),
            "disorder_delta": disorder_delta,
            "frame_preserved": int(frame_preserved),
            "strand": strand,
        })

    return results


def _process_gene_chunk(args: Tuple) -> List[Dict[str, Any]]:
    """Worker function for parallel processing of a chunk of transcripts."""
    transcript_chunk, feat_lookup = args
    all_results = []
    for tid, gene_id, strand, exon_coords, protein_id in transcript_chunk:
        feat_row = feat_lookup.get(protein_id) or feat_lookup.get(tid)
        results = _process_transcript(tid, gene_id, strand, exon_coords, protein_id, feat_row)
        all_results.extend(results)
    return all_results


class SyntheticExonSkipGenerator:
    """
    Generates synthetic exon-skip pairs for SSL pretraining.

    Reads GENCODE GTF and Daedalus Phase A UniProt tables to enumerate
    all possible single-exon-skip variants from multi-exon protein-coding
    transcripts. Computes rule-based functional change labels.

    Usage
    -----
    gen = SyntheticExonSkipGenerator()
    gen.generate(
        gtf_path=Path("data/raw/gencode/gencode.v44.annotation.gtf.gz"),
        daedalus_interim=Path("/path/to/daedalus/phase_a/data/interim"),
        output_path=Path("data/interim/pretrain_events/synthetic_exon_skips.tsv"),
    )
    """

    def generate(
        self,
        gtf_path: Path,
        daedalus_interim: Path,
        output_path: Path,
        max_genes: Optional[int] = None,
        n_workers: int = 4,
    ) -> pd.DataFrame:
        """
        Main generation method.

        Parameters
        ----------
        gtf_path : Path
            Path to GENCODE GTF file (.gtf or .gtf.gz).
        daedalus_interim : Path
            Path to Daedalus Phase A interim directory.
        output_path : Path
            Output TSV path for synthetic exon-skip events.
        max_genes : int | None
            Limit to first N genes (for debugging).
        n_workers : int
            Number of parallel worker processes.

        Returns
        -------
        pd.DataFrame
            Generated events dataframe.
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        print(f"[SyntheticExonSkipGenerator] Parsing GTF: {gtf_path}")
        exon_df = _parse_gtf(gtf_path, feature_type="exon")
        exon_df = exon_df[exon_df["transcript_type"] == "protein_coding"].copy()

        print(f"  Parsed {len(exon_df)} exon records for protein-coding transcripts.")

        # Load UniProt features if available
        # uniprot_features.tsv is long-format (one row per feature per protein).
        # uniprot_gencode_map.tsv maps primary_accession -> gencode_protein_id/transcript_id.
        feat_lookup: Dict[str, Any] = {}
        feat_path = Path(daedalus_interim) / "uniprot_features.tsv"
        map_path  = Path(daedalus_interim) / "uniprot_gencode_map.tsv"
        if feat_path.exists():
            feat_df = pd.read_csv(feat_path, sep="\t", low_memory=False)
            # Aggregate per accession: count TM helices, domains, active sites
            agg: Dict[str, Dict[str, float]] = {}
            for _, row in feat_df.iterrows():
                acc = str(row.get("primary_accession", ""))
                if not acc:
                    continue
                if acc not in agg:
                    agg[acc] = {"n_tm_helices": 0.0, "n_domains": 0.0, "n_active_sites": 0.0}
                ftype = str(row.get("feature_type", "")).lower()
                if "transmembrane" in ftype:
                    agg[acc]["n_tm_helices"] += 1
                elif ftype == "domain":
                    agg[acc]["n_domains"] += 1
                elif "active site" in ftype:
                    agg[acc]["n_active_sites"] += 1

            # Build lookup keyed by ENSP/ENST ID using the gencode map
            if map_path.exists():
                map_df = pd.read_csv(map_path, sep="\t", low_memory=False)
                for _, mrow in map_df.iterrows():
                    acc = str(mrow.get("primary_accession", ""))
                    if acc not in agg:
                        continue
                    for id_col in ("gencode_protein_id", "gencode_transcript_id"):
                        gid = str(mrow.get(id_col, ""))
                        if gid and gid != "nan":
                            feat_lookup[gid] = agg[acc]
                print(f"  Loaded {len(feat_lookup)} ENSP/ENST → UniProt feature entries.")
            else:
                # Fall back to accession-keyed lookup
                feat_lookup = agg  # type: ignore
                print(f"  Loaded {len(feat_lookup)} accession → UniProt feature entries (no gencode map).")
        else:
            print(f"  uniprot_features.tsv not found; domain_disrupted/tm_lost will be 0.")

        # Build per-transcript exon lists
        transcript_data: Dict[str, Dict] = {}
        for _, row in exon_df.iterrows():
            tid = row["transcript_id"]
            if tid not in transcript_data:
                transcript_data[tid] = {
                    "gene_id": row["gene_id"],
                    "strand": row["strand"],
                    "exons": [],
                    "protein_id": "",
                }
            transcript_data[tid]["exons"].append(
                (int(row["start"]) - 1, int(row["end"]))  # convert to 0-based half-open
            )

        # Sort exons by genomic position (strand-aware)
        transcripts_list = []
        gene_ids_seen: Dict[str, int] = {}

        for tid, data in transcript_data.items():
            if len(data["exons"]) < 3:
                continue  # need >= 3 exons

            gid = data["gene_id"]
            if max_genes is not None:
                if gid not in gene_ids_seen:
                    gene_ids_seen[gid] = len(gene_ids_seen)
                if gene_ids_seen[gid] >= max_genes:
                    continue

            strand = data["strand"]
            exons = sorted(data["exons"], key=lambda x: x[0])
            if strand == "-":
                exons = list(reversed(exons))

            # Use transcript_id as protein_id proxy (Ensembl ENSP)
            protein_id = data.get("protein_id", "") or tid.replace("ENST", "ENSP")

            transcripts_list.append((tid, gid, strand, exons, protein_id))

        print(
            f"  Processing {len(transcripts_list)} multi-exon transcripts "
            f"with {n_workers} workers..."
        )

        # Chunk for parallel processing
        chunk_size = max(100, len(transcripts_list) // (n_workers * 4))
        chunks = [
            transcripts_list[i : i + chunk_size]
            for i in range(0, len(transcripts_list), chunk_size)
        ]

        all_results = []
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(_process_gene_chunk, (chunk, feat_lookup)): i
                for i, chunk in enumerate(chunks)
            }
            completed = 0
            for future in as_completed(futures):
                try:
                    results = future.result()
                    all_results.extend(results)
                except Exception as exc:
                    warnings.warn(f"Worker failed: {exc}", stacklevel=2)
                completed += 1
                if completed % 10 == 0:
                    print(
                        f"  Progress: {completed}/{len(chunks)} chunks, "
                        f"{len(all_results)} events so far"
                    )

        df = pd.DataFrame(all_results)
        df.to_csv(output_path, sep="\t", index=False)

        print(
            f"[SyntheticExonSkipGenerator] Saved {len(df)} synthetic exon-skip events "
            f"to {output_path}"
        )
        print(f"  domain_disrupted: {df['domain_disrupted'].sum()}")
        print(f"  nmd_triggered:    {df['nmd_triggered'].sum()}")
        print(f"  tm_lost:          {df['tm_lost'].sum()}")

        return df
