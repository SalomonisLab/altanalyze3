"""Consolidate every isoform interaction into one flat TXT, keyed to the source (Vidal lab 2016/2025 /
UniProt) isoform ID harmonized to the GFF "final name", its best-mapped ENST, its Ens91 structure, and
the genomic exon-coordinate tuples of that structure.

Columns (tab-delimited, with header):
  1. Symbol             gene symbol
  2. ENSG               Ensembl gene id
  3. isoform_structure  Ens91 E/I structure string ('' if the isoform has only a coordinate/UNK locus)
  4. exon_coordinates   genomic exon tuples 'chrom:strand:(s1,e1),(s2,e2),...' for the structure
  5. interaction_type   PPI or PDI
  6. target             partner gene symbol (PPI) or DNA bait id (PDI)
  7. source             Vidal2025_Y2H | Vidal2025_eY1H | Vidal2016_Y2H | UniProt
  8. source_isoform_id  Vidal2025 -> GFF final name (SYMBOL|i/n|well); Vidal2016 -> Yang Isoform_ID;
                          UniProt -> UniProt id
  9. best_ENST          best-mapped Ensembl transcript id

Vidal2025 (atlas) structures + IDs come from the GFF transcript_associations (the authoritative 2025
catalogue annotation: every clone has a structure, E/I or coordinate); exon coordinates from the
c_6k GFF exon lines (exact, all clones). Yang/UniProt structures from the high-accuracy project
mappings, with structureless ones recovered from the MDS-AML full-length catalog
(structureless_resolved.tsv); their exon coordinates from structure_coords.sqlite. Detected/positive
interactions only; de-duplicated on all columns; non-ENSG / non-E-I placeholders sanitized.
"""

import os
import re
import logging

import pandas as pd

from . import config
from .crosswalk import reference_structures

logger = logging.getLogger(__name__)

OUTPUT_NAME = "isoform_interactions.txt"
COLUMNS = ["Symbol", "ENSG", "isoform_structure", "exon_coordinates", "interaction_type", "target",
           "source", "source_isoform_id", "best_ENST", "activity_score"]
_ENSG_RE = re.compile(r"^ENSG[0-9]+")
_EI_RE = re.compile(r"[EI][0-9]")
_TID_RE = re.compile(r'transcript_id[ =]"?([^";]+)"?')


def _read(path, what):
    p = config.resolve_shared(path)
    if not os.path.exists(p):
        logger.warning("[iso2function.export] missing %s (%s)", path, what)
        return None
    return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, na_values=[""]).fillna("")


def _is_true(v):
    return str(v).strip().lower() in config.TRUE_TOKENS


def _clean_structure(s):
    """Keep only valid E#/I# tokens; drop coordinate/UNK pseudo-structures and any embedded garbage
    token (e.g. a stray 'ENSG...' that leaked into a session structure string)."""
    s = str(s or "").strip()
    if not s:
        return ""
    toks = [t for t in s.split("|") if re.match(r"^[EI][0-9]", t)]
    return "|".join(toks)


def _clean_ensg(g):
    g = str(g or "").strip()
    return g if _ENSG_RE.match(g) else ""


def _gencode_enst_to_ensg():
    """{base_ENST -> base_ENSG} from BOTH the GENCODE pc_translations FASTA (>ENSP|ENST|ENSG|...) and the
    Ensembl116 pep FASTA ('gene:ENSG... transcript:ENST...'). Together they resolve ENST->ENSG for any
    transcript (GENCODE v49 + the Ensembl116-only ids the session matched against)."""
    import gzip
    out = {}
    pg = config.resolve_shared(config.GENCODE_PROTEIN_FASTA)
    if os.path.exists(pg):
        opener = gzip.open if pg.endswith(".gz") else open
        with opener(pg, "rt") as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                parts = line[1:].split("|")
                enst = next((x.split(".")[0] for x in parts if x.startswith("ENST")), "")
                ensg = next((x.split(".")[0] for x in parts if x.startswith("ENSG")), "")
                if enst and ensg:
                    out.setdefault(enst, ensg)
    pe = config.resolve_shared(config.ENSEMBL_PEP_FASTA)
    if os.path.exists(pe):
        gre = re.compile(r"gene:(ENSG[0-9]+)"); tre = re.compile(r"transcript:(ENST[0-9]+)")
        with gzip.open(pe, "rt") as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                g, t = gre.search(line), tre.search(line)
                if g and t:
                    out.setdefault(t.group(1), g.group(1))
    logger.info("[iso2function.export] ENST->ENSG map: %d (GENCODE + Ensembl116 pep)", len(out))
    return out


def _enst_utr_variants():
    """{base_ENST -> [protein-identical, same-gene base_ENSTs]} — the UTR variants of each transcript.
    Two transcripts of the same gene that encode an identical protein differ only in 5'/3' UTR (or CDS
    start), so one cloned ORF legitimately maps to all of them. Built from BOTH the GENCODE
    pc_translations FASTA and the Ensembl116 pep FASTA, grouped by (ENSG, protein sequence)."""
    import gzip
    import collections

    groups = collections.defaultdict(set)   # (ENSG, protein) -> {ENST}

    def add(enst, ensg, prot):
        if enst and ensg and prot:
            groups[(ensg, prot)].add(enst)

    def parse(path, pipe_header):
        if not path or not os.path.exists(path):
            return
        enst = ensg = None
        seq = []
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    add(enst, ensg, "".join(seq))
                    seq = []
                    if pipe_header:
                        parts = line[1:].split("|")
                        enst = next((x.split(".")[0] for x in parts if x.startswith("ENST")), None)
                        ensg = next((x.split(".")[0] for x in parts if x.startswith("ENSG")), None)
                    else:
                        g = re.search(r"gene:(ENSG[0-9]+)", line)
                        t = re.search(r"transcript:(ENST[0-9]+)", line)
                        enst, ensg = (t.group(1) if t else None), (g.group(1) if g else None)
                else:
                    seq.append(line.strip())
            add(enst, ensg, "".join(seq))

    parse(config.resolve_shared(config.GENCODE_PROTEIN_FASTA), pipe_header=True)
    parse(config.resolve_shared(config.ENSEMBL_PEP_FASTA), pipe_header=False)
    out = {}
    multi = 0
    for ensts in groups.values():
        u = sorted(ensts)
        if len(u) > 1:
            multi += 1
        for e in u:
            out.setdefault(e, set()).update(u)
    out = {e: sorted(v) for e, v in out.items()}
    logger.info("[iso2function.export] UTR-variant map: %d transcripts (%d protein groups with >1 isoform)",
                len(out), multi)
    return out


def _load_affinity(data_dir):
    """Per-edge affinity (0-1) to simulate kinetics, from the manuscript's quantitative measurements:
    PPI binding affinity = atlas N2H log2 NLR; PDI/TF-target affinity = atlas Y1H-luciferase Log2(FC).
    Each is robust min-maxed to 0-1 (clipped to the 1st-99th percentile). Returns (nlr_map, ppi_median,
    pdival_map, pdi_median) keyed by (atlas SYMBOL-N clone, partner/bait); the median is the imputed
    value for edges with no measurement (binary-only edges and UniProt)."""
    import numpy as np

    def build(fname, k1, k2, valcol):
        df = _read(os.path.join(data_dir, fname), fname)
        if df is None:
            return {}, 0.5
        v = pd.to_numeric(df[valcol], errors="coerce")
        sub = df[v.notna()].copy()
        sub["_v"] = v[v.notna()].to_numpy()
        if not len(sub):
            return {}, 0.5
        lo, hi = np.percentile(sub["_v"], 1), np.percentile(sub["_v"], 99)
        rng = (hi - lo) or 1.0
        m = {(r[k1], r[k2]): round(min(1.0, max(0.0, (r["_v"] - lo) / rng)), 4)
             for _, r in sub.iterrows()}
        med = round(float(np.median(list(m.values()))), 4) if m else 0.5
        return m, med

    nlr, ppi_med = build("ppi_n2h.tsv", "clone_id", "gene_symbol_partner", "log2_nlr")
    pdival, pdi_med = build("pdi_validation.tsv", "clone_id", "Bait", "log2_fc")
    logger.info("[iso2function.affinity] N2H binding affinities=%d (median=%.3f); PDI Log2FC affinities="
                "%d (median=%.3f)", len(nlr), ppi_med, len(pdival), pdi_med)
    return nlr, ppi_med, pdival, pdi_med


def _review(df):
    """Per-row consistency review (logged every run). An interaction is internally consistent when, if
    it has an ENST, it also has an ENSG and (where derivable) a structure + exon coords. Returns the
    issue counts so the caller can iterate until clean."""
    # FIXABLE: if a row has an ENST it must have an ENSG (and a structure must have an ENSG); these are
    # derivable and must be zero.
    fixable = {
        "ENST_but_no_ENSG": int(((df["best_ENST"] != "") & (df["ENSG"] == "")).sum()),
        "structure_but_no_ENSG": int(((df["isoform_structure"] != "") & (df["ENSG"] == "")).sum()),
        "no_id_at_all": int((df["source_isoform_id"] == "").sum()),
    }
    # EXPECTED edge cases (not bugs): UNK atlas clones have genomic exon coords but no Ens91 E/I
    # structure; Ensembl116-only genes have an Ens91 structure string but no Ens91 exon coordinates
    # (the gene postdates the Ens91 exon reference).
    expected = {
        "exon_coords_but_no_structure (UNK atlas clones)":
            int(((df["exon_coordinates"] != "") & (df["isoform_structure"] == "")).sum()),
        "structure_but_no_exon_coords (Ensembl116-only genes)":
            int(((df["isoform_structure"] != "") & (df["exon_coordinates"] == "")).sum()),
    }
    bad = {k: v for k, v in fixable.items() if v}
    if bad:
        logger.warning("[iso2function.review] FIXABLE CONSISTENCY ISSUES: %s", bad)
    else:
        logger.info("[iso2function.review] OK — 0 fixable consistency issues")
    logger.info("[iso2function.review] expected edge cases: %s", {k: v for k, v in expected.items() if v})
    return {"fixable": fixable, "expected": expected}


# --------------------------------------------------------------------------- GFF sources (Vidal2025)
def _load_gff_structures():
    """{GFF clone name -> E/I structure} from the 2025-catalogue GFF transcript_associations."""
    p = config.resolve_shared(config.TFISO_GFF_STRUCTURES)
    out = {}
    if not os.path.exists(p):
        logger.warning("[iso2function.export] GFF structures missing: %s", config.TFISO_GFF_STRUCTURES)
        return out
    with open(p) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                out[parts[3]] = _clean_structure(parts[2])     # E/I only; coord -> ''
    logger.info("[iso2function.export] GFF structures: %d clones (%d with E/I)", len(out),
                sum(1 for v in out.values() if v))
    return out


def _load_gff_exons():
    """{GFF clone name -> 'chrom:strand:(s,e),...'} from the c_6k GFF exon lines (exact, all clones)."""
    p = config.resolve_shared(config.TFISO_GFF)
    by_clone = {}
    if not os.path.exists(p):
        logger.warning("[iso2function.export] c_6k GFF missing: %s", config.TFISO_GFF)
        return {}
    with open(p) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = _TID_RE.search(f[8])
            if not m or "|" not in m.group(1):
                continue
            try:
                s, e = int(f[3]), int(f[4])
            except ValueError:
                continue
            by_clone.setdefault(m.group(1), {"chrom": f[0], "strand": f[6], "ex": []})["ex"].append(
                (min(s, e), max(s, e)))
    out = {}
    for cid, d in by_clone.items():
        ex = sorted(set(d["ex"]))
        out[cid] = f"{d['chrom']}:{d['strand']}:" + ",".join(f"({s},{e})" for s, e in ex)
    logger.info("[iso2function.export] c_6k GFF exon coords: %d clones", len(out))
    return out


def _fmt_exons(got):
    """Format reference_structures.StructureCoords.exons() output as 'chrom:strand:(s,e),...'."""
    if not got:
        return ""
    chrom, strand, exons = got
    return f"{chrom}:{strand}:" + ",".join(f"({s},{e})" for s, e in exons)


def build_interactions_txt(data_dir=None, out_path=None):
    data_dir = data_dir or config.DATA_DIR
    out_path = out_path or os.path.join(data_dir, OUTPUT_NAME)

    # ---- atlas (Vidal2025) maps: clone_id -> gff name / ensg / enst, + GFF structures & exon coords
    cts = _read(os.path.join(data_dir, "clone_to_structure.tsv"), "clone_to_structure")
    atlas = {}
    if cts is not None:
        for _, r in cts.iterrows():
            atlas[r["clone_id"]] = {"gff": r.get("gff_clone_label", ""), "ensg": r.get("ensg", ""),
                                    "enst": r.get("enst", ""), "symbol": r.get("gene_symbol", "")}
    gff_struct = _load_gff_structures()
    gff_exons = _load_gff_exons()

    # ---- Yang / UniProt maps (project mappings) keyed by source isoform id
    def _kv(df, kcol, *cols):
        return {r[kcol]: {c: r.get(c, "") for c in cols} for _, r in df.iterrows()} if df is not None else {}
    ye = _read(config.YANG_ID_TO_ENST, "Yang ENST"); ys = _read(config.YANG_STRUCTURES, "Yang struct")
    yang = {}
    for d, sc, ec in ((ye, "matched_ENST", "Ensembl_gene"), (ys, "matched_ENST", "Ensembl_gene")):
        if d is None:
            continue
        for _, r in d.iterrows():
            y = yang.setdefault(r["Isoform_ID"], {"enst": "", "ensg": "", "structure": "", "symbol": r.get("Gene_Symbol", "")})
            y["enst"] = y["enst"] or r.get("matched_ENST", ""); y["ensg"] = y["ensg"] or r.get("Ensembl_gene", "")
            if "structure" in r:
                y["structure"] = y["structure"] or r.get("structure", "")
    ue = _read(config.UNIPROT_ID_TO_ENST, "UP ENST"); us = _read(config.UNIPROT_STRUCTURES_MAP, "UP struct")
    up = {}
    for d in (ue, us):
        if d is None:
            continue
        for _, r in d.iterrows():
            u = up.setdefault(r["UniProt_ID"], {"enst": "", "ensg": "", "structure": "", "symbol": r.get("Gene_Symbol", "")})
            u["enst"] = u["enst"] or r.get("matched_ENST", ""); u["ensg"] = u["ensg"] or r.get("aa_gene", "")
            if "structure" in r:
                u["structure"] = u["structure"] or r.get("structure", "")

    # ---- structureless isoforms recovered from the MDS-AML full-length catalog
    sr = _read(os.path.join(data_dir, "structureless_resolved.tsv"), "structureless_resolved")
    resolved = {}
    if sr is not None:
        for _, r in sr.iterrows():
            resolved[(r.get("source", ""), r.get("source_isoform_id", ""))] = (r.get("structure", ""), r.get("ensg", ""))

    # universal exon-coordinate source: Ens91 exon reference (token -> genomic coords), works for any
    # Ens91 structure (Yang/UniProt). Atlas uses the exact c_6k GFF clone exons (gff_exons).
    try:
        exon_ref = reference_structures.load_exon_reference()
    except FileNotFoundError as e:
        logger.warning("[iso2function.export] Ens91 exon reference unavailable: %s", e)
        exon_ref = {}
    exon_cache = {}

    def exon_coords_for(ensg, structure):
        if not (ensg and structure):
            return ""
        key = (ensg, structure)
        if key not in exon_cache:
            exon_cache[key] = reference_structures.exon_coords_for_structure(ensg, structure, exon_ref)
        return exon_cache[key]

    nlr, ppi_med, pdival, pdi_med = _load_affinity(data_dir)   # per-edge binding/TF-target affinity (0-1)
    rows = []

    def emit(symbol, ensg, structure, exon, itype, target, source, sid, enst, activity):
        rows.append([symbol, ensg, structure, exon, itype, target, source, sid, enst, activity])

    # ---- Vidal2025 atlas: Y2H + eY1H ----
    for fname, itype, src, clone_col, gene_col, tgt_cols in (
            ("ppi_y2h.tsv", "PPI", "Vidal2025_Y2H", "ad_clone_id", "ad_gene_symbol", ("db_gene_symbol", "db_orf_id")),
            ("pdi_ey1h.tsv", "PDI", "Vidal2025_eY1H", "clone_id", "gene_symbol", ("bait_id",))):
        df = _read(os.path.join(data_dir, fname), fname)
        if df is None:
            continue
        for _, r in df.iterrows():
            if not _is_true(r.get("detected")):
                continue
            cid = r.get(clone_col, ""); a = atlas.get(cid, {})
            gff = a.get("gff", "")
            structure = gff_struct.get(gff, "")
            ensg = _clean_ensg(a.get("ensg", ""))
            tgt = next((r.get(c) for c in tgt_cols if r.get(c)), "")
            act = (nlr.get((cid, tgt), ppi_med) if itype == "PPI"
                   else pdival.get((cid, tgt), pdi_med))
            emit(a.get("symbol") or r.get(gene_col, ""), ensg, structure, gff_exons.get(gff, ""),
                 itype, tgt, src, gff or cid, a.get("enst", ""), act)

    # ---- Vidal2016 (Yang) PPI ----
    p1 = _read(os.path.join(data_dir, "paper1_ppi.tsv"), "paper1_ppi")
    if p1 is not None:
        for _, r in p1.iterrows():
            if not _is_true(r.get("detected")):
                continue
            iid = r.get("Isoform_ID", ""); y = yang.get(iid, {})
            structure, ensg = y.get("structure", ""), y.get("ensg", "")
            if not structure:                                  # recover from cohort full-length match
                structure, ensg2 = resolved.get(("Vidal2016_Y2H", iid), ("", ""))
                ensg = ensg or ensg2
            emit(y.get("symbol") or r.get("Gene_Symbol", ""), _clean_ensg(ensg), _clean_structure(structure),
                 exon_coords_for(_clean_ensg(ensg), _clean_structure(structure)), "PPI",
                 r.get("Interactor_Symbol", ""), "Vidal2016_Y2H", iid, y.get("enst", ""), ppi_med)

    # ---- UniProt curated partners ----
    upf = _read(os.path.join(data_dir, "uniprot_isoform_function.tsv"), "uniprot_isoform_function")
    if upf is not None:
        for _, r in upf.iterrows():
            partners = [p for p in str(r.get("uniprot_partners", "")).split(";") if p]
            if not partners:
                continue
            iid, acc = r.get("isoform_id", ""), r.get("accession", "")
            u = up.get(iid) or up.get(acc) or {}
            sid = iid if iid in up else (acc if acc in up else iid)
            structure, ensg = u.get("structure", ""), u.get("ensg", "")
            if not structure:
                structure, ensg2 = resolved.get(("UniProt", sid), ("", ""))
                ensg = ensg or ensg2
            structure, ensg = _clean_structure(structure), _clean_ensg(ensg)
            exon = exon_coords_for(ensg, structure)
            for tgt in partners:
                emit(u.get("symbol") or r.get("gene_symbol", ""), ensg, structure, exon, "PPI", tgt,
                     "UniProt", sid, u.get("enst", ""), ppi_med)

    df = pd.DataFrame(rows, columns=COLUMNS)
    for c in COLUMNS:
        df[c] = df[c].astype(str).str.replace(r"[\t\r\n]+", " ", regex=True).str.strip()
    df["isoform_structure"] = df["isoform_structure"].map(_clean_structure)
    df["ENSG"] = df["ENSG"].map(_clean_ensg)

    # ---- MANE Select fallback: ANY row still missing an ENST gets the gene's standard reference
    # transcript ("main select isoform as the standard ref"), so every interaction with a known gene
    # carries an ENST. ----
    from .crosswalk import clonelist_resolve
    _mane_ensg, _mane_sym = clonelist_resolve._load_mane()
    no_enst = df["best_ENST"].astype(str).str.strip() == ""
    df.loc[no_enst, "best_ENST"] = df.loc[no_enst, "Symbol"].map(lambda s: _mane_sym.get(s, "")).fillna("")

    # ---- consistency back-fill: every row with an ENST must carry its ENSG, structure, exon coords ----
    # ENST -> (ENSG, structure) from the full Ens91 reference; ENSG also from GENCODE for the rest.
    e2s = reference_structures.load_enst_to_structure()
    enst2ensg = _gencode_enst_to_ensg()
    base = df["best_ENST"].astype(str).str.split(".").str[0]
    ref_gene = base.map(lambda e: e2s.get(e, {}).get("gene", "") if e else "")
    ref_struct = base.map(lambda e: e2s.get(e, {}).get("structure", "") if e else "")
    gc_gene = base.map(lambda e: enst2ensg.get(e, "") if e else "")
    df["ENSG"] = df["ENSG"].where(df["ENSG"] != "", ref_gene.where(ref_gene != "", gc_gene)).map(_clean_ensg)
    df["isoform_structure"] = df["isoform_structure"].where(df["isoform_structure"] != "",
                                                            ref_struct.map(_clean_structure))
    need = (df["isoform_structure"] != "") & (df["exon_coordinates"] == "")
    df.loc[need, "exon_coordinates"] = [exon_coords_for(g, s) for g, s in
                                        zip(df.loc[need, "ENSG"], df.loc[need, "isoform_structure"])]

    # structure fallback: rows still lacking an Ens91 structure (their own transcript isn't in the Ens91
    # reference) inherit the gene's MANE-Select structure (the "main select isoform as the standard ref").
    still = df["isoform_structure"].astype(str).str.strip() == ""
    if still.any():
        mane_struct = df.loc[still, "ENSG"].astype(str).str.split(".").str[0].map(
            lambda g: e2s.get(_mane_ensg.get(g, ""), {}).get("structure", ""))
        df.loc[still, "isoform_structure"] = mane_struct.map(_clean_structure)
        need2 = (df["isoform_structure"] != "") & (df["exon_coordinates"] == "")
        df.loc[need2, "exon_coordinates"] = [exon_coords_for(g, s) for g, s in
                                             zip(df.loc[need2, "ENSG"], df.loc[need2, "isoform_structure"])]

    # ---- ADD additional transcripts (gap-fill only): one ORF maps to several protein-identical
    # (UTR-variant) ENSTs of the same gene. PRESERVE every existing row exactly as built; only APPEND a
    # new row per ADDITIONAL ENST, carrying that transcript's own Ens91 structure. Never modify or
    # replace an existing row's ENST/structure. ----
    utr = _enst_utr_variants()
    extra = []
    for row in df.itertuples(index=False):
        b = str(row.best_ENST).split(".")[0]
        if not b:
            continue
        for v in utr.get(b, []):
            if v == b:
                continue                                   # the original ENST stays as its original row
            # the variant is an ENST ANNOTATION of the same ORF; it carries the ORF's Ens91 structure
            # (NOT a structure pulled from a mixed-build table). Structure stays Ens91; only ENST varies.
            extra.append([row.Symbol, row.ENSG, row.isoform_structure, row.exon_coordinates,
                          row.interaction_type, row.target, row.source, row.source_isoform_id, v,
                          row.activity_score])
    if extra:
        df = pd.concat([df, pd.DataFrame(extra, columns=COLUMNS)], ignore_index=True)

    # ---- ENFORCE Ens91-only structures: isoform_structure must be pure Ens91 (every non novel-site
    # token is a real Ens91 exon of that gene). Structures cannot be mixed across Ensembl builds, so any
    # structure that fails (e.g. v116 numbering that leaked from a reference table, or an Ensembl116-only
    # gene) is BLANKED here, while the best_ENST annotation (any build) is kept. ----
    valid_pairs = set(exon_ref.keys())

    def _is_ens91(ensg, struct):
        if not struct:
            return True
        for t in struct.split("|"):
            if not t or "_" in t:                          # novel splice-site token (gff_process) -> allowed
                continue
            if (ensg, t) not in valid_pairs:
                return False
        return True

    bad = ~df.apply(lambda r: _is_ens91(r["ENSG"], r["isoform_structure"]), axis=1)
    nbad = int(bad.sum())
    if nbad:
        df.loc[bad, "isoform_structure"] = ""
        df.loc[bad, "exon_coordinates"] = ""
        logger.info("[iso2function.export] Ens91 enforcement: blanked %d non-Ens91 structures (ENST kept)", nbad)

    df = df[df["target"] != ""].drop_duplicates().reset_index(drop=True)
    df.to_csv(out_path, sep="\t", index=False)
    _review(df)

    by = df.groupby("source").agg(n=("target", "size"),
                                  struct=("isoform_structure", lambda s: (s != "").sum()),
                                  exon=("exon_coordinates", lambda s: (s != "").sum()),
                                  enst=("best_ENST", lambda s: (s != "").sum())).to_dict("index")
    logger.info("[iso2function.export] %s: %d interactions", out_path, len(df))
    for k, v in sorted(by.items()):
        logger.info("    %-16s n=%-6d ENST=%-6d structure=%-6d exon_coords=%-6d", k, v["n"], v["enst"],
                    v["struct"], v["exon"])
    return df


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config.ensure_dirs()
    build_interactions_txt()
