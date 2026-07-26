#!/usr/bin/env python3
"""variant_extraction.py — supervised per-cell variant -> cell-barcode caller (SNV / indel / insertion
/ FLT3-ITD) for single-cell (KINNEX/Iso-Seq) or bulk RNA-seq BAMs.

Given a mutation-position file + a BAM, report, per cell barcode, whether the cell carries the
reference (WT) or variant (MUT) allele at each locus, plus a per-cell genotype matrix.

CLI (unchanged) + one new flag:
  python3 variant_extraction.py --bam BAM --sample NAME --mutations FILE --reference genome.fa \
                                --output-dir DIR [--bulk] [--min-mapq 20]

Mutation file (tab/space-delimited, one per line; '#' comments allowed):
  chr17   31159086  31159086  NF1c.281T>G   SNV
  chr17   31163273  31163273  NF1c.377delA  Indel
  chr5    171410542 171410542 NPM1_ins      INS      # NPM1c-style small insertion (dominant length)
  chr13   28033939  28034211  FLT3_ITD      ITD      # FLT3-ITD region (start end); soft-clip + I-op
Columns: chrom  start  [end]  label  type. `type` in {SNV(default), Indel, INS, ITD}. For ITD the
3rd column is the region END (JM/TKD1 window); FLT3 exon edges are auto-loaded when the label/gene is FLT3.

Calling ENGINE = the validated altanalyze3.components.bam.variant_scan primitives (reference-free
=/X SNV composition with FASTA fallback, two-major-allele indel, dominant-length insertion,
dual-mechanism FLT3-ITD [long-read insertion-op + short-read soft-clip breakpoint]). This reproduces
the reference genotypes (5801M_pre RUNX1 W279* 256 MUT/435 WT; SRSF2 P95R 2053/1286).

Outputs (per sample; filenames preserved for backward compatibility):
  <sample>_complete_analysis.tsv       per (cell, variant): allele status + read-level evidence
  <sample>_mutation_matrix.csv         per cell: genotype, mutations_present, <label>_WT/_MUT/_GENOTYPE
  <sample>_mutation_matrix_summary.txt  WT/MUT cell counts per locus and overall
Barcodes in the matrix are reverse-complemented (isoform-BAM -> junction/h5ad orientation), as before.

Historical note: the previous version emitted a degenerate strand-bias Fisher test (comparing
identical rows -> p≈1) and several hard-coded placeholder PCR-artifact columns; those are removed.
The SNV rule changed from "any non-reference base = variant" to "dominant alternate base = variant"
(rejects sequencing-error bases); on the validated loci the two agree.
"""
import os
import time
import gzip
import argparse
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("variant_extraction")

try:  # packaged
    from altanalyze3.components.bam.variant_scan import (
        scan_snv, scan_indel, scan_insertion, scan_itd, scan_itd_softclip, rc, core_bc,
        expected_indel_class)
    from altanalyze3.components.bam.codon_mapper import parse_gff_cds
except Exception:  # standalone (run from this directory)
    from variant_scan import (
        scan_snv, scan_indel, scan_insertion, scan_itd, scan_itd_softclip, rc, core_bc,
        expected_indel_class)
    try:
        from codon_mapper import parse_gff_cds
    except Exception:
        parse_gff_cds = None

import pysam


_VTYPES = {"snv": "SNV", "indel": "Indel", "ins": "INS", "itd": "ITD"}


def _norm_chrom(c):
    return c if c.startswith("chr") else "chr" + c


def load_mutations(mutation_files):
    """Parse one or more mutation/panel files into a single locus list (order preserved).

    Accepts a single path or a list, so the panels can be passed directly:
        --mutations panels/positive_control_panel.tsv panels/beataml_hotspot_panel.tsv
    """
    if isinstance(mutation_files, (str, bytes, os.PathLike)):
        mutation_files = [mutation_files]
    muts = []
    for f in mutation_files:
        got = parse_mutation_file(f)
        log.info("  %s: %d loci", os.path.basename(str(f)), len(got))
        muts += got
    return muts


def parse_mutation_file(mutation_file):
    """Return list of dicts: chrom, pos, end, label, vtype, gene, source, expected_uids, edges.

    Two input layouts are accepted and auto-detected:

    1. Simple (whitespace, headerless) -- the historical format:
           chr17  31159086  31159086  NF1c.281T>G  SNV
       3rd column is an optional end coordinate (ITD region).

    2. Panel (tab-delimited, with a header line beginning `chrom<TAB>pos`):
           chrom  pos  vtype  gene  label  source  expected_uids  ref_aa  codon  notes
       `notes` may carry `end=<int>` and `edges=<int,int,...>` for a structural (ITD) target;
       when present these override the annotation-derived exon edges.
    """
    with open(mutation_file) as fh:
        lines = [ln.rstrip("\r\n") for ln in fh]
    body = [ln for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]
    if not body:
        return []
    head = body[0].split("\t")
    if len(head) >= 5 and head[0].strip().lower() == "chrom" and head[1].strip().lower() == "pos":
        return _parse_panel(body[0], body[1:])
    return _parse_simple(body)


def _parse_simple(body):
    muts = []
    for line in body:
        p = line.split()
        if len(p) < 2:
            continue
        chrom = _norm_chrom(p[0])
        pos = int(p[1])
        end = None
        label_i = 2
        if len(p) > 2 and p[2].isdigit():
            end = int(p[2]); label_i = 3
        label = p[label_i] if len(p) > label_i else f"{chrom}_{pos}"
        vtype = "SNV"
        for tok in p[label_i + 1:]:
            if tok.lower() in _VTYPES:
                vtype = _VTYPES[tok.lower()]
                break
        gene = label.split("_")[0].split("p.")[0].split("c.")[0]
        muts.append(dict(chrom=chrom, pos=pos, end=(end or pos), label=label, vtype=vtype,
                         gene=gene, source="", expected_uids="", edges=[]))
    return muts


def _parse_panel(header_line, rows):
    cols = [c.strip() for c in header_line.split("\t")]
    idx = {c: i for i, c in enumerate(cols)}
    muts = []
    for line in rows:
        p = line.split("\t")

        def get(name, default=""):
            i = idx.get(name)
            return p[i].strip() if (i is not None and i < len(p)) else default

        if not get("chrom") or not get("pos"):
            continue
        chrom = _norm_chrom(get("chrom"))
        pos = int(get("pos"))
        vtype = _VTYPES.get(get("vtype").lower(), "SNV")
        gene = get("gene")
        label = get("label") or f"{chrom}_{pos}"
        # structural targets carry their region end and exon edges in `notes`
        end, edges = pos, []
        for kv in get("notes").split(";"):
            kv = kv.strip()
            if kv.startswith("end="):
                try:
                    end = int(kv[4:])
                except ValueError:
                    pass
            elif kv.startswith("edges="):
                edges = [int(x) for x in kv[6:].split(",") if x.strip().isdigit()]
        if not gene:
            gene = label.split("_")[0].split("p.")[0].split("c.")[0]
        muts.append(dict(chrom=chrom, pos=pos, end=end, label=label, vtype=vtype, gene=gene,
                         source=get("source"), expected_uids=get("expected_uids"), edges=edges))
    return muts


# Exon boundaries are needed by the ITD soft-clip caller to reject ordinary splice-junction clips
# (they pile up at exon edges in normal samples too). They are read from the AltAnalyze exon
# annotation shipped with the package, so no flag is required; --exon_annot overrides the file and
# --gff falls back to a GENCODE GFF3. No genomic coordinate is hard-coded anywhere.
_RESOURCE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             "long_read", "resources")
DEFAULT_EXON_ANNOT = os.path.join(_RESOURCE_DIR, "Hs_Ensembl_exon.txt.gz")
DEFAULT_GENE_ANNOT = os.path.join(_RESOURCE_DIR, "Hs_Ensembl-annotations.txt.gz")


def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _symbol_to_ensg(symbols, gene_annot=None):
    """{symbol -> ENSG} for the requested symbols, from Hs_Ensembl-annotations.txt[.gz].

    The table also carries non-Ensembl accessions (e.g. LRG_457 for FLT3); only ENSG rows are kept.
    """
    path = gene_annot or DEFAULT_GENE_ANNOT
    want = set(symbols)
    out = {}
    if not os.path.exists(path):
        log.warning("gene annotation not found: %s", path)
        return out
    with _open_text(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or not p[0].startswith("ENSG"):
                continue
            if p[1] in want and p[1] not in out:
                out[p[1]] = p[0]
    return out


def load_exon_edges(genes, exon_annot=None, gene_annot=None):
    """{gene_symbol -> sorted [exon boundary coordinates]} from the AltAnalyze exon annotation.

    One streaming pass over the exon table for all requested genes. Boundaries are the exon
    region start/stop of every E-block of the gene (I-blocks are introns and are skipped).
    """
    genes = sorted(set(g for g in genes if g))
    if not genes:
        return {}
    path = exon_annot or DEFAULT_EXON_ANNOT
    if not os.path.exists(path):
        log.warning("exon annotation not found: %s -- ITD splice-edge filtering disabled", path)
        return {}
    sym2ensg = _symbol_to_ensg(genes, gene_annot)
    ensg2sym = {v: k for k, v in sym2ensg.items()}
    missing = [g for g in genes if g not in sym2ensg]
    if missing:
        log.warning("no Ensembl gene id for: %s", ",".join(missing))
    edges = {}
    with _open_text(path) as fh:
        header = fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            sym = ensg2sym.get(p[0])
            if sym is None or not p[1].startswith("E"):
                continue
            try:
                edges.setdefault(sym, set()).update((int(p[4]), int(p[5])))
            except ValueError:
                continue
    return {g: sorted(v) for g, v in edges.items()}


def _gff_exon_edges(gff_path, genes):
    """Fallback: exon (CDS) boundaries per gene from a GENCODE GFF3 via codon_mapper."""
    if not (gff_path and parse_gff_cds and os.path.exists(gff_path)):
        return {}
    try:
        parsed = parse_gff_cds(gff_path)
    except Exception as e:
        log.warning("could not parse --gff %s: %s", gff_path, e)
        return {}
    out = {}
    for g in set(genes):
        t = next((v for v in parsed.get(g, {}).values() if v["mane"]), None)
        if t:
            out[g] = sorted(set(b for s, e, _ in t["cds"] for b in (s, e)))
    return out


def call_locus(bam, fa, m, min_mapq, bulk, gene_edges=None):
    """Return (per_cell {core: {'mut_reads','wt_reads'}}, summary dict) for one mutation."""
    chrom, pos, vtype, label = m["chrom"], m["pos"], m["vtype"], m["label"]
    cell = {}
    summ = dict(label=label, gene=m["gene"], chrom=chrom, pos=pos, vtype=vtype)
    if vtype == "SNV":
        s = scan_snv(bam, fa, chrom, pos, min_mapq=min_mapq, bulk=bulk)
        ref_b, alt_b = s["ref_base"], s["alt_base"]
        for c, bc in s["cell_bases"].items():
            mut = bc.get(alt_b, 0) if alt_b else 0
            wt = bc.get(ref_b, 0) if ref_b else 0
            cell[c] = {"mut_reads": mut, "wt_reads": wt}
        note = ""
        if s["fasta_ref"] and ref_b and s["fasta_ref"] != ref_b:
            note = f"fasta_ref={s['fasta_ref']}!=reads_ref={ref_b}"
        summ.update(depth=s["depth"], ref_base=ref_b, alt_base=alt_b, alt_reads=s["alt_reads"],
                    ref_reads=s["ref_reads"], other_reads=s["other_reads"], note=note,
                    vaf=round(s["alt_reads"] / s["depth"], 4) if s["depth"] else 0.0)
    elif vtype == "Indel":
        # supervised: detect the allele class the panel names (ins/del/any), not the most popular one
        s = scan_indel(bam, chrom, pos, min_mapq=min_mapq, bulk=bulk,
                       expected=expected_indel_class(label))
        for c, d in s["cell_alleles"].items():
            cell[c] = {"mut_reads": d["indel"], "wt_reads": d["clean"]}
        summ.update(depth=s["depth"], ref_base="ref", alt_base=s["major_indel"] or "NA",
                    alt_reads=s["n_indel_reads"], ref_reads=s["depth"] - s["n_indel_reads"],
                    other_reads=0,
                    # sorted so the rendered spectrum is byte-identical run to run (a dict's
                    # insertion order here follows set iteration, which is hash-seed dependent)
                    note=str(dict(sorted(s["allele_counts"].items(),
                                         key=lambda kv: (-kv[1], kv[0][0], kv[0][1]))))[:80],
                    vaf=round(s["n_indel_reads"] / s["depth"], 4) if s["depth"] else 0.0)
    elif vtype == "INS":
        s = scan_insertion(bam, chrom, pos, min_len=3, max_len=60, min_mapq=min_mapq,
                           require_dominant=True, dom_tol=1, bulk=bulk)
        for c, d in s["cell"].items():
            cell[c] = {"mut_reads": d["ins"], "wt_reads": d["clean"]}
        summ.update(depth=s["depth"], ref_base="ref",
                    alt_base=(f"ins{s['dom_len']}bp" if s["dom_len"] else "NA"),
                    alt_reads=s["n_ins_reads"], ref_reads=s["depth"] - s["n_ins_reads"],
                    other_reads=0,
                    note=(f"dom_len={s['dom_len']} dom_frac={s['dom_count']}/{s['n_ins_total']} "
                          f"lens={str(s['len_counts'])[:60]}"),
                    vaf=round(s["n_ins_reads"] / s["depth"], 4) if s["depth"] else 0.0)
    elif vtype == "ITD":
        end = m["end"] if m["end"] > pos else pos + 300
        # exon edges of THIS gene (not a single hard-coded gene's), so the soft-clip caller can
        # distinguish an ITD breakpoint from a normal splice-junction clip.
        edges = m.get("edges") or (gene_edges or {}).get(m["gene"]) or []
        if not edges:
            log.warning("%s: no exon boundaries for %s -- ITD soft-clip splice filtering is OFF "
                        "(breakpoints at exon edges may be miscalled)", label, m["gene"])
        s_ins = scan_itd(bam, chrom, pos, end, min_mapq=min_mapq, bulk=bulk)
        s_sc = scan_itd_softclip(bam, fa, chrom, pos, end, edges,
                                 min_mapq=max(1, min_mapq // 20), bulk=bulk)
        use_sc = s_sc["n_itd_reads"] > s_ins["n_itd_reads"]
        cm = s_sc["cell"] if use_sc else s_ins["cell"]
        for c, d in cm.items():
            cell[c] = {"mut_reads": d["itd"], "wt_reads": d["clean"]}
        depth = s_sc["n_span"] if use_sc else s_ins["n_span"]
        n_itd = s_sc["n_itd_reads"] if use_sc else s_ins["n_itd_reads"]
        summ.update(depth=depth, ref_base="ref",
                    alt_base=(f"ITD_softclip@{s_sc['bp']}({'dup' if s_sc['tandem_dup'] else 'ins'})"
                              if use_sc
                              else (f"ITD{s_ins['dom_len']}bp" if s_ins["dom_len"] else "NA")),
                    alt_reads=n_itd, ref_reads=depth - n_itd, other_reads=0,
                    note=(f"insop={s_ins['n_itd_reads']} softclip={s_sc['n_itd_reads']}"
                          f"@{s_sc['bp']} dup={s_sc['tandem_dup']} edges={len(edges)}"),
                    vaf=round(n_itd / depth, 4) if depth else 0.0)
    else:
        summ.update(depth=0, ref_base="NA", alt_base="NA", alt_reads=0, ref_reads=0, vaf=0.0)
    return cell, summ


def genotype_of(mut_reads, wt_reads, mut_min=1, wt_min=1):
    if mut_reads >= mut_min:
        return "MUT"
    if wt_reads >= wt_min:
        return "WT"
    return "UNK"


def run(sample, bam_path, mutation_file, reference, out_dir, min_mapq=20, bulk=False,
        gff_path=None, exon_annot=None, gene_annot=None):
    os.makedirs(out_dir, exist_ok=True)
    muts = load_mutations(mutation_file)
    log.info("Loaded %d mutation positions", len(muts))
    bam = pysam.AlignmentFile(bam_path, "rb")
    fa = pysam.FastaFile(reference) if reference and os.path.exists(reference) else None
    # Exon boundaries for every gene carrying an ITD target. Default source is the annotation
    # shipped with the package; --exon_annot replaces it and --gff is a GENCODE fallback.
    itd_genes = set(m["gene"] for m in muts if m["vtype"] == "ITD")
    gene_edges = {}
    if itd_genes:
        gene_edges = _gff_exon_edges(gff_path, itd_genes) if gff_path else {}
        if not gene_edges:
            gene_edges = load_exon_edges(itd_genes, exon_annot, gene_annot)
        for g in sorted(itd_genes):
            log.info("ITD exon boundaries for %s: %d", g, len(gene_edges.get(g, [])))

    # A panel may list the same variant once per expecting donor (same chrom/pos/label, different
    # expected_uids), and different labels may alias one position. Each distinct locus is therefore
    # scanned ONCE and the result reused: one summary row per input row (so per-donor expectations
    # are preserved), but per-cell rows emitted once per label -- otherwise cells are double counted.
    per_locus_cells = {}   # label -> {core: {'mut_reads','wt_reads'}}
    summaries = []
    read_rows = []
    scan_cache = {}        # (chrom,pos,end,vtype) -> (cells, summary)
    label_pos = {}
    for m in muts:
        key = (m["chrom"], m["pos"], m["end"], m["vtype"])
        if key in scan_cache:
            cells, base_summ = scan_cache[key]
        else:
            cells, base_summ = call_locus(bam, fa, m, min_mapq, bulk, gene_edges)
            scan_cache[key] = (cells, base_summ)
            log.info("  %s %s:%d (%s) ref=%s alt=%s depth=%s vaf=%s -> %d cells",
                     m["label"], m["chrom"], m["pos"], m["vtype"], base_summ.get("ref_base"),
                     base_summ.get("alt_base"), base_summ.get("depth"), base_summ.get("vaf"),
                     len(cells))
        summ = dict(base_summ)
        summ.update(label=m["label"], gene=m["gene"], source=m.get("source", ""),
                    expected_uids=m.get("expected_uids", ""))
        summaries.append(summ)

        if m["label"] in per_locus_cells:
            if label_pos[m["label"]] != (m["chrom"], m["pos"]):
                log.warning("label %s maps to two positions (%s vs %s); per-cell rows keep the first",
                            m["label"], label_pos[m["label"]], (m["chrom"], m["pos"]))
            continue
        per_locus_cells[m["label"]] = cells
        label_pos[m["label"]] = (m["chrom"], m["pos"])
        for core, d in cells.items():
            gt = genotype_of(d["mut_reads"], d["wt_reads"])
            read_rows.append(dict(cell_barcode=rc(core), bam_barcode=core, chr=m["chrom"],
                                  position=m["pos"], label=m["label"], gene=m["gene"], vtype=m["vtype"],
                                  ref_base=summ.get("ref_base") or "", alt_base=summ.get("alt_base") or "",
                                  mut_reads=d["mut_reads"], wt_reads=d["wt_reads"],
                                  allele_status=("V" if gt == "MUT" else ("R" if gt == "WT" else "UNK")),
                                  genotype=gt, locus_depth=summ.get("depth"),
                                  locus_alt_reads=summ.get("alt_reads"), locus_vaf=summ.get("vaf")))
    bam.close()
    if fa:
        fa.close()

    # --- per-read/per-cell forensic table ---
    tsv = os.path.join(out_dir, f"{sample}_complete_analysis.tsv")
    _write_tsv(read_rows, tsv)
    log.info("wrote %s", tsv)

    # --- per-locus read counts (machine readable; one row per input row) ---
    # This is the cross-sample aggregation input: column set is a contract, do not reorder.
    rc_cols = ["uid", "chrom", "pos", "gene", "label", "vtype", "source", "expected",
               "depth", "ref_base", "ref_reads", "alt_base", "alt_reads", "other_reads",
               "vaf", "n_cells_mut", "n_cells_wt", "note"]
    rc_rows = []
    for s in summaries:
        cells = per_locus_cells.get(s["label"], {})
        n_mut = sum(1 for d in cells.values() if genotype_of(d["mut_reads"], d["wt_reads"]) == "MUT")
        n_wt = sum(1 for d in cells.values() if genotype_of(d["mut_reads"], d["wt_reads"]) == "WT")
        exp_uids = s.get("expected_uids", "")
        expected = "yes" if (exp_uids == "ALL" or sample in exp_uids.split(";")) else "no"
        rc_rows.append({"uid": sample, "chrom": s["chrom"], "pos": s["pos"], "gene": s["gene"],
                        "label": s["label"], "vtype": s["vtype"], "source": s.get("source", ""),
                        "expected": expected, "depth": s.get("depth", 0),
                        "ref_base": s.get("ref_base") or "", "ref_reads": s.get("ref_reads", 0),
                        "alt_base": s.get("alt_base") or "", "alt_reads": s.get("alt_reads", 0),
                        "other_reads": s.get("other_reads", 0), "vaf": s.get("vaf", 0.0),
                        "n_cells_mut": n_mut, "n_cells_wt": n_wt, "note": s.get("note", "")})
    rc_path = os.path.join(out_dir, f"{sample}_variant_readcounts.tsv")
    _write_tsv(rc_rows, rc_path, fieldnames=rc_cols)
    log.info("wrote %s", rc_path)

    # --- per-cell genotype matrix (RC barcodes) ---
    labels = list(dict.fromkeys(m["label"] for m in muts))   # ordered-unique: one column set per label
    all_cores = sorted(set(c for cells in per_locus_cells.values() for c in cells))
    matrix_rows = []
    for core in all_cores:
        row = {"cell_barcode": rc(core)}
        gts = []
        present = []
        for lab in labels:
            d = per_locus_cells[lab].get(core, {"mut_reads": 0, "wt_reads": 0})
            gt = genotype_of(d["mut_reads"], d["wt_reads"])
            row[f"{lab}_WT"] = d["wt_reads"]
            row[f"{lab}_MUT"] = d["mut_reads"]
            row[f"{lab}_GENOTYPE"] = gt
            gts.append(gt)
            if gt == "MUT":
                present.append(lab)
        row["genotype"] = "MUT" if "MUT" in gts else ("WT" if "WT" in gts else "UNK")
        row["mutations_present"] = ";".join(present) if present else "None"
        matrix_rows.append(row)
    # column order: cell_barcode, genotype, mutations_present, then per-label triples
    cols = ["cell_barcode", "genotype", "mutations_present"]
    for lab in labels:
        cols += [f"{lab}_WT", f"{lab}_MUT", f"{lab}_GENOTYPE"]
    mcsv = os.path.join(out_dir, f"{sample}_mutation_matrix.csv")
    _write_csv(matrix_rows, mcsv, cols)
    log.info("wrote %s", mcsv)

    # --- summary ---
    summ_path = os.path.join(out_dir, f"{sample}_mutation_matrix_summary.txt")
    with open(summ_path, "w") as fh:
        fh.write(f"Mutation Matrix Summary - {sample}\n" + "=" * 50 + "\n\n")
        fh.write(f"Total genotyped cells: {len(matrix_rows)}\n")
        fh.write(f"Total loci: {len(labels)}\n\n")
        overall = {"MUT": 0, "WT": 0, "UNK": 0}
        for r in matrix_rows:
            overall[r["genotype"]] += 1
        fh.write("Final genotype (per cell, any-locus MUT dominates):\n")
        for g in ("MUT", "WT", "UNK"):
            n = overall[g]
            fh.write(f"  {g}: {n} cells ({100*n/len(matrix_rows):.1f}%)\n" if matrix_rows else f"  {g}: 0\n")
        fh.write("\nPer-locus (read-genotyped cells):\n")
        for lab in labels:
            cells = per_locus_cells[lab]
            mut = sum(1 for d in cells.values() if genotype_of(d["mut_reads"], d["wt_reads"]) == "MUT")
            wt = sum(1 for d in cells.values() if genotype_of(d["mut_reads"], d["wt_reads"]) == "WT")
            s = next(x for x in summaries if x["label"] == lab)
            fh.write(f"  {lab} ({s['chrom']}:{s['pos']} {s['vtype']}): "
                     f"{mut+wt} cells  WT={wt}  MUT={mut}  depth={s.get('depth')}  "
                     f"ref={s.get('ref_base')} alt={s.get('alt_base')} VAF={s.get('vaf')}\n")
    log.info("wrote %s", summ_path)
    return matrix_rows, summaries


_MOUNT_ALIASES = [("/data/salomonis-archive/", "/Volumes/salomonis-archive/"),
                  ("/data/salomonis2/", "/Volumes/salomonis2/"),
                  ("/data/saljh8/", "/Users/saljh8/")]


def resolve_bam(path):
    """Return the path if it exists, else the first existing cluster<->local mount alias of it."""
    if os.path.exists(path):
        return path
    for a, b in _MOUNT_ALIASES:
        for src, dst in ((a, b), (b, a)):
            if path.startswith(src):
                alt = dst + path[len(src):]
                if os.path.exists(alt):
                    return alt
    return None


def run_from_metadata(metadata, mutation_files, reference, out_dir, min_mapq=20, bulk=False,
                      gff_path=None, exon_annot=None, gene_annot=None, samples=None):
    """Scan every sample in an sclr-style metadata file (columns: uid, bam, ...).

    One BAM per uid is scanned (the first that exists; alternates are counted in the manifest).
    Writes the same per-sample outputs as run(), plus input_manifest.tsv recording exactly which
    BAM produced each sample's results.
    """
    import csv as _csv
    os.makedirs(out_dir, exist_ok=True)
    order, uid_bams = [], {}
    with open(metadata) as fh:
        for r in _csv.DictReader(fh, delimiter="\t"):
            uid = (r.get("uid") or "").strip()
            bam = (r.get("bam") or "").strip()
            if not uid or not bam:
                continue
            if uid not in uid_bams:
                uid_bams[uid] = []; order.append(uid)
            uid_bams[uid].append(bam)
    if samples:
        want = set(samples)
        order = [u for u in order if u in want]
    log.info("metadata %s -> %d samples", metadata, len(order))

    manifest, failed = [], []
    for i, uid in enumerate(order, 1):
        found = [p for p in (resolve_bam(b) for b in uid_bams[uid]) if p]
        if not found:
            log.warning("[%d/%d] %s: no BAM on disk -- skipped", i, len(order), uid)
            failed.append((uid, "no BAM on disk")); continue
        t0 = time.time()
        try:
            run(uid, found[0], mutation_files, reference, out_dir, min_mapq=min_mapq, bulk=bulk,
                gff_path=gff_path, exon_annot=exon_annot, gene_annot=gene_annot)
            manifest.append({"uid": uid, "bam": found[0], "n_alt_bams": len(found) - 1,
                             "seconds": round(time.time() - t0, 1)})
            log.info("[%d/%d] %s done in %.0fs", i, len(order), uid, time.time() - t0)
        except Exception as e:
            log.exception("[%d/%d] %s FAILED: %s", i, len(order), uid, e)
            failed.append((uid, repr(e)))
    _write_tsv(manifest, os.path.join(out_dir, "input_manifest.tsv"),
               fieldnames=["uid", "bam", "n_alt_bams", "seconds"])
    log.info("cohort done: %d/%d samples -> %s", len(manifest), len(order), out_dir)
    for uid, why in failed:
        log.warning("  FAILED %s: %s", uid, why)
    return manifest, failed


def _write_tsv(rows, path, fieldnames=None):
    import csv
    if not rows:
        open(path, "w").close(); return
    with open(path, "w", newline="") as fh:
        # lineterminator="\n": csv defaults to "\r\n"; with newline="" that carriage return is
        # written literally and every awk/cut/R parse of the LAST column silently fails.
        w = csv.DictWriter(fh, fieldnames=fieldnames or list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def _write_csv(rows, path, cols):
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def parse_arguments():
    ap = argparse.ArgumentParser(description="Supervised variant/indel/insertion/ITD -> cell barcode reporting")
    ap.add_argument("--sample", "--s", default=None, help="Sample name (single-sample mode)")
    ap.add_argument("--bam", "--b", default=None, help="Path to BAM file (index required)")
    ap.add_argument("--metadata", default=None,
                    help="sclr metadata TSV (uid, bam, ...): scan EVERY sample in it. "
                         "Use instead of --bam/--sample.")
    ap.add_argument("--samples", nargs="*", default=None,
                    help="With --metadata: restrict to these uids (default: all)")
    ap.add_argument("--mutations", "--m", required=True, nargs="+",
                    help="One or more mutation/panel files (simple or headered panel format)")
    ap.add_argument("--reference", "--r", required=True, help="Reference FASTA (.fai required)")
    ap.add_argument("--output-dir", "--o", required=True, help="Output directory")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--bulk", action="store_true",
                    help="bulk BAM without cell barcodes (e.g. STAR RNA-seq): count all reads as one pseudo-cell")
    ap.add_argument("--exon_annot", dest="exon_annot", default=None,
                    help=f"Exon annotation for ITD splice-edge filtering "
                         f"(default: bundled {os.path.basename(DEFAULT_EXON_ANNOT)})")
    ap.add_argument("--gene_annot", dest="gene_annot", default=None,
                    help=f"Gene symbol -> Ensembl id table "
                         f"(default: bundled {os.path.basename(DEFAULT_GENE_ANNOT)})")
    ap.add_argument("--gff", default=None,
                    help="GENCODE GFF3 alternative to --exon_annot for exon edges in ITD mode")
    return ap.parse_args()


def main():
    args = parse_arguments()
    if args.metadata:
        run_from_metadata(args.metadata, args.mutations, args.reference, args.output_dir,
                          min_mapq=args.min_mapq, bulk=args.bulk, gff_path=args.gff,
                          exon_annot=args.exon_annot, gene_annot=args.gene_annot,
                          samples=args.samples)
        return
    if not (args.sample and args.bam):
        raise SystemExit("provide --metadata (cohort) OR both --sample and --bam (one sample)")
    run(args.sample, args.bam, args.mutations, args.reference, args.output_dir,
        min_mapq=args.min_mapq, bulk=args.bulk, gff_path=args.gff,
        exon_annot=args.exon_annot, gene_annot=args.gene_annot)


if __name__ == "__main__":
    main()
