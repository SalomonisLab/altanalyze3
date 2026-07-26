#!/usr/bin/env python3.11
"""variant_scan.py — supervised per-cell variant scanner for long-read single-cell KINNEX BAMs.

TEST COPY of the variant-calling tool (developed here per instructions before any change to the
canonical altanalyze3/components/bam/variant_extraction.py). Built on the VALIDATED reference-free
=/X CIGAR genotyper (altanalyze3.variant_impact.genotype_from_bam, which reproduced the documented
5801M_pre RUNX1 256/435 and SRSF2 2053/1286).

What it adds over the two existing callers, for THIS task:
  * SNV : full per-base composition + depth + VAF at each locus (specificity signal), not just the
          two-major-allele reduction. ref = consensus of '=' (match) reads; alt = dominant non-ref.
  * Indel: major indel (type,len) within a window vs clean spanning reads, per cell.
  * INS  : targeted small-insertion scan (NPM1 exon-12 4bp dup) — insertion CIGAR within a window.
  * Per-cell barcodes reported in BOTH orientations (bam + reverse-complement) so they join to the
    isoform h5ad (RC) or the raw BAM as needed.
  * One BAM opened once; all panel loci fetched by index (tiny region reads) -> fast.

FASTA (genome.fa) here is NOT chr-prefixed; BAM/panel ARE. We strip 'chr' for FASTA fetches only.
"""
import os, sys, csv, argparse, logging
from collections import defaultdict, Counter
import pysam

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("variant_scan")
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
BARCODE_TAGS = ("CB", "BC", "XC", "UB")


def rc(seq):
    return "".join(_COMP.get(b, b) for b in reversed(str(seq)))


def core_bc(barcode):
    import re
    bc = str(barcode).strip().split("-")[0].split(".")[0]
    m = re.match(r"^[ACGTN]+", bc.upper())
    return m.group(0) if m else bc.upper()


def get_bc(read, bulk=False):
    for tg in BARCODE_TAGS:
        try:
            return read.get_tag(tg)
        except KeyError:
            continue
    # bulk (no single-cell barcode, e.g. BEAT-AML STAR RNA-seq): one pseudo-cell so every read
    # is counted; VAF/depth/indel/ITD logic is identical, only the per-cell dimension collapses.
    return "BULK" if bulk else None


def base_and_class_at(read, target0):
    """(query_base, '=' | 'X' | 'M') at reference position target0 (0-based), or (None,None)."""
    refpos = read.reference_start
    qpos = 0
    seq = read.query_sequence
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            if refpos <= target0 < refpos + length:
                qi = qpos + (target0 - refpos)
                base = seq[qi].upper() if seq and qi < len(seq) else None
                return base, {0: "M", 7: "=", 8: "X"}[op]
            refpos += length; qpos += length
        elif op == 1:
            qpos += length
        elif op in (2, 3):
            refpos += length
        elif op == 4:
            qpos += length
        if refpos > target0:
            break
    return None, None


def indels_in_window(read, target0, window):
    """[(type,len,refpos)] indels whose ref position is within `window` of target0."""
    out = []
    refpos = read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            refpos += length
        elif op == 1:  # insertion
            if abs(refpos - target0) <= window:
                out.append(("ins", length, refpos))
        elif op == 2:  # deletion
            if refpos - window <= target0 <= refpos + length + window:
                out.append(("del", length, refpos))
            refpos += length
        elif op == 3:
            refpos += length
    return out


def scan_snv(bam, fa, chrom, pos1, min_mapq=20, bulk=False):
    """Return dict: depth, ref_base, ref_reads, alt_base, alt_reads, other_reads, fasta_ref,
    and per-cell {core: Counter(base)} covering the locus."""
    target0 = pos1 - 1
    base_counts = Counter()
    cls_counts = Counter()          # '=', 'X', 'M'
    cell_bases = defaultdict(Counter)
    for read in bam.fetch(chrom, target0, target0 + 1):
        if read.is_unmapped or read.mapping_quality < min_mapq or read.cigartuples is None:
            continue
        base, cls = base_and_class_at(read, target0)
        if base is None:
            continue
        bc = get_bc(read, bulk)
        if not bc:
            continue
        base_counts[base] += 1
        cls_counts[cls] += 1
        cell_bases[core_bc(bc)][base] += 1
    fasta_ref = None
    if fa is not None:
        fc = chrom[3:] if chrom.startswith("chr") else chrom
        try:
            fasta_ref = fa.fetch(fc, target0, target0 + 1).upper()
        except Exception:
            fasta_ref = None
    depth = sum(base_counts.values())
    # ref = consensus of '=' reads; if BAM had them. Fall back to FASTA ref, then majority base.
    ref_base = fasta_ref
    if ref_base is None and depth:
        ref_base = base_counts.most_common(1)[0][0]
    alt_candidates = [(b, c) for b, c in base_counts.items() if b != ref_base]
    alt_candidates.sort(key=lambda x: -x[1])
    alt_base = alt_candidates[0][0] if alt_candidates else None
    ref_reads = base_counts.get(ref_base, 0)
    alt_reads = base_counts.get(alt_base, 0) if alt_base else 0
    other_reads = depth - ref_reads - alt_reads
    return {"depth": depth, "ref_base": ref_base, "ref_reads": ref_reads,
            "alt_base": alt_base, "alt_reads": alt_reads, "other_reads": other_reads,
            "fasta_ref": fasta_ref, "n_match": cls_counts.get("=", 0), "n_mismatch": cls_counts.get("X", 0),
            "n_plainM": cls_counts.get("M", 0), "cell_bases": cell_bases}


def expected_indel_class(label):
    """('ins' | 'del' | 'any' | None) — the indel class this supervised locus is looking for.

    Taken from the HGVS in the panel label, because a supervised scan already knows the variant:
      c.487_488ins / c.2257_2258insG / dup -> 'ins'
      c.377delA / c.385delG               -> 'del'
      c.1992_1994delATTinsGGA (delins)    -> 'any' (both events are the variant)
      p.G643fs (frameshift, class unstated) -> 'any' (any indel frameshifts)
      no indel HGVS (e.g. c.516A>C)       -> None (fall back to the major-allele model)
    """
    if not label:
        return None
    s = str(label).lower()
    if "delins" in s or "fs" in s:
        return "any"
    has_ins = "ins" in s or "dup" in s
    has_del = "del" in s
    if has_ins and has_del:
        return "any"
    if has_ins:
        return "ins"
    if has_del:
        return "del"
    return None


def scan_indel(bam, chrom, pos1, window=5, min_mapq=20, bulk=False, expected=None):
    """MAJOR indel allele (type,len) within window vs clean spanning reads, per cell.

    Two-major-allele model (matches the validated genotype_from_bam): a cell is MUT only if it
    carries the SINGLE most-common indel allele; a cell with a *different* (minor/error) indel is
    neither MUT nor WT. This suppresses PacBio homopolymer indel errors that would otherwise inflate
    MUT at sites like ASXL1 c.1934dupG."""
    target0 = pos1 - 1
    allele_counts = Counter()
    reads = []  # (core, set_of_alleles, spans_clean_bool)
    fs = max(0, target0 - window); fe = target0 + window + 1
    for read in bam.fetch(chrom, fs, fe):
        if read.is_unmapped or read.mapping_quality < min_mapq or read.cigartuples is None:
            continue
        bc = get_bc(read, bulk)
        if not bc:
            continue
        inds = set((t, l) for t, l, _ in indels_in_window(read, target0, window))
        base, cls = base_and_class_at(read, target0)
        for a in inds:
            allele_counts[a] += 1
        reads.append((core_bc(bc), inds, base is not None))
    # SUPERVISED rule (default when the locus states its variant class): a read is mutant if it
    # carries an indel of the class we were asked to detect. A cell with the variant is MUT; a cell
    # with only clean spanning reads is WT; a cell carrying only some OTHER indel is neither, and is
    # not reported. This replaces a popularity contest that discarded cells carrying the very variant
    # under test whenever a competing allele tied or outnumbered it.
    #
    # UNSUPERVISED fallback (expected is None): single major allele, deterministically tie-broken.
    # `inds` is a set of (type,len) tuples containing strings, so its iteration order -- hence Counter
    # insertion order -- varies with per-process string hash randomization, and most_common() breaks
    # ties by insertion order. Sort explicitly: highest count, then type, then shortest length.
    def _matches(inds):
        if expected == "any":
            return bool(inds)
        return any(t == expected for t, _l in inds)

    major = (min(allele_counts.items(), key=lambda kv: (-kv[1], kv[0][0], kv[0][1]))[0]
             if allele_counts else None)
    cell_alleles = defaultdict(lambda: {"indel": 0, "clean": 0})
    n_mut_reads = 0
    n_span = 0
    for c, inds, clean in reads:
        if inds or clean:
            n_span += 1
        hit = _matches(inds) if expected else (major is not None and major in inds)
        if hit:
            cell_alleles[c]["indel"] += 1
            n_mut_reads += 1
        elif not inds and clean:
            cell_alleles[c]["clean"] += 1
    if expected:
        obs = sorted((kv for kv in allele_counts.items()
                      if expected == "any" or kv[0][0] == expected),
                     key=lambda kv: (-kv[1], kv[0][0], kv[0][1]))
        called = f"{obs[0][0][0]}{obs[0][0][1]}" if obs else f"{expected}0"
    else:
        called = f"{major[0]}{major[1]}" if major else None
    return {"depth": n_span, "major_indel": called, "expected": expected,
            "n_indel_reads": n_mut_reads, "allele_counts": dict(allele_counts),
            "cell_alleles": cell_alleles}


def scan_insertion(bam, chrom, pos1, window=30, min_len=2, max_len=60, min_mapq=20,
                   require_dominant=True, dom_tol=0, bulk=False):
    """Targeted insertion scan (e.g. NPM1 exon-12 4bp dup, FLT3-ITD).

    A real recurrent insertion event has a DOMINANT length (NPM1c = 4 bp; an ITD = one specific
    duplication size). Scattered lengths at low fraction are PacBio insertion errors (severe at very
    highly-expressed genes like NPM1). With require_dominant, a cell is INS only if it carries an
    insertion whose length is within dom_tol of the single most-common insertion length. We report
    the dominant length + its fraction so the caller can gate on a clean event vs an error smear."""
    target0 = pos1 - 1
    len_counts = Counter()
    read_ins = []   # (core, best_ins_len_or_None, spans_anchor_bool)
    fs = max(0, target0 - window); fe = target0 + window + 1
    for read in bam.fetch(chrom, fs, fe):
        if read.is_unmapped or read.mapping_quality < min_mapq or read.cigartuples is None:
            continue
        bc = get_bc(read, bulk)
        if not bc:
            continue
        c = core_bc(bc)
        refpos = read.reference_start
        found = None
        for op, length in read.cigartuples:
            if op in (0, 7, 8):
                refpos += length
            elif op == 1:
                if min_len <= length <= max_len and abs(refpos - target0) <= window:
                    found = length if found is None else max(found, length)  # bulk-aware via get_bc
            elif op in (2, 3):
                refpos += length
        base, _ = base_and_class_at(read, target0)
        if found:
            len_counts[found] += 1
        read_ins.append((c, found, base is not None))
    dom_len = len_counts.most_common(1)[0][0] if len_counts else None
    dom_count = len_counts.get(dom_len, 0) if dom_len else 0
    cell = defaultdict(lambda: {"ins": 0, "clean": 0})
    n_ins_reads = 0
    for c, found, clean in read_ins:
        is_event = (found is not None and (not require_dominant or
                    (dom_len is not None and abs(found - dom_len) <= dom_tol)))
        if is_event:
            cell[c]["ins"] += 1
            n_ins_reads += 1
        elif found is None and clean:
            cell[c]["clean"] += 1
    depth = n_ins_reads + sum(v["clean"] for v in cell.values())
    return {"len_counts": dict(len_counts), "n_ins_reads": n_ins_reads, "dom_len": dom_len,
            "dom_count": dom_count, "n_ins_total": sum(len_counts.values()),
            "cell": cell, "depth": depth}


def scan_itd(bam, chrom, start1, end1, min_len=6, max_len=400, min_mapq=20, min_dom_frac=0.0, bulk=False):
    """FLT3-ITD-style scan: an internal tandem duplication appears as an INSERTION (CIGAR I) of
    length min_len..max_len anywhere in the [start1,end1] CDS window (the FLT3 juxtamembrane / TKD1
    region). Reports the dominant insertion length and per-cell ITD membership (length within 3bp of
    the dominant, to absorb alignment jitter of the same duplication)."""
    s0 = start1 - 1; e0 = end1
    len_counts = Counter()
    read_ins = []  # (core, ins_len_or_None, spans_bool)
    for read in bam.fetch(chrom, s0, e0):
        if read.is_unmapped or read.mapping_quality < min_mapq or read.cigartuples is None:
            continue
        bc = get_bc(read, bulk)
        if not bc:
            continue
        c = core_bc(bc)
        refpos = read.reference_start
        best = None
        spans = False
        for op, length in read.cigartuples:
            if op in (0, 7, 8):
                if refpos < e0 and refpos + length > s0:
                    spans = True
                refpos += length
            elif op == 1:
                if min_len <= length <= max_len and s0 <= refpos <= e0:
                    best = length if best is None else max(best, length)
            elif op in (2, 3):
                refpos += length
        if best:
            len_counts[best] += 1
        read_ins.append((c, best, spans))
    dom_len = len_counts.most_common(1)[0][0] if len_counts else None
    cell = defaultdict(lambda: {"itd": 0, "clean": 0})
    n_itd = 0
    for c, best, spans in read_ins:
        if best is not None and dom_len is not None and abs(best - dom_len) <= 3:
            cell[c]["itd"] += 1; n_itd += 1
        elif best is None and spans:
            cell[c]["clean"] += 1
    n_span = sum(1 for _, _, sp in read_ins if sp) + n_itd
    return {"len_counts": dict(len_counts), "dom_len": dom_len, "n_itd_reads": n_itd,
            "n_span": n_span, "cell": cell}


def scan_itd_softclip(bam, fa, chrom, start1, end1, exon_edges, min_clip=10, seed=12,
                      min_off=3, max_off=250, min_mapq=1, min_bp=8, edge_tol=6, bulk=False):
    """FLT3-ITD from SHORT-READ (STAR) BAMs, where an ITD is NOT an insertion op but a pile-up of
    SOFT-CLIPPED reads at the duplication breakpoint (verified on BEAT-AML: ins>=10bp = 0 in every
    ITD+ sample; the signal is soft-clips). Logic:
      * collect soft-clips (>= min_clip) whose breakpoint is inside [start1,end1] (the JM/TKD1 window)
      * cluster breakpoints (+/-3 bp); EXCLUDE any within edge_tol of an exon boundary (those are
        normal splice soft-clips -- they pile up at the exon edge in controls too)
      * the dominant remaining (mid-exon) breakpoint with >= min_bp reads is the ITD; a tandem-
        duplication check (clip re-matches local reference at offset min_off..max_off) is reported as
        supporting evidence but NOT required (non-templated-insertion ITDs still pile up).
    Returns dominant breakpoint, supporting reads, VAF (reads/spanning), tandem_dup flag, per cell."""
    s0 = start1 - 1; e0 = end1
    fc = chrom[3:] if chrom.startswith("chr") else chrom
    pad = max_off + 50
    try:
        ref = fa.fetch(fc, s0 - pad, e0 + pad).upper() if fa is not None else ""
    except Exception:
        ref = ""
    roff = s0 - pad
    from collections import defaultdict as _dd
    spans = 0
    bp = _dd(list)   # breakpoint pos -> [(core, side, clip_seq)]
    for read in bam.fetch(chrom, s0, e0):
        if read.is_unmapped or read.mapping_quality < min_mapq or read.cigartuples is None:
            continue
        seq = read.query_sequence
        if not seq:
            continue
        bc = get_bc(read, bulk)
        if not bc:
            continue
        c = core_bc(bc)
        spans += 1
        ct = read.cigartuples
        if ct[0][0] == 4 and ct[0][1] >= min_clip and read.reference_start is not None \
                and s0 <= read.reference_start <= e0:
            bp[read.reference_start].append((c, "L", seq[:ct[0][1]]))
        if ct[-1][0] == 4 and ct[-1][1] >= min_clip and read.reference_end is not None \
                and s0 <= read.reference_end <= e0:
            bp[read.reference_end].append((c, "R", seq[-ct[-1][1]:]))

    def near_edge(p):
        return any(abs(p - ed) <= edge_tol for ed in exon_edges)

    best = None
    for center in sorted(bp, key=lambda p: -len(bp[p])):
        if near_edge(center):
            continue
        clips = [x for p in bp if abs(p - center) <= 3 for x in bp[p]]
        if len(clips) < min_bp:
            continue
        dup = 0
        for _c, side, cseq in clips:
            s = cseq[:seed] if side == "R" else cseq[-seed:]
            j = ref.find(s) if ref else -1
            if j >= 0 and min_off <= abs((roff + j) - center) <= max_off:
                dup += 1
        if best is None or len(clips) > len(best[1]):
            best = (center, clips, dup)
    cell = defaultdict(lambda: {"itd": 0, "clean": 0})
    if best is None:
        return {"bp": None, "n_itd_reads": 0, "n_span": spans, "vaf": 0.0, "tandem_dup": False,
                "cell": cell}
    center, clips, dup = best
    itd_cores = Counter(c for c, _s, _q in clips)
    for c, n in itd_cores.items():
        cell[c]["itd"] += n
    n_itd = len(clips)
    return {"bp": center, "n_itd_reads": n_itd, "n_span": spans,
            "vaf": round(n_itd / spans, 4) if spans else 0.0,
            "tandem_dup": dup >= 0.4 * n_itd, "cell": cell}


def load_panel(path):
    rows = []
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for x in r:
            rows.append(x)
    return rows


def load_cell_states(path):
    if not path or not os.path.exists(path):
        return None
    m = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or "\t" not in line:
                continue
            p = line.split("\t")
            m[core_bc(p[0])] = p[1]
            m[core_bc(rc(p[0]))] = p[1]   # store both orientations
    return m


def scan_bam(uid, bam_path, panel, fa, out_dir, states=None, min_mapq=20,
             mut_min=1, wt_min=1, bulk=False):
    os.makedirs(out_dir, exist_ok=True)
    bam = pysam.AlignmentFile(bam_path, "rb")
    refset = set(bam.references)
    summary_rows = []
    barcode_rows = []
    for v in panel:
        chrom = v["chrom"]; pos = int(v["pos"]); vtype = v["vtype"]
        gene = v["gene"]; label = v["label"]; source = v["source"]
        exp = v.get("expected_uids", "")
        expected = "yes" if (exp == "ALL" or uid in exp.split(";")) else "no"
        if chrom not in refset:
            summary_rows.append(dict(uid=uid, chrom=chrom, pos=pos, gene=gene, label=label,
                                     vtype=vtype, source=source, expected=expected, depth=0,
                                     ref_base="NA", ref_reads=0, alt_base="NA", alt_reads=0,
                                     other_reads=0, vaf=0.0, n_cells_mut=0, n_cells_wt=0,
                                     note="contig_absent"))
            continue
        if vtype == "SNV":
            s = scan_snv(bam, fa, chrom, pos, min_mapq=min_mapq, bulk=bulk)
            depth = s["depth"]; ref_b = s["ref_base"]; alt_b = s["alt_base"]
            n_mut = n_wt = 0
            for cbc, bc_counts in s["cell_bases"].items():
                a = bc_counts.get(alt_b, 0) if alt_b else 0
                r_ = bc_counts.get(ref_b, 0) if ref_b else 0
                if a >= mut_min:
                    n_mut += 1
                    barcode_rows.append(dict(uid=uid, gene=gene, label=label, vtype=vtype,
                                             source=source, expected=expected, bam_barcode=cbc,
                                             rc_barcode=rc(cbc),
                                             cell_state=(states.get(cbc) if states else ""),
                                             mut_reads=a, wt_reads=r_, alt_base=alt_b, ref_base=ref_b))
                elif r_ >= wt_min:
                    n_wt += 1
            vaf = s["alt_reads"] / depth if depth else 0.0
            note = ""
            if s["fasta_ref"] and ref_b and s["fasta_ref"] != ref_b:
                note = f"fasta_ref={s['fasta_ref']}!=reads_ref={ref_b}"
            summary_rows.append(dict(uid=uid, chrom=chrom, pos=pos, gene=gene, label=label,
                                     vtype=vtype, source=source, expected=expected, depth=depth,
                                     ref_base=ref_b, ref_reads=s["ref_reads"], alt_base=alt_b,
                                     alt_reads=s["alt_reads"], other_reads=s["other_reads"],
                                     vaf=round(vaf, 4), n_cells_mut=n_mut, n_cells_wt=n_wt, note=note))
        elif vtype in ("Indel",):
            s = scan_indel(bam, chrom, pos, window=5, min_mapq=min_mapq, bulk=bulk)
            n_mut = n_wt = 0
            for cbc, d in s["cell_alleles"].items():
                if d["indel"] >= mut_min:
                    n_mut += 1
                    barcode_rows.append(dict(uid=uid, gene=gene, label=label, vtype=vtype,
                                             source=source, expected=expected, bam_barcode=cbc,
                                             rc_barcode=rc(cbc),
                                             cell_state=(states.get(cbc) if states else ""),
                                             mut_reads=d["indel"], wt_reads=d["clean"],
                                             alt_base=s["major_indel"] or "", ref_base="ref"))
                elif d["clean"] >= wt_min:
                    n_wt += 1
            depth = s["depth"]
            vaf = s["n_indel_reads"] / depth if depth else 0.0
            summary_rows.append(dict(uid=uid, chrom=chrom, pos=pos, gene=gene, label=label,
                                     vtype=vtype, source=source, expected=expected, depth=depth,
                                     ref_base="ref", ref_reads=depth - s["n_indel_reads"],
                                     alt_base=s["major_indel"] or "NA", alt_reads=s["n_indel_reads"],
                                     other_reads=0, vaf=round(vaf, 4), n_cells_mut=n_mut,
                                     n_cells_wt=n_wt, note=str(s["allele_counts"])[:80]))
        elif vtype in ("INS",):
            # NPM1c-style: a real event has a dominant insertion length (4bp). require_dominant
            # restricts INS cells to the dominant length, killing the highly-expressed-gene error smear.
            s = scan_insertion(bam, chrom, pos, window=30, min_len=3, max_len=60, min_mapq=min_mapq,
                               require_dominant=True, dom_tol=1, bulk=bulk)
            n_mut = n_wt = 0
            for cbc, d in s["cell"].items():
                if d["ins"] >= mut_min:
                    n_mut += 1
                    barcode_rows.append(dict(uid=uid, gene=gene, label=label, vtype=vtype,
                                             source=source, expected=expected, bam_barcode=cbc,
                                             rc_barcode=rc(cbc),
                                             cell_state=(states.get(cbc) if states else ""),
                                             mut_reads=d["ins"], wt_reads=d["clean"],
                                             alt_base=f"ins{s['dom_len']}bp", ref_base="ref"))
                elif d["clean"] >= wt_min:
                    n_wt += 1
            depth = s["depth"]
            vaf = s["n_ins_reads"] / depth if depth else 0.0
            summary_rows.append(dict(uid=uid, chrom=chrom, pos=pos, gene=gene, label=label,
                                     vtype=vtype, source=source, expected=expected, depth=depth,
                                     ref_base="ref", ref_reads=depth - s["n_ins_reads"],
                                     alt_base=(f"ins{s['dom_len']}bp" if s["dom_len"] else "NA"),
                                     alt_reads=s["n_ins_reads"], other_reads=0, vaf=round(vaf, 4),
                                     n_cells_mut=n_mut, n_cells_wt=n_wt,
                                     note=f"dom_len={s['dom_len']} dom_frac={s['dom_count']}/{s['n_ins_total']} lens={str(s['len_counts'])[:60]}"))
        elif vtype in ("ITD",):
            # FLT3-ITD: region encoded as pos=start, notes carries end=<e> and edges=<exon boundaries>.
            # Two mechanisms are tried and the stronger is reported:
            #   (1) scan_itd       -- ITD as an INSERTION op (long-read KINNEX; full-length reads)
            #   (2) scan_itd_softclip -- ITD as a SOFT-CLIP breakpoint pile-up (short-read STAR;
            #       verified on BEAT-AML where ins ops are absent), excluding exon-boundary splice clips.
            end = pos + 300
            edges = []
            for kv in v.get("notes", "").split(";"):
                if kv.startswith("end="):
                    end = int(kv.split("=")[1])
                elif kv.startswith("edges="):
                    edges = [int(x) for x in kv.split("=")[1].split(",") if x]
            s_ins = scan_itd(bam, chrom, pos, end, min_len=6, max_len=400, min_mapq=min_mapq, bulk=bulk)
            s_sc = scan_itd_softclip(bam, fa, chrom, pos, end, edges, min_mapq=max(1, min_mapq // 20),
                                     bulk=bulk)
            use_sc = s_sc["n_itd_reads"] > s_ins["n_itd_reads"]
            cellmap = s_sc["cell"] if use_sc else s_ins["cell"]
            alt = (f"ITD_softclip@{s_sc['bp']}({'dup' if s_sc['tandem_dup'] else 'ins'})" if use_sc
                   else (f"ITD{s_ins['dom_len']}bp" if s_ins["dom_len"] else "NA"))
            n_itd_reads = s_sc["n_itd_reads"] if use_sc else s_ins["n_itd_reads"]
            depth = s_sc["n_span"] if use_sc else s_ins["n_span"]
            n_mut = n_wt = 0
            for cbc, d in cellmap.items():
                if d["itd"] >= mut_min:
                    n_mut += 1
                    barcode_rows.append(dict(uid=uid, gene=gene, label=label, vtype=vtype,
                                             source=source, expected=expected, bam_barcode=cbc,
                                             rc_barcode=rc(cbc),
                                             cell_state=(states.get(cbc) if states else ""),
                                             mut_reads=d["itd"], wt_reads=d["clean"],
                                             alt_base=alt, ref_base="ref"))
                elif d["clean"] >= wt_min:
                    n_wt += 1
            vaf = n_itd_reads / depth if depth else 0.0
            summary_rows.append(dict(uid=uid, chrom=chrom, pos=pos, gene=gene, label=label,
                                     vtype=vtype, source=source, expected=expected, depth=depth,
                                     ref_base="ref", ref_reads=depth - n_itd_reads,
                                     alt_base=alt, alt_reads=n_itd_reads, other_reads=0,
                                     vaf=round(vaf, 4), n_cells_mut=n_mut, n_cells_wt=n_wt,
                                     note=f"insop={s_ins['n_itd_reads']} softclip={s_sc['n_itd_reads']}@{s_sc['bp']} dup={s_sc['tandem_dup']}"))
    bam.close()
    return summary_rows, barcode_rows


def write_tsv(rows, path, fieldnames=None):
    if not rows:
        return
    fieldnames = fieldnames or list(rows[0].keys())
    with open(path, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--uid", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--panel", required=True, action="append", help="panel TSV (repeatable)")
    ap.add_argument("--fasta", default="/Users/saljh8/Dropbox/Revio/Other/Variants/SNV/genome.fa")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cell-annot", default=None)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--bulk", action="store_true",
                    help="bulk BAM without cell barcodes (e.g. BEAT-AML STAR RNA-seq): count all reads as one pseudo-cell")
    args = ap.parse_args()
    panel = []
    for p in args.panel:
        panel += load_panel(p)
    fa = pysam.FastaFile(args.fasta) if os.path.exists(args.fasta) else None
    states = load_cell_states(args.cell_annot)
    sm, bc = scan_bam(args.uid, args.bam, panel, fa, args.out, states=states, min_mapq=args.min_mapq,
                      bulk=args.bulk)
    write_tsv(sm, os.path.join(args.out, f"{args.uid}_variant_readcounts.tsv"))
    write_tsv(bc, os.path.join(args.out, f"{args.uid}_mut_barcodes.tsv"))
    log.info("%s: %d panel loci scanned, %d MUT cell rows", args.uid, len(sm), len(bc))
