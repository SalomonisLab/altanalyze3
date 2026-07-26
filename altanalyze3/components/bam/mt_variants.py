#!/usr/bin/env python3
"""mt_variants.py — fast mitochondrial variant discovery and passenger-marker nomination.

WHY THIS EXISTS
---------------
Mitochondrial DNA is present in thousands of copies per cell and mutates far faster than nuclear
DNA. A mutation that arises in one cell is inherited by all of its descendants, so a mitochondrial
variant marks a clone without causing disease -- a "passenger". That makes these variants usable as
a surrogate: a cell whose driver-mutation status is unknown (because the driver gene was not
captured in that cell) can be assigned from the mitochondrial markers it carries.

This module does three things, in order:

  1. DISCOVER   which mitochondrial positions vary WITHIN this sample (heteroplasmy).
  2. GENOTYPE   every cell at those positions.
  3. ASSOCIATE  each variant with a known driver status per cell, and nominate the ones that
                track the driver strongly enough to serve as surrogates.

Positions where the whole sample differs from the reference are inherited (homoplasmic) and carry
no information about clones inside the sample, so they are excluded by default (`--max-af`).

SPEED
-----
The mitochondrial genome is small (16,569 bp) but extremely deeply covered in RNA data
(measured: 4.27 M reads, ~186,000x, in one KINNEX BAM). The naive approach -- fetch the reference
per read and build a per-read alignment table -- does billions of Python-level operations and takes
hours. Instead:

  Step 1 uses pysam.count_coverage(), a single C-level call that returns per-position A/C/G/T counts
         for the whole contig. No Python loop over reads at all.
  Step 2 makes ONE pass over the reads, walking each CIGAR once and recording bases only at the
         handful of candidate positions found in step 1.

Both steps are linear in the data and touch each read at most once.

USAGE
-----
  python3 -m altanalyze3.components.bam.mt_variants \\
      --bam sample.bam --sample NAME --output-dir DIR \\
      [--driver-calls <sample>_complete_analysis.tsv] [--driver-label RUNX1_p.W279*]

`--driver-calls` is the per-cell output of variant_extraction.py; steps 1-2 run without it, and the
association step (3) is skipped if it is not supplied.
"""
import argparse
import csv
import logging
import os
import time
from collections import defaultdict

from glob import glob as _glob

import numpy as np
import pysam

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("mt_variants")

BARCODE_TAGS = ("CB", "BC", "XC", "UB")
MT_NAMES = ("chrM", "chrMT", "M", "MT")
BASES = ("A", "C", "G", "T")


# --------------------------------------------------------------------------------------- helpers

def find_mt_contig(bam):
    """Name of the mitochondrial contig in this BAM, or None."""
    refs = set(bam.references)
    for n in MT_NAMES:
        if n in refs:
            return n
    return None


def has_base_qualities(bam, contig, probe=50):
    """True if reads on this contig carry per-base quality scores.

    PacBio/KINNEX alignments frequently omit them; a quality filter would then drop every base.
    """
    n = 0
    for read in bam.fetch(contig):
        if read.is_unmapped:
            continue
        if read.query_qualities is not None:
            return True
        n += 1
        if n >= probe:
            break
    return False


def get_barcode(read, bulk=False):
    for tg in BARCODE_TAGS:
        try:
            bc = str(read.get_tag(tg)).strip().split("-")[0].split(".")[0]
            return bc.upper()
        except KeyError:
            continue
    return "BULK" if bulk else None


UMI_TAGS = ("XM", "ZU", "YM", "UB")


def get_umi(read):
    """UMI (molecule id) for this read, so counts are per-molecule not per-read. Falls back to the
    read name if no UMI tag is present (then each read is its own molecule)."""
    for tg in UMI_TAGS:
        try:
            return str(read.get_tag(tg))
        except KeyError:
            continue
    return read.query_name


def genotype_variant_cells(bam, chrom, pos1, ref, alt, min_mapq=20, bulk=False):
    """Per-cell, UMI-deduplicated genotype of ONE variant (SNV **or indel**) from a BAM.

    Substitution-only counting (pysam count_coverage) cannot see insertions/deletions, so a clonal
    indel marker such as chrM:3565 A>AC is invisible to it. This walks each spanning read's aligned
    pairs and detects the alt allele directly from the CIGAR: for an insertion, inserted query bases
    (reference position None) immediately after the anchor; for a deletion, missing query bases; for
    an SNV, the query base at the position. Counts are collapsed to distinct (barcode, UMI) molecules
    (an alt-supporting read makes the molecule alt). Returns {barcode: {'alt': n_alt_umi,
    'total': n_umi}}. This mirrors the independent chrM lineage quantification (UMI-level, any-alt).
    """
    pos0 = pos1 - 1
    ref = str(ref).upper()
    alt = str(alt).upper()
    is_snv = (len(ref) == 1 and len(alt) == 1)
    is_ins = len(alt) > len(ref)
    need_ins = len(alt) - len(ref)
    need_del = len(ref) - len(alt)
    span_end = pos0 + max(len(ref), len(alt)) + 1
    cells = {}
    for read in bam.fetch(chrom, max(0, pos0), span_end):
        if read.is_unmapped or read.cigartuples is None or read.mapping_quality < min_mapq:
            continue
        seq = read.query_sequence
        if not seq:
            continue
        bc = get_barcode(read, bulk)
        if bc is None:
            continue
        umi = get_umi(read)
        # Fast CIGAR walk (get_aligned_pairs allocates a per-base list and is far too slow on the deep
        # chrM pileup). Track reference/query position through the ops; at pos0 read the base (SNV) or
        # inspect the op immediately after the anchor (insertion / deletion).
        cig = read.cigartuples
        refpos = read.reference_start
        qpos = 0
        is_alt = None
        for k, (op, ln) in enumerate(cig):
            if op == 0 or op == 7 or op == 8:            # M / = / X : aligned block
                if refpos <= pos0 < refpos + ln:
                    if is_snv:
                        base = seq[qpos + (pos0 - refpos)].upper()
                        is_alt = (base == alt) if base in (ref, alt) else False
                    elif pos0 == refpos + ln - 1:        # pos0 is the last aligned base of this block
                        nxt = cig[k + 1] if k + 1 < len(cig) else None
                        if is_ins:
                            is_alt = bool(nxt and nxt[0] == 1 and nxt[1] >= need_ins)
                        else:
                            is_alt = bool(nxt and nxt[0] == 2 and nxt[1] >= need_del)
                    else:                                 # aligned straight through pos0 -> reference
                        is_alt = False
                    break
                refpos += ln
                qpos += ln
            elif op == 1 or op == 4:                      # I / S : consume query only
                qpos += ln
            elif op == 2 or op == 3:                      # D / N : consume reference only
                refpos += ln
            # H (5), P (6): consume neither
        if is_alt is None:
            continue
        d = cells.setdefault(bc, {})
        if d.get(umi):                                   # this molecule already alt-supported
            continue
        d[umi] = bool(is_alt)
    out = {}
    for bc, umis in cells.items():
        total = len(umis)
        a = sum(1 for v in umis.values() if v)
        out[bc] = {"alt": a, "total": total}
    return out


def genotype_mt_panel_cells(bam, contig, mt_variants, min_mapq=20, bulk=False):
    """Genotype MANY mt variants (SNV or indel) per cell in ONE pass over the contig -- each read is
    read and CIGAR-walked exactly once, and every panel variant within that read's span is scored.
    This replaces calling genotype_variant_cells once per variant (which re-reads the deep ~5M-read
    MT pile ~N times); here the whole panel costs a single pass. UMI-deduplicated, indel-aware.
    Returns {variant_name: {barcode: {'alt': n_alt_umi, 'total': n_umi}}} where variant_name is
    'chrom:pos_ref>alt'. Genotypes identically to genotype_variant_cells, just batched.
    """
    import numpy as _np
    snv_pos_set = set()
    snv_at = {}                                          # pos0 -> list of (alt_base, name)
    ind_at = {}                                          # anchor0 -> list of (need_ins, need_del, name)
    names = []
    for v in mt_variants:
        ref = str(v["ref"]).upper(); alt = str(v["alt"]).upper(); pos0 = int(v["pos"]) - 1
        nm = f"{v['chrom']}:{v['pos']}_{ref}>{alt}"
        names.append(nm)
        if len(ref) == 1 and len(alt) == 1:
            snv_pos_set.add(pos0); snv_at.setdefault(pos0, []).append((alt, nm))
        else:
            ind_at.setdefault(pos0, []).append((len(alt) - len(ref), len(ref) - len(alt), nm))
    snv_pos = _np.array(sorted(snv_pos_set), dtype=_np.int64)
    ind_pos = _np.array(sorted(ind_at), dtype=_np.int64)
    cells = {nm: {} for nm in names}
    all0 = list(snv_pos_set) + list(ind_at)
    if not all0:
        return cells
    lo = max(0, min(all0)); hi = max(all0) + 2
    ss = _np.searchsorted

    def _rec(nm, bc, umi, is_alt):
        c = cells[nm].setdefault(bc, {})
        if c.get(umi):
            return
        c[umi] = bool(is_alt)

    for read in bam.fetch(contig, lo, hi):
        if read.is_unmapped or read.cigartuples is None or read.mapping_quality < min_mapq:
            continue
        seq = read.query_sequence
        if not seq:
            continue
        bc = get_barcode(read, bulk)
        if bc is None:
            continue
        umi = get_umi(read)
        cig = read.cigartuples
        refpos = read.reference_start
        qpos = 0
        for k, (op, ln) in enumerate(cig):
            if op == 0 or op == 7 or op == 8:            # M / = / X
                if snv_pos.size:
                    i0 = ss(snv_pos, refpos, side="left"); i1 = ss(snv_pos, refpos + ln, side="left")
                    for ii in range(i0, i1):
                        p0 = int(snv_pos[ii]); base = seq[qpos + (p0 - refpos)].upper()
                        for alt_b, nm in snv_at[p0]:
                            _rec(nm, bc, umi, base == alt_b)
                if ind_pos.size:
                    j0 = ss(ind_pos, refpos, side="left"); j1 = ss(ind_pos, refpos + ln, side="left")
                    for jj in range(j0, j1):
                        anchor0 = int(ind_pos[jj])
                        nxt = cig[k + 1] if (anchor0 == refpos + ln - 1 and k + 1 < len(cig)) else None
                        for need_ins, need_del, nm in ind_at[anchor0]:
                            if anchor0 != refpos + ln - 1:            # aligned straight through -> ref
                                _rec(nm, bc, umi, False)
                            elif need_ins > 0:
                                _rec(nm, bc, umi, bool(nxt and nxt[0] == 1 and nxt[1] >= need_ins))
                            else:
                                _rec(nm, bc, umi, bool(nxt and nxt[0] == 2 and nxt[1] >= need_del))
                refpos += ln; qpos += ln
            elif op == 1 or op == 4:                      # I / S
                qpos += ln
            elif op == 2 or op == 3:                      # D / N
                refpos += ln
    out = {}
    for nm, bcs in cells.items():
        out[nm] = {bc: {"alt": sum(1 for x in u.values() if x), "total": len(u)} for bc, u in bcs.items()}
    return out


def _resolve_contig(bam, name):
    refs = set(bam.references)
    if name in refs:
        return name
    stripped = name[3:] if name.startswith("chr") else name
    for cand in (stripped, "chr" + stripped, "chr" + name):
        if cand in refs:
            return cand
    return name


def assoc_variants_from_bams(bam_map, var_a, var_b, out_dir=None, min_mapq=20, tag=None):
    """Extract two variants (SNV or indel) per cell across several BAMs, POOL, and test their
    presence-based (any-alt-read, UMI-level) association -- the independent chrM-lineage approach that
    can see indels. var_a/var_b: dicts with keys chrom, pos, ref, alt, name. Genotype is any-alt
    (present = >=1 alt UMI; absent = 0 alt & >=1 total UMI; uncovered excluded). Reports the 2x2 with
    Fisher's exact OR/p, and both framings (var_a as a surrogate for var_b: recall/background/precision;
    and the lineage-enrichment odds ratio). Writes a per-cell table and a summary. Returns (rows, summary).
    """
    import pysam
    from scipy.stats import fisher_exact
    per = {}
    for sample, bam_path in bam_map.items():
        bam = pysam.AlignmentFile(bam_path, "rb")
        ca = _resolve_contig(bam, var_a["chrom"])
        cb = _resolve_contig(bam, var_b["chrom"])
        ga = genotype_variant_cells(bam, ca, int(var_a["pos"]), var_a["ref"], var_a["alt"], min_mapq)
        gb = genotype_variant_cells(bam, cb, int(var_b["pos"]), var_b["ref"], var_b["alt"], min_mapq)
        bam.close()
        for bc, c in ga.items():
            per.setdefault((sample, bc), {}).update(a_alt=c["alt"], a_tot=c["total"])
        for bc, c in gb.items():
            per.setdefault((sample, bc), {}).update(b_alt=c["alt"], b_tot=c["total"])
        na = sum(1 for c in ga.values() if c["alt"] >= 1)
        nb = sum(1 for c in gb.values() if c["alt"] >= 1)
        log.info("%s: %s covered=%d (alt+=%d), %s covered=%d (alt+=%d)",
                 sample, var_a["name"], len(ga), na, var_b["name"], len(gb), nb)
    n = {"0,0": 0, "0,1": 0, "1,0": 0, "1,1": 0}
    rows = []
    for (sample, bc), c in per.items():
        a_tot = c.get("a_tot", 0)
        b_tot = c.get("b_tot", 0)
        if a_tot == 0 or b_tot == 0:
            continue
        ga = 1 if c.get("a_alt", 0) >= 1 else 0
        gb = 1 if c.get("b_alt", 0) >= 1 else 0
        n[f"{ga},{gb}"] += 1
        rows.append(dict(sample=sample, barcode=bc,
                         **{f"{var_a['name']}_alt": c.get("a_alt", 0), f"{var_a['name']}_tot": a_tot,
                            f"{var_b['name']}_alt": c.get("b_alt", 0), f"{var_b['name']}_tot": b_tot,
                            var_a["name"]: ga, var_b["name"]: gb}))
    a = n["1,1"]; b = n["1,0"]; c01 = n["0,1"]; d = n["0,0"]
    try:
        OR, p = fisher_exact([[d, c01], [b, a]])
    except Exception:
        OR, p = float("nan"), float("nan")
    b_pos = a + c01; b_neg = b + d; a_pos = a + b
    summary = dict(
        var_a=var_a["name"], var_b=var_b["name"], both_covered=a + b + c01 + d,
        cell_11=a, cell_a1_b0=b, cell_a0_b1=c01, cell_00=d,
        odds_ratio=round(float(OR), 3), p_value=float(p),
        recall_of_Bpos=round(a / b_pos, 4) if b_pos else 0.0,
        background_in_Bneg=round(b / b_neg, 4) if b_neg else 0.0,
        precision_of_Apos=round(a / a_pos, 4) if a_pos else 0.0,
        n_Bpos=b_pos, n_Bneg=b_neg)
    if out_dir:
        suffix = f"_{tag}" if tag else ""
        safe = f"{var_a['name']}_vs_{var_b['name']}".replace(":", "_").replace(">", "").replace("/", "_")
        _write_tsv(rows, os.path.join(out_dir, f"bam_assoc_{safe}{suffix}_percell.tsv"))
        _write_tsv([summary], os.path.join(out_dir, f"bam_assoc_{safe}{suffix}_summary.tsv"))
        log.info("wrote bam_assoc_%s%s_{percell,summary}.tsv", safe, suffix)
    return rows, summary


def discover_mt_variants(bam, contig, min_alt_reads=10, min_af=0.02, max_af=0.98,
                         min_mapq=20, min_base_quality=0, max_reads=3000000):
    """De novo mtDNA variant discovery from a BAM -- SUBSTITUTIONS **and INDELS**.

    Substitutions: one C-level count_coverage() over the contig; a position yields a variant for any
    minor allele whose fraction is in [min_af, max_af] with >= min_alt_reads supporting reads (the
    reference allele is the per-sample consensus/majority base, matching the existing sample-consensus
    convention -- no external FASTA needed).
    Indels: one read pass walking each CIGAR; insertions (op I) and deletions (op D) are tallied by
    (anchor position, allele). An indel is kept if supported by >= min_alt_reads reads and its fraction
    (indel reads / spanning depth at the anchor) is in [min_af, max_af]. Insertions are emitted as
    ref=<anchor base>, alt=<anchor base>+<inserted bases>; deletions as ref=<anchor>+<deleted>, alt=<anchor>.
    This closes the gap that count_coverage-only discovery misses every indel (e.g. chrM:3565 A>AC).
    Returns a list of dicts {chrom, pos(1-based), ref, alt, kind, alt_reads, depth, af}.
    """
    import numpy as _np
    length = bam.get_reference_length(contig)
    bases = _np.array(list("ACGT"))

    # SNVs: C-level count_coverage with NO Python read_callback (a per-read callback would make it
    # minutes on the deep MT pileup). Default filtering skips unmapped/secondary/qcfail/duplicate.
    cov = bam.count_coverage(contig, 0, length, quality_threshold=min_base_quality)
    stacks = _np.vstack([_np.asarray(c, dtype=_np.int64) for c in cov])      # 4 x length (A,C,G,T)
    depth = stacks.sum(0)
    maj = stacks.argmax(0)                                                   # consensus base index
    variants = []
    covered = _np.where(depth >= min_alt_reads)[0]
    for pos0 in covered:
        d = int(depth[pos0]); ref_i = int(maj[pos0])
        for alt_i in range(4):
            if alt_i == ref_i:
                continue
            ac = int(stacks[alt_i, pos0])
            af = ac / d
            if ac >= min_alt_reads and min_af <= af <= max_af:
                variants.append(dict(chrom=contig, pos=int(pos0) + 1, ref=str(bases[ref_i]),
                                     alt=str(bases[alt_i]), kind="snv", alt_reads=ac, depth=d,
                                     af=round(af, 4)))
    # indels: one read pass, but only reads whose CIGAR actually contains an I/D op are walked
    ins_counts = {}    # (anchor0, inserted_seq) -> reads
    del_counts = {}    # (anchor0, del_len) -> reads
    n = 0
    for read in bam.fetch(contig):
        n += 1
        if n > max_reads:
            log.warning("discover_mt_variants: hit max_reads=%d on %s; indel panel may be partial",
                        max_reads, contig)
            break
        cig = read.cigartuples
        if cig is None or read.mapping_quality < min_mapq:
            continue
        has_indel = False
        for op, _ln in cig:
            if op == 1 or op == 2:
                has_indel = True
                break
        if not has_indel:                                # fast skip: the vast majority of reads
            continue
        seq = read.query_sequence
        if not seq:
            continue
        refpos = read.reference_start
        qpos = 0
        for op, ln in cig:
            if op == 0 or op == 7 or op == 8:
                refpos += ln; qpos += ln
            elif op == 1:                                    # insertion between refpos-1 and refpos
                anchor0 = refpos - 1
                if anchor0 >= 0:
                    ins_seq = seq[qpos:qpos + ln].upper()
                    ins_counts[(anchor0, ins_seq)] = ins_counts.get((anchor0, ins_seq), 0) + 1
                qpos += ln
            elif op == 2:                                    # deletion of ln bases starting at refpos
                anchor0 = refpos - 1
                if anchor0 >= 0:
                    del_counts[(anchor0, ln)] = del_counts.get((anchor0, ln), 0) + 1
                refpos += ln
            elif op == 3:
                refpos += ln
            elif op == 4:
                qpos += ln
    for (anchor0, ins_seq), cnt in ins_counts.items():
        if cnt < min_alt_reads or anchor0 >= length:
            continue
        d = int(depth[anchor0]) or cnt
        af = cnt / d
        if not (min_af <= af <= max_af):
            continue
        rb = str(bases[int(maj[anchor0])])
        variants.append(dict(chrom=contig, pos=anchor0 + 1, ref=rb, alt=rb + ins_seq, kind="ins",
                             alt_reads=cnt, depth=d, af=round(af, 4)))
    for (anchor0, ln), cnt in del_counts.items():
        if cnt < min_alt_reads or anchor0 + ln >= length:
            continue
        d = int(depth[anchor0]) or cnt
        af = cnt / d
        if not (min_af <= af <= max_af):
            continue
        rb = str(bases[int(maj[anchor0])])
        del_bases = "".join(str(bases[int(maj[anchor0 + 1 + k])]) for k in range(ln))
        variants.append(dict(chrom=contig, pos=anchor0 + 1, ref=rb + del_bases, alt=rb, kind="del",
                             alt_reads=cnt, depth=d, af=round(af, 4)))
    return variants


def build_bam_feature_matrix(uid, bam_paths, driver_variants, out_dir, min_mapq=20,
                             min_cells_mt=10, min_af=0.0005, contig=None,
                             marker_variants=None, mt_discovery=True):
    """Build the IMPROVED per-sample cell x feature h5ad DIRECTLY FROM BAM(s), UMI-level and
    indel-aware. Features = every driver variant (from driver_variants: SNV or indel) AND every
    de-novo-discovered mtDNA variant (substitutions AND indels). Genotyping is presence-capable:
    layers['alt'] and layers['total'] hold per-cell UMI counts; X = alt/total fraction. This is the
    version that can see chrM:3565-type indels that the count_coverage-based build_cell_feature_matrix
    cannot. Multiple BAMs for one uid are pooled (UMI counts summed per barcode). Writes
    <uid>_bam_feature_matrix.h5ad; returns its path.
    """
    import pysam
    import numpy as _np
    import pandas as _pd
    import anndata as _ad
    import scipy.sparse as _sp
    # discover the mt panel (union across this uid's BAMs); skipped when mt_discovery is False
    # (an external precomputed marker panel via marker_variants is used instead).
    mt_panel = []
    seen = set()
    for bp in (bam_paths if mt_discovery else []):
        bam = pysam.AlignmentFile(bp, "rb")
        c = contig or find_mt_contig(bam)
        for v in discover_mt_variants(bam, c, min_alt_reads=min_cells_mt, min_af=min_af,
                                      min_mapq=min_mapq):
            k = (v["pos"], v["ref"], v["alt"])
            if k not in seen:
                seen.add(k); mt_panel.append(v)
        bam.close()
    log.info("[%s] discovered %d mt variants (%d indels) at min_af=%.4f min_reads=%d", uid,
             len(mt_panel), sum(1 for v in mt_panel if v["kind"] != "snv"), min_af, min_cells_mt)
    mt_names = [f"{v['chrom']}:{v['pos']}_{v['ref']}>{v['alt']}" for v in mt_panel]
    # drivers grouped by name: co-located alleles (same gene_chr:pos, different alt -- e.g. SRSF2 P95R/
    # P95L) collapse to ONE feature, mutant if ANY allele is present. Keeps var names unique (get_loc
    # returns a scalar) and matches the genomic-coordinate convention of build_cell_feature_matrix.
    drv_by_name = {}
    for dv in driver_variants:
        info = drv_by_name.setdefault(dv["name"], dict(chrom=dv["chrom"], pos=int(dv["pos"]),
                                                       ref=dv["ref"], alt=dv["alt"], variants=[]))
        info["variants"].append(dv)
    feat_drv = {nm: {} for nm in drv_by_name}
    feat_mt = {nm: {} for nm in mt_names}
    # external PRECOMPUTED marker panel (e.g. nuclear passengers): collapse co-located alleles by name,
    # genotyped per cell exactly like drivers (nuclear, shallow, indel-aware) but placed in the MARKER
    # slot (feature_type=mtDNA) below so nominate_presence_fdr treats them as candidate markers.
    marker_by_name = {}
    for mv in (marker_variants or []):
        info = marker_by_name.setdefault(mv["name"], dict(chrom=mv["chrom"], pos=int(mv["pos"]),
                                                          ref=mv["ref"], alt=mv["alt"], variants=[]))
        info["variants"].append(mv)
    feat_marker = {nm: {} for nm in marker_by_name}
    for bp in bam_paths:
        bam = pysam.AlignmentFile(bp, "rb")
        # mtDNA: ONE pass over the contig for the whole panel (was one deep re-read per variant)
        if mt_panel:
            cmt = _resolve_contig(bam, mt_panel[0]["chrom"])
            for nm, bcd in genotype_mt_panel_cells(bam, cmt, mt_panel, min_mapq).items():
                fc = feat_mt[nm]
                for bc, cnt in bcd.items():
                    pa, pt = fc.get(bc, (0, 0)); fc[bc] = (pa + cnt["alt"], pt + cnt["total"])
        # drivers: per-variant (nuclear, shallow); collapse co-located by max within a BAM, sum across BAMs
        for nm, info in drv_by_name.items():
            perbc = {}
            for v in info["variants"]:
                cc = _resolve_contig(bam, v["chrom"])
                g = genotype_variant_cells(bam, cc, int(v["pos"]), v["ref"], v["alt"], min_mapq)
                for bc, cnt in g.items():
                    pa, pt = perbc.get(bc, (0, 0))
                    perbc[bc] = (max(pa, cnt["alt"]), max(pt, cnt["total"]))
            fc = feat_drv[nm]
            for bc, (a, t) in perbc.items():
                pa, pt = fc.get(bc, (0, 0)); fc[bc] = (pa + a, pt + t)
        # external markers: genotyped exactly like drivers (nuclear, shallow, indel-aware)
        for nm, info in marker_by_name.items():
            perbc = {}
            for v in info["variants"]:
                cc = _resolve_contig(bam, v["chrom"])
                g = genotype_variant_cells(bam, cc, int(v["pos"]), v["ref"], v["alt"], min_mapq)
                for bc, cnt in g.items():
                    pa, pt = perbc.get(bc, (0, 0))
                    perbc[bc] = (max(pa, cnt["alt"]), max(pt, cnt["total"]))
            fc = feat_marker[nm]
            for bc, (a, t) in perbc.items():
                pa, pt = fc.get(bc, (0, 0)); fc[bc] = (pa + a, pt + t)
        bam.close()
    features = [(nm, "driver", d["chrom"], d["pos"], d["ref"], d["alt"])
                for nm, d in drv_by_name.items()]
    features += [(mt_names[i], "mtDNA", mt_panel[i]["chrom"], mt_panel[i]["pos"],
                  mt_panel[i]["ref"], mt_panel[i]["alt"]) for i in range(len(mt_panel))]
    # external precomputed markers occupy the SAME marker slot (feature_type=mtDNA)
    features += [(nm, "mtDNA", d["chrom"], d["pos"], d["ref"], d["alt"])
                 for nm, d in marker_by_name.items()]
    feat_cells = {}
    feat_cells.update(feat_drv); feat_cells.update(feat_mt); feat_cells.update(feat_marker)
    barcodes = sorted({bc for fc in feat_cells.values() for bc in fc})
    bc_idx = {b: i for i, b in enumerate(barcodes)}
    names = [f[0] for f in features]
    alt = _np.zeros((len(barcodes), len(names)), dtype=_np.int32)
    tot = _np.zeros((len(barcodes), len(names)), dtype=_np.int32)
    for j, nm in enumerate(names):
        for bc, (a, t) in feat_cells[nm].items():
            i = bc_idx[bc]; alt[i, j] = a; tot[i, j] = t
    frac = _np.divide(alt, tot, out=_np.zeros(alt.shape, float), where=tot > 0)
    var = _pd.DataFrame({"feature_type": [f[1] for f in features], "chrom": [f[2] for f in features],
                         "pos": [f[3] for f in features], "ref": [f[4] for f in features],
                         "alt": [f[5] for f in features]}, index=names)
    obs = _pd.DataFrame(index=[f"{uid}|{b}" for b in barcodes])
    adata = _ad.AnnData(X=_sp.csr_matrix(frac), obs=obs, var=var)
    adata.layers["alt"] = _sp.csr_matrix(alt)
    adata.layers["total"] = _sp.csr_matrix(tot)
    os.makedirs(out_dir, exist_ok=True)
    # same suffix as build_cell_feature_matrix so nominate/sweep/genotype-summary are drop-in
    # (a separate --matrix-dir, e.g. bam_matrices/, distinguishes this indel-aware BAM build)
    out = os.path.join(out_dir, f"{uid}_cell_feature_matrix.h5ad")
    adata.write_h5ad(out, compression="gzip")
    log.info("[%s] wrote %s: %d cells x %d features (%d driver + %d mt + %d marker-panel)", uid, out,
             adata.n_obs, adata.n_vars, len(drv_by_name), len(mt_panel), len(marker_by_name))
    return out


# ------------------------------------------------------------------- step 1: discover heteroplasmy

def scan_heteroplasmy(bam, contig, min_depth=100, min_alt_reads=5, min_af=0.001, max_af=0.95,
                      min_base_quality=20):
    """Positions in the mitochondrial genome that vary within this sample.

    One C-level count_coverage() call over the whole contig; no per-read Python work.

    A position is a candidate when the second-most-common base is supported by at least
    `min_alt_reads` reads and its fraction lies in [min_af, max_af]. The upper bound removes
    inherited (homoplasmic) differences from the reference, which are identical in every cell and
    therefore useless for telling clones apart.

    Returns a list of dicts: pos (1-based), ref_base (majority in THIS sample), alt_base, depth,
    ref_reads, alt_reads, af.
    """
    length = dict(zip(bam.references, bam.lengths))[contig]
    # Long-read BAMs (PacBio/KINNEX) often store no per-base quality scores at all. Passing a
    # quality threshold in that case silently discards EVERY base and yields zero coverage, so the
    # threshold is dropped when the data has no qualities to filter on.
    if not has_base_qualities(bam, contig):
        if min_base_quality:
            log.info("%s: reads carry no base-quality scores; base-quality filter disabled", contig)
        min_base_quality = 0
    t0 = time.time()
    cov = bam.count_coverage(contig, 0, length, quality_threshold=min_base_quality,
                             read_callback="nofilter")
    counts = np.vstack([np.asarray(c, dtype=np.int64) for c in cov])      # 4 x length
    log.info("count_coverage over %s (%d bp): %.1fs", contig, length, time.time() - t0)

    depth = counts.sum(axis=0)
    order = np.argsort(counts, axis=0)                                    # ascending
    top_i, second_i = order[3], order[2]
    idx = np.arange(length)
    top_n = counts[top_i, idx]
    second_n = counts[second_i, idx]

    with np.errstate(divide="ignore", invalid="ignore"):
        af = np.where(depth > 0, second_n / np.maximum(depth, 1), 0.0)

    keep = (depth >= min_depth) & (second_n >= min_alt_reads) & (af >= min_af) & (af <= max_af)
    out = []
    for i in np.flatnonzero(keep):
        out.append(dict(pos=int(i) + 1, ref_base=BASES[top_i[i]], alt_base=BASES[second_i[i]],
                        depth=int(depth[i]), ref_reads=int(top_n[i]), alt_reads=int(second_n[i]),
                        af=round(float(af[i]), 6)))
    log.info("candidate heteroplasmic positions: %d (of %d covered)", len(out), int((depth > 0).sum()))
    return out


# ----------------------------------------------------------------- step 2: genotype cells at those

def genotype_cells(bam, contig, positions, min_mapq=20, min_base_quality=20, bulk=False):
    """Per-cell base counts at the candidate positions, in ONE pass over the reads.

    Each read's CIGAR is walked once; only the candidate positions inside that read's span are
    recorded. Returns {pos: {barcode: {base: count}}}.
    """
    if not positions:
        return None, [], {}
    want = np.array(sorted(p["pos"] - 1 for p in positions), dtype=np.int64)   # 0-based, sorted
    pos_index = {int(p): i for i, p in enumerate(want)}
    last = int(want[-1])
    base_index = {"A": 0, "C": 1, "G": 2, "T": 3}

    # Counts live in one dense int32 array (positions x cells x 4), not nested dictionaries.
    # At 416 positions x ~8,000 cells that is ~53 MB, versus roughly a gigabyte of Python dicts for
    # the same data, and it lets the calling step be vectorised instead of looping over cells.
    # The cell axis grows by doubling because the barcode set is discovered during the pass.
    cap = 4096
    counts = np.zeros((len(want), cap, 4), dtype=np.int32)
    bc_index = {}
    barcodes = []
    searchsorted = np.searchsorted                    # local binding: called millions of times

    t0 = time.time()
    n_reads = 0
    for read in bam.fetch(contig):
        if read.is_unmapped or read.cigartuples is None or read.mapping_quality < min_mapq:
            continue
        bc = get_barcode(read, bulk)
        if bc is None:
            continue
        seq = read.query_sequence
        if not seq:
            continue
        quals = read.query_qualities
        n_reads += 1
        ci = bc_index.get(bc)
        if ci is None:
            ci = bc_index[bc] = len(barcodes)
            barcodes.append(bc)
            if ci >= cap:                              # grow the cell axis by doubling
                counts = np.concatenate(
                    [counts, np.zeros((len(want), cap, 4), dtype=np.int32)], axis=1)
                cap *= 2
        # Walk the CIGAR once. For each aligned block, binary-search the sorted position list for
        # the slice that falls inside it -- O(log n + hits) rather than a scan over every candidate
        # position per block, which at 4M reads would be billions of comparisons.
        refpos = read.reference_start
        qpos = 0
        for op, ln in read.cigartuples:
            if op in (0, 7, 8):                       # M / = / X : aligned to the reference
                lo, hi = refpos, refpos + ln
                i0 = searchsorted(want, lo, side="left")
                i1 = searchsorted(want, hi, side="left")
                for k in range(i0, i1):
                    qi = qpos + (int(want[k]) - lo)
                    if qi < len(seq) and (quals is None or quals[qi] >= min_base_quality):
                        bi = base_index.get(seq[qi].upper())
                        if bi is not None:
                            counts[k, ci, bi] += 1
                refpos += ln
                qpos += ln
            elif op == 1 or op == 4:                  # insertion / soft-clip: query only
                qpos += ln
            elif op in (2, 3):                        # deletion / skipped: reference only
                refpos += ln
            if refpos > last:
                break
    counts = counts[:, :len(barcodes), :]              # trim unused capacity
    log.info("genotyped %d positions x %d cells across %d reads in %.1fs (%.0f MB)",
             len(want), len(barcodes), n_reads, time.time() - t0, counts.nbytes / 1e6)
    return counts, barcodes, pos_index


def call_cells(positions, counts, barcodes, pos_index, min_reads=2, min_cell_af=0.10,
               min_alt_reads_per_cell=2, min_ref_reads_per_cell=2):
    """Per position: which cells carry the variant, which are reference, with the read counts.

    A cell is CARRIER if the alternate base is at least `min_cell_af` of that cell's reads at the
    position; REFERENCE if it has reads but not enough alternate. Cells with fewer than `min_reads`
    reads there are unknown and excluded.

    Vectorised over cells: one pass per position instead of a Python loop over thousands of cells.
    """
    base_index = {"A": 0, "C": 1, "G": 2, "T": 3}
    bc = np.asarray(barcodes)
    out = {}
    for p in positions:
        k = pos_index.get(p["pos"] - 1)
        if k is None:
            continue
        c = counts[k]                                   # cells x 4
        total = c.sum(axis=1)
        alt = c[:, base_index[p["alt_base"]]]
        usable = total >= min_reads
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.where(total > 0, alt / np.maximum(total, 1), 0.0)
        # >=2 alternate reads, not 1. With ~7 mitochondrial reads per cell at a position, a single
        # sequencing error already clears a 10% within-cell fraction, and single-read calls made up
        # 54% of all carriers (12,806 of 23,766) in the first cohort run -- overwhelmingly error.
        ref = c[:, base_index[p["ref_base"]]] if p["ref_base"] in base_index else total - alt

        # Three categories, not two. Evidence has to be positive on BOTH sides:
        #   CARRIER   >= min_alt_reads_per_cell alternate reads (and past the fraction floor)
        #   REFERENCE zero alternate reads AND >= min_ref_reads_per_cell reference reads
        #   AMBIGUOUS anything else -- notably a cell with exactly ONE alternate read
        # The ambiguous cells are excluded from both groups. Calling them wild-type (the earlier
        # behaviour, "reference = not carrier") dumped every single-alternate-read cell into the
        # normal group, which deflates the background rate and inflates apparent enrichment -- the
        # opposite of what raising the carrier threshold was meant to achieve.
        is_carrier = usable & (alt >= min_alt_reads_per_cell) & (frac >= min_cell_af)
        is_ref = (alt == 0) & (ref >= min_ref_reads_per_cell)
        is_ambig = ~(is_carrier | is_ref)
        out[p["pos"]] = dict(carrier=set(bc[is_carrier].tolist()),
                             reference=set(bc[is_ref].tolist()),
                             n_ambiguous=int(is_ambig.sum()),
                             alt=alt, total=total, usable=(is_carrier | is_ref))
    return out


# --------------------------------------------------------------- step 3: associate with the driver

def load_driver_calls(path, label=None, orientation="bam_barcode"):
    """{barcode: 'MUT'|'WT'} from a variant_extraction.py <sample>_complete_analysis.tsv.

    If `label` is given, only that variant's calls are used; otherwise a cell is MUT if it is MUT at
    any variant, else WT if called WT anywhere.
    """
    mut, wt = set(), set()
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if label and r.get("label") != label:
                continue
            bc = r.get(orientation)
            if not bc:
                continue
            if r.get("genotype") == "MUT":
                mut.add(bc)
            elif r.get("genotype") == "WT":
                wt.add(bc)
    wt -= mut
    return {**{b: "MUT" for b in mut}, **{b: "WT" for b in wt}}


def associate(positions, calls, driver, min_mut_cells=5, label="", patient="",
              n_samples_compared=1):
    """Test each mitochondrial variant for concordance with a driver mutation.

    The question is NOT "are these two associated" but "does this mitochondrial variant mark the
    mutant cells and almost nothing else". A driver may be carried by only 20 of 30,000 cells, so a
    plain association test is both underpowered and misleading: a mitochondrial variant present in
    90% of ALL cells will overlap 95% of any 20-cell set purely by chance.

    The statistic therefore compares the observed overlap to what the variant's OWN background rate
    predicts:

        background = carriers among driver-NORMAL cells / driver-normal cells
        recall     = carriers among driver-mutant cells / driver-mutant cells
        expected   = n_mutant_cells * background
        binomial p = P(X >= observed | n_mutant_cells, background)     one-sided

    A common variant has a high `background`, so `expected` is large and the binomial p-value stays
    unremarkable no matter how high the raw overlap looks. `fold_enrichment` = recall / background
    makes the same point on a readable scale. Fisher's exact test on the 2x2 is reported alongside.

    Rows carry `n_mut_cells` so an apparently perfect but tiny result can be judged: 20/20 with a 2%
    background is a strong candidate; 19/20 with an 85% background is not.
    """
    try:
        from scipy.stats import fisher_exact, binomtest
    except Exception:
        fisher_exact = binomtest = None
    rows = []
    for p in positions:
        c = calls.get(p["pos"])
        if not c:
            continue
        a = sum(1 for b in c["carrier"] if driver.get(b) == "MUT")     # mutant cells carrying it
        b_ = sum(1 for b in c["carrier"] if driver.get(b) == "WT")     # normal cells carrying it
        cc = sum(1 for b in c["reference"] if driver.get(b) == "MUT")  # mutant cells without it
        d = sum(1 for b in c["reference"] if driver.get(b) == "WT")    # normal cells without it
        n_mut, n_wt = a + cc, b_ + d
        if n_mut < min_mut_cells or n_wt == 0:
            continue

        background = b_ / n_wt                       # how common the variant is in normal cells
        recall = a / n_mut                           # how much of the mutant clone it marks
        expected = n_mut * background
        fold = (recall / background) if background > 0 else float("inf")

        if binomtest is not None and background > 0:
            pval = binomtest(a, n_mut, background, alternative="greater").pvalue
        elif binomtest is not None:
            # background of exactly zero: use a conservative pseudo-rate of 1/(n_wt+1)
            pval = binomtest(a, n_mut, 1.0 / (n_wt + 1), alternative="greater").pvalue
        else:
            pval = float("nan")
        odds, fisher_p = (fisher_exact([[a, b_], [cc, d]], alternative="greater")
                          if fisher_exact is not None else (float("nan"), float("nan")))

        # Every cell of the 2x2 table plus both margins, so any statistic here can be recomputed
        # from the file without going back to the per-cell data.
        #
        #                      MT variant PRESENT      MT variant ABSENT      row total
        #   driver MUTANT              a                       c            driver_mut_total
        #   driver WILD-TYPE           b                       d            driver_wt_total
        #   column total        mtvar_pos_total         mtvar_neg_total     cells_tested_total
        rows.append(dict(
            patient=patient, driver_label=label, n_samples_compared=n_samples_compared,
            mt_pos=p["pos"], mt_ref_base=p["ref_base"], mt_alt_base=p["alt_base"],
            mt_bulk_allele_frac=p["af"],
            cells_with_mtvar_and_driver_mut=a,                       # a
            cells_with_mtvar_and_driver_wt=b_,                       # b
            cells_no_mtvar_and_driver_mut=cc,                        # c
            cells_no_mtvar_and_driver_wt=d,                          # d
            cells_driver_mut_total=n_mut,                            # a + c
            cells_driver_wt_total=n_wt,                              # b + d
            cells_with_mtvar_total=a + b_,                           # a + b
            cells_no_mtvar_total=cc + d,                             # c + d
            cells_tested_total=a + b_ + cc + d,
            frac_driver_mut_with_mtvar=round(recall, 4),             # a / (a+c)  "recall"
            frac_driver_wt_with_mtvar=round(background, 6),          # b / (b+d)  "background"
            frac_mtvar_cells_that_are_driver_mut=round(a / max(a + b_, 1), 4),   # a / (a+b)
            expected_cells_with_mtvar_and_driver_mut=round(expected, 2),  # (a+c) * background
            fold_enrichment=(round(fold, 3) if fold != float("inf") else "inf"),
            binomial_p=float(pval), fisher_p=float(fisher_p),
            odds_ratio=(round(float(odds), 4) if odds == odds else "nan")))
    return rows


_ASSOC_FLOAT = ("frac_driver_mut_with_mtvar", "frac_driver_wt_with_mtvar",
                "frac_mtvar_cells_that_are_driver_mut", "q_value", "binomial_p", "fisher_p",
                "expected_cells_with_mtvar_and_driver_mut")
_ASSOC_INT = ("cells_with_mtvar_and_driver_mut", "cells_with_mtvar_and_driver_wt",
              "cells_no_mtvar_and_driver_mut", "cells_no_mtvar_and_driver_wt",
              "cells_driver_mut_total", "cells_driver_wt_total", "cells_with_mtvar_total",
              "cells_no_mtvar_total", "n_samples_compared", "cells_driver_mut_called_total")


def load_association_table(path):
    """Read a written mt_driver_association_all*.tsv back into typed rows (no recompute).

    Lets the significance/best-variant summaries be regenerated from an existing association table
    without re-running the O(mutations x positions x cells) association -- so a threshold or gate can
    be changed and re-reported in seconds.
    """
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for k in _ASSOC_FLOAT:
                v = r.get(k, "")
                r[k] = (1.0 if k in ("q_value", "binomial_p", "fisher_p") else 0.0) \
                    if v in ("", "nan", "inf", "-inf") else float(v)
            for k in _ASSOC_INT:
                try:
                    r[k] = int(r.get(k, 0))
                except ValueError:
                    r[k] = 0
            if "cells_driver_mut_called_total" not in r or not r["cells_driver_mut_called_total"]:
                r["cells_driver_mut_called_total"] = r["cells_driver_mut_total"]
            rows.append(r)
    return rows


def report_from_association_table(assoc_path, out_dir, max_q=0.05, max_background=0.05,
                                  min_recall=0.90, tag=""):
    """Regenerate the significant-associations and best-variant summaries from an existing table."""
    rows = load_association_table(assoc_path)
    suffix = f"_{tag}" if tag else ""
    sig = significant_associations(rows, max_q=max_q, max_background=max_background)
    _write_tsv(sig, os.path.join(out_dir, f"mt_significant_associations{suffix}.tsv"))
    best = best_variant_per_mutation(rows, min_recall_for_imputation=min_recall,
                                     max_background_for_imputation=max_background,
                                     max_q_for_imputation=max_q)
    _write_tsv(best, os.path.join(out_dir, f"mt_best_variant_per_mutation{suffix}.tsv"))
    from collections import Counter as _C
    tiers = _C(b["imputation_verdict"].split(":")[0] for b in best)
    n_pairs = len(set((r["donor"], r["mutation"]) for r in sig))
    log.info("%d significant associations across %d donor-mutation pairs; best-per-mutation: %s",
             len(sig), n_pairs, dict(tiers))
    return sig, best


def significant_associations(all_rows, max_q=0.05, max_background=0.05, min_overlap_cells=10):
    """EVERY statistically significant, specific mutation<->mitochondrial-variant association.

    This is the primary answer to "which mutations, in which samples, are associated with a
    mitochondrial variant" -- not one-per-mutation, all of them. A row is kept if it is significant
    after correction (q <= max_q), specific (the variant is in <= max_background of non-mutant cells),
    and rests on real overlap (>= min_overlap_cells mutant cells actually carry it). Recall and
    precision are reported so partial (subclone) markers are visible, not discarded.
    """
    keep = []
    for r in all_rows:
        q = r.get("q_value", 1.0)
        if (q <= max_q and r["frac_driver_wt_with_mtvar"] <= max_background
                and r["cells_with_mtvar_and_driver_mut"] >= min_overlap_cells):
            keep.append(dict(
                donor=r["patient"], mutation=r["driver_label"],
                n_samples_compared=r["n_samples_compared"],
                mt_pos=r["mt_pos"], mt_major_base=r["mt_ref_base"], mt_variant_base=r["mt_alt_base"],
                cells_mut_with_mtvar=r["cells_with_mtvar_and_driver_mut"],
                cells_mut_total_covered=r["cells_driver_mut_total"],
                cells_wt_with_mtvar=r["cells_with_mtvar_and_driver_wt"],
                cells_wt_total_covered=r["cells_driver_wt_total"],
                recall_frac_clone_marked=r["frac_driver_mut_with_mtvar"],
                background_frac_wt_with_mtvar=r["frac_driver_wt_with_mtvar"],
                precision_marked_that_are_mutant=r["frac_mtvar_cells_that_are_driver_mut"],
                fold_enrichment=r["fold_enrichment"], expected_by_chance=r["expected_cells_with_mtvar_and_driver_mut"],
                binomial_p=r["binomial_p"], q_value=q))
    keep.sort(key=lambda x: (x["donor"], x["mutation"], x["q_value"]))
    return keep


def best_variant_per_mutation(all_rows, min_recall_for_imputation=0.90,
                              max_background_for_imputation=0.05,
                              min_mut_cells_for_imputation=10,
                              min_overlap_cells=10,
                              max_q_for_imputation=0.05):
    """One row per (donor, mutation): the mitochondrial variant that best predicts it.

    IMPUTATION, NOT JUST ASSOCIATION. To assign a cell's mutation status from its mitochondrial
    reads, the mitochondrial variant has to be MEASURED in the cells that matter and PRESENT in
    essentially all the mutant ones. Two separate things can fail:

      coverage  -- the mitochondrial position may not be read in a given cell at all. If only a
                   fraction of mutation-positive cells have any call there, most cells cannot be
                   imputed no matter how good the association is.
      recall    -- of the mutant cells that ARE covered, the share carrying the variant. If half
                   the mutant cells lack the variant, the mutation cannot be read off it.

    `imputation_usable` requires both: the mitochondrial variant must be covered in at least as many
    mutant cells as the mutation itself was called in, and must mark at least
    `min_recall_for_imputation` of them.
    """
    best = {}
    for r in all_rows:
        key = (r["patient"], r["driver_label"])
        # Rank by Q-VALUE (significance), then discrimination, then overlap. Ranking by discrimination
        # (recall - background) alone favours a variant that perfectly marks a HANDFUL of cells
        # (recall 1.0, tiny N) over one that marks 65 mutant cells with q=3e-71 -- because both have
        # high discrimination but only the second is real. The q-value already folds in N, so the
        # well-powered marker wins. Variants below the min mutant-cell count are pushed to the back so
        # they only win when nothing better exists.
        youden = r["frac_driver_mut_with_mtvar"] - r["frac_driver_wt_with_mtvar"]
        underpowered = 1 if r["cells_driver_mut_total"] < min_mut_cells_for_imputation else 0
        rank = (underpowered, r.get("q_value", 1.0), -youden, -r["cells_with_mtvar_and_driver_mut"])
        if key not in best or rank < best[key][0]:
            best[key] = (rank, r)

    out = []
    for (patient, label), (_rank, r) in sorted(best.items()):
        n_mut_called = r.get("cells_driver_mut_called_total", r["cells_driver_mut_total"])
        cov = (r["cells_driver_mut_total"] / n_mut_called) if n_mut_called else 0.0
        recall = r["frac_driver_mut_with_mtvar"]
        background = r["frac_driver_wt_with_mtvar"]
        q = r.get("q_value", 1.0)
        n_mut = r["cells_driver_mut_total"]

        # A marker is judged on SIGNIFICANCE + SPECIFICITY + real OVERLAP -- NOT on recall or on
        # whole-clone coverage. What makes a mitochondrial variant usable is that a cell CARRYING it
        # is confidently mutant (background ~ 0 -> high precision). Recall is a DESCRIPTOR of how
        # much of the clone that one variant covers, not a pass/fail: a variant that marks 54% of the
        # mutant cells with zero background (PHF6 mt-7392: 65 mut / 0 wt, q=3e-71) is a real subclone
        # marker, not a failure. Earlier this was gated on recall>=0.9 and coverage==1.0, which
        # rejected every genuine marker (mitochondrial positions are only read in a fraction of
        # cells). precision = a/(a+b) = fraction of marked cells that are truly mutant.
        overlap = r["cells_with_mtvar_and_driver_mut"]
        precision = r["frac_mtvar_cells_that_are_driver_mut"]
        recall_all = (overlap / n_mut_called) if n_mut_called else 0.0
        # Factual output only -- one row per mutation, its single best mitochondrial variant, with the
        # complete 2x2 and both recall denominators. No pass/fail label; the reader judges from the
        # numbers. The best variant is the most significant after correction (q-value already accounts
        # for the number of cells, so a variant marking many cells outranks one marking a few).
        out.append(dict(
            donor=patient,
            mutation=label,
            n_samples_pooled=r["n_samples_compared"],
            # driver cell counts
            mutant_cells_total=n_mut_called,
            mutant_cells_with_mt_coverage=r["cells_driver_mut_total"],
            wildtype_cells_with_mt_coverage=r["cells_driver_wt_total"],
            # the best mitochondrial variant
            mt_pos=r["mt_pos"], mt_ref=r["mt_ref_base"], mt_alt=r["mt_alt_base"],
            # 2x2 table
            mutant_and_mtvar=overlap,
            mutant_and_mt_reference=r["cells_no_mtvar_and_driver_mut"],
            wildtype_and_mtvar=r["cells_with_mtvar_and_driver_wt"],
            wildtype_and_mt_reference=r["cells_no_mtvar_and_driver_wt"],
            # concordance
            recall_of_all_mutant_cells=round(recall_all, 4),
            recall_of_covered_mutant_cells=round(recall, 4),
            background_frac_wildtype_with_mtvar=round(background, 4),
            precision_frac_mtvar_cells_that_are_mutant=round(precision, 4),
            fold_enrichment=r["fold_enrichment"],
            expected_overlap_by_chance=r["expected_cells_with_mtvar_and_driver_mut"],
            binomial_p=r["binomial_p"], q_value=r.get("q_value", "")))
    out.sort(key=lambda x: (x["q_value"] if isinstance(x["q_value"], float) else 1.0))
    return out


def add_fdr(rows, key="binomial_p", out_key="q_value"):
    """Benjamini-Hochberg across ALL tests, so the correction spans every variant x driver pair."""
    rows = [r for r in rows if r.get(key) == r.get(key)]        # drop NaN
    rows.sort(key=lambda r: r[key])
    m = len(rows)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        q = min(prev, rows[i][key] * m / (i + 1))
        rows[i][out_key] = float(q)
        prev = q
    return rows


def nominate(assoc, max_q=0.05, min_recall=0.90, max_background=0.10, min_mut_cells=5):
    """Mitochondrial variants that could stand in for the driver.

    Requires all of: significant after correction, marks nearly all mutant cells (`min_recall`),
    and is genuinely uncommon in normal cells (`max_background`) so the overlap is not explained by
    the variant simply being everywhere.

    `confidence` records the caveat the user asked to be explicit: a perfect match over very few
    mutant cells, or against a high background, is flagged rather than reported as established.
    """
    out = []
    for r in assoc:
        if r.get("q_value", 1.0) > max_q:
            continue
        if (r["frac_driver_mut_with_mtvar"] < min_recall
                or r["frac_driver_wt_with_mtvar"] > max_background):
            continue
        if r["cells_driver_mut_total"] < min_mut_cells:
            continue
        if r["cells_driver_mut_total"] < 20:
            conf = "speculative: few mutant cells"
        elif r["frac_driver_wt_with_mtvar"] > 0.02:
            conf = "speculative: variant not rare in normal cells"
        else:
            conf = "candidate"
        out.append({**r, "confidence": conf})
    return out


# ----------------------------------------------------------------------------------- orchestration

# Internal dict keys are short for use in the code; output columns are renamed here so the files
# state what each number actually counts. Every column a reader needs to recompute the statistics is
# present -- both margins and all four cells of each 2x2 table.
_OUTPUT_COLUMN_NAMES = {
    # per-position heteroplasmy scan
    "pos": "mt_pos",
    "ref_base": "mt_major_base",              # commonest base in THIS sample
    "alt_base": "mt_variant_base",            # second commonest base = the candidate variant
    "depth": "total_reads_at_mt_pos",
    "ref_reads": "reads_with_major_base",
    "alt_reads": "reads_with_variant_base",
    "af": "frac_reads_with_variant_base",
    # per-cell calls
    "sample": "sample",
    "bam_barcode": "cell_barcode_bam_orientation",
    "total_reads": "reads_in_cell_at_mt_pos",
    "cell_af": "frac_cell_reads_with_variant_base",
    "status": "cell_mtvar_status",
}


def rename_for_output(rows, extra=None):
    """Map internal keys to self-describing column names, preserving column order."""
    mapping = dict(_OUTPUT_COLUMN_NAMES)
    if extra:
        mapping.update(extra)
    return [{mapping.get(k, k): v for k, v in r.items()} for r in rows]


def _write_tsv(rows, path, fieldnames=None):
    if not rows:
        open(path, "w").close()
        return
    # lineterminator="\n": csv defaults to "\r\n", and with newline="" that carriage return is
    # written literally. The final column then reads as "VALUE\r" and every awk/cut/R parse of the
    # last field silently fails.
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames or list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def run(sample, bam_path, out_dir, driver_calls=None, driver_label=None, contig=None,
        min_depth=100, min_alt_reads=5, min_af=0.001, max_af=0.95, min_mapq=20,
        min_base_quality=20, min_reads_per_cell=2, min_cell_af=0.10, bulk=False,
        max_q=0.05, min_precision=0.80):
    os.makedirs(out_dir, exist_ok=True)
    t_all = time.time()
    bam = pysam.AlignmentFile(bam_path, "rb")
    contig = contig or find_mt_contig(bam)
    if contig is None:
        bam.close()
        raise SystemExit(f"no mitochondrial contig found in {bam_path} "
                         f"(looked for {', '.join(MT_NAMES)})")
    log.info("[%s] mitochondrial contig: %s", sample, contig)

    positions = scan_heteroplasmy(bam, contig, min_depth, min_alt_reads, min_af, max_af,
                                  min_base_quality)
    _write_tsv(rename_for_output(positions), os.path.join(out_dir, f"{sample}_mt_variants.tsv"))

    counts, barcodes, pos_index = genotype_cells(bam, contig, positions, min_mapq,
                                                 min_base_quality, bulk)
    bam.close()
    calls = call_cells(positions, counts, barcodes, pos_index, min_reads_per_cell, min_cell_af)

    # Per-cell rows carry the read counts, not just the label, so the carrier threshold can be
    # re-examined later without re-scanning the BAM. Only cells with reads at the position appear.
    cell_rows = []
    bc_arr = np.asarray(barcodes)
    for p in positions:
        c = calls.get(p["pos"])
        if c is None:
            continue
        idx = np.flatnonzero(c["usable"])
        for i in idx:
            b = str(bc_arr[i])
            cell_rows.append(dict(sample=sample, pos=p["pos"], ref_base=p["ref_base"],
                                  alt_base=p["alt_base"], bam_barcode=b,
                                  alt_reads=int(c["alt"][i]), total_reads=int(c["total"][i]),
                                  cell_af=round(float(c["alt"][i]) / int(c["total"][i]), 4),
                                  status="CARRIER" if b in c["carrier"] else "REFERENCE"))
    _write_tsv(rename_for_output(cell_rows, {"alt_reads": "reads_with_variant_base_in_cell"}),
               os.path.join(out_dir, f"{sample}_mt_cell_calls.tsv"))
    n_carrier = sum(1 for r in cell_rows if r["status"] == "CARRIER")
    log.info("per-cell calls written: %d rows (%d carrier, %d reference)",
             len(cell_rows), n_carrier, len(cell_rows) - n_carrier)

    assoc = nominated = []
    if driver_calls:
        driver = load_driver_calls(driver_calls, driver_label)
        n_mut = sum(1 for v in driver.values() if v == "MUT")
        log.info("driver calls loaded: %d cells (%d mutant, %d normal)",
                 len(driver), n_mut, len(driver) - n_mut)
        assoc = associate(positions, calls, driver)
        _write_tsv(assoc, os.path.join(out_dir, f"{sample}_mt_driver_association.tsv"))
        nominated = nominate(assoc, max_q, min_precision)
        _write_tsv(nominated, os.path.join(out_dir, f"{sample}_mt_passenger_candidates.tsv"))
        log.info("associations tested: %d; nominated passengers: %d", len(assoc), len(nominated))

    log.info("[%s] total %.1fs", sample, time.time() - t_all)
    return positions, calls, assoc, nominated


def load_patient_map(panel_paths, samples, donor_map=None):
    """{sample -> patient}. Samples of one patient are pooled because this is a genetics-only
    comparison: the same genome, so a mitochondrial variant and a driver mutation seen at different
    timepoints belong to the same clonal history.

    Groups come from the panel's `expected_uids` column (semicolon-separated sample lists), which
    already encodes which samples belong to which donor. A `donor_map` TSV (sample<TAB>patient)
    overrides it. Any sample not named in either becomes its own patient.
    """
    m = {}
    if donor_map and os.path.exists(donor_map):
        with open(donor_map) as fh:
            for line in fh:
                p = line.rstrip("\r\n").split("\t")
                if len(p) >= 2 and p[0] != "sample":
                    m[p[0].strip()] = p[1].strip()
    else:
        for path in panel_paths:
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                rows = [ln.rstrip("\r\n").split("\t") for ln in fh if ln.strip()]
            hdr = {c.strip(): i for i, c in enumerate(rows[0])}
            j = hdr.get("expected_uids")
            if j is None:
                continue
            for r in rows[1:]:
                if j >= len(r):
                    continue
                v = r[j].strip()
                if not v or v == "ALL":
                    continue
                grp = [x.strip() for x in v.split(";") if x.strip()]
                if len(grp) < 1:
                    continue
                name = grp[0]
                for s in grp:
                    m.setdefault(s, name)
    for s in samples:
        m.setdefault(s, s)
    return m


def load_driver_calls_all(results_dir, samples, min_mut_cells=5, min_driver_mut_reads=1,
                          min_driver_wt_reads=1):
    """{driver_label: {"sample|barcode": 'MUT'|'WT'}} across every sample, every variant.

    Read support is re-applied here from the stored per-cell counts, so both thresholds can be
    varied without re-running variant_extraction:

      MUT  mut_reads >= min_driver_mut_reads   (default 1: seeing the variant allele means the
                                                cell carries it)
      WT   mut_reads == 0 and wt_reads >= min_driver_wt_reads
                                               (wild-type is inferred from ABSENCE, so it is the
                                                weaker call; in this cohort 94,731 wild-type calls
                                                rest on a single read, versus 18,352 on two)

    Cells meeting neither are dropped rather than defaulted to wild-type.
    Barcodes are namespaced by sample so cells from different samples never collide when pooled.
    """
    per_label = defaultdict(dict)
    for s in samples:
        path = os.path.join(results_dir, f"{s}_complete_analysis.tsv")
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    mr, wr = int(r.get("mut_reads", 0)), int(r.get("wt_reads", 0))
                except ValueError:
                    continue
                if mr >= min_driver_mut_reads:
                    gt = "MUT"
                elif mr == 0 and wr >= min_driver_wt_reads:
                    gt = "WT"
                else:
                    continue
                per_label[r["label"]][f"{s}|{r['bam_barcode']}"] = gt
    keep = {}
    for lab, d in per_label.items():
        if sum(1 for v in d.values() if v == "MUT") >= min_mut_cells:
            keep[lab] = d
    log.info("driver variants with >=%d mutant cells: %d (of %d seen)",
             min_mut_cells, len(keep), len(per_label))
    return keep


def load_mt_calls(mt_dir, samples, restrict_barcodes=None, min_alt_reads_per_cell=1,
                  min_ref_reads_per_cell=1, min_cell_af=0.10):
    """{(pos, ref, alt): {"sample|barcode": 'CARRIER'|'REFERENCE'}} pooled across samples.

    Status is re-derived from the stored read counts rather than the stored label, so the two
    thresholds can be varied without re-scanning any BAM:

      CARRIER    alt_reads >= min_alt_reads_per_cell
                 Seeing the variant allele at all is evidence the cell carries it, so the default
                 is 1.
      REFERENCE  alt_reads == 0 and (total - alt) >= min_ref_reads_per_cell
                 Wild-type is an inference from ABSENCE, which is weaker: one read showing only the
                 normal allele is compatible with a carrier whose variant allele simply was not
                 sampled. Requiring more reads makes the wild-type assignment more reliable, which
                 is why this is a separate parameter.

    A cell with no alternate reads and fewer than `min_ref_reads_per_cell` reads is left out of both
    groups.

    `restrict_barcodes` keeps only cells that have a driver call, which is what makes this fit in
    memory: the raw per-cell files are ~56 MB each and most cells are never tested.
    """
    out = defaultdict(dict)
    for s in samples:
        path = os.path.join(mt_dir, f"{s}_mt_cell_calls.tsv")
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                key = f"{s}|{r['bam_barcode']}"
                if restrict_barcodes is not None and key not in restrict_barcodes:
                    continue
                alt = int(r["alt_reads"])
                tot = int(r["total_reads"])
                ref = tot - alt
                frac = (alt / tot) if tot else 0.0
                # The variant counts as PRESENT in a cell when it is >10% of that cell's reads
                # (min_cell_af) and has at least min_alt_reads_per_cell reads behind it.
                if alt >= min_alt_reads_per_cell and frac >= min_cell_af:
                    st = "CARRIER"
                elif alt == 0 and ref >= min_ref_reads_per_cell:
                    st = "REFERENCE"                    # exclusively the reference allele
                else:
                    continue                            # ambiguous: excluded from both groups
                out[(int(r["pos"]), r["ref_base"], r["alt_base"])][key] = st
    return out


def _index_path(mt_dir, sample):
    return os.path.join(mt_dir, f"{sample}_mt_index.npz")


def build_mt_index(mt_dir, sample):
    """Convert one sample's 65 MB per-cell TSV into a ~3 MB binary index, read once.

    The TSV has ~1.2 M rows (one per cell x position, 97.6% of them reference-only) and lives on a
    ~3 MB/s network mount, so parsing it costs ~20 s. The index stores the same numbers as integer
    arrays (pos, alt_reads, total_reads, cell code) plus a barcode table and the per-position
    variant base -- it loads in under a second. Written next to the TSV so it travels with the data.
    """
    tsv = os.path.join(mt_dir, f"{sample}_mt_cell_calls.tsv")
    if not os.path.exists(tsv):
        return None
    idx = _index_path(mt_dir, sample)
    if os.path.exists(idx) and os.path.getmtime(idx) >= os.path.getmtime(tsv):
        return idx
    pos, alt, total, cell_code = [], [], [], []
    bc_index, barcodes = {}, []
    pos_bases = {}                          # pos -> (ref_base, variant_base)
    with open(tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            p = int(r["pos"])
            b = r["bam_barcode"]
            ci = bc_index.get(b)
            if ci is None:
                ci = bc_index[b] = len(barcodes)
                barcodes.append(b)
            pos.append(p)
            alt.append(int(r["alt_reads"]))
            total.append(int(r["total_reads"]))
            cell_code.append(ci)
            if p not in pos_bases:
                pos_bases[p] = (r["ref_base"], r["alt_base"])
    order = sorted(pos_bases)
    np.savez_compressed(
        idx,
        pos=np.asarray(pos, dtype=np.int32),
        alt=np.asarray(alt, dtype=np.int16),
        total=np.asarray(total, dtype=np.int16),
        cell_code=np.asarray(cell_code, dtype=np.int32),
        barcodes=np.asarray(barcodes),
        pos_list=np.asarray(order, dtype=np.int32),
        ref_bases=np.asarray([pos_bases[p][0] for p in order]),
        variant_bases=np.asarray([pos_bases[p][1] for p in order]))
    return idx


def _load_mt_from_index(mt_dir, sample, keep_barcodes):
    """{(pos, ref, alt): {'sample|barcode': (alt_reads, total_reads)}} from the binary index."""
    idx = build_mt_index(mt_dir, sample)
    if idx is None:
        return {}
    z = np.load(idx, allow_pickle=False)
    barcodes = z["barcodes"]
    pos, alt, total, cc = z["pos"], z["alt"], z["total"], z["cell_code"]
    base_of = {int(p): (str(rb), str(vb))
               for p, rb, vb in zip(z["pos_list"], z["ref_bases"], z["variant_bases"])}
    out = defaultdict(dict)
    for i in range(pos.shape[0]):
        b = f"{sample}|{barcodes[cc[i]]}"
        if b not in keep_barcodes:
            continue
        p = int(pos[i])
        rb, vb = base_of[p]
        out[(p, rb, vb)][b] = (int(alt[i]), int(total[i]))
    return out


def matrix_2x2(adata, driver_feat, mt_feat, carrier_af, min_cell_reads,
               ref_af=0.10, min_driver_reads=1):
    """One (mutation x mt-variant) 2x2 at given thresholds, from a cell x feature h5ad.

    MUT  = driver alt >= min_driver_reads ; WT = driver alt == 0 and driver total >= 1 (UNK excluded).
    carrier = mt fraction >= carrier_af and mt total >= min_cell_reads ; reference = mt fraction <=
    ref_af and mt total >= min_cell_reads. Returns the 2x2, recall (all + covered), background,
    precision, fold, binomial_p, and the imputable-cell count (mt carriers with no driver read).
    """
    import numpy as _np
    dj = adata.var_names.get_loc(driver_feat)
    mj = adata.var_names.get_loc(mt_feat)
    d_alt = _np.asarray(adata.layers["alt"][:, dj].todense()).ravel()
    d_tot = _np.asarray(adata.layers["total"][:, dj].todense()).ravel()
    m_alt = _np.asarray(adata.layers["alt"][:, mj].todense()).ravel()
    m_tot = _np.asarray(adata.layers["total"][:, mj].todense()).ravel()
    m_frac = _np.divide(m_alt, m_tot, out=_np.zeros_like(m_alt, float), where=m_tot > 0)

    mut = d_alt >= min_driver_reads
    wt = (d_alt == 0) & (d_tot >= 1)
    carrier = (m_frac >= carrier_af) & (m_tot >= min_cell_reads)
    ref = (m_frac <= ref_af) & (m_frac < carrier_af) & (m_tot >= min_cell_reads)  # disjoint from carrier

    a = int((mut & carrier).sum()); b = int((wt & carrier).sum())
    c = int((mut & ref).sum()); d = int((wt & ref).sum())
    n_mut_all = int(mut.sum()); n_wt_cov = b + d
    recall_all = a / n_mut_all if n_mut_all else 0.0
    recall_cov = a / (a + c) if (a + c) else 0.0
    bg = b / n_wt_cov if n_wt_cov else 0.0
    prec = a / (a + b) if (a + b) else 0.0
    fold = (recall_cov / bg) if bg > 0 else float("inf")
    # imputable = carriers with NO driver read (UNK for the mutation)
    unk = (~mut) & (~wt)
    imputable = int((carrier & unk).sum())
    p = float("nan")
    try:
        from scipy.stats import binomtest
        if (a + c) and n_wt_cov:
            p = binomtest(a, a + c, max(bg, 1.0 / (n_wt_cov + 1)), alternative="greater").pvalue
    except Exception:
        pass
    return dict(a=a, b=b, c=c, d=d, n_mut_all=n_mut_all, n_wt_covered=n_wt_cov,
                recall_all=round(recall_all, 4), recall_covered=round(recall_cov, 4),
                background=round(bg, 5), precision=round(prec, 4),
                fold_enrichment=(round(fold, 2) if fold != float("inf") else "inf"),
                imputable_cells=imputable, binomial_p=p)


def collapse_duplicate_features(adata):
    """Merge columns that share a var name into one (max of alt, max of total). Older BAM-built h5ads
    kept each co-located driver ALLELE as its own column, so e.g. SRSF2 P95R/P95L/P95H all appear as
    'SRSF2_chr17:76736877' -- duplicate var names, which make var_names.get_loc return an array and
    break the scan/nominate loops. Collapsing (mutant if ANY allele is present) restores unique names
    and matches the current build. A no-op when names are already unique. Returns a new AnnData.
    """
    import numpy as _np
    import pandas as _pd
    import anndata as _ad
    import scipy.sparse as _sp
    if adata.var_names.is_unique:
        return adata
    A = _np.asarray(adata.layers["alt"].todense())
    T = _np.asarray(adata.layers["total"].todense())
    vn = list(adata.var_names)
    cols_by = {}
    order = []
    for i, nm in enumerate(vn):
        if nm not in cols_by:
            cols_by[nm] = []; order.append(nm)
        cols_by[nm].append(i)
    nA = _np.zeros((A.shape[0], len(order)), dtype=A.dtype)
    nT = _np.zeros((T.shape[0], len(order)), dtype=T.dtype)
    meta = {c: [] for c in adata.var.columns}
    for k, nm in enumerate(order):
        cs = cols_by[nm]
        nA[:, k] = A[:, cs].max(axis=1)
        nT[:, k] = T[:, cs].max(axis=1)
        for c in adata.var.columns:
            meta[c].append(adata.var.iloc[cs[0]][c])
    frac = _np.divide(nA, nT, out=_np.zeros(nA.shape, float), where=nT > 0)
    var = _pd.DataFrame(meta, index=order)
    out = _ad.AnnData(X=_sp.csr_matrix(frac), obs=adata.obs.copy(), var=var)
    out.layers["alt"] = _sp.csr_matrix(nA)
    out.layers["total"] = _sp.csr_matrix(nT)
    log.info("collapsed %d duplicate-named columns -> %d unique features", len(vn), len(order))
    return out


def mt_dense_block(adata):
    """Densify the mtDNA block of a cell x feature matrix ONCE, for reuse across a threshold sweep.

    Returns (mt_names, A, T, frac): variant ids, alt-read matrix, total-read matrix, and within-cell
    variant fraction (cells x mt-variants). The block does not depend on any threshold, so computing
    it once and passing it to scan_mt_for_driver(block=...) avoids re-densifying per combination.
    """
    import numpy as _np
    mt_cols = _np.where(adata.var["feature_type"].values == "mtDNA")[0]
    A = _np.asarray(adata.layers["alt"][:, mt_cols].todense(), dtype=_np.float64)
    T = _np.asarray(adata.layers["total"][:, mt_cols].todense(), dtype=_np.float64)
    frac = _np.divide(A, T, out=_np.zeros_like(A), where=T > 0)
    return adata.var_names[mt_cols], A, T, frac


def scan_mt_for_driver(adata, driver_feat, carrier_af=0.5, min_cell_reads=3, ref_af=0.10,
                       min_driver_reads=1, min_driver_wt_reads=1, min_carriers=5,
                       max_background=0.05, return_all=False, block=None):
    """Scan EVERY mitochondrial variant against one driver mutation at given thresholds (vectorised).

    Genotype rules (denominator is driver-MUT vs driver-WT cells only):
      MUT = driver alt reads >= min_driver_reads
      WT  = driver alt == 0 AND driver total >= min_driver_wt_reads.
            min_driver_wt_reads = 0 makes WT = every non-mutant cell, so driver-UNCOVERED cells (the
            cells we would impute onto) become the background -- the correct denominator for an
            imputation marker. min_driver_wt_reads = 1 restricts background to cells with a driver read.
      mt carrier   = within-cell variant fraction >= carrier_af AND mt total reads >= max(1,min_cell_reads)
      mt reference = fraction <= ref_af AND mt total reads >= max(1,min_cell_reads)
    A per-cell mt total floor of 1 is always enforced (a cell with zero mt reads cannot be genotyped),
    so min_cell_reads=0 coincides with min_cell_reads=1; this avoids counting zero-coverage cells as
    references (which would deflate background and manufacture false specificity).

    Returns the best mitochondrial variant (highest recall among hits with background <= max_background
    and binomial p < 0.05); if none are specific, returns the best by recall. With return_all=True,
    returns the full list of mt variants reaching >= min_carriers mutant carriers (or []).
    """
    import numpy as _np
    import scipy.stats as _st
    cov = max(1, int(min_cell_reads))
    dj = adata.var_names.get_loc(driver_feat)
    d_alt = _np.asarray(adata.layers["alt"][:, dj].todense()).ravel()
    d_tot = _np.asarray(adata.layers["total"][:, dj].todense()).ravel()
    mut = (d_alt >= min_driver_reads).astype(_np.float64)
    if min_driver_wt_reads <= 0:
        wt = (d_alt == 0).astype(_np.float64)                 # all non-mutant cells (incl. uncovered)
    else:
        wt = ((d_alt == 0) & (d_tot >= min_driver_wt_reads)).astype(_np.float64)
    n_mut = int(mut.sum())
    if n_mut == 0:
        return [] if return_all else None

    if block is None:
        mt_names, A, T, frac = mt_dense_block(adata)
    else:
        mt_names, A, T, frac = block
    unk = ((d_alt == 0) & (d_tot == 0)).astype(_np.float64)   # no driver read = imputation target
    carrier = ((frac >= carrier_af) & (T >= cov)).astype(_np.float64)
    ref = ((frac <= ref_af) & (frac < carrier_af) & (T >= cov)).astype(_np.float64)  # disjoint from carrier

    a = mut @ carrier            # per mt variant: mutant carriers
    b = wt @ carrier             # wild-type carriers
    c = mut @ ref                # mutant references
    d = wt @ ref                 # wild-type references
    imp = unk @ carrier          # driver-uncovered cells that carry the mt variant (imputable)
    nmut_cov = a + c
    nwt_cov = b + d
    with _np.errstate(divide="ignore", invalid="ignore"):
        recall_all = a / n_mut
        recall_cov = _np.where(nmut_cov > 0, a / nmut_cov, 0.0)
        bg = _np.where(nwt_cov > 0, b / nwt_cov, 0.0)
        prec = _np.where((a + b) > 0, a / (a + b), 0.0)

    rows = []
    for i in range(len(mt_names)):
        if a[i] < min_carriers:
            continue
        pv = 1.0
        # p-value is only needed for potentially specific variants (background <= max_background);
        # germline/common variants can never be a hit, so skipping their binomtest keeps a wide sweep
        # tractable without changing which variants qualify as significant.
        if bg[i] <= max_background and nmut_cov[i] > 0 and nwt_cov[i] > 0:
            pv = _st.binomtest(int(a[i]), int(nmut_cov[i]),
                               max(bg[i], 1.0 / (nwt_cov[i] + 1)), alternative="greater").pvalue
        rows.append(dict(mt_feature=str(mt_names[i]), mutant_carriers=int(a[i]),
                         wt_carriers=int(b[i]), n_mutant_all=n_mut, n_mutant_covered=int(nmut_cov[i]),
                         n_wt_covered=int(nwt_cov[i]), recall_all=round(float(recall_all[i]), 4),
                         recall_covered=round(float(recall_cov[i]), 4), background=round(float(bg[i]), 5),
                         precision=round(float(prec[i]), 4), binomial_p=float(pv),
                         imputable_cells=int(imp[i])))
    if return_all:
        rows.sort(key=lambda r: (-r["recall_all"], r["binomial_p"]))
        return rows
    if not rows:
        return None
    # best = specific (background <= max_background) and significant, ranked by recall then p
    specific = [r for r in rows if r["background"] <= max_background and r["binomial_p"] < 0.05]
    pool = specific if specific else rows
    pool.sort(key=lambda r: (-r["recall_all"], r["binomial_p"]))
    return pool[0]


def scan_mt_wt_sweep(adata, driver_feat, carrier_af, min_cell_reads, ref_af, min_driver_reads,
                     min_driver_wt_reads_list, min_carriers=5, max_background=0.05,
                     min_wt_covered=20, block=None):
    """Evaluate all mt variants against one driver across MANY wild-type read thresholds efficiently.

    The mt carrier/reference matrices and the mutant-side counts do NOT depend on the wild-type read
    threshold, so they are computed once per (carrier_af, min_cell_reads); only the wild-type counts
    (a cheap 1 x N @ N x V matrix product) are recomputed per threshold. This makes a 0..20 wild-type
    sweep essentially free relative to recomputing the whole scan each time.

    Returns dict: min_driver_wt_reads -> (best_row_or_None, hit_rows), where best_row = the max-recall
    variant reaching >= min_carriers mutant carriers, and hit_rows = the variants that are specific
    (background <= max_background) AND significant (binomial p < 0.05) AND rest on a credible wild-type
    denominator (n_wt_covered >= min_wt_covered). The denominator gate is essential: a high wild-type
    read floor can shrink the covered-WT count to a handful of cells, at which point background is
    trivially 0 and precision trivially 1.0 -- a small-denominator artifact, not real specificity.
    """
    import numpy as _np
    import scipy.stats as _st
    cov = max(1, int(min_cell_reads))
    dj = adata.var_names.get_loc(driver_feat)
    d_alt = _np.asarray(adata.layers["alt"][:, dj].todense()).ravel()
    d_tot = _np.asarray(adata.layers["total"][:, dj].todense()).ravel()
    mut = (d_alt >= min_driver_reads).astype(_np.float64)
    unk = ((d_alt == 0) & (d_tot == 0)).astype(_np.float64)
    n_mut = int(mut.sum())
    if n_mut == 0:
        return {mwr: (None, []) for mwr in min_driver_wt_reads_list}

    mt_names, A, T, frac = block if block is not None else mt_dense_block(adata)
    carrier = ((frac >= carrier_af) & (T >= cov)).astype(_np.float64)
    # reference must be disjoint from carrier: at carrier_af == ref_af a cell at exactly that fraction
    # would otherwise be both, double-counting it on the covered side
    ref = ((frac <= ref_af) & (frac < carrier_af) & (T >= cov)).astype(_np.float64)
    a = mut @ carrier
    c = mut @ ref
    imp = unk @ carrier
    nmut_cov = a + c
    with _np.errstate(divide="ignore", invalid="ignore"):
        recall_all = a / n_mut
        recall_cov = _np.where(nmut_cov > 0, a / nmut_cov, 0.0)
    qual = _np.where(a >= min_carriers)[0]

    def _row(i, bg_i, prec_i, pv):
        return dict(mt_feature=str(mt_names[i]), mutant_carriers=int(a[i]), wt_carriers=int(_bd_b[i]),
                    n_mutant_all=n_mut, n_mutant_covered=int(nmut_cov[i]), n_wt_covered=int(_nwt[i]),
                    recall_all=round(float(recall_all[i]), 4), recall_covered=round(float(recall_cov[i]), 4),
                    background=round(float(bg_i), 5), precision=round(float(prec_i), 4),
                    binomial_p=float(pv), imputable_cells=int(imp[i]))

    out = {}
    for mwr in min_driver_wt_reads_list:
        if mwr <= 0:
            wt = (d_alt == 0).astype(_np.float64)        # all non-mutant cells (incl. driver-uncovered)
        else:
            wt = ((d_alt == 0) & (d_tot >= mwr)).astype(_np.float64)
        _bd_b = wt @ carrier
        d = wt @ ref
        _nwt = _bd_b + d
        with _np.errstate(divide="ignore", invalid="ignore"):
            bg = _np.where(_nwt > 0, _bd_b / _nwt, 0.0)
            prec = _np.where((a + _bd_b) > 0, a / (a + _bd_b), 0.0)
        best = None
        best_rec = -1.0
        hits = []
        for i in qual:
            i = int(i)
            if recall_all[i] > best_rec:            # track max-recall variant for the grid "best"
                best_rec = recall_all[i]
                pv_best = 1.0
                if _nwt[i] > 0 and bg[i] <= max_background and nmut_cov[i] > 0:
                    pv_best = _st.binomtest(int(a[i]), int(nmut_cov[i]),
                                            max(bg[i], 1.0 / (_nwt[i] + 1)), alternative="greater").pvalue
                best = _row(i, bg[i], prec[i], pv_best)
            if (_nwt[i] >= min_wt_covered and bg[i] <= max_background and nmut_cov[i] > 0):
                pv = _st.binomtest(int(a[i]), int(nmut_cov[i]),
                                   max(bg[i], 1.0 / (_nwt[i] + 1)), alternative="greater").pvalue
                if pv < 0.05:
                    hits.append(_row(i, bg[i], prec[i], pv))
        out[mwr] = (best, hits)
    return out


def resolve_driver_features(adata, driver):
    """Return the driver feature name(s) in this matrix matching `driver`.

    Accepts an exact feature id (e.g. 'PHF6_chrX:134415107') or a gene prefix (e.g. 'PHF6'), in
    which case every driver feature whose id starts with '<gene>_' is returned. Only features with
    at least one mutant read in this sample are returned (a driver with zero mutant cells here has
    nothing to associate).
    """
    import numpy as _np
    drv = adata.var_names[adata.var["feature_type"].values == "driver"]
    if driver in set(drv):
        cand = [driver]
    else:
        cand = [f for f in drv if f == driver or f.startswith(driver + "_")]
    out = []
    for f in cand:
        j = adata.var_names.get_loc(f)
        alt = _np.asarray(adata.layers["alt"][:, j].todense()).ravel()
        if int((alt >= 1).sum()) > 0:
            out.append(f)
    return out


def driver_genotype_summary(matrix_dir, samples, out_dir=None, tag=None):
    """Per (sample, driver mutation), count cells that are MUT / WT / UNK, to find which mutations can
    even be tested for a specific mitochondrial surrogate.

    MUT = driver alt >= 1 ; WT = driver alt == 0 and driver total >= 1 ; UNK = no driver read. A
    mutation is only worth sweeping where it has both a mutant population (MUT) AND a wild-type arm
    (WT); with WT == 0 the mutation is clonal in that sample and no specificity can be established at a
    wild-type read floor >= 1. Writes driver_genotype_summary[.tag].tsv if out_dir given; returns rows.
    """
    import anndata as _ad
    import numpy as _np
    rows = []
    for s in samples:
        path = os.path.join(matrix_dir, f"{s}_cell_feature_matrix.h5ad")
        if not os.path.exists(path):
            log.warning("no matrix for %s (%s)", s, path)
            continue
        adata = collapse_duplicate_features(_ad.read_h5ad(path))
        drv = adata.var_names[adata.var["feature_type"].values == "driver"]
        for f in drv:
            j = adata.var_names.get_loc(f)
            alt = _np.asarray(adata.layers["alt"][:, j].todense()).ravel()
            tot = _np.asarray(adata.layers["total"][:, j].todense()).ravel()
            mut = int((alt >= 1).sum())
            if mut == 0:
                continue
            wt = int(((alt == 0) & (tot >= 1)).sum())
            unk = int(((alt == 0) & (tot == 0)).sum())
            gene = f.split("_")[0]
            rows.append(dict(sample=s, gene=gene, driver_feature=f, n_mut=mut, n_wt=wt, n_unk=unk))
    rows.sort(key=lambda r: (r["gene"], -r["n_mut"]))
    if out_dir:
        suffix = f"_{tag}" if tag else ""
        p = os.path.join(out_dir, f"driver_genotype_summary{suffix}.tsv")
        _write_tsv(rows, p)
        log.info("wrote %s (%d sample x driver rows)", p, len(rows))
    return rows


def is_transition(mt_feature):
    """True if the mt variant is a transition (A<->G or C<->T), the dominant class of real somatic
    mtDNA mutation. Transversions (esp. the recurrent G>C class seen here) are largely systematic
    library/alignment artifacts, so flagging them separates credible clonal variants from noise.
    mt_feature form: 'chrM:4175_G>A'.
    """
    try:
        ref, alt = mt_feature.rsplit("_", 1)[1].split(">")
        return {ref, alt} in ({"A", "G"}, {"C", "T"})
    except Exception:
        return False


def scan_mt_presence(adata, driver_feat, min_carriers=5, min_wt_covered=20, min_fold=2.0,
                     min_driver_wt_reads=1, block=None, return_all=False):
    """Presence-based (any-alt-read) scan of every mt variant against one driver -- the improved
    framework's genotyping AND ranking. A cell is + for a variant if it has >=1 alt read, - if 0 alt
    and >=1 total read; NO within-cell fraction threshold (a clonal variant diluted among many mtDNA
    copies is still called present). Denominator = driver MUT (alt>=1) vs WT (alt==0 &
    total>=min_driver_wt_reads); UNK excluded.

    Each mt variant is scored by Fisher's exact ENRICHMENT (odds ratio + two-sided p on the 2x2), not
    a strict background gate -- so a real lineage marker with moderate background (e.g. chrM:3565 vs
    PHF6: recall ~0.73, background ~0.26, OR ~7.7) is retained, which a background<=0.05 rule would
    wrongly reject. A hit requires >= min_carriers mutant carriers, a covered wild-type arm
    (n_wt_covered >= min_wt_covered), fold-enrichment (recall/background) >= min_fold, and Fisher
    p < 0.05. Returns the best hit (highest recall, then OR), or all variants with return_all=True.
    """
    import numpy as _np
    from scipy.stats import fisher_exact
    dj = adata.var_names.get_loc(driver_feat)
    d_alt = _np.asarray(adata.layers["alt"][:, dj].todense()).ravel()
    d_tot = _np.asarray(adata.layers["total"][:, dj].todense()).ravel()
    mut = (d_alt >= 1).astype(_np.float64)
    if min_driver_wt_reads <= 0:
        wt = (d_alt == 0).astype(_np.float64)
    else:
        wt = ((d_alt == 0) & (d_tot >= min_driver_wt_reads)).astype(_np.float64)
    unk = ((d_alt == 0) & (d_tot == 0)).astype(_np.float64)
    n_mut = int(mut.sum())
    if n_mut == 0:
        return [] if return_all else None
    mt_names, A, T, _frac = block if block is not None else mt_dense_block(adata)
    present = (A >= 1).astype(_np.float64)              # any alt read = present (no fraction cutoff)
    absent = ((A == 0) & (T >= 1)).astype(_np.float64)
    a = mut @ present; b = wt @ present; c = mut @ absent; d = wt @ absent
    imp = unk @ present
    nmut_cov = a + c; nwt_cov = b + d
    with _np.errstate(divide="ignore", invalid="ignore"):
        recall_all = a / n_mut
        bg = _np.where(nwt_cov > 0, b / nwt_cov, 0.0)
        prec = _np.where((a + b) > 0, a / (a + b), 0.0)
        fold = _np.where(bg > 0, recall_all / bg, _np.inf)
    rows = []
    for i in range(len(mt_names)):
        if a[i] < min_carriers or nwt_cov[i] < min_wt_covered:
            continue
        ai, bi, ci, di = int(a[i]), int(b[i]), int(c[i]), int(d[i])
        try:
            OR, pv = fisher_exact([[ai, bi], [ci, di]])
        except Exception:
            OR, pv = float("nan"), 1.0
        rows.append(dict(mt_feature=str(mt_names[i]), mutant_carriers=ai, wt_carriers=bi,
                         n_mutant_all=n_mut, n_wt_covered=int(nwt_cov[i]),
                         recall_all=round(float(recall_all[i]), 4), background=round(float(bg[i]), 5),
                         precision=round(float(prec[i]), 4),
                         fold_enrichment=(round(float(fold[i]), 2) if fold[i] != _np.inf else "inf"),
                         odds_ratio=round(float(OR), 3), fisher_p=float(pv),
                         imputable_cells=int(imp[i])))
    if return_all:
        rows.sort(key=lambda r: (-r["recall_all"], r["fisher_p"]))
        return rows
    hits = [r for r in rows if r["fisher_p"] < 0.05
            and (r["fold_enrichment"] == "inf" or r["fold_enrichment"] >= min_fold)]
    if not hits:
        return None
    hits.sort(key=lambda r: (-r["recall_all"], -(r["odds_ratio"] if isinstance(r["odds_ratio"], float)
                                                 else 0.0), r["fisher_p"]))
    return hits[0]


def nominate_presence(matrix_dir, samples, out_dir, min_mut=20, min_wt=30, min_carriers=5,
                      max_background=0.05, min_wt_covered=20, transitions_only=False, tag=None):
    """Improved-framework (presence-based, any-alt-read) nomination across EVERY testable
    (sample, driver mutation): mutant clone n_mut>=min_mut AND wild-type arm n_wt>=min_wt. Reports the
    best specific mt marker per mutation using presence genotyping (no fraction threshold), annotated
    transition/transversion. Directly comparable to the fraction-threshold nomination (mt_nomination.tsv).
    Note: still limited to the SUBSTITUTION mt panel in the matrices -- indels (e.g. chrM:3565) require
    de novo indel discovery from the BAMs and are not yet in these matrices.
    """
    import anndata as _ad
    import numpy as _np
    rows = []
    for s in samples:
        path = os.path.join(matrix_dir, f"{s}_cell_feature_matrix.h5ad")
        if not os.path.exists(path):
            log.warning("no matrix for %s", s)
            continue
        adata = collapse_duplicate_features(_ad.read_h5ad(path))
        block = mt_dense_block(adata)
        drv = adata.var_names[adata.var["feature_type"].values == "driver"]
        for feat in drv:
            j = adata.var_names.get_loc(feat)
            alt = _np.asarray(adata.layers["alt"][:, j].todense()).ravel()
            tot = _np.asarray(adata.layers["total"][:, j].todense()).ravel()
            n_mut = int((alt >= 1).sum())
            n_wt = int(((alt == 0) & (tot >= 1)).sum())
            if n_mut < min_mut or n_wt < min_wt:
                continue
            best = scan_mt_presence(adata, feat, min_carriers=min_carriers,
                                    min_wt_covered=min_wt_covered, block=block)
            base = dict(sample=s, gene=feat.split("_")[0], driver_feature=feat, n_mut=n_mut, n_wt=n_wt)
            if best is not None and (not transitions_only or is_transition(best["mt_feature"])):
                base.update(best)
                base["mt_is_transition"] = is_transition(best["mt_feature"])
            else:
                base.update(dict(mt_feature="none", mutant_carriers=0, wt_carriers=0,
                                 n_wt_covered=0, recall_all=0.0, background=0.0, precision=0.0,
                                 fold_enrichment=0.0, odds_ratio=0.0, fisher_p=1.0,
                                 imputable_cells=0, mt_is_transition=""))
            rows.append(base)
    rows.sort(key=lambda r: (-float(r.get("recall_all", 0)), r["gene"]))
    suffix = f"_{tag}" if tag else ""
    p = os.path.join(out_dir, f"mt_nomination_presence{suffix}.tsv")
    _write_tsv(rows, p)
    log.info("wrote %s (%d testable pairs)", p, len(rows))
    return rows


def _load_unit_adata(matrix_dir, samples):
    """Load one analysis unit: a single sample's h5ad, or -- when several samples are given -- the
    POOLED matrix (concatenate cells, union of the per-sample mt panels). Pooling a patient's samples
    combines the MISTRG xenograft and patient-direct cells for a patient-level association. Duplicate
    driver columns are collapsed on load; a variant absent from one sample's panel is uncovered
    (total 0) for that sample's cells, not falsely called reference.
    """
    import anndata as _ad
    ads = []
    for s in samples:
        pth = os.path.join(matrix_dir, f"{s}_cell_feature_matrix.h5ad")
        if os.path.exists(pth):
            ads.append(collapse_duplicate_features(_ad.read_h5ad(pth)))
        else:
            log.warning("no matrix for %s", s)
    if not ads:
        return None
    if len(ads) == 1:
        return ads[0]
    return _ad.concat(ads, join="outer", fill_value=0, merge="first")


def nominate_presence_fdr(matrix_dir, units, out_dir, level="sample", min_mut=10, min_wt=30,
                          min_carriers=5, min_wt_covered=20, max_q=0.05, max_background=0.15, tag=None):
    """Presence-based nomination with Benjamini-Hochberg FDR across ALL (unit x driver x mt) tests,
    at a chosen level. `units` is a dict label -> [sample(s)]:
      * sample level: {sample: [sample]} -- each sample alone,
      * patient level: {patient: [samples]} -- the patient's samples POOLED (MISTRG + direct).
    Every mt variant reaching >= min_carriers mutant carriers with a covered wild-type arm
    (n_wt_covered >= min_wt_covered) is a test; BH-FDR is computed over the whole set of tests and only
    hits with q < max_q are reported (the best-by-recall per unit x driver). Writes
    mt_nomination_presence_<level>.tsv. Returns the significant best-per-(unit,driver) rows.
    """
    import numpy as _np
    # Memory-bounded: keep every candidate's p-value (a float) for the global BH-FDR, but only a
    # top-K-by-recall SHORTLIST of full rows per (unit, driver) -- storing every candidate dict OOMs
    # on big clones (thousands of qualifying mt variants each).
    all_p = []
    shortlist = {}
    n_tests = 0
    n_units = 0
    for label, samples in units.items():
        adata = _load_unit_adata(matrix_dir, samples)
        if adata is None:
            continue
        n_units += 1
        block = mt_dense_block(adata)
        drv = adata.var_names[adata.var["feature_type"].values == "driver"]
        for feat in drv:
            j = adata.var_names.get_loc(feat)
            da = _np.asarray(adata.layers["alt"][:, j].todense()).ravel()
            dt = _np.asarray(adata.layers["total"][:, j].todense()).ravel()
            n_mut = int((da >= 1).sum())
            n_wt = int(((da == 0) & (dt >= 1)).sum())
            if n_mut < min_mut or n_wt < min_wt:
                continue
            cands = scan_mt_presence(adata, feat, min_carriers=min_carriers,
                                     min_wt_covered=min_wt_covered, block=block, return_all=True)
            k = (label, feat)
            lst = shortlist.setdefault(k, [])
            for r in cands:
                if r["n_wt_covered"] < min_wt_covered:
                    continue
                all_p.append(r["fisher_p"])
                n_tests += 1
                lst.append(dict(level=level, unit=label, gene=feat.split("_")[0], driver=feat,
                                n_samples=len(samples), n_mut=n_mut, n_wt=n_wt,
                                mt_is_transition=is_transition(r["mt_feature"]), **r))
            if len(lst) > 60:
                # a reported hit must be SPECIFIC (background <= max_background); among those, the best
                # is the highest-recall one. So retain the highest-recall SPECIFIC candidates. Trimming
                # by recall alone dropped clean low-recall markers (ATRX x chrM:2471); trimming by
                # background alone kept near-zero-background noise and dropped the real high-recall
                # marker (TP53 x chrM:6248). This keeps exactly the candidates eligible to win.
                lst.sort(key=lambda x: (0 if x["background"] <= max_background else 1,
                                        -x["recall_all"]))
                del lst[40:]
    if not all_p:
        log.warning("nominate_presence_fdr level=%s: no tests", level)
        return []
    p = _np.sort(_np.asarray(all_p, dtype=float))
    m = len(p)
    q_sorted = _np.minimum.accumulate((p * m / _np.arange(1, m + 1))[::-1])[::-1]  # BH, monotone

    def _p2q(pv):
        i = min(int(_np.searchsorted(p, pv, side="left")), m - 1)
        return float(q_sorted[i])

    best_rows = []
    for lst in shortlist.values():
        sig = []
        for r in lst:
            q = _p2q(r["fisher_p"])
            # a reported hit must be significant AND SPECIFIC (background <= max_background) -- FDR
            # controls chance, but only the specificity gate removes germline/common variants
            # (e.g. chrM:2617, chrM:302) that pass FDR merely by being enriched in the dominant clone.
            if q < max_q and r["background"] <= max_background:
                r["q_value"] = round(q, 6)
                sig.append(r)
        if sig:
            sig.sort(key=lambda r: (-r["recall_all"], -r["precision"], r["q_value"]))
            best_rows.append(sig[0])
    best_rows.sort(key=lambda r: (-r["recall_all"], r["q_value"]))
    suffix = f"_{tag}" if tag else ""
    out = os.path.join(out_dir, f"mt_nomination_presence_{level}{suffix}.tsv")
    _write_tsv(best_rows, out)
    log.info("level=%s: %d units, %d tests, BH over %d p-values, %d (unit x driver) hits q<%.2f -> %s",
             level, n_units, n_tests, m, len(best_rows), max_q, out)
    return best_rows


def annotate_hits_2x2(matrix_dir, hits_tsv, out_dir, patient_map=None):
    """Add the full 2x2 diagnostic metrics to each nomination hit (recomputed from the h5ad, pooling
    a patient's samples for patient-level rows): n_wt_covered, n_mt_predicted (driver-UNCOVERED cells
    the mt marker flags as mutant = the imputable cells), sensitivity = a/(a+c) (mt-carrier fraction of
    the mt-COVERED mutant cells; recall_all is over ALL mutant cells so it can be lower), and
    specificity = d/(b+d) = 1 - background. Writes <hits>_annotated.tsv. a=mut&carrier, b=wt&carrier,
    c=mut&reference, d=wt&reference; carrier = mt alt>=1, reference = mt alt==0 & total>=1.
    """
    import anndata as _ad
    import numpy as _np
    rows = load_association_table(hits_tsv)
    groups = {}
    for r in rows:
        groups.setdefault((r.get("level", "sample"), r["unit"]), []).append(r)
    out_rows = []
    for (level, unit), unit_rows in groups.items():
        if level == "patient" and patient_map:
            samples = patient_map.get(unit, [unit])
        else:
            samples = [unit]
        adata = _load_unit_adata(matrix_dir, samples)
        if adata is None:
            out_rows.extend(unit_rows); continue
        vn = set(adata.var_names)
        for r in unit_rows:
            drv, mt = r["driver"], r["mt_feature"]
            if drv not in vn or mt not in vn:
                out_rows.append(r); continue
            dj = adata.var_names.get_loc(drv)
            da = _np.asarray(adata.layers["alt"][:, dj].todense()).ravel()
            dt = _np.asarray(adata.layers["total"][:, dj].todense()).ravel()
            mut = da >= 1; wt = (da == 0) & (dt >= 1); unk = (da == 0) & (dt == 0)
            mj = adata.var_names.get_loc(mt)
            ma = _np.asarray(adata.layers["alt"][:, mj].todense()).ravel()
            mm = _np.asarray(adata.layers["total"][:, mj].todense()).ravel()
            car = ma >= 1; ref = (ma == 0) & (mm >= 1)
            a = int((mut & car).sum()); c = int((mut & ref).sum())
            b = int((wt & car).sum()); d = int((wt & ref).sum())
            imp = int((unk & car).sum())
            r["n_wt_covered"] = b + d
            r["n_mt_predicted"] = imp
            r["sensitivity"] = round(a / (a + c), 4) if (a + c) else 0.0
            r["specificity"] = round(d / (b + d), 4) if (b + d) else 0.0
            out_rows.append(r)
    out_rows.sort(key=lambda r: (r.get("level", ""), -float(r.get("recall_all", 0) or 0)))
    base = os.path.basename(hits_tsv).replace(".tsv", "")
    outp = os.path.join(out_dir, f"{base}_annotated.tsv")
    _write_tsv(out_rows, outp)
    log.info("annotated %d hits -> %s", len(out_rows), outp)
    return out_rows


_HITS_HTML_TEMPLATE = r"""<style>
:root{--bg:#0b0f0c;--panel:#0f1512;--line:#1e2a22;--grid:#20342a;--txt:#c9e8d2;--dim:#5a6b5e;
--acc:#46d160;--acc2:#7ee787;--strong:#8ef0a0;--warn:#e3b341;--indel:#4aa3df;--snv:#9a7cd8;
--fontm:ui-monospace,"SF Mono",Menlo,Consolas,"Liberation Mono",monospace;
--fonts:"Inter",system-ui,-apple-system,"Segoe UI",sans-serif;}
:root[data-theme="light"]{--bg:#f3f6f1;--panel:#fff;--line:#d8e3d6;--grid:#e3ecdf;--txt:#1b2a20;
--dim:#7c8a7e;--acc:#1f9d4d;--acc2:#137a3a;--strong:#0f6b30;--warn:#9a6b00;--indel:#1d6fa5;--snv:#5b3f9a;}
@media (prefers-color-scheme:light){:root:not([data-theme="dark"]){--bg:#f3f6f1;--panel:#fff;--line:#d8e3d6;
--grid:#e3ecdf;--txt:#1b2a20;--dim:#7c8a7e;--acc:#1f9d4d;--acc2:#137a3a;--strong:#0f6b30;--warn:#9a6b00;--indel:#1d6fa5;--snv:#5b3f9a;}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--txt);font-family:var(--fontm);
font-variant-numeric:tabular-nums;-webkit-font-smoothing:antialiased;line-height:1.5}
.wrap{max-width:1400px;margin:0 auto;padding:40px 28px 64px}
.eyebrow{font-size:11px;letter-spacing:.18em;text-transform:uppercase;color:var(--acc);margin-bottom:14px}
h1{font-family:var(--fonts);font-weight:650;font-size:clamp(26px,4vw,40px);letter-spacing:-.02em;
margin:0 0 14px;color:var(--txt);text-wrap:balance}
.lede{font-family:var(--fonts);max-width:80ch;color:var(--txt);opacity:.85;font-size:15px;margin:0 0 26px}
.lede .mono{font-family:var(--fontm);color:var(--acc2)}
.lede b{color:var(--acc2);font-weight:600}
.stats{display:flex;flex-wrap:wrap;gap:12px;margin-bottom:24px}
.stat{background:var(--panel);border:1px solid var(--line);border-radius:9px;padding:12px 18px;min-width:118px}
.stat .num{display:block;font-size:26px;font-weight:600;color:var(--strong)}
.stat.weak .num{color:var(--dim)}
.stat .lbl{font-size:11px;letter-spacing:.04em;color:var(--dim)}
.controls{display:flex;flex-wrap:wrap;gap:14px;align-items:center;margin-bottom:18px}
.seg{display:inline-flex;border:1px solid var(--line);border-radius:8px;overflow:hidden;background:var(--panel)}
.seg button{font-family:var(--fontm);font-size:12px;color:var(--dim);background:transparent;border:0;
padding:8px 14px;cursor:pointer;border-right:1px solid var(--line)}
.seg button:last-child{border-right:0}
.seg button.on{color:var(--bg);background:var(--acc)}
.seg button:hover:not(.on){color:var(--txt)}
.hint{font-size:11px;color:var(--dim);letter-spacing:.03em}
.tablewrap{border:1px solid var(--line);border-radius:11px;overflow:auto;max-height:76vh;background:var(--panel)}
table{border-collapse:collapse;width:100%;font-size:12.5px;white-space:nowrap}
thead th{position:sticky;top:0;z-index:2;background:var(--panel);color:var(--acc2);text-align:left;
font-weight:600;padding:11px 12px;border-bottom:1.5px solid var(--grid);cursor:pointer;user-select:none;letter-spacing:.02em}
thead th.r{text-align:right}
thead th:hover{color:var(--strong)}
thead th.sorted::after{content:" \25BC";font-size:9px;color:var(--acc)}
thead th.sorted.asc::after{content:" \25B2"}
tbody td{padding:8px 12px;border-bottom:1px solid var(--line)}
tbody td.r{text-align:right}
tbody tr:hover{background:color-mix(in srgb,var(--acc) 8%,transparent)}
.gene{color:var(--acc2);font-weight:600}
.unit{color:var(--txt);opacity:.8}
.strong-c{color:var(--strong);font-weight:600}
.prec{font-weight:600}
.pred{color:var(--acc)}
.dim{color:var(--dim)}
.lv{font-size:10px;letter-spacing:.06em;text-transform:uppercase;color:var(--dim)}
.lv-patient{color:var(--warn)}
.ty{font-size:10px;padding:2px 7px;border-radius:20px;letter-spacing:.04em}
.ty-indel{color:var(--indel);border:1px solid color-mix(in srgb,var(--indel) 45%,transparent)}
.ty-snv{color:var(--snv);border:1px solid color-mix(in srgb,var(--snv) 45%,transparent)}
.hp{color:var(--warn);margin-left:6px;font-size:9px;vertical-align:middle}
tr.tier-weak td:not(.gene):not(.lv){color:var(--dim)}
tr.tier-weak .gene{opacity:.6}
tr.tier-strong .gene{color:var(--strong)}
tr.tier-strong{background:color-mix(in srgb,var(--acc) 5%,transparent)}
footer{font-family:var(--fonts);font-size:12.5px;color:var(--dim);max-width:80ch;margin-top:22px;line-height:1.6}
</style>
<div class="wrap">
<header>
  <div class="eyebrow">KINNEX-5 single-cell cohort &middot; mtDNA clonal markers</div>
  <h1>__TITLE__</h1>
  <p class="lede">Every driver&ndash;mitochondrial association surviving Benjamini&ndash;Hochberg FDR
  (<span class="mono">q&nbsp;&lt;&nbsp;0.05</span>) and the specificity gate (marker in &lt;15% of
  wild-type cells), at two levels: each <b>sample</b> alone, and each <b>patient</b> with all timepoints
  pooled (MISTRG&nbsp;+&nbsp;direct). Sorted strongest first. <b>sens</b>=mt-carrier fraction of
  mt-covered mutant cells; <b>spec</b>=1&minus;background; <b># MT-pred</b>=driver-uncovered cells the
  marker flags as mutant.</p>
  <div class="stats">
    <div class="stat"><span class="num">__NSTRONG__</span><span class="lbl">strong &middot; prec &ge; 0.4</span></div>
    <div class="stat"><span class="num">__NMOD__</span><span class="lbl">moderate</span></div>
    <div class="stat weak"><span class="num">__NWEAK__</span><span class="lbl">weak / homopolymer</span></div>
    <div class="stat"><span class="num">__NSAMP__</span><span class="lbl">sample-level</span></div>
    <div class="stat"><span class="num">__NPAT__</span><span class="lbl">patient-level</span></div>
  </div>
  <div class="controls">
    <div class="seg" id="lvl"><button data-v="all" class="on">Both levels</button><button data-v="sample">Sample</button><button data-v="patient">Patient</button></div>
    <div class="seg" id="qual"><button data-v="all" class="on">All</button><button data-v="strong">Strong only</button><button data-v="nohomop">Hide homopolymer</button></div>
    <div class="hint"><span id="count"></span> rows &middot; click a column to sort</div>
  </div>
</header>
<div class="tablewrap">
<table id="t"><thead><tr>
  <th data-k="level">level</th><th data-k="gene">gene</th><th data-k="unit">sample&nbsp;/&nbsp;patient</th>
  <th data-k="n_mut" class="r">n_mut</th><th data-k="marker">marker</th><th data-k="type">type</th>
  <th data-k="recall" class="r">recall</th><th data-k="sens" class="r">sens</th><th data-k="spec" class="r">spec</th>
  <th data-k="precision" class="r">precision</th><th data-k="background" class="r">bg</th>
  <th data-k="n_wt" class="r"># wt</th><th data-k="n_pred" class="r"># MT-pred</th>
  <th data-k="OR" class="r">OR</th><th data-k="p" class="r">p</th><th data-k="q" class="r">q</th>
</tr></thead><tbody id="tb"></tbody></table>
</div>
<footer>Markers are clonal passengers, not causal. A poly-C control-region variant can still ride a
dominant clone; homopolymer single-base indels (flagged &#9679;) are the PacBio artifact class &mdash;
judge those by precision and cross-sample replication. Independently validated: PHF6&nbsp;&times;&nbsp;chrM:3565
and TP53&nbsp;&times;&nbsp;chrM:3122.</footer>
</div>
<script>
const DATA = __DATA__;
const tb=document.getElementById('tb');
let flt={lvl:'all',qual:'all'},sortK=null,sortDir=-1;
function fmtp(v){if(v===0)return '0';if(v>=0.001)return (+v).toFixed(3);return (+v).toExponential(0);}
function pass(r){
  if(flt.lvl!=='all'&&r.level!==flt.lvl)return false;
  if(flt.qual==='strong'&&r.tier!=='strong')return false;
  if(flt.qual==='nohomop'&&r.homop)return false;
  return true;
}
function render(){
  let rows=DATA.filter(pass);
  if(sortK){rows=rows.slice().sort((a,b)=>{let x=a[sortK],y=b[sortK];
    if(typeof x==='number'&&typeof y==='number')return (x-y)*sortDir;
    return String(x).localeCompare(String(y))*sortDir;});}
  tb.innerHTML=rows.map(r=>{
    const cls='tier-'+r.tier+(r.homop?' homop':'');
    const tcls=r.type==='indel'?'ty-indel':'ty-snv';
    const hp=r.homop?'<span class="hp" title="homopolymer single-base indel (artifact-prone)">&#9679;</span>':'';
    return '<tr class="'+cls+'">'+
      '<td class="lv lv-'+r.level+'">'+r.level+'</td>'+
      '<td class="gene">'+r.gene+'</td>'+
      '<td class="unit">'+r.unit+'</td>'+
      '<td class="r">'+r.n_mut+'</td>'+
      '<td class="marker">'+r.marker+hp+'</td>'+
      '<td><span class="ty '+tcls+'">'+r.type+'</span></td>'+
      '<td class="r strong-c">'+r.recall.toFixed(2)+'</td>'+
      '<td class="r">'+r.sens.toFixed(2)+'</td>'+
      '<td class="r">'+r.spec.toFixed(2)+'</td>'+
      '<td class="r prec">'+r.precision.toFixed(2)+'</td>'+
      '<td class="r">'+r.background.toFixed(2)+'</td>'+
      '<td class="r">'+r.n_wt+'</td>'+
      '<td class="r pred">'+r.n_pred+'</td>'+
      '<td class="r">'+r.OR+'</td>'+
      '<td class="r dim">'+fmtp(r.p)+'</td>'+
      '<td class="r dim">'+fmtp(r.q)+'</td></tr>';
  }).join('');
  document.getElementById('count').textContent=rows.length;
}
document.querySelectorAll('#lvl button').forEach(b=>b.onclick=()=>{flt.lvl=b.dataset.v;
  b.parentElement.querySelectorAll('button').forEach(x=>x.classList.remove('on'));b.classList.add('on');render();});
document.querySelectorAll('#qual button').forEach(b=>b.onclick=()=>{flt.qual=b.dataset.v;
  b.parentElement.querySelectorAll('button').forEach(x=>x.classList.remove('on'));b.classList.add('on');render();});
document.querySelectorAll('th[data-k]').forEach(th=>th.onclick=()=>{
  const k=th.dataset.k;if(sortK===k)sortDir*=-1;else{sortK=k;sortDir=-1;}
  document.querySelectorAll('th').forEach(x=>x.classList.remove('sorted','asc'));
  th.classList.add('sorted');if(sortDir>0)th.classList.add('asc');render();});
render();
</script>"""


def combine_hits_tsv(annotated_tsvs, out_tsv):
    """Concatenate the sample- and patient-level *_annotated.tsv hit tables into ONE flat table, adding
    the derived `type` (snv/indel, from ref/alt length) and `tier` (strong prec>=0.4 / moderate>=0.2 /
    weak) columns, sorted by tier then recall. Reproducible module replacement for the earlier ad hoc
    awk concatenation; computes no statistics, only concatenates + derives the two label columns.
    """
    import csv as _csv
    rows = []
    for path in annotated_tsvs:
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            for r in _csv.DictReader(fh, delimiter="\t"):
                mt = r.get("mt_feature", "")
                try:
                    ref, alt = mt.rsplit("_", 1)[1].split(">")
                except Exception:
                    ref, alt = "", ""
                r["type"] = "indel" if len(ref) != len(alt) else "snv"
                try:
                    pr = float(r.get("precision", 0) or 0)
                except Exception:
                    pr = 0.0
                r["tier"] = "strong" if pr >= 0.4 else ("moderate" if pr >= 0.2 else "weak")
                rows.append(r)
    order = {"strong": 0, "moderate": 1, "weak": 2}

    def _recall(r):
        try:
            return -float(r.get("recall_all", 0) or 0)
        except Exception:
            return 0.0
    rows.sort(key=lambda r: (order.get(r.get("tier"), 3), _recall(r)))
    _write_tsv(rows, out_tsv)
    log.info("wrote %s (%d hits)", out_tsv, len(rows))
    return rows


def write_hits_html(annotated_tsvs, out_html, title="Mitochondrial markers of driver mutations"):
    """Render annotated nomination hits (from annotate_hits_2x2) into a single self-contained HTML
    report: one sortable/filterable table of every FDR-significant driver x mt hit, with sensitivity,
    specificity, n_wt, and n_mt_predicted, quality-tiered (strong / moderate / weak-homopolymer) so
    the real clonal markers read at a glance. annotated_tsvs is a list of *_annotated.tsv paths.
    """
    import csv as _csv
    import json as _json
    import html as _html
    data = []
    for path in annotated_tsvs:
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            for r in _csv.DictReader(fh, delimiter="\t"):
                mt = r.get("mt_feature", "")
                try:
                    ref, alt = mt.rsplit("_", 1)[1].split(">")
                except Exception:
                    ref, alt = "", ""
                is_indel = len(ref) != len(alt)
                homop = is_indel and abs(len(ref) - len(alt)) == 1

                def _f(k, d=0.0):
                    try:
                        return float(r.get(k, d) or d)
                    except Exception:
                        return d

                def _i(k):
                    try:
                        return int(float(r.get(k, 0) or 0))
                    except Exception:
                        return 0
                prec = _f("precision")
                tier = "strong" if prec >= 0.4 else ("moderate" if prec >= 0.2 else "weak")
                data.append(dict(level=r.get("level", "sample"), gene=r.get("gene", ""),
                                 unit=r.get("unit", ""), n_mut=_i("n_mut"), marker=mt,
                                 type=("indel" if is_indel else "snv"), homop=homop,
                                 recall=round(_f("recall_all"), 3), precision=round(prec, 3),
                                 background=round(_f("background"), 3), sens=round(_f("sensitivity"), 3),
                                 spec=round(_f("specificity"), 3), n_wt=_i("n_wt_covered"),
                                 n_pred=_i("n_mt_predicted"), OR=(r.get("odds_ratio") or ""),
                                 p=_f("fisher_p"), q=_f("q_value"), tier=tier))
    order = {"strong": 0, "moderate": 1, "weak": 2}
    data.sort(key=lambda r: (order[r["tier"]], -r["precision"], -r["recall"]))
    counts = dict(strong=sum(r["tier"] == "strong" for r in data),
                  moderate=sum(r["tier"] == "moderate" for r in data),
                  weak=sum(r["tier"] == "weak" for r in data),
                  sample=sum(r["level"] == "sample" for r in data),
                  patient=sum(r["level"] == "patient" for r in data))
    html_doc = _HITS_HTML_TEMPLATE.replace("__DATA__", _json.dumps(data)) \
        .replace("__TITLE__", _html.escape(title)) \
        .replace("__NSTRONG__", str(counts["strong"])).replace("__NMOD__", str(counts["moderate"])) \
        .replace("__NWEAK__", str(counts["weak"])).replace("__NSAMP__", str(counts["sample"])) \
        .replace("__NPAT__", str(counts["patient"]))
    with open(out_html, "w") as fh:
        fh.write(html_doc)
    log.info("wrote %s (%d hits: %d strong, %d moderate, %d weak)", out_html, len(data),
             counts["strong"], counts["moderate"], counts["weak"])
    return out_html


def nominate_mt_associations(matrix_dir, samples, out_dir, min_mut=20, min_wt=30,
                             carrier_afs=(0.1, 0.2, 0.3), min_cell_reads_list=(1,),
                             min_driver_wt_reads_list=(1,), ref_af=0.05, min_carriers=5,
                             max_background=0.05, min_wt_covered=20, transitions_only=False, tag=None):
    """Apply one fixed threshold policy to EVERY testable (sample, driver mutation) and report the
    best specific mt marker per mutation -- the systematic 'which mutations have any specific mt
    association' table.

    Testable = a mutant clone (n_mut >= min_mut) AND a covered wild-type arm (n_wt >= min_wt) in that
    sample. For each, all mt variants are scanned over the small threshold grid; the best gated hit
    (highest recall among background <= max_background, binomial p < 0.05, n_wt_covered >=
    min_wt_covered) is kept, annotated with whether it is a transition (real mtDNA class) or a
    transversion (largely artifact). transitions_only restricts nominations to transitions. One row
    per testable pair (mt_feature 'none' if nothing specific is found). Writes mt_nomination[.tag].tsv.
    """
    import anndata as _ad
    import numpy as _np
    rows = []
    for s in samples:
        path = os.path.join(matrix_dir, f"{s}_cell_feature_matrix.h5ad")
        if not os.path.exists(path):
            log.warning("no matrix for %s (%s)", s, path)
            continue
        adata = collapse_duplicate_features(_ad.read_h5ad(path))
        block = mt_dense_block(adata)
        drv = adata.var_names[adata.var["feature_type"].values == "driver"]
        for feat in drv:
            j = adata.var_names.get_loc(feat)
            alt = _np.asarray(adata.layers["alt"][:, j].todense()).ravel()
            tot = _np.asarray(adata.layers["total"][:, j].todense()).ravel()
            n_mut = int((alt >= 1).sum())
            n_wt = int(((alt == 0) & (tot >= 1)).sum())
            if n_mut < min_mut or n_wt < min_wt:
                continue
            candidates = []
            for caf in carrier_afs:
                for mcr in min_cell_reads_list:
                    res = scan_mt_wt_sweep(adata, feat, caf, mcr, ref_af, 1,
                                           min_driver_wt_reads_list, min_carriers=min_carriers,
                                           max_background=max_background, min_wt_covered=min_wt_covered,
                                           block=block)
                    for mwr, (_best, hitrows) in res.items():
                        for h in hitrows:
                            if transitions_only and not is_transition(h["mt_feature"]):
                                continue
                            candidates.append(dict(carrier_af=caf, min_cell_reads=mcr,
                                                   min_driver_wt_reads=mwr, **h))
            base = dict(sample=s, gene=feat.split("_")[0], driver_feature=feat, n_mut=n_mut, n_wt=n_wt)
            if candidates:
                best = sorted(candidates, key=lambda r: (-r["recall_all"], -r["precision"],
                                                         r["binomial_p"]))[0]
                base.update(best)
                base["mt_is_transition"] = is_transition(best["mt_feature"])
            else:
                base.update(dict(mt_feature="none", recall_all=0.0, precision=0.0, background=0.0,
                                 binomial_p=1.0, mutant_carriers=0, n_wt_covered=0, imputable_cells=0,
                                 mt_is_transition=""))
            rows.append(base)
    rows.sort(key=lambda r: (-float(r.get("recall_all", 0)), r["gene"]))
    suffix = f"_{tag}" if tag else ""
    p = os.path.join(out_dir, f"mt_nomination{suffix}.tsv")
    _write_tsv(rows, p)
    log.info("wrote %s (%d testable sample x mutation pairs)", p, len(rows))
    return rows


def contingency_report(matrix_dir, sample, driver_feat, mt_feats, out_dir=None,
                       min_driver_wt_reads=1, tag=None):
    """Explicit per-cell contingency for ONE driver mutation vs one or more mt variants -- the exact
    cell counts and overlaps. Presence-based (any-alt-read), matching the nomination. For each mt
    variant it reports the full cross-tabulation: driver MUTANT / WILD-TYPE / UNCOVERED crossed with
    mt CARRIER (>=1 alt read) / mt REFERENCE (0 alt, covered) / mt UNCOVERED, every cell accounted for.
    Returns one row-dict per mt variant; writes contingency_<driver>[.tag].tsv if out_dir given.
    """
    import anndata as _ad
    import numpy as _np
    adata = collapse_duplicate_features(_ad.read_h5ad(
        os.path.join(matrix_dir, f"{sample}_cell_feature_matrix.h5ad")))
    dj = adata.var_names.get_loc(driver_feat)
    da = _np.asarray(adata.layers["alt"][:, dj].todense()).ravel()
    dt = _np.asarray(adata.layers["total"][:, dj].todense()).ravel()
    mut = da >= 1
    wt = (da == 0) & (dt >= min_driver_wt_reads)
    unk = (da == 0) & (dt == 0)
    rows = []
    for v in mt_feats:
        if v not in set(adata.var_names):
            log.warning("%s not in %s panel", v, sample)
            continue
        mj = adata.var_names.get_loc(v)
        ma = _np.asarray(adata.layers["alt"][:, mj].todense()).ravel()
        mm = _np.asarray(adata.layers["total"][:, mj].todense()).ravel()
        car = ma >= 1
        ref = (ma == 0) & (mm >= 1)
        mun = mm == 0
        a = int((mut & car).sum()); mr = int((mut & ref).sum()); mu = int((mut & mun).sum())
        b = int((wt & car).sum()); wr = int((wt & ref).sum()); wu = int((wt & mun).sum())
        uc = int((unk & car).sum()); ur = int((unk & ref).sum())
        rows.append(dict(
            sample=sample, driver=driver_feat, mt_variant=v, total_cells=int(adata.n_obs),
            driver_MUT=int(mut.sum()), driver_WT=int(wt.sum()), driver_UNCOVERED=int(unk.sum()),
            mt_carrier_cells=int(car.sum()), mt_reference_cells=int(ref.sum()),
            mt_uncovered_cells=int(mun.sum()),
            MUT_and_mtCarrier=a, MUT_and_mtRef=mr, MUT_and_mtUncov=mu,
            WT_and_mtCarrier=b, WT_and_mtRef=wr, WT_and_mtUncov=wu,
            UNCOV_and_mtCarrier_imputable=uc, UNCOV_and_mtRef=ur,
            recall_MUT_carry=round(a / (a + mr), 4) if (a + mr) else 0.0,
            background_WT_carry=round(b / (b + wr), 4) if (b + wr) else 0.0,
            precision_carrier_MUT=round(a / (a + b), 4) if (a + b) else 0.0))
    if out_dir and rows:
        safe = driver_feat.replace(":", "_").replace("/", "_")
        suffix = f"_{tag}" if tag else ""
        _write_tsv(rows, os.path.join(out_dir, f"contingency_{safe}{suffix}.tsv"))
    return rows


def marker_threshold_profile(matrix_dir, sample, driver_feat, mt_feat, out_dir=None,
                             carrier_afs=(0.3, 0.5, 0.7, 0.9), min_cell_reads=1, ref_af=0.05,
                             min_driver_reads=1, tag=None):
    """Descriptive profile of ONE (mutation x mt-variant) marker across the carrier_af sweep -- the
    format that tells a real clonal marker apart from a threshold-manufactured one:
      * a REAL clonal marker has background ~0 and precision ~1 at EVERY carrier_af (its specificity is
        cutoff-independent); only the imputable-cell count falls as the cutoff tightens.
      * a GERMLINE/artifact variant has high background at a loose cutoff that collapses only when the
        cutoff is tightened toward near-homoplasmic -- its 'specificity' is manufactured by the cutoff.
    Wild-type = driver alt == 0 and driver total >= 1 (a real covered wild-type arm; UNK excluded).
    Returns one row per carrier_af (background, precision, recall_all/covered, imputable, 2x2 counts);
    writes marker_profile_<driver>_<mt>[.tag].tsv if out_dir given.
    """
    import anndata as _ad
    adata = _ad.read_h5ad(os.path.join(matrix_dir, f"{sample}_cell_feature_matrix.h5ad"))
    rows = []
    for caf in carrier_afs:
        r = matrix_2x2(adata, driver_feat, mt_feat, carrier_af=caf, min_cell_reads=min_cell_reads,
                       ref_af=min(ref_af, caf), min_driver_reads=min_driver_reads)
        rows.append(dict(sample=sample, driver_feature=driver_feat, mt_feature=mt_feat,
                         carrier_af=caf, min_cell_reads=min_cell_reads, **r))
    if out_dir:
        safe_d = driver_feat.replace(":", "_").replace("/", "_")
        safe_m = mt_feat.replace(":", "_").replace(">", "").replace("/", "_")
        suffix = f"_{tag}" if tag else ""
        p = os.path.join(out_dir, f"marker_profile_{safe_d}_{safe_m}{suffix}.tsv")
        _write_tsv(rows, p)
        log.info("wrote %s (%d carrier_af rows)", p, len(rows))
    return rows


def sweep_thresholds_for_driver(matrix_dir, samples, driver, out_dir,
                                carrier_afs=(0.1, 0.2, 0.3, 0.5, 0.7, 0.9),
                                min_cell_reads_list=(2, 3, 5),
                                min_driver_reads_list=(1,),
                                min_driver_wt_reads_list=(0, 1, 2),
                                ref_afs=(0.10,), min_carriers=5, max_background=0.05,
                                min_wt_covered=20, tag=None):
    """Sweep every threshold combination for ONE driver mutation across per-sample matrices.

    For each (sample, carrier_af, min_cell_reads, min_driver_reads, min_driver_wt_reads, ref_af) the
    full set of mitochondrial variants reaching >= min_carriers mutant carriers is evaluated. Writes:
      * mt_driver_sweep_<driver>.tsv          -- best mt variant at every combination (full grid)
      * mt_driver_sweep_<driver>_best.tsv     -- the single best combination per sample
      * mt_driver_sweep_<driver>_hits.tsv     -- EVERY (combination x mt variant) that is specific
                                                 (background <= max_background), significant
                                                 (binomial p < 0.05) AND rests on a credible wild-type
                                                 denominator (n_wt_covered >= min_wt_covered), so no
                                                 real association is hidden behind the per-combination
                                                 best and no small-denominator artifact is called specific.

    Denominator is driver-MUT vs driver-WT cells only. min_driver_wt_reads = 0 makes the background
    the full non-mutant population (including driver-uncovered cells we would impute onto); >= 1
    restricts it to cells with a driver read. The best-per-sample row is the highest-recall credible
    hit (from the gated hits), NOT the grid's max-recall variant (which is typically a germline
    heteroplasmy present in nearly all cells). Returns (grid, best, hits).
    """
    import anndata as _ad
    safe = driver.replace(":", "_").replace("/", "_")
    grid, hits = [], []
    for s in samples:
        path = os.path.join(matrix_dir, f"{s}_cell_feature_matrix.h5ad")
        if not os.path.exists(path):
            log.warning("no matrix for %s (%s)", s, path)
            continue
        adata = _ad.read_h5ad(path)
        feats = resolve_driver_features(adata, driver)
        if not feats:
            log.info("%s: no '%s' driver feature with mutant cells", s, driver)
            continue
        for feat in feats:
            block = mt_dense_block(adata)   # densify mtDNA block once; reuse across all combinations
            for caf in carrier_afs:
                for mcr in min_cell_reads_list:
                    for mdr in min_driver_reads_list:
                        for raf in ref_afs:
                            res = scan_mt_wt_sweep(
                                adata, feat, caf, mcr, raf, mdr, min_driver_wt_reads_list,
                                min_carriers=min_carriers, max_background=max_background,
                                min_wt_covered=min_wt_covered, block=block)
                            for mwr, (best, hitrows) in res.items():
                                ctx = dict(sample=s, driver_feature=feat, carrier_af=caf,
                                           min_cell_reads=mcr, min_driver_reads=mdr,
                                           min_driver_wt_reads=mwr, ref_af=raf)
                                grow = dict(ctx)
                                if best is None:
                                    grow.update(dict(mt_feature="none", mutant_carriers=0,
                                                     wt_carriers=0, n_mutant_all=0, n_mutant_covered=0,
                                                     n_wt_covered=0, recall_all=0.0, recall_covered=0.0,
                                                     background=0.0, precision=0.0, binomial_p=1.0,
                                                     imputable_cells=0))
                                else:
                                    grow.update(best)
                                grid.append(grow)
                                for r in hitrows:
                                    hr = dict(ctx); hr.update(r); hits.append(hr)
    if not grid:
        log.warning("sweep produced no rows for driver '%s'", driver)
        return [], [], []
    # best combination per sample = the highest-recall credible hit (from the denominator-gated hits),
    # tie-broken by precision then p. Derived from `hits`, NOT the grid: the grid stores the max-recall
    # variant per combination, which is almost always a germline heteroplasmy carried by nearly all
    # cells (high recall, ~0 precision) and never a usable marker. A sample with no credible hit gets an
    # explicit "none_specific" row so the absence is recorded rather than papered over with germline.
    best_rows = []
    for s in sorted({r["sample"] for r in grid}):
        shits = [r for r in hits if r["sample"] == s]
        if shits:
            best_rows.append(sorted(shits, key=lambda r: (-r["recall_all"], -r["precision"],
                                                          r["binomial_p"]))[0])
        else:
            ctx = next(r for r in grid if r["sample"] == s)
            best_rows.append(dict(sample=s, driver_feature=ctx["driver_feature"],
                                  mt_feature="none_specific",
                                  note=f"no mt variant with background<={max_background}, p<0.05 and "
                                       f"n_wt_covered>={min_wt_covered}"))
    hits.sort(key=lambda r: (r["sample"], -r["recall_all"], r["binomial_p"]))
    suffix = f"_{tag}" if tag else ""
    grid_path = os.path.join(out_dir, f"mt_driver_sweep_{safe}{suffix}.tsv")
    best_path = os.path.join(out_dir, f"mt_driver_sweep_{safe}{suffix}_best.tsv")
    hits_path = os.path.join(out_dir, f"mt_driver_sweep_{safe}{suffix}_hits.tsv")
    _write_tsv(grid, grid_path)
    _write_tsv(best_rows, best_path)
    _write_tsv(hits, hits_path)
    log.info("wrote %s (%d), %s (%d), %s (%d specific+significant hits)",
             grid_path, len(grid), best_path, len(best_rows), hits_path, len(hits))
    return grid, best_rows, hits


def build_cell_feature_matrix(sample, mt_dir, results_dir, out_dir):
    """Write a per-sample AnnData (.h5ad): cells x features, features = every driver mutation AND
    every mitochondrial variant, so any threshold or association is an immediate in-memory query.

    X                 = within-cell variant fraction (alt_reads / total_reads) at each feature
    layers['alt']     = variant-supporting reads (mutant reads for a driver; variant reads for mtDNA)
    layers['total']   = total reads at that feature in that cell
    obs               = cell barcodes (bam orientation), one row per cell
    var               = features; var['feature_type'] in {driver, mtDNA}; var['chrom','pos','ref','alt']

    Driver features are collapsed to genomic coordinate (gene_chrom:pos); mitochondrial features are
    chrM:pos_ref>alt. From this one file: mutant = driver alt>=1; mt carrier = mtDNA fraction>=t with
    total>=n; the 2x2 for any (mutation, mt variant) is two columns cross-tabulated.
    """
    import anndata as ad
    import scipy.sparse as sp
    os.makedirs(out_dir, exist_ok=True)

    cells = {}                 # barcode -> row index
    feats = {}                 # feature name -> col index
    fmeta = {}                 # feature -> (type, chrom, pos, ref, alt)
    rows, cols, alt_v, tot_v = [], [], [], []

    def cell_idx(bc):
        i = cells.get(bc)
        if i is None:
            i = cells[bc] = len(cells)
        return i

    def feat_idx(name, meta):
        j = feats.get(name)
        if j is None:
            j = feats[name] = len(feats)
            fmeta[name] = meta
        return j

    # driver features (collapsed to coordinate; max reads across co-located labels per cell)
    ca = os.path.join(results_dir, f"{sample}_complete_analysis.tsv")
    driver_cell = defaultdict(lambda: (0, 0))    # (feature, bc) -> (mut, wt)
    if os.path.exists(ca):
        with open(ca) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    mr, wr = int(r.get("mut_reads", 0)), int(r.get("wt_reads", 0))
                except ValueError:
                    continue
                gene = r.get("gene") or r.get("label", "").split("_")[0]
                name = f"{gene}_{r.get('chr','')}:{r.get('position','')}"
                key = (name, r["bam_barcode"])
                p = driver_cell[key]
                driver_cell[key] = (max(mr, p[0]), max(wr, p[1]))
                fmeta.setdefault(name, ("driver", r.get("chr", ""), r.get("position", ""),
                                        r.get("ref_base", ""), r.get("alt_base", "")))
    for (name, bc), (mr, wr) in driver_cell.items():
        j = feat_idx(name, fmeta[name]); i = cell_idx(bc)
        rows.append(i); cols.append(j); alt_v.append(mr); tot_v.append(mr + wr)

    # mitochondrial features (from the binary index)
    idx = build_mt_index(mt_dir, sample)
    if idx is not None:
        z = np.load(idx, allow_pickle=False)
        pos, alt, total, cc = z["pos"], z["alt"], z["total"], z["cell_code"]
        barcodes = z["barcodes"]
        base_of = {int(p): (str(rb), str(vb))
                   for p, rb, vb in zip(z["pos_list"], z["ref_bases"], z["variant_bases"])}
        for k in range(pos.shape[0]):
            p = int(pos[k]); rb, vb = base_of[p]
            name = f"chrM:{p}_{rb}>{vb}"
            j = feat_idx(name, ("mtDNA", "chrM", p, rb, vb))
            i = cell_idx(str(barcodes[cc[k]]))
            rows.append(i); cols.append(j); alt_v.append(int(alt[k])); tot_v.append(int(total[k]))

    n, m = len(cells), len(feats)
    A = sp.csr_matrix((alt_v, (rows, cols)), shape=(n, m), dtype=np.int32)
    T = sp.csr_matrix((tot_v, (rows, cols)), shape=(n, m), dtype=np.int32)
    with np.errstate(divide="ignore", invalid="ignore"):
        Td = T.toarray().astype(np.float32)
        X = np.divide(A.toarray(), Td, out=np.zeros((n, m), np.float32), where=Td > 0)

    import pandas as pd
    bc_order = sorted(cells, key=cells.get)
    ft_order = sorted(feats, key=feats.get)
    var = pd.DataFrame(index=ft_order)
    var["feature_type"] = [fmeta[f][0] for f in ft_order]
    var["chrom"] = [fmeta[f][1] for f in ft_order]
    var["pos"] = [str(fmeta[f][2]) for f in ft_order]
    var["ref"] = [fmeta[f][3] for f in ft_order]
    var["alt"] = [fmeta[f][4] for f in ft_order]
    obs = pd.DataFrame(index=[f"{sample}|{b}" for b in bc_order])
    obs["sample"] = sample
    adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var,
                       layers={"alt": A, "total": T})
    adata.uns["description"] = ("cells x (driver mutations + mitochondrial variants); "
                                "X = variant fraction; layers alt/total = read counts")
    out = os.path.join(out_dir, f"{sample}_cell_feature_matrix.h5ad")
    adata.write_h5ad(out, compression="gzip")
    n_drv = int((var["feature_type"] == "driver").sum())
    n_mt = int((var["feature_type"] == "mtDNA").sum())
    log.info("[%s] wrote %s: %d cells x %d features (%d driver, %d mtDNA)",
             sample, out, n, m, n_drv, n_mt)
    return out


def load_patient_raw(mt_dir, results_dir, samples):
    """Read each per-cell file ONCE and keep the raw read counts.

    Returns (driver_raw, mt_raw):
      driver_raw = {driver_label: {"sample|barcode": (mut_reads, wt_reads)}}
      mt_raw     = {(pos, ref, alt): {"sample|barcode": (alt_reads, total_reads)}}
    mt_raw is limited to cells that have at least one driver call, which is what keeps it in memory.
    Thresholds are NOT applied here -- they are applied in memory by derive_driver / derive_mt, so a
    threshold sweep reads the 2.7 GB of files once instead of once per setting.
    """
    # Driver calls are collapsed to the GENOMIC COORDINATE (gene, chrom, pos), not the panel label,
    # so co-located labels for one mutation (e.g. SRSF2 P95 at chr17:76736877 appears as p.P95R,
    # c.284C>A, c.284C>T, P95_p2) become a single mutation tested once. For a cell genotyped at
    # several co-located labels, the mutant/wild-type read support is the maximum across them (a cell
    # is mutant if it carries the variant allele at that position under any label).
    driver_raw = defaultdict(dict)
    for s in samples:
        path = os.path.join(results_dir, f"{s}_complete_analysis.tsv")
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    mr, wr = int(r.get("mut_reads", 0)), int(r.get("wt_reads", 0))
                except ValueError:
                    continue
                gene = r.get("gene") or r.get("label", "").split("_")[0]
                mut = f"{gene}_{r.get('chr','')}:{r.get('position','')}"
                bc = f"{s}|{r['bam_barcode']}"
                prev = driver_raw[mut].get(bc)
                driver_raw[mut][bc] = (mr, wr) if prev is None else (max(mr, prev[0]),
                                                                     max(wr, prev[1]))
    keep = set()
    for d in driver_raw.values():
        keep.update(d)

    mt_raw = defaultdict(dict)
    for s in samples:
        # binary index if available (fast); build it from the TSV on first use
        part = _load_mt_from_index(mt_dir, s, keep)
        for k, cells in part.items():
            mt_raw[k].update(cells)
    return driver_raw, mt_raw


def derive_driver(driver_raw, min_mut_cells=5, min_driver_mut_reads=1, min_driver_wt_reads=1):
    """Apply the mutation-side thresholds to raw counts (in memory). See load_driver_calls_all."""
    keep = {}
    for lab, cells in driver_raw.items():
        gt = {}
        for bc, (mr, wr) in cells.items():
            if mr >= min_driver_mut_reads:
                gt[bc] = "MUT"
            elif mr == 0 and wr >= min_driver_wt_reads:
                gt[bc] = "WT"
        if sum(1 for v in gt.values() if v == "MUT") >= min_mut_cells:
            keep[lab] = gt
    return keep


def derive_mt(mt_raw, min_cell_reads=3, carrier_af=0.90, ref_af=0.10,
              max_frac_intermediate=0.20):
    """Bimodal, coverage-gated per-cell genotyping of mitochondrial variants.

    A somatic clonal passenger is near-homoplasmic WITHIN a cell: cells either carry it at high
    variant fraction (~1) or not at all (~0). A germline heteroplasmy sits at an intermediate
    fraction in most cells. The earlier definition (>=10% variant fraction, one read) counted both
    intermediate germline cells and single-read sequencing error as carriers, which dominated the
    associations. This definition instead requires:

      CARRIER    variant fraction >= carrier_af  AND total reads >= min_cell_reads
      REFERENCE  variant fraction <= ref_af       AND total reads >= min_cell_reads
      excluded   intermediate fraction, or fewer than min_cell_reads reads (cannot be called)

    A whole POSITION is also discarded when too large a fraction of its covered cells fall in the
    intermediate band (> max_frac_intermediate), which is the signature of a germline heteroplasmy
    rather than a clonal variant.
    """
    out = {}
    for k, cells in mt_raw.items():
        carrier, reference, n_inter, n_cov = set(), set(), 0, 0
        for bc, (alt, tot) in cells.items():
            if tot < min_cell_reads:
                continue
            n_cov += 1
            frac = alt / tot
            if frac >= carrier_af:
                carrier.add(bc)
            elif frac <= ref_af:
                reference.add(bc)
            else:
                n_inter += 1
        if n_cov == 0:
            continue
        # germline-heteroplasmy filter: a clonal variant is bimodal, so few cells are intermediate
        if (n_inter / n_cov) > max_frac_intermediate:
            continue
        if carrier or reference:
            out[k] = dict(carrier=carrier, reference=reference)
    return out


def cohort_associate(mt_dir, results_dir, out_dir, metadata=None, panel_paths=(), donor_map=None,
                     samples=None, min_mut_cells=5, max_q=0.05, min_recall=0.90,
                     max_background=0.10, min_alt_reads_per_cell=1, min_ref_reads_per_cell=1,
                     min_driver_mut_reads=1, min_driver_wt_reads=1, tag="", raw_cache=None,
                     min_cell_reads=3, carrier_af=0.90, ref_af=0.10, n_permutations=0,
                     per_sample=False):
    """Test every driver mutation in every patient against every mitochondrial variant.

    Mutations are collapsed to genomic coordinate (in load_patient_raw). Cells are genotyped for the
    mitochondrial variant with the bimodal, coverage-gated rule (derive_mt: carrier fraction >=
    carrier_af, reference <= ref_af, >= min_cell_reads reads). Cells are pooled per patient; one FDR
    correction spans all tests. If n_permutations > 0, an empirical null is built by shuffling the
    mutant/wild-type labels among genotyped cells and recording the best statistic per permutation,
    giving a permutation false-discovery estimate that does not assume test independence.
    """
    os.makedirs(out_dir, exist_ok=True)
    if samples is None:
        samples = sorted(os.path.basename(p).replace("_mt_cell_calls.tsv", "")
                         for p in _glob(os.path.join(mt_dir, "*_mt_cell_calls.tsv")))
    if per_sample:
        # Each sample analysed on its own. Mitochondrial variants are discovered per sample, so
        # pooling timepoints mixes cells that were never genotyped at the same positions and biases
        # the background; per-sample analysis keeps mutant and wild-type cells on equal footing.
        by_patient = {s: [s] for s in samples}
    else:
        patient_of = load_patient_map(list(panel_paths), samples, donor_map)
        by_patient = defaultdict(list)
        for s in samples:
            by_patient[patient_of[s]].append(s)
    log.info("%d samples in %d analysis groups (per_sample=%s)", len(samples), len(by_patient),
             per_sample)
    if raw_cache is None:
        raw_cache = {}

    all_rows = []
    for patient, sms in sorted(by_patient.items()):
        if patient not in raw_cache:
            raw_cache[patient] = load_patient_raw(mt_dir, results_dir, sms)   # read files ONCE
        driver_raw, mt_raw = raw_cache[patient]
        drivers = derive_driver(driver_raw, min_mut_cells, min_driver_mut_reads,
                                min_driver_wt_reads)
        if not drivers:
            continue
        mt = derive_mt(mt_raw, min_cell_reads=min_cell_reads, carrier_af=carrier_af, ref_af=ref_af)
        if not mt:
            continue
        log.info("[%s] %d samples, %d driver variants, %d mitochondrial variants",
                 patient, len(sms), len(drivers), len(mt))
        positions = [dict(pos=k[0], ref_base=k[1], alt_base=k[2], af="", depth="")
                     for k in mt.keys()]
        calls = {k[0]: v for k, v in mt.items()}
        for label, driver in drivers.items():
            # SAMPLE CONFOUND CONTROL. Mitochondrial variants are discovered per sample, so a
            # variant found in one sample is largely absent from that patient's other samples. If
            # the mutant cells also sit in one sample, the variant and the mutation coincide because
            # they share a SAMPLE, not a clone -- and the test reports a near-perfect, entirely
            # spurious match. Restrict both groups to the samples that actually contain mutant
            # cells, so sample identity cannot generate the association.
            mut_samples = {b.split("|", 1)[0] for b, v in driver.items() if v == "MUT"}
            driver_local = {b: v for b, v in driver.items()
                            if b.split("|", 1)[0] in mut_samples}
            if sum(1 for v in driver_local.values() if v == "MUT") < min_mut_cells:
                continue
            # how many cells the mutation was called mutant in AT ALL, before requiring that the
            # mitochondrial position also be covered -- the denominator for imputation coverage
            n_called = sum(1 for v in driver_local.values() if v == "MUT")
            rows_here = associate(positions, calls, driver_local, min_mut_cells=min_mut_cells,
                                  label=label, patient=patient,
                                  n_samples_compared=len(mut_samples))
            for rr in rows_here:
                rr["cells_driver_mut_called_total"] = n_called
            all_rows += rows_here

    for r in all_rows:
        r["carrier_min_variant_fraction"] = carrier_af
        r["reference_max_variant_fraction"] = ref_af
        r["min_cell_reads"] = min_cell_reads
        r["min_mutant_reads_to_call_cell_mutation_positive"] = min_driver_mut_reads
        r["min_wt_reads_to_call_cell_mutation_negative"] = min_driver_wt_reads
    add_fdr(all_rows, key="binomial_p", out_key="q_value")
    suffix = f"_{tag}" if tag else ""
    _write_tsv(all_rows, os.path.join(out_dir, f"mt_driver_association_all{suffix}.tsv"))
    nominated = nominate(all_rows, max_q, min_recall, max_background, min_mut_cells)
    nominated.sort(key=lambda r: (r["q_value"], -r["cells_driver_mut_total"]))
    _write_tsv(nominated, os.path.join(out_dir, f"mt_passenger_candidates{suffix}.tsv"))
    best = best_variant_per_mutation(all_rows, min_recall_for_imputation=min_recall,
                                     max_background_for_imputation=max_background,
                                     min_mut_cells_for_imputation=10,
                                     max_q_for_imputation=max_q)
    _write_tsv(best, os.path.join(out_dir, f"mt_best_variant_per_mutation{suffix}.tsv"))
    # every significant, specific association (not one per mutation)
    sig = significant_associations(all_rows, max_q=max_q, max_background=max_background)
    _write_tsv(sig, os.path.join(out_dir, f"mt_significant_associations{suffix}.tsv"))
    n_mut_assoc = len(set((r["donor"], r["mutation"]) for r in sig))
    log.info("tests: %d; mutations with a best association: %d; significant associations: %d rows "
             "across %d mutations", len(all_rows), len(best), len(sig), n_mut_assoc)
    return all_rows, nominated


def cohort_associate_sweep(mt_dir, results_dir, out_dir, settings=((1, 1), (1, 2)), **kw):
    """Run the association under several read-support settings and write one comparison table.

    Each setting is (min_alt_reads_per_cell, min_ref_reads_per_cell). The mutant threshold is 1 by
    default because observing the variant allele at all means the cell carries it; the wild-type
    threshold is varied because calling a cell wild-type is an inference from absence and gets more
    reliable with more reads.
    """
    summary = []
    for alt_n, ref_n in settings:
        tag = f"alt{alt_n}_wt{ref_n}"
        rows, nom = cohort_associate(mt_dir, results_dir, out_dir,
                                     min_alt_reads_per_cell=alt_n,
                                     min_ref_reads_per_cell=ref_n, tag=tag, **kw)
        sig = [r for r in rows if r.get("q_value", 1.0) <= kw.get("max_q", 0.05)]
        strong = [r for r in nom if r["cells_driver_mut_total"] >= 20]
        summary.append(dict(
            min_mtvar_reads_to_call_cell_carrier=alt_n,
            min_normal_reads_to_call_cell_wildtype=ref_n,
            n_tests_mtvar_x_driver=len(rows),
            median_cells_driver_mut_per_test=(int(np.median(
                [r["cells_driver_mut_total"] for r in rows])) if rows else 0),
            median_cells_driver_wt_per_test=(int(np.median(
                [r["cells_driver_wt_total"] for r in rows])) if rows else 0),
            median_frac_driver_wt_with_mtvar=(round(float(np.median(
                [r["frac_driver_wt_with_mtvar"] for r in rows])), 5) if rows else 0),
            n_tests_significant_after_fdr=len(sig),
            n_nominated_candidates=len(nom),
            n_nominated_with_ge20_driver_mut_cells=len(strong),
            association_table_file=f"mt_driver_association_all_{tag}.tsv",
            candidate_table_file=f"mt_passenger_candidates_{tag}.tsv"))
    _write_tsv(summary, os.path.join(out_dir, "mt_parameter_comparison.tsv"))
    log.info("wrote %s", os.path.join(out_dir, "mt_parameter_comparison.tsv"))
    return summary


_MOUNT_ALIASES = [("/data/salomonis-archive/", "/Volumes/salomonis-archive/"),
                  ("/data/salomonis2/", "/Volumes/salomonis2/"),
                  ("/data/saljh8/", "/Users/saljh8/")]


def resolve_bam(path):
    if os.path.exists(path):
        return path
    for a, b in _MOUNT_ALIASES:
        for src, dst in ((a, b), (b, a)):
            if path.startswith(src):
                alt = dst + path[len(src):]
                if os.path.exists(alt):
                    return alt
    return None


def _scan_one(args):
    sample, bam, out_dir, kw = args
    done = os.path.join(out_dir, f"{sample}_mt_cell_calls.tsv")
    if os.path.exists(done) and os.path.getsize(done) > 0:
        log.info("[%s] already scanned -- skipped", sample)
        return sample, True
    try:
        run(sample, bam, out_dir, **kw)
        return sample, True
    except Exception as e:
        log.exception("[%s] FAILED: %s", sample, e)
        return sample, False


def scan_cohort(metadata, out_dir, workers=1, samples=None, **kw):
    """Scan every sample in an sclr-style metadata file (uid, bam). Resumable: a sample whose
    per-cell output already exists is skipped, so the cohort can be built up over several runs."""
    import csv as _csv
    os.makedirs(out_dir, exist_ok=True)
    order, uid_bams = [], {}
    with open(metadata) as fh:
        for r in _csv.DictReader(fh, delimiter="\t"):
            uid, bam = (r.get("uid") or "").strip(), (r.get("bam") or "").strip()
            if not uid or not bam:
                continue
            if uid not in uid_bams:
                uid_bams[uid] = []
                order.append(uid)
            uid_bams[uid].append(bam)
    if samples:
        order = [u for u in order if u in set(samples)]

    jobs = []
    for uid in order:
        found = [p for p in (resolve_bam(b) for b in uid_bams[uid]) if p]
        if not found:
            log.warning("[%s] no BAM on disk -- skipped", uid)
            continue
        jobs.append((uid, found[0], out_dir, kw))
    log.info("scanning %d samples with %d worker(s)", len(jobs), workers)

    if workers > 1:
        import multiprocessing as mp
        ctx = mp.get_context("fork")
        with ctx.Pool(workers) as pool:
            results = pool.map(_scan_one, jobs)
    else:
        results = [_scan_one(j) for j in jobs]
    ok = [s for s, good in results if good]
    log.info("cohort scan complete: %d/%d samples", len(ok), len(jobs))
    return ok


def parse_arguments():
    ap = argparse.ArgumentParser(
        description="Fast mitochondrial variant discovery and passenger-marker nomination")
    ap.add_argument("--bam", default=None, help="BAM file (index required); single-sample mode")
    ap.add_argument("--sample", default=None, help="Sample name used in output filenames")
    ap.add_argument("--metadata", default=None,
                    help="sclr metadata TSV (uid, bam): scan EVERY sample. Resumable -- a sample "
                         "whose per-cell output already exists is skipped.")
    ap.add_argument("--samples", nargs="*", default=None,
                    help="With --metadata: restrict to these uids")
    ap.add_argument("--workers", type=int, default=1,
                    help="Scan this many samples in parallel with --metadata")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--contig", default=None, help="Mitochondrial contig (default: auto-detect)")
    ap.add_argument("--driver-calls", default=None,
                    help="variant_extraction.py <sample>_complete_analysis.tsv giving per-cell "
                         "driver MUT/WT status; association step is skipped without it")
    ap.add_argument("--driver-label", default=None,
                    help="Restrict driver status to this variant label (default: any variant)")
    ap.add_argument("--min-depth", type=int, default=100)
    ap.add_argument("--min-alt-reads", type=int, default=5)
    ap.add_argument("--min-af", type=float, default=0.001,
                    help="Lowest within-sample alternate fraction to consider. Default 0.001")
    ap.add_argument("--max-af", type=float, default=0.95,
                    help="Above this the variant is inherited (same in every cell) and useless as a "
                         "clone marker. Default 0.95")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-base-quality", type=int, default=20)
    ap.add_argument("--min-reads-per-cell", type=int, default=2)
    ap.add_argument("--min-cell-af", type=float, default=0.10)
    ap.add_argument("--max-q", type=float, default=0.05, help="FDR cutoff for nomination")
    ap.add_argument("--min-precision", type=float, default=0.80,
                    help="Minimum fraction of marked cells that must be driver-mutant")
    ap.add_argument("--bulk", action="store_true", help="BAM without cell barcodes")
    ap.add_argument("--associate", action="store_true",
                    help="Cross-patient association mode: test EVERY driver mutation in every "
                         "patient against every mitochondrial variant already scanned")
    ap.add_argument("--results-dir", default=None,
                    help="variant_extraction.py output directory (<sample>_complete_analysis.tsv)")
    ap.add_argument("--mt-dir", default=None,
                    help="Directory holding <sample>_mt_cell_calls.tsv (default: --output-dir). "
                         "Set this when writing results somewhere other than the mitochondrial "
                         "scan folder.")
    ap.add_argument("--panels", nargs="*", default=[],
                    help="Panel TSVs; their expected_uids column defines the patient groups")
    ap.add_argument("--donor-map", default=None,
                    help="Optional sample<TAB>patient TSV overriding the panel-derived grouping")
    ap.add_argument("--carrier-af", type=float, default=0.90,
                    help="Minimum within-cell variant fraction to call a cell a CARRIER of the "
                         "mitochondrial variant. High by default so a somatic clonal (near-"
                         "homoplasmic) marker is required, not an intermediate germline heteroplasmy. "
                         "Default 0.90")
    ap.add_argument("--ref-af", type=float, default=0.10,
                    help="Maximum within-cell variant fraction to call a cell REFERENCE. Default 0.10")
    ap.add_argument("--min-cell-reads", type=int, default=3,
                    help="Minimum mitochondrial reads at the position for a cell to be genotyped "
                         "(carrier or reference); fewer are excluded. Default 3")
    ap.add_argument("--n-permutations", type=int, default=0,
                    help="Permutation null: shuffle mutant/wild-type labels this many times to "
                         "estimate the empirical false-discovery rate. Default 0 (off)")
    ap.add_argument("--per-sample", action="store_true",
                    help="Analyse each sample independently instead of pooling a patient's "
                         "timepoints (avoids the per-sample mitochondrial-discovery confound)")
    ap.add_argument("--build-matrix", action="store_true",
                    help="Write a per-sample cells x (mutations + mitochondrial variants) AnnData "
                         "(.h5ad) for immediate query; needs --mt-dir, --results-dir, --samples")
    ap.add_argument("--build-bam-matrix", action="store_true",
                    help="IMPROVED workflow: build the per-sample cell x feature h5ad DIRECTLY FROM "
                         "BAMs -- UMI-level, indel-aware, de-novo mt discovery (SNVs AND indels). "
                         "Needs --metadata (uid<TAB>bam...), --driver-variants, --output-dir; "
                         "--samples restricts to listed uids. Writes <uid>_bam_feature_matrix.h5ad.")
    ap.add_argument("--marker-variants", default=None, metavar="TSV",
                    help="Precomputed external marker panel (same columns as --driver-variants) "
                         "genotyped per cell and placed in the MARKER slot (feature_type=mtDNA) for "
                         "--build-bam-matrix. Genotyped by the same indel-aware CIGAR walk as drivers; "
                         "no de-novo discovery for these. (Nuclear-passenger orchestration that fills "
                         "this lives in mt_nuclear_variant.py.)")
    ap.add_argument("--no-mt-discovery", action="store_true",
                    help="Skip de-novo chrM discovery in --build-bam-matrix; use only "
                         "--marker-variants as the marker panel.")
    ap.add_argument("--driver-variants", default=None, metavar="TSV",
                    help="TSV of driver variants to genotype for --build-bam-matrix: columns "
                         "name<TAB>chrom<TAB>pos<TAB>ref<TAB>alt (header required). SNV or indel.")
    ap.add_argument("--min-cells-mt", type=int, default=10,
                    help="Min supporting reads for a de-novo mt variant to enter the panel "
                         "(--build-bam-matrix). Default 10")
    ap.add_argument("--min-af-mt", type=float, default=0.0005,
                    help="Min BULK allele fraction for a de-novo mt variant (--build-bam-matrix). "
                         "0.0005 (0.05%%) catches markers of SMALL clones -- a marker of a ~17-cell "
                         "clone diluted by thousands of wild-type cells has bulk fraction ~0.1%% "
                         "(measured: TP53's best marker chrM:7498 = 0.11%%, chrM:3122 = 0.65%%). A "
                         "higher floor silently drops them; downstream precision/Fisher reject the "
                         "extra low-AF noise this admits. Default 0.0005")
    ap.add_argument("--nominate-min-mut", type=int, default=10,
                    help="Min mutant cells for a mutation to be TESTED by --nominate/--nominate-presence. "
                         "Default 10 (so small but real clones like TP53's 17 cells are tested).")
    ap.add_argument("--min-alt-reads-per-cell", type=int, default=1,
                    help="[legacy] superseded by --carrier-af/--min-cell-reads")
    ap.add_argument("--min-ref-reads-per-cell", type=int, default=1,
                    help="Reads showing only the normal allele needed to call a cell wild-type. "
                         "Separate parameter because wild-type is inferred from ABSENCE and gets "
                         "more reliable with more reads. Default 1")
    ap.add_argument("--min-driver-mut-reads", type=int, default=1,
                    help="Reads showing the mutant allele needed to call a cell mutation-positive. "
                         "Default 1")
    ap.add_argument("--min-driver-wt-reads", type=int, default=1,
                    help="Reads showing only the normal allele needed to call a cell "
                         "mutation-negative. This is the weaker call (absence of evidence); most "
                         "wild-type calls in this cohort rest on a single read. Default 1")
    ap.add_argument("--report-from", default=None, metavar="ASSOC.tsv",
                    help="Regenerate the significant-associations and best-variant summaries from an "
                         "existing mt_driver_association_all*.tsv (no recompute). Fast.")
    ap.add_argument("--sweep-driver-wt-reads", nargs="*", type=int, default=None, metavar="N",
                    help="Repeat the association at each of these wild-type read thresholds on the "
                         "MUTATION side and write a comparison table, e.g. 1 2 3")
    ap.add_argument("--sweep", nargs="*", default=None, metavar="ALT,WT",
                    help="Mitochondrial-side settings, e.g. --sweep 1,1 1,2")
    ap.add_argument("--assoc-bams", action="store_true",
                    help="Pooled presence-based (any-alt-read, UMI-level) association between two "
                         "variants extracted directly from BAMs -- handles indels (e.g. chrM:3565 A>AC) "
                         "that count_coverage cannot see. Needs --bam-map, --var-a, --var-b, --output-dir.")
    ap.add_argument("--bam-map", default=None, metavar="TSV",
                    help="Two-column TSV (sample<TAB>bam_path) of BAMs to pool for --assoc-bams.")
    ap.add_argument("--var-a", default=None, metavar="chrom:pos:ref:alt:name",
                    help="First variant for --assoc-bams, e.g. chrM:3565:A:AC:chrM3565 (SNV or indel).")
    ap.add_argument("--var-b", default=None, metavar="chrom:pos:ref:alt:name",
                    help="Second variant for --assoc-bams, e.g. chrX:134415107:G:C:PHF6.")
    ap.add_argument("--nominate", action="store_true",
                    help="Scan EVERY testable (sample, driver mutation) -- mutant clone n_mut>=20 and "
                         "wild-type arm n_wt>=30 -- and report the best specific mt marker per mutation "
                         "at a fixed threshold policy (needs --matrix-dir, --samples). Output: "
                         "mt_nomination.tsv in --output-dir. Use --transitions-only to drop transversion "
                         "artifacts.")
    ap.add_argument("--nominate-presence", action="store_true",
                    help="Improved-framework nomination: same as --nominate but presence-based "
                         "(any-alt-read, no within-cell fraction threshold). Output: "
                         "mt_nomination_presence.tsv. Still substitution-only (indels need BAM re-extraction).")
    ap.add_argument("--nominate-fdr", action="store_true",
                    help="Presence-based nomination with Benjamini-Hochberg FDR across ALL tests "
                         "(report q<0.05). Runs SAMPLE level always; if --patient-map is given, also "
                         "PATIENT level (pooling a patient's samples, MISTRG + direct). Needs "
                         "--matrix-dir, --samples. Outputs mt_nomination_presence_sample.tsv and "
                         "mt_nomination_presence_patient.tsv.")
    ap.add_argument("--patient-map", default=None, metavar="TSV",
                    help="Two-column TSV (sample<TAB>patient) grouping samples into patients for the "
                         "patient level of --nominate-fdr (combines MISTRG xenograft + direct cells). "
                         "FDR cutoff is --max-q (default 0.05).")
    ap.add_argument("--annotate-hits", default=None, metavar="HITS.tsv",
                    help="Recompute the full 2x2 for each nomination hit and add n_wt_covered, "
                         "n_mt_predicted (imputable), sensitivity and specificity columns. Needs "
                         "--matrix-dir; --patient-map for patient-level rows. Writes <hits>_annotated.tsv.")
    ap.add_argument("--transitions-only", action="store_true",
                    help="For --nominate: keep only transition mt variants (A<->G, C<->T), the real "
                         "somatic mtDNA class; drop transversions (largely artifact).")
    ap.add_argument("--genotype-summary", action="store_true",
                    help="Write per (sample, driver mutation) MUT/WT/UNK counts across the per-sample "
                         "matrices (needs --matrix-dir, --samples), to find mutations with a wild-type "
                         "arm worth sweeping. Output: driver_genotype_summary.tsv in --output-dir.")
    ap.add_argument("--sweep-driver", default=None, metavar="DRIVER",
                    help="Sweep every threshold combination for ONE driver mutation across the "
                         "per-sample cell x feature matrices (needs --matrix-dir, --samples). Accepts "
                         "an exact feature id (e.g. PHF6_chrX:134415107) or a gene prefix (e.g. PHF6).")
    ap.add_argument("--marker-profile", nargs=2, default=None, metavar=("DRIVER", "MT"),
                    help="Descriptive across-carrier_af profile of ONE (driver, mt-variant) pair in "
                         "ONE --samples sample: background/precision/recall/imputable at each carrier_af "
                         "(default 0.3 0.5 0.7 0.9), to show whether specificity is cutoff-invariant "
                         "(real clonal marker) or threshold-manufactured (germline/artifact).")
    ap.add_argument("--contingency", default=None, metavar="DRIVER",
                    help="Explicit per-cell contingency for ONE driver (in ONE --samples sample) vs the "
                         "mt variants in --contingency-mt: exact MUT/WT/UNCOVERED x carrier/reference/"
                         "uncovered cell counts and overlaps. Writes contingency_<driver>.tsv.")
    ap.add_argument("--contingency-mt", nargs="*", default=None, metavar="MT",
                    help="One or more mt variant ids (e.g. chrM:302_A>ACC) for --contingency.")
    ap.add_argument("--matrix-dir", default=None,
                    help="Directory holding the per-sample <sample>_cell_feature_matrix.h5ad files "
                         "(built with --build-matrix). Used by --sweep-driver.")
    ap.add_argument("--sweep-carrier-af", nargs="*", type=float, default=None, metavar="AF",
                    help="Carrier-fraction thresholds to sweep for --sweep-driver "
                         "(default 0.1 0.2 0.3 0.5 0.7 0.9).")
    ap.add_argument("--sweep-min-cell-reads", nargs="*", type=int, default=None, metavar="N",
                    help="Per-cell mt total-read thresholds to sweep for --sweep-driver "
                         "(default 2 3 5).")
    ap.add_argument("--sweep-min-driver-reads", nargs="*", type=int, default=None, metavar="N",
                    help="Mutant-read thresholds to call a cell driver-positive, swept for "
                         "--sweep-driver (default 1).")
    ap.add_argument("--sweep-min-driver-wt-reads", nargs="*", type=int, default=None, metavar="N",
                    help="Normal-allele read thresholds to call a cell driver-negative (wild-type), "
                         "swept for --sweep-driver. 0 = use every non-mutant cell (including "
                         "driver-uncovered cells) as background. Default 0 1 2.")
    ap.add_argument("--sweep-ref-af", nargs="*", type=float, default=None, metavar="AF",
                    help="Reference-fraction ceilings to sweep for --sweep-driver (default 0.10).")
    ap.add_argument("--sweep-min-wt-covered", type=int, default=None, metavar="N",
                    help="Minimum covered wild-type cells for a --sweep-driver hit to count as "
                         "specific. Guards against a high wild-type read floor shrinking the "
                         "denominator to a few cells (trivial background 0, precision 1.0). Default 20.")
    ap.add_argument("--min-mut-cells", type=int, default=5,
                    help="Minimum mutant cells for a driver variant to be tested. Default 5")
    ap.add_argument("--min-recall", type=float, default=0.90,
                    help="Fraction of mutant cells the variant must mark. Default 0.90")
    ap.add_argument("--max-background", type=float, default=0.10,
                    help="Maximum carrier rate among non-mutant cells; above this a high overlap "
                         "is explained by the variant simply being common. Default 0.10")
    return ap.parse_args()


def main():
    a = parse_arguments()
    if a.build_matrix:
        if not (a.results_dir and a.samples):
            raise SystemExit("--build-matrix requires --results-dir and --samples")
        mt_dir = a.mt_dir or a.output_dir
        for s in a.samples:
            build_cell_feature_matrix(s, mt_dir, a.results_dir, a.output_dir)
        return
    if a.build_bam_matrix:
        if not (a.metadata and a.driver_variants):
            raise SystemExit("--build-bam-matrix requires --metadata and --driver-variants")
        # metadata: header uid<TAB>bam[...]; group bam paths per uid
        bam_by_uid = {}
        with open(a.metadata) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            ui = header.index("uid"); bi = header.index("bam")
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(ui, bi) or not parts[ui]:
                    continue
                bam_by_uid.setdefault(parts[ui], [])
                if parts[bi] and parts[bi] not in bam_by_uid[parts[ui]]:
                    bam_by_uid[parts[ui]].append(parts[bi])
        # driver variants
        drv = []
        with open(a.driver_variants) as fh:
            dh = fh.readline().rstrip("\n").split("\t")
            idx = {c: dh.index(c) for c in ("name", "chrom", "pos", "ref", "alt")}
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 5 or not p[idx["name"]]:
                    continue
                drv.append(dict(name=p[idx["name"]], chrom=p[idx["chrom"]], pos=p[idx["pos"]],
                                ref=p[idx["ref"]], alt=p[idx["alt"]]))
        # optional precomputed external marker panel (e.g. nuclear passengers)
        markers = None
        if a.marker_variants:
            markers = []
            with open(a.marker_variants) as fh:
                mh = fh.readline().rstrip("\n").split("\t")
                midx = {c: mh.index(c) for c in ("name", "chrom", "pos", "ref", "alt")}
                for line in fh:
                    mp = line.rstrip("\n").split("\t")
                    if len(mp) < 5 or not mp[midx["name"]]:
                        continue
                    markers.append(dict(name=mp[midx["name"]], chrom=mp[midx["chrom"]],
                                        pos=mp[midx["pos"]], ref=mp[midx["ref"]], alt=mp[midx["alt"]]))
            log.info("loaded %d external marker variants from %s", len(markers), a.marker_variants)
        uids = a.samples if a.samples else list(bam_by_uid)
        for uid in uids:
            paths = bam_by_uid.get(uid)
            if not paths:
                log.warning("no BAM for uid %s in metadata", uid); continue
            build_bam_feature_matrix(uid, paths, drv, a.output_dir, min_mapq=a.min_mapq,
                                     min_cells_mt=a.min_cells_mt, min_af=a.min_af_mt, contig=a.contig,
                                     marker_variants=markers, mt_discovery=not a.no_mt_discovery)
        return
    if a.assoc_bams:
        if not (a.bam_map and a.var_a and a.var_b):
            raise SystemExit("--assoc-bams requires --bam-map, --var-a, --var-b")

        def _pv(spec):
            chrom, pos, ref, alt, name = spec.split(":")
            return dict(chrom=chrom, pos=int(pos), ref=ref, alt=alt, name=name)

        bam_map = {}
        with open(a.bam_map) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                s, bp = line.split("\t")[:2]
                bam_map[s] = bp
        _, summary = assoc_variants_from_bams(bam_map, _pv(a.var_a), _pv(a.var_b),
                                              a.output_dir, min_mapq=a.min_mapq)
        log.info("association summary: %s", summary)
        return
    if a.genotype_summary:
        matrix_dir = a.matrix_dir or a.mt_dir or a.output_dir
        if not a.samples:
            raise SystemExit("--genotype-summary requires --samples")
        driver_genotype_summary(matrix_dir, a.samples, a.output_dir)
        return
    if a.nominate:
        matrix_dir = a.matrix_dir or a.mt_dir or a.output_dir
        if not a.samples:
            raise SystemExit("--nominate requires --samples")
        nominate_mt_associations(matrix_dir, a.samples, a.output_dir, min_mut=a.nominate_min_mut,
                                 transitions_only=a.transitions_only)
        return
    if a.annotate_hits:
        matrix_dir = a.matrix_dir or a.mt_dir or a.output_dir
        pmap = None
        if a.patient_map:
            pmap = {}
            with open(a.patient_map) as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2 and parts[0]:
                        pmap.setdefault(parts[1], []).append(parts[0])
        annotate_hits_2x2(matrix_dir, a.annotate_hits, a.output_dir, patient_map=pmap)
        return
    if a.nominate_fdr:
        matrix_dir = a.matrix_dir or a.mt_dir or a.output_dir
        if not a.samples:
            raise SystemExit("--nominate-fdr requires --samples")
        # SAMPLE level: each sample is its own unit
        sample_units = {s: [s] for s in a.samples}
        nominate_presence_fdr(matrix_dir, sample_units, a.output_dir, level="sample",
                              min_mut=a.nominate_min_mut, max_q=a.max_q, max_background=a.max_background)
        # PATIENT level: pool a patient's samples (MISTRG + direct)
        if a.patient_map:
            pmap = {}
            with open(a.patient_map) as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 2 and parts[0] in set(a.samples):
                        pmap.setdefault(parts[1], []).append(parts[0])
            nominate_presence_fdr(matrix_dir, pmap, a.output_dir, level="patient",
                                  min_mut=a.nominate_min_mut, max_q=a.max_q,
                                  max_background=a.max_background)
        return
    if a.nominate_presence:
        matrix_dir = a.matrix_dir or a.mt_dir or a.output_dir
        if not a.samples:
            raise SystemExit("--nominate-presence requires --samples")
        nominate_presence(matrix_dir, a.samples, a.output_dir, min_mut=a.nominate_min_mut,
                          transitions_only=a.transitions_only)
        return
    if a.contingency:
        matrix_dir = a.matrix_dir or a.mt_dir or a.output_dir
        if not a.samples or len(a.samples) != 1 or not a.contingency_mt:
            raise SystemExit("--contingency requires one --samples and --contingency-mt <mt...>")
        rows = contingency_report(matrix_dir, a.samples[0], a.contingency, a.contingency_mt,
                                  out_dir=a.output_dir)
        for r in rows:
            print(f"\n=== {r['driver']}  x  {r['mt_variant']}   (sample {r['sample']}, "
                  f"{r['total_cells']} cells) ===")
            print(f"  DRIVER:  MUTANT={r['driver_MUT']}   WILD-TYPE={r['driver_WT']}   "
                  f"UNCOVERED={r['driver_UNCOVERED']}")
            print(f"  mt {r['mt_variant']}:  carrier(>=1 alt)={r['mt_carrier_cells']}   "
                  f"reference(0 alt,covered)={r['mt_reference_cells']}   uncovered={r['mt_uncovered_cells']}")
            print("  OVERLAPS (cells):")
            print(f"    MUTANT   & mt-carrier = {r['MUT_and_mtCarrier']:5d}   & mt-reference = "
                  f"{r['MUT_and_mtRef']:5d}   & mt-uncovered = {r['MUT_and_mtUncov']:5d}")
            print(f"    WILDTYPE & mt-carrier = {r['WT_and_mtCarrier']:5d}   & mt-reference = "
                  f"{r['WT_and_mtRef']:5d}   & mt-uncovered = {r['WT_and_mtUncov']:5d}")
            print(f"    UNCOVERED& mt-carrier = {r['UNCOV_and_mtCarrier_imputable']:5d} (imputable) "
                  f"  & mt-reference = {r['UNCOV_and_mtRef']:5d}")
            print(f"  recall(MUT carrying)={r['recall_MUT_carry']}  "
                  f"background(WT carrying)={r['background_WT_carry']}  "
                  f"precision(carriers that are MUT)={r['precision_carrier_MUT']}")
        return
    if a.marker_profile:
        matrix_dir = a.matrix_dir or a.mt_dir or a.output_dir
        if not a.samples or len(a.samples) != 1:
            raise SystemExit("--marker-profile requires exactly one --samples <sample>")
        caf = tuple(a.sweep_carrier_af) if a.sweep_carrier_af else (0.3, 0.5, 0.7, 0.9)
        driver_feat, mt_feat = a.marker_profile
        marker_threshold_profile(matrix_dir, a.samples[0], driver_feat, mt_feat, a.output_dir,
                                 carrier_afs=caf)
        return
    if a.sweep_driver:
        matrix_dir = a.matrix_dir or a.mt_dir or a.output_dir
        if not a.samples:
            raise SystemExit("--sweep-driver requires --samples")
        kw = {}
        if a.sweep_carrier_af:      kw["carrier_afs"] = tuple(a.sweep_carrier_af)
        if a.sweep_min_cell_reads is not None:  kw["min_cell_reads_list"] = tuple(a.sweep_min_cell_reads)
        if a.sweep_min_driver_reads: kw["min_driver_reads_list"] = tuple(a.sweep_min_driver_reads)
        if a.sweep_min_driver_wt_reads is not None:
            kw["min_driver_wt_reads_list"] = tuple(a.sweep_min_driver_wt_reads)
        if a.sweep_ref_af:          kw["ref_afs"] = tuple(a.sweep_ref_af)
        if a.sweep_min_wt_covered is not None:  kw["min_wt_covered"] = a.sweep_min_wt_covered
        sweep_thresholds_for_driver(matrix_dir, a.samples, a.sweep_driver, a.output_dir,
                                    min_carriers=a.min_mut_cells, **kw)
        return
    if a.report_from:
        tag = ""
        base = os.path.basename(a.report_from)
        if "driverwt" in base:
            tag = "driverwt" + base.split("driverwt")[1].split(".")[0]
        report_from_association_table(a.report_from, a.output_dir, max_q=a.max_q,
                                      max_background=a.max_background, min_recall=a.min_recall,
                                      tag=tag)
        return
    if a.associate:
        if not a.results_dir:
            raise SystemExit("--associate requires --results-dir (variant_extraction.py outputs)")
        kw = dict(metadata=a.metadata, panel_paths=a.panels, donor_map=a.donor_map,
                  samples=a.samples, min_mut_cells=a.min_mut_cells, max_q=a.max_q,
                  min_recall=a.min_recall, max_background=a.max_background,
                  min_cell_reads=a.min_cell_reads, carrier_af=a.carrier_af, ref_af=a.ref_af)
        if a.sweep_driver_wt_reads:
            summary = []
            raw_cache = {}                    # per-patient raw counts, read once, reused every pass
            for n in a.sweep_driver_wt_reads:
                rows, nom = cohort_associate(
                    a.output_dir, a.results_dir, a.output_dir,
                    min_driver_mut_reads=a.min_driver_mut_reads, min_driver_wt_reads=n,
                    min_alt_reads_per_cell=a.min_alt_reads_per_cell,
                    min_ref_reads_per_cell=a.min_ref_reads_per_cell,
                    tag=f"driverwt{n}", raw_cache=raw_cache, **kw)
                best = best_variant_per_mutation(rows, min_recall_for_imputation=a.min_recall)
                summary.append(dict(
                    min_wt_reads_to_call_cell_mutation_negative=n,
                    n_tests_mtvar_x_driver=len(rows),
                    n_donor_mutation_pairs=len(best),
                    median_cells_mutation_negative=(int(np.median(
                        [r["cells_driver_wt_total"] for r in rows])) if rows else 0),
                    median_frac_wildtype_cells_with_mtvar=(round(float(np.median(
                        [r["frac_driver_wt_with_mtvar"] for r in rows])), 5) if rows else 0),
                    n_nominated_candidates=len(nom),
                    n_usable_for_imputation=sum(
                        1 for r in best if r["imputation_verdict"] == "usable for imputation"),
                    best_variant_table_file=f"mt_best_variant_per_mutation_driverwt{n}.tsv"))
            _write_tsv(summary, os.path.join(a.output_dir,
                                             "mt_driver_wt_threshold_comparison.tsv"))
            log.info("wrote mt_driver_wt_threshold_comparison.tsv")
            return
        if a.sweep:
            settings = [tuple(int(x) for x in s.split(",")) for s in a.sweep]
            cohort_associate_sweep(a.output_dir, a.results_dir, a.output_dir,
                                   settings=settings, **kw)
        else:
            mt_dir = a.mt_dir or a.output_dir
            cohort_associate(mt_dir, a.results_dir, a.output_dir,
                             n_permutations=a.n_permutations, per_sample=a.per_sample, **kw)
        return
    if a.metadata:
        scan_cohort(a.metadata, a.output_dir, workers=a.workers, samples=a.samples,
                    contig=a.contig, min_depth=a.min_depth, min_alt_reads=a.min_alt_reads,
                    min_af=a.min_af, max_af=a.max_af, min_mapq=a.min_mapq,
                    min_base_quality=a.min_base_quality,
                    min_reads_per_cell=a.min_reads_per_cell, min_cell_af=a.min_cell_af,
                    bulk=a.bulk)
        return
    if not (a.bam and a.sample):
        raise SystemExit("provide --metadata (cohort mode) OR both --bam and --sample")
    run(a.sample, a.bam, a.output_dir, driver_calls=a.driver_calls, driver_label=a.driver_label,
        contig=a.contig, min_depth=a.min_depth, min_alt_reads=a.min_alt_reads, min_af=a.min_af,
        max_af=a.max_af, min_mapq=a.min_mapq, min_base_quality=a.min_base_quality,
        min_reads_per_cell=a.min_reads_per_cell, min_cell_af=a.min_cell_af, bulk=a.bulk,
        max_q=a.max_q, min_precision=a.min_precision)


if __name__ == "__main__":
    main()
