#!/usr/bin/env python3
"""mt_nuclear_variant.py -- NUCLEAR passenger-variant surrogate discovery for driver mutations.

Companion to mt_variants.py. Where mt_variants uses de-novo MITOCHONDRIAL heteroplasmy as the clonal
passenger markers, this module uses NUCLEAR passenger variants -- the non-germline and high-confidence
SNVs AND indels from the unsupervised genomic discovery -- as the markers, and asks, per driver
mutation, which passenger best marks the mutant clone so it can impute driver status in cells where
the driver gene was not captured.

It REUSES the mt_variants engine unchanged, so the statistics are IDENTICAL to the mitochondrial
framework (see mt_variants/METHODS_mt_variant_association.md):
  * genotyping + matrix build : mt_variants.build_bam_feature_matrix(marker_variants=..., mt_discovery=False)
  * testable-pair summary      : mt_variants.driver_genotype_summary
  * nomination (Fisher / BH-FDR + specificity gate, sample & patient) : mt_variants.nominate_presence_fdr
The only new code here is (1) discovering the nuclear passenger panel from the discovery output TSVs
and (2) the two-stage CLI.

CLI
---
  # STAGE 1 (per uid; compute nodes): build cell x (driver + passenger) matrices
  python3 -m altanalyze3.components.bam.mt_nuclear_variant --build-passenger-matrix \
      --metadata META --driver-variants DV --discovery-manifest DISC_MANIFEST \
      --samples UID --output-dir OUTDIR

  # STAGE 2 (after all uids finish): nominate the best passenger surrogate per driver
  python3 -m altanalyze3.components.bam.mt_nuclear_variant --nominate \
      --matrix-dir OUTDIR --output-dir OUTDIR --patient-map PMAP --samples UID1 UID2 ...
"""
import argparse
import csv
import logging
import os
from glob import glob

from altanalyze3.components.bam import mt_variants as mtv

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("mt_nuclear_variant")

DEFAULT_SUBPATHS = ("results/merged/high_confidence_final.tsv",     # >=3 callers
                    "results/merged/medium_confidence_final.tsv",   # 2 callers
                    "results/merged/low_confidence_final.tsv",      # 1 caller
                    "annotation/gatk/non_germline_final.tsv")       # GATK-only somatic


def discovery_outdirs_by_uid(manifest_path):
    """uid -> [genome_variant_detection outdirs] from a genomic-discovery manifest (cols incl uid,
    outdir). Multiple discovery samples (BAMs) per uid are collected so their passengers UNION."""
    by_uid = {}
    with open(manifest_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ui = header.index("uid"); oi = header.index("outdir")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(ui, oi) or not p[ui]:
                continue
            by_uid.setdefault(p[ui], [])
            if p[oi] and p[oi] not in by_uid[p[ui]]:
                by_uid[p[ui]].append(p[oi])
    return by_uid


def discover_passenger_variants(outdirs, exclude_keys=(), subpaths=DEFAULT_SUBPATHS):
    """Build a nuclear passenger marker panel from precomputed genomic-discovery output TSVs -- the
    nuclear analog of mt_variants.discover_mt_variants (chrM heteroplasmy from a BAM).

    Reads each uid outdir's high-confidence and GATK non-germline final tables (columns
    CHROM POS REF ALT GENE ...), UNIONs their variants (SNV AND indel), drops any whose
    CHROM:POS:REF:ALT is in exclude_keys (the driver panel), and returns a deduplicated list of
    marker dicts {name=GENE_chrom:pos, chrom, pos, ref, alt}.
    """
    exclude = set(exclude_keys)
    seen = set(); markers = []
    for od in outdirs:
        for sp in subpaths:
            path = os.path.join(od, sp)
            if not os.path.isfile(path):
                continue
            with open(path) as fh:
                header = fh.readline().rstrip("\n").split("\t")
                try:
                    ci = {c: header.index(c) for c in ("CHROM", "POS", "REF", "ALT", "GENE")}
                except ValueError:
                    log.warning("passenger TSV missing expected columns, skipped: %s", path); continue
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) <= max(ci.values()):
                        continue
                    chrom, pos, ref, alt, gene = (p[ci["CHROM"]], p[ci["POS"]], p[ci["REF"]],
                                                  p[ci["ALT"]], p[ci["GENE"]])
                    key = f"{chrom}:{pos}:{ref}:{alt}"
                    if key in exclude or key in seen:
                        continue
                    seen.add(key)
                    markers.append(dict(name=f"{gene}_{chrom}:{pos}", chrom=chrom, pos=pos,
                                        ref=ref, alt=alt))
    return markers


def _read_driver_panel(path):
    """Read a name<TAB>chrom<TAB>pos<TAB>ref<TAB>alt driver TSV -> (list of dicts, set of coord keys)."""
    drv = []; keys = set()
    with open(path) as fh:
        dh = fh.readline().rstrip("\n").split("\t")
        idx = {c: dh.index(c) for c in ("name", "chrom", "pos", "ref", "alt")}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5 or not p[idx["name"]]:
                continue
            drv.append(dict(name=p[idx["name"]], chrom=p[idx["chrom"]], pos=p[idx["pos"]],
                            ref=p[idx["ref"]], alt=p[idx["alt"]]))
            keys.add(f'{p[idx["chrom"]]}:{p[idx["pos"]]}:{p[idx["ref"]]}:{p[idx["alt"]]}')
    return drv, keys


def _bam_paths_by_uid(metadata_path):
    """uid -> [bam paths] from the metadata (multiple BAMs per uid pooled)."""
    by_uid = {}
    with open(metadata_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ui = header.index("uid"); bi = header.index("bam")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(ui, bi) or not p[ui]:
                continue
            by_uid.setdefault(p[ui], [])
            if p[bi] and p[bi] not in by_uid[p[ui]]:
                by_uid[p[ui]].append(p[bi])
    return by_uid


def build_passenger_matrix(uid, bam_paths, driver_variants, driver_keys, outdirs, out_dir,
                           subpaths=DEFAULT_SUBPATHS, min_mapq=20, contig=None):
    """Discover this uid's nuclear passenger panel, write it for provenance, then build the cell x
    (driver + passenger) matrix by reusing mt_variants.build_bam_feature_matrix with the passenger
    panel in the marker slot and NO chrM discovery. Returns the h5ad path."""
    markers = discover_passenger_variants(outdirs, exclude_keys=driver_keys, subpaths=subpaths)
    panels_dir = os.path.join(out_dir, "panels"); os.makedirs(panels_dir, exist_ok=True)
    panel_path = os.path.join(panels_dir, f"{uid}_passenger_panel.tsv")
    with open(panel_path, "w") as ph:
        ph.write("name\tchrom\tpos\tref\talt\n")
        for m in markers:
            ph.write(f"{m['name']}\t{m['chrom']}\t{m['pos']}\t{m['ref']}\t{m['alt']}\n")
    n_ind = sum(1 for m in markers if len(m["ref"]) != 1 or len(m["alt"]) != 1)
    log.info("[%s] passenger panel: %d variants (%d indel, %d snv) -> %s", uid, len(markers),
             n_ind, len(markers) - n_ind, panel_path)
    return mtv.build_bam_feature_matrix(uid, bam_paths, driver_variants, out_dir, min_mapq=min_mapq,
                                        contig=contig, marker_variants=markers, mt_discovery=False)


def _panel_ref_alt(panels_dir):
    """marker name (GENE_chrom:pos) -> (ref, alt) from the per-uid passenger panels, so type and the
    full variant can be recovered (the marker name alone does not encode ref/alt)."""
    ra = {}
    for p in glob(os.path.join(panels_dir, "*_passenger_panel.tsv")):
        with open(p) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                ra.setdefault(row["name"], (row["ref"], row["alt"]))
    return ra


def write_passenger_hits_table(annotated_paths, panels_dir, out_path):
    """Emit the best passenger->driver associations in the SAME layout as mt_marker_hits_ALL.tsv but
    with the marker columns renamed MT->Var (#mut-Var, #wt-Var, #Var-pred) and the CORRECT snv/indel
    `type` recovered from the panel ref/alt. Columns: level, gene, unit, n_mut, marker, type, recall,
    sens, spec, prec, bg, #mut-Var, #wt-Var, #wt, #Var-pred, OR, p; sorted by recall descending.
    Sample and patient rows are combined (both levels present in annotated_paths)."""
    ra = _panel_ref_alt(panels_dir)
    cols = ["level", "gene", "driver", "unit", "n_mut", "marker", "type", "recall", "sens", "spec",
            "prec", "bg", "#mut-Var", "#wt-Var", "#wt", "#Var-pred", "OR", "p"]
    rows = []
    for path in annotated_paths:
        if not os.path.isfile(path):
            continue
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                name = r.get("mt_feature", "")
                ref, alt = ra.get(name, ("", ""))
                typ = "indel" if (ref and alt and len(ref) != len(alt)) else "snv"
                marker = f"{name}_{ref}>{alt}" if ref else name
                rows.append({"level": r.get("level", ""), "gene": r.get("gene", ""),
                             "driver": r.get("driver", ""),
                             "unit": r.get("unit", ""), "n_mut": r.get("n_mut", ""), "marker": marker,
                             "type": typ, "recall": r.get("recall_all", ""),
                             "sens": r.get("sensitivity", ""), "spec": r.get("specificity", ""),
                             "prec": r.get("precision", ""), "bg": r.get("background", ""),
                             "#mut-Var": r.get("mutant_carriers", ""), "#wt-Var": r.get("wt_carriers", ""),
                             "#wt": r.get("n_wt_covered", ""), "#Var-pred": r.get("n_mt_predicted", ""),
                             "OR": r.get("odds_ratio", ""), "p": r.get("fisher_p", "")})

    def _rk(x):
        try:
            return -float(x["recall"] or 0)
        except Exception:
            return 0.0
    rows.sort(key=_rk)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    log.info("wrote %s (%d hits)", out_path, len(rows))
    return rows


def driver_passenger_full_scan(adata, driver_feat, min_carriers=3, min_wt_covered=20,
                               min_mut_covered=10):
    """For ONE driver coordinate, compute the full 2x2 for EVERY marker in the matrix -- NO specificity
    or FDR gate, so surrogates can be selected by sensitivity -- but REQUIRE a valid comparison arm:
    n_wt_covered >= min_wt_covered AND n_mut_covered >= min_mut_covered, else sensitivity/specificity
    are undefined (a passenger with no covered wild-type cells shows a meaningless sens=1/spec=0/p=1).
    presence = alt>=1; a=mut&carrier, b=wt&carrier, c=mut&reference, d=wt&reference. Returns rows
    (marker, n_mut, recall, sensitivity=a/(a+c), specificity=d/(b+d), precision=a/(a+b),
    background=b/(b+d), #mut-Var, #wt-Var, #wt, #Var-pred, OR, p) sorted by sensitivity."""
    import numpy as _np
    from scipy.stats import fisher_exact
    j = adata.var_names.get_loc(driver_feat)
    da = _np.asarray(adata.layers["alt"][:, j].todense()).ravel()
    dt = _np.asarray(adata.layers["total"][:, j].todense()).ravel()
    mut = (da >= 1).astype(float); wt = ((da == 0) & (dt >= 1)).astype(float)
    unk = ((da == 0) & (dt == 0)).astype(float)
    n_mut = int(mut.sum())
    names, A, T, _ = mtv.mt_dense_block(adata)
    present = (A >= 1).astype(float); absent = ((A == 0) & (T >= 1)).astype(float)
    a = mut @ present; b = wt @ present; c = mut @ absent; d = wt @ absent
    imp = unk @ present
    rows = []
    for i in range(len(names)):
        ai, bi, ci, di = int(a[i]), int(b[i]), int(c[i]), int(d[i])
        nmut_cov = ai + ci; nwt_cov = bi + di
        if ai < min_carriers or nwt_cov < min_wt_covered or nmut_cov < min_mut_covered:
            continue
        try:
            ORv, pv = fisher_exact([[ai, bi], [ci, di]])
        except Exception:
            ORv, pv = float("nan"), 1.0
        rows.append(dict(marker=str(names[i]), n_mut=n_mut,
                         recall=round(ai / n_mut, 4) if n_mut else 0.0,
                         sensitivity=round(ai / nmut_cov, 4) if nmut_cov else 0.0,
                         specificity=round(di / nwt_cov, 4) if nwt_cov else 0.0,
                         precision=round(ai / (ai + bi), 4) if (ai + bi) else 0.0,
                         background=round(bi / nwt_cov, 4) if nwt_cov else 0.0,
                         mut_Var=ai, wt_Var=bi, n_wt=nwt_cov, Var_pred=int(imp[i]),
                         OR=(round(float(ORv), 3) if ORv == ORv else "inf"), p=float(pv)))
    rows.sort(key=lambda r: (-r["sensitivity"], -r["precision"]))
    return rows


def nominate_all_significant(matrix_dir, units, out_dir, level="sample", min_mut=10, min_wt=30,
                             min_carriers=5, min_wt_covered=20, max_q=0.05, max_background=0.15):
    """Like mt_variants.nominate_presence_fdr but emits EVERY significant (unit x driver x marker)
    association (q<max_q AND background<=max_background), NOT just the best per (unit,driver). Reuses
    mt_variants.scan_mt_presence / mt_dense_block / _load_unit_adata. Writes
    passenger_all_significant_<level>.tsv; returns the rows."""
    import numpy as _np
    all_p = []; shortlist = {}
    for label, samples in units.items():
        adata = mtv._load_unit_adata(matrix_dir, samples)
        if adata is None:
            continue
        block = mtv.mt_dense_block(adata)
        drv = adata.var_names[adata.var["feature_type"].values == "driver"]
        for feat in drv:
            j = adata.var_names.get_loc(feat)
            da = _np.asarray(adata.layers["alt"][:, j].todense()).ravel()
            dt = _np.asarray(adata.layers["total"][:, j].todense()).ravel()
            n_mut = int((da >= 1).sum()); n_wt = int(((da == 0) & (dt >= 1)).sum())
            if n_mut < min_mut or n_wt < min_wt:
                continue
            cands = mtv.scan_mt_presence(adata, feat, min_carriers=min_carriers,
                                         min_wt_covered=min_wt_covered, block=block, return_all=True)
            lst = shortlist.setdefault((label, feat), [])
            for r in cands:
                if r["n_wt_covered"] < min_wt_covered:
                    continue
                all_p.append(r["fisher_p"])
                lst.append(dict(level=level, unit=label, gene=feat.split("_")[0], driver=feat,
                                n_samples=len(samples), n_mut=n_mut, n_wt=n_wt, **r))
    if not all_p:
        log.warning("nominate_all_significant level=%s: no tests", level); return []
    p = _np.sort(_np.asarray(all_p, dtype=float)); m = len(p)
    q_sorted = _np.minimum.accumulate((p * m / _np.arange(1, m + 1))[::-1])[::-1]

    def _p2q(pv):
        i = min(int(_np.searchsorted(p, pv, side="left")), m - 1)
        return float(q_sorted[i])
    out_rows = []
    for lst in shortlist.values():
        for r in lst:
            q = _p2q(r["fisher_p"])
            if q < max_q and r["background"] <= max_background:
                r["q_value"] = round(q, 6); out_rows.append(r)
    out_rows.sort(key=lambda r: (-r["recall_all"], r["q_value"]))
    out = os.path.join(out_dir, f"passenger_all_significant_{level}.tsv")
    mtv._write_tsv(out_rows, out)
    log.info("level=%s: %d significant associations -> %s", level, len(out_rows), out)
    return out_rows


def all_significant_passenger_table(matrix_dir, samples, out_dir, panels_dir, patient_map=None,
                                    min_mut=10, max_q=0.05, max_background=0.15):
    """EVERYTHING: every significant passenger->driver association (uncollapsed), sample + patient,
    annotated, in the MT->Var layout -> passenger_all_significant_marker_hits.tsv."""
    nominate_all_significant(matrix_dir, {s: [s] for s in samples}, out_dir, level="sample",
                             min_mut=min_mut, max_q=max_q, max_background=max_background)
    pmap = None
    if patient_map:
        smp = set(samples); pmap = {}
        with open(patient_map) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2 and parts[0] in smp:
                    pmap.setdefault(parts[1], []).append(parts[0])
        nominate_all_significant(matrix_dir, pmap, out_dir, level="patient",
                                 min_mut=min_mut, max_q=max_q, max_background=max_background)
    annotated = []
    for lvl in ("sample", "patient"):
        src = os.path.join(out_dir, f"passenger_all_significant_{lvl}.tsv")
        if os.path.isfile(src):
            mtv.annotate_hits_2x2(matrix_dir, src, out_dir,
                                  patient_map=(pmap if lvl == "patient" else None))
            ann = os.path.join(out_dir, f"passenger_all_significant_{lvl}_annotated.tsv")
            if os.path.isfile(ann):
                annotated.append(ann)
    return write_passenger_hits_table(annotated, panels_dir,
                                      os.path.join(out_dir, "passenger_all_significant_marker_hits.tsv"))


def nominate_passenger_surrogates(matrix_dir, samples, out_dir, patient_map=None,
                                  min_mut=10, max_q=0.05, max_background=0.15):
    """Run the mt_variants engine on the passenger matrices END-TO-END -- the SAME four stages as the
    mitochondrial framework, all reused unchanged: driver_genotype_summary + nominate_presence_fdr
    (sample + patient) + annotate_hits_2x2 (adds sensitivity / specificity / n_mt_predicted) +
    write_hits_html. Then rename mt_* -> passenger_* and write the combined passenger_marker_hits_ALL.tsv
    (the nuclear analog of mt_marker_hits_ALL.tsv). Identical statistics; only the marker panel differs."""
    mtv.driver_genotype_summary(matrix_dir, samples, out_dir)
    sample_units = {s: [s] for s in samples}
    mtv.nominate_presence_fdr(matrix_dir, sample_units, out_dir, level="sample",
                              min_mut=min_mut, max_q=max_q, max_background=max_background)
    pmap = None
    if patient_map:
        smp = set(samples); pmap = {}
        with open(patient_map) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2 and parts[0] in smp:
                    pmap.setdefault(parts[1], []).append(parts[0])
        mtv.nominate_presence_fdr(matrix_dir, pmap, out_dir, level="patient",
                                  min_mut=min_mut, max_q=max_q, max_background=max_background)
    annotated = []
    for lvl in ("sample", "patient"):
        src = os.path.join(out_dir, f"mt_nomination_presence_{lvl}.tsv")
        if not os.path.isfile(src):
            continue
        dst = os.path.join(out_dir, f"passenger_nomination_presence_{lvl}.tsv")
        os.replace(src, dst)
        log.info("wrote %s", dst)
        # STAGE 3: annotate (recompute the 2x2 -> sensitivity / specificity / n_mt_predicted); patient
        # rows pool via pmap. annotate_hits_2x2 writes <hits>_annotated.tsv into out_dir.
        mtv.annotate_hits_2x2(matrix_dir, dst, out_dir, patient_map=(pmap if lvl == "patient" else None))
        ann = os.path.join(out_dir, f"passenger_nomination_presence_{lvl}_annotated.tsv")
        if os.path.isfile(ann):
            annotated.append(ann)
    if annotated:
        # STAGE 4: HTML report + combined flat table
        mtv.write_hits_html(annotated, os.path.join(out_dir, "passenger_marker_hits.html"),
                            title="Nuclear passenger markers of driver mutations")
        # passenger flat table: mt layout, but marker columns renamed MT->Var and the correct
        # snv/indel type recovered from the panel ref/alt (the marker name does not encode it).
        write_passenger_hits_table(annotated, os.path.join(out_dir, "panels"),
                                   os.path.join(out_dir, "passenger_marker_hits_ALL.tsv"))
    # Also emit the COMPLETE significant table: every driver x marker association (mt AND nuclear),
    # uncollapsed and NOT gated on specificity, so no significant hit is hidden by the best-per-driver
    # nomination or by the surrogate specificity gate. The bg (background) column is kept so specificity
    # can be filtered afterward. This is the single all-markers view; mt rows are just the chrM markers.
    all_significant_passenger_table(matrix_dir, samples, out_dir, os.path.join(out_dir, "panels"),
                                    patient_map=patient_map, min_mut=min_mut, max_q=max_q,
                                    max_background=1.0)


def _read_discovery_variants_with_truth(outdirs, subpaths=DEFAULT_SUBPATHS):
    """Every discovery variant with the callers' own BAM evidence, for detection validation.
    Returns list of dict(chrom,pos,ref,alt,callers,caller_count,bam_alt,bam_dp,vtype). vtype in
    {snv, ins, del, mnv, complex}. Deduplicated by (chrom,pos,ref,alt) across the given TSVs."""
    def _num(x):
        try:
            return float(x)
        except Exception:
            return None
    seen = set()
    out = []
    for od in outdirs:
        for sp in subpaths:
            p = os.path.join(od, sp)
            if not os.path.isfile(p):
                continue
            with open(p) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    ch = row.get("CHROM"); pos = row.get("POS")
                    ref = row.get("REF"); alt = row.get("ALT")
                    if not (ch and pos and ref and alt):
                        continue
                    k = (ch, pos, ref, alt)
                    if k in seen:
                        continue
                    seen.add(k)
                    lr, la = len(ref), len(alt)
                    if lr == 1 and la == 1:
                        vt = "snv"
                    elif la > lr:
                        vt = "ins"
                    elif lr > la:
                        vt = "del"
                    elif lr == la:
                        vt = "mnv"
                    else:
                        vt = "complex"
                    try:
                        posi = int(pos)
                    except Exception:
                        continue
                    out.append(dict(chrom=ch, pos=posi, ref=ref, alt=alt,
                                    callers=row.get("CALLERS", ""), caller_count=row.get("CALLER_COUNT", ""),
                                    bam_alt=_num(row.get("BAM_ALT")), bam_dp=_num(row.get("BAM_DP")),
                                    vtype=vt))
    return out


def validate_detection(bam_path, variants, min_mapq=20):
    """Compare the genotyper's detection against the callers' BAM_ALT, using the SAME detection code the
    per-cell matrices use (mt_variants.genotype_variant_cells), with a TARGETED fetch at each variant
    (conventional BAM seek -- a few bp window per variant, NOT a whole-contig scan) and pooled over all
    cell barcodes (bulk). detected = pooled alt reads > 0. Returns per-variant result dicts."""
    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb")
    results = []
    for v in variants:
        cc = mtv._resolve_contig(bam, v["chrom"])
        g = mtv.genotype_variant_cells(bam, cc, int(v["pos"]), v["ref"], v["alt"],
                                       min_mapq=min_mapq, bulk=True)
        ga = sum(x.get("alt", 0) for x in g.values())
        gt = sum(x.get("total", 0) for x in g.values())
        results.append(dict(chrom=v["chrom"], pos=v["pos"], ref=v["ref"], alt=v["alt"],
                            vtype=v["vtype"], callers=v["callers"], caller_count=v["caller_count"],
                            bam_alt=v["bam_alt"], bam_dp=v["bam_dp"],
                            geno_alt=int(ga), geno_total=int(gt), detected=int(ga > 0)))
    bam.close()
    return results


def run_detection_validation(metadata, discovery_manifest, samples, out_dir, subpaths=DEFAULT_SUBPATHS,
                             min_mapq=20):
    """Multi-BAM detection validation: for each uid, genotype every discovery variant and compare to the
    callers' BAM_ALT. Writes per-variant TSV per sample + a pooled per-stratum (vtype x caller_count)
    detection summary. Detection rate is measured ONLY on truly-present variants (BAM_ALT > 0)."""
    import collections
    os.makedirs(out_dir, exist_ok=True)
    bam_by_uid = _bam_paths_by_uid(metadata)
    disc = discovery_outdirs_by_uid(discovery_manifest)
    pooled = collections.defaultdict(lambda: [0, 0])   # (vtype, caller_count) -> [detected, present]
    per_sample = {}
    for uid in samples:
        paths = bam_by_uid.get(uid)
        od = disc.get(uid, [])
        if not paths or not od:
            log.warning("skip %s (no bam or no discovery outdir)", uid); continue
        variants = _read_discovery_variants_with_truth(od, subpaths)
        res = validate_detection(paths[0], variants, min_mapq=min_mapq)
        outp = os.path.join(out_dir, f"{uid}_detection_validation.tsv")
        cols = ["chrom", "pos", "ref", "alt", "vtype", "caller_count", "callers",
                "bam_alt", "bam_dp", "geno_alt", "geno_total", "detected"]
        with open(outp, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            for r in res:
                w.writerow({k: r.get(k) for k in cols})
        det = present = 0
        for r in res:
            if r["bam_alt"] and r["bam_alt"] > 0:
                present += 1; det += r["detected"]
                pooled[(r["vtype"], r["caller_count"])][0] += r["detected"]
                pooled[(r["vtype"], r["caller_count"])][1] += 1
        per_sample[uid] = (det, present)
        log.info("[%s] detection on BAM_ALT>0 variants: %d/%d = %.4f -> %s",
                 uid, det, present, (det / present if present else 0.0), outp)
    # pooled summary
    sump = os.path.join(out_dir, "detection_summary.tsv")
    with open(sump, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["stratum_vtype", "caller_count", "detected", "present", "detection_rate"])
        by_vtype = collections.defaultdict(lambda: [0, 0])
        for (vt, cc), (d, n) in sorted(pooled.items()):
            w.writerow([vt, cc, d, n, round(d / n, 4) if n else ""])
            by_vtype[vt][0] += d; by_vtype[vt][1] += n
        w.writerow([])
        for vt, (d, n) in sorted(by_vtype.items()):
            w.writerow([vt, "ALL", d, n, round(d / n, 4) if n else ""])
        alld = sum(v[0] for v in by_vtype.values()); alln = sum(v[1] for v in by_vtype.values())
        w.writerow(["ALL", "ALL", alld, alln, round(alld / alln, 4) if alln else ""])
    log.info("wrote pooled detection summary -> %s (overall %d/%d = %.4f)",
             sump, alld, alln, (alld / alln if alln else 0.0))
    return sump


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-passenger-matrix", action="store_true",
                    help="STAGE 1: per uid, build cell x (driver + nuclear-passenger) matrix.")
    ap.add_argument("--nominate", action="store_true",
                    help="STAGE 2: nominate the best passenger surrogate per driver (sample + patient).")
    ap.add_argument("--validate-detection", action="store_true",
                    help="Validate that the genotyper detects the callers' variants (multi-BAM, all "
                         "variants): per-variant TSV + per-stratum detection_summary.tsv vs BAM_ALT.")
    ap.add_argument("--metadata", default=None, help="uid<TAB>bam[...] table (pools BAMs per uid).")
    ap.add_argument("--driver-variants", default=None, help="name chrom pos ref alt driver panel TSV.")
    ap.add_argument("--discovery-manifest", default=None,
                    help="Genomic-discovery manifest (cols incl uid, outdir).")
    ap.add_argument("--passenger-subpaths", nargs="*", default=None,
                    help="Per-outdir passenger TSV subpaths (default: high_confidence_final + non_germline_final).")
    ap.add_argument("--samples", nargs="*", default=None, help="uids to process.")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--matrix-dir", default=None, help="Where the h5ads are (--nominate; default output-dir).")
    ap.add_argument("--patient-map", default=None, help="sample<TAB>patient TSV for the patient level.")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--contig", default=None)
    ap.add_argument("--nominate-min-mut", type=int, default=10)
    ap.add_argument("--max-q", type=float, default=0.05)
    ap.add_argument("--max-background", type=float, default=0.15)
    a = ap.parse_args()

    subpaths = tuple(a.passenger_subpaths) if a.passenger_subpaths else DEFAULT_SUBPATHS

    if a.build_passenger_matrix:
        if not (a.metadata and a.driver_variants and a.discovery_manifest):
            raise SystemExit("--build-passenger-matrix requires --metadata, --driver-variants, --discovery-manifest")
        bam_by_uid = _bam_paths_by_uid(a.metadata)
        drv, driver_keys = _read_driver_panel(a.driver_variants)
        disc = discovery_outdirs_by_uid(a.discovery_manifest)
        uids = a.samples if a.samples else list(bam_by_uid)
        for uid in uids:
            paths = bam_by_uid.get(uid)
            if not paths:
                log.warning("no BAM for uid %s in metadata", uid); continue
            outdirs = disc.get(uid, [])
            if not outdirs:
                log.warning("no discovery outdir for uid %s in %s", uid, a.discovery_manifest); continue
            build_passenger_matrix(uid, paths, drv, driver_keys, outdirs, a.output_dir,
                                   subpaths=subpaths, min_mapq=a.min_mapq, contig=a.contig)
        return

    if a.nominate:
        matrix_dir = a.matrix_dir or a.output_dir
        if not a.samples:
            raise SystemExit("--nominate requires --samples")
        nominate_passenger_surrogates(matrix_dir, a.samples, a.output_dir, patient_map=a.patient_map,
                                      min_mut=a.nominate_min_mut, max_q=a.max_q,
                                      max_background=a.max_background)
        return

    if a.validate_detection:
        if not (a.metadata and a.discovery_manifest and a.samples):
            raise SystemExit("--validate-detection requires --metadata, --discovery-manifest, --samples")
        run_detection_validation(a.metadata, a.discovery_manifest, a.samples, a.output_dir,
                                 subpaths=subpaths, min_mapq=a.min_mapq)
        return

    raise SystemExit("choose one action: --build-passenger-matrix, --nominate, or --validate-detection")


if __name__ == "__main__":
    main()
