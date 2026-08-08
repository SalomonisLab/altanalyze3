#!/usr/bin/env python3
"""altanalyze3 LOH pileup -- build the per-cell chr x germline-het-SNV ALT/REF matrices that feed
loh_classifier, for one multi-timepoint patient. Faithful port of Xuan's S1+S2+S3
(2026-03-05_RUNX1_LOH_classifier), generalized to ANY chromosome and driven by the GATK-germline
calls produced per sample by the discovery pipeline.

Stages (all in one patient run):
  S1  union het-SNVs: from EACH timepoint's GATK germline VCF, keep chrom biallelic SNVs with GT 0/1,
      DP >= min_dp, GQ >= min_gq, AD-VAF in [vaf_lo, vaf_hi]; take the UNION of positions across timepoints.
  S2  per-cell pileup: for each timepoint BAM, pileup every union position (CB barcode, XM UMI, MAPQ>=20),
      dedup per (barcode, UMI), count ALT/REF UMIs per (barcode, SNV).  MAS-seq reads have no base
      qualities, so base-quality filtering is a no-op (matches the source pipeline).
  S3  matrices: concatenate timepoints -> combined_alt.txt / combined_ref.txt (cells x SNV, cell =
      'timepoint::barcode') + cell_metadata.tsv (cell, sample, stage). Append the anchor + driver
      per-cell genotypes (WT->REF, MUT->ALT) from results_variant_extraction so loh_classifier can
      select SNVs against the driver trajectory and label cells by the anchor clone.

Run: python -m altanalyze3.components.bam.loh_pileup --patient 5801 --chrom chr21 \
        --timepoints preHMA:5801M_pre postHMA:5801M_post sAML:5801A_pre,5801A_pre_2 ... \
        --manifest MANIFEST --extraction-dir results_variant_extraction --drivers RUNX1p.W279* SRSF2p.P95R \
        --output OUT
"""
import argparse
import csv
import gzip
import os
import numpy as np
import pandas as pd


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _rc(barcode):
    """Reverse-complement a cell barcode to cellHarmony/annotation orientation.

    variant_extraction writes per-cell driver genotypes with barcodes in this orientation
    (isoform-BAM CB -> junction/h5ad), so the het-SNV pileup barcodes MUST be reverse-complemented
    to match; otherwise the SNV features and the driver/anchor genotypes land on disjoint rows and
    every labeled cell has all-zero SNV features (the classifier then predicts the class prior at
    chance AUC). Any 10x '-N' gem suffix is dropped to match the bare-16mer keys."""
    seq = barcode.split("-", 1)[0]
    return seq.translate(_COMPLEMENT)[::-1]


# ── S1: union of germline het-SNVs on a chromosome ─────────────────────────────
def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def extract_het_union(germline_vcfs, chrom, min_dp=10, min_gq=20, vaf_lo=0.25, vaf_hi=0.75):
    """UNION of chrom biallelic het SNV positions across the given germline VCFs. Returns a DataFrame
    (chrom,pos,ref,alt,snv_id,n_samples,max_dp)."""
    reg = {}
    for vcf in germline_vcfs:
        if not os.path.exists(vcf):
            continue
        with _open(vcf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if f[0] != chrom or "," in f[4] or len(f[3]) != 1 or len(f[4]) != 1:
                    continue
                fmt = dict(zip(f[8].split(":"), f[9].split(":")))
                gt = fmt.get("GT", "./.").replace("|", "/").split("/")
                if set(gt) != {"0", "1"}:
                    continue
                try:
                    dp = int(fmt.get("DP", "0")); gq = int(fmt.get("GQ", "0"))
                except ValueError:
                    continue
                if dp < min_dp or gq < min_gq:
                    continue
                ad = fmt.get("AD", "")
                if "," not in ad:
                    continue
                try:
                    rc, ac = int(ad.split(",")[0]), int(ad.split(",")[1])
                except (ValueError, IndexError):
                    continue
                tot = rc + ac
                if tot == 0 or not (vaf_lo <= ac / tot <= vaf_hi):
                    continue
                key = (f[0], int(f[1]), f[3], f[4])
                r = reg.setdefault(key, {"n": 0, "max_dp": 0})
                r["n"] += 1; r["max_dp"] = max(r["max_dp"], tot)
    rows = [{"chrom": k[0], "pos": k[1], "ref": k[2], "alt": k[3],
             "snv_id": f"{k[0]}:{k[1]}:{k[2]}>{k[3]}", "n_samples": v["n"], "max_dp": v["max_dp"]}
            for k, v in reg.items()]
    return pd.DataFrame(sorted(rows, key=lambda r: r["pos"]))


# ── S2: per-cell pileup at the union positions ─────────────────────────────────
def _umi(read):
    for t in ("XM", "ZU", "YM", "UB"):
        if read.has_tag(t):
            return read.get_tag(t)
    return read.query_name


def pileup_bam(bam_path, cands, min_mapq=20, max_depth=100000):
    """{(barcode, snv_id): [ref_umi, alt_umi]} from one BAM; dedup per (barcode, UMI)."""
    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb")
    refset = set(bam.references)
    out = {}
    for _, row in cands.iterrows():
        chrom, pos0 = row["chrom"], int(row["pos"]) - 1
        if chrom not in refset:
            continue
        ref_b, alt_b, snv = row["ref"].upper(), row["alt"].upper(), row["snv_id"]
        seen = {}                                            # (bc,umi) -> base (within-run dedup)
        for col in bam.pileup(chrom, pos0, pos0 + 1, truncate=True, min_base_quality=0,
                              max_depth=max_depth):
            if col.reference_pos != pos0:
                continue
            for pr in col.pileups:
                rd = pr.alignment
                if rd.mapping_quality < min_mapq or pr.is_del or pr.is_refskip or not rd.has_tag("CB"):
                    continue
                bc = _rc(rd.get_tag("CB")); key = (bc, _umi(rd))   # annotation orientation (see _rc)
                if key in seen:
                    continue
                b = rd.query_sequence[pr.query_position].upper()
                seen[key] = b
                cell = out.setdefault((bc, snv), [0, 0])
                if b == ref_b:
                    cell[0] += 1
                elif b == alt_b:
                    cell[1] += 1
    bam.close()
    return out


# ── S3: build combined matrices + append driver genotypes ──────────────────────
def _driver_cols(extraction_dir, sample, drivers):
    """{barcode: {driver: (ref_umi_like, alt_umi_like)}} from <sample>_mutation_matrix.csv (WT/MUT counts)."""
    p = os.path.join(extraction_dir, f"{sample}_mutation_matrix.csv")
    out = {}
    if not os.path.exists(p):
        return out
    with open(p) as fh:
        rd = csv.DictReader(fh)
        cols = rd.fieldnames or []
        want = {d: (f"{d}_WT", f"{d}_MUT") for d in drivers if f"{d}_WT" in cols and f"{d}_MUT" in cols}
        for r in rd:
            bc = r["cell_barcode"]
            out[bc] = {d: (int(r[w] or 0), int(r[m] or 0)) for d, (w, m) in want.items()}
    return out


def build_patient_matrices(patient, chrom, timepoints, manifest_bams, germline_vcfs, extraction_dir,
                           drivers, out_dir, min_dp=10, min_gq=20, min_mapq=20, limit=None,
                           sample_library=None, annotations=None, celltypes=None):
    """timepoints: list of (stage, [sample1, sample2...]) in chronological order (multi-run timepoints merged).
    When `annotations` is given, cells are RESTRICTED to barcodes annotated for each sample's library
    (sample_library[sample]) — the required single-cell filter — and a celltype column is written."""
    os.makedirs(out_dir, exist_ok=True)
    all_vcfs = [v for ss in timepoints for s in ss[1] for v in [germline_vcfs.get(s)] if v]
    cands = extract_het_union(all_vcfs, chrom, min_dp, min_gq)
    if limit:
        cands = cands.head(int(limit))                       # test scope only
    cands.to_csv(os.path.join(out_dir, f"{chrom}_het_candidates.tsv"), sep="\t", index=False)
    print(f"[{patient}] {chrom}: {len(cands)} union germline het-SNVs")
    if cands.empty:
        raise SystemExit("no het-SNVs -- check chrom/germline VCFs")

    alt_frames, ref_frames, meta_rows = [], [], []
    for stage, samples in timepoints:
        # merge all run-BAMs of this timepoint (sum ALT/REF UMIs per barcode)
        merged = {}
        drv = {}
        for s in samples:
            bam = manifest_bams.get(s)
            if not bam or not os.path.exists(bam):
                print(f"  [warn] {s}: bam missing"); continue
            valid = annotations.get((sample_library or {}).get(s, "")) if annotations is not None else None
            pu = pileup_bam(bam, cands, min_mapq)
            for (bc, snv), (rc, ac) in pu.items():
                if valid is not None and bc not in valid:
                    continue                                         # drop unannotated cells (required filter)
                m = merged.setdefault((bc, snv), [0, 0]); m[0] += rc; m[1] += ac
            for bc, dd in _driver_cols(extraction_dir, s, drivers).items():
                if valid is not None and bc not in valid:
                    continue
                d0 = drv.setdefault(bc, {})
                for dname, (w, mt) in dd.items():
                    cur = d0.setdefault(dname, [0, 0]); cur[0] += w; cur[1] += mt
        # assemble cell x SNV for this timepoint
        bcs = sorted({bc for (bc, _snv) in merged} | set(drv))
        idx = [f"{stage}::{bc}" for bc in bcs]
        a = pd.DataFrame(0, index=idx, columns=cands["snv_id"].tolist(), dtype=np.int32)
        r = pd.DataFrame(0, index=idx, columns=cands["snv_id"].tolist(), dtype=np.int32)
        for (bc, snv), (rc, ac) in merged.items():
            cell = f"{stage}::{bc}"; r.at[cell, snv] = rc; a.at[cell, snv] = ac
        for bc, dd in drv.items():
            cell = f"{stage}::{bc}"
            for dname, (w, mt) in dd.items():
                r.at[cell, dname] = w; a.at[cell, dname] = mt
        alt_frames.append(a); ref_frames.append(r)
        lib0 = (sample_library or {}).get(samples[0], "")
        meta_rows += [{"cell": c, "sample": samples[0], "stage": stage,
                       "celltype": (celltypes or {}).get((lib0, c.split("::", 1)[1]), "NA")} for c in idx]
        print(f"  {stage}: {len(idx)} cells")

    feats = cands["snv_id"].tolist() + list(drivers)
    A = pd.concat(alt_frames).reindex(columns=feats, fill_value=0)
    R = pd.concat(ref_frames).reindex(columns=feats, fill_value=0)
    A.to_csv(os.path.join(out_dir, "combined_alt.txt"), sep="\t")
    R.to_csv(os.path.join(out_dir, "combined_ref.txt"), sep="\t")
    pd.DataFrame(meta_rows).to_csv(os.path.join(out_dir, "cell_metadata.tsv"), sep="\t", index=False)
    print(f"[{patient}] wrote matrices: {A.shape[0]} cells x {A.shape[1]} features -> {out_dir}")


def _load_manifest_bams(manifest):
    bams = {}
    with open(manifest) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            bams[r["sample"]] = r["bam"]
    return bams


def _germline_vcfs(manifest):
    vc = {}
    with open(manifest) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            vc[r["sample"]] = os.path.join(r["outdir"], "annotation", "gatk", "germline.vcf.gz")
    return vc


def _load_manifest_libraries(manifest):
    """{sample -> cellHarmony library name} from the discovery manifest (the authoritative bridge
    between a KINNEX BAM sample and its short-read cellHarmony annotation library)."""
    lib = {}
    with open(manifest) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lib[r["sample"]] = r.get("library", "")
    return lib


def _load_annotations(path, libraries):
    """Parse the cellHarmony unique-barcode file ('BARCODE-1.LIBRARY<TAB>celltype' per annotated cell).
    Returns {library -> {bare_barcode}} and {(library, bare_barcode) -> celltype}, restricted to the
    requested libraries. Barcodes are stored bare (no '-1') and are already in annotation orientation,
    so they match _rc(BAM CB) and the variant_extraction driver barcodes."""
    want = set(libraries)
    ann, ct = {}, {}
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            cell, _, celltype = ln.partition("\t")
            bc, _, lib = cell.partition(".")
            if lib not in want:
                continue
            if bc.endswith("-1"):
                bc = bc[:-2]
            ann.setdefault(lib, set()).add(bc)
            ct[(lib, bc)] = celltype
    return ann, ct


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--patient", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--timepoints", nargs="+", required=True,
                    help="ordered stage:sample[,sample2...] tokens, e.g. preHMA:5801M_pre sAML:5801A_pre_1,5801A_pre_2")
    ap.add_argument("--manifest", required=True, help="discovery manifest (sample,uid,library,group,bam,outdir,job_script)")
    ap.add_argument("--extraction-dir", dest="extraction_dir", required=True, help="results_variant_extraction (driver genotypes)")
    ap.add_argument("--drivers", nargs="+", required=True, help="driver columns to append (anchor + LOH target), e.g. SRSF2p.P95R RUNX1p.W279*")
    ap.add_argument("--output", required=True)
    ap.add_argument("--min-dp", dest="min_dp", type=int, default=10)
    ap.add_argument("--min-gq", dest="min_gq", type=int, default=20)
    ap.add_argument("--min-mapq", dest="min_mapq", type=int, default=20)
    ap.add_argument("--annotations", default=None,
                    help="cellHarmony unique-barcode file (BARCODE-1.LIBRARY<TAB>celltype). When given, "
                         "cells are RESTRICTED to annotated barcodes for each sample's manifest library "
                         "(the required single-cell filter) and a celltype column is added to cell_metadata.")
    ap.add_argument("--limit", type=int, default=None, help="test only: first N het-SNVs")
    args = ap.parse_args()
    tps = []
    for tok in args.timepoints:
        stage, samples = tok.split(":", 1)
        tps.append((stage, samples.split(",")))
    sample_library = _load_manifest_libraries(args.manifest)
    annotations = celltypes = None
    if args.annotations:
        needed = {sample_library.get(s, "") for _st, ss in tps for s in ss}
        annotations, celltypes = _load_annotations(args.annotations, needed)
        n = sum(len(v) for v in annotations.values())
        print(f"[{args.patient}] annotation filter: {len(annotations)} libraries, {n} annotated barcodes")
    build_patient_matrices(args.patient, args.chrom, tps, _load_manifest_bams(args.manifest),
                           _germline_vcfs(args.manifest), args.extraction_dir, args.drivers, args.output,
                           min_dp=args.min_dp, min_gq=args.min_gq, min_mapq=args.min_mapq, limit=args.limit,
                           sample_library=sample_library, annotations=annotations, celltypes=celltypes)


if __name__ == "__main__":
    main()
