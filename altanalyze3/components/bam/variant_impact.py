"""altanalyze3 variant-impact

Cell-state-resolved expression/isoform impact of a called variant, end to end in one command:

  1. GENOTYPE each cell straight from the BAM's extended CIGAR (= match -> WT read, X mismatch ->
     MUT read) at the variant position. No reference FASTA -- the BAM encodes match/mismatch.
  2. For each requested matrix level (--level gene|isoform|both), per CELL STATE that has enough
     called MUT and WT cells, run the benchmarked MarkerFinder (cellHarmony.find_markers_from_adata)
     MUT-vs-WT, then IMPUTE the undetected (0-read) cells of that state from the signature and EXPAND
     the genotype, and render the MarkerFinder heatmap (labeled MUT/WT, PDF) on the expanded set.
  3. A COMBINED (all cell states pooled) run using the confident calls (MUT>=--mut-min, WT>=--wt-min).
  4. A cross-cell-state marker CONSISTENCY table (which genes/isoforms recur across cell states).

Only orchestration + genotyping + the spec'd imputation live here; markers and heatmaps come from the
benchmarked cellHarmony code.
"""

import os
import re
import logging

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import issparse

from altanalyze3.components.cellHarmony.cellHarmony_lite import normalize_adata
from altanalyze3.components.visualization.marker_heatmap_h5ad import generate_marker_heatmap_from_adata

logger = logging.getLogger(__name__)
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


# ---- barcode helpers (orientation-aware matching) ----
def _rc(seq):
    return "".join(_COMP.get(b, b) for b in reversed(seq))


def _core(barcode):
    bc = str(barcode).strip().split("-")[0].split(".")[0]
    m = re.match(r"^[ACGTN]+", bc.upper())
    return m.group(0) if m else bc.upper()


def _pick_orientation(ref_cores, other_cores):
    other = list(other_cores)
    direct = sum(1 for c in other if c in ref_cores)
    rc = sum(1 for c in other if _rc(c) in ref_cores)
    return ("rc", rc) if rc > direct else ("direct", direct)


def _safe(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name))


# ---- genotyping straight from the BAM (= / X extended CIGAR; no reference FASTA) ----
def _base_and_class_at(read, target0):
    """(query_base, '=' | 'X' | 'M') at the reference position, or (None, None) (deletion/skip)."""
    refpos = read.reference_start
    qpos = 0
    seq = read.query_sequence
    for op, length in read.cigartuples:
        if op in (0, 7, 8):                         # M, =, X -> consume query + ref
            if refpos <= target0 < refpos + length:
                qi = qpos + (target0 - refpos)
                base = seq[qi].upper() if seq and qi < len(seq) else None
                return base, {0: "M", 7: "=", 8: "X"}[op]
            refpos += length; qpos += length
        elif op == 1:                                # I -> query only
            qpos += length
        elif op in (2, 3):                           # D, N -> ref only
            refpos += length
        elif op in (4,):                             # S -> query only
            qpos += length
        if refpos > target0:
            break
    return None, None


def _indels_in_window(read, target0, window):
    """List of (type, length) indels whose reference position is within `window` bp of the target.
    The window tolerates left/right-alignment of the same short indel (CIGAR only, no FASTA)."""
    out = []
    refpos = read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            refpos += length
        elif op == 1:                                # insertion
            if abs(refpos - target0) <= window:
                out.append(("ins", length))
        elif op == 2:                                # deletion
            if refpos - window <= target0 <= refpos + length + window:
                out.append(("del", length))
            refpos += length
        elif op == 3:
            refpos += length
    return out


def genotype_from_bam(bam_path, chrom, pos1, vtype="SNV",
                      barcode_tags=("CB", "BC", "XC", "UB"), min_mapq=20, indel_window=5):
    """Per-cell allele read counts at a 1-based position, straight from the BAM (no reference FASTA).
    Only the TWO MAJOR alleles are counted (minor/error alleles ignored).

    SNV  : reference allele = consensus base of '=' (match) reads (== the genome reference); variant
           allele = major base of 'X' (mismatch) reads. mut = variant-allele reads, wt = reference reads.
    Indel: variant allele = the MAJOR indel (type,length carried by the most reads, within
           `indel_window` bp to absorb left/right-alignment); reference = reads spanning with no indel.
           mut = major-indel reads, wt = clean reads.
    Returns DataFrame[core_barcode] -> (mut, wt); the resolved ref/alt allele are in .attrs.
    """
    import pysam
    from collections import Counter
    target = pos1 - 1
    is_indel = str(vtype).lower().startswith("indel")
    bam = pysam.AlignmentFile(bam_path, "rb")
    # for indels, fetch a window so reads whose (left/right-aligned) indel sits a few bp away are seen
    fs = max(0, target - indel_window) if is_indel else target
    fe = (target + indel_window + 1) if is_indel else (target + 1)
    recs = []                                        # (core, base, cls, indels[list of (type,len)])
    for read in bam.fetch(chrom, fs, fe):
        if read.is_unmapped or read.mapping_quality < min_mapq or read.cigartuples is None:
            continue
        bc = None
        for tg in barcode_tags:
            try:
                bc = read.get_tag(tg); break
            except KeyError:
                continue
        if not bc:
            continue
        base, cls = _base_and_class_at(read, target)
        indels = _indels_in_window(read, target, indel_window) if is_indel else []
        recs.append((_core(bc), base, cls, indels))
    bam.close()

    mut, wt = {}, {}
    if is_indel:
        # two major alleles: the MAJOR indel (type,length carried by the most reads) vs reference
        # (reads spanning the locus with no indel). Minor indels (errors / a different event) are
        # excluded from both, exactly like the SNV major-alt-vs-reference logic.
        allele_counts = Counter(a for _, _, _, inds in recs for a in set(inds))
        major = allele_counts.most_common(1)[0][0] if allele_counts else None
        ref_allele = "ref(no-indel)"
        alt_allele = (f"{major[0]}{major[1]}" if major else None)
        for core, base, cls, inds in recs:
            if major and major in inds:
                mut[core] = mut.get(core, 0) + 1
            elif not inds and base is not None:
                wt[core] = wt.get(core, 0) + 1
    else:
        ref_allele = (Counter(b for _, b, c, i in recs if c == "=" and b).most_common(1) or [(None,)])[0][0]
        alt_allele = (Counter(b for _, b, c, i in recs if c == "X" and b).most_common(1) or [(None,)])[0][0]
        n_M = sum(1 for _, _, c, _ in recs if c == "M")
        if n_M and ref_allele is None:
            logger.warning("genotype_from_bam %s:%d uses standard 'M' CIGAR (no =/X) -- a reference "
                           "FASTA would be needed to call this locus.", chrom, pos1)
        for core, base, cls, ind in recs:
            if cls == "X" and base == alt_allele:
                mut[core] = mut.get(core, 0) + 1
            elif base == ref_allele:
                wt[core] = wt.get(core, 0) + 1
    bcs = sorted(set(mut) | set(wt))
    df = pd.DataFrame({"mut": [mut.get(b, 0) for b in bcs], "wt": [wt.get(b, 0) for b in bcs]}, index=bcs)
    df.attrs["ref_allele"] = ref_allele
    df.attrs["alt_allele"] = alt_allele
    logger.info("genotype %s:%d (%s) ref=%s alt=%s -> %d cells", chrom, pos1, vtype,
                ref_allele, alt_allele, len(df))
    return df


def _call_genotype(counts, mut_min, wt_min):
    """counts: DataFrame(mut, wt) per core barcode -> Series core->{MUT,WT,UNK}."""
    g = np.where(counts["mut"].to_numpy() >= mut_min, "MUT",
                 np.where(counts["wt"].to_numpy() >= wt_min, "WT", "UNK"))
    return pd.Series(g, index=counts.index)


# ---- cell states ----
def load_cell_states(cell_annot_path):
    mapping = {}
    with open(cell_annot_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or "\t" not in line:
                continue
            bc, clust = line.split("\t")[:2]
            mapping[_core(bc)] = clust
    return mapping


def _load_symbol_map(path):
    if not path or not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path, sep="\t", header=None, dtype=str)
        return dict(zip(df[0], df[1]))
    except Exception as e:
        logger.warning("gene_symbol map not loaded (%s)", e)
        return None


# ---- imputation (spec'd: median-normalized signature centroids, threshold from true positives) ----
def _corr_rows(V, centroid):
    Vc = V - V.mean(axis=1, keepdims=True)
    cc = centroid - centroid.mean()
    num = Vc @ cc
    den = np.sqrt((Vc ** 2).sum(axis=1) * (cc ** 2).sum())
    den[den == 0] = np.nan
    return num / den


def impute_unknown(expr_all, geno_called, geno_unknown, signature, min_called=10):
    """Predict MUT/WT for the undetected cells from the MUT/WT signature (per cell type).
    expr_all: cells x features (normalized). geno_called: Series barcode->MUT/WT. Threshold set
    from the called (true-positive) cells: 10th percentile of correctly-scored called-cell margins."""
    sig = [f for f in signature if f in expr_all.columns]
    called = geno_called[geno_called.isin(["MUT", "WT"])]
    if len(sig) < 3 or (called == "MUT").sum() < min_called or (called == "WT").sum() < min_called:
        return None
    mut_c = np.median(expr_all.loc[called.index[called == "MUT"], sig].to_numpy(dtype=float), axis=0)
    wt_c = np.median(expr_all.loc[called.index[called == "WT"], sig].to_numpy(dtype=float), axis=0)
    ref = np.median(np.vstack([mut_c, wt_c]), axis=0)   # method B (validated best): median of centroids

    def score(bcs):
        # median-center by the median-of-centroids, then nearest centroid (Euclidean). >0 -> nearer MUT.
        V = expr_all.loc[bcs, sig].to_numpy(dtype=float) - ref
        dM = np.linalg.norm(V - (mut_c - ref), axis=1)
        dW = np.linalg.norm(V - (wt_c - ref), axis=1)
        return dW - dM

    cs = pd.Series(score(called.index), index=called.index)
    correct = ((cs > 0) & (called == "MUT")) | ((cs < 0) & (called == "WT"))
    margin = np.abs(cs[correct])
    thr = float(np.percentile(margin, 10)) if len(margin) else 0.0
    unk = [b for b in geno_unknown if b in expr_all.index]
    if unk:
        us = pd.Series(score(unk), index=unk)
        pred = np.where(us >= thr, "MUT", np.where(us <= -thr, "WT", "UNK"))
        out = pd.DataFrame({"barcode": unk, "imputed": pred, "score": us.values, "threshold": thr})
    else:
        out = pd.DataFrame(columns=["barcode", "imputed", "score", "threshold"])
    # VERIFICATION: the classifier must re-predict the ORIGINAL called cells it was built from
    # (resubstitution). concordance = fraction whose sign(score) matches the called genotype. A low
    # value means the signature does not separate the genuine calls -> the expansion is not trustworthy.
    sign = pd.Series(np.where(cs > 0, "MUT", "WT"), index=called.index)
    out.attrs["concordance"] = float((sign == called).mean())
    out.attrs["mut_recovery"] = float((sign[called == "MUT"] == "MUT").mean())
    out.attrs["wt_recovery"] = float((sign[called == "WT"] == "WT").mean())
    out.attrs["n_called_mut"] = int((called == "MUT").sum())
    out.attrs["n_called_wt"] = int((called == "WT").sum())
    return out


# ---- expression helpers ----
def _expr_subset(norm, cells, features):
    """Dense (cells x features) frame from the already-normalized adata, for a SMALL feature set
    (the MarkerFinder signature) -- never densifies the full (isoform) matrix."""
    fset = set(norm.var_names)
    feats = [f for f in features if f in fset]
    sub = norm[list(cells), feats]
    X = sub.X.toarray() if issparse(sub.X) else np.asarray(sub.X)
    return pd.DataFrame(X, index=list(sub.obs_names), columns=feats)


def _markerfinder(base, cells, labels, out_dir, tag, n_markers, species=None):
    """Production MarkerFinder + raster heatmap (generate_marker_heatmap_from_adata, marker_method
    'markerfinder', imshow body). Returns the top markers as (feature, direction MUT/WT)."""
    sub = base[base.obs_names.isin(list(cells))].copy()
    sub.obs["genotype"] = pd.Series(labels, index=list(cells)).loc[sub.obs_names].astype(str).values
    normalize_adata(sub)
    try:
        res = generate_marker_heatmap_from_adata(
            sub, out=os.path.join(out_dir, f"{tag}_marker_heatmap.pdf"),
            cluster_key="genotype", top_n=n_markers, marker_method="markerfinder",
            cells_per_cluster=100000, seed=0, export_networks=False, network_top_n=0, network_jobs=1,
            species=species, write_heatmap_tsv=True, write_expression_tsv=False,
            write_heatmap_cache=False, pval_threshold=None)
    except Exception as e:
        logger.warning("MarkerFinder heatmap failed for %s: %s", tag, e)
        return None
    try:
        m = pd.read_csv(res["markers_tsv"], sep="\t")
        gcol = "gene" if "gene" in m.columns else m.columns[0]
        ccol = "cluster" if "cluster" in m.columns else m.columns[1]
        return pd.DataFrame({"feature": m[gcol].astype(str), "direction": m[ccol].astype(str)})
    except Exception:
        return None


# ---- accuracy gate (>=85% re-classification of the ORIGINAL calls) ----
ACCURACY_MIN = 0.85


def _write_verification(res, out_dir, tag):
    """Write how well the classifier re-predicts the ORIGINAL called cells (resubstitution; the
    classifier is centroid-based so train==test is acceptable per spec)."""
    a = res.attrs
    with open(os.path.join(out_dir, f"{tag}_classifier_verification.txt"), "w") as fh:
        fh.write("classifier\tmedian-centered nearest-centroid (MUT vs WT)\n")
        fh.write(f"called_MUT\t{a.get('n_called_mut')}\ncalled_WT\t{a.get('n_called_wt')}\n")
        fh.write(f"pct_correct_MUT\t{100 * a.get('mut_recovery', 0):.1f}\n")
        fh.write(f"pct_correct_WT\t{100 * a.get('wt_recovery', 0):.1f}\n")
        fh.write(f"accuracy\t{100 * a.get('concordance', 0):.1f}\nrequired\t{100 * ACCURACY_MIN:.0f}\n")
        fh.write(f"status\t{'PASS' if a.get('concordance', 0) >= ACCURACY_MIN else 'FAIL (<85%; not deployed)'}\n")


# ---- per (sample, variant, level): MarkerFinder + classify + impute, per cell state AND combined ----
def process_level(base, norm, states, counts, sample, variant, level, out_dir,
                  mut_min, wt_min, min_cells, n_markers, species=None):
    core = {b: _core(b) for b in base.obs_names}
    h5c = {_core(b): b for b in base.obs_names}
    vor, _ = _pick_orientation(set(h5c), counts.index)
    cnt = counts.copy(); cnt.index = [_core(_rc(c)) if vor == "rc" else c for c in cnt.index]
    cnt = cnt[~cnt.index.duplicated()]
    cor, _ = _pick_orientation(set(h5c), states.keys())
    stmap = {(_core(_rc(k)) if cor == "rc" else k): v for k, v in states.items()}
    st = pd.Series({b: stmap.get(core[b]) for b in base.obs_names})
    mco = pd.Series({b: cnt["mut"].get(core[b], 0) for b in base.obs_names})
    wco = pd.Series({b: cnt["wt"].get(core[b], 0) for b in base.obs_names})
    g = pd.Series(np.where(mco >= mut_min, "MUT", np.where(wco >= wt_min, "WT", "UNK")), index=base.obs_names)

    summary, deployed = [], []
    # every read-genotyped cell is in the final sample-level table as 'called'
    for b in g.index[g.isin(["MUT", "WT"])]:
        deployed.append({"sample": sample, "variant": variant, "level": level, "barcode": b,
                         "cell_state": st.get(b), "genotype": g[b], "source": "called"})

    # one unit per cell type with both groups present, plus the all-cell-types-combined unit
    units = [("combined", g.index[g.isin(["MUT", "WT"])])]
    for s in sorted(st.dropna().unique()):
        units.append((s, st.index[(st == s) & g.isin(["MUT", "WT"])]))

    for unit, cells in units:
        gl = g.loc[cells]; mut, wt = cells[gl == "MUT"], cells[gl == "WT"]
        if len(mut) < min_cells or len(wt) < min_cells:
            continue
        is_comb = (unit == "combined")
        udir = os.path.join(out_dir, "combined") if is_comb else os.path.join(out_dir, "per_celltype", _safe(unit))
        os.makedirs(udir, exist_ok=True)
        tag = "combined" if is_comb else _safe(unit)
        lab = ["MUT"] * len(mut) + ["WT"] * len(wt)
        # MarkerFinder on ALL the called cells of this unit (no subset) + labeled raster heatmap
        mk = _markerfinder(base, list(mut) + list(wt), lab, udir, tag, n_markers, species=species)
        if mk is None or not len(mk):
            continue
        sig = mk["feature"].astype(str).tolist()
        unk = g.index[g == "UNK"] if is_comb else st.index[(st == unit) & (g == "UNK")]
        res = impute_unknown(_expr_subset(norm, list(mut) + list(wt) + list(unk), sig),
                             g.loc[list(mut) + list(wt)], pd.Index(unk), sig, min_called=min_cells)
        if res is None:
            continue
        acc = float(res.attrs.get("concordance", 0.0))
        _write_verification(res, udir, tag)
        imp = res[res["imputed"].isin(["MUT", "WT"])] if len(res) else res
        if len(imp):
            imp.to_csv(os.path.join(udir, f"{tag}_imputed.csv"), index=False)
        n_im = int((imp["imputed"] == "MUT").sum()) if len(imp) else 0
        n_iw = int((imp["imputed"] == "WT").sum()) if len(imp) else 0
        deploy = (not is_comb) and acc >= ACCURACY_MIN          # deploy per-cell-state calls only when >=85%
        if deploy:
            for _, r in imp.iterrows():
                deployed.append({"sample": sample, "variant": variant, "level": level,
                                 "barcode": r["barcode"], "cell_state": unit,
                                 "genotype": r["imputed"], "source": "imputed"})
        summary.append({"sample": sample, "variant": variant, "level": level, "unit": unit,
                        "n_MUT": len(mut), "n_WT": len(wt),
                        "pct_correct_MUT": round(100 * res.attrs.get("mut_recovery", 0), 1),
                        "pct_correct_WT": round(100 * res.attrs.get("wt_recovery", 0), 1),
                        "accuracy_pct": round(100 * acc, 1), "pass_85": acc >= ACCURACY_MIN,
                        "n_imputed_MUT": n_im, "n_imputed_WT": n_iw, "deployed": deploy})
        logger.info("[%s|%s|%s] %s: MUT=%d WT=%d acc=%.1f%% imputed=%d/%d %s", sample, variant, level,
                    unit, len(mut), len(wt), 100 * acc, n_im, n_iw,
                    "DEPLOYED" if deploy else ("(combined)" if is_comb else "below 85%"))
    return summary, deployed


def _parse_variants(spec):
    """Accept 'chr:pos:label[:type]' (semicolon-separated) OR a file of 'chrom start end label type'.
    Returns (chrom, pos, label, vtype) with vtype SNV (default) or Indel."""
    out = []
    if spec and os.path.exists(spec):
        with open(spec) as fh:
            for line in fh:
                p = line.split()
                if len(p) >= 2 and not line.startswith("#"):
                    chrom = p[0] if p[0].startswith("chr") else "chr" + p[0]
                    out.append((chrom, int(p[1]), p[3] if len(p) > 3 else f"{chrom}_{p[1]}",
                                p[4] if len(p) > 4 else "SNV"))
    else:
        for tok in (spec or "").split(";"):
            tok = tok.strip()
            if not tok:
                continue
            parts = tok.split(":")
            chrom = parts[0] if parts[0].startswith("chr") else "chr" + parts[0]
            out.append((chrom, int(parts[1]), parts[2] if len(parts) > 2 else f"{chrom}_{parts[1]}",
                        parts[3] if len(parts) > 3 else "SNV"))
    return out


def run_variant_impact(args):
    metadata = pd.read_csv(args.metadata, sep="\t")
    out_root = os.path.abspath(getattr(args, "output", None) or os.path.join(os.getcwd(), "variant_impact"))
    os.makedirs(out_root, exist_ok=True)
    levels = ["gene", "isoform"] if args.level == "both" else [args.level]
    variants = _parse_variants(args.variants)
    if not variants:
        raise ValueError("--variants required: 'chr:pos:label[;...]' or a mutations file path")
    cell_annot = getattr(args, "cell_annot", None)
    sym = _load_symbol_map(getattr(args, "gene_symbol", None))
    only = set(getattr(args, "samples", None) or [])
    states = load_cell_states(cell_annot) if cell_annot and os.path.exists(cell_annot) else None
    all_summary, all_calls = [], []

    for _, r in metadata.iterrows():
        uid = str(r.get("uid", "")).strip()
        library = str(r.get("library", uid)).strip()
        bam = str(r.get("bam", ""))
        if not uid or (only and uid not in only and library not in only):
            continue
        if not os.path.exists(bam):
            logger.warning("[%s] bam not found: %s -- skipped", uid, bam); continue
        bam_dir = os.path.dirname(bam)
        for (chrom, pos, label, vtype) in variants:
            counts = genotype_from_bam(bam, chrom, pos, vtype=vtype,
                                       min_mapq=getattr(args, "min_mapq", 20),
                                       indel_window=getattr(args, "indel_window", 5))
            logger.info("[%s] %s @ %s:%d (%s) ref=%s alt=%s -> %d cells genotyped",
                        uid, label, chrom, pos, vtype, counts.attrs.get("ref_allele"),
                        counts.attrs.get("alt_allele"), len(counts))
            for level in levels:
                h5 = os.path.join(bam_dir, f"{library}-{level}.h5ad")
                if not os.path.exists(h5):
                    logger.warning("[%s] %s h5ad not found: %s -- skipped", uid, level, h5); continue
                base = ad.read_h5ad(h5)
                if level == "gene" and sym:
                    base.var_names = [sym.get(g, g) for g in base.var_names]
                    base.var_names_make_unique()
                st_map = states
                if st_map is None and "cluster" in base.obs.columns:
                    st_map = {_core(b): c for b, c in zip(base.obs_names, base.obs["cluster"])}
                if st_map is None:
                    logger.warning("[%s] no --cell_annot and no obs['cluster']; skipped", uid); continue
                norm = base.copy(); normalize_adata(norm)
                out_dir = os.path.join(out_root, _safe(uid), _safe(label), level)
                os.makedirs(out_dir, exist_ok=True)
                try:
                    summ, dep = process_level(base, norm, st_map, counts, uid, label, level, out_dir,
                                              mut_min=args.mut_min, wt_min=args.wt_min,
                                              min_cells=args.min_cells, n_markers=args.top_n)
                    all_summary.extend(summ); all_calls.extend(dep)
                except Exception as e:
                    logger.exception("[%s|%s|%s] failed: %s", uid, label, level, e)
    if all_summary:
        pd.DataFrame(all_summary).to_csv(os.path.join(out_root, "variant_impact_summary.tsv"),
                                         sep="\t", index=False)
    if all_calls:
        # FINAL sample-level variant result: read-called cells + imputed cells from the cell states
        # whose classifier re-predicted the originals at >= 85%.
        pd.DataFrame(all_calls).to_csv(os.path.join(out_root, "sample_variant_calls.tsv"),
                                       sep="\t", index=False)
    logger.info("variant-impact done -> %s  (per-unit summary: variant_impact_summary.tsv ; "
                "final calls: sample_variant_calls.tsv)", out_root)
