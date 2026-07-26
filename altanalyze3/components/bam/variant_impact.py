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
from altanalyze3.components.visualization.marker_heatmap_h5ad import (
    generate_marker_heatmap_from_adata, _build_heatmap, _plot_heatmap)

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
            elif ref_allele is not None and base is not None and base == ref_allele:
                # `base` is None for a read carrying a DELETION / splice-skip over the locus, and
                # `ref_allele` is None at a ~100%-VAF locus (no '=' reads at all). Without these two
                # guards `None == None` counted every no-base read as WILD-TYPE, inventing a WT
                # population for clonal/hemizygous variants (e.g. GATA1 c.212A>G, VAF 0.999: 54
                # phantom WT cells where variant_extraction correctly reports 0).
                wt[core] = wt.get(core, 0) + 1
    bcs = sorted(set(mut) | set(wt))
    df = pd.DataFrame({"mut": [mut.get(b, 0) for b in bcs], "wt": [wt.get(b, 0) for b in bcs]}, index=bcs)
    df.attrs["ref_allele"] = ref_allele
    df.attrs["alt_allele"] = alt_allele
    logger.info("genotype %s:%d (%s) ref=%s alt=%s -> %d cells", chrom, pos1, vtype,
                ref_allele, alt_allele, len(df))
    return df


def _inserted_seqs_in_window(read, target0, window):
    """(refpos, inserted_query_seq) for each insertion (I) op whose reference position is within
    `window` bp of the target. Lets a supervised query confirm the EXACT inserted bases (not just an
    insertion of matching length), while tolerating left/right-alignment via the window."""
    out = []
    refpos = read.reference_start
    qpos = 0
    seq = read.query_sequence or ""
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            refpos += length; qpos += length
        elif op == 1:                                    # insertion (query only)
            if abs(refpos - target0) <= window:
                out.append((refpos, seq[qpos:qpos + length].upper()))
            qpos += length
        elif op in (2, 3):                               # deletion / skip (ref only)
            refpos += length
        elif op == 4:                                    # soft clip (query only)
            qpos += length
    return out


def genotype_from_bam_supervised(bam_path, chrom, pos1, ref, alt,
                                 barcode_tags=("CB", "BC", "XC", "UB"), min_mapq=1, indel_window=5):
    """FULLY SUPERVISED per-cell genotype: confirm the EXACT pre-defined REF>ALT allele from the caller,
    with NO allele discovery (unlike genotype_from_bam, which takes the major mismatch/indel in the
    reads). Same (mut, wt) DataFrame + .attrs contract as genotype_from_bam, so it is a drop-in.

    min_mapq defaults to 1 -- the SAME threshold the discovery callers used (kinnex_variant_pipeline.sh:
    LongcallR --min-mapq 1, Clair3 --min_mq 5). The callers accept low-MAPQ reads in repeat/homopolymer
    regions; a MAPQ>=20 filter would silently drop the very reads that support those indels and lose
    them (measured: an A-insertion at chr21:5116990 has 573 supporting reads, ALL MAPQ 0-14).

      SNV       : mut = reads whose base == ALT ; wt = reads whose base == REF (others ignored).
      Insertion : mut = reads carrying an insertion within indel_window whose inserted sequence ==
                  ALT[len(REF):] (exact bases) ; wt = reads spanning cleanly with no such insertion.
      Deletion  : mut = reads carrying a deletion within indel_window of length len(REF)-len(ALT) ;
                  wt = reads spanning cleanly with no such deletion.
    Left/right-alignment is absorbed by the window; the allele identity is checked by sequence (ins)
    or exact length at the locus (del), not by "the major event in the pile".
    """
    import pysam
    ref = str(ref).upper(); alt = str(alt).upper()
    target = pos1 - 1
    lr, la = len(ref), len(alt)
    is_snv = (lr == 1 and la == 1)
    is_ins = la > lr
    is_del = lr > la
    exp_ins = alt[lr:] if (is_ins and ref == alt[:lr]) else (alt if is_ins else "")
    del_len = lr - la
    bam = pysam.AlignmentFile(bam_path, "rb")
    fs = target if is_snv else max(0, target - indel_window)
    fe = (target + 1) if is_snv else (target + max(lr, la) + indel_window + 1)
    mut, wt = {}, {}
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
        core = _core(bc)
        base, _cls = _base_and_class_at(read, target)
        if is_snv:
            if base == alt:
                mut[core] = mut.get(core, 0) + 1
            elif base == ref:
                wt[core] = wt.get(core, 0) + 1
        elif is_ins:
            ins = _inserted_seqs_in_window(read, target, indel_window)
            if any(s == exp_ins for _, s in ins):
                mut[core] = mut.get(core, 0) + 1
            elif base is not None and not ins:
                wt[core] = wt.get(core, 0) + 1
        elif is_del:
            dels = [ln for (t, ln) in _indels_in_window(read, target, indel_window) if t == "del"]
            if del_len in dels:
                mut[core] = mut.get(core, 0) + 1
            elif base is not None and not dels:
                wt[core] = wt.get(core, 0) + 1
        # MNV / complex (len(ref)==len(alt)>1) are neither snv/ins/del here; the caller set has none
        # (verified), so they are left uncounted rather than silently mis-called.
    bam.close()
    bcs = sorted(set(mut) | set(wt))
    df = pd.DataFrame({"mut": [mut.get(b, 0) for b in bcs], "wt": [wt.get(b, 0) for b in bcs]}, index=bcs)
    df.attrs["ref_allele"] = ref; df.attrs["alt_allele"] = alt
    logger.info("genotype(supervised) %s:%d %s>%s -> %d cells (MUT>=1: %d)", chrom, pos1, ref, alt,
                len(df), int((df["mut"] > 0).sum()))
    return df


def genotype_from_variant_extraction(matrix_csv, label):
    """Per-cell allele read counts read from a `variant_extraction.py` output
    `<sample>_mutation_matrix.csv` -- the same (mut, wt) contract as `genotype_from_bam`, so the two
    are interchangeable as the genotype source.

    That matrix carries, per cell, `<label>_MUT`, `<label>_WT` and `<label>_GENOTYPE` for every scanned
    locus; `cell_barcode` is written reverse-complemented (junction/h5ad orientation). Use this when the
    variants have already been extracted, instead of re-scanning the BAM.
    Returns DataFrame[barcode] -> (mut, wt), or None if the locus is absent.
    """
    m = pd.read_csv(matrix_csv, low_memory=False)
    mut_col, wt_col = f"{label}_MUT", f"{label}_WT"
    if mut_col not in m.columns or wt_col not in m.columns:
        logger.warning("variant_extraction matrix %s has no locus %s", os.path.basename(matrix_csv), label)
        return None
    df = pd.DataFrame(
        {"mut": pd.to_numeric(m[mut_col], errors="coerce").fillna(0).astype(int).to_numpy(),
         "wt": pd.to_numeric(m[wt_col], errors="coerce").fillna(0).astype(int).to_numpy()},
        index=[str(b) for b in m["cell_barcode"]])
    df = df[~df.index.duplicated()]
    df.attrs["ref_allele"] = "ref"
    df.attrs["alt_allele"] = label
    logger.info("genotype_from_variant_extraction %s: %d cells (MUT>=1: %d, WT>=1: %d)", label,
                len(df), int((df["mut"] >= 1).sum()), int((df["wt"] >= 1).sum()))
    return df


def variant_extraction_loci(readcounts_tsv):
    """{(chrom, pos): label} for every locus in a `variant_extraction.py` <sample>_variant_readcounts.tsv,
    so a caller holding genomic coordinates can resolve the matrix column label."""
    if not os.path.exists(readcounts_tsv):
        return {}
    d = pd.read_csv(readcounts_tsv, sep="\t")
    return {(str(r.chrom), int(r.pos)): str(r.label) for r in d.itertuples()}


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


# ---- 3-panel MarkerFinder output + PCA -------------------------------------------------------
# Per analysis unit we emit, using the SAME benchmarked MarkerFinder heatmap code throughout:
#   1  <tag>_marker_heatmap.pdf              TRAINING cells (read-genotyped MUT/WT), top marker genes
#   2  <tag>_2_predicted_marker_heatmap.pdf  NEWLY IMPUTED cells with their assigned MUT/WT labels,
#                                            markers recomputed FROM those cells
#   3  <tag>_3_combined_marker_heatmap.pdf   training + predicted as FOUR annotated groups, restricted
#                                            to the TRAINING genes in the TRAINING gene order
#   4  <tag>_4_combined_predictedGenes_marker_heatmap.pdf
#                                            same four-group combined plot as (3), but using the
#                                            PREDICTED-cell genes from (2), in the predicted gene order
#   5  <tag>_5_combined_pca.pdf              PCA of the plot-3 matrix (training vs imputed cells)
#   6  <tag>_6_training_volcano.pdf          volcano of the TRAINING MUT-vs-WT DEGs (cellHarmony
#                                            moderated t-test; fold>=1.2 non-log & raw p<0.05), top 10 up and
#                                            top 10 down labelled by fold  (+ _6_training_DEGs.tsv)
#   7  <tag>_7_predicted_volcano.pdf         same test on the PREDICTED cells, labelling ONLY the
#                                            training DEGs                (+ _7_predicted_DEGs.tsv)
G_TRAIN_MUT, G_TRAIN_WT = "MUT-training", "WT-training"
G_PRED_MUT, G_PRED_WT = "MUT-predicted", "WT-predicted"
GROUP_ORDER4 = [G_TRAIN_MUT, G_TRAIN_WT, G_PRED_MUT, G_PRED_WT]
# explicit RGB (editable PDF; no named colors)
_PCA_COLORS = {G_TRAIN_MUT: (0.706, 0.157, 0.129), G_TRAIN_WT: (0.122, 0.353, 0.616),
               G_PRED_MUT: (0.937, 0.604, 0.553), G_PRED_WT: (0.604, 0.761, 0.898)}


def _pca_from_heatmap(heatmap_df, column_clusters, cluster_order, out_path, title=""):
    """PCA derived from the combined heatmap matrix (genes x cells, row z-scored): cells projected on
    PC1/PC2, coloured by the four training/predicted groups."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    X = np.nan_to_num(heatmap_df.to_numpy(dtype=float).T, nan=0.0, posinf=0.0, neginf=0.0)  # cells x genes
    if X.shape[0] < 3 or X.shape[1] < 2:
        return None
    X = X - X.mean(axis=0, keepdims=True)
    U, S, _ = np.linalg.svd(X, full_matrices=False)          # PCA via SVD (numpy only)
    pc = U[:, :2] * S[:2]
    ev = (S ** 2) / max(float((S ** 2).sum()), 1e-12) * 100.0

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    cc = np.asarray([str(c) for c in column_clusters])
    for grp in cluster_order:
        m = cc == grp
        if not m.any():
            continue
        ax.scatter(pc[m, 0], pc[m, 1], s=14, linewidths=0, c=[_PCA_COLORS.get(grp, (0.5, 0.5, 0.5))],
                   alpha=0.85 if grp.endswith("training") else 0.55,
                   label=f"{grp} (n={int(m.sum())})")
    ax.set_xlabel(f"PC1 ({ev[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)")
    ax.set_title(title or "Training + imputed cells (training marker genes)", fontsize=10)
    ax.legend(fontsize=7, frameon=False, loc="best")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


# ---- differential expression + volcano (cellHarmony_differential moderated t-test) ----------
VOLCANO_FOLD = 1.2                     # NON-LOG fold-change threshold (user spec): 1.2x
VOLCANO_LOG2FC = float(np.log2(VOLCANO_FOLD))  # = 0.263 on the log2 axis (|log2FC| >= log2(1.2))
VOLCANO_ADJP = 0.05                    # significance p cutoff -- RAW p by default (see _volcano use_rawp)
_V_UP, _V_DOWN, _V_NS = (0.706, 0.157, 0.129), (0.122, 0.353, 0.616), (0.72, 0.74, 0.78)
_V_HI = (0.90, 0.49, 0.06)             # training DEGs carried onto the predicted volcano


def _moderated_de(base, case_cells, ctrl_cells, case="MUT", ctrl="WT"):
    """Limma-like moderated t-test with empirical-Bayes shrinkage -- the cellHarmony-differential
    method (cellHarmony_differential._moderated_t_test), run on log2(CP10k+1) single cells so the
    returned `log2fc` is a true log2 fold change."""
    from altanalyze3.components.cellHarmony.cellHarmony_differential import _moderated_t_test
    cells = list(case_cells) + list(ctrl_cells)
    if len(case_cells) < 2 or len(ctrl_cells) < 2:
        return pd.DataFrame(columns=["gene", "log2fc", "t", "pval", "fdr"])
    sub = base[cells].copy()
    X = sub.X.toarray() if issparse(sub.X) else np.asarray(sub.X, dtype=float)
    tot = X.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1.0
    sub.X = np.log2(X / tot * 1e4 + 1.0)          # cellHarmony log2CP10k convention
    sub.obs["vi_de"] = [case] * len(case_cells) + [ctrl] * len(ctrl_cells)
    try:
        df, _ = _moderated_t_test(sub, "vi_de", case, ctrl, "variant_impact")
    except Exception as e:
        logger.warning("moderated t-test failed: %s", e)
        return pd.DataFrame(columns=["gene", "log2fc", "t", "pval", "fdr"])
    return df


def _label_points(ax, rows, colors, y_of, x_of):
    """Annotate without piling labels on top of each other: labels go outward (left of left-hand
    points, right of right-hand points) and are fanned vertically in offset points."""
    rows = list(rows)
    left = sorted([i for i, r in enumerate(rows) if x_of(r) < 0], key=lambda i: y_of(rows[i]))
    right = sorted([i for i, r in enumerate(rows) if x_of(r) >= 0], key=lambda i: y_of(rows[i]))
    for side, idxs in (("left", left), ("right", right)):
        for k, i in enumerate(idxs):
            r = rows[i]
            ax.annotate(str(r["gene"]), (x_of(r), y_of(r)),
                        xytext=(-4 if side == "left" else 4, (k % 6) * 7 - 17),
                        textcoords="offset points", fontsize=5,
                        ha="right" if side == "left" else "left", va="center",
                        color=colors[i], annotation_clip=False)


def _volcano(df, out_path, title, label_genes=None, top_n=10,
             lfc_min=VOLCANO_LOG2FC, alpha=VOLCANO_ADJP, use_rawp=True):
    """Volcano on the log2 axis. Significant = |log2FC| >= lfc_min (= log2(1.2), a NON-LOG 1.2-fold
    cutoff) AND p < alpha, where p is the RAW p by default (`use_rawp=True`, the cellHarmony-differential
    `use_rawp` convention -- BH across ~16k genes on a few dozen single cells leaves nothing significant),
    or the BH-adjusted p when `use_rawp=False`.
    label_genes=None  -> label the top `top_n` up and down by log2 fold change. If nothing is
                         significant these are drawn in GREY and marked n.s., never in the DEG colours.
    label_genes given -> the TRAINING DEGs, in their own colour AND labelled, over this data set's
                         own up/down colouring. Returns this data set's significant genes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    pcol = "pval" if use_rawp else "fdr"
    plab = "raw p" if use_rawp else "adjusted p"
    d = df.dropna(subset=["log2fc", pcol]).copy()
    if d.empty:
        return []
    d["sig"] = (d[pcol] < alpha) & (d["log2fc"].abs() >= lfc_min)
    y = -np.log10(d[pcol].clip(lower=1e-300).to_numpy(dtype=float))
    x = d["log2fc"].to_numpy(dtype=float)
    up = (d["sig"] & (d["log2fc"] > 0)).to_numpy()
    dn = (d["sig"] & (d["log2fc"] < 0)).to_numpy()
    ns = ~(up | dn)
    yv = dict(zip(d["gene"].astype(str), y))

    fig, ax = plt.subplots(figsize=(6.0, 5.4))
    ax.scatter(x[ns], y[ns], s=5, c=[_V_NS], linewidths=0, alpha=0.55)
    ax.scatter(x[up], y[up], s=11, c=[_V_UP], linewidths=0, label=f"up (n={int(up.sum())})")
    ax.scatter(x[dn], y[dn], s=11, c=[_V_DOWN], linewidths=0, label=f"down (n={int(dn.sum())})")
    for v in (-lfc_min, lfc_min):                       # single 2-point lines
        ax.axvline(v, linestyle="--", linewidth=0.7, color=(0, 0, 0))
    ax.axhline(-np.log10(alpha), linestyle="--", linewidth=0.7, color=(0, 0, 0))

    note = ""
    if label_genes is None:
        pool, sig_pool = d[d["sig"]], True
        if pool.empty:                                  # nothing significant -> label top by fold, GREY
            pool, sig_pool = d[d["log2fc"].abs() >= lfc_min], False
            if pool.empty:
                pool = d
            note = "  (no gene significant; top by fold shown in grey, n.s.)"
        lab = pd.concat([pool.nlargest(top_n, "log2fc"), pool.nsmallest(top_n, "log2fc")])
        lab_colors = [(_V_UP if v > 0 else _V_DOWN) if sig_pool else _V_NS for v in lab["log2fc"]]
    else:
        lab = d[d["gene"].astype(str).isin(set(map(str, label_genes)))]
        if len(lab):                                    # training DEGs: own colour AND labels
            ax.scatter(lab["log2fc"].to_numpy(dtype=float),
                       np.array([yv[str(g)] for g in lab["gene"]]),
                       s=26, c=[_V_HI], linewidths=0, label=f"training DEG (n={len(lab)})")
        else:
            note = "  (no training DEGs to carry over)"
        lab_colors = [_V_HI] * len(lab)
    if len(lab):
        _label_points(ax, [r for _, r in lab.iterrows()], lab_colors,
                      y_of=lambda r: yv[str(r["gene"])], x_of=lambda r: float(r["log2fc"]))
    ax.set_xlabel("log2 fold change (MUT vs WT)")
    ax.set_ylabel(f"-log10 {plab}")
    ax.set_title(f"{title}\nfold >= {2 ** lfc_min:.2f} (non-log; |log2FC| >= {lfc_min:.2f}), "
                 f"{plab} < {alpha}{note}", fontsize=8)
    ax.margins(x=0.18)
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return d.loc[d["sig"], "gene"].astype(str).tolist()


def _quiet_render_loggers():
    """Silence the fontTools/matplotlib per-glyph subset chatter emitted while writing many PDFs so a
    batch re-render leaves a readable progress log."""
    for n in ("fontTools", "fontTools.subset", "fontTools.ttLib", "matplotlib",
              "matplotlib.font_manager", "PIL"):
        logging.getLogger(n).setLevel(logging.ERROR)


def rerender_volcanoes(root, use_rawp=True, recurse=True, summary_csv=None, quiet_font_logs=True):
    """Re-render the training (panel 6) and predicted (panel 7) volcano PDFs IN PLACE from the
    already-saved `*_6_training_DEGs.tsv` / `*_7_predicted_DEGs.tsv` tables, using the CURRENT
    thresholds (VOLCANO_LOG2FC = log2(1.2), a non-log 1.2-fold cutoff) and p convention. The DEG
    statistics (log2fc/pval/fdr) are read straight from disk -- no re-genotyping or re-classification
    -- so this is the fast path to refresh volcanoes after a threshold / p-value change, reproducing
    the exact panels emitted by process_level (training top-10 up/down by fold; predicted labels the
    training DEGs). The predicted table's `training_DEG` flag is refreshed to the recomputed set.

    `recurse=True` walks every `*_6_training_DEGs.tsv` under `root`, so pointing it at a whole project
    tree refreshes every sample/variant/cell-state volcano in one call. When `summary_csv` is given, a
    per-PDF summary (relative path, panel, genes tested, genes significant) is written there.
    Returns a list of dicts (one per PDF written)."""
    import glob
    if quiet_font_logs:
        _quiet_render_loggers()
    root = os.path.abspath(root)
    suffix = "_6_training_DEGs.tsv"
    pat = os.path.join(root, "**", "*" + suffix) if recurse else os.path.join(root, "*" + suffix)
    written = []
    for train_tsv in sorted(glob.glob(pat, recursive=recurse)):
        d_dir = os.path.dirname(train_tsv)
        tag = os.path.basename(train_tsv)[:-len(suffix)]
        de_train = pd.read_csv(train_tsv, sep="\t")
        out6 = os.path.join(d_dir, f"{tag}_6_training_volcano.pdf")
        train_degs = _volcano(de_train, out6, f"{tag}: training MUT vs WT", use_rawp=use_rawp)
        written.append({"tag": tag, "panel": "6_training", "pdf": out6,
                        "pdf_rel": os.path.relpath(out6, root),
                        "n_tested": int(len(de_train)), "n_sig": int(len(train_degs))})
        logger.info("[revolcano] %s training -> %s (%d significant)", tag, out6, len(train_degs))
        pred_tsv = os.path.join(d_dir, f"{tag}_7_predicted_DEGs.tsv")
        if os.path.exists(pred_tsv):
            de_pred = pd.read_csv(pred_tsv, sep="\t")
            de_pred["training_DEG"] = de_pred["gene"].astype(str).isin(set(train_degs))
            de_pred.to_csv(pred_tsv, sep="\t", index=False)   # keep the flag consistent with the new DEG set
            out7 = os.path.join(d_dir, f"{tag}_7_predicted_volcano.pdf")
            _volcano(de_pred, out7, f"{tag}: predicted MUT vs WT  (labels = training DEGs)",
                     label_genes=train_degs, use_rawp=use_rawp)
            written.append({"tag": tag, "panel": "7_predicted", "pdf": out7,
                            "pdf_rel": os.path.relpath(out7, root),
                            "n_tested": int(len(de_pred)), "n_sig": int(len(train_degs))})
            logger.info("[revolcano] %s predicted -> %s", tag, out7)
    if summary_csv and written:
        os.makedirs(os.path.dirname(os.path.abspath(summary_csv)), exist_ok=True)
        pd.DataFrame(written).to_csv(summary_csv, index=False)
        logger.info("[revolcano] wrote summary of %d PDFs -> %s", len(written), summary_csv)
    return written


def _render_combined(sub, grp, genes, gene_dirs, dir_groups, out_pdf, out_tsv):
    """Four-group combined heatmap over a FIXED gene list, preserving the caller's gene ORDER, rendered
    with the same MarkerFinder heatmap code (_build_heatmap/_plot_heatmap). `dir_groups` = (group for
    MUT-marker genes, group for WT-marker genes) used for the row colour bar."""
    try:
        hm, _, counts, ordered = _build_heatmap(sub, "vi_group", genes, GROUP_ORDER4, False, None)
    except Exception as e:
        logger.warning("combined heatmap failed (%s): %s", os.path.basename(out_pdf), e)
        return None, None
    if hm is None or hm.empty:
        return None, None
    col_clusters = [grp[b] for b in ordered]
    row_clusters = [dir_groups[0] if gene_dirs.get(str(gene)) == "MUT" else dir_groups[1]
                    for gene in hm.index]
    try:
        _plot_heatmap(hm, out_pdf, counts, GROUP_ORDER4, col_clusters, row_clusters)
    except Exception as e:
        logger.warning("combined heatmap render failed (%s): %s", os.path.basename(out_pdf), e)
        return None, None
    hm.to_csv(out_tsv, sep="\t")
    return hm, col_clusters


def _predicted_and_combined_panels(base, mut_cells, wt_cells, imp, markers, out_dir, tag,
                                   n_markers, species=None):
    """Emit panels 2-5 (see above). `markers` is the TRAINING MarkerFinder table (feature, direction),
    whose row order defines the gene order for panel 3; panel 4 repeats the plot with the PREDICTED-cell
    genes from panel 2, in the predicted gene order."""
    pred_mut = [b for b in imp.loc[imp["imputed"] == "MUT", "barcode"].tolist() if b in base.obs_names]
    pred_wt = [b for b in imp.loc[imp["imputed"] == "WT", "barcode"].tolist() if b in base.obs_names]

    # --- panel 2: MarkerFinder on the NEW cells, markers derived from those cells ---
    pred_markers = None
    if len(pred_mut) >= 3 and len(pred_wt) >= 3:
        pred_markers = _markerfinder(base, pred_mut + pred_wt,
                                     ["MUT"] * len(pred_mut) + ["WT"] * len(pred_wt),
                                     out_dir, f"{tag}_2_predicted", n_markers, species=species)

    # --- shared four-group cell set ---
    grp = {}
    for b in mut_cells:
        grp[b] = G_TRAIN_MUT
    for b in wt_cells:
        grp[b] = G_TRAIN_WT
    for b in pred_mut:
        grp[b] = G_PRED_MUT
    for b in pred_wt:
        grp[b] = G_PRED_WT
    if not grp:
        return
    sub = base[base.obs_names.isin(list(grp))].copy()
    sub.obs["vi_group"] = pd.Series(grp).reindex(sub.obs_names).values
    normalize_adata(sub)

    # --- panel 3: combined, TRAINING genes in TRAINING order ---
    hm3, col3 = _render_combined(
        sub, grp, [str(f) for f in markers["feature"]],
        dict(zip(markers["feature"].astype(str), markers["direction"].astype(str))),
        (G_TRAIN_MUT, G_TRAIN_WT),
        os.path.join(out_dir, f"{tag}_3_combined_marker_heatmap.pdf"),
        os.path.join(out_dir, f"{tag}_3_combined_heatmap_matrix.tsv"))

    # --- panel 4: same plot, PREDICTED-cell genes from panel 2, in the predicted gene order ---
    if pred_markers is not None and len(pred_markers):
        _render_combined(
            sub, grp, [str(f) for f in pred_markers["feature"]],
            dict(zip(pred_markers["feature"].astype(str), pred_markers["direction"].astype(str))),
            (G_PRED_MUT, G_PRED_WT),
            os.path.join(out_dir, f"{tag}_4_combined_predictedGenes_marker_heatmap.pdf"),
            os.path.join(out_dir, f"{tag}_4_combined_predictedGenes_heatmap_matrix.tsv"))

    # --- panel 5: PCA of the panel-3 matrix ---
    if hm3 is not None:
        _pca_from_heatmap(hm3, col3, GROUP_ORDER4,
                          os.path.join(out_dir, f"{tag}_5_combined_pca.pdf"),
                          title=f"{tag}: training + imputed (training marker genes)")

    # --- panel 6: TRAINING volcano (cellHarmony moderated t-test), top 10 up/down labelled by fold ---
    de_train = _moderated_de(base, mut_cells, wt_cells)
    train_degs = []
    if len(de_train):
        de_train.to_csv(os.path.join(out_dir, f"{tag}_6_training_DEGs.tsv"), sep="\t", index=False)
        train_degs = _volcano(de_train, os.path.join(out_dir, f"{tag}_6_training_volcano.pdf"),
                              f"{tag}: training MUT vs WT")

    # --- panel 7: PREDICTED volcano, labelling ONLY the training DEGs ---
    de_pred = _moderated_de(base, pred_mut, pred_wt)
    if len(de_pred):
        de_pred["training_DEG"] = de_pred["gene"].astype(str).isin(set(train_degs))
        de_pred.to_csv(os.path.join(out_dir, f"{tag}_7_predicted_DEGs.tsv"), sep="\t", index=False)
        _volcano(de_pred, os.path.join(out_dir, f"{tag}_7_predicted_volcano.pdf"),
                 f"{tag}: predicted MUT vs WT  (labels = training DEGs)",
                 label_genes=train_degs)


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
            # panels 2-4: predicted-cell MarkerFinder, 4-group combined heatmap on the training genes
            # (training order), and the PCA derived from that combined matrix
            _predicted_and_combined_panels(base, list(mut), list(wt), imp, mk, udir, tag,
                                           n_markers, species=species)
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


def subset_h5ad_cells(h5ad_path, rows, layer="counts"):
    """Raw-counts AnnData for the given global row indices of a large .h5ad, read with h5py in
    CONTIGUOUS RUNS so a multi-GB matrix is never loaded into memory. This is what lets variant
    genotypes be paired with an expression atlas that does not live beside the BAM (e.g. a combined
    short-read atlas), instead of the `<bam_dir>/<library>-<level>.h5ad` convention.
    obs carries every column of the source obs for the selected cells."""
    import h5py
    from scipy.sparse import csr_matrix
    rows = np.sort(np.asarray(rows))
    if rows.size == 0:
        raise ValueError("subset_h5ad_cells: no rows requested")
    runs, s = [], 0
    for i in range(1, len(rows) + 1):
        if i == len(rows) or rows[i] != rows[i - 1] + 1:
            runs.append((rows[s], rows[i - 1])); s = i
    with h5py.File(h5ad_path, "r") as f:
        X = f["layers"][layer] if (layer and "layers" in f and layer in f["layers"]) else f["X"]
        indptr = X["indptr"][:]
        var_names = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]])
        obs_index = np.array([x.decode() if isinstance(x, bytes) else x for x in f["obs"]["_index"][:]])
        data_l, ind_l, rowlens = [], [], []
        for r0, r1 in runs:
            a, b = int(indptr[r0]), int(indptr[r1 + 1])
            data_l.append(X["data"][a:b]); ind_l.append(X["indices"][a:b])
            rowlens.append(np.diff(indptr[r0:r1 + 2]))
        obs = {}
        for k in f["obs"].keys():
            if k in ("_index",):
                continue
            g = f["obs"][k]
            try:
                if isinstance(g, h5py.Group) and "categories" in g and "codes" in g:
                    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]])
                    codes = g["codes"][:][rows]
                    obs[k] = np.where(codes >= 0, cats[np.clip(codes, 0, None)], "NA")
                else:
                    v = g[:][rows]
                    obs[k] = np.array([x.decode() if isinstance(x, bytes) else x for x in v])
            except Exception:
                continue
    M = csr_matrix((np.concatenate(data_l), np.concatenate(ind_l),
                    np.concatenate([[0], np.cumsum(np.concatenate(rowlens))])),
                   shape=(len(rows), len(var_names)))
    A = ad.AnnData(X=M, obs=pd.DataFrame(obs, index=obs_index[rows]),
                   var=pd.DataFrame(index=var_names))
    A.layers["counts"] = M.copy()
    A.var_names_make_unique()
    return A


def pooled_genotype_table(base, member_counts, group_col, mut_min=1, wt_min=1):
    """MUT/WT/UNK per FULL barcode for an AnnData that POOLS several libraries (e.g. a donor's
    timepoints). `member_counts` maps a group value (obs[group_col], one per library) to that
    library's per-cell (mut, wt) table. Full-barcode keying is required because 16 bp cell cores
    collide across libraries -- keying on the core alone would leak genotypes between samples."""
    core_arr = pd.Series([_core(b) for b in base.obs_names], index=base.obs_names)
    grp = pd.Series(base.obs[group_col].astype(str).values, index=base.obs_names)
    geno = pd.Series("UNK", index=base.obs_names)
    for member, cnt in member_counts.items():
        if cnt is None or not len(cnt):
            continue
        m = (grp == str(member)).to_numpy()
        if not m.any():
            continue
        cells = base.obs_names[m]
        cores = core_arr[m]
        # genotype tables may be keyed on either orientation; align to the h5ad cores
        idx = set(cnt.index)
        keys = cores if sum(1 for c in cores[:200] if c in idx) >= sum(
            1 for c in cores[:200] if _core(_rc(c)) in idx) else cores.map(lambda c: _core(_rc(c)))
        mv = keys.map(cnt["mut"]).fillna(0).to_numpy()
        wv = keys.map(cnt["wt"]).fillna(0).to_numpy()
        geno.loc[cells] = np.where(mv >= mut_min, "MUT", np.where(wv >= wt_min, "WT", "UNK"))
    return geno


def process_pooled(base, norm, geno, states, sample, variant, level, out_dir,
                   min_cells=15, n_markers=50, species=None):
    """`process_level` for a POOLED (multi-library) AnnData whose genotypes are already resolved per
    full barcode (see `pooled_genotype_table`). Same units, gate, imputation and 7-panel output."""
    st = pd.Series({b: states.get(b, states.get(_core(b))) for b in base.obs_names})
    summary, deployed = [], []
    for b in geno.index[geno.isin(["MUT", "WT"])]:
        deployed.append({"sample": sample, "variant": variant, "level": level, "barcode": b,
                         "cell_state": st.get(b), "genotype": geno[b], "source": "called"})
    units = [("combined", geno.index[geno.isin(["MUT", "WT"])])]
    for s in sorted(pd.Series(st).dropna().unique()):
        units.append((s, geno.index[(st == s).values & geno.isin(["MUT", "WT"]).values]))
    for unit, cells in units:
        gl = geno.loc[cells]; mut, wt = cells[gl == "MUT"], cells[gl == "WT"]
        if len(mut) < min_cells or len(wt) < min_cells:
            continue
        is_comb = unit == "combined"
        udir = os.path.join(out_dir, "combined") if is_comb else \
            os.path.join(out_dir, "per_celltype", _safe(unit))
        os.makedirs(udir, exist_ok=True)
        tag = "combined" if is_comb else _safe(unit)
        mk = _markerfinder(base, list(mut) + list(wt), ["MUT"] * len(mut) + ["WT"] * len(wt),
                           udir, tag, n_markers, species=species)
        if mk is None or not len(mk):
            continue
        sig = mk["feature"].astype(str).tolist()
        unk = geno.index[geno == "UNK"] if is_comb else \
            geno.index[(st == unit).values & (geno == "UNK").values]
        res = impute_unknown(_expr_subset(norm, list(mut) + list(wt) + list(unk), sig),
                             geno.loc[list(mut) + list(wt)], pd.Index(unk), sig, min_called=min_cells)
        if res is None:
            continue
        acc = float(res.attrs.get("concordance", 0.0))
        _write_verification(res, udir, tag)
        imp = res[res["imputed"].isin(["MUT", "WT"])] if len(res) else res
        if len(imp):
            imp.to_csv(os.path.join(udir, f"{tag}_imputed.csv"), index=False)
            _predicted_and_combined_panels(base, list(mut), list(wt), imp, mk, udir, tag,
                                           n_markers, species=species)
        n_im = int((imp["imputed"] == "MUT").sum()) if len(imp) else 0
        n_iw = int((imp["imputed"] == "WT").sum()) if len(imp) else 0
        deploy = (not is_comb) and acc >= ACCURACY_MIN
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
    return summary, deployed


def zscore_genes(expr):
    """Gene-wise z-score of a genes x samples matrix (log scale). Compute ONCE and reuse -- scoring
    thousands of permuted/random signatures against a re-z-scored matrix is the dominant cost."""
    z = expr.copy()
    z.index = pd.Index([str(g) for g in z.index])
    sd = z.std(axis=1).replace(0, np.nan)
    return z.sub(z.mean(axis=1), axis=0).div(sd, axis=0).dropna(how="all")


def score_signature(z, up_genes, down_genes=None):
    """Per-sample MUT-vs-WT signature score on a PRE-Z-SCORED matrix (see `zscore_genes`):
    mean(z[up]) - mean(z[down]). Returns a Series indexed by sample."""
    zidx = set(z.index)
    # NB: `down_genes or []` would raise on a numpy array (ambiguous truth value)
    down_genes = [] if down_genes is None else list(down_genes)
    up = [g for g in dict.fromkeys(map(str, up_genes)) if g in zidx]
    dn = [g for g in dict.fromkeys(map(str, down_genes)) if g in zidx]
    if not up and not dn:
        return None
    s = pd.Series(0.0, index=z.columns)
    if up:
        s = s + z.loc[up].mean(axis=0)
    if dn:
        s = s - z.loc[dn].mean(axis=0)
    s.attrs["n_up"], s.attrs["n_down"] = len(up), len(dn)
    return s


def auroc(score, is_positive):
    """AUROC from the Mann-Whitney U statistic (no sklearn dependency)."""
    s = np.asarray(score, dtype=float)
    y = np.asarray(is_positive).astype(bool)
    n1, n0 = int(y.sum()), int((~y).sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    r = pd.Series(s).rank().to_numpy()
    return float((r[y].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def external_signature_validation(z, is_mutant, up_genes, down_genes=None,
                                  n_perm=500, n_random=100, seed=0):
    """Test whether a single-cell-derived signature separates MUTANT from WILD-TYPE samples in an
    INDEPENDENT bulk cohort. `z` is a PRE-Z-SCORED genes x samples matrix (`zscore_genes`).

    Reports AUROC plus three nulls, because a bulk cohort will happily produce a high AUROC for the
    wrong reason (cell-type composition, blast %, batch):
      p_mwu         Mann-Whitney on the score
      p_perm        LABEL permutation -- the honest null for "does this signature track genotype"
      p_randomgenes RANDOM gene sets of the same size -- guards against any gene set scoring
    """
    from scipy import stats
    s = score_signature(z, up_genes, down_genes)
    if s is None:
        return None
    y = np.asarray(is_mutant).astype(bool)
    sv = s.to_numpy()
    a = auroc(sv, y)
    try:
        p = float(stats.mannwhitneyu(sv[y], sv[~y], alternative="two-sided").pvalue)
    except Exception:
        p = np.nan
    rng = np.random.RandomState(seed)
    perm = np.array([auroc(sv, rng.permutation(y)) for _ in range(n_perm)])
    p_perm = float((np.abs(perm - 0.5) >= abs(a - 0.5)).mean())
    n_up, n_dn = s.attrs.get("n_up", 0), s.attrs.get("n_down", 0)
    genes = np.asarray(z.index)
    rnd = []
    for _ in range(n_random):
        pick = rng.choice(genes, size=max(1, n_up + n_dn), replace=False)
        rs = score_signature(z, pick[:n_up], pick[n_up:])
        if rs is not None:
            rnd.append(auroc(rs.to_numpy(), y))
    rnd = np.array([v for v in rnd if np.isfinite(v)])
    p_rand = float((np.abs(rnd - 0.5) >= abs(a - 0.5)).mean()) if len(rnd) else np.nan
    return {"auroc": round(a, 3), "p_mwu": p, "p_perm": round(p_perm, 4),
            "p_randomgenes": round(p_rand, 4) if np.isfinite(p_rand) else np.nan,
            "n_up": n_up, "n_down": n_dn, "n_mut": int(y.sum()), "n_wt": int((~y).sum())}


def donor_map_from_panel(panel_tsv):
    """{sample_uid: donor} parsed from a variant panel whose `notes` carry `donor=<ID>` and whose
    `expected_uids` lists that donor's samples. Collapsing samples to donors is what keeps serial
    timepoints and MISTRG xenografts of the SAME patient from being counted as separate patients."""
    panel = pd.read_csv(panel_tsv, sep="\t")
    out = {}
    for _, r in panel.iterrows():
        m = re.search(r"donor=([^;]+)", str(r.get("notes", "")))
        if not m:
            continue
        for uid in str(r.get("expected_uids", "")).split(";"):
            uid = uid.strip()
            if uid:
                out[uid] = m.group(1)
    return out


def mutation_recurrence(calls, donor_map=None, uid_col="uid"):
    """How many distinct PATIENTS carry each mutation.

    Returns (per_locus, per_gene): DataFrames counting distinct donors for each (gene, chrom, pos)
    and for each gene. Samples absent from `donor_map` are treated as their own donor. Recurrence is
    what licenses a signature learned in one patient to be tested in another -- a mutation seen in a
    single donor is confounded with that donor's background and cannot be cross-validated.
    """
    df = calls.copy()
    dm = donor_map or {}
    df["donor"] = df[uid_col].map(lambda u: dm.get(str(u), str(u)))
    per_locus = (df.groupby(["gene", "chrom", "pos"])["donor"]
                 .agg(n_donors="nunique", donors=lambda s: ";".join(sorted(set(s))))
                 .reset_index().sort_values("n_donors", ascending=False))
    per_gene = (df.groupby("gene")["donor"]
                .agg(n_donors="nunique", donors=lambda s: ";".join(sorted(set(s))))
                .reset_index().sort_values("n_donors", ascending=False))
    return per_locus, per_gene


def genotype_series(base, counts, mut_min=1, wt_min=1):
    """MUT/WT/UNK per cell of `base`, resolving barcode ORIENTATION against the genotype table
    (long-read isoform BAM barcodes are the reverse-complement of short-read/h5ad barcodes)."""
    h5c = {_core(b): b for b in base.obs_names}
    vor, _ = _pick_orientation(set(h5c), counts.index)
    cnt = counts.copy()
    cnt.index = [_core(_rc(c)) if vor == "rc" else c for c in cnt.index]
    cnt = cnt[~cnt.index.duplicated()]
    core = {b: _core(b) for b in base.obs_names}
    mco = pd.Series({b: cnt["mut"].get(core[b], 0) for b in base.obs_names})
    wco = pd.Series({b: cnt["wt"].get(core[b], 0) for b in base.obs_names})
    return pd.Series(np.where(mco >= mut_min, "MUT", np.where(wco >= wt_min, "WT", "UNK")),
                     index=base.obs_names)


def markers_from_cells(base, cells, labels, top_n, mode="pearson", rho_threshold=0.0):
    """Render-free MUT/WT signature for the many cross-validation folds. Two modes,
    matching the two `generate_marker_heatmap_from_adata` marker_method options:
      mode='pearson'  (== marker_method 'markerfinder'): Pearson correlation of each
                      gene to the MUT/WT on-off indicator (cellHarmony MarkerFinder).
      mode='wilcoxon' (== marker_method 'scanpy'): scanpy rank_genes_groups Wilcoxon
                      rank-sum, top_n up-genes per group.
    `rho_threshold` (pearson mode) defaults to 0.0 so the TOP top_n markers per group
    are always returned -- at high n with a weak per-cell signal, the MarkerFinder
    default 0.2 can return ZERO markers and silently drop the whole unit.
    Returns the marker feature list (top_n per group), or None."""
    sub = base[list(cells)].copy()
    sub.obs["genotype"] = pd.Series(labels, index=list(cells)).loc[sub.obs_names].astype(str).values
    normalize_adata(sub)
    if mode in ("wilcoxon", "scanpy"):
        try:
            import scanpy as sc
            sc.tl.rank_genes_groups(sub, "genotype", method="wilcoxon", n_genes=top_n)
            names = sub.uns["rank_genes_groups"]["names"]
            feats = []
            for grp in names.dtype.names:
                feats.extend(str(g) for g in list(names[grp])[:top_n])
            return list(dict.fromkeys(feats)) or None
        except Exception as e:
            logger.warning("wilcoxon markers failed: %s", e)
            return None
    from altanalyze3.components.cellHarmony.markerFinder import find_markers_from_adata
    try:
        mo = find_markers_from_adata(sub, "genotype", n_markers=top_n, direction="up",
                                     rho_threshold=rho_threshold, write_outputs=False)
    except Exception:
        return None
    return mo.markers["marker"].astype(str).tolist()


def _balanced_acc(true, pred):
    true = np.asarray(true); pred = np.asarray(pred)
    accs = [float((pred[true == c] == c).mean()) for c in ("MUT", "WT") if (true == c).any()]
    return float(np.mean(accs)) if accs else np.nan


def heldout_cv(base, norm, mut_cells, wt_cells, top_n=50, k=3, repeats=3, min_train=6, seed=0):
    """HONEST held-out accuracy of the deployed classifier. The pipeline's 85% gate is
    RESUBSTITUTION (train == test), an optimistic ceiling; this refits the ENTIRE classifier
    (MarkerFinder signature + median centroids + threshold) on training folds only and scores the
    held-out cells, so a state that passes here carries genuinely recoverable within-state signal.
    Returns balanced accuracy pooled over held-out predictions (chance = 50)."""
    mut_cells, wt_cells = list(mut_cells), list(wt_cells)
    if min(len(mut_cells), len(wt_cells)) < k * 2:
        return {"cv_balanced_acc": np.nan, "cv_mut_recall": np.nan, "cv_wt_recall": np.nan,
                "cv_n_test": 0, "cv_note": "too_few_for_cv"}
    rng = np.random.RandomState(seed)
    true_all, pred_all = [], []
    for rep in range(repeats):
        mo, wo = rng.permutation(mut_cells), rng.permutation(wt_cells)
        mf, wf = np.array_split(mo, k), np.array_split(wo, k)
        for fi in range(k):
            te = list(mf[fi]) + list(wf[fi])
            tr_m = [c for j in range(k) if j != fi for c in mf[j]]
            tr_w = [c for j in range(k) if j != fi for c in wf[j]]
            if min(len(tr_m), len(tr_w)) < min_train:
                continue
            sig = markers_from_cells(base, tr_m + tr_w, ["MUT"] * len(tr_m) + ["WT"] * len(tr_w), top_n)
            if not sig:
                continue
            geno_tr = pd.Series(["MUT"] * len(tr_m) + ["WT"] * len(tr_w), index=tr_m + tr_w)
            res = impute_unknown(_expr_subset(norm, tr_m + tr_w + te, sig), geno_tr,
                                 pd.Index(te), sig, min_called=min_train)
            if res is None or not len(res):
                continue
            sc = dict(zip(res["barcode"], res["score"]))
            truth = {c: "MUT" for c in mf[fi]}; truth.update({c: "WT" for c in wf[fi]})
            for c in te:
                if c in sc and not np.isnan(sc[c]):
                    true_all.append(truth[c]); pred_all.append("MUT" if sc[c] > 0 else "WT")
    if not true_all:
        return {"cv_balanced_acc": np.nan, "cv_mut_recall": np.nan, "cv_wt_recall": np.nan,
                "cv_n_test": 0, "cv_note": "no_folds"}
    ta, pa = np.array(true_all), np.array(pred_all)
    return {"cv_balanced_acc": round(100 * _balanced_acc(ta, pa), 1),
            "cv_mut_recall": round(100 * float((pa[ta == "MUT"] == "MUT").mean()), 1) if (ta == "MUT").any() else np.nan,
            "cv_wt_recall": round(100 * float((pa[ta == "WT"] == "WT").mean()), 1) if (ta == "WT").any() else np.nan,
            "cv_n_test": int(len(ta)), "cv_note": f"{repeats}x{k}fold"}


def cross_state_transfer(base, norm, geno, states, train_states, test_states, top_n=50, min_cells=15):
    """Train the MUT/WT classifier on ALL called cells of each train state and apply it to each other
    state's called cells: balanced accuracy vs the read genotype. Off-diagonal = honest cross-cell-state
    generalization (the test cells were never seen in training)."""
    st = pd.Series({b: states.get(_core(b), states.get(b)) for b in base.obs_names})

    def called(s):
        cells = geno.index[(st == s).values & geno.isin(["MUT", "WT"]).values]
        return cells[geno.loc[cells] == "MUT"], cells[geno.loc[cells] == "WT"]
    rows = []
    for a in train_states:
        am, aw = called(a)
        if min(len(am), len(aw)) < min_cells:
            continue
        sig = markers_from_cells(base, list(am) + list(aw), ["MUT"] * len(am) + ["WT"] * len(aw), top_n)
        if not sig:
            continue
        geno_tr = pd.Series(["MUT"] * len(am) + ["WT"] * len(aw), index=list(am) + list(aw))
        for b in test_states:
            if b == a:
                continue
            bm, bw = called(b)
            if min(len(bm), len(bw)) < min_cells:
                continue
            te = list(bm) + list(bw)
            res = impute_unknown(_expr_subset(norm, list(am) + list(aw) + te, sig), geno_tr,
                                 pd.Index(te), sig, min_called=min_cells)
            if res is None or not len(res):
                continue
            sc = dict(zip(res["barcode"], res["score"]))
            truth = {c: "MUT" for c in bm}; truth.update({c: "WT" for c in bw})
            tr = np.array([truth[c] for c in te if c in sc])
            pr = np.array(["MUT" if sc[c] > 0 else "WT" for c in te if c in sc])
            if not len(tr):
                continue
            rows.append({"train_state": a, "test_state": b,
                         "n_train_MUT": len(am), "n_train_WT": len(aw),
                         "n_test_MUT": int((tr == "MUT").sum()), "n_test_WT": int((tr == "WT").sum()),
                         "balanced_acc": round(100 * _balanced_acc(tr, pr), 1)})
    return rows


def run_impact_for_matrix(base, states, counts, sample, variant, out_dir, level="gene",
                          mut_min=1, wt_min=1, min_cells=15, top_n=50, compute_cv=False,
                          species=None):
    """Run the full per-cell-state variant-impact analysis for ONE library using an expression matrix
    supplied by the caller -- i.e. an atlas that does NOT live beside the BAM (`subset_h5ad_cells`),
    with genotypes from either `genotype_from_bam` or `genotype_from_variant_extraction`.
    Emits the 7 standard panels per unit and optionally the honest held-out CV."""
    os.makedirs(out_dir, exist_ok=True)
    norm = base.copy(); normalize_adata(norm)
    summ, dep = process_level(base, norm, states, counts, sample, variant, level, out_dir,
                              mut_min=mut_min, wt_min=wt_min, min_cells=min_cells, n_markers=top_n,
                              species=species)
    if compute_cv and summ:
        g = genotype_series(base, counts, mut_min, wt_min)
        st = pd.Series({b: states.get(_core(b)) for b in base.obs_names})
        for row in summ:
            unit = row["unit"]
            cells = g.index[g.isin(["MUT", "WT"])] if unit == "combined" else \
                g.index[(st == unit).values & g.isin(["MUT", "WT"]).values]
            gl = g.loc[cells]
            row.update(heldout_cv(base, norm, cells[gl == "MUT"], cells[gl == "WT"], top_n))
    return summ, dep


def _parse_variants(spec):
    """Accept 'chr:pos:label[:type[:ref:alt]]' (semicolon-separated) OR a file of
    'chrom start end label type [ref] [alt]'. Returns (chrom, pos, label, vtype, ref, alt); ref/alt are
    None unless supplied. They are REQUIRED for --supervised (which confirms the exact caller allele)
    and IGNORED by the default discovery genotyper."""
    out = []
    if spec and os.path.exists(spec):
        with open(spec) as fh:
            for line in fh:
                p = line.split()
                if len(p) >= 2 and not line.startswith("#"):
                    chrom = p[0] if p[0].startswith("chr") else "chr" + p[0]
                    out.append((chrom, int(p[1]), p[3] if len(p) > 3 else f"{chrom}_{p[1]}",
                                p[4] if len(p) > 4 else "SNV",
                                p[5] if len(p) > 5 else None, p[6] if len(p) > 6 else None))
    else:
        for tok in (spec or "").split(";"):
            tok = tok.strip()
            if not tok:
                continue
            parts = tok.split(":")
            chrom = parts[0] if parts[0].startswith("chr") else "chr" + parts[0]
            out.append((chrom, int(parts[1]), parts[2] if len(parts) > 2 else f"{chrom}_{parts[1]}",
                        parts[3] if len(parts) > 3 else "SNV",
                        parts[4] if len(parts) > 4 else None, parts[5] if len(parts) > 5 else None))
    return out


_DISC_SUBPATHS = ("results/merged/high_confidence_final.tsv",
                  "results/merged/medium_confidence_final.tsv",
                  "results/merged/low_confidence_final.tsv",
                  "annotation/gatk/non_germline_final.tsv")


def _read_discovery_variants(discovery_outdir, subpaths=_DISC_SUBPATHS):
    """Every variant the four callers reported for ONE sample, with their BAM evidence:
    dict(chrom,pos,ref,alt,vtype,callers,caller_count,bam_alt). Deduplicated by (chrom,pos,ref,alt)."""
    import csv

    def _num(x):
        try:
            return float(x)
        except Exception:
            return None
    seen = set()
    out = []
    for sp in subpaths:
        p = os.path.join(discovery_outdir, sp)
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
                vt = "SNV" if (lr == 1 and la == 1) else "Indel"
                out.append(dict(chrom=ch, pos=int(pos), ref=ref, alt=alt, vtype=vt,
                                callers=row.get("CALLERS", ""), caller_count=row.get("CALLER_COUNT", ""),
                                bam_alt=_num(row.get("BAM_ALT"))))
    return out


def validate_supervised(bam_path, discovery_outdir, out_dir, min_mapq=1, indel_window=5, vtypes=None,
                        limit=None):
    """Validate the FULLY SUPERVISED genotyper against the callers' own BAM evidence, for EVERY variant
    the four algorithms reported in THIS ONE SAMPLE (no other sample, no discovery). For each variant it
    runs genotype_from_bam_supervised and compares MUT-read detection to BAM_ALT. Writes:
      supervised_detection_per_variant.tsv  (all variants, detected flag)
      supervised_detection_failures.tsv     (ONLY BAM_ALT>0 variants the genotyper misses -- the point)
      supervised_detection_summary.tsv      (per-type detection rate on BAM_ALT>0 variants)
    Detection rate is measured only on truly-present variants (BAM_ALT>0)."""
    import csv
    import collections
    os.makedirs(out_dir, exist_ok=True)
    variants = _read_discovery_variants(discovery_outdir)
    if vtypes:
        vset = {v.lower() for v in vtypes}
        variants = [v for v in variants if v["vtype"].lower() in vset]
    if limit and len(variants) > limit:
        step = max(1, len(variants) // limit)           # even, deterministic spread across the panel
        variants = variants[::step][:limit]
    rows = []
    pooled = collections.defaultdict(lambda: [0, 0])       # vtype -> [detected, present]
    for v in variants:
        df = genotype_from_bam_supervised(bam_path, v["chrom"], v["pos"], v["ref"], v["alt"],
                                          min_mapq=min_mapq, indel_window=indel_window)
        gm = int(df["mut"].sum()) if len(df) else 0
        present = bool(v["bam_alt"] and v["bam_alt"] > 0)
        det = int(gm > 0)
        rows.append(dict(chrom=v["chrom"], pos=v["pos"], ref=v["ref"], alt=v["alt"], vtype=v["vtype"],
                         caller_count=v["caller_count"], callers=v["callers"], bam_alt=v["bam_alt"],
                         geno_mut_reads=gm, detected=det, present=int(present)))
        if present:
            pooled[v["vtype"]][0] += det; pooled[v["vtype"]][1] += 1
    cols = ["chrom", "pos", "ref", "alt", "vtype", "caller_count", "callers", "bam_alt",
            "geno_mut_reads", "detected", "present"]
    with open(os.path.join(out_dir, "supervised_detection_per_variant.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader()
        for r in rows:
            w.writerow(r)
    fails = [r for r in rows if r["present"] and not r["detected"]]
    fcols = ["chrom", "pos", "ref", "alt", "vtype", "caller_count", "callers", "bam_alt", "geno_mut_reads"]
    with open(os.path.join(out_dir, "supervised_detection_failures.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fcols, delimiter="\t"); w.writeheader()
        for r in fails:
            w.writerow({k: r[k] for k in fcols})
    sp = os.path.join(out_dir, "supervised_detection_summary.tsv")
    with open(sp, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["vtype", "detected", "present", "detection_rate"])
        td = tp = 0
        for vt, (d, n) in sorted(pooled.items()):
            w.writerow([vt, d, n, round(d / n, 4) if n else ""]); td += d; tp += n
        w.writerow(["ALL", td, tp, round(td / tp, 4) if tp else ""])
    logger.info("supervised validation: %d variants; present(BAM_ALT>0)=%d detected=%d (%.4f); "
                "failures=%d -> %s", len(rows), tp, td, (td / tp if tp else 0.0), len(fails), out_dir)
    return sp


# ============================================================================
# PASSENGER -> DRIVER ASSOCIATION  (supervised genotyping here + mt_variants statistics)
# ----------------------------------------------------------------------------
def _passenger_panel_from_discovery(discovery_outdir, exclude_keys=(), subpaths=_DISC_SUBPATHS):
    """Nuclear + mitochondrial passenger markers = UNION of the four-caller tiers (high+medium+low) and
    GATK non_germline for ONE sample. Returns list of dict(name, chrom, pos, ref, alt). Driver
    coordinates (exclude_keys, 'chrom:pos:ref:alt') are removed. No discovery -- reads existing TSVs."""
    import csv
    exclude = set(exclude_keys)
    seen = set()
    out = []
    for sp in subpaths:
        p = os.path.join(discovery_outdir, sp)
        if not os.path.isfile(p):
            continue
        with open(p) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                ch = row.get("CHROM"); pos = row.get("POS")
                ref = row.get("REF"); alt = row.get("ALT"); gene = row.get("GENE", "") or ""
                if not (ch and pos and ref and alt):
                    continue
                key = f"{ch}:{pos}:{ref}:{alt}"
                if key in exclude or key in seen:
                    continue
                seen.add(key)
                out.append(dict(name=f"{gene}_{ch}:{pos}", chrom=ch, pos=int(pos), ref=ref, alt=alt))
    return out


def _driver_panel(driver_tsv):
    """Driver (pathogenic) panel: name<TAB>chrom<TAB>pos<TAB>ref<TAB>alt. Returns (list of dicts, set of
    'chrom:pos:ref:alt' coord keys to exclude from the passenger panel)."""
    drv = []
    keys = set()
    with open(driver_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: header.index(c) for c in ("name", "chrom", "pos", "ref", "alt")}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(idx.values()):
                continue
            drv.append(dict(name=p[idx["name"]], chrom=p[idx["chrom"]], pos=int(p[idx["pos"]]),
                            ref=p[idx["ref"]], alt=p[idx["alt"]]))
            keys.add(f'{p[idx["chrom"]]}:{p[idx["pos"]]}:{p[idx["ref"]]}:{p[idx["alt"]]}')
    return drv, keys


def build_supervised_feature_matrix(bam_path, uid, driver_variants, marker_variants, out_dir,
                                    min_mapq=1, indel_window=5):
    """cell x (driver + passenger-marker) matrix, genotyped with the SUPERVISED exact-allele genotyper
    (genotype_from_bam_supervised), in the EXACT h5ad layout mt_variants.nominate_presence_fdr consumes:
    layers['alt']/['total'], var['feature_type'] in {driver, mtDNA}, obs index 'uid|barcode'. Co-located
    driver alleles (same name) collapse to ONE feature (mut if ANY allele present). Writes
    {uid}_cell_feature_matrix.h5ad and returns its path."""
    import anndata as _ad
    import scipy.sparse as _sp

    def _group(vs):
        g = {}
        for v in vs:
            g.setdefault(v["name"], dict(chrom=v["chrom"], pos=int(v["pos"]), ref=v["ref"],
                                         alt=v["alt"], variants=[]))["variants"].append(v)
        return g
    drv = _group(driver_variants)
    mrk = _group(marker_variants)
    feat_cells = {}
    for name, info in list(drv.items()) + list(mrk.items()):
        perbc = {}
        for v in info["variants"]:
            df = genotype_from_bam_supervised(bam_path, v["chrom"], int(v["pos"]), v["ref"], v["alt"],
                                              min_mapq=min_mapq, indel_window=indel_window)
            for core, row in df.iterrows():
                a = int(row["mut"]); t = int(row["mut"] + row["wt"])
                pa, pt = perbc.get(core, (0, 0))
                perbc[core] = (max(pa, a), max(pt, t))   # collapse co-located by max within a cell
        feat_cells[name] = perbc
    features = [(nm, "driver", d["chrom"], d["pos"], d["ref"], d["alt"]) for nm, d in drv.items()]
    features += [(nm, "mtDNA", d["chrom"], d["pos"], d["ref"], d["alt"]) for nm, d in mrk.items()]
    names = [f[0] for f in features]
    barcodes = sorted({c for nm in names for c in feat_cells.get(nm, {})})
    bc_idx = {b: i for i, b in enumerate(barcodes)}
    alt = np.zeros((len(barcodes), len(names)), dtype=np.int32)
    tot = np.zeros((len(barcodes), len(names)), dtype=np.int32)
    for j, nm in enumerate(names):
        for c, (a, t) in feat_cells.get(nm, {}).items():
            i = bc_idx[c]; alt[i, j] = a; tot[i, j] = t
    frac = np.divide(alt, tot, out=np.zeros(alt.shape, float), where=tot > 0)
    var = pd.DataFrame({"feature_type": [f[1] for f in features], "chrom": [f[2] for f in features],
                        "pos": [f[3] for f in features], "ref": [f[4] for f in features],
                        "alt": [f[5] for f in features]}, index=names)
    obs = pd.DataFrame(index=[f"{uid}|{b}" for b in barcodes])
    adata = _ad.AnnData(X=_sp.csr_matrix(frac), obs=obs, var=var)
    adata.layers["alt"] = _sp.csr_matrix(alt)
    adata.layers["total"] = _sp.csr_matrix(tot)
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"{uid}_cell_feature_matrix.h5ad")
    adata.write_h5ad(out, compression="gzip")
    logger.info("[%s] supervised matrix: %d cells x %d features (%d driver + %d passenger) -> %s",
                uid, adata.n_obs, adata.n_vars, len(drv), len(mrk), out)
    return out


def _matrix_marker_ref_alt(matrix_dir):
    """feature name -> (ref, alt) read from the built matrices' var, for snv/indel type + full display."""
    import anndata as _ad
    import glob as _glob
    ra = {}
    for h in _glob.glob(os.path.join(matrix_dir, "*_cell_feature_matrix.h5ad")):
        v = _ad.read_h5ad(h, backed="r").var
        for nm, rf, al in zip(v.index, v["ref"], v["alt"]):
            ra.setdefault(str(nm), (str(rf), str(al)))
    return ra


def _write_supervised_hits_table(annotated_paths, matrix_dir, out_path):
    """Image-format hits table (the fields you specified): level, gene, driver, unit, n_mut, marker,
    type, recall, sens, spec, prec, bg, #mut-Var, #wt-Var, #wt, #Var-pred, OR, p. marker gets ref>alt
    appended; type is snv/indel from the exact allele. Sample + patient rows combined, recall-sorted."""
    import csv
    ra = _matrix_marker_ref_alt(matrix_dir)
    cols = ["level", "gene", "driver", "unit", "n_mut", "marker", "type", "recall", "sens", "spec",
            "prec", "bg", "#mut-Var", "#wt-Var", "#wt", "#Var-pred", "OR", "p"]
    rows = []
    for path in annotated_paths:
        if not os.path.isfile(path):
            continue
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                name = r.get("mt_feature", "")
                rf, al = ra.get(name, ("", ""))
                typ = "indel" if (rf and al and len(rf) != len(al)) else "snv"
                marker = f"{name}_{rf}>{al}" if rf else name
                rows.append({"level": r.get("level", ""), "gene": r.get("gene", ""),
                             "driver": r.get("driver", ""), "unit": r.get("unit", ""),
                             "n_mut": r.get("n_mut", ""), "marker": marker, "type": typ,
                             "recall": r.get("recall_all", ""), "sens": r.get("sensitivity", ""),
                             "spec": r.get("specificity", ""), "prec": r.get("precision", ""),
                             "bg": r.get("background", ""), "#mut-Var": r.get("mutant_carriers", ""),
                             "#wt-Var": r.get("wt_carriers", ""), "#wt": r.get("n_wt_covered", ""),
                             "#Var-pred": r.get("n_mt_predicted", ""), "OR": r.get("odds_ratio", ""),
                             "p": r.get("fisher_p", "")})

    def _rk(x):
        # rank by SENSITIVITY (primary). Specificity/background is unreliable in single-cell data: an
        # undetected mutation (allelic dropout) makes a mutant cell look wild-type, so it inflates the
        # "wt-carrying-marker" count -- a high-sensitivity passenger is a good surrogate regardless.
        try:
            return (-float(x["sens"] or 0), -float(x["recall"] or 0))
        except Exception:
            return (0.0, 0.0)
    rows.sort(key=_rk)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    logger.info("wrote %s (%d hits)", out_path, len(rows))
    return rows


def assess_extraction_recovery(discovery_outdir, extraction_readcounts_tsv, out_path,
                               subpaths=_DISC_SUBPATHS):
    """Step 7: how many / what percent of the four-algorithm variants (union) for this sample does
    variant_extract.py (results_variant_extraction) capture. Compared by chrom:pos. Writes a summary."""
    import csv
    disc = set()
    for sp in subpaths:
        p = os.path.join(discovery_outdir, sp)
        if not os.path.isfile(p):
            continue
        with open(p) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("CHROM") and row.get("POS"):
                    disc.add((row["CHROM"], str(row["POS"])))
    ext = set()
    if os.path.isfile(extraction_readcounts_tsv):
        with open(extraction_readcounts_tsv) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                ch = row.get("chrom"); ps = row.get("pos")
                if ch and ps:
                    ext.add((ch, str(ps)))
    total = len(disc)
    captured = len(disc & ext)
    missed = total - captured
    with open(out_path, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"four_caller_variants_union\t{total}\n")
        fh.write(f"captured_by_variant_extract\t{captured}\n")
        fh.write(f"missed_by_variant_extract\t{missed}\n")
        fh.write(f"pct_missed\t{round(100 * missed / total, 2) if total else ''}\n")
    logger.info("[recovery] %d/%d four-caller variants missed by variant_extract (%.2f%%) -> %s",
                missed, total, (100 * missed / total if total else 0.0), out_path)
    return out_path


def run_supervised_passenger_association(args):
    """END TO END, supervised: for each sample build the cell x (driver + four-caller-union passenger)
    matrix with the EXACT-allele genotyper (this module), then run the mt_variants.py statistics
    (nominate_presence_fdr + annotate_hits_2x2, sample + patient) and write the image-format hits table.
    Also writes the per-sample variant_extract recovery (step 7). Re-runs NONE of the four callers."""
    from altanalyze3.components.bam import mt_variants as mtv
    metadata = pd.read_csv(args.metadata, sep="\t")
    out_root = os.path.abspath(getattr(args, "output", None) or os.path.join(os.getcwd(),
                                                                             "passenger_association"))
    matrix_dir = os.path.join(out_root, "matrices")
    os.makedirs(matrix_dir, exist_ok=True)
    drv, drv_keys = _driver_panel(args.driver_variants)
    only = set(getattr(args, "samples", None) or [])
    mapq = getattr(args, "min_mapq", None) or 1
    iwin = getattr(args, "indel_window", 5)
    rve = getattr(args, "results_extraction", None)
    built = []
    if not getattr(args, "nominate_only", False):
        for _, r in metadata.iterrows():
            uid = str(r.get("uid", "")).strip()
            bam = str(r.get("bam", ""))
            if not uid or (only and uid not in only):
                continue
            if not os.path.exists(bam):
                logger.warning("[%s] bam missing: %s", uid, bam); continue
            disc = os.path.join(os.path.dirname(bam), "genome_variant_detection")
            if not os.path.isdir(disc):
                logger.warning("[%s] no discovery outdir: %s", uid, disc); continue
            markers = _passenger_panel_from_discovery(disc, exclude_keys=drv_keys)
            lim = getattr(args, "limit_markers", None)
            if lim:
                markers = markers[:int(lim)]          # test-scope only (leave unset for production)
            build_supervised_feature_matrix(bam, uid, drv, markers, matrix_dir, min_mapq=mapq,
                                            indel_window=iwin)
            built.append(uid)
            if rve:
                assess_extraction_recovery(disc, os.path.join(rve, f"{uid}_variant_readcounts.tsv"),
                                           os.path.join(out_root, f"{uid}_extraction_recovery.tsv"))
    if getattr(args, "build_only", False):
        logger.info("build-only: %d matrices -> %s", len(built), matrix_dir); return
    import glob as _glob
    samples = [os.path.basename(h).replace("_cell_feature_matrix.h5ad", "")
               for h in _glob.glob(os.path.join(matrix_dir, "*_cell_feature_matrix.h5ad"))]
    if not samples:
        raise SystemExit("no matrices to nominate on: %s" % matrix_dir)
    min_mut = getattr(args, "nominate_min_mut", None) or 10
    max_q = getattr(args, "max_q", 0.05)
    max_bg = getattr(args, "max_background", 0.15)
    mtv.driver_genotype_summary(matrix_dir, samples, out_root)
    mtv.nominate_presence_fdr(matrix_dir, {s: [s] for s in samples}, out_root, level="sample",
                              min_mut=min_mut, max_q=max_q, max_background=max_bg)
    pmap = None
    if getattr(args, "patient_map", None):
        smp = set(samples); pmap = {}
        with open(args.patient_map) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2 and parts[0] in smp:
                    pmap.setdefault(parts[1], []).append(parts[0])
        mtv.nominate_presence_fdr(matrix_dir, pmap, out_root, level="patient",
                                  min_mut=min_mut, max_q=max_q, max_background=max_bg)
    annotated = []
    for lvl in ("sample", "patient"):
        src = os.path.join(out_root, f"mt_nomination_presence_{lvl}.tsv")
        if not os.path.isfile(src):
            continue
        dst = os.path.join(out_root, f"passenger_nomination_presence_{lvl}.tsv")
        os.replace(src, dst)
        mtv.annotate_hits_2x2(matrix_dir, dst, out_root, patient_map=(pmap if lvl == "patient" else None))
        ann = os.path.join(out_root, f"passenger_nomination_presence_{lvl}_annotated.tsv")
        if os.path.isfile(ann):
            annotated.append(ann)
    _write_supervised_hits_table(annotated, matrix_dir,
                                 os.path.join(out_root, "passenger_marker_hits_ALL.tsv"))
    logger.info("passenger association done -> %s", out_root)


def _discovery_outdirs_by_uid(manifest_path):
    """uid -> [discovery outdirs] from the genomic-discovery manifest (columns include uid, outdir)."""
    import csv as _csv
    by_uid = {}
    with open(manifest_path) as fh:
        for row in _csv.DictReader(fh, delimiter="\t"):
            u = row.get("uid"); od = row.get("outdir")
            if not u:
                continue
            by_uid.setdefault(u, [])
            if od and od not in by_uid[u]:
                by_uid[u].append(od)
    return by_uid


_FUNCTIONAL = ("missense_variant", "frameshift_variant", "stop_gained", "stop_lost", "start_lost",
               "splice_acceptor_variant", "splice_donor_variant", "inframe_insertion",
               "inframe_deletion", "protein_altering_variant")

# artifact-prone gene classes (KINNEX-RNA multi-mapping / hypervariable / ubiquitously recurrent):
# mitochondrial, HLA/MHC, immunoglobulin/T-cell-receptor, ribosomal proteins + rRNA, histones, snRNA.
_ARTIFACT_GENE_PREFIX = ("MT-", "HLA-", "IGH", "IGK", "IGL", "IGLV", "TRAV", "TRBV", "TRGV", "TRDV",
                         "RPL", "RPS", "MRPL", "MRPS", "RNA5", "RNA18", "RNA28", "RNA45", "Y_RNA",
                         "HIST", "H1-", "H2A", "H2B", "H3-", "H4-", "RNU", "U1", "U2", "U3", "U4",
                         "U5", "U6", "U7", "U8", "SNOR", "SCARNA")


def nominate_functional_drivers(discovery_manifest, out_path, subpaths=_DISC_SUBPATHS, max_gnomad=0.001,
                                gene_whitelist=None, require_cosmic_clinvar=False, priority_genes=None):
    """Nominate ADDITIONAL candidate driver mutations from ALL four-caller variants (no caller re-run):
    read each sample's discovery tier tables, keep variants with a protein-altering CONSEQUENCE that are
    somatic (SOMATIC_CLASS not germline) and rare (gnomAD_AF_max <= max_gnomad), and rank by RECURRENCE
    across samples then COSMIC/ClinVar support. Writes candidate_drivers.tsv (gene, chrom, pos, ref, alt,
    consequence, HGVSp, cosmic, clinvar, gnomAD, n_samples, samples)."""
    import csv as _csv
    by_uid = _discovery_outdirs_by_uid(discovery_manifest)

    def _num(x):
        try:
            return float(x)
        except Exception:
            return None
    agg = {}                                            # (chrom,pos,ref,alt) -> record + sample set
    for uid, ods in by_uid.items():
        for od in ods:
            for sp in subpaths:
                p = os.path.join(od, sp)
                if not os.path.isfile(p):
                    continue
                with open(p) as fh:
                    for row in _csv.DictReader(fh, delimiter="\t"):
                        cons = row.get("CONSEQUENCE", "") or ""
                        if not any(f in cons for f in _FUNCTIONAL):
                            continue
                        gene = row.get("GENE", "") or ""
                        if any(gene.startswith(pfx) for pfx in _ARTIFACT_GENE_PREFIX):
                            continue                    # drop artifact-prone gene classes
                        clinvar_path = "pathogenic" in (row.get("CLINVAR", "") or "").lower()
                        in_prio = bool(priority_genes) and (gene in priority_genes)
                        in_wl = (gene_whitelist is None) or (gene in gene_whitelist)
                        if not (in_wl or clinvar_path):
                            continue                    # KEEP: cancer gene (whitelist) OR ClinVar-pathogenic
                        # source of support (for annotation + ranking): AML/MDS panel > cancer gene > ClinVar-only
                        gene_source = "aml_mds" if in_prio else ("cancer_gene" if in_wl else "clinvar_pathogenic")
                        sc = (row.get("SOMATIC_CLASS", "") or "").upper()
                        if "GERMLINE" in sc:
                            continue
                        g = _num(row.get("gnomAD_AF_max"))
                        if g is not None and g > max_gnomad:
                            continue
                        if require_cosmic_clinvar and not (row.get("COSMIC") or row.get("CLINVAR")):
                            continue                    # cancer-database (COSMIC/ClinVar) gating
                        k = (row.get("CHROM"), row.get("POS"), row.get("REF"), row.get("ALT"))
                        if None in k:
                            continue
                        rec = agg.get(k)
                        if rec is None:
                            rec = agg[k] = dict(gene=row.get("GENE", ""), consequence=cons,
                                                hgvsp=row.get("HGVSp", ""), cosmic=row.get("COSMIC", ""),
                                                clinvar=row.get("CLINVAR", ""), gnomad=row.get("gnomAD_AF_max", ""),
                                                source=gene_source, samples=set())
                        else:
                            # keep the strongest source seen for this variant across tiers/samples
                            _rank = {"aml_mds": 0, "cancer_gene": 1, "clinvar_pathogenic": 2}
                            if _rank[gene_source] < _rank[rec["source"]]:
                                rec["source"] = gene_source
                        rec["samples"].add(uid)
    _srank = {"aml_mds": 0, "cancer_gene": 1, "clinvar_pathogenic": 2}
    rows = []
    for (ch, pos, ref, alt), r in agg.items():
        rows.append((ch, pos, ref, alt, r["gene"], r["consequence"], r["hgvsp"], r["cosmic"],
                     r["clinvar"], r["gnomad"], r["source"], len(r["samples"]),
                     ",".join(sorted(r["samples"]))))
    # rank: AML/MDS panel first, then cancer gene, then ClinVar-only; within each, recurrence desc
    rows.sort(key=lambda x: (_srank.get(x[10], 9), -x[11], 0 if (x[7] or x[8]) else 1))
    with open(out_path, "w", newline="") as fh:
        w = _csv.writer(fh, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt", "gene", "consequence", "HGVSp", "cosmic", "clinvar",
                    "gnomAD_AF_max", "gene_source", "n_samples", "samples"])
        for row in rows:
            w.writerow(row)
    logger.info("nominate_functional_drivers: %d candidate driver mutations -> %s", len(rows), out_path)
    return out_path


def assess_matrix_quantification(matrix_dir, out_path):
    """Cohort performance: per sample, of the four-caller-union passenger variants genotyped into the
    matrix, how many were QUANTIFIED at the cell-barcode level (>=1 cell with read coverage) vs all-zero
    (no single-cell reads). Reads the built matrices only (no BAM). Writes a per-sample + cohort report."""
    import csv
    import glob as _glob
    import anndata as _ad
    import scipy.sparse as _sp
    rows = []
    tot_q = tot_p = 0
    for h in sorted(_glob.glob(os.path.join(matrix_dir, "*_cell_feature_matrix.h5ad"))):
        uid = os.path.basename(h).replace("_cell_feature_matrix.h5ad", "")
        a = _ad.read_h5ad(h)
        pmask = (a.var["feature_type"].values == "mtDNA")
        total = a.layers["total"]
        total = total.toarray() if _sp.issparse(total) else np.asarray(total)
        covered = (total[:, pmask].sum(axis=0) > 0)
        npass = int(pmask.sum()); nq = int(covered.sum())
        rows.append((uid, npass, nq, round(nq / npass, 4) if npass else ""))
        tot_p += npass; tot_q += nq
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "four_caller_passengers", "quantified_at_cell_level", "frac_quantified"])
        for r in rows:
            w.writerow(r)
        w.writerow(["ALL", tot_p, tot_q, round(tot_q / tot_p, 4) if tot_p else ""])
    logger.info("matrix quantification: %d/%d passengers quantified at cell level (%.4f) -> %s",
                tot_q, tot_p, (tot_q / tot_p if tot_p else 0.0), out_path)
    return out_path


def export_patient_variants(matrix_dir, patient_map_path, hits_table, out_h5ad_dir, out_tsv_dir):
    """Per PATIENT (all its samples combined): (1) write a patient h5ad of driver + passenger genotypes
    -- obs index is 'uid|barcode' so the SAMPLE NAME is carried in every cell id, var is the UNION of
    features across the patient's samples -- to out_h5ad_dir; (2) for each PATIENT-level significant
    (driver, best-passenger) hit, write a TSV of the cells that CARRY that passenger (barcode, sample) --
    the per-cell prediction of the driver mutation -- to out_tsv_dir. e.g. for RUNX1_chr21:34799432 in
    patient 5801, every cell carrying LINC02256_chr15:32580386 is written as a RUNX1 prediction."""
    import csv
    import anndata as _ad
    import scipy.sparse as _sp
    os.makedirs(out_h5ad_dir, exist_ok=True)
    os.makedirs(out_tsv_dir, exist_ok=True)
    pmap = {}
    with open(patient_map_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) >= 2:
                pmap.setdefault(p[1], []).append(p[0])
    hits = []
    if hits_table and os.path.isfile(hits_table):
        with open(hits_table) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("level") == "patient":
                    hits.append(r)
    for patient, uids in pmap.items():
        mats = [_ad.read_h5ad(os.path.join(matrix_dir, f"{u}_cell_feature_matrix.h5ad"))
                for u in uids if os.path.isfile(os.path.join(matrix_dir, f"{u}_cell_feature_matrix.h5ad"))]
        if not mats:
            continue
        comb = mats[0] if len(mats) == 1 else _ad.concat(mats, join="outer", merge="first")
        combp = os.path.join(out_h5ad_dir, f"{patient}_variants.h5ad")
        comb.write_h5ad(combp, compression="gzip")
        logger.info("[%s] patient h5ad: %d cells x %d features -> %s", patient, comb.n_obs, comb.n_vars, combp)
        alt = comb.layers["alt"]
        alt = alt.toarray() if _sp.issparse(alt) else np.asarray(alt)
        vnames = list(comb.var_names)
        for r in [h for h in hits if h.get("unit") == patient]:
            driver = r.get("driver", "")
            marker_full = r.get("marker", "")
            marker = marker_full.rsplit("_", 1)[0] if (">" in marker_full) else marker_full
            if marker not in vnames:
                continue
            j = vnames.index(marker)
            tp = os.path.join(out_tsv_dir, f"{patient}__{_safe(driver)}__pred_by__{_safe(marker)}.tsv")
            with open(tp, "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["barcode", "sample", "predicted_driver", "passenger_marker", "sens", "spec", "p"])
                n = 0
                for i in range(comb.n_obs):
                    if alt[i, j] >= 1:
                        bc = str(comb.obs_names[i])
                        w.writerow([bc, bc.split("|")[0], driver, marker_full,
                                    r.get("sens", ""), r.get("spec", ""), r.get("p", "")])
                        n += 1
            logger.info("[%s] %s predicted by %s: %d cells -> %s", patient, driver, marker, n, tp)


def run_variant_impact(args):
    metadata = pd.read_csv(args.metadata, sep="\t")
    out_root = os.path.abspath(getattr(args, "output", None) or os.path.join(os.getcwd(), "variant_impact"))
    os.makedirs(out_root, exist_ok=True)
    if getattr(args, "passenger_association", False):
        run_supervised_passenger_association(args)
        return
    if getattr(args, "nominate_functional_drivers", False):
        wl = None
        gw = getattr(args, "gene_whitelist", None)
        if gw and os.path.isfile(gw):
            wl = {ln.strip() for ln in open(gw) if ln.strip() and not ln.startswith("#")}
        prio = None
        pg = getattr(args, "priority_genes", None)
        if pg and os.path.isfile(pg):
            prio = {ln.strip() for ln in open(pg) if ln.strip() and not ln.startswith("#")}
        nominate_functional_drivers(args.discovery_manifest,
                                    os.path.join(out_root, "candidate_drivers.tsv"),
                                    gene_whitelist=wl, priority_genes=prio,
                                    require_cosmic_clinvar=getattr(args, "require_cosmic_clinvar", False))
        return
    if getattr(args, "assess_quantification", False):
        mdir = getattr(args, "matrix_dir", None) or os.path.join(out_root, "matrices")
        assess_matrix_quantification(mdir, os.path.join(out_root, "matrix_quantification.tsv"))
        return
    if getattr(args, "export_patient", False):
        mdir = getattr(args, "matrix_dir", None) or os.path.join(out_root, "matrices")
        export_patient_variants(mdir, args.patient_map, getattr(args, "hits_table", None),
                                getattr(args, "h5ad_dir", None) or os.path.join(out_root, "patient_h5ad"),
                                getattr(args, "tsv_dir", None) or os.path.join(out_root, "patient_predictions"))
        return
    # MAPQ: supervised paths default to the callers' threshold (1, matching kinnex_variant_pipeline.sh);
    # the discovery genotyper keeps 20. An explicit --min-mapq overrides either.
    _sup = getattr(args, "supervised", False) or getattr(args, "validate_supervised", False)
    _mapq = getattr(args, "min_mapq", None)
    if _mapq is None:
        _mapq = 1 if _sup else 20
    if getattr(args, "validate_supervised", False):
        # Supervised-genotyper validation: for the requested sample(s), genotype EVERY four-caller
        # variant of THAT sample (discovery outdir = <bam_dir>/genome_variant_detection) and report the
        # specific variants the supervised genotyper fails to detect. No impact/MarkerFinder pipeline.
        only = set(getattr(args, "samples", None) or [])
        for _, r in metadata.iterrows():
            uid = str(r.get("uid", "")).strip()
            library = str(r.get("library", uid)).strip()
            bam = str(r.get("bam", ""))
            if not uid or (only and uid not in only and library not in only):
                continue
            if not os.path.exists(bam):
                logger.warning("[%s] bam not found: %s -- skipped", uid, bam); continue
            disc = os.path.join(os.path.dirname(bam), "genome_variant_detection")
            if not os.path.isdir(disc):
                logger.warning("[%s] discovery outdir not found: %s -- skipped", uid, disc); continue
            validate_supervised(bam, disc, os.path.join(out_root, _safe(uid), "supervised_validation"),
                                min_mapq=_mapq, indel_window=getattr(args, "indel_window", 5))
        logger.info("supervised validation done -> %s", out_root)
        return
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
        for (chrom, pos, label, vtype, ref, alt) in variants:
            if getattr(args, "supervised", False):
                if not (ref and alt):
                    raise ValueError("--supervised requires REF and ALT per variant "
                                     f"('chr:pos:label:type:REF:ALT' or file cols 6,7); missing for "
                                     f"{chrom}:{pos} ({label})")
                counts = genotype_from_bam_supervised(bam, chrom, pos, ref, alt,
                                                      min_mapq=_mapq,
                                                      indel_window=getattr(args, "indel_window", 5))
            else:
                counts = genotype_from_bam(bam, chrom, pos, vtype=vtype,
                                           min_mapq=_mapq,
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


# ============================================================================
# PER-CELL CLASSIFIER BENCHMARK  (method assessment / improvement)
# ----------------------------------------------------------------------------
# Question 1 (viability / "where does the sensitivity lie"): for a given
#   mutation x sample x cell-state, is there a REAL per-cell expression signal
#   separating read-genotyped MUT from WT cells? Answered by held-out AUROC vs a
#   LABEL-PERMUTATION null, which is valid at small n.
# Question 2 (can we beat the deployed classifier): head-to-head, on IDENTICAL
#   cross-validation folds, of the deployed median-centered nearest-centroid
#   ("centroid") against L2 logistic regression ("logreg") and shrinkage LDA
#   ("lda"). A challenger is adopted only if it beats the centroid by a paired
#   test across the real-signal units; otherwise the centroid is kept.
# Everything is PER CELL. Feature (marker) selection is redone on the TRAINING
# cells of every fold, so no test cell ever informs its own signature.
# ============================================================================

BENCH_METHODS = ("centroid", "logreg", "lda")


def _bench_folds(mut, wt, k=5, repeats=5, seed=0):
    """Repeated stratified k-fold: list of (rep, train_mut, train_wt, test_mut, test_wt).
    In each repeat every cell is held out exactly once, so a per-repeat AUROC uses
    independent (once-held-out) predictions."""
    rng = np.random.RandomState(seed)
    mut, wt = list(mut), list(wt)
    out = []
    for rep in range(repeats):
        mf = np.array_split(rng.permutation(mut), k)
        wf = np.array_split(rng.permutation(wt), k)
        for fi in range(k):
            te_m, te_w = list(mf[fi]), list(wf[fi])
            tr_m = [c for j in range(k) if j != fi for c in mf[j]]
            tr_w = [c for j in range(k) if j != fi for c in wf[j]]
            if min(len(tr_m), len(tr_w)) < 2 or not te_m or not te_w:
                continue
            out.append((rep, tr_m, tr_w, te_m, te_w))
    return out


def _centroid_scores(Xall, tr_m, tr_w, te):
    """Deployed classifier's per-cell score on a precomputed signature matrix Xall
    (rows = cells, cols = signature features). Identical rule to `impute_unknown`:
    median-of-centroids reference, then nearest-centroid margin (dW - dM); >0 -> MUT."""
    mut_c = np.median(Xall.loc[list(tr_m)].to_numpy(dtype=float), axis=0)
    wt_c = np.median(Xall.loc[list(tr_w)].to_numpy(dtype=float), axis=0)
    ref = np.median(np.vstack([mut_c, wt_c]), axis=0)
    V = Xall.loc[list(te)].to_numpy(dtype=float) - ref
    dM = np.linalg.norm(V - (mut_c - ref), axis=1)
    dW = np.linalg.norm(V - (wt_c - ref), axis=1)
    return dict(zip([str(b) for b in te], (dW - dM).astype(float)))


def _sk_scores(Xall, tr_m, tr_w, te, method):
    """Standardized linear classifier (logreg / lda) on the same signature matrix;
    returns the signed decision function per test cell (higher -> MUT)."""
    from sklearn.preprocessing import StandardScaler
    ytr = np.array([1] * len(tr_m) + [0] * len(tr_w))
    Xtr = Xall.loc[list(tr_m) + list(tr_w)].to_numpy(dtype=float)
    Xte = Xall.loc[list(te)].to_numpy(dtype=float)
    if len(set(ytr)) < 2 or Xtr.shape[0] < 4:
        return None
    sc = StandardScaler().fit(Xtr)
    Atr, Ate = sc.transform(Xtr), sc.transform(Xte)
    try:
        if method == "logreg":
            from sklearn.linear_model import LogisticRegression
            clf = LogisticRegression(penalty="l2", C=1.0, max_iter=2000, class_weight="balanced")
        elif method == "lda":
            from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
            clf = LinearDiscriminantAnalysis(solver="lsqr", shrinkage="auto")
        else:
            return None
        clf.fit(Atr, ytr)
        p = clf.decision_function(Ate)
    except Exception as e:
        logger.warning("bench sklearn %s failed: %s", method, e)
        return None
    return dict(zip([str(b) for b in te], np.asarray(p, dtype=float)))


def _en_full_scores(norm, tr_m, tr_w, te, n_hvg=2000):
    """Genome-wide elastic-net logistic regression: NO MarkerFinder pre-selection --
    the top n_hvg highest-variance genes (chosen on TRAIN only) fed to an L1/L2
    (elastic-net) logistic model with embedded selection. Tests whether the mutation's
    program lives in genes the top-50 correlation markers miss. Signed decision fn per
    test cell (>0 -> MUT)."""
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    cells = list(tr_m) + list(tr_w) + list(te)
    sub = norm[cells]
    X = sub.X.toarray() if issparse(sub.X) else np.asarray(sub.X, dtype=float)
    ntr = len(tr_m) + len(tr_w)
    Xtr_all, Xte_all = X[:ntr], X[ntr:]
    if Xtr_all.shape[0] < 4 or Xte_all.shape[0] == 0:
        return None
    hv = np.argsort(Xtr_all.var(axis=0))[::-1][:n_hvg]
    Xtr, Xte = Xtr_all[:, hv], Xte_all[:, hv]
    y = np.array([1] * len(tr_m) + [0] * len(tr_w))
    try:
        sc = StandardScaler().fit(Xtr)
        clf = LogisticRegression(penalty="elasticnet", l1_ratio=0.5, C=0.5, solver="saga",
                                 max_iter=3000, class_weight="balanced")
        clf.fit(sc.transform(Xtr), y)
        p = clf.decision_function(sc.transform(Xte))
    except Exception as e:
        logger.warning("en_full failed: %s", e)
        return None
    return dict(zip([str(b) for b in te], np.asarray(p, dtype=float)))


def _fold_scores(base, norm, tr_m, tr_w, te, top_n, methods, mode="pearson"):
    """Select the MUT/WT signature on the TRAIN cells once (marker mode = pearson or
    wilcoxon), pull the dense signature matrix for train+test once, then score every
    method on it. 'en_full' bypasses marker selection (genome-wide elastic-net)."""
    if list(methods) == ["en_full"]:                    # marker-free; skip MarkerFinder
        s = _en_full_scores(norm, tr_m, tr_w, te)
        return {"en_full": s} if s else {}
    sig = markers_from_cells(base, list(tr_m) + list(tr_w),
                             ["MUT"] * len(tr_m) + ["WT"] * len(tr_w), top_n, mode=mode)
    if not sig:
        return {}
    Xall = _expr_subset(norm, list(tr_m) + list(tr_w) + list(te), sig)
    if Xall.empty or not set(map(str, te)).issubset(set(map(str, Xall.index))):
        # keep only test cells actually present
        pass
    have = set(map(str, Xall.index))
    te2 = [c for c in te if str(c) in have]
    trm2 = [c for c in tr_m if str(c) in have]
    trw2 = [c for c in tr_w if str(c) in have]
    if min(len(trm2), len(trw2)) < 2 or not te2:
        return {}
    out = {}
    for m in methods:
        if m == "centroid":
            out[m] = _centroid_scores(Xall, trm2, trw2, te2)
        elif m == "en_full":
            out[m] = _en_full_scores(norm, trm2, trw2, te2)
        else:
            out[m] = _sk_scores(Xall, trm2, trw2, te2, m)
    return out


def benchmark_unit(base, norm, mut_cells, wt_cells, top_n=50, k=5, repeats=5,
                   methods=BENCH_METHODS, n_perm=0, seed=0, label="", marker_mode="pearson"):
    """Held-out AUROC of each method on ONE unit (read-genotyped MUT/WT cells).
    Per-repeat AUROC (once-held-out predictions) -> mean +/- SD across repeats.
    n_perm>0 adds a label-permutation null for the reference 'centroid' method ->
    two-sided p for 'is the per-cell signal real'. `marker_mode` selects the
    signature method (pearson vs wilcoxon). All predictions are per cell."""
    mut_cells = [str(c) for c in mut_cells]
    wt_cells = [str(c) for c in wt_cells]
    folds = _bench_folds(mut_cells, wt_cells, k, repeats, seed)
    if not folds:
        return None
    truth = {c: 1 for c in mut_cells}
    truth.update({c: 0 for c in wt_cells})
    per = {m: {} for m in methods}                       # method -> rep -> {cell: score}
    for (rep, tr_m, tr_w, te_m, te_w) in folds:
        fs = _fold_scores(base, norm, tr_m, tr_w, list(te_m) + list(te_w), top_n, methods, mode=marker_mode)
        for m in methods:
            if fs.get(m):
                per[m].setdefault(rep, {}).update(fs[m])

    def repeat_aurocs(scmap):
        vals = []
        for _rep, sc in scmap.items():
            cells = [c for c in sc if c in truth]
            if len({truth[c] for c in cells}) == 2:
                a = auroc([sc[c] for c in cells], [truth[c] == 1 for c in cells])
                if np.isfinite(a):
                    vals.append(a)
        return vals

    def confusion(scmap):
        # pooled over all held-out repeats; predict MUT iff score>0 (deployed sign rule).
        tp = fp = tn = fn = 0
        for _rep, sc in scmap.items():
            for c, s in sc.items():
                if c not in truth or not np.isfinite(s):
                    continue
                pred_mut = s > 0
                if truth[c] == 1:
                    tp += int(pred_mut); fn += int(not pred_mut)
                else:
                    fp += int(pred_mut); tn += int(not pred_mut)
        sens = tp / (tp + fn) if (tp + fn) else np.nan
        spec = tn / (tn + fp) if (tn + fp) else np.nan
        return tp, fp, tn, fn, sens, spec

    out = {"label": label, "marker_mode": marker_mode, "n_MUT": len(mut_cells), "n_WT": len(wt_cells)}
    for m in methods:
        av = repeat_aurocs(per[m])
        out[f"{m}_auroc"] = round(float(np.mean(av)), 4) if av else np.nan
        out[f"{m}_auroc_sd"] = round(float(np.std(av)), 4) if av else np.nan
        out[f"{m}_n_rep"] = len(av)
        tp, fp, tn, fn, sens, spec = confusion(per[m])
        out[f"{m}_TP"], out[f"{m}_FP"], out[f"{m}_TN"], out[f"{m}_FN"] = tp, fp, tn, fn
        out[f"{m}_sens"] = round(float(sens), 4) if np.isfinite(sens) else np.nan
        out[f"{m}_spec"] = round(float(spec), 4) if np.isfinite(spec) else np.nan
    if n_perm and "centroid" in methods:
        obs = out.get("centroid_auroc", np.nan)
        rng = np.random.RandomState(seed + 7)
        allc = mut_cells + wt_cells
        nm = len(mut_cells)
        null = []
        for _ in range(int(n_perm)):
            perm = list(rng.permutation(allc))
            pm, pw = perm[:nm], perm[nm:]
            ptruth = {c: 1 for c in pm}
            ptruth.update({c: 0 for c in pw})
            sc_all = {}
            for (_r, a, b, cc, dd) in _bench_folds(pm, pw, k, 1, seed=int(rng.randint(1, 10 ** 6))):
                fs = _fold_scores(base, norm, a, b, list(cc) + list(dd), top_n, ("centroid",), mode=marker_mode)
                if fs.get("centroid"):
                    sc_all.update(fs["centroid"])
            cells = [c for c in sc_all if c in ptruth]
            if len({ptruth[c] for c in cells}) == 2:
                na = auroc([sc_all[c] for c in cells], [ptruth[c] == 1 for c in cells])
                if np.isfinite(na):
                    null.append(na)
        if null and np.isfinite(obs):
            null = np.asarray(null)
            out["perm_p_centroid"] = round(float((np.sum(np.abs(null - 0.5) >= abs(obs - 0.5)) + 1)
                                                  / (len(null) + 1)), 4)
        else:
            out["perm_p_centroid"] = np.nan
        out["perm_n"] = len(null)
    return out


def benchmark_method_comparison(df, methods=BENCH_METHODS, reference="centroid"):
    """Paired comparison across benchmarked units: mean/median held-out AUROC per
    method and a Wilcoxon signed-rank test of each challenger vs the reference
    (`centroid`, the deployed classifier), paired by unit within each marker_mode.
    Answers 'does any classifier systematically beat the current one'."""
    from scipy import stats
    rows = []
    for mode, sub in df.groupby("marker_mode"):
        for m in methods:
            paired = sub.dropna(subset=[f"{m}_auroc", f"{reference}_auroc"])
            a = paired[f"{m}_auroc"].to_numpy(dtype=float)
            b = paired[f"{reference}_auroc"].to_numpy(dtype=float)
            p = np.nan
            if m != reference and len(paired) >= 5 and np.any(a != b):
                try:
                    p = float(stats.wilcoxon(a, b, zero_method="wilcox").pvalue)
                except Exception:
                    p = np.nan
            rows.append({"marker_mode": mode, "method": m, "n_units": int(len(paired)),
                         "mean_auroc": round(float(np.mean(a)), 4) if len(a) else np.nan,
                         "median_auroc": round(float(np.median(a)), 4) if len(a) else np.nan,
                         "wins_vs_ref": int((a > b).sum()), "losses_vs_ref": int((a < b).sum()),
                         "wilcoxon_p_vs_ref": round(p, 4) if np.isfinite(p) else np.nan})
    return pd.DataFrame(rows)


def benchmark_units(base, norm, units, top_n=50, k=5, repeats=5, methods=BENCH_METHODS,
                    n_perm=0, min_cells=15, seed=0, marker_modes=("pearson",)):
    """Run `benchmark_unit` over many units and marker modes. `units` = list of
    (label, mut_cells, wt_cells). Skips units with < min_cells in either class.
    Returns a long DataFrame, one row per (unit, marker_mode)."""
    rows = []
    for mode in marker_modes:
        for (label, mut_cells, wt_cells) in units:
            if len(mut_cells) < min_cells or len(wt_cells) < min_cells:
                continue
            r = benchmark_unit(base, norm, mut_cells, wt_cells, top_n=top_n, k=k, repeats=repeats,
                               methods=methods, n_perm=n_perm, seed=seed, label=label, marker_mode=mode)
            if r:
                rows.append(r)
                logger.info("bench [%s] %s: " % (mode, label)
                            + " ".join(f"{m}={r.get(m + '_auroc')}" for m in methods))
    return pd.DataFrame(rows)


# ============================================================================
# MARKER RECURRENCE ACROSS WILD-TYPE REFERENCES  (marker-stability test)
# ----------------------------------------------------------------------------
# Hold the MUTANT set fixed and swap the WILD-TYPE reference across several cell
# states; extract MarkerFinder markers for each reference; ask how much the marker
# genes recur across references and whether that recurrence exceeds chance. If the
# mutant expression program is real, the MUT-direction markers recur across every
# WT reference; if it is noise, the marker sets differ each time (overlap ~ chance).
# ============================================================================

def markers_with_direction(base, cells, labels, top_n, mode="pearson", rho_threshold=0.0):
    """Marker features WITH the group each marks: DataFrame[feature, group] where group
    is 'MUT' or 'WT' (the class the gene is UP in). pearson = cellHarmony MarkerFinder;
    wilcoxon = scanpy rank_genes_groups."""
    sub = base[list(cells)].copy()
    sub.obs["genotype"] = pd.Series(labels, index=list(cells)).loc[sub.obs_names].astype(str).values
    normalize_adata(sub)
    if mode in ("wilcoxon", "scanpy"):
        try:
            import scanpy as sc
            sc.tl.rank_genes_groups(sub, "genotype", method="wilcoxon", n_genes=top_n)
            names = sub.uns["rank_genes_groups"]["names"]
            rows = [(str(g), grp) for grp in names.dtype.names for g in list(names[grp])[:top_n]]
            return pd.DataFrame(rows, columns=["feature", "group"]).drop_duplicates("feature")
        except Exception as e:
            logger.warning("wilcoxon markers_with_direction failed: %s", e)
            return pd.DataFrame(columns=["feature", "group"])
    from altanalyze3.components.cellHarmony.markerFinder import find_markers_from_adata
    try:
        mo = find_markers_from_adata(sub, "genotype", n_markers=top_n, direction="up",
                                     rho_threshold=rho_threshold, write_outputs=False)
    except Exception as e:
        logger.warning("pearson markers_with_direction failed: %s", e)
        return pd.DataFrame(columns=["feature", "group"])
    mk = mo.markers
    if mk is None or not len(mk):
        return pd.DataFrame(columns=["feature", "group"])
    return pd.DataFrame({"feature": mk["marker"].astype(str),
                         "group": mk["cluster"].astype(str)}).drop_duplicates("feature")


def _expressed_universe(base, cells):
    """Genes with any expression across `cells` (the pool MarkerFinder can draw from)."""
    sub = base[list(dict.fromkeys(cells))]
    X = sub.X
    nz = np.asarray((X != 0).sum(axis=0)).ravel() if issparse(X) else (np.asarray(X) != 0).sum(axis=0)
    return [str(g) for g, e in zip(sub.var_names, nz) if e > 0]


def marker_recurrence(base, mut_cells, wt_references, top_n=50, mode="pearson",
                      n_random=2000, seed=0):
    """Fixed MUT set vs several WT references (`wt_references` = {name: [cells]}).
    Extract markers per reference, then measure recurrence of the marker genes across
    references and whether it exceeds a random-gene-set null. Returns the per-reference
    marker sets, pairwise/n-way overlaps, expected overlap and permutation p -- computed
    both for ALL markers and for MUT-direction markers only (the mutant program)."""
    import itertools
    universe = _expressed_universe(base, list(mut_cells) + [c for ws in wt_references.values() for c in ws])
    U = len(universe)
    sets_all, sets_mut, sizes = {}, {}, {}
    for name, wt in wt_references.items():
        md = markers_with_direction(base, list(mut_cells) + list(wt),
                                    ["MUT"] * len(mut_cells) + ["WT"] * len(wt), top_n, mode)
        sets_all[name] = set(md["feature"])
        sets_mut[name] = set(md.loc[md["group"] == "MUT", "feature"])
        sizes[name] = {"n_wt": len(wt), "n_markers": len(sets_all[name]), "n_mut_markers": len(sets_mut[name])}

    def stats(sets):
        names = list(sets)
        vals = [sets[n] for n in names]
        obs = len(set.intersection(*vals)) if all(len(v) for v in vals) else 0
        pair = {f"{a}&{b}": len(sets[a] & sets[b]) for a, b in itertools.combinations(names, 2)}
        if not all(len(v) for v in vals) or U == 0:
            return {"n_way_overlap": obs, "expected": np.nan, "p_random": np.nan, "pairwise": pair}
        rng = np.random.RandomState(seed)
        uni = np.array(universe)
        ss = [len(v) for v in vals]
        null = np.array([len(set.intersection(*[set(rng.choice(uni, s, replace=False)) for s in ss]))
                         for _ in range(n_random)])
        return {"n_way_overlap": obs, "expected": round(float(null.mean()), 3),
                "p_random": round(float((null >= obs).mean()), 4), "pairwise": pair}

    return {"universe": U, "top_n": top_n, "mode": mode, "sizes": sizes,
            "sets_all": sets_all, "sets_mut": sets_mut,
            "all_markers": stats(sets_all), "mut_markers": stats(sets_mut)}


def marker_overlap_random_mut_null(base, pool_cells, wt_references, n_mut, top_n=50,
                                   mode="pearson", n_perm=200, seed=0,
                                   pooled_up=None, pooled_down=None):
    """Negative control for `marker_recurrence`: is the marker recurrence driven by the
    MUTATION, or would ANY equally-sized cell group replicate it against the same
    wild-type references? Draw random `n_mut`-cell pseudo-mutant sets from `pool_cells`
    (excluding the wild-type reference cells) and recompute markers vs each reference.
    Tracks, as null distributions:
      - pairwise & n-way overlaps among references, for UP-in-mutant and DOWN-in-mutant;
      - each reference's overlap with a fixed pooled UP / pooled DOWN marker set
        (`pooled_up`/`pooled_down`), so the observed per-state-vs-pooled overlap can be
        judged mutation-specific vs cell-identity."""
    import itertools
    rng = np.random.RandomState(seed)
    wt_all = set(c for ws in wt_references.values() for c in ws)
    pool = [c for c in dict.fromkeys(pool_cells) if c not in wt_all]
    names = list(wt_references)
    pair_up = {f"{a}&{b}": [] for a, b in itertools.combinations(names, 2)}
    pair_down = {f"{a}&{b}": [] for a, b in itertools.combinations(names, 2)}
    nway_up, nway_down = [], []
    up_vs_pool = {n: [] for n in names}
    down_vs_pool = {n: [] for n in names}
    for _ in range(int(n_perm)):
        pmut = list(rng.choice(pool, size=min(n_mut, len(pool)), replace=False))
        up, down = {}, {}
        for name, wt in wt_references.items():
            md = markers_with_direction(base, pmut + list(wt),
                                        ["MUT"] * len(pmut) + ["WT"] * len(wt), top_n, mode)
            up[name] = set(md.loc[md["group"] == "MUT", "feature"])
            down[name] = set(md.loc[md["group"] == "WT", "feature"])
            if pooled_up is not None:
                up_vs_pool[name].append(len(up[name] & pooled_up))
            if pooled_down is not None:
                down_vs_pool[name].append(len(down[name] & pooled_down))
        for a, b in itertools.combinations(names, 2):
            pair_up[f"{a}&{b}"].append(len(up[a] & up[b]))
            pair_down[f"{a}&{b}"].append(len(down[a] & down[b]))
        nway_up.append(len(set.intersection(*[up[n] for n in names])) if all(len(up[n]) for n in names) else 0)
        nway_down.append(len(set.intersection(*[down[n] for n in names])) if all(len(down[n]) for n in names) else 0)
    arr = lambda d: {k: np.array(v) for k, v in d.items()}
    return {"pairwise_up": arr(pair_up), "pairwise_down": arr(pair_down),
            "nway_up": np.array(nway_up), "nway_down": np.array(nway_down),
            "up_vs_pooled": arr(up_vs_pool), "down_vs_pooled": arr(down_vs_pool),
            # back-compat aliases (UP direction was the original output)
            "pairwise_null": arr(pair_up), "nway_null": np.array(nway_up),
            "n_perm": int(n_perm), "n_mut": int(n_mut)}


def plot_venn3(sets, out_pdf, title="", subtitle=""):
    """3-set Venn (rendered manually; editable PDF, Arial, explicit RGB). `sets` =
    ordered {name: set} with exactly three entries. Region counts are exact."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    names = list(sets)
    A, B, C = (sets[names[0]], sets[names[1]], sets[names[2]])
    vals = {"A": len(A - B - C), "B": len(B - A - C), "C": len(C - A - B),
            "AB": len((A & B) - C), "AC": len((A & C) - B), "BC": len((B & C) - A),
            "ABC": len(A & B & C)}
    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    cen = {"A": (-0.33, 0.22), "B": (0.33, 0.22), "C": (0.0, -0.38)}
    col = {"A": (0.85, 0.37, 0.34), "B": (0.30, 0.55, 0.80), "C": (0.55, 0.75, 0.45)}
    for k in "ABC":
        ax.add_patch(Circle(cen[k], 0.62, facecolor=col[k], edgecolor=(0.2, 0.2, 0.2), lw=1.2, alpha=0.45))
    pos = {"A": (-0.62, 0.5), "B": (0.62, 0.5), "C": (0.0, -0.85), "AB": (0.0, 0.52),
           "AC": (-0.42, -0.28), "BC": (0.42, -0.28), "ABC": (0.0, 0.0)}
    for k, (x, y) in pos.items():
        ax.text(x, y, str(vals[k]), ha="center", va="center", fontsize=12, fontweight="bold")
    for k, lbl in zip(("A", "B", "C"), ((-0.7, 1.0), (0.7, 1.0), (0.0, -1.18))):
        ax.text(lbl[0], lbl[1], f"{names['ABC'.index(k)]}\n(n={len(sets[names['ABC'.index(k)]])})",
                ha="center", va="center", fontsize=8)
    ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.45, 1.35); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title((title + ("\n" + subtitle if subtitle else "")), fontsize=9)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    return out_pdf


def main(argv=None):
    """Direct module entry: python -m altanalyze3.components.bam.variant_impact ...  Provided because a
    partial altanalyze3 checkout may not import the top-level CLI (parser.py); these flags mirror the
    `variant-impact` subcommand for the supervised + passenger-association paths."""
    import argparse
    ap = argparse.ArgumentParser("variant_impact")
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--output", default=None)
    ap.add_argument("--samples", nargs="*", default=None)
    ap.add_argument("--variants", default=None)
    ap.add_argument("--level", default="gene")
    ap.add_argument("--cell_annot", default=None)
    ap.add_argument("--gene_symbol", default=None)
    ap.add_argument("--supervised", action="store_true")
    ap.add_argument("--validate-supervised", dest="validate_supervised", action="store_true")
    ap.add_argument("--passenger-association", dest="passenger_association", action="store_true")
    ap.add_argument("--build-only", dest="build_only", action="store_true")
    ap.add_argument("--nominate-only", dest="nominate_only", action="store_true")
    ap.add_argument("--driver-variants", dest="driver_variants", default=None)
    ap.add_argument("--results-extraction", dest="results_extraction", default=None)
    ap.add_argument("--patient-map", dest="patient_map", default=None)
    ap.add_argument("--min-mapq", dest="min_mapq", default=None, type=int)
    ap.add_argument("--indel-window", dest="indel_window", default=5, type=int)
    ap.add_argument("--max-q", dest="max_q", default=0.05, type=float)
    ap.add_argument("--max-background", dest="max_background", default=0.15, type=float)
    ap.add_argument("--nominate-min-mut", dest="nominate_min_mut", default=10, type=int)
    ap.add_argument("--limit-markers", dest="limit_markers", default=None, type=int, help="test only: genotype the first N passengers per sample")
    ap.add_argument("--export-patient", dest="export_patient", action="store_true", help="Per patient: combined driver+passenger h5ad + per-prediction cell-barcode TSVs (needs --matrix-dir, --patient-map, --hits-table, --h5ad-dir, --tsv-dir).")
    ap.add_argument("--matrix-dir", dest="matrix_dir", default=None)
    ap.add_argument("--hits-table", dest="hits_table", default=None)
    ap.add_argument("--h5ad-dir", dest="h5ad_dir", default=None)
    ap.add_argument("--tsv-dir", dest="tsv_dir", default=None)
    ap.add_argument("--nominate-functional-drivers", dest="nominate_functional_drivers", action="store_true", help="Nominate candidate driver mutations from all 4-caller variants by functional impact + recurrence (needs --discovery-manifest).")
    ap.add_argument("--assess-quantification", dest="assess_quantification", action="store_true", help="Cohort report: per sample, how many 4-caller passengers were quantified at the cell-barcode level (needs --matrix-dir).")
    ap.add_argument("--discovery-manifest", dest="discovery_manifest", default=None)
    ap.add_argument("--gene-whitelist", dest="gene_whitelist", default=None, help="one gene per line; keep functional-driver candidates whose GENE is in this cancer-gene list (OncoKB + AML/MDS panel + MEIS1). ClinVar-pathogenic variants are kept regardless.")
    ap.add_argument("--priority-genes", dest="priority_genes", default=None, help="one gene per line; the AML/MDS panel subset. Candidates in these genes get gene_source=aml_mds and sort to the top.")
    ap.add_argument("--require-cosmic-clinvar", dest="require_cosmic_clinvar", action="store_true", help="functional-driver nomination: keep only variants with COSMIC or ClinVar support.")
    ap.add_argument("--mut-min", dest="mut_min", default=1, type=int)
    ap.add_argument("--wt-min", dest="wt_min", default=2, type=int)
    ap.add_argument("--min-cells", dest="min_cells", default=20, type=int)
    ap.add_argument("--top-n", dest="top_n", default=50, type=int)
    ap.add_argument("--impute", action="store_true")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_variant_impact(a)


if __name__ == "__main__":
    main()
