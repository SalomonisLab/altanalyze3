#!/usr/bin/env python3
"""
UDON pseudobulk CLI — run the UDON disease-heterogeneity protocol directly on a
PRE-COMPUTED pseudobulk h5ad, with optional SAMPLE-SPECIFIC matched-control
normalization, GO-Elite enrichment, and a results folder + log.

Example:
  python3 run_udon_pseudobulk.py \
    --pseudobulk  /path/pseudobulk_scaled_log2_hashed.h5ad \
    --counts-pseudobulk /path/pseudobulk_counts_hashed.h5ad \
    --control-mode matched --species Hs

Outputs (in <pseudobulk dir>/UDON_results/ by default):
  udon_run.log, control_annotation.tsv, matched_controls_long.tsv, sample_sex.tsv,
  udon_core/{marker_genes_all.txt, marker_heatmap.txt, udon_clusters.txt},
  udon_result.h5ad, marker_heatmap.pdf, goelite/GOElite_UDON*.tsv(.pdf), RUN_SUMMARY.md
"""
import os, sys, argparse, logging, re
import numpy as np, pandas as pd, anndata as ad

UDON_DIR = os.path.dirname(os.path.abspath(__file__))
if UDON_DIR not in sys.path:
    sys.path.insert(0, UDON_DIR)

import pseudobulk_protocol as P


def setup_logger(logpath):
    lg = logging.getLogger("udon"); lg.setLevel(logging.INFO); lg.handlers = []
    fmt = logging.Formatter("%(asctime)s %(levelname)s | %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(logpath); fh.setFormatter(fmt); lg.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout); sh.setFormatter(fmt); lg.addHandler(sh)
    return lg


SUSPECT_RE = re.compile(r"TotalSeq|empty", re.IGNORECASE)
def flag_suspect(samples, annots):
    return np.array([bool(SUSPECT_RE.search(str(s))) or str(a) == "0"
                     for s, a in zip(samples, annots)])


def main():
    ap = argparse.ArgumentParser(description="UDON pseudobulk protocol")
    ap.add_argument("--pseudobulk", required=True, help="pseudobulk h5ad (X = log2(CP10K+1))")
    ap.add_argument("--counts-pseudobulk", default=None, help="raw-count pseudobulk h5ad (sex calling)")
    ap.add_argument("--output-dir", default=None)
    ap.add_argument("--species", default="Hs")
    ap.add_argument("--sample-col", default="Sample")
    ap.add_argument("--celltype-col", default="Hs-BM-titrated-reference-centroid")
    ap.add_argument("--annot-col", default="Annotation")
    ap.add_argument("--ncells-col", default="n_cells")
    ap.add_argument("--control-label", default="Control")
    ap.add_argument("--control-mode", choices=["matched", "collective", "annotation"], default="matched")
    ap.add_argument("--control-annotation", default=None, help="2-col TSV for --control-mode annotation")
    ap.add_argument("--sample-metadata", default=None,
                    help="TSV (index=Sample) with a 'chemistry' column; enables sex+platform matching")
    ap.add_argument("--match-platform", action="store_true", default=True,
                    help="require matched controls to share platform/chemistry (default on)")
    ap.add_argument("--no-match-platform", dest="match_platform", action="store_false")
    ap.add_argument("--hvg-n", type=int, default=250)
    ap.add_argument("--topk", type=int, default=4)
    ap.add_argument("--score-floor", type=float, default=0.40)
    ap.add_argument("--score-margin", type=float, default=0.05)
    ap.add_argument("--rank", type=int, default=None)
    ap.add_argument("--modality", choices=["rna", "adt"], default="rna",
                    help="adt: skip RNA-specific protein-coding/non-coding gene filters and use all features")
    ap.add_argument("--fast", action="store_true",
                    help="use validated-fast feature selection (vectorized) + SVD rank instead of the "
                         "slow per-gene loop + full-eig rank (nimfa clustering unchanged)")
    ap.add_argument("--no-goelite", action="store_true")
    ap.add_argument("--keep-suspect", action="store_true", help="do not drop TotalSeq/empty/'0' pseudobulks")
    ap.add_argument("--subset-celltypes", default=None, help="comma list to restrict cell types (smoke test)")
    ap.add_argument("--cell-type", default=None,
                    help="restrict the FULL analysis to ONE cell type (e.g. 'MPP-MEP'); output dir gets a "
                         "_<celltype> suffix unless --output-dir is given")
    ap.add_argument("--no-gene-filter", action="store_true",
                    help="skip the up-front UDON gene filter (protein-coding only; drop RPL/RPS ribosomal + "
                         "XIST/TSIX; keep Y) that reduces batch effects (RNA modality only)")
    args = ap.parse_args()

    pb_path = os.path.abspath(args.pseudobulk)
    if args.output_dir:
        outdir = os.path.abspath(args.output_dir)
    else:
        base = os.path.join(os.path.dirname(pb_path), "UDON_results")
        outdir = f"{base}_{args.cell_type.replace('/', '-')}" if args.cell_type else base
    os.makedirs(outdir, exist_ok=True)
    lg = setup_logger(os.path.join(outdir, "udon_run.log"))
    log = lambda m: lg.info(m)
    log(f"=== UDON pseudobulk protocol ===")
    log(f"input: {pb_path}")
    log(f"output: {outdir}")
    log(f"control-mode: {args.control_mode} | species: {args.species}")

    adata = ad.read_h5ad(pb_path)
    log(f"loaded pseudobulk: {adata.shape} (obs cols: {list(adata.obs.columns)})")
    for c in (args.sample_col, args.celltype_col, args.annot_col):
        if c not in adata.obs.columns:
            log(f"FATAL: obs column '{c}' missing"); sys.exit(2)

    # drop suspect non-patient pseudobulks
    if not args.keep_suspect:
        susp = flag_suspect(adata.obs[args.sample_col].astype(str).values,
                            adata.obs[args.annot_col].astype(str).values)
        if susp.any():
            bad = sorted(set(adata.obs[args.sample_col].astype(str).values[susp]))
            log(f"dropping {int(susp.sum())} suspect pseudobulks from samples: {bad}")
            adata = adata[~susp].copy()
    if args.cell_type and args.cell_type not in set(adata.obs[args.celltype_col].astype(str)):
        # validate now; the actual restriction is applied to the FOLDS (control matching needs all
        # cell types to score sample similarity, so we must NOT subset the AnnData before matching).
        log(f"FATAL: --cell-type '{args.cell_type}' not found. Available (first 40): "
            f"{sorted(set(adata.obs[args.celltype_col].astype(str)))[:40]}"); sys.exit(2)
    if args.subset_celltypes:
        cts = set(args.subset_celltypes.split(","))
        adata = adata[adata.obs[args.celltype_col].astype(str).isin(cts)].copy()
        log(f"SUBSET to cell types {cts}: {adata.shape}")

    annots = adata.obs[args.annot_col].astype(str).values
    n_ctrl = int((annots == args.control_label).sum()); n_dis = int((annots != args.control_label).sum())
    log(f"control pseudobulks: {n_ctrl} | disease pseudobulks: {n_dis}")

    # ---- control selection ----
    selected = None
    if args.control_mode == "matched":
        counts_ad = ad.read_h5ad(os.path.abspath(args.counts_pseudobulk)) if args.counts_pseudobulk else None
        log("predicting sample sex (chrY-vs-XIST dosage)...")
        sex_of = P.predict_sample_sex(adata, args.sample_col, counts_adata=counts_ad, log_input=True)
        pd.Series(sex_of, name="predicted_sex").rename_axis("Sample").to_csv(os.path.join(outdir, "sample_sex.tsv"), sep="\t")
        log(f"sex: male={sum(v=='male' for v in sex_of.values())} female={sum(v=='female' for v in sex_of.values())}")
        # platform (chemistry) per sample, for sex+platform matching
        chem_of = None
        if args.match_platform:
            if not args.sample_metadata or not os.path.exists(args.sample_metadata):
                log("FATAL: --match-platform needs --sample-metadata TSV with a 'chemistry' column "
                    "(use --no-match-platform to match by sex only)"); sys.exit(2)
            smeta = pd.read_csv(args.sample_metadata, sep="\t", index_col=0)
            chem_of = smeta["chemistry"].astype(str).to_dict()
            log(f"matching by sex + platform (chemistry); {len(set(chem_of.values()))} platforms")
        log("selecting matched controls per disease sample...")
        selected, long_df = P.select_matched_controls(
            adata, args.sample_col, args.celltype_col, args.annot_col, args.ncells_col,
            sex_of, chemistry_of=chem_of, control_label=args.control_label, hvg_n=args.hvg_n,
            topk=args.topk, score_floor=args.score_floor, score_margin=args.score_margin, logger=log)
        long_df.to_csv(os.path.join(outdir, "matched_controls_long.tsv"), sep="\t", index=False)
        P.write_control_annotation(selected, os.path.join(outdir, "control_annotation.tsv"))
        fold_mode = "matched"
    elif args.control_mode == "annotation":
        if not args.control_annotation:
            log("FATAL: --control-annotation required for control-mode=annotation"); sys.exit(2)
        selected = P.read_control_annotation(os.path.abspath(args.control_annotation))
        log(f"read control annotation: {len(selected)} disease samples")
        fold_mode = "matched"
    else:
        fold_mode = "collective"

    # ---- fold matrix ----
    log("computing pseudobulk folds (disease - control baseline, log2 space)...")
    folds, fold_obs = P.build_fold_matrix(
        adata, args.sample_col, args.celltype_col, args.annot_col,
        control_label=args.control_label, mode=fold_mode, selected_controls=selected, logger=log)
    if args.cell_type:                                    # restrict to ONE cell type AFTER matching+folds
        keep = fold_obs[args.celltype_col].astype(str) == args.cell_type
        folds = folds.loc[:, fold_obs.index[keep]]; fold_obs = fold_obs[keep]
        log(f"RESTRICTED folds to cell type '{args.cell_type}': {folds.shape[1]} pseudobulks")
    if args.modality == "rna" and not args.no_gene_filter:   # up-front batch-reduction gene filter
        folds = folds.loc[P.filter_udon_genes(list(folds.index), species=args.species, logger=log)]
    udon = P.assemble_udon_adata(folds, fold_obs)
    log(f"UDON adata assembled: {udon.shape} (varm pseudobulk_folds: {folds.shape})")

    # ---- UDON feature selection + clustering (run from UDON_DIR for relative resources) ----
    prev_cwd = os.getcwd(); os.chdir(UDON_DIR)
    try:
        from feature_selection import feature_selection_wrapper
        from clustering_wrapper import clustering_wrapper
        rank = args.rank
        if args.modality == "adt":
            # ADT: RNA gene filters (protein-coding, HLA/Y/dot strips) would wrongly drop
            # surface markers (CD*, HLA.ABC, ...). Use ALL features for NMF.
            log(f"modality=adt: skipping RNA gene filters; using all {udon.n_vars} ADT features")
            udon.var["correlated_genes"] = True
            if args.fast and rank is None:
                import fast_feature_selection as FFS
                rank = FFS.fast_rank(folds)
                log(f"fast rank (randomized SVD): {rank}")
        elif args.fast:
            import fast_feature_selection as FFS
            log("FAST feature selection (vectorized; validated-identical to original)...")
            pc = None
            pcf = os.path.join(UDON_DIR, "ProteinCoding-Hs-Mm.txt")
            if os.path.exists(pcf):
                gcol = pd.read_csv(pcf, sep="\t", names=["g"])["g"]
                pc = {x for x in gcol if (str(x).isupper() if args.species == "Hs" else True)}
            fs = FFS.fast_feature_selection(folds, pc)
            # subset to the same gene universe the original markerFinder sees (name+protein filtered)
            udon = P.assemble_udon_adata(folds.loc[fs["after_name_filter"]], fold_obs)
            udon.var["correlated_genes"] = udon.var_names.isin(fs["correlated_genes"])
            log(f"fast FS: {len(fs['after_name_filter'])} genes after filter, "
                f"correlated_genes={int(udon.var['correlated_genes'].sum())}")
            if rank is None:
                rank = FFS.fast_rank(folds.loc[fs["correlated_genes"]])
                log(f"fast rank (randomized SVD): {rank}")
        else:
            log("UDON feature selection (original)...")
            udon = feature_selection_wrapper(udon, species=args.species)
        core_out = os.path.join(outdir, "udon_core")
        if os.path.exists(core_out):
            import shutil; shutil.rmtree(core_out)
        log("UDON NMF clustering + marker finding...")
        udon = clustering_wrapper(udon, output_filename=core_out, rank=rank)
    finally:
        os.chdir(prev_cwd)

    clusters = udon.uns.get("udon_clusters")
    n_clusters = clusters["cluster"].nunique() if clusters is not None and len(clusters) else 0
    log(f"UDON clusters found: {n_clusters}")

    # ---- save result h5ad. uns holds pandas DataFrames (clusters/markers/heatmap) that
    #      can trip h5ad serialization; fall back to a slim object (X/obs/var/varm) since the
    #      key tables are already written as TSVs in udon_core/. ----
    h5_path = os.path.join(outdir, "udon_result.h5ad")
    try:
        u = udon.copy()
        u.obs.index = u.obs.index.astype(str)
        u.var.index = u.var.index.astype(str)
        for c in u.obs.columns:
            if u.obs[c].dtype == object:
                u.obs[c] = u.obs[c].astype(str)
        u.write_h5ad(h5_path)
        log("wrote udon_result.h5ad")
    except Exception as e:
        try:
            slim = ad.AnnData(X=udon.X, obs=udon.obs.astype(str), var=udon.var.copy())
            slim.write_h5ad(h5_path)
            log(f"wrote slim udon_result.h5ad (uns tables in udon_core/; full write failed: {e})")
        except Exception as e2:
            log(f"h5ad write skipped ({e2}); key tables already in udon_core/")

    # ---- heatmap: canonical MarkerFinder layout/colors (visualizations.plot_markers_df),
    #      now with the dense quadmesh rasterized -> writes marker_heatmap.png + .pdf
    #      (the layout is unchanged; only the >200MB vector blow-up is fixed). ----
    try:
        os.chdir(UDON_DIR)
        from visualizations import plot_markers_df
        plot_markers_df(udon.uns["marker_heatmap"], udon.uns["udon_marker_genes_top_n"],
                        udon.uns["udon_clusters"], os.path.join(outdir, "marker_heatmap.pdf"))
        log("wrote marker_heatmap.png + marker_heatmap.pdf (MarkerFinder layout, rasterized)")
        # standard UDON binary heatmaps: donor x cluster + cell-type x cluster
        from udon_binary_heatmaps import make_udon_binary_heatmaps
        study_of = None
        if args.sample_metadata and os.path.exists(args.sample_metadata):
            _sm = pd.read_csv(args.sample_metadata, sep="\t", index_col=0)
            if "Dataset" in _sm.columns:
                study_of = _sm["Dataset"].astype(str).to_dict()
        donor_of = None
        _pbdir = os.path.dirname(os.path.abspath(pb_path))
        for _c in (os.path.join(_pbdir, "UDON", "AML_harmonized_metadata.xlsx"),
                   os.path.join(_pbdir, "AML_harmonized_metadata.xlsx")):
            if os.path.exists(_c):
                _cl = pd.read_excel(_c, sheet_name="Clinical_Metadata")
                donor_of = dict(zip(_cl["Sample"].astype(str), _cl["Donor_ID"].astype(str))); break
        make_udon_binary_heatmaps(udon.uns["udon_clusters"], outdir, donor_of=donor_of, study_of=study_of)
        log("wrote celltype_cluster_heatmap" + (" + donor_cluster_heatmap" if donor_of else ""))
    except Exception as e:
        log(f"heatmap skipped: {e}")
    finally:
        os.chdir(prev_cwd)

    # ---- GO-Elite ----
    if not args.no_goelite and n_clusters > 0:
        try:
            from goelite_enrichment import run_goelite_on_udon
            run_goelite_on_udon(udon, os.path.join(outdir, "goelite"), species=args.species, logger=log)
        except Exception as e:
            log(f"GO-Elite skipped: {e}")

    # ---- run summary ----
    with open(os.path.join(outdir, "RUN_SUMMARY.md"), "w") as f:
        f.write(f"# UDON pseudobulk run\n\n")
        f.write(f"- input: `{pb_path}`\n- control mode: **{args.control_mode}**\n")
        f.write(f"- control pseudobulks: {n_ctrl}; disease pseudobulks: {n_dis}\n")
        f.write(f"- fold matrix: {folds.shape[0]} genes x {folds.shape[1]} disease pseudobulks\n")
        f.write(f"- UDON clusters: {n_clusters}\n")
        if selected is not None:
            f.write(f"- matched controls/sample (mean): {np.mean([len(v) for v in selected.values()]):.2f}\n")
    log("=== DONE ===")


if __name__ == "__main__":
    main()
