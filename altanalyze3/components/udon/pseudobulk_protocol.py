#!/usr/bin/env python3
"""
UDON pseudobulk protocol (optional): run UDON directly on a PRE-COMPUTED pseudobulk
h5ad (one row per Sample x cell-type), with control normalization against either
ALL controls (collective, legacy behavior) or SAMPLE-SPECIFIC matched controls.

Matched-control selection ports the validated control-matching method
(evaluation/RESULTS.md): per cell-type variable genes selected on CONTROLS only,
sex-matched Pearson correlation, depth-weighted aggregation across shared cell
types; each disease sample is normalized to its best-scoring sex-matched controls
(>= threshold, minimum 1). This captures batch/normal variation rather than AML
disease variance.

Input pseudobulk h5ad expectations (obs):
  - SAMPLE_COL   : patient/sample id            (default 'Sample')
  - CELLTYPE_COL : cell population/cluster       (default 'Hs-BM-titrated-reference-centroid')
  - ANNOT_COL    : disease/condition; controls == 'Control'   (default 'Annotation')
  - N_CELLS_COL  : cells per pseudobulk (weight) (default 'n_cells')
X is expected to be log2(CP10K+1) (UDON folds = log-difference of pseudobulk means).
"""
import os, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

CHRY_MARKERS = ["RPS4Y1","DDX3Y","EIF1AY","KDM5D","UTY","USP9Y","NLGN4Y","TMSB4Y",
                "PRKY","ZFY","TSPY1","BCORP1","PRY","TBL1Y","TXLNGY","AMELY"]
SEX_GENES = set(CHRY_MARKERS + ["XIST"])


def filter_udon_genes(genes, species="Hs", pc_path=None, min_pc_frac=0.5, logger=print):
    """Restrict UDON features to minimize batch effects, done up front on the feature set:
      - protein-coding only (requires a species annotation; applied ONLY if the majority of
        detected features match it, else the namespace likely differs and we don't restrict),
      - drop ribosomal genes (symbol begins with RPL/RPS, e.g. RPL*/RPS*; mouse Rpl/Rps too),
      - drop XIST and TSIX (X-inactivation; strong sex-batch drivers).
    Y-chromosome genes are ALWAYS kept (loss-of-Y in confirmed males is a real disease signal),
    which also exempts the Y-linked RPS4Y1/RPS4Y2 from the ribosomal rule. Returns kept genes."""
    import os
    if pc_path is None:
        pc_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ProteinCoding-Hs-Mm.txt")
    pc = set(pd.read_csv(pc_path, sep="\t", names=["g"])["g"].astype(str)) if os.path.exists(pc_path) else set()
    sym = [str(g).split("_")[0] for g in genes]
    frac = float(np.mean([s in pc for s in sym])) if pc else 0.0
    restrict_pc = bool(pc) and frac >= min_pc_frac
    if pc and not restrict_pc:
        logger(f"[gene-filter] only {frac:.0%} of features in the {species} protein-coding annotation "
               f"(<{min_pc_frac:.0%}); NOT restricting to protein-coding (namespace mismatch?)")
    elif restrict_pc:
        logger(f"[gene-filter] protein-coding ({species}): {frac:.0%} of features matched -> restricting")
    Y = {g.upper() for g in CHRY_MARKERS}
    keep, n_nc, n_rib, n_x = [], 0, 0, 0
    for g, s in zip(genes, sym):
        u = s.upper()
        if u in Y:                                       # always keep Y (loss-of-Y indicator)
            keep.append(g); continue
        if restrict_pc and s not in pc:
            n_nc += 1; continue
        if u in ("XIST", "TSIX"):
            n_x += 1; continue
        if u.startswith("RPL") or u.startswith("RPS"):   # ribosomal (RPS4Y* already kept above)
            n_rib += 1; continue
        keep.append(g)
    logger(f"[gene-filter] {len(genes)} -> {len(keep)} features (dropped {n_nc} non-coding, "
           f"{n_rib} ribosomal RPL/RPS, {n_x} XIST/TSIX; Y-chromosome genes kept)")
    return keep


def udon_restriction():
    """Shared UDON precursor config (so the non-CLI study-aware scripts honour --cell-type too).
    Reads env vars: UDON_CELL_TYPE (restrict to one cell type), UDON_GENE_FILTER (default on),
    UDON_SPECIES (default Hs), UDON_TAG (arbitrary output-folder suffix, e.g. to keep a gene-filtered
    run separate from an old one). Returns (cell_type|None, gene_filter:bool, species, path_suffix)."""
    import os
    ct = (os.environ.get("UDON_CELL_TYPE") or "").strip() or None
    gf = (os.environ.get("UDON_GENE_FILTER", "1").strip().lower() not in ("0", "false", "no", "off"))
    sp = os.environ.get("UDON_SPECIES", "Hs")
    tag = (os.environ.get("UDON_TAG") or "").strip()
    sfx = (("_" + ct.replace("/", "-")) if ct else "") + (("_" + tag) if tag else "")
    return ct, gf, sp, sfx


def restrict_folds(folds, fobs, celltype_col, logger=print):
    """Apply the env-configured cell-type restriction (pseudobulk columns) + gene filter (rows)
    to a fold matrix. fobs is indexed by pseudobulk. Returns (folds, fobs)."""
    ct, gf, sp, _ = udon_restriction()
    if ct:
        if celltype_col in fobs.columns:
            keep = fobs[celltype_col].astype(str).values == ct
        else:                                            # fall back to parsing 'celltype__Sample'
            keep = np.array([str(c).split("__", 1)[0] == ct for c in folds.columns])
        folds = folds.loc[:, folds.columns[keep]]; fobs = fobs.loc[folds.columns]
        logger(f"[restriction] cell type '{ct}': {folds.shape[1]} pseudobulks")
    if gf:
        folds = folds.loc[filter_udon_genes(list(folds.index), species=sp, logger=logger)]
    return folds, fobs


# --------------------------------------------------------------------------- #
# Sex prediction (chrY-vs-XIST dosage; validated 27/27 vs name-encoded sex)
# --------------------------------------------------------------------------- #
def predict_sample_sex(adata, sample_col, counts_adata=None, log_input=True):
    """Per-sample sex. Uses raw-count pseudobulk if given (preferred), else the
    main (log2) matrix. male if chrY_cp10k > XIST_cp10k."""
    src = counts_adata if counts_adata is not None else adata
    is_log = (counts_adata is None) and log_input
    samples = src.obs[sample_col].astype(str).values
    uS = pd.unique(samples)
    codes = pd.Categorical(samples, categories=uS).codes
    ind = sp.csr_matrix((np.ones(len(samples)), (codes, np.arange(len(samples)))),
                        shape=(len(uS), len(samples)))
    X = src.X
    if is_log:                      # undo log2(x+1) so aggregation is in linear CP10K-ish space
        X = X.copy()
        X = X.toarray() if sp.issparse(X) else np.asarray(X)
        X = np.power(2.0, X) - 1.0
        X = sp.csr_matrix(X)
    else:
        X = sp.csr_matrix(X)
    Xs = sp.csr_matrix(ind @ X)                       # samples x genes (summed)
    tot = np.asarray(Xs.sum(axis=1)).ravel()
    cp10k = Xs.multiply(1e4 / np.maximum(tot, 1)[:, None]).tocsr()
    yi = [src.var_names.get_loc(g) for g in CHRY_MARKERS if g in src.var_names]
    chry = np.asarray(cp10k[:, yi].sum(axis=1)).ravel() if yi else np.zeros(len(uS))
    xist = (np.asarray(cp10k[:, src.var_names.get_loc("XIST")].todense()).ravel()
            if "XIST" in src.var_names else np.zeros(len(uS)))
    eps = 0.1
    sexes = np.where(np.log2(chry+eps) - np.log2(xist+eps) > 0, "male", "female")
    return dict(zip(uS, sexes))


# --------------------------------------------------------------------------- #
# Matched-control selection (cell-type control-HVG, sex-matched, depth-weighted)
# --------------------------------------------------------------------------- #
def select_matched_controls(adata, sample_col, celltype_col, annot_col, n_cells_col,
                            sex_of, chemistry_of=None, control_label="Control", hvg_n=250,
                            min_controls_ct=4, min_shared_ct=2, topk=4, score_floor=0.40,
                            score_margin=0.05, logger=print):
    """Return {disease_sample: [(control_sample, score), ...]} and a long-form
    DataFrame. Controls are matched by **sex AND platform (chemistry)** — both hard
    constraints — and ranked by expression similarity (Pearson on control-HVG per
    cell type, depth-weighted across shared cell types). If chemistry_of is None,
    only sex is enforced (legacy)."""
    X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
    var_names = adata.var_names
    keep_gene = np.array([g not in SEX_GENES for g in var_names])
    pop = adata.obs[celltype_col].astype(str).values
    samp = adata.obs[sample_col].astype(str).values
    annot = adata.obs[annot_col].astype(str).values
    ncell = adata.obs[n_cells_col].astype(float).values if n_cells_col in adata.obs else np.ones(len(samp))
    is_ctrl = annot == control_label

    acc = {}
    for ct in pd.unique(pop):
        m = pop == ct
        ci = np.where(m & is_ctrl)[0]; di = np.where(m & ~is_ctrl)[0]
        if len(ci) < min_controls_ct or len(di) == 0:
            continue
        sub = X[np.where(m)[0]][:, keep_gene]
        order = np.where(m)[0]
        cmask = is_ctrl[order]
        cvar = sub[cmask].var(axis=0)
        hvg = np.argsort(cvar)[::-1][:hvg_n]
        sub = sub[:, hvg]
        # row-standardize for Pearson
        Z = sub - sub.mean(1, keepdims=True)
        Z /= (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-12)
        cZ = Z[cmask]; dZ = Z[~cmask]
        c_samp = samp[order][cmask]; d_samp = samp[order][~cmask]
        c_nc = ncell[order][cmask]; d_nc = ncell[order][~cmask]
        S = dZ @ cZ.T                                  # disease x control Pearson
        for i, ds in enumerate(d_samp):
            for j, cs in enumerate(c_samp):
                if sex_of.get(cs) != sex_of.get(ds):          # hard constraint: sex
                    continue
                if chemistry_of is not None and chemistry_of.get(cs) != chemistry_of.get(ds):
                    continue                                   # hard constraint: platform (chemistry)
                w = min(float(d_nc[i]), float(c_nc[j]))
                acc.setdefault((ds, cs), []).append((S[i, j], w))
    by_d = {}
    for (ds, cs), lst in acc.items():
        if len(lst) < min_shared_ct:
            continue
        sims = np.array([x[0] for x in lst]); ws = np.array([x[1] for x in lst])
        score = float(np.average(sims, weights=ws)) if ws.sum() > 0 else float(sims.mean())
        by_d.setdefault(ds, []).append((cs, score, len(lst)))
    # rank + threshold: keep controls >= max(floor, top-margin), min 1, cap topk
    selected = {}
    rows = []
    for ds, lst in by_d.items():
        lst.sort(key=lambda t: t[1], reverse=True)
        top = lst[0][1]
        thr = max(score_floor, top - score_margin)
        chosen = [t for t in lst if t[1] >= thr][:topk]
        if not chosen:
            chosen = [lst[0]]
        selected[ds] = [(c, s) for c, s, _ in chosen]
        for rank, (c, s, n) in enumerate(chosen, 1):
            rows.append({"disease_sample": ds, "control_sample": c, "score": round(s, 4),
                         "rank": rank, "n_shared_celltypes": n})
    logger(f"[match] matched {len(selected)} disease samples; "
           f"mean controls/sample = {np.mean([len(v) for v in selected.values()]):.2f}")
    return selected, pd.DataFrame(rows)


def write_control_annotation(selected, path):
    """2-column TSV: disease_sample <tab> ctrl (score),ctrl (score),... (matches the
    evaluation deliverable format; also the file UDON can read back via --control-annotation)."""
    with open(path, "w") as f:
        f.write("disease_sample\tcontrols\n")
        for ds in sorted(selected):
            cell = ",".join(f"{c} ({s:.3f})" for c, s in selected[ds])
            f.write(f"{ds}\t{cell}\n")


def read_control_annotation(path):
    """Parse the 2-column control annotation file -> {disease_sample: [control,...]}."""
    sel = {}
    with open(path) as f:
        next(f)
        for line in f:
            ds, _, rest = line.rstrip("\n").partition("\t")
            ctrls = []
            for tok in rest.split(","):
                tok = tok.strip()
                if not tok:
                    continue
                ctrls.append(tok.split(" (")[0].strip())
            if ds:
                sel[ds] = ctrls
    return sel


# --------------------------------------------------------------------------- #
# Fold matrix from pre-made pseudobulks
# --------------------------------------------------------------------------- #
def build_fold_matrix(adata, sample_col, celltype_col, annot_col, control_label="Control",
                      mode="matched", selected_controls=None, drop_unmatched=True, logger=print):
    """Return (folds_df genes x disease-pseudobulks, obs_df for the columns).
    mode='matched' uses selected_controls[disease_sample]; mode='collective' uses all
    controls of the cell type. fold = disease_pb - mean(control_pbs) in log2 space.

    drop_unmatched=True (matched mode): if a disease pseudobulk's matched controls have
    NO pseudobulk in its cell type, the pseudobulk is DROPPED (not analyzed) rather than
    falling back to the collective control mean (collective fallback reintroduces batch)."""
    X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
    pop = adata.obs[celltype_col].astype(str).values
    samp = adata.obs[sample_col].astype(str).values
    annot = adata.obs[annot_col].astype(str).values
    is_ctrl = annot == control_label
    genes = list(adata.var_names)

    # index control rows by (celltype, sample)
    ctrl_idx = {}
    for r in np.where(is_ctrl)[0]:
        ctrl_idx.setdefault((pop[r], samp[r]), r)
    ctrl_by_ct = {}
    for (ct, s), r in ctrl_idx.items():
        ctrl_by_ct.setdefault(ct, []).append((s, r))

    fold_cols = {}; col_meta = []
    n_skipped = n_matched = n_fallback = 0
    for r in np.where(~is_ctrl)[0]:
        ct, s = pop[r], samp[r]
        baseline_kind = "collective"
        if mode == "collective":
            base_rows = [rr for _, rr in ctrl_by_ct.get(ct, [])]
        else:
            # selected_controls[s] may be a list of control names OR of (name, score) tuples
            raw = (selected_controls or {}).get(s, [])
            wanted = set(c[0] if isinstance(c, (tuple, list)) else c for c in raw)
            base_rows = [rr for cs, rr in ctrl_by_ct.get(ct, []) if cs in wanted]
            if base_rows:
                baseline_kind = "matched"; n_matched += 1
            elif drop_unmatched:                     # no matched control in this cell type -> DROP (not analyzed)
                n_fallback += 1
                continue
            else:                                    # legacy: fall back to collective control mean
                base_rows = [rr for _, rr in ctrl_by_ct.get(ct, [])]
                n_fallback += 1
        if not base_rows:
            n_skipped += 1
            continue
        baseline = X[base_rows].mean(axis=0)
        col = f"{ct}__{s}"
        fold_cols[col] = X[r] - baseline
        col_meta.append({"pseudobulk": col, sample_col: s, celltype_col: ct,
                         annot_col: annot[r], "n_controls_used": len(base_rows),
                         "baseline_kind": baseline_kind})
    folds = pd.DataFrame(fold_cols, index=genes)
    # remove_negatives (shift each gene so min >= 0), drop NA columns (legacy parity)
    folds = folds.sub(folds.min(axis=1), axis=0)
    folds = folds.dropna(axis=1)
    obs_df = pd.DataFrame(col_meta).set_index("pseudobulk").loc[folds.columns]
    logger(f"[folds] {folds.shape[0]} genes x {folds.shape[1]} disease pseudobulks "
           f"(mode={mode}); matched-baseline={n_matched}, collective-fallback={n_fallback}, "
           f"skipped={n_skipped}")
    return folds, obs_df


def assemble_udon_adata(folds, obs_df):
    """Build the AnnData UDON expects: var=genes, varm['pseudobulk_folds']=folds,
    X=folds.T (pseudobulks x genes), obs=column metadata."""
    a = ad.AnnData(X=folds.T.values.astype(np.float32),
                   obs=obs_df.copy(),
                   var=pd.DataFrame(index=folds.index))
    a.varm["pseudobulk_folds"] = folds
    return a
