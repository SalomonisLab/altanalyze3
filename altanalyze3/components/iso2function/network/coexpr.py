"""Expression- and assay-informed condition-specific interaction networks.

Computes, for a donor-paired contrast between two groups (cell-state x covariate), the per-edge
interaction "activity" EA_ij = w_ij * sqrt(P_i * P_j), where P is isoform/gene abundance and w_ij is the
assay affinity (isoform_interactions.activity_score). It then tests which edges differ between the two
groups (donor-level, moderated paired statistics), evaluates PDI regulatory consequence via target-gene
DE coherence (recovering the unencoded activation/repression direction), and ranks the TF isoforms that
drive the network differences.

Defensibility: the replication unit is the biological donor (NOT the cell). Each group carries one
isoform/gene pseudobulk per donor; statistics are paired across donors with empirical-Bayes variance
moderation (limma-style) appropriate for small n. Effect sizes accompany every p-value. The assay
affinity is a soft prior (measured for ~190 edges, median elsewhere) so expression carries the signal.

Inputs (AltAnalyze3 long-read workflow): isoform_combined_pseudo_cluster_tpm-filtered.h5ad and
gene_combined_pseudo_cluster_tpm-filtered.h5ad, whose obs are '{cell_state}.{donor_library}-isoform/-gene'
and whose var are 'ENSG:ENST' (isoform) / 'ENSG' (gene). Edges from data/isoform_interactions.txt.
"""

import os
import gzip
import logging

import numpy as np
import pandas as pd

from .. import config

logger = logging.getLogger(__name__)
EPS = 1e-3                       # pseudo-count for log of activity / CPM
DETECT_TPM = 1.0                # a molecule counts as "expressed" in a donor at CPM >= this

# Dataset-specific input names (overridable per deployment, e.g. cohorts that label files cpm vs tpm).
ISO_H5AD = "isoform_combined_pseudo_cluster_tpm-filtered.h5ad"
GENE_H5AD = "gene_combined_pseudo_cluster_tpm-filtered.h5ad"
META = "mds_metadata_bam.txt"
EDGE_TYPE = None                 # None = all interactions; "PDI" or "PPI" to restrict per deployment


# --------------------------------------------------------------------------- loaders
def _read_h5ad(path):
    import anndata as ad
    return ad.read_h5ad(path)


def load_symbol_to_ensg(edges):
    """symbol -> ENSG, from the MANE summary (authoritative) plus the edge table's own Symbol/ENSG."""
    s2g = {}
    p = config.MANE_SUMMARY
    if os.path.exists(p):
        with gzip.open(p, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            h = {c: i for i, c in enumerate(hdr)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > max(h["symbol"], h["Ensembl_Gene"]):
                    s2g.setdefault(f[h["symbol"]], f[h["Ensembl_Gene"]].split(".")[0])
    for sym, g in zip(edges["Symbol"], edges["ENSG"]):
        if sym and g and g.startswith("ENSG"):
            s2g.setdefault(sym, g)
    return s2g


def load_mane_ensg2enst():
    """{base_ENSG -> base MANE-Select ENST} (the gene's productive reference transcript)."""
    out = {}
    p = config.MANE_SUMMARY
    if os.path.exists(p):
        with gzip.open(p, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t"); h = {c: i for i, c in enumerate(hdr)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > max(h["Ensembl_Gene"], h["Ensembl_nuc"], h["MANE_status"]) \
                        and f[h["MANE_status"]] == "MANE Select":
                    out.setdefault(f[h["Ensembl_Gene"]].split(".")[0], f[h["Ensembl_nuc"]].split(".")[0])
    return out


def reference_isoform_expr(iso_tpm, mane=None):
    """Score the TARGET at ISOFORM resolution: a (state,donor)-by-ENSG matrix holding each gene's
    reference (MANE-Select) transcript TPM, taken from the isoform matrix. The MANE isoform is the
    productive, Not-NMD reference (so non-productive/NMD/intron-retention transcripts never contribute
    to the partner's abundance). Genes whose MANE isoform is absent from the isoform quantification get
    NO column -> their interactions are not scored (only isoform-resolved-on-both-sides edges remain)."""
    mane = mane or load_mane_ensg2enst()
    cols = {ensg: iso_tpm[enst].to_numpy() for ensg, enst in mane.items() if enst in iso_tpm.columns}
    ref = pd.DataFrame(cols, index=iso_tpm.index)
    logger.info("[coexpr] reference-isoform target expression: %d genes have an expressed MANE isoform",
                ref.shape[1])
    return ref


def load_groups(dataset_dir, states, iso_h5ad=None, gene_h5ad=None, keep_donors=None):
    iso_h5ad = iso_h5ad or ISO_H5AD
    gene_h5ad = gene_h5ad or GENE_H5AD
    """Return (iso_tpm, gene_tpm, donors_by_state). Each matrix is a DataFrame indexed by (state, donor)
    with columns = ENST (isoform TPM) or ENSG (gene TPM), restricted to the requested states."""
    iso = _read_h5ad(os.path.join(dataset_dir, iso_h5ad))
    gene = _read_h5ad(os.path.join(dataset_dir, gene_h5ad))

    def parse(adata, suffix):
        rows = {}
        for i, name in enumerate(adata.obs_names):
            if not name.endswith(suffix):
                continue
            state, _, donor = name[: -len(suffix)].partition(".")
            if (states is None or state in states) and (keep_donors is None or donor in keep_donors):
                rows[(state, donor)] = i      # states/donors=None -> all; else subset BEFORE densifying
        return rows

    iso_rows, gene_rows = parse(iso, "-isoform"), parse(gene, "-gene")
    # isoform columns -> base ENST (keep first occurrence; novel 'ENSG:coord.donor' cols are ignored)
    enst_col = {}
    for j, v in enumerate(iso.var_names):
        if ":" in v:
            tail = v.split(":", 1)[1]
            if tail.startswith("ENST"):
                enst_col.setdefault(tail.split(".")[0], j)
    ensg_col = {str(v).split(".")[0]: j for j, v in enumerate(gene.var_names)}

    def densify(X, row_ids, col_ids):         # subset sparse rows+cols, THEN densify (memory-safe)
        sub = X[row_ids][:, col_ids]
        return sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)

    ikeys = sorted(iso_rows)
    iso_tpm = pd.DataFrame(densify(iso.X, [iso_rows[k] for k in ikeys], list(enst_col.values())),
                           index=pd.MultiIndex.from_tuples(ikeys, names=["state", "donor"]),
                           columns=list(enst_col.keys()))
    gkeys = sorted(gene_rows)
    gene_tpm = pd.DataFrame(densify(gene.X, [gene_rows[k] for k in gkeys], list(ensg_col.values())),
                            index=pd.MultiIndex.from_tuples(gkeys, names=["state", "donor"]),
                            columns=list(ensg_col.keys()))

    all_states = states if states is not None else sorted({st for (st, d) in iso_rows})
    donors_by_state = {s: sorted({d for (st, d) in iso_rows if st == s}) for s in all_states}
    logger.info("[coexpr] loaded states=%s donors=%s | isoform ENSTs=%d gene ENSGs=%d",
                states, {s: len(d) for s, d in donors_by_state.items()}, iso_tpm.shape[1], gene_tpm.shape[1])
    return iso_tpm, gene_tpm, donors_by_state


def load_max_isoform_expr(dataset_dir, states, iso_h5ad=None, keep_donors=None):
    """Per-gene MAX isoform CPM: for each ENSG, the maximum over ALL its expressed isoforms (known +
    novel) per (state, donor) — the dominant isoform's abundance. Used as the abundance of a target /
    interacting gene when it has several isoforms (values are CPM, not TPM)."""
    import collections
    iso_h5ad = iso_h5ad or ISO_H5AD
    iso = _read_h5ad(os.path.join(dataset_dir, iso_h5ad))
    rows = {}
    for i, name in enumerate(iso.obs_names):
        if not name.endswith("-isoform"):
            continue
        state, _, donor = name[: -len("-isoform")].partition(".")
        if (states is None or state in states) and (keep_donors is None or donor in keep_donors):
            rows[(state, donor)] = i
    gene_cols = collections.defaultdict(list)
    for j, v in enumerate(iso.var_names):
        gene_cols[str(v).split(":")[0]].append(j)
    keys = sorted(rows)
    sub = iso.X[[rows[k] for k in keys]]                  # subset rows BEFORE densifying
    Xsub = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
    data = {g: Xsub[:, cols].max(axis=1) for g, cols in gene_cols.items()}
    out = pd.DataFrame(data, index=pd.MultiIndex.from_tuples(keys, names=["state", "donor"]))
    logger.info("[coexpr] per-gene max-isoform CPM: %d genes x %d (state,donor)", out.shape[1], out.shape[0])
    return out


def load_edges(interactions_txt=None):
    edges = pd.read_csv(interactions_txt or config.INTERACTIONS_TXT, sep="\t", dtype=str).fillna("")
    edges["activity_score"] = pd.to_numeric(edges["activity_score"], errors="coerce").fillna(0.0)
    if EDGE_TYPE:                # restrict to a single interaction type (e.g. PDI = regulatory only)
        edges = edges[edges["interaction_type"] == EDGE_TYPE].reset_index(drop=True)
    return edges


# --------------------------------------------------------------------------- edge activity
def edge_activity(edges, iso_tpm, gene_tpm, sym2ensg, states, donors):
    """Per-edge interaction activity EA = w * sqrt(P_i * P_j) for every (state, donor).
    P_i = isoform TPM of the source isoform (best_ENST); P_j = gene TPM of the target gene.
    Returns a long DataFrame: one row per (edge, state, donor) with EA and the two potentials.
    Edges whose source isoform or target gene is absent from the matrices are dropped."""
    keys = [(s, d) for s in states for d in donors[s]]
    iso_lookup = {k: iso_tpm.loc[k] for k in keys if k in iso_tpm.index}
    gene_lookup = {k: gene_tpm.loc[k] for k in keys if k in gene_tpm.index}

    # unique (source ENST, target ENSG, type, w) edges; carry Symbol + source_isoform_id for reporting
    e = edges.copy()
    e["tgt_ensg"] = e["target"].map(sym2ensg).fillna("")
    e = e[(e["best_ENST"] != "") & (e["tgt_ensg"] != "")]
    e = e.drop_duplicates(["best_ENST", "tgt_ensg", "interaction_type"])

    iso_ok = set(iso_tpm.columns)
    gene_ok = set(gene_tpm.columns)
    e = e[e["best_ENST"].isin(iso_ok) & e["tgt_ensg"].isin(gene_ok)]

    rows = []
    for r in e.itertuples(index=False):
        for (s, d) in keys:
            il, gl = iso_lookup.get((s, d)), gene_lookup.get((s, d))
            if il is None or gl is None:
                continue
            pi = float(il.get(r.best_ENST, 0.0))
            pj = float(gl.get(r.tgt_ensg, 0.0))
            ea = r.activity_score * np.sqrt(max(pi, 0.0) * max(pj, 0.0))
            rows.append((r.Symbol, r.source_isoform_id, r.best_ENST, r.target, r.tgt_ensg,
                         r.interaction_type, r.activity_score, s, d, pi, pj, ea))
    cols = ["Symbol", "source_isoform_id", "source_ENST", "target", "target_ENSG", "interaction_type",
            "affinity", "state", "donor", "P_source", "P_target", "EA"]
    ea_long = pd.DataFrame(rows, columns=cols)
    logger.info("[coexpr] edge activity: %d edges x %d (state,donor) = %d rows",
                e.shape[0], len(keys), len(ea_long))
    return ea_long


# --------------------------------------------------------------------------- moderated paired stats
def _moderated_paired(diff_matrix, prior_df=4.0):
    """limma-style empirical-Bayes moderated paired t-test. diff_matrix: edges x donors of paired
    differences (log2 group A - group B). Returns (mean_diff, t, df, p) arrays."""
    from scipy import stats
    n = diff_matrix.shape[1]
    mean_d = np.nanmean(diff_matrix, axis=1)
    var_d = np.nanvar(diff_matrix, axis=1, ddof=1)
    var_d = np.where(np.isfinite(var_d), var_d, np.nan)
    s2_prior = np.nanmedian(var_d[var_d > 0]) if np.any(var_d > 0) else 1.0
    s2_post = (prior_df * s2_prior + (n - 1) * np.nan_to_num(var_d)) / (prior_df + (n - 1))
    se = np.sqrt(s2_post / n)
    t = np.where(se > 0, mean_d / se, 0.0)
    dfree = prior_df + (n - 1)
    p = 2 * stats.t.sf(np.abs(t), dfree)
    return mean_d, t, dfree, p


def _bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    idx = np.where(ok)[0]
    if not len(idx):
        return q
    order = idx[np.argsort(p[idx])]
    m = len(order)
    ranked = p[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.clip(ranked, 0, 1)
    return q


# --------------------------------------------------------------------------- differential edges
def differential_edges(ea_long, state_a, state_b, donors_common, min_detect=2, p_alpha=0.05):
    """Donor-paired differential activity of each edge between state_a and state_b. Returns a per-edge
    table with log2 fold-change, moderated p, BH-FDR, detection, attribution (source vs target driven),
    and a class (gained/lost/rewired/quantitative/unchanged)."""
    key = ["source_isoform_id", "source_ENST", "target", "target_ENSG", "interaction_type", "affinity",
           "Symbol"]
    ea_long = ea_long[ea_long["donor"].isin(donors_common)].copy()
    piv = ea_long.pivot_table(index=key, columns=["state", "donor"], values="EA")
    psrc = ea_long.pivot_table(index=key, columns=["state", "donor"], values="P_source")
    ptgt = ea_long.pivot_table(index=key, columns=["state", "donor"], values="P_target")

    A = np.log2(piv[state_a].reindex(columns=donors_common).to_numpy() + EPS)
    B = np.log2(piv[state_b].reindex(columns=donors_common).to_numpy() + EPS)
    diff = A - B
    mean_d, t, dfree, p = _moderated_paired(diff)
    fdr = _bh(p)

    eaA, eaB = piv[state_a].mean(axis=1).to_numpy(), piv[state_b].mean(axis=1).to_numpy()
    detA = (piv[state_a].reindex(columns=donors_common).to_numpy() > 0).sum(axis=1)
    detB = (piv[state_b].reindex(columns=donors_common).to_numpy() > 0).sum(axis=1)
    # attribution: how much of the activity change is source-isoform vs target driven (log2 fold-changes)
    d_src = (np.log2(psrc[state_a].mean(axis=1).to_numpy() + EPS)
             - np.log2(psrc[state_b].mean(axis=1).to_numpy() + EPS))
    d_tgt = (np.log2(ptgt[state_a].mean(axis=1).to_numpy() + EPS)
             - np.log2(ptgt[state_b].mean(axis=1).to_numpy() + EPS))

    out = piv.index.to_frame(index=False)
    out["mean_EA_%s" % state_a] = eaA
    out["mean_EA_%s" % state_b] = eaB
    out["log2FC"] = mean_d
    out["t"] = t
    out["p"] = p
    out["FDR"] = fdr
    out["detect_%s" % state_a] = detA
    out["detect_%s" % state_b] = detB
    out["dlog2_source"] = d_src
    out["dlog2_target"] = d_tgt
    out["driver"] = np.where(np.abs(d_src) >= np.abs(d_tgt), "source_isoform", "target_gene")

    keep = (detA >= min_detect) | (detB >= min_detect)
    out = out[keep].reset_index(drop=True)
    out["significant_FDR"] = out["FDR"] < 0.1
    out["significant"] = out["p"] < p_alpha          # lax raw-p call (test/exploratory)

    # Effect-size + detection driven classification. gained/lost require a ROBUST detection difference
    # (majority of donors on vs minority off) so sparse 2-vs-1 flips don't masquerade as on/off; both-
    # detected edges with a large fold-change are quantitative/rewired. FDR is a confidence annotation
    # (n is typically small), not the gate; the contrast leads with effect sizes.
    n = len(donors_common)
    strict_on = max(min_detect + 1, int(np.ceil(0.75 * n)))
    strict_off = int(np.floor(0.25 * n))
    ca, cb = "detect_%s" % state_a, "detect_%s" % state_b

    def classify(row, lfc=1.0):
        a, b = row[ca], row[cb]
        if a >= strict_on and b <= strict_off:
            return "gained_%s" % state_a
        if b >= strict_on and a <= strict_off:
            return "lost_%s" % state_a
        if a >= min_detect and b >= min_detect:
            if abs(row["log2FC"]) >= lfc:
                return "rewired" if row["driver"] == "source_isoform" else "quantitative"
            return "unchanged"
        return "state_biased_lowconf"     # detection differs but underpowered (e.g. 2-vs-1 of 4 donors)

    out["class"] = out.apply(classify, axis=1)
    out["differential"] = out["class"].isin(
        ["gained_%s" % state_a, "lost_%s" % state_a, "rewired", "quantitative"])
    out = out.reindex(out["log2FC"].abs().sort_values(ascending=False).index).reset_index(drop=True)
    logger.info("[coexpr] differential edges: %d tested | effect-size diff: %d | raw-p<%.2g: %d | FDR<0.1: %d",
                len(out), int(out["differential"].sum()), p_alpha, int(out["significant"].sum()),
                int(out["significant_FDR"].sum()))
    return out


# --------------------------------------------------------------------------- gene DE (for PDI coherence)
def gene_differential(gene_tpm, state_a, state_b, donors_common):
    """Per-gene donor-paired differential expression (log2FC + moderated p), used for PDI target
    coherence. Returns a DataFrame indexed by ENSG."""
    a = np.log2(gene_tpm.loc[[(state_a, d) for d in donors_common]].to_numpy() + EPS)
    b = np.log2(gene_tpm.loc[[(state_b, d) for d in donors_common]].to_numpy() + EPS)
    diff = (a - b).T                      # genes x donors
    mean_d, t, dfree, p = _moderated_paired(diff)
    g = pd.DataFrame({"ENSG": gene_tpm.columns, "log2FC": mean_d, "p": p}, ).set_index("ENSG")
    g["FDR"] = _bh(g["p"].to_numpy())
    g["signed_score"] = -np.log10(np.clip(g["p"], 1e-300, 1)) * np.sign(g["log2FC"])
    return g


# --------------------------------------------------------------------------- PDI regulon coherence
def pdi_coherence(edge_diff, gene_de, state_a, state_b, m1h=None):
    """For each TF isoform with PDI targets, test whether its target genes coherently respond
    (differential expression vs background) and infer activation/repression direction; combine the
    isoform's own edge differential with the regulon-response p (Stouffer) to boost confident,
    consequential PDIs. Returns a per-TF-isoform table."""
    from scipy import stats
    pdi = edge_diff[edge_diff["interaction_type"] == "PDI"].copy()
    if pdi.empty:
        return pd.DataFrame()
    bg = gene_de["log2FC"].dropna()
    abs_bg = bg.abs()
    rows = []
    for (sid, sym, enst), grp in pdi.groupby(["source_isoform_id", "Symbol", "source_ENST"]):
        tgt = [g for g in grp["target_ENSG"].unique() if g in gene_de.index]
        if len(tgt) < 3:
            continue
        d = gene_de.loc[tgt, "log2FC"].dropna()
        if len(d) < 3:
            continue
        # regulon RESPONSE (either direction): are targets more DE than background?
        u, p_resp = stats.mannwhitneyu(d.abs(), abs_bg, alternative="greater")
        # direction: net shift of targets + sign consistency
        up, down = int((d > 0).sum()), int((d < 0).sum())
        net = float(np.median(d))
        p_sign = stats.binomtest(max(up, down), up + down, 0.5).pvalue if (up + down) else 1.0
        # isoform's own activity change (median log2FC of its PDI edges) -> couple to direction
        iso_fc = float(np.nanmedian(grp["log2FC"]))
        if up + down and abs(up - down) / (up + down) >= 0.2:
            regulation = "activation" if (net > 0) == (iso_fc > 0) else "repression"
        else:
            regulation = "mixed"
        # combine the isoform's best PDI-edge p with the regulon-response p (Stouffer)
        edge_p = float(np.nanmin(grp["p"])) if grp["p"].notna().any() else 1.0
        z = np.array([stats.norm.isf(min(max(edge_p, 1e-300), 1 - 1e-16)),
                      stats.norm.isf(min(max(p_resp, 1e-300), 1 - 1e-16))])
        z_comb = z.sum() / np.sqrt(2)
        p_comb = float(stats.norm.sf(z_comb))
        m1h_call = (m1h.get(sid) or m1h.get(sym, "")) if m1h else ""
        rows.append((sid, sym, enst, len(tgt), iso_fc, net, up, down, p_resp, p_sign, regulation,
                     m1h_call, edge_p, p_comb))
    out = pd.DataFrame(rows, columns=[
        "source_isoform_id", "Symbol", "source_ENST", "n_targets", "isoform_edge_log2FC",
        "target_median_log2FC", "targets_up", "targets_down", "p_regulon_response", "p_direction",
        "inferred_regulation", "m1h_call", "p_edge", "p_combined"])
    if not out.empty:
        out["FDR_combined"] = _bh(out["p_combined"].to_numpy())
        out = out.sort_values("p_combined").reset_index(drop=True)
    logger.info("[coexpr] PDI coherence: %d TF isoforms with >=3 expressed targets | regulon FDR<0.1: %d",
                len(out), int((out.get("FDR_combined", pd.Series(dtype=float)) < 0.1).sum()))
    return out


# --------------------------------------------------------------------------- driver TF isoforms
def driver_isoforms(edge_diff, pdi_coh):
    """Rank TF isoforms by network impact: PPI rewiring (sum of -log10FDR*|log2FC| over significant
    edges) plus PDI regulatory consequence (-log10 of the combined regulon p)."""
    sig = edge_diff[edge_diff["differential"]].copy()
    sig["impact"] = -np.log10(np.clip(sig["p"], 1e-300, 1)) * sig["log2FC"].abs()
    ppi = (sig.groupby(["source_isoform_id", "Symbol"])
           .agg(n_sig_edges=("impact", "size"), ppi_impact=("impact", "sum")).reset_index())
    if pdi_coh is not None and not pdi_coh.empty:
        reg = pdi_coh.assign(reg_impact=-np.log10(np.clip(pdi_coh["p_combined"], 1e-300, 1)))[
            ["source_isoform_id", "Symbol", "reg_impact", "inferred_regulation", "n_targets"]]
        drivers = ppi.merge(reg, on=["source_isoform_id", "Symbol"], how="outer")
    else:
        drivers = ppi.assign(reg_impact=0.0, inferred_regulation="", n_targets=0)
    for c in ["n_sig_edges", "ppi_impact", "reg_impact", "n_targets"]:
        drivers[c] = pd.to_numeric(drivers[c], errors="coerce").fillna(0)
    drivers["network_impact"] = drivers["ppi_impact"] + drivers["reg_impact"]
    drivers = drivers.sort_values("network_impact", ascending=False).reset_index(drop=True)
    return drivers


def load_m1h():
    p = os.path.join(config.DATA_DIR, "activation_m1h.tsv")
    if not os.path.exists(p):
        return {}
    df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
    col = next((c for c in df.columns if "call" in c.lower() or "m1h" in c.lower()), None)
    cid = next((c for c in df.columns if "clone" in c.lower() or "isoform" in c.lower()), None)
    return dict(zip(df[cid], df[col])) if (col and cid) else {}
