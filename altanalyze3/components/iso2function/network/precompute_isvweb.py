#!/usr/bin/env python3
"""Pre-compute the compact, queryable artifact the ISV-web network tab loads at startup.

Produces ONE JSON (`isvweb_network.json`) in the dataset/run dir:
  - edges:    isoform-resolved interactions [src_iso, src_enst, src_gene, tgt_gene, tgt_ensg, type(PDI/PPI), score]
  - iso_meta: per source isoform -> {gene, structure(E/I), enst}
  - iso_cpm:  source-isoform CPM  : ENST  -> {"<state>||<group>": cpm}   (only nonzero stored)
  - gene_cpm: target-gene max-isoform CPM : ENSG -> {"<state>||<group>": cpm}
  - contexts: sorted list of "<state>||<group>"
Expression is the per-(cell-state x covariate group) donor-mean from the AltAnalyze3 long-read pseudobulks
(reusing iso2function.network.coexpr). Meant to run as a precompute default in the LR workflow.

Usage: python -m altanalyze3.components.iso2function.network.precompute_isvweb --dataset_dir DIR [--out DIR]
"""
import os, json, argparse, logging
import numpy as np, pandas as pd
from . import coexpr

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger(__name__)


def _load_protein_meta(dataset_dir):
    """gff-output/protein_summary.txt -> {versionless transcript id: (protein_length_aa, nmd_status)}.
    Keyed by ENST (reference) and bare molecule id (novel); used to annotate each source isoform with its
    amino-acid length and to drop NMD isoforms (no stable protein -> no PPI/PDI) from the network."""
    path = os.path.join(dataset_dir, "gff-output", "protein_summary.txt")
    if not os.path.exists(path):
        log.info("[precompute] no protein_summary.txt -> no AA/NMD annotation"); return {}
    pm = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            base = p[1].strip().split(".")[0]
            try:
                aa = int(p[2])
            except (ValueError, IndexError):
                aa = None
            pm.setdefault(base, (aa, p[3].strip()))
    log.info("[precompute] protein_summary: %d transcripts (AA + NMD status)", len(pm))
    return pm


def _gene_to_isoforms(dataset_dir):
    """ENSG -> [base ENST, ...] from the isoform h5ad var names ('ENSG:ENST'). Backed read (var names only,
    not the matrix). Used to list the known isoforms underlying each target/partner gene's expression."""
    import anndata as ad
    gene2ensts = {}
    try:
        iso = ad.read_h5ad(os.path.join(dataset_dir, coexpr.ISO_H5AD), backed="r")
    except Exception as e:
        log.info("[precompute] could not read isoform var names for gene-isoform map: %s", e); return gene2ensts
    for v in iso.var_names:
        s = str(v)
        if ":" not in s:
            continue
        ensg = s.split(":")[0].split(".")[0]; tail = s.split(":", 1)[1]
        if tail.startswith("ENST"):
            gene2ensts.setdefault(ensg, []).append(tail.split(".")[0])
    try:
        iso.file.close()
    except Exception:
        pass
    return gene2ensts


def _group_donors(dataset_dir, meta):
    m = pd.read_csv(os.path.join(dataset_dir, meta), sep="\t", dtype=str).fillna("")
    groups = [g for g in sorted(set(m["groups"])) if g]
    gd = {}
    for g in groups:
        sel = m[m["groups"] == g]
        gd[g] = set(sel["library"]) | (set(sel["uid"]) if "uid" in sel.columns else set())
    return groups, gd


def precompute(dataset_dir, out_dir=None, meta=None):
    meta = meta or coexpr.META
    out_dir = out_dir or dataset_dir
    edges = coexpr.load_edges()
    sym2ensg = coexpr.load_symbol_to_ensg(edges)
    edges = edges.copy()
    edges["target_ENSG"] = edges["target"].map(sym2ensg).fillna("")
    keep = edges[(edges["best_ENST"] != "") & (edges["target_ENSG"] != "")].drop_duplicates(
        ["source_isoform_id", "best_ENST", "target", "target_ENSG", "interaction_type"])
    src_ensts = set(keep["best_ENST"]); tgt_ensgs = set(keep["target_ENSG"])
    log.info("[precompute] %d edges (dedup), %d source isoforms, %d target genes", len(keep), len(src_ensts), len(tgt_ensgs))

    # target-gene -> its known ENST isoforms (the isoforms underlying each partner/target gene's expression)
    gene2ensts = _gene_to_isoforms(dataset_dir)
    tgt_iso_ensts = set(e for g in tgt_ensgs for e in gene2ensts.get(g, []))
    log.info("[precompute] %d target-gene isoforms (known ENSTs) to break down", len(tgt_iso_ensts))

    groups, gd = _group_donors(dataset_dir, meta)
    log.info("[precompute] covariate groups: %s", groups)
    iso_cpm, gene_cpm, tgt_iso_cpm, contexts = {}, {}, {}, set()
    for g in groups:
        iso_tpm, _, _ = coexpr.load_groups(dataset_dir, None, keep_donors=gd[g])
        mx = coexpr.load_max_isoform_expr(dataset_dir, None, keep_donors=gd[g])
        im = iso_tpm.groupby(level="state").mean(); gm = mx.groupby(level="state").mean()
        scols = [e for e in src_ensts if e in im.columns]; gcols = [e for e in tgt_ensgs if e in gm.columns]
        ticols = [e for e in tgt_iso_ensts if e in im.columns]
        for st in im.index:
            ctx = f"{st}||{g}"; contexts.add(ctx)
            s = im.loc[st, scols]
            for ens, v in s[s > 0].items(): iso_cpm.setdefault(ens, {})[ctx] = round(float(v), 3)
            ts = im.loc[st, ticols]
            for ens, v in ts[ts > 0].items(): tgt_iso_cpm.setdefault(ens, {})[ctx] = round(float(v), 3)
            if st in gm.index:
                gs = gm.loc[st, gcols]
                for ens, v in gs[gs > 0].items(): gene_cpm.setdefault(ens, {})[ctx] = round(float(v), 3)
        log.info("[precompute]   group %s: %d states", g, im.shape[0])

    prot = _load_protein_meta(dataset_dir)
    # gene -> [ENST...] for isoforms actually expressed somewhere; + AA length per target isoform
    gene_isoforms = {}
    for g in tgt_ensgs:
        es = [e for e in gene2ensts.get(g, []) if e in tgt_iso_cpm]
        if es:
            gene_isoforms[g] = es
    tgt_iso_aa = {e: prot.get(e.split(".")[0], (None, ""))[0] for e in tgt_iso_cpm}
    log.info("[precompute] gene-isoform breakdown: %d genes, %d isoforms w/ CPM", len(gene_isoforms), len(tgt_iso_cpm))
    iso_meta = {}
    for r in keep.itertuples():
        if r.source_isoform_id not in iso_meta:
            aa, nmd = prot.get((r.best_ENST or "").split(".")[0], (None, ""))
            iso_meta[r.source_isoform_id] = {"gene": r.Symbol, "enst": r.best_ENST,
                                             "structure": r.isoform_structure, "aa": aa, "nmd": nmd}
    n_nmd = sum(1 for v in iso_meta.values() if v["nmd"] in ("NMD", "Potential-NMD"))
    log.info("[precompute] AA-annotated %d/%d source isoforms; %d flagged NMD/Potential-NMD (excluded at query)",
             sum(1 for v in iso_meta.values() if v["aa"] is not None), len(iso_meta), n_nmd)
    out_edges = [{"src_iso": r.source_isoform_id, "src_enst": r.best_ENST, "src_gene": r.Symbol,
                  "tgt_gene": r.target, "tgt_ensg": r.target_ENSG, "type": r.interaction_type,
                  "score": (round(float(r.activity_score), 4) if pd.notna(r.activity_score) else None)}
                 for r in keep.itertuples()]
    obj = {"edges": out_edges, "iso_meta": iso_meta, "iso_cpm": iso_cpm, "gene_cpm": gene_cpm,
           "gene_isoforms": gene_isoforms, "tgt_iso_cpm": tgt_iso_cpm, "tgt_iso_aa": tgt_iso_aa,
           "contexts": sorted(contexts), "groups": groups}
    out_path = os.path.join(out_dir, "isvweb_network.json")
    with open(out_path, "w") as fh: json.dump(obj, fh)
    log.info("[precompute] wrote %s  (%d edges, %d iso nodes w/ CPM, %d gene nodes w/ CPM, %d contexts)",
             out_path, len(out_edges), len(iso_cpm), len(gene_cpm), len(contexts))
    return out_path


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset_dir", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--meta", default=None)
    a = ap.parse_args()
    precompute(a.dataset_dir, a.out, a.meta)
