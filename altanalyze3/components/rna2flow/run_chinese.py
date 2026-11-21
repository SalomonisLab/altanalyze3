#!/usr/bin/env python3
"""Run every CITE-seq -> flow transfer method on the Chinese thymus data and score it."""
import json, os, sys, time
import numpy as np, pandas as pd, anndata as ad

sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
sys.path.insert(0, "/Users/saljh8/Dropbox/Code/pyInfinityFlow-main")
from altanalyze3.components.rna2flow.io import read_flowjo_rds
from altanalyze3.components.rna2flow.crosswalk import build_crosswalk
from altanalyze3.components.rna2flow.transfer import METHODS
from altanalyze3.components.rna2flow.evaluate import score_transfer, adversarial_controls
from altanalyze3.components.rna2flow.marker_concordance import concordance

CN = "/Users/saljh8/Dropbox/Collaborations/Grimes/Thymus/Chinese"
FLOW_RDS = "/Users/saljh8/Dropbox/Collaborations/Grimes/Thymus/FlowData/Ungated_J8DW/FlowSOM_RData_J8DW.RDS"
FLOW_LAB = "/Users/saljh8/Dropbox/Collaborations/Grimes/Thymus/FlowData/Ungated_J8DW/FlowSOM_Results_J8DW.csv"
OUT = f"{CN}/rna2flow"
os.makedirs(OUT, exist_ok=True)

def log(m): print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)

log("reading flow")
fl = read_flowjo_rds(FLOW_RDS)
lab = pd.read_csv(FLOW_LAB)
flowsom = lab.iloc[:, 1].astype(str).values
assert len(flowsom) == fl.X.shape[0], (len(flowsom), fl.X.shape)
log("flow: %d events x %d channels, %d FlowSOM metaclusters"
    % (fl.X.shape[0], fl.X.shape[1], len(set(flowsom))))

log("reading CITE-seq ADT + labels")
adt = ad.read_h5ad(f"{CN}/PRJCA039526_thymus_ADT_dsb_annotated.h5ad")
sct = ad.read_h5ad(f"{CN}/scTriangulate/run/sctriangulate.h5ad", backed="r")
sct_lab = pd.Series(sct.obs["pruned"].astype(str).values, index=sct.obs_names)
ch = ad.read_h5ad(f"{CN}/combined_with_umap_and_markers.h5ad", backed="r")
bm = pd.Series(ch.obs["Mm-MarrowAtlas-L4"].astype(str).values,
               index=[f"{n.rsplit('.',1)[1]}_{n.rsplit('.',1)[0]}" for n in ch.obs_names])

label_sets = {"Author_celltype": adt.obs["celltype"].astype(str).values}
label_sets["Ferchen_L4_BM"] = bm.reindex(adt.obs_names).fillna("unassigned").values
sct_idx = pd.Index([f"{n.rsplit('.',1)[1]}_{n.rsplit('.',1)[0]}" for n in sct.obs_names])
label_sets["scTriangulate"] = pd.Series(sct_lab.values, index=sct_idx).reindex(adt.obs_names).fillna("unassigned").values
for k, v in label_sets.items():
    log("  %-18s %d labels, %d unassigned" % (k, len(set(v)), int((v == "unassigned").sum())))

log("crosswalk")
pairs, unf, _ = build_crosswalk(fl.antibodies, list(adt.var_names), extra_aliases={"CD8": "CD8a"})
pairs = pairs.drop_duplicates(subset=["cite"], keep="first").reset_index(drop=True)
log("  %d shared markers; unmatched flow: %s" % (len(pairs), ", ".join(unf)))
pairs.to_csv(f"{OUT}/crosswalk_chinese.tsv", sep="\t", index=False)

fi = {c: i for i, c in enumerate(fl.channels)}
ci = {c: i for i, c in enumerate(adt.var_names)}
flow_cols = []
for f in pairs["flow"]:
    hit = [c for c in fl.channels if (c.split("::", 1)[1].strip() if "::" in c else c.strip()) == f]
    flow_cols.append(fi[hit[0]])
Fx = np.asarray(fl.X[:, flow_cols], dtype=np.float64)
Cx = np.asarray(adt.X, dtype=np.float64)[:, [ci[c] for c in pairs["cite"]]]
log("  matrices: CITE %s  FLOW %s" % (Cx.shape, Fx.shape))

rows = []
for lname, labels in label_sets.items():
    keep = labels != "unassigned"
    C, L = Cx[keep], labels[keep]
    log("== %s: %d CITE cells, %d labels" % (lname, keep.sum(), len(set(L))))
    for mname, fn in METHODS.items():
        t = time.time()
        try:
            pred = fn(C, Fx, L)
            s = score_transfer(pred, flowsom)
            s.update(adversarial_controls(fn, C, Fx, L, flowsom))
            cdf, csum = concordance(C, L, Fx, pred, pairs)
            s.update(csum)
            cdf.to_csv(f"{OUT}/concordance_{lname}_{mname}.tsv", sep="\t", index=False)
            s.update(label_set=lname, method=mname, seconds=round(time.time() - t, 1),
                     n_markers=len(pairs), n_cite=int(keep.sum()))
            rows.append(s)
            log("   %-16s ARI=%.4f | MarkerFinder pass %d/%d rho=%.3f %s | top=%.2f shuffled=%+.3f %.0fs"
                % (mname, s["ARI"], s["populations_passing"], s["populations_compared"],
                   s.get("median_rank_spearman", float("nan")),
                   "MEETS TARGET" if s.get("meets_target") else "", s["largest_label_fraction"],
                   s["shuffled_labels_ARI"], s["seconds"]))
            np.save(f"{OUT}/pred_{lname}_{mname}.npy", pred.astype(str))
        except Exception as e:
            log("   %-16s FAILED: %s" % (mname, str(e)[:160]))
            rows.append(dict(label_set=lname, method=mname, error=str(e)[:300]))
        # Write after EVERY method. A hard crash in one must not cost the rows already earned.
        pd.DataFrame(rows).to_csv(f"{OUT}/transfer_benchmark.tsv", sep="\t", index=False)

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/transfer_benchmark.tsv", sep="\t", index=False)
log("wrote %s/transfer_benchmark.tsv" % OUT)
