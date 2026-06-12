#!/usr/bin/env python3
"""Validate satay_udon_core against the ORIGINAL pyudon satay_udon.py on synthetic data:
  1. make_mean_based_binary numerically identical to pyudon's.
  2. fishers_clinical_feats (donor_of=None == per-pseudobulk) reproduces pyudon's p-values on the
     cells both test (a >= floor).
  3. cmh_clinical_feats (batch-stratified) runs and returns finite p-values.
  4. load_metadata parses a standard table (binary + categorical + opt-in mean-binarize).
Reproducible: `python validate_satay_core.py`. No network, no real data."""
import os, sys, importlib.util, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import satay_udon_core as C

# load the ORIGINAL pyudon satay_udon.py directly (light deps: pandas/numpy/scipy/statsmodels)
_spec = importlib.util.spec_from_file_location("py_satay", "/Users/saljh8/Downloads/pyudon-main/pyudon/satay_udon.py")
PY = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(PY)

rng = np.random.RandomState(0)
cell_types = ["CT1", "CT2", "CT3"]
donors = [f"D{i:02d}" for i in range(1, 25)]
# build udon_clusters: index celltype__donor, column 'cluster'; plant signal: mut+ donors -> U1
mut = {d: int(i < 12) for i, d in enumerate(donors)}          # first 12 donors mutated
rows = {}
for ct in cell_types:
    for d in donors:
        if ct == "CT1":
            cl = "U1" if mut[d] and rng.rand() < 0.8 else rng.choice(["U2", "U3"])
        else:
            cl = rng.choice(["U1", "U2", "U3"])
        rows[f"{ct}__{d}"] = cl
udon = pd.DataFrame({"cluster": rows})
donor_meta = pd.DataFrame({"mut": [mut[d] for d in donors],
                           "batch": rng.choice(["b1", "b2"], len(donors)),
                           "sex": rng.choice(["M", "F"], len(donors)),
                           "age": rng.randint(30, 80, len(donors))}, index=donors)

ok = True

# --- 1. make_mean_based_binary ---
mine = C.make_mean_based_binary(donor_meta, ["age"]).values.ravel()
theirs = PY.make_mean_based_binary(donor_meta, ["age"]).values.ravel()
m1 = np.array_equal(mine, theirs)
print(f"[1] make_mean_based_binary identical to pyudon: {m1}"); ok &= m1

# --- 2. fishers_clinical_feats vs pyudon (per-pseudobulk == donor_of=None) ---
mine_f = C.fishers_clinical_feats(donor_meta, "mut", udon.copy(), donor_of=None, n_units=3, min_pb=1)["p_val"]
theirs_f = PY.fishers_clinical_feats(donor_meta.copy(), "mut", udon_clusters=udon.copy(), p_val=0.1, n_samples=3)["p_val"]
mine_f = mine_f.reindex(index=theirs_f.index, columns=theirs_f.columns)
both = (~mine_f.isna()) & (~theirs_f.isna())
diff = (mine_f.values[both.values] - theirs_f.values[both.values])
m2 = both.values.sum() > 0 and np.allclose(diff, 0, atol=1e-9)
print(f"[2] fishers p-values match pyudon on {int(both.values.sum())} co-tested cells (max|Δ|="
      f"{np.abs(diff).max() if diff.size else 0:.2e}): {m2}"); ok &= m2
# show the planted CT1/U1 signal
print(f"    CT1/U1 fisher p = {theirs_f.loc['CT1','U1']:.2e} (planted enrichment)")

# --- 3. cmh_clinical_feats runs + finite ---
try:
    cmh = C.cmh_clinical_feats(donor_meta, "mut", "batch", udon.copy(), donor_of=None, n_units=3, min_pb=1)
    m3 = np.isfinite(cmh.values).any()
    print(f"[3] cmh_clinical_feats runs; finite p in CT1 row: {dict(cmh.loc['CT1'].dropna().round(3))}; ok={m3}")
except Exception as e:
    m3 = False; print(f"[3] cmh_clinical_feats FAILED: {e}")
ok &= m3

# --- 4. load_metadata on a standard sample-level table ---
tmp = "/tmp/_satay_meta_test.tsv"
samp_meta = pd.DataFrame({"Sample": donors, "Donor_ID": donors,
                          "mut": [mut[d] for d in donors],
                          "sex": donor_meta["sex"].values, "age": donor_meta["age"].values})
samp_meta.to_csv(tmp, sep="\t", index=False)
try:
    C.load_metadata(tmp)                       # should ERROR: 'age' is continuous, not in mean_binarize
    m4a = False
except ValueError as e:
    m4a = "age" in str(e); print(f"[4a] load_metadata strict-errors on raw numeric 'age': {m4a}")
binary, d_of = C.load_metadata(tmp, mean_binarize=["age"])
m4b = ("mut" in binary.columns and any(c.startswith("sex=") for c in binary.columns)
       and "age" in binary.columns and d_of is not None)
print(f"[4b] load_metadata(mean_binarize=['age']) -> covariates {list(binary.columns)}; donor map set={d_of is not None}: {m4b}")
ok &= (m4a and m4b)

print("\n=== VALIDATION", "PASSED" if ok else "FAILED", "===")
sys.exit(0 if ok else 1)
