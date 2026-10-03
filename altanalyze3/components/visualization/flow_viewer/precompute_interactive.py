#!/usr/bin/env python3
"""Assemble the interactive rna2flow bundle: embeddings, labels, channels, gates.

Everything the browser needs is written as compact binary (float32 / int16 codes) so a
100,000-event panel loads without a server round trip per interaction.
"""
import argparse, glob, json, os, re, sys
import numpy as np, pandas as pd

sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
sys.path.insert(0, "/Users/saljh8/Dropbox/Code/pyInfinityFlow-main")
from altanalyze3.components.rna2flow.io import read_fcs

FLOW_ROOT = ("/Users/saljh8/Dropbox/Collaborations/Grimes/Thymus/FlowData/full_dataset/"
             "UMAP and FlowSOM Anlaysis")
FCS = os.path.join(FLOW_ROOT, "Concat_MultiLin_Thymus_8_24_26.fcs")
SKIP = {"SSC-W", "SSC-H", "FSC-W", "FSC-H", "SSC-B-W", "SSC-B-H", "SSC-B-A",
        "Time", "AF-A", "SampleID", "TUBENAME"}


def codes(series):
    """Categorical -> (int16 codes, level list). int16 keeps the payload small."""
    c = pd.Categorical(pd.Series(series).astype(str))
    return c.codes.astype(np.int16), list(map(str, c.categories))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    ap.add_argument("--transfer-dir", action="append", default=[],
                    help="directory holding pred_*.npy from a benchmark run")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    arr = os.path.join(a.out, "arrays"); os.makedirs(arr, exist_ok=True)

    from altanalyze3.components.rna2flow.transform import read_fcs_anndata, logicle_anndata
    ad_fcs = read_fcs_anndata(FCS); logicle_anndata(ad_fcs, in_place=True)
    raw = read_fcs(FCS)
    n = ad_fcs.shape[0]
    X = np.asarray(ad_fcs.X, dtype=np.float32)
    chan = [(det, ab) for det, ab in zip(map(str, ad_fcs.var_names), raw.antibodies)]
    keep = [i for i, (d, ab) in enumerate(chan) if ab not in SKIP]
    channels = [chan[i][1] for i in keep]
    X = X[:, keep]
    X.tofile(os.path.join(arr, "channels.f32"))
    display_transform = {channels[j]: {'detector': str(ad_fcs.var.index[i]),
        **{k: float(ad_fcs.var.iloc[i][k]) for k in ['LOGICLE_T','LOGICLE_W','LOGICLE_M','LOGICLE_A']},
        'USE_LOGICLE': bool(ad_fcs.var.iloc[i]['USE_LOGICLE']),
        'IMPUTED': bool(ad_fcs.var.iloc[i]['IMPUTED'])} for j,i in enumerate(keep)}
    print("channels: %d x %d events" % (len(channels), n))

    # ---- embeddings: every flow UMAP run, joined on EventNumberDP ----
    embeddings = {}
    for f in sorted(glob.glob(os.path.join(FLOW_ROOT, "**", "UMAP_Results_*.csv"), recursive=True)):
        rid = re.search(r"UMAP_Results_(\w+)\.csv", f).group(1)
        d = pd.read_csv(f)
        if len(d) != n:
            print("  skip UMAP %s: %d rows vs %d events" % (rid, len(d), n)); continue
        if d.iloc[:,0].duplicated().any() or set(d.iloc[:,0])!=set(range(1,n+1)):
            raise ValueError('UMAP event identifiers do not match FCS: '+f)
        d=d.set_index(d.columns[0]).reindex(range(1,n+1))
        xy = d.iloc[:, :2].to_numpy(np.float32)
        xy.tofile(os.path.join(arr, "emb_flow_%s.f32" % rid))
        embeddings["flow_UMAP_%s" % rid] = {"file": "emb_flow_%s.f32" % rid, "n": n,
                                            "space": "flow"}
        print("  embedding flow_UMAP_%s" % rid)

    # ---- labels on flow events: every FlowSOM run + every transferred label set ----
    labels = {}
    for f in sorted(glob.glob(os.path.join(FLOW_ROOT, "**", "FlowSOM_Results_*.csv"), recursive=True)):
        rid = re.search(r"FlowSOM_Results_(\w+)\.csv", f).group(1)
        d = pd.read_csv(f)
        if len(d) != n: continue
        if d.iloc[:,0].duplicated().any() or set(d.iloc[:,0])!=set(range(1,n+1)):
            raise ValueError('FlowSOM event identifiers do not match FCS: '+f)
        d=d.set_index(d.columns[0]).reindex(range(1,n+1))
        c, lev = codes(d.iloc[:, 0])
        c.tofile(os.path.join(arr, "lab_FlowSOM_%s.i16" % rid))
        labels["FlowSOM_%s" % rid] = {"file": "lab_FlowSOM_%s.i16" % rid, "levels": lev,
                                      "space": "flow", "source": "FlowSOM"}
        print("  labels FlowSOM_%s (%d clusters)" % (rid, len(lev)))

    for td in a.transfer_dir:
        for f in sorted(glob.glob(os.path.join(td, "pred_*.npy"))):
            name = os.path.basename(f)[5:-4]
            v = np.load(f, allow_pickle=True)
            if len(v) != n: continue
            c, lev = codes(v)
            key = "transfer_%s" % name
            c.tofile(os.path.join(arr, "lab_%s.i16" % key))
            labels[key] = {"file": "lab_%s.i16" % key, "levels": lev, "space": "flow",
                           "source": "transferred", "origin": os.path.basename(td.rstrip("/"))}
    print("  label sets: %d" % len(labels))

    json.dump({"n_events": int(n), "channels": channels,
               "embeddings": embeddings, "labels": labels,
               "channels_file": "channels.f32", 'display_transform': display_transform,
               'provenance': {'fcs': FCS, 'transform': 'per-channel pyInfinityFlow Logicle'}},
              open(os.path.join(a.out, "manifest.json"), "w"), indent=1)
    print("bundle -> %s" % a.out)


if __name__ == "__main__":
    main()
