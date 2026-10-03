#!/usr/bin/env python3
"""Add CITE-seq spaces to a flow_viewer bundle.

A bundle holds several SPACES. Each space has its own event/cell count, its own embeddings,
its own label sets and its own feature matrix:

  flow   100,000 FCS events, logicle channels, FlowSOM runs, transferred annotations
  cite_* CITE-seq cells, ADT features, author / bone-marrow / scTriangulate annotations

Keeping them as separate spaces is deliberate: flow events and CITE cells are different
objects, and pretending one index serves both is how a viewer starts lying.
"""
import argparse, json, os, sys
import numpy as np, pandas as pd, anndata as ad

sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")


def codes(v):
    c = pd.Categorical(pd.Series(v).astype(str))
    return c.codes.astype(np.int16), list(map(str, c.categories))


def add_space(man, arr, name, X, features, embeddings, label_cols, prefix):
    n = X.shape[0]
    fX = os.path.join(arr, "%s_features.f32" % prefix)
    np.asarray(X, dtype=np.float32).tofile(fX)
    sp = {"n": int(n), "features": list(features),
          "features_file": os.path.basename(fX), "embeddings": {}, "labels": {}}
    for ename, xy in embeddings.items():
        f = "%s_emb_%s.f32" % (prefix, ename)
        np.asarray(xy, dtype=np.float32).tofile(os.path.join(arr, f))
        sp["embeddings"][ename] = {"file": f}
    for lname, vals in label_cols.items():
        c, lev = codes(vals)
        f = "%s_lab_%s.i16" % (prefix, lname)
        c.tofile(os.path.join(arr, f))
        sp["labels"][lname] = {"file": f, "levels": lev}
    man["spaces"][name] = sp
    print("  space %-14s %6d cells | %2d features | %d embeddings | %d label sets"
          % (name, n, len(features), len(sp["embeddings"]), len(sp["labels"])))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--sct-h5ad", required=True)
    ap.add_argument("--adt-h5ad", required=True)
    ap.add_argument("--rna-umap", required=True)
    ap.add_argument("--adt-umap", required=True)
    ap.add_argument("--prefix", default="chinese")
    a = ap.parse_args()
    root = os.path.abspath(a.bundle); arr = os.path.join(root, "arrays")
    man = json.load(open(os.path.join(root, "manifest.json")))

    # Fold the existing flat flow layout into the new spaces layout, once.
    if "spaces" not in man:
        man["spaces"] = {"flow": {**{k:man[k] for k in ("display_transform","provenance") if k in man}, "n": man["n_events"], "features": man["channels"],
                                  "features_file": man["channels_file"],
                                  "embeddings": man["embeddings"], "labels": man["labels"]}}
        print("  space flow           %6d events | %2d channels | %d embeddings | %d label sets"
              % (man["n_events"], len(man["channels"]), len(man["embeddings"]), len(man["labels"])))

    s = ad.read_h5ad(a.sct_h5ad)
    emb = {"scTriangulate_UMAP": np.asarray(s.obsm["X_umap"], np.float32)}
    labs = {c: s.obs[c].astype(str).values for c in
            ["pruned", "Author_celltype", "Mm-MarrowAtlas-L4", "Leiden_RNA_r2",
             "Leiden_ADT_r2", "ICGS3_RNA"] if c in s.obs.columns}
    ab = [v for v in s.var_names if str(v).startswith("AB_")]
    Xab = np.asarray(s[:, ab].X.todense() if hasattr(s[:, ab].X, "todense") else s[:, ab].X,
                     dtype=np.float32)
    add_space(man, arr, "cite_scTriangulate_%s" % a.prefix, Xab,
              [v[3:] for v in ab], emb, labs, "%s_sct" % a.prefix)

    adt = ad.read_h5ad(a.adt_h5ad)
    rna_xy = pd.read_csv(a.rna_umap, sep="\t", index_col=0)
    adt_xy = pd.read_csv(a.adt_umap, sep="\t", index_col=0)
    common = adt.obs_names.intersection(rna_xy.index).intersection(adt_xy.index)
    sub = adt[common]
    X = np.asarray(sub.X.todense() if hasattr(sub.X, "todense") else sub.X, dtype=np.float32)
    add_space(man, arr, "cite_markerUMAP_%s" % a.prefix, X, list(sub.var_names),
              {"RNA_marker_UMAP": rna_xy.loc[common, ["UMAP-X", "UMAP-Y"]].to_numpy(np.float32),
               "ADT_marker_UMAP": adt_xy.loc[common, ["UMAP-X", "UMAP-Y"]].to_numpy(np.float32)},
              {c: sub.obs[c].astype(str).values for c in ["celltype", "library"]
               if c in sub.obs.columns},
              "%s_mk" % a.prefix)

    json.dump(man, open(os.path.join(root, "manifest.json"), "w"), indent=1)
    print("spaces now: %s" % ", ".join(man["spaces"]))


if __name__ == "__main__":
    main()
