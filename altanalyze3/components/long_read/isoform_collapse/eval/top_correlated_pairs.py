#!/usr/bin/env python3
"""Export the most-correlated isoform<->junction pairs (best within-structure junction per isoform,
Spearman, per cell) for BOTH WTA and EM, side by side, to a TSV. Reuses the per-method best_corr_*.tsv
already written by em_vs_wta_benchmark.py (isoform, gene, best_junction_corr, best_junction, iso_reads)."""
import pandas as pd, os
B="/Users/saljh8/Dropbox/Revio/test/Isoform-Workflow-BAM"; D=f"{B}/analysis/em_vs_wta"
ANN="/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart100/ensembl/Hs/Hs_Ensembl-annotations.txt"
sym={}
for line in open(ANN):
    t=line.rstrip("\n").split("\t")
    if len(t)>=2: sym[t[0]]=t[1]

w=pd.read_csv(f"{D}/best_corr_wta.tsv",sep="\t")
e=pd.read_csv(f"{D}/best_corr_em.tsv",sep="\t")
w["gene_symbol"]=w["gene"].map(sym); e["gene_symbol"]=e["gene"].map(sym)

# side-by-side on shared isoforms
m=w.merge(e,on=["isoform","gene"],suffixes=("_wta","_em"))
m["gene_symbol"]=m["gene"].map(sym)
m=m[["isoform","gene","gene_symbol","iso_reads_wta",
     "best_junction_wta","best_junction_corr_wta",
     "best_junction_em","best_junction_corr_em"]]
m["max_corr"]=m[["best_junction_corr_wta","best_junction_corr_em"]].max(axis=1)
m=m.sort_values("max_corr",ascending=False)
m.to_csv(f"{D}/top_correlated_pairs_wta_vs_em.tsv",sep="\t",index=False)
print(f"wrote {D}/top_correlated_pairs_wta_vs_em.tsv  ({len(m):,} shared isoform-junction pairs)")

# also each method's own ranked best pairs
for lab,df in (("wta",w),("em",e)):
    d=df.sort_values("best_junction_corr",ascending=False)[
        ["isoform","gene","gene_symbol","best_junction","best_junction_corr","iso_reads"]]
    d.to_csv(f"{D}/top_correlated_pairs_{lab}.tsv",sep="\t",index=False)
    print(f"wrote {D}/top_correlated_pairs_{lab}.tsv  ({len(d):,})")
print("\n=== top 10 by max corr (shared) ===")
print(m.head(10).to_string(index=False))
print("DONE")
