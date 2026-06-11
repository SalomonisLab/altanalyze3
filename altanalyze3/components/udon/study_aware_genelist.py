#!/usr/bin/env python3
"""Conserved vs unique UDON clusters by DIRECT top-50 marker GENE-LIST overlap
across studies (no GO/ontology), with a SHARED-CORE FLOOR so conserved programs
are crisp (not transitive chains):

  * pairwise overlap of top-50 gene lists (+ hypergeometric significance);
  * COMPLETE-linkage grouping on Jaccard distance (every member mutually similar,
    no transitive hubs);
  * a group is a CONSERVED program only if it spans >= MIN_STUDIES studies AND its
    members share a COMMON CORE of >= CORE_FLOOR genes (full intersection);
  * all other clusters are UNIQUE (study/subtype-specific).

Writes inspectable intermediates:
  conserved_programs.tsv   (program, #clusters, #studies, members, CORE genes)
  unique_clusters.tsv      (study-specific clusters + their top markers)
  conserved_vs_unique.tsv  (per-cluster call; consumed by study_aware_integrate.py)
  genelist_jaccard.tsv     (pairwise top-50 Jaccard)
"""
import os, sys, numpy as np, pandas as pd
from scipy.stats import hypergeom
from sklearn.cluster import AgglomerativeClustering
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pseudobulk_protocol as P

_CELLTYPE, _GENEFILT, _SPECIES, _SFX = P.udon_restriction()   # honours --cell-type via UDON_CELL_TYPE
SA = "/Users/saljh8/Dropbox/Collaborations/Grimes/UDON/cellHarmony-datasets/final/pseudobulk/UDON/study_aware" + _SFX
N_UNIVERSE = 18451
JACCARD_LINK = 0.08      # complete-linkage: every pair in a group shares Jaccard >= this (~7/50 genes)
CORE_FLOOR = 3           # a conserved program's members must share >= this many genes (full intersection)
MIN_STUDIES = 2


def main():
    m = pd.read_csv(os.path.join(SA, "per_cluster_markers.tsv"), sep="\t")
    keys = m["cluster_key"].tolist(); study = np.asarray(m["study"].astype(str))
    sets = [set(str(x).split(",")) for x in m["markers"]]
    n = len(keys)
    print(f"{n} per-study clusters; complete-linkage on top-50 gene-list Jaccard "
          f"(link>={JACCARD_LINK}); conserved needs >={MIN_STUDIES} studies AND >={CORE_FLOOR} shared core genes")

    Jm = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            k = len(sets[i] & sets[j]); u = len(sets[i] | sets[j])
            Jm[i, j] = Jm[j, i] = k / u if u else 0.0
    pd.DataFrame(Jm, index=keys, columns=keys).to_csv(os.path.join(SA, "genelist_jaccard.tsv"), sep="\t")

    D = 1.0 - Jm
    lab = AgglomerativeClustering(n_clusters=None, distance_threshold=1.0 - JACCARD_LINK,
                                  metric="precomputed", linkage="complete").fit_predict(D)

    # classify each group: conserved if >=MIN_STUDIES studies AND full-intersection core >= CORE_FLOOR
    group_of = np.full(n, -1)              # final conserved program id (-1 = unique)
    prog_rows = []; pid = 0
    for g in np.unique(lab):
        idx = np.where(lab == g)[0]
        studies = set(study[idx])
        core = set.intersection(*[sets[i] for i in idx]) if len(idx) > 1 else set()
        if len(studies) >= MIN_STUDIES and len(core) >= CORE_FLOOR:
            for i in idx:
                group_of[i] = pid
            prog_rows.append({"program": f"C{pid}", "n_clusters": len(idx), "n_studies": len(studies),
                              "studies": ",".join(sorted(studies)),
                              "members": ",".join(keys[i] for i in idx),
                              "n_core_genes": len(core), "core_genes": ",".join(sorted(core))})
            pid += 1

    status = ["conserved" if group_of[i] >= 0 else "unique" for i in range(n)]
    cvu = pd.DataFrame({"cluster_key": keys, "study": study, "genelist_group": group_of,
                        "n_studies_in_group": [len(set(study[group_of == group_of[i]])) if group_of[i] >= 0 else 1
                                               for i in range(n)],
                        "status": status, "n_markers": [len(s) for s in sets]})
    cvu.to_csv(os.path.join(SA, "conserved_vs_unique.tsv"), sep="\t", index=False)
    prog = pd.DataFrame(prog_rows)
    prog.to_csv(os.path.join(SA, "conserved_programs.tsv"), sep="\t", index=False)

    uniq = cvu[cvu["status"] == "unique"][["cluster_key", "study"]].copy()
    uniq["top_markers"] = [",".join(sorted(sets[keys.index(k)])[:15]) for k in uniq["cluster_key"]]
    uniq.to_csv(os.path.join(SA, "unique_clusters.tsv"), sep="\t", index=False)

    ncons = status.count("conserved"); nuniq = status.count("unique")
    print(f"\n=> {pid} CONSERVED programs ({ncons} clusters); {nuniq} UNIQUE clusters")
    print("\n--- CONSERVED programs (intermediate result, with shared core signature) ---")
    for r in prog_rows:
        print(f"  {r['program']}: {r['n_studies']} studies, {r['n_clusters']} clusters | "
              f"core({r['n_core_genes']})={r['core_genes'][:80]}")
        print(f"       members: {r['members']}")
    print(f"\nwrote: conserved_programs.tsv, unique_clusters.tsv, conserved_vs_unique.tsv, genelist_jaccard.tsv")
    print(f"UNIQUE clusters per study:")
    print(uniq["study"].value_counts().to_string())


if __name__ == "__main__":
    main()
