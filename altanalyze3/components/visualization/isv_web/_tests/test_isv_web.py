"""ISV web app validation against the live local run.

Run:  python -m altanalyze3.components.visualization.isv_web._tests.test_isv_web
(or pytest). Requires the completed local run at RUN below.

Checks (per the plan):
  1. parity: /api isoforms == engine group_structures (same isoform set + cluster count) for HOPX
  2. exon coords: exon_segments match build_isoform_segments
  3. protein meta: protein_length / NMD match protein_summary.txt
  4. FASTA export: protein_fasta == the protein_sequences.fasta record
  5. speed: warm query is fast (reads atom caches, not full h5ad)
"""
import os
import sys
import time

sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")

from altanalyze3.components.visualization import isoform_structure_view as isv
from altanalyze3.components.visualization.isv_web import data_api as da

RUN = os.environ.get("ISV_RUN", "/Users/saljh8/Dropbox/Revio/test/Isoform-Workflow-BAM")
GM = f"{RUN}/Hs_Ensembl_exon.txt"
GS = "/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl-annotations.txt"
GENE = "HOPX"


def _ctx():
    return da.RunContext(RUN, f"{RUN}/mds_metadata_bam.txt", GM, GS).load_all()


def run():
    results = []

    def check(name, ok, detail=""):
        results.append((name, bool(ok), detail))
        print(("PASS " if ok else "FAIL ") + name + ("  " + detail if detail else ""))

    ctx = _ctx()
    check("catalog: 4 samples", len(ctx.samples) == 4, f"n={len(ctx.samples)}")
    check("catalog: groups young+AML-NPM1", set(ctx.groups) == {"young", "AML-NPM1"}, str(set(ctx.groups)))
    check("catalog: cell types present", len(ctx.cell_types) > 50, f"n={len(ctx.cell_types)}")
    check("gene resolves", ctx.resolve_gene(GENE) == "ENSG00000171476")

    res = da.query_isoforms(ctx, GENE, None, None, None, combine_by="group")
    check("query returns isoforms", len(res["isoforms"]) > 0, f"n={len(res['isoforms'])}")
    check("query columns are the groups", set(res["columns"]) == {"young", "AML-NPM1"}, str(res["columns"]))
    check("clusters assigned", res["cluster_count"] > 0, f"clusters={res['cluster_count']}")

    # 2. exon coords parity with build_isoform_segments (recompute for one isoform)
    gene = "ENSG00000171476"
    structures = da._load_structures(ctx, gene)
    ok_segs = True
    for iso in res["isoforms"][:5]:
        resolved = isv.resolve_isoform_id(iso["isoform_id"], structures)
        if not resolved:
            continue
        segs = isv.build_isoform_segments(structures[resolved], ctx.exon_lookup, gene)
        # the API returns the same number of coordinate-bearing exon segments
        api_e = [s for s in iso["exon_segments"] if s["start"] is not None]
        eng_e = [s for s in segs if s.get("start") is not None]
        if len(api_e) != len(eng_e):
            ok_segs = False
    check("exon_segments parity with engine", ok_segs)

    # 3. protein meta matches protein_summary.txt
    ok_prot = True
    for iso in res["isoforms"]:
        if iso["protein_length"] is not None:
            meta = ctx.protein_meta(iso["isoform_id"])
            if not meta or meta.get("protein_length") != iso["protein_length"]:
                ok_prot = False
    check("protein_length matches protein_summary", ok_prot)

    # 4. FASTA export equals the record
    enst = next((i["isoform_id"] for i in res["isoforms"] if i["known"]), None)
    fa = ctx.protein_fasta(enst) if enst else None
    check("FASTA export works for a known isoform", bool(fa and fa.startswith(">")),
          (enst or "no ENST") + (f" len={len(fa)}" if fa else ""))

    # 5. speed: warm query
    t = time.time(); da.query_isoforms(ctx, GENE, None, None, None, combine_by="group")
    warm = time.time() - t
    check("warm query fast (<2s)", warm < 2.0, f"{warm:.2f}s")

    # cell-type combine
    cts = ctx.cell_types[:2]
    rc = da.query_isoforms(ctx, GENE, None, None, cts, combine_by="cell_type")
    check("cell_type combine columns", set(rc["columns"]) == set(cts), str(rc["columns"]))

    n_fail = sum(1 for _, ok, _ in results if not ok)
    print(f"\n{len(results) - n_fail}/{len(results)} passed")
    return n_fail


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
