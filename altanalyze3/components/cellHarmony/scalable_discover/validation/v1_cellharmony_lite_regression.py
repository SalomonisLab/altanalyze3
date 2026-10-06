"""V1: the QC-only edit to cellHarmony_lite leaves the reference path unchanged, and the
QC-only mode returns exactly the cells and values the reference path holds before alignment.

  python v1_cellharmony_lite_regression.py --orig <HEAD copy of cellHarmony_lite.py> \
      --h5 <a.h5>:<sample> --h5 <b.h5>:<sample> --reference <states .txt> --work <scratch> --out <json>

Both runs use the arguments flask/pipeline.py:run_cellharmony_pipeline passes, ambient on.
"""
import argparse, importlib.util, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import adata_digest, file_record, matrix_digest, parse_h5

ap = argparse.ArgumentParser()
ap.add_argument("--orig", required=True); ap.add_argument("--h5", action="append", required=True)
ap.add_argument("--reference", required=True); ap.add_argument("--work", required=True); ap.add_argument("--out", required=True)
args = ap.parse_args()
H5 = parse_h5(args.h5); WORK = Path(args.work)

spec = importlib.util.spec_from_file_location("cellHarmony_lite_orig", args.orig)
orig = importlib.util.module_from_spec(spec); sys.modules["cellHarmony_lite_orig"] = orig; spec.loader.exec_module(orig)
from altanalyze3.components.cellHarmony import cellHarmony_lite as new


def run(module, ref, outdir):
    kw = dict(h5_files=H5, h5ad_file=None, cellharmony_ref=ref, output_dir=str(outdir), export_cptt=False,
              export_h5ad=False, min_genes=500, min_cells=0, min_counts=1000, mit_percent=15,
              generate_umap=False, save_adata=False, unsupervised_cluster=False, gene_translation_file=None,
              metacell_align=False, ambient_correct_cutoff="auto", ambient_memory_efficient=True,
              concat_on_disk=True, concat_batch_size=1, stream_10x_inputs=True, return_adata=True)
    if ref is not None:
        kw.update(alignment_mode="cosine", min_alignment_score=0.4)
    return module.combine_and_align_h5(**kw)


df_o, ad_o = run(orig, args.reference, WORK / "orig")
df_n, ad_n = run(new, args.reference, WORK / "new")
s_o, s_n = adata_digest(ad_o), adata_digest(ad_n)
assign_o = file_record(WORK / "orig" / "cellHarmony_lite_assignments.txt")
assign_n = file_record(WORK / "new" / "cellHarmony_lite_assignments.txt")
none_df, ad_q = run(new, None, WORK / "qconly")
s_q = adata_digest(ad_q)
sub_q = ad_q[ad_n.obs_names]
per_matrix = {"X": matrix_digest(sub_q.X) == s_n["X"]}
for k in ad_n.layers:
    per_matrix[f"layers/{k}"] = (k in sub_q.layers) and matrix_digest(sub_q.layers[k]) == s_n["layers"][k]
result = {
    "inputs": {"h5": [dict(file_record(p), sample=s) for p, s in H5], "reference": file_record(args.reference),
               "original_module": file_record(args.orig)},
    "reference_path_identical": {
        "assignments_file_sha256_equal": assign_o["sha256"] == assign_n["sha256"],
        "assignments_rows": int(len(df_n)), "adata_digest_equal": s_o == s_n, "orig": s_o, "new": s_n},
    "qc_only_mode": {
        "returned_frame_is_none": none_df is None, "qc_only_shape": s_q["shape"], "aligned_shape": s_n["shape"],
        "aligned_cells_contained_in_qc_only": int(len(ad_n.obs_names.intersection(ad_q.obs_names))),
        "aligned_cells": int(ad_n.n_obs), "matrices_equal_on_aligned_cells": per_matrix,
        "var_names_equal": s_q["var_names"] == s_n["var_names"],
        "first_50_barcodes_identical_order": list(ad_n.obs_names[:50]) == list(sub_q.obs_names[:50])},
}
Path(args.out).write_text(json.dumps(result, indent=2))
print(json.dumps({k: v for k, v in result.items() if k != "inputs"}, indent=2, default=str)[:3000])
