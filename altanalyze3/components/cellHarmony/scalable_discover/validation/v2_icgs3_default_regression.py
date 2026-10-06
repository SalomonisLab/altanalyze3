"""V2: the ICGS3 edits leave its default path unchanged.

  python v2_icgs3_default_regression.py --orig <HEAD copy of ICGS.py> \
      --h5 <a.h5>:<sample> --h5 <b.h5>:<sample> --work <scratch> --out <json>

The HEAD copy must sit beside ProteinCoding-Hs-Mm.txt, which it reads relative to itself.
Three runs on one counts h5ad (QC-only cellHarmony-lite output, ambient on):
  orig      HEAD ICGS.py, ICGS3Config defaults; reads BioMarkers from its hard-coded path
  new       edited ICGS.py, ICGS3Config defaults; reads the bundled BioMarkers copy
  networks  edited ICGS.py, defaults + export_marker_networks=True
orig == new on every file proves the default path, the bundled BioMarkers file included, is
unchanged. new == networks proves the option only adds the network files.
"""
import argparse, importlib.util, json, sys
from pathlib import Path
import anndata as ad
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import file_record, parse_h5

ap = argparse.ArgumentParser()
ap.add_argument("--orig", required=True); ap.add_argument("--h5", action="append", required=True)
ap.add_argument("--work", required=True); ap.add_argument("--out", required=True)
args = ap.parse_args()
WORK = Path(args.work); WORK.mkdir(parents=True, exist_ok=True)
from altanalyze3.components.cellHarmony import cellHarmony_lite
from altanalyze3.components.clustering import ICGS as new

counts_path = WORK / "icgs3_input_counts.h5ad"
if not counts_path.exists():
    _, q = cellHarmony_lite.combine_and_align_h5(
        h5_files=parse_h5(args.h5), cellharmony_ref=None, output_dir=str(WORK / "qc"), min_genes=500,
        min_cells=0, min_counts=1000, mit_percent=15, ambient_correct_cutoff="auto", ambient_memory_efficient=True,
        concat_on_disk=True, concat_batch_size=1, stream_10x_inputs=True, return_adata=True)
    ad.AnnData(X=q.layers["counts"], obs=q.obs.copy(), var=q.var.copy()).write_h5ad(counts_path, compression="lzf")

spec = importlib.util.spec_from_file_location("ICGS_orig", args.orig)
orig = importlib.util.module_from_spec(spec); sys.modules["ICGS_orig"] = orig; spec.loader.exec_module(orig)


def run(module, name, **extra):
    out = WORK / name
    module.run_icgs3(module.ICGS3Config(input_paths=[str(counts_path)], output_dir=str(out), species="Hs", **extra))
    return out


FILES = ["icgs3_clusters.tsv", "icgs3_cell_barcode_clusters.tsv", "sNMF/icgs3_pagerank_downsampling.tsv",
         "MarkerFinder/icgs3_marker_heatmap_markers.tsv", "MarkerFinder/icgs3_marker_heatmap_redundant_markers.tsv",
         "MarkerFinder/icgs3_marker_heatmap_fold_matrix.centroids.tsv", "UMAPs/icgs3_umap.tsv",
         "GO-Elite/icgs3_biomarker_enrichment.tsv", "GO-Elite/icgs3_cell_state_predictions.tsv"]
runs = {"orig": run(orig, "orig"), "new": run(new, "new"), "networks": run(new, "networks", export_marker_networks=True)}
table = {f: {k: file_record(v / f)["sha256"] for k, v in runs.items()} for f in FILES}
net_dir = runs["networks"] / "MarkerFinder" / "icgs3_marker_heatmap_networks"
n_clusters = sum(1 for _ in open(runs["new"] / "GO-Elite" / "icgs3_cell_state_predictions.tsv")) - 1
result = {
    "inputs": {"h5": [dict(file_record(p), sample=s) for p, s in parse_h5(args.h5)],
               "icgs3_input_counts": file_record(counts_path), "original_module": file_record(args.orig)},
    "files": table,
    "orig_equals_new_all_files": all(t["orig"] is not None and t["orig"] == t["new"] for t in table.values()),
    "new_equals_networks_all_files": all(t["new"] is not None and t["new"] == t["networks"] for t in table.values()),
    "network_tsv_files": len(list(net_dir.glob("*.tsv"))) if net_dir.exists() else 0,
    "clusters_with_biomarker_prediction": n_clusters,
    "default_runs_wrote_networks": any((runs[k] / "MarkerFinder" / "icgs3_marker_heatmap_networks").exists()
                                       for k in ("orig", "new")),
}
Path(args.out).write_text(json.dumps(result, indent=2))
print(json.dumps({k: v for k, v in result.items() if k not in ("files", "inputs")}, indent=2))
