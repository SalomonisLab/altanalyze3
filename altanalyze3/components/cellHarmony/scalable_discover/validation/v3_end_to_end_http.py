"""V3: one scALABLE-discover job, end to end, through the running server's HTTP API only.

  python v3_end_to_end_http.py --base-url http://127.0.0.1:8010 --species human \
      --h5 <a.h5>:<sample> --h5 <b.h5>:<sample> --ambient yes --out <json>

The same calls the browser makes: upload, save QC, run, poll status, then every Explore
view, every PDF, every download, Chat, and the removed differential routes. A call counts
as passing only on HTTP 200 with a non-empty, well-formed body; each check records why.
"""
import argparse, json, sys, time
from pathlib import Path
from urllib.parse import quote

import httpx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import file_record, parse_h5

ap = argparse.ArgumentParser()
ap.add_argument("--base-url", required=True); ap.add_argument("--species", required=True)
ap.add_argument("--h5", action="append", required=True); ap.add_argument("--ambient", default="no")
ap.add_argument("--out", required=True); ap.add_argument("--timeout-min", type=float, default=60)
args = ap.parse_args()
base = args.base_url.rstrip("/")
client = httpx.Client(timeout=600)
checks = []


def check(name, response, test=lambda r: True, detail=""):
    ok = response.status_code == 200 and test(response)
    checks.append({"check": name, "status": response.status_code, "ok": bool(ok),
                   "bytes": len(response.content), "detail": detail})
    return ok


inputs = parse_h5(args.h5)
files = [("files", (Path(p).name, open(p, "rb"), "application/octet-stream")) for p, _ in inputs]
data = {"species": args.species, "reference": "icgs3", "sample_names": [s for _, s in inputs]}
r = client.post(f"{base}/api/jobs", data=data, files=files)
check("upload", r, lambda r: "job_id" in r.json())
job = r.json()["job_id"]
r = client.post(f"{base}/api/jobs/{job}/qc", json={"min_genes": 500, "min_counts": 1000, "min_cells": 0,
                                                  "mit_percent": 15, "ambient_correction": args.ambient,
                                                  "impute_modalities": []})
check("save_qc", r, lambda r: r.json()["qc"]["impute_modalities"] == [])
r = client.post(f"{base}/api/jobs/{job}/qc", json={"impute_modalities": ["adt"]})
checks.append({"check": "qc_refuses_imputation", "status": r.status_code, "ok": r.status_code == 400,
               "bytes": len(r.content), "detail": r.text[:200]})
client.post(f"{base}/api/jobs/{job}/qc", json={"min_genes": 500, "min_counts": 1000, "min_cells": 0,
                                               "mit_percent": 15, "ambient_correction": args.ambient})
started = time.time()
r = client.post(f"{base}/api/jobs/{job}/run")
check("run", r, lambda r: r.json().get("status") == "queued")
status = {}
while time.time() - started < args.timeout_min * 60:
    status = client.get(f"{base}/api/jobs/{job}/status").json()
    if status.get("status") in {"completed", "failed"}:
        break
    time.sleep(5)
elapsed = round(time.time() - started, 1)
checks.append({"check": "pipeline_completed", "status": 200, "ok": status.get("status") == "completed",
               "bytes": 0, "detail": f"{status.get('status')}: {status.get('message')}"})
if status.get("status") != "completed":
    Path(args.out).write_text(json.dumps({"job_id": job, "checks": checks, "status": status}, indent=2, default=str))
    sys.exit(f"pipeline did not complete: {status.get('message')}")

icgs3 = status["icgs3_analysis"]
clusters = icgs3["clusters"]
first = clusters[0]
gene = status.get("default_gene")
umap = client.get(f"{base}/api/jobs/{job}/umap")
check("umap_clusters", umap, lambda r: len(r.json()["query"]) == icgs3["clustered_cells"] and r.json()["reference"] == [],
      "query points equal ICGS3-clustered cells; no reference points")
check("umap_color_by_prediction", client.get(f"{base}/api/jobs/{job}/umap?color_by=ICGS3_cell_state_prediction"),
      lambda r: r.json()["color_by"] == "ICGS3_cell_state_prediction")
check("expression_umap", client.get(f"{base}/api/jobs/{job}/expression?gene={quote(gene)}"),
      lambda r: bool(r.json()), f"gene {gene}")
check("genes", client.get(f"{base}/api/jobs/{job}/genes"), lambda r: len(r.json().get("genes", r.json())) > 0)
check("dotplot_default_markers", client.get(f"{base}/api/jobs/{job}/dotplot"), lambda r: bool(r.json()))
check("combplot", client.get(f"{base}/api/jobs/{job}/combplot?genes={quote(gene)}"), lambda r: bool(r.json()))
check("plot_variables", client.get(f"{base}/api/jobs/{job}/plot-variables"), lambda r: bool(r.json()))
check("display_filters", client.get(f"{base}/api/jobs/{job}/display-filters"), lambda r: bool(r.json()))
heat = client.get(f"{base}/api/jobs/{job}/marker/heatmap.tsv?cells_per_sample=10")
heat_rows = heat.text.count("\n") - 1 if heat.status_code == 200 else 0
check("marker_heatmap_tsv_sampled", heat, lambda r: heat_rows > 0, f"{heat_rows} marker rows")
check("marker_heatmap_tsv_all_columns", client.get(f"{base}/api/jobs/{job}/marker/heatmap.tsv"), lambda r: r.text.count("\n") > 1)
check("marker_heatmap_viewer", client.get(f"{base}/jobs/{job}/marker/heatmap/viewer?cells_per_sample=10"),
      lambda r: "morpheus" in r.text.lower())
check("marker_heatmap_pdf", client.get(f"{base}/api/jobs/{job}/marker/heatmap.pdf"), lambda r: r.content[:4] == b"%PDF")
networks = status["marker_analysis"]["networks"]
net_states = [n["population"] for n in networks]
check("marker_network", client.get(f"{base}/api/jobs/{job}/marker/network?population={quote(net_states[0])}"),
      lambda r: len(r.json()["elements"]) > 0, f"{len(networks)} of {len(clusters)} clusters have a network")
check("marker_network_pdf", client.get(f"{base}/api/jobs/{job}/marker/network/pdf?population={quote(net_states[0])}"),
      lambda r: r.content[:4] == b"%PDF")
check("umap_pdf_clusters", client.get(f"{base}/api/jobs/{job}/umap/pdf?mode=cluster"), lambda r: r.content[:4] == b"%PDF")
check("expression_pdf_umap", client.get(f"{base}/api/jobs/{job}/expression/pdf?gene={quote(gene)}&mode=umap"),
      lambda r: r.content[:4] == b"%PDF")
check("expression_pdf_violin", client.get(f"{base}/api/jobs/{job}/expression/pdf?gene={quote(gene)}&mode=violin"),
      lambda r: r.content[:4] == b"%PDF")
fastcomm = status.get("fastcomm_analysis") or {}
def _counts_only(summary):
    """The fastComm summary without its file paths: a result file must hold no local path."""
    return {k: v for k, v in (summary or {}).items() if isinstance(v, (int, float)) and not isinstance(v, bool)}


checks.append({"check": "fastcomm_enabled", "status": 200, "ok": bool(fastcomm.get("enabled")), "bytes": 0,
               "detail": fastcomm.get("message") or json.dumps(_counts_only(fastcomm.get("summary")))[:200]})
if fastcomm.get("enabled"):
    focus = fastcomm["populations"][0]
    for plot in ("focused_incoming", "focused_outgoing", "cell_state_network", "state_heatmap", "lr_dotplot",
                 "top_table", "per_sample"):
        check(f"fastcomm_{plot}", client.get(f"{base}/api/jobs/{job}/fastcomm/plot?population={quote(focus)}&plot_type={plot}"),
              lambda r: bool(r.json()))
regnet = client.get(f"{base}/api/jobs/{job}/integrated/network?cell_state={quote(first)}&source=marker")
regnet_body = regnet.json() if regnet.status_code == 200 else {}
check("regulatory_network_marker", regnet, lambda r: isinstance(regnet_body, dict),
      json.dumps({k: regnet_body.get(k) for k in ("available", "note", "status")}
                 | {"edges": len(regnet_body.get("edges", []))})[:300])
for key in status.get("artifacts", {}):
    check(f"download_{key}", client.get(f"{base}/api/jobs/{job}/download/{key}"), lambda r: len(r.content) > 0)
check("download_log", client.get(f"{base}/api/jobs/{job}/log"), lambda r: b"Job completed." in r.content)
examples = client.get(f"{base}/api/jobs/{job}/chat-examples")
check("chat_examples", examples, lambda r: len(r.json()["examples"]) > 0 and not r.json()["has_contrast"])
chat = []
for question in examples.json()["examples"][:6]:
    a = client.post(f"{base}/api/jobs/{job}/chat", json={"question": question})
    body = a.json() if a.headers.get("content-type", "").startswith("application/json") else {}
    chat.append({"question": question, "status": a.status_code, "intent": body.get("intent"),
                 "answer": str(body.get("answer", body.get("detail", "")))[:240]})
checks.append({"check": "chat_answers", "status": 200, "ok": all(c["status"] == 200 for c in chat),
               "bytes": 0, "detail": f"{sum(c['status'] == 200 for c in chat)} of {len(chat)} answered 200"})
for path in ("/api/jobs/{job}/differential/status", "/api/jobs/{job}/differential/interactive/summary",
             "/api/meta/reference-preview?species=human&reference=icgs3"):
    r = client.get(base + path.format(job=job))
    checks.append({"check": f"removed {path.split('?')[0]}", "status": r.status_code, "ok": r.status_code == 404,
                   "bytes": len(r.content), "detail": "removed route answers 404"})
r = client.post(f"{base}/api/jobs/{job}/differential", json={"population_col": "ICGS3_cluster",
                                                             "group1_samples": ["a"], "group2_samples": ["b"]})
checks.append({"check": "removed POST differential", "status": r.status_code, "ok": r.status_code in {404, 405},
               "bytes": len(r.content), "detail": "removed route"})

result = {
    "inputs": [dict(file_record(p), sample=s) for p, s in inputs], "species": args.species,
    "ambient_correction": args.ambient, "job_id": job, "pipeline_seconds": elapsed,
    "icgs3": {k: icgs3[k] for k in ("qc_cells", "qc_genes", "clustered_cells", "unassigned_cells", "n_clusters", "clusters")},
    "cell_state_predictions": sorted(set(p["population"] for p in client.get(
        f"{base}/api/jobs/{job}/umap?color_by=ICGS3_cell_state_prediction").json()["query"])),
    "marker_heatmap_rows_sampled_view": heat_rows,
    "marker_networks": len(networks),
    "fastcomm": {"status": fastcomm.get("status"), "sample_key": fastcomm.get("sample_key"),
                 "summary_counts": _counts_only(fastcomm.get("summary"))},
    "bundle": {k: v for k, v in (status.get("bundle") or {}).items() if k in ("status", "n_cells", "threshold", "reason", "seconds", "returncode")},
    "checks_passed": sum(c["ok"] for c in checks), "checks_total": len(checks),
    "failed_checks": [c for c in checks if not c["ok"]],
    "checks": checks, "chat": chat,
}
Path(args.out).write_text(json.dumps(result, indent=2, default=str))
print(json.dumps({k: v for k, v in result.items() if k not in ("checks", "chat", "cell_state_predictions")},
                 indent=2, default=str)[:4000])
