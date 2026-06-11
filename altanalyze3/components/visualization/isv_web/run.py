"""CLI launcher for the ISV web app.

  python -m altanalyze3.components.visualization.isv_web.run \
      --metadata /path/to/metadata.txt \
      --gene_model /path/to/Hs_Ensembl_exon.txt \
      --gene_symbol /path/to/Hs_Ensembl-annotations.txt \
      [--run_dir <dir containing per-sample h5ads + gff-output>] [--port 8050]

Files are assumed LOCAL. run_dir defaults to the metadata's directory.
"""
from __future__ import annotations

import argparse
import os
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(description="Launch the interactive ISV web viewer.")
    ap.add_argument("--metadata", required=True, help="metadata TSV (uid/library/reverse/groups)")
    ap.add_argument("--gene_model", required=True, help="Hs_Ensembl_exon.txt (exon coords)")
    ap.add_argument("--gene_symbol", required=True, help="Hs_Ensembl-annotations.txt (ENSG->symbol)")
    ap.add_argument("--run_dir", default=None, help="dir with per-sample *-isoform.h5ad + gff-output (default: metadata dir)")
    ap.add_argument("--gff_output", default=None, help="gff-output dir (default: <run_dir>/gff-output)")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8050)
    ap.add_argument("--build_reads_index", action="store_true",
                    help="build/refresh the per-gene molecule retrieval index at startup, then serve "
                         "(makes the read-level view retrieve in ~ms; idempotent)")
    a = ap.parse_args(argv)

    run_dir = a.run_dir or os.path.dirname(os.path.abspath(a.metadata))

    # Import here so the engine import cost is paid once, after arg parsing.
    from . import data_api as da
    from .server import create_app
    import uvicorn

    print(f"[isv_web] loading run context from {run_dir} ...", flush=True)
    ctx = da.RunContext(run_dir, a.metadata, a.gene_model, a.gene_symbol,
                        gff_output_dir=a.gff_output).load_all()

    if a.build_reads_index:
        from . import precompute_reads_index as pri
        print("[isv_web] building per-gene molecule retrieval index ...", flush=True)
        pri.build_reads_index(ctx.viewer_sample_dict, ctx.barcode_sample_dict, ctx.reads_index_dir)

    ctx.prewarm_reads()   # make the first read-level click instant
    print(f"[isv_web] ready: {len(ctx.samples)} samples, {len(ctx.groups)} groups, "
          f"{len(ctx.cell_types)} cell types, {len(ctx.all_genes)} genes.", flush=True)
    print(f"[isv_web] serving on http://{a.host}:{a.port}", flush=True)

    app = create_app(ctx)
    uvicorn.run(app, host=a.host, port=a.port, log_level="info")


if __name__ == "__main__":
    main()
