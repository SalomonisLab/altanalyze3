"""Interactive ISV (isoform structure viewer) web app.

A self-contained FastAPI + D3/SVG web interface over the AltAnalyze3 isoform structure viewer
engine and its fast precompute caches (gene_indexes_v2/). Independent of the cellHarmony web app.

Run:  python -m altanalyze3.components.visualization.isv_web.run --metadata <metadata.txt> --port 8050
"""
