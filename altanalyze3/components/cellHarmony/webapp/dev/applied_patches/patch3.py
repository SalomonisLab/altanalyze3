P = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/app.py"
src = open(P).read()

# Read-only GET handlers that do synchronous numpy / pandas / file work. As
# `async def` they ran ON the event loop, so one slow figure blocked every other
# request in the app. FastAPI runs a plain `def` handler in its threadpool.
NAMES = [
    "reference_preview", "umap", "expression", "marker_network", "grn_network",
    "fastcomm_network", "fastcomm_plot", "job_genes", "dotplot", "combplot",
    "plot_variables", "job_display_filters", "umap_pdf", "expression_pdf",
    "marker_heatmap_tsv", "marker_heatmap_viewer", "marker_heatmap_pdf",
    "marker_network_pdf", "differential_heatmap_data", "differential_volcano_data",
    "differential_go_data", "differential_network_data", "differential_table_data",
    "differential_gene_data", "differential_gene_pdf", "differential_rendered_pdf",
]
for name in NAMES:
    old = f"    async def {name}("
    new = f"    def {name}("
    assert src.count(old) == 1, f"{name}: {src.count(old)} matches"
    assert f"await {name}(" not in src, f"{name} is awaited somewhere"
    src = src.replace(old, new)

open(P, "w").write(src)
print(f"patch3 applied: {len(NAMES)} handlers moved to the threadpool")
