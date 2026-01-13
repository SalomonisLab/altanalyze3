import argparse
import os
import scanpy as sc
import matplotlib.pyplot as plt
import math

plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['figure.facecolor'] = 'white'


parser = argparse.ArgumentParser(description="Generate UMAP plots from a preprocessed h5ad file")
parser.add_argument('--h5ad', type=str, required=True, help='Input h5ad file with UMAP coordinates and annotations')
parser.add_argument('--outdir', type=str, default='figures', help='Directory to save output plots')
parser.add_argument('--genes', type=str, default='', help='Space-separated list of genes to plot')
args = parser.parse_args()

# Parse genes
genes = args.genes.split() if args.genes else []

# Load data
print(f"Reading AnnData from: {args.h5ad}")
adata = sc.read_h5ad(args.h5ad)

# Collect all plots
plot_entries = []

if "leiden" in adata.obs:
    plot_entries.append(("leiden", True))  # categorical
if "AssignedCellState" in adata.obs:
    plot_entries.append(("AssignedCellState", True))

available_genes = [gene for gene in genes if gene in adata.var_names]
missing_genes = [gene for gene in genes if gene not in adata.var_names]
if missing_genes:
    print(f"Warning: These genes were not found in adata.var_names and will be skipped: {missing_genes}")

plot_entries.extend([(gene, False) for gene in available_genes])  # False = continuous

n_plots = len(plot_entries)
if n_plots == 0:
    print("No valid UMAP plots to generate. Exiting.")
    exit()

# Compute grid size (wrap to new row after 5 columns)
cols = 5
rows = math.ceil(n_plots / cols)

fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 5 * rows))
axes = axes.flatten() if n_plots > 1 else [axes]

# Plot each UMAP in the appropriate subplot
for i, (feature, is_categorical) in enumerate(plot_entries):
    sc.pl.embedding(
        adata,
        basis='umap',
        color=feature,
        ax=axes[i],
        show=False,
        legend_loc='on data' if is_categorical else 'right margin',
        legend_fontsize=6,
        legend_fontweight='normal' if is_categorical else 'bold',
        title=feature
    )

# Hide unused subplots
for ax in axes[n_plots:]:
    ax.set_visible(False)

# Save to PDF
output_path = os.path.join(args.outdir, "umap_all_panels.pdf")
plt.tight_layout()
plt.savefig(output_path, bbox_inches="tight")
print(f"Saved combined UMAP plot to: {output_path}")

# === Second plot: vertically stacked violin plots across leiden clusters ===
if "leiden" in adata.obs and available_genes:
    print("Generating vertically stacked violin plots for genes across 'leiden' clusters...")

    import matplotlib.pyplot as plt

    n_genes = len(available_genes)
    fig, axes = plt.subplots(nrows=n_genes, ncols=1, figsize=(6, 2.5 * n_genes), constrained_layout=True)

    if n_genes == 1:
        axes = [axes]  # ensure iterable if only one gene

    for i, gene in enumerate(available_genes):
        sc.pl.violin(
            adata,
            keys=gene,
            groupby="leiden",
            stripplot=False,
            ax=axes[i],
            show=False
        )

    violin_path = os.path.join(args.outdir, "violin_leiden_stacked.pdf")
    plt.savefig(violin_path)
    plt.close()
    print(f"Saved violin plot to: {violin_path}")
else:
    print("Skipping violin plot: 'leiden' clustering or valid genes not found.")
