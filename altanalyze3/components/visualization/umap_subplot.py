import argparse
import os
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import math

plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['figure.facecolor'] = 'white'


_BUILTIN_GRADIENTS = {
    "blue_yellow_red": ["#1f4ac0", "#f5e616", "#d4181c"],
    "blue_white_red": ["#1f4ac0", "#ffffff", "#d4181c"],
}


def _resolve_cmap(name: str):
    key = name.lower()
    if key in _BUILTIN_GRADIENTS:
        return LinearSegmentedColormap.from_list(key, _BUILTIN_GRADIENTS[key])
    return plt.get_cmap(name)


parser = argparse.ArgumentParser(description="Generate UMAP plots from a preprocessed h5ad file")
parser.add_argument('--h5ad', type=str, required=True, help='Input h5ad file with UMAP coordinates and annotations')
parser.add_argument('--outdir', type=str, default='figures', help='Directory to save output plots')
parser.add_argument('--genes', type=str, default='', help='Space-separated list of genes to plot')
parser.add_argument('--cmap', type=str, default='blue_yellow_red',
                    help='Colormap for continuous features. Built-ins: blue_yellow_red, blue_white_red. '
                         'Any matplotlib name (e.g. RdYlBu_r, viridis, coolwarm) also works.')
parser.add_argument('--vmin', type=str, default='p2',
                    help="Lower clip for continuous features. 'pN' = N-th percentile (e.g. p2). "
                         "Or a number. Use 'auto' for scanpy default.")
parser.add_argument('--vmax', type=str, default='p98',
                    help="Upper clip for continuous features. 'pN' = N-th percentile. Or a number. 'auto' for default.")
parser.add_argument('--contrast', type=float, default=1.0,
                    help="Scale factor on the displayed range, computed per-gene from --vmin/--vmax percentiles. "
                         ">1 narrows the range (brighter); <1 widens it (muted). 1.0 leaves vmin/vmax untouched.")
parser.add_argument('--out-name', type=str, default='umap_all_panels',
                    help="Base filename for the output PDF/PNG (no extension). Defaults to 'umap_all_panels'.")
args = parser.parse_args()


def _parse_clip(val):
    if val is None or str(val).lower() == 'auto':
        return None
    return val  # scanpy accepts 'pN' strings and floats directly

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

cmap_obj = _resolve_cmap(args.cmap)
vmin_val = _parse_clip(args.vmin)
vmax_val = _parse_clip(args.vmax)


def _resolve_clip_value(val, gene_values):
    if val is None:
        return None
    if isinstance(val, str) and val.startswith('p'):
        import numpy as np
        return float(np.percentile(gene_values, float(val[1:])))
    return float(val)


def _scaled_clip(feature, base_vmin, base_vmax, contrast):
    if contrast == 1.0 and not (isinstance(base_vmin, str) and base_vmin.startswith('p')) \
       and not (isinstance(base_vmax, str) and base_vmax.startswith('p')):
        return base_vmin, base_vmax
    import numpy as np
    col = adata[:, feature].X
    vals = col.toarray().ravel() if hasattr(col, 'toarray') else np.asarray(col).ravel()
    lo = _resolve_clip_value(base_vmin, vals)
    hi = _resolve_clip_value(base_vmax, vals)
    if lo is None or hi is None or contrast == 1.0:
        return lo, hi
    center = (lo + hi) / 2.0
    half = (hi - lo) / 2.0 / max(contrast, 1e-6)
    return center - half, center + half


# Plot each UMAP in the appropriate subplot
for i, (feature, is_categorical) in enumerate(plot_entries):
    extra = {}
    if not is_categorical:
        extra['cmap'] = cmap_obj
        lo, hi = _scaled_clip(feature, vmin_val, vmax_val, args.contrast)
        if lo is not None:
            extra['vmin'] = lo
        if hi is not None:
            extra['vmax'] = hi
    sc.pl.embedding(
        adata,
        basis='umap',
        color=feature,
        ax=axes[i],
        show=False,
        legend_loc='on data' if is_categorical else 'right margin',
        legend_fontsize=6,
        legend_fontweight='normal' if is_categorical else 'bold',
        title=feature,
        **extra,
    )

# Hide unused subplots
for ax in axes[n_plots:]:
    ax.set_visible(False)

# Save to PDF
output_path = os.path.join(args.outdir, f"{args.out_name}.pdf")
png_path = os.path.join(args.outdir, f"{args.out_name}.png")
plt.tight_layout()
plt.savefig(output_path, bbox_inches="tight")
plt.savefig(png_path, bbox_inches="tight", dpi=150)
print(f"Saved combined UMAP plot to: {output_path}")
print(f"Saved combined UMAP plot (PNG) to: {png_path}")

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
