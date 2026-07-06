import os, sys
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import re

plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['pdf.fonttype'] = 42  # Embed fonts as TrueType (editable text)
plt.rcParams['ps.fonttype'] = 42   # Same for PostScript
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['figure.facecolor'] = 'white'


def parse_isoform_length(isoform_label):
    """Parse the protein length from an isoform token like 'ENST00000681979|973' or
    'ENST00000681979|973(NMD)'. Returns the leading integer, or None if unparseable.
    (The plain int(split('|')[1]) form throws on the '(NMD)' suffix and drops the row.)"""
    try:
        tail = isoform_label.split("|", 1)[1]
    except Exception:
        return None
    m = re.match(r"(\d+)", tail)
    return int(m.group(1)) if m else None


def extract_cluster(comp, cluster_order=None):
    """Return the cell-state/cluster label embedded in a comparison name.

    Handles both dPSI naming conventions:
      cluster-events  : <group>-<cluster>-Others-<cluster>  (e.g. young-MEP-2-Others-MEP-2)
      covariate-events: <group1>-<group2>-<cluster>         (e.g. young-aged-cDC2-1)

    With cluster_order (the known label vocabulary), the cluster is the longest known
    label that is a hyphen/space-bounded suffix of the name; falls back to the legacy
    'Others-' split when no vocabulary is supplied.
    """
    if cluster_order:
        cands = [cl for cl in cluster_order
                 if comp == cl or comp.endswith("-" + cl) or comp.endswith(" " + cl)]
        if cands:
            return max(cands, key=len)
    if "Others-" in comp:
        return comp.split("Others-")[-1]
    return comp


def resolve_pval_column(col_indices, pval_type):
    """Pick the p-value column for the requested type ('raw'|'fdr'), auto-detecting the
    statistic family present (ebayes or mwu). Returns None if no match is found."""
    for family in ("ebayes", "mwu"):
        col = f"{family}AdjPval" if pval_type == "fdr" else f"{family}Pval"
        if col in col_indices:
            return col
    return None


def analyze_splicing_events(
    stats_files,
    output_dir='output_plots',
    pval_threshold=0.05,
    pval_type='raw',
    diff_means_threshold=1.0,
    custom_events=None,
    sample_min=3,
    cluster_order=[],
    write_significant=True
):
    event_data = {}  
    overall_event_data = {}  
    cond1_dict = {}
    cond2_dict = {}
    diff_dict = {}

    event_direction_map = {
        "inclusion": "inclusion_proximal",
        "proximal": "inclusion_proximal",
        "exclusion": "exclusion_distal",
        "distal": "exclusion_distal",
        "null": None
    }

    for stats_file in stats_files:
        #if ("HSC" not in stats_file) and ("MPP" not in stats_file): continue
        comparison_name = os.path.basename(stats_file).replace("_stats-annotated.txt", "")
        significant_events = 0
        with open(stats_file, 'r') as f:
            header = f.readline().strip().split("\t")
            col_indices = {col: i for i, col in enumerate(header)}
            condition_cols = [col for col in header if col.startswith("N_")]
            condition1_name = condition_cols[0].replace("N_", "")
            condition2_name = condition_cols[1].replace("N_", "")

            pval_col = resolve_pval_column(col_indices, pval_type)
            required_cols = {"Feature", "Isoform_1|Length", "Isoform_2|Length", "Event_Type", "Event_Direction"}
            if pval_col is None or not required_cols.issubset(set(col_indices.keys())):
                print(f"Skipping {stats_file}: missing required columns or no '{pval_type}' p-value column.")
                continue

            sig_out = None
            if write_significant:
                os.makedirs(output_dir, exist_ok=True)
                sig_out = open(os.path.join(output_dir, comparison_name + "_significant-results.tsv"), "w")
                sig_out.write("SourceFile\t" + "\t".join(header) + "\n")

            diff_means_col = [col for col in header if col.startswith("DiffMeans_")]
            for line in f:
                cols = line.strip().split("\t")

                feature = cols[col_indices["Feature"]]
                isoform_1 = cols[col_indices["Isoform_1|Length"]]
                isoform_2 = cols[col_indices["Isoform_2|Length"]]
                event_type = cols[col_indices["Event_Type"]]
                event_direction = cols[col_indices["Event_Direction"]]
                isoform_1_length = parse_isoform_length(isoform_1)
                isoform_2_length = parse_isoform_length(isoform_2)
                if isoform_1_length is None or isoform_2_length is None:
                    continue

                event_direction = event_direction_map.get(event_direction, None)
                if event_direction is None:
                    continue
                try:
                    pval = float(cols[col_indices[pval_col]])
                except:
                    pval = 1

                n_count1 = int(cols[col_indices[condition_cols[0]]])
                n_count2 = int(cols[col_indices[condition_cols[1]]])

                if not diff_means_col:
                    print(f"Skipping {stats_file}: No 'DiffMeans_*' column found.")
                    continue

                try: 
                    diff_means = abs(float(cols[col_indices[diff_means_col[0]]]))
                except:
                    continue

                event_tuple = tuple(sorted([isoform_1, isoform_2]))

                if custom_events:
                    if feature not in custom_events:
                        continue
                else:
                    if (pval > pval_threshold) or (diff_means < diff_means_threshold) or (n_count1 < sample_min) or (n_count2 < sample_min):
                        continue

                if sig_out is not None:
                    sig_out.write(comparison_name + "\t" + line.strip() + "\n")

                if event_direction not in event_data:
                    event_data[event_direction] = {}
                if comparison_name not in event_data[event_direction]:
                    event_data[event_direction][comparison_name] = {}
                if event_type not in event_data[event_direction][comparison_name]:
                    event_data[event_direction][comparison_name][event_type] = set()
                significant_events += 1

                event_data[event_direction][comparison_name][event_type].add(event_tuple)

                if comparison_name not in overall_event_data:
                    overall_event_data[comparison_name] = {}
                if event_type not in overall_event_data[comparison_name]:
                    overall_event_data[comparison_name][event_type] = set()

                overall_event_data[comparison_name][event_type].add(event_tuple)

                if event_direction == "exclusion_distal":
                    isoform_2_length, isoform_1_length = isoform_1_length, isoform_2_length

                cond1_dict.setdefault(comparison_name, []).append(isoform_1_length)
                cond2_dict.setdefault(comparison_name, []).append(isoform_2_length)
                diff_dict.setdefault(comparison_name, []).append(isoform_2_length - isoform_1_length)

        if sig_out is not None:
            sig_out.close()
        print(comparison_name, significant_events, 'unique events')

    #print("Comparisons seen:", len(overall_event_data))
    #print("Sample of comparison keys:", list(overall_event_data.keys())[:5])
    #print("Cluster order:", cluster_order)

    if cond1_dict and cond2_dict:
        plot_protein_length_distributions(
            cond1_dict, cond2_dict, diff_dict, output_dir, cluster_order=cluster_order
        )

    def plot_stacked_bar(data_dict, title, output_filename, cluster_order=[]):
        comparisons = list(data_dict.keys())
        event_types = sorted({etype for comp in data_dict.values() for etype in comp})

        if cluster_order:
            cluster_order = cluster_order[::-1]
            comparison_to_cluster = {comp: extract_cluster(comp, cluster_order) for comp in comparisons}
            comparisons = sorted(
                comparisons,
                key=lambda c: cluster_order.index(comparison_to_cluster[c]) if comparison_to_cluster[c] in cluster_order else len(cluster_order)
            )

        event_counts = {etype: [len(data_dict.get(comp, {}).get(etype, set())) for comp in comparisons] for etype in event_types}

        # Fixed event-type -> RGB colour map so a given event type is the same colour in
        # every plot (positional colouring would shift as the set of event types changes).
        event_color_map = {
            'Alt C-Terminal Exon':     '#808080',  # grey
            "Alt3' SS":                '#87CEEB',  # skyblue
            "Alt5' SS":                '#FA8072',  # salmon
            'AltPolyA':                '#800080',  # purple
            'AltPromoter':             '#FFD700',  # gold
            'Cassette Exon':           '#6495ED',  # cornflowerblue
            'Intron Retention':        '#00FF00',  # lime
            'Mutually-Exclusive Exon': '#0000FF',  # blue
            'Trans-Splicing':          '#FF0000',  # red
        }
        default_color = '#BBBBBB'

        fig, ax = plt.subplots(figsize=(10, 6))
        bottom = [0] * len(comparisons)

        for etype in event_types:
            if etype not in event_color_map:
                print(f"Warning: no colour defined for event type '{etype}'; using {default_color}")
            ax.barh(comparisons, event_counts[etype], left=bottom, label=etype,
                    color=event_color_map.get(etype, default_color))
            bottom = [bottom[j] + event_counts[etype][j] for j in range(len(comparisons))]

        ax.set_xlabel("Count of Splicing Events")
        ax.set_ylabel("Comparisons (Stats Files)")
        ax.set_title(title)
        fontsize = min(8, max(3, int(300 / len(comparisons))))
        ax.tick_params(axis='y', labelsize=fontsize)
        ax.set_ylim(-0.5, len(comparisons) - 0.5)
        ax.legend(title="Event Type", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()

        os.makedirs(output_dir, exist_ok=True)
        output_pdf_path = os.path.join(output_dir, output_filename.replace("/", "_"))
        plt.savefig(output_pdf_path)
        plt.close()
        print(f"Saved plot: {output_pdf_path}")

    plot_stacked_bar(overall_event_data, "Splicing Events Across All Comparisons", "splicing_events_all.pdf", cluster_order=cluster_order)

    if not event_data:
        print("No events passed the filtering criteria. No PDFs generated.")
        return

    for event_direction, comparison_data in event_data.items():
        plot_stacked_bar(comparison_data, f"Splicing Events for {event_direction}", f"splicing_events_{event_direction}.pdf", cluster_order=cluster_order)

def plot_protein_length_distributions(cond1_dict, cond2_dict, diff_dict, output_dir, cluster_order=[]):
    """Generate and save HORIZONTAL box plots for protein length distributions and differences across all comparisons using IQR. 
       Includes statistical significance tests and annotations for each comparison.
    """
    import pandas as pd
    cluster_order = cluster_order[::-1]
    
    def format_p_value(p):
        return "< 1e-100" if p < 1e-100 else f"{p:.2e}"

    def get_iqr_limits(data, scale=2):
        q1, q3 = np.percentile(data, [25, 75])
        iqr = q3 - q1
        lower = max(q1 - (scale * iqr), min(data))
        upper = min(q3 + (scale * iqr), max(data))
        return lower, upper

    # Order comparisons by cluster_order, keeping ALL comparisons (unknown clusters last).
    # Uses the shared extract_cluster so covariate-events names (no 'Others-') also match.
    all_comparisons = sorted(set(cond1_dict) & set(cond2_dict) & set(diff_dict))
    if cluster_order:
        comparisons = sorted(
            all_comparisons,
            key=lambda c: cluster_order.index(extract_cluster(c, cluster_order))
            if extract_cluster(c, cluster_order) in cluster_order else len(cluster_order)
        )
    else:
        comparisons = all_comparisons

    ### --- Protein Length Distribution Plot (Horizontal) --- ###
    fig1, ax1 = plt.subplots(figsize=(6, max(6, len(comparisons) * 0.5)))
    box_data = []
    positions = []
    bracket_positions = []
    p_texts = []

    for i, comp in enumerate(comparisons):
        c1 = np.array(cond1_dict[comp])  # Cluster
        c2 = np.array(cond2_dict[comp])  # Others
        box_data.extend([c2, c1])
        positions.extend([2 * i + 1, 2 * i + 2])

        try:
            # Paired test: c1 (isoform_1) and c2 (isoform_2) are matched per event
            _, p = stats.wilcoxon(c1, c2, alternative='two-sided')
        except:
            p = 1
        p_texts.append(format_p_value(p) if p < 0.05 else None)
        ymax = np.percentile(np.concatenate([c1, c2]), 95)
        bracket_positions.append((2 * i + 1, 2 * i + 2, ymax))

    bp = ax1.boxplot(box_data, positions=positions, vert=False, widths=0.6, patch_artist=True,
                     boxprops=dict(color="black"), medianprops=dict(color="black"), showfliers=False)

    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor("grey" if i % 2 == 0 else "lightblue")

    import matplotlib.patches as mpatches
    legend_handles = [
        mpatches.Patch(color='grey', label='Others'),
        mpatches.Patch(color='lightblue', label='Cluster')
    ]
    ax1.legend(handles=legend_handles, loc='lower right', fontsize=8)

    ax1.set_yticks([np.mean([positions[2*i], positions[2*i+1]]) for i in range(len(comparisons))])
    ax1.set_yticklabels([extract_cluster(c, cluster_order) for c in comparisons])
    ax1.set_xlabel("Protein Length (Amino Acids)")
    ax1.set_title("Protein Length Distributions")

    for (y1, y2, x), ptxt in zip(bracket_positions, p_texts):
        if ptxt:
            ax1.plot([x + 5, x + 10, x + 10, x + 5], [y1, y1, y2, y2], color='black')
            ax1.text(x + 12, (y1 + y2) / 2, f"p = {ptxt}", va='center', fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "protein_length_distribution_all_comparisons.pdf"))
    plt.close()

    ### --- Protein Length Difference Plot (Horizontal) --- ###
    fig2, ax2 = plt.subplots(figsize=(8, len(comparisons) * 0.5)) # tighter height
    diff_data = []
    p_texts_diff = []

    for comp in comparisons:
        c1 = np.array(cond1_dict[comp])  # Cluster
        c2 = np.array(cond2_dict[comp])  # Others
        d = c1 - c2  # Cluster - Others
        diff_data.append(d)
        try:
            _, p = stats.wilcoxon(d, alternative='two-sided')
        except:
            p = 1
        p_texts_diff.append(format_p_value(p) if p < 0.05 else None)

    positions = np.arange(1, len(comparisons) + 1)

    bp = ax2.boxplot(diff_data, positions=positions, vert=False, widths=0.6, patch_artist=True,
                     boxprops=dict(facecolor="purple"), medianprops=dict(color="black"), showfliers=False)

    # Get fixed right margin from axis limits after plotting
    x_lim = ax2.get_xlim()
    x_text = x_lim[1] + 5  # Fixed offset outside the plot area

    for i, (d, pos) in enumerate(zip(diff_data, positions)):
        direction = "↑ Cluster" if np.median(d) > 0 else "↑ Others"
        color = "blue" if np.median(d) > 0 else "grey"
        ax2.text(x_text, pos, direction, va='center', fontsize=7, color=color)

        if p_texts_diff[i]:
            q3 = np.percentile(d, 75)
            xlim = ax2.get_xlim()
            pad = (xlim[1] - xlim[0]) * 0.01  # 1% of plot width
            ax2.text(q3 + pad, pos - 0.25, f"p = {p_texts_diff[i]}", va='center', fontsize=8, color="black", clip_on=True)


    ax2.set_yticks(positions)
    ax2.set_yticklabels([extract_cluster(c, cluster_order) for c in comparisons], fontsize=8)
    ax2.set_xlabel("Difference in Protein Length (Amino Acids)")
    ax2.set_title("Difference in Protein Length (Cluster - Others)")
    ax2.axvline(0, color='red', linestyle='dashed', linewidth=1)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "protein_length_difference_all_comparisons.pdf"))
    plt.close()

    print("Saved updated horizontal multi-comparison protein length plots with IQR and significance annotations.")
def import_feature_ids(file_path):
    with open(file_path, 'r') as f:
        feature_ids = [line.strip() for line in f if line.strip()]
    print('features to filter for imported...')
    return feature_ids

def derive_cluster_order(input_dir):
    """Without a cellHarmony labels file, build the cluster vocabulary from cluster-events
    filenames (the authoritative <cluster>-Others-<cluster> labels) in the input folder and
    a sibling dPSI-cluster-events folder, ordered alphabetically."""
    vocab = set()
    sibling = os.path.join(os.path.dirname(os.path.normpath(input_dir)), "dPSI-cluster-events")
    for d in {input_dir, sibling}:
        for f in glob.glob(os.path.join(d, "*_stats-annotated.txt")):
            comp = os.path.basename(f).replace("_stats-annotated.txt", "")
            if "Others-" in comp:
                vocab.add(comp.split("Others-")[-1])
    return sorted(vocab)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description="Splice-event bar charts and protein-length plots from AltAnalyze3 dPSI "
                    "*_stats-annotated.txt files (dPSI-cluster-events or dPSI-covariate-events)."
    )
    parser.add_argument('--input-dir', required=True,
                        help="Folder of *_stats-annotated.txt files (e.g. .../dPSI-covariate-events).")
    parser.add_argument('--group1', default=None,
                        help="First condition/group to select (e.g. young). For cluster-events this "
                             "is the group compared against Others.")
    parser.add_argument('--group2', default=None,
                        help="Second condition/group for covariate-events (e.g. aged). Either "
                             "filename order (group1-group2 or group2-group1) is matched.")
    parser.add_argument('--dpsi', type=float, default=0.1,
                        help="Minimum |DiffMeans| (dPSI) threshold. Default 0.1.")
    parser.add_argument('--pval', type=float, default=0.05,
                        help="p-value threshold for filtering. Default 0.05.")
    parser.add_argument('--pval-type', choices=['raw', 'fdr'], default='raw',
                        help="Filter on raw p-values (ebayesPval/mwuPval) or FDR-adjusted "
                             "(ebayesAdjPval/mwuAdjPval). Default raw.")
    parser.add_argument('--sample-min', type=int, default=3,
                        help="Minimum N per condition. Default 3.")
    parser.add_argument('--barcode-dir', default=None,
                        help="Optional cellHarmony labels file; sets cluster (row) order via "
                             "isoform_matrix.return_cluster_order. If omitted, cluster order is "
                             "derived alphabetically from cluster-events filenames.")
    parser.add_argument('--custom-events', default=None,
                        help="Optional file of Feature IDs to restrict to (bypasses thresholds).")
    parser.add_argument('--output-dir', default=None,
                        help="Output folder. Default: <input-dir>/<comparison-name>.")
    parser.add_argument('--exclude-clusters', default=None,
                        help="Comma-separated cluster labels to drop from the row order "
                             "(e.g. 'ERP-7,MEP-Eryth-2').")
    parser.add_argument('--write-significant', action=argparse.BooleanOptionalAction, default=True,
                        help="Write <comparison>_significant-results.tsv per file. Default on "
                             "(--no-write-significant to disable).")
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        sys.exit(f"Input dir not found: {args.input_dir}")

    def comp_name(path):
        return os.path.basename(path).replace("_stats-annotated.txt", "")

    def selected(path):
        c = comp_name(path)
        if args.group1 and args.group2:
            return (c.startswith(f"{args.group1}-{args.group2}-") or
                    c.startswith(f"{args.group2}-{args.group1}-"))
        if args.group1:
            return c.startswith(f"{args.group1}-")
        return True

    all_files = sorted(glob.glob(os.path.join(args.input_dir, "*_stats-annotated.txt")))
    input_files = [f for f in all_files if selected(f)]
    if not input_files:
        sys.exit(f"No *_stats-annotated.txt files in {args.input_dir} match the requested comparison.")
    print(f"Selected {len(input_files)} of {len(all_files)} files in {args.input_dir}")

    if args.group1 and args.group2:
        comparison_label = f"{args.group1}-{args.group2}"
    elif args.group1:
        comparison_label = f"{args.group1}-vs-Others"
    else:
        comparison_label = "all-comparisons"
    output_dir = args.output_dir or os.path.join(args.input_dir, comparison_label)
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output -> {output_dir}")

    if args.barcode_dir:
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
        sys.path.insert(0, repo_root)
        import altanalyze3.components.long_read.isoform_matrix as iso
        cluster_order = iso.return_cluster_order([args.barcode_dir])
    else:
        cluster_order = derive_cluster_order(args.input_dir)
        print(f"No --barcode-dir given; derived {len(cluster_order)} cluster labels "
              f"(alphabetical) from cluster-events filenames.")

    if args.exclude_clusters:
        drop = {c.strip() for c in args.exclude_clusters.split(",") if c.strip()}
        n_before = len(input_files)
        # resolve each file's cluster against the full order BEFORE removing the excluded labels
        input_files = [f for f in input_files
                       if extract_cluster(comp_name(f), cluster_order) not in drop]
        cluster_order = [c for c in cluster_order if c not in drop]
        print(f"Excluded clusters {sorted(drop)}: dropped {n_before - len(input_files)} file(s)")

    custom_events = import_feature_ids(args.custom_events) if args.custom_events else None

    analyze_splicing_events(
        input_files,
        output_dir=output_dir,
        pval_threshold=args.pval,
        pval_type=args.pval_type,
        diff_means_threshold=args.dpsi,
        custom_events=custom_events,
        sample_min=args.sample_min,
        cluster_order=cluster_order,
        write_significant=args.write_significant,
    )