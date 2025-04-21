import os, sys
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['figure.facecolor'] = 'white'

def analyze_splicing_events(
    stats_files,
    output_dir='output_plots', 
    mwu_pval_threshold=0.05, 
    diff_means_threshold=1.0,
    custom_events=None,
    sample_min=3,
    cluster_order=[]
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

            required_cols = {"Feature", "Isoform_1|Length", "Isoform_2|Length", "Event_Type", "Event_Direction", "mwuPval"}
            if not required_cols.issubset(set(col_indices.keys())):
                print(f"Skipping {stats_file}: Missing required columns.")
                continue

            diff_means_col = [col for col in header if col.startswith("DiffMeans_")]
            for line in f:
                cols = line.strip().split("\t")

                feature = cols[col_indices["Feature"]]
                isoform_1 = cols[col_indices["Isoform_1|Length"]]
                isoform_2 = cols[col_indices["Isoform_2|Length"]]
                event_type = cols[col_indices["Event_Type"]]
                event_direction = cols[col_indices["Event_Direction"]]
                try:
                    isoform_1_length = int(isoform_1.split("|")[1]) 
                    isoform_2_length = int(isoform_2.split("|")[1]) 
                except:
                    continue

                event_direction = event_direction_map.get(event_direction, None)
                if event_direction is None:
                    continue
                try:
                    mwu_pval = float(cols[col_indices["mwuPval"]])
                except:
                    mwu_pval = 1

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
                    if (mwu_pval > mwu_pval_threshold) or (diff_means < diff_means_threshold) or (n_count1 < sample_min) or (n_count2 < sample_min):
                        continue

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

        def extract_cluster(comp):
            if "Others-" in comp:
                return comp.split("Others-")[-1]
            return comp

        if cluster_order:
            cluster_order = cluster_order[::-1]
            comparison_to_cluster = {comp: extract_cluster(comp) for comp in comparisons}
            comparisons = sorted(
                comparisons,
                key=lambda c: cluster_order.index(comparison_to_cluster[c]) if comparison_to_cluster[c] in cluster_order else len(cluster_order)
            )

        event_counts = {etype: [len(data_dict.get(comp, {}).get(etype, set())) for comp in comparisons] for etype in event_types}

        fig, ax = plt.subplots(figsize=(10, 6))
        bottom = [0] * len(comparisons)
        colors = ['grey','skyblue','salmon','purple', 'gold','cornflowerblue','lime','blue','red']

        for i, etype in enumerate(event_types):
            ax.barh(comparisons, event_counts[etype], left=bottom, label=etype, color=colors[i % len(colors)])
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

    def extract_cluster(comp):
        if "Others-" in comp:
            return comp.split("Others-")[-1]
        return comp

    # Filter and order comparisons based on cluster_order
    all_comparisons = sorted(set(cond1_dict) & set(cond2_dict) & set(diff_dict))
    if cluster_order:
        cluster_to_comp = {extract_cluster(c): c for c in all_comparisons if extract_cluster(c) in cluster_order}
        comparisons = [cluster_to_comp[cluster] for cluster in cluster_order if cluster in cluster_to_comp]
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
            _, p = stats.mannwhitneyu(c2, c1, alternative='two-sided')
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
    ax1.set_yticklabels([extract_cluster(c) for c in comparisons])
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
    ax2.set_yticklabels([extract_cluster(c) for c in comparisons], fontsize=8)
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

if __name__ == '__main__':
    stats_folder = "dPSI-cluster-events"
    output_dir = "dPSI-cluster-events/output_plots"
    os.makedirs(output_dir, exist_ok=True)

    sys.path.insert(0, "/data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3")
    import altanalyze3.components.long_read.isoform_matrix as iso
    bc_dir = '/data/salomonis-archive/LabFiles/Nathan/Revio/Young-Healthy-Revio/Illumina/SoupX-hg38/cellHarmony-labels-young-old.txt'
    cluster_order = iso.return_cluster_order([bc_dir])

    mwu_pval_threshold = 0.05
    diff_means_threshold = 0.1

    custom_events = None 
    if custom_events is not None:
        custom_events = import_feature_ids(custom_events)

    input_files = glob.glob(os.path.join(stats_folder, "*stats-annotated.txt"))
    if not input_files:
        print(f"No input files found in {stats_folder}. Ensure the folder contains valid files.")
    else:
        analyze_splicing_events(
            input_files,
            output_dir=output_dir, 
            custom_events=custom_events,
            mwu_pval_threshold=mwu_pval_threshold,
            diff_means_threshold=diff_means_threshold,
            sample_min=3,
            cluster_order=cluster_order,
        )