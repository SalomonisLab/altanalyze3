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

# Fixed event-type -> RGB colour map so a given event type is the same colour in
# every plot (positional colouring would shift as the set of event types changes).
EVENT_COLOR_MAP = {
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
DEFAULT_EVENT_COLOR = '#BBBBBB'


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
    write_significant=True,
    states_axis='y',
    mirror=False,
    mirror_min_samples=0,
    mirror_min_events=0,
    strict_thresholds=False,
    additional_fdr_threshold=None,
    figure_caption=None,
):
    if additional_fdr_threshold is not None and (not np.isfinite(additional_fdr_threshold)
                                                or not 0 <= additional_fdr_threshold <= 1):
        raise ValueError('additional_fdr_threshold must be finite and within [0,1]')
    event_data = {}
    overall_event_data = {}
    comparison_samples = {}  # mirror only: comparison -> {condition: largest per-event N}
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
            fdr_col = resolve_pval_column(col_indices, 'fdr')
            if additional_fdr_threshold is not None and fdr_col is None:
                raise ValueError(f'{stats_file}: FDR column required for the requested joint filter')
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
            comp_n = [0, 0]; n_unparsed = 0
            for line in f:
                cols = line.strip().split("\t")
                if mirror:
                    # comparison-level sample count = the largest per-event N of each condition, read from every
                    # row before any event filter
                    for k, ccol in enumerate(condition_cols[:2]):
                        try:
                            comp_n[k] = max(comp_n[k], int(cols[col_indices[ccol]]))
                        except (ValueError, IndexError):
                            n_unparsed += 1

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
                if not np.isfinite(pval) or not np.isfinite(diff_means):
                    continue

                event_tuple = tuple(sorted([isoform_1, isoform_2]))

                if custom_events:
                    if feature not in custom_events:
                        continue
                else:
                    outside = ((pval >= pval_threshold) or (diff_means <= diff_means_threshold)
                               if strict_thresholds else
                               (pval > pval_threshold) or (diff_means < diff_means_threshold))
                    if outside or (n_count1 < sample_min) or (n_count2 < sample_min):
                        continue
                    if additional_fdr_threshold is not None:
                        try: qval = float(cols[col_indices[fdr_col]])
                        except (ValueError, IndexError): continue
                        if not np.isfinite(qval) or (qval >= additional_fdr_threshold if strict_thresholds
                                                    else qval > additional_fdr_threshold):
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
        if mirror:
            comparison_samples[comparison_name] = {condition1_name: comp_n[0], condition2_name: comp_n[1]}
            if n_unparsed:
                print(f"Warning: {comparison_name}: {n_unparsed} N_ values could not be read as integers "
                      f"(not counted in the comparison sample count)")
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
            # states_axis 'y' lists the first state at the top (reversed order); 'x' runs left to right
            cluster_order = cluster_order[::-1] if states_axis == 'y' else list(cluster_order)
            comparison_to_cluster = {comp: extract_cluster(comp, cluster_order) for comp in comparisons}
            comparisons = sorted(
                comparisons,
                key=lambda c: cluster_order.index(comparison_to_cluster[c]) if comparison_to_cluster[c] in cluster_order else len(cluster_order)
            )

        event_counts = {etype: [len(data_dict.get(comp, {}).get(etype, set())) for comp in comparisons] for etype in event_types}

        event_color_map = EVENT_COLOR_MAP
        default_color = DEFAULT_EVENT_COLOR

        if states_axis == 'x':
            _plot_states_on_x(comparisons, event_types, event_counts, event_color_map, default_color,
                              cluster_order, output_dir, output_filename)
            return

        fig, ax = plt.subplots(figsize=(10, 6))
        bottom = [0] * len(comparisons)

        for etype in event_types:
            if etype not in event_color_map:
                print(f"Warning: no colour defined for event type '{etype}'; using {default_color}")
            ax.barh(comparisons, event_counts[etype], left=bottom, label=etype,
                    color=event_color_map.get(etype, default_color))
            bottom = [bottom[j] + event_counts[etype][j] for j in range(len(comparisons))]

        # y-labels = the clean cell-state / cluster name (ASDC), not the full comparison
        # (healthy-ASDC-Others-ASDC). bars were placed at 0..n-1 in `comparisons` order.
        ax.set_yticks(range(len(comparisons)))
        ax.set_yticklabels([extract_cluster(c, cluster_order) for c in comparisons])
        ax.set_xlabel("Count of Splicing Events")
        ax.set_ylabel("Cell state")
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

    if mirror:
        _plot_mirror_states_on_x(event_data, overall_event_data, comparison_samples, cluster_order, output_dir,
                                 mirror_min_samples, mirror_min_events, caption=figure_caption)

    if not event_data:
        print("No events passed the filtering criteria. No PDFs generated.")
        return

    for event_direction, comparison_data in event_data.items():
        plot_stacked_bar(comparison_data, f"Splicing Events for {event_direction}", f"splicing_events_{event_direction}.pdf", cluster_order=cluster_order)

def _plot_mirror_states_on_x(event_data, overall_event_data, comparison_samples, cluster_order, output_dir,
                             min_samples=0, min_events=0, output_stem="splicing_events_inclusion_exclusion_mirror",caption=None):
    """Mirror panel: one bar per cell state along the x axis, inclusion/proximal events stacked upward and
    exclusion/distal events stacked downward (same event-type colours), state labels below the axes.

    A comparison is drawn only when every condition has >= min_samples samples (the largest per-event N of that
    condition in its stats file) and it has >= min_events events (unique isoform pairs summed over event types,
    the height of the splicing_events_all bar). Every comparison, with its counts and any exclusion reason, is
    written to <output_stem>.tsv."""
    from matplotlib.patches import Patch
    from matplotlib.ticker import FuncFormatter
    inc = event_data.get("inclusion_proximal", {})
    exc = event_data.get("exclusion_distal", {})
    comps = list(comparison_samples)
    if cluster_order:
        rank = {cl: i for i, cl in enumerate(cluster_order)}
        comps = sorted(comps, key=lambda c: rank.get(extract_cluster(c, cluster_order), len(cluster_order)))
    count = lambda d, c: sum(len(s) for s in d.get(c, {}).values())
    rows, kept = [], []
    for c in comps:
        reasons = [f"{g} samples {n} < {min_samples}" for g, n in comparison_samples[c].items() if n < min_samples]
        total = count(overall_event_data, c)
        if total < min_events:
            reasons.append(f"events {total} < {min_events}")
        row = dict(comparison=c, cell_state=extract_cluster(c, cluster_order))
        row.update({f"samples_{g}": n for g, n in comparison_samples[c].items()})
        row.update(total_events=total, inclusion_proximal_events=count(inc, c), exclusion_distal_events=count(exc, c),
                   kept=not reasons, excluded_because="; ".join(reasons))
        rows.append(row)
        if not reasons:
            kept.append(c)
    os.makedirs(output_dir, exist_ok=True)
    table_path = os.path.join(output_dir, output_stem + ".tsv")
    pd.DataFrame(rows).to_csv(table_path, sep="\t", index=False)
    print(f"Mirror plot: {len(kept)} of {len(comps)} comparisons kept (min samples per condition {min_samples}, "
          f"min events {min_events}); table: {table_path}")
    if not kept:
        print("No comparison passed the mirror filters. No mirror PDF generated.")
        return

    event_types = sorted({et for d in (inc, exc) for c in kept for et in d.get(c, {})})
    n = len(kept)
    xs = np.arange(n)
    fig, ax = plt.subplots(figsize=(max(6.0, 0.13 * n + 1.4), 4.2))
    up, down = np.zeros(n), np.zeros(n)
    handles = []
    for etype in event_types:
        if etype not in EVENT_COLOR_MAP:
            print(f"Warning: no colour defined for event type '{etype}'; using {DEFAULT_EVENT_COLOR}")
        color = EVENT_COLOR_MAP.get(etype, DEFAULT_EVENT_COLOR)
        vi = np.array([len(inc.get(c, {}).get(etype, ())) for c in kept], dtype=float)
        ve = np.array([len(exc.get(c, {}).get(etype, ())) for c in kept], dtype=float)
        nz = vi > 0
        ax.bar(xs[nz], vi[nz], bottom=up[nz], width=0.8, linewidth=0, color=color)
        nz = ve > 0
        ax.bar(xs[nz], -ve[nz], bottom=-down[nz], width=0.8, linewidth=0, color=color)
        up += vi
        down += ve
        handles.append(Patch(facecolor=color, edgecolor='none', label=f"{etype} ({int(vi.sum()):,} / {int(ve.sum()):,})"))
    ax.axhline(0, color='#000000', linewidth=0.5)
    ax.set_xticks(xs)
    ax.set_xticklabels([extract_cluster(c, cluster_order) for c in kept], rotation=90,
                       fontsize=min(7, max(4, int(400 / max(n, 1)))))
    ax.set_xlim(-0.6, n - 0.4)
    top, bot = max(up.max(), 1.0), max(down.max(), 1.0)
    ax.set_ylim(-1.05 * bot, 1.05 * top)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v):g}"))
    ax.set_ylabel("Splicing events (unique isoform pairs)", fontsize=7)
    ax.tick_params(axis='y', labelsize=6)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.text(1.01, top / (top + bot) * 0.5 + bot / (top + bot), "Inclusion / proximal", transform=ax.transAxes,
            rotation=90, ha='left', va='center', fontsize=6)
    ax.text(1.01, bot / (top + bot) * 0.5, "Exclusion / distal", transform=ax.transAxes,
            rotation=90, ha='left', va='center', fontsize=6)
    ax.legend(handles=handles, title=(f"Event type (inclusion / exclusion; totals {int(up.sum()):,} / "
                                      f"{int(down.sum()):,} in {n} states)"), fontsize=5.5, title_fontsize=6, ncol=5,
              frameon=False, loc='lower left', bbox_to_anchor=(0.0, 1.01), handlelength=0.9, columnspacing=0.8)
    # Illustrator-ready: no white background boxes and no per-artist clip groups (the y limits cover the data)
    fig.patch.set_visible(False)
    ax.patch.set_visible(False)
    for art in ax.get_children():
        if hasattr(art, "set_clip_on"):
            art.set_clip_on(False)
    pdf_path = os.path.join(output_dir, output_stem + ".pdf")
    if caption:
        import textwrap
        fig.text(.05,-.39,'\n'.join(textwrap.wrap(str(caption),width=125)),ha='left',va='top',fontsize=6)
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot: {pdf_path}")

def _plot_states_on_x(comparisons, event_types, event_counts, event_color_map, default_color, cluster_order,
                      output_dir, output_filename):
    """Figure-panel layout: one stacked bar per cell state along the x axis, the event-type legend above the
    axes and no title. Zero-height segments are not drawn (they become empty paths in Illustrator)."""
    n = len(comparisons)
    fig, ax = plt.subplots(figsize=(max(6.0, 0.13 * n + 1.4), 3.0))
    xs = np.arange(n)
    bottom = np.zeros(n)
    for etype in event_types:
        if etype not in event_color_map:
            print(f"Warning: no colour defined for event type '{etype}'; using {default_color}")
        v = np.asarray(event_counts[etype], dtype=float)
        nz = v > 0
        ax.bar(xs[nz], v[nz], bottom=bottom[nz], width=0.8, linewidth=0,
               color=event_color_map.get(etype, default_color), label=f"{etype} ({int(v.sum()):,})")
        bottom += v
    ax.set_xticks(xs)
    ax.set_xticklabels([extract_cluster(c, cluster_order) for c in comparisons], rotation=90,
                       fontsize=min(7, max(4, int(400 / max(n, 1)))))
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylabel("Splicing events (unique isoform pairs)", fontsize=7)
    ax.tick_params(axis='y', labelsize=6)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.legend(title=f"Event type (total {int(bottom.sum()):,} in {int((bottom > 0).sum())} states)", fontsize=5.5,
              title_fontsize=6, ncol=5, frameon=False, loc='lower left', bbox_to_anchor=(0.0, 1.01),
              handlelength=0.9, columnspacing=0.8)
    os.makedirs(output_dir, exist_ok=True)
    output_pdf_path = os.path.join(output_dir, output_filename.replace("/", "_"))
    fig.savefig(output_pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot: {output_pdf_path}")


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

# Canonical hematopoietic lineage cell-state order (top->bottom of the maturation hierarchy). Same order
# used by the lineage/cell-state stacked plots. Used to order the bars by lineage when no --barcode-dir
# cellHarmony labels file is supplied. A cellHarmony order via --barcode-dir still overrides this.
CANONICAL_CELLSTATE_ORDER = [
    # HSPC
    "HSC-1", "HSC-2", "MPP-1", "MPP-2", "MPP-MEP", "LMPP-1", "LMPP-1-cycling", "LMPP-2", "Multilin-1", "Multilin-2",
    # Erythroid
    "MEP-1", "MEP-2", "MEP-Eryth-1", "MEP-Eryth-2", "ERP-1", "ERP-2", "ERP-3", "ERP-4", "ERP-5", "ERP-6", "ERP-7",
    "ERP-8", "Erythroblast-1", "Erythroblast-2", "Erythroblast-3",
    # Megakaryocyte
    "MKP-early", "MKP-late", "MK-Platelet",
    # Mast-Basophil
    "BMCP-1", "BMCP-2",
    # Granulocyte
    "preNeu", "immNeu-1", "immNeu-2",
    # GMP-Myeloid
    "MultiLin-GMP-1", "MultiLin-GMP-2", "MultiLin-GMP-3",
    "Myeloid intermediate 1", "Myeloid intermediate 2", "Myeloid intermediate 3",
    # Monocyte
    "cMOP", "Mono-1", "Mono-2", "Intermediate Mono-1", "Intermediate Mono-2", "Intermediate Mono-3",
    "Classical-Mono", "Non-Classical Mono-1", "Non-Classical Mono-2", "Mac",
    # Dendritic
    "MDP-1", "MDP-2", "MDP-3", "MDP-4", "MDP-5", "MDP-6", "pre-DC-1", "pre-DC-2", "pre-DC-3",
    "cDC1", "cDC2-1", "cDC2-2", "pDC", "ASDC",
    # B-lymphoid
    "CLP", "Pro-B-Early", "Pro-B-Early-cycling", "Pro-B-1", "Pro-B-2", "Pro-B-3",
    "Pro-B-cycling-1", "Pro-B-cycling-2", "Transitional-B-1", "Transitional-B-2", "pre-B",
    "B Memory-1", "B Memory-2", "Plasma Cell",
    # T-NK
    "CD4 TCM", "T CD4 Activated", "T CD4 Naive-1", "T CD4 Naive-2", "T CD8 Naive", "CD8 TEM", "MAIT", "NK-Mature",
]


def derive_cluster_order(input_dir):
    """Without a cellHarmony labels file, build the cluster vocabulary from cluster-events
    filenames (the authoritative <cluster>-Others-<cluster> labels) in the input folder and
    a sibling dPSI-cluster-events folder, ordered by the canonical hematopoietic LINEAGE order
    (clusters not in the canonical list are appended alphabetically at the end)."""
    vocab = set()
    sibling = os.path.join(os.path.dirname(os.path.normpath(input_dir)), "dPSI-cluster-events")
    for d in {input_dir, sibling}:
        for f in glob.glob(os.path.join(d, "*_stats-annotated.txt")):
            comp = os.path.basename(f).replace("_stats-annotated.txt", "")
            if "Others-" in comp:
                vocab.add(comp.split("Others-")[-1])
    known = [c for c in CANONICAL_CELLSTATE_ORDER if c in vocab]
    unknown = sorted(vocab - set(known))
    return known + unknown


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
    parser.add_argument('--cluster-order', default=None,
                        help='Explicit native cell-state order, one state per line; avoids requiring barcode files')
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
    parser.add_argument('--states-axis', choices=['y', 'x'], default='y',
                        help="Axis that lists the cell states in the event-count plots. y (default): horizontal "
                             "bars, one row per state. x: figure-panel layout, one vertical stacked bar per state "
                             "along the x axis, legend above, no title.")
    parser.add_argument('--mirror', action=argparse.BooleanOptionalAction, default=False,
                        help="Also write splicing_events_inclusion_exclusion_mirror.pdf/.tsv: one bar per cell state "
                             "along x, inclusion/proximal events stacked up, exclusion/distal events stacked down, "
                             "state labels at the bottom. Default off.")
    parser.add_argument('--mirror-min-samples', type=int, default=0,
                        help="Mirror plot only: keep a comparison when every condition has at least this many "
                             "samples (largest per-event N_<condition> in its stats file). Default 0 (no filter).")
    parser.add_argument('--mirror-min-events', type=int, default=0,
                        help="Mirror plot only: drop a comparison with fewer than this many splicing events "
                             "(unique isoform pairs over all event types and directions). Default 0 (no filter).")
    parser.add_argument('--strict-thresholds', action='store_true',
                        help='Use p < cutoff and |dPSI| > cutoff; default retains historical inclusive boundaries')
    parser.add_argument('--additional-fdr', type=float,
                        help='Also require native FDR below this threshold, in addition to --pval selection')
    parser.add_argument('--caption-file', help='Optional scientific figure legend for the mirror panel')
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

    if args.cluster_order:
        with open(args.cluster_order) as handle:
            cluster_order = [line.strip() for line in handle if line.strip()]
        if len(cluster_order) != len(set(cluster_order)):
            raise ValueError('--cluster-order contains duplicate state identities')
    elif args.barcode_dir:
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
        states_axis=args.states_axis,
        mirror=args.mirror,
        mirror_min_samples=args.mirror_min_samples,
        mirror_min_events=args.mirror_min_events,
        strict_thresholds=args.strict_thresholds,
        additional_fdr_threshold=args.additional_fdr,
        figure_caption=open(args.caption_file).read() if args.caption_file else None,
    )
