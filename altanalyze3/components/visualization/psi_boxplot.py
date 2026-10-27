"""Plot donor-level inclusion/(inclusion+exclusion) from complete saved count tables.

The supplied roster is independent of observed nonzero counts. Missing donor/isoform
records are errors; observed zero support produces undefined PSI, never zero PSI.
This module does not query, filter, reassign or normalize native RNA measurements.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import textwrap


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def add_arguments(parser):
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--counts', help='Native long TSV: cohort, donor, isoform_id, assigned_RNA')
    source.add_argument('--psi-table', help='Complete native event PSI table: cohort, donor, event_uid, PSI')
    parser.add_argument('--counts-sha256', help='Require the exact source table hash')
    parser.add_argument('--psi-sha256', help='Require the exact native PSI table hash')
    parser.add_argument('--event-uid', help='Literal event identity required with --psi-table')
    parser.add_argument('--donor-roster', required=True, help='JSON mapping each condition to its complete donor list')
    parser.add_argument('--inclusion-isoforms', nargs='+')
    parser.add_argument('--exclusion-isoforms', nargs='+')
    parser.add_argument('--cohort-order', nargs='+', required=True)
    parser.add_argument('--labels', nargs='+', help=r'One display label per condition; literal \n wraps a label')
    parser.add_argument('--condition-column', default='cohort')
    parser.add_argument('--donor-column', default='donor')
    parser.add_argument('--stat-donor-column', help='Native PSI table donor ID: tests use each donor mean once; plotted values stay unchanged')
    parser.add_argument('--unit-label', default='donors', help='Unit for displayed circles/n (e.g. pseudobulks)')
    parser.add_argument('--isoform-column', default='isoform_id')
    parser.add_argument('--count-column', default='assigned_RNA')
    parser.add_argument('--title', default='Inclusion PSI in donor pseudobulks')
    parser.add_argument('--ylabel', default='PSI (inclusion / [inclusion + exclusion])')
    parser.add_argument('--caption-file', help='Complete UTF-8 scientific figure caption')
    parser.add_argument('--reference-condition', help='Compare every other condition with this reference')
    parser.add_argument('--annotations-table', help='Existing statistics TSV: reference, comparison, p, test, correction; no test is rerun')
    parser.add_argument('--marker-table', help='Complete donor-style TSV: cohort, donor, marker, label; changes point shapes only')
    parser.add_argument('--stat-correction', choices=['holm', 'bh', 'none'], default='holm',
                        help='Correction across reference comparisons (default: Holm)')
    parser.add_argument('--y-min', type=float, default=0.0)
    parser.add_argument('--y-max', type=float, default=1.03)
    parser.add_argument('--show-clipped-points', action='store_true', help='Mark values outside the axis at its boundary; retain true values in TSV')
    parser.add_argument('--width', type=float, default=8.4, help='Figure width in inches')
    parser.add_argument('--height', type=float, default=6.1, help='Figure height in inches')
    parser.add_argument('--caption-wrap', type=int, default=116, help='Characters per caption line')
    parser.add_argument('--annotation-fontsize', type=float, default=12, help='Statistical annotation font size')
    parser.add_argument('--ylabel-pad', type=float, default=4, help='Y-axis title padding in points')
    parser.add_argument('--out', required=True, help='New PDF; existing outputs are never overwritten')


def donor_psi_records(path, roster, cohorts, inclusion, exclusion,
                      columns=('cohort', 'donor', 'isoform_id', 'assigned_RNA')):
    """Retain the exact roster and demand an explicit observed count for every role."""
    if len(cohorts) != len(set(cohorts)) or set(cohorts) != set(roster):
        raise ValueError('Condition order must cover the complete independent donor roster exactly')
    ids = list(inclusion) + list(exclusion)
    if not inclusion or not exclusion or len(ids) != len(set(ids)):
        raise ValueError('Inclusion/exclusion IDs must be nonempty, unique and disjoint')
    for condition in cohorts:
        donors = roster[condition]
        if not isinstance(donors, list) or not donors or len(donors) != len(set(donors)):
            raise ValueError('Each condition requires a nonempty, unique donor roster')
    expected = {(c, d, i) for c in cohorts for d in roster[c] for i in ids}
    counts = {}
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not set(columns) <= set(reader.fieldnames or []):
            raise ValueError('Count table lacks required columns: ' + str(columns))
        for row in reader:
            key = tuple(row[c] for c in columns[:3])
            if key in counts:
                raise ValueError('Duplicate donor/isoform count record: ' + str(key))
            value = float(row[columns[3]])
            if not math.isfinite(value) or value < 0:
                raise ValueError('Counts must be finite and nonnegative: ' + str(key))
            counts[key] = value
    if set(counts) != expected:
        raise ValueError('Count/roster identities differ: missing=' + str(sorted(expected-set(counts))) +
                         '; unexpected=' + str(sorted(set(counts)-expected)))
    records = []
    for condition in cohorts:
        for donor in roster[condition]:
            inc = sum(counts[condition, donor, i] for i in inclusion)
            exc = sum(counts[condition, donor, i] for i in exclusion)
            total = inc + exc
            records.append(dict(condition=condition, donor=donor, inclusion_reads=inc,
                                exclusion_reads=exc, total_informative_reads=total,
                                PSI=inc/total if total else None,
                                status='OBSERVED_INCLUSION_FRACTION' if total else
                                'UNDEFINED_ZERO_INFORMATIVE_READS_RETAINED'))
    return records


def native_psi_records(path, roster, cohorts, event_uid, identity_column='donor', biological_donor_column=None):
    """Consume literal native PSI; missing values and every roster role are retained."""
    if len(cohorts) != len(set(cohorts)) or set(cohorts) != set(roster):
        raise ValueError('Condition order must cover the complete independent donor roster exactly')
    if any(not donors or len(donors) != len(set(donors)) for donors in roster.values()):
        raise ValueError('Donor rosters must be nonempty and unique')
    expected = {(c,d) for c in cohorts for d in roster[c]}
    records = {}
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'cohort',identity_column,'event_uid','PSI'} | ({biological_donor_column} if biological_donor_column else set())
        if not required <= set(reader.fieldnames or []):
            raise ValueError('Native PSI table lacks required identity/value columns')
        for row in reader:
            key = row['cohort'],row[identity_column]
            if key in records: raise ValueError('Duplicate native PSI donor role: '+str(key))
            if row['event_uid'] != event_uid: raise ValueError('Native PSI event identity mismatch')
            value = float(row['PSI']) if row['PSI'].strip().lower() not in ('','na','nan') else None
            if value is not None and (not math.isfinite(value) or not 0 <= value <= 1):
                raise ValueError('PSI must be missing or finite within [0,1]')
            records[key] = dict(condition=key[0], donor=key[1], event_uid=event_uid, PSI=value,
                                status='NATIVE_EVENT_PSI' if value is not None else 'NATIVE_PSI_UNDEFINED_RETAINED')
            if biological_donor_column:
                if not row[biological_donor_column].strip(): raise ValueError('Missing biological donor identity')
                records[key]['biological_donor'] = row[biological_donor_column]
    if set(records) != expected:
        raise ValueError('Native PSI table does not retain the complete donor roster: '+str(expected-set(records)))
    return [records[c,d] for c in cohorts for d in roster[c]]


def run(args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    # Axis headroom is for annotations; native_psi_records still requires PSI in [0, 1].
    if not (0 <= args.y_min < args.y_max and args.y_max <= 1.2):
        raise ValueError('Require 0 <= y-min < y-max <= 1.2')
    if args.width <= 0 or args.height <= 0:
        raise ValueError('Figure dimensions must be positive')
    roster = json.loads(Path(args.donor_roster).read_text())
    source_hash = sha256(args.psi_table or args.counts)
    expected_hash = args.psi_sha256 if args.psi_table else args.counts_sha256
    if expected_hash and source_hash != expected_hash:
        raise ValueError('PSI source SHA256 mismatch' if args.psi_table else 'Count source SHA256 mismatch')
    if args.psi_table:
        if not args.event_uid or args.inclusion_isoforms or args.exclusion_isoforms:
            raise ValueError('Native PSI input requires --event-uid and no isoform count roles')
        records = native_psi_records(args.psi_table, roster, args.cohort_order, args.event_uid,
                                     args.donor_column, args.stat_donor_column)
    else:
        records = donor_psi_records(args.counts, roster, args.cohort_order,
                                   args.inclusion_isoforms, args.exclusion_isoforms,
                                   (args.condition_column, args.donor_column, args.isoform_column, args.count_column))
    marker_table = getattr(args, 'marker_table', None)
    marker_legend = {}
    if marker_table:
        from matplotlib.markers import MarkerStyle
        with Path(marker_table).open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not {'cohort', 'donor', 'marker', 'label'} <= set(reader.fieldnames or []):
                raise ValueError('Marker table requires cohort, donor, marker and label')
            styles = {}
            for row in reader:
                key = row['cohort'], row['donor']
                if key in styles or not row['label'].strip():
                    raise ValueError('Duplicate donor style or missing marker label')
                MarkerStyle(row['marker'])
                if row['label'] in marker_legend and marker_legend[row['label']] != row['marker']:
                    raise ValueError('Each marker label must have one consistent shape')
                marker_legend[row['label']] = row['marker']
                styles[key] = row
        if set(styles) != {(r['condition'], r['donor']) for r in records}:
            raise ValueError('Marker table must retain every donor role exactly')
        for row in records:
            style = styles[row['condition'], row['donor']]
            row['point_marker'], row['point_marker_label'] = style['marker'], style['label']
    labels = args.labels or args.cohort_order
    if len(labels) != len(args.cohort_order):
        raise ValueError('Provide one label per condition')
    caption = (Path(args.caption_file).read_text().strip() if args.caption_file else
               ('Each circle is a donor; values are literal native junction-event PSI. ' if args.psi_table else
                'Each circle is a donor; PSI = inclusion/(inclusion+exclusion) within the explicit isoform universe. ')+
               'Boxes show median/IQR; whiskers extend to observations within 1.5 IQR. All defined points are shown. '
               'Zero-support PSI remains undefined; n reports defined/roster donors. No statistical test is applied.')
    out = Path(args.out)
    outputs = [out, out.with_suffix('.donor_PSI.tsv'), out.with_suffix('.manifest.json')]
    if args.reference_condition or args.annotations_table:
        outputs.append(out.with_suffix('.statistics.tsv'))
        if args.stat_donor_column:
            outputs.append(out.with_suffix('.statistical_units.tsv'))
    if any(p.exists() for p in outputs):
        raise FileExistsError('Existing figure/source outputs retained: ' + str(outputs))
    data = [[r['PSI'] for r in records if r['condition'] == c and r['PSI'] is not None]
            for c in args.cohort_order]
    stats = []; statistical_units = []
    if args.reference_condition:
        from altanalyze3.components.visualization.violin import _pairwise_group_tests, _significance_stars
        reference = args.reference_condition
        if reference not in roster:
            raise ValueError('Reference condition is absent from the complete roster')
        stat_roster = {c:set(roster[c]) for c in args.cohort_order}
        stat_data = data
        if args.stat_donor_column:
            stat_data = []
            for c in args.cohort_order:
                by_donor = {r['biological_donor']:[] for r in records if r['condition'] == c}
                for row in records:
                    if row['condition'] == c and row['PSI'] is not None:
                        by_donor.setdefault(row['biological_donor'], []).append(row['PSI'])
                stat_data.append([float(np.mean(values)) for values in by_donor.values() if values])
                statistical_units.extend(dict(condition=c, biological_donor=donor,
                                               mean_saved_PSI=float(np.mean(values)) if values else None,
                                               defined_pseudobulks=len(values),
                                               status='DEFINED_DONOR_MEAN' if values else 'NO_DEFINED_NATIVE_PSI_RETAINED')
                                         for donor,values in by_donor.items())
                stat_roster[c] = {r['biological_donor'] for r in records if r['condition'] == c}
        for condition in args.cohort_order:
            if condition != reference and stat_roster[reference] & stat_roster[condition]:
                raise ValueError('Independent-group test requires disjoint reference/comparison donors')
        stats = _pairwise_group_tests(
            [dict(label=c, values=v) for c,v in zip(args.cohort_order, stat_data)],
            test='mannwhitney', correction=args.stat_correction, reference=reference)
        for row in stats:
            row['reference'] = reference
            row['comparison'] = row['group_2'] if row['group_1'] == reference else row['group_1']
            row['label'] = _significance_stars(row['p_adjusted'])
        # A supplied caption must not retain a contradictory no-test statement.
        if args.caption_file and 'no statistical test' in caption.lower():
            raise ValueError('Caption contradicts the requested statistical comparison')
        if not args.caption_file:
            caption = caption.replace('No statistical test is applied.', '').strip()
        caption += (' Two-sided Mann–Whitney U tests versus '+reference+'; '+
                    ('Holm-adjusted' if args.stat_correction == 'holm' else
                     'BH-adjusted' if args.stat_correction == 'bh' else 'unadjusted')+
                    ' P values across '+str(len(stats))+' comparisons. '
                    'ns: P>0.05; *: P≤0.05; **: P≤0.01; ***: P≤0.001. '
                    'Only defined donor PSI enters each test; zero-support donors remain in the roster.')
        if args.stat_donor_column:
            caption += ' Tests use the mean of available saved PSI values within each donor; circles show individual pseudobulks.'
    if args.annotations_table:
        if args.reference_condition:
            raise ValueError('Use existing annotations or newly computed reference tests, not both')
        from altanalyze3.components.visualization.violin import _significance_stars
        with Path(args.annotations_table).open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not {'reference','comparison','p','test','correction'} <= set(reader.fieldnames or []):
                raise ValueError('Existing statistics table lacks required comparison/method fields')
            stats = list(reader)
        seen = set(); references = set()
        for row in stats:
            reference, comparison = row['reference'], row['comparison']
            if reference not in roster or comparison not in roster or reference == comparison or comparison in seen:
                raise ValueError('Existing statistics do not match unique displayed comparisons')
            raw_p = row['p'].strip()
            p = None if raw_p.lower() in ('', 'na', 'nan') else float(raw_p)
            if p is not None and (not math.isfinite(p) or not 0 <= p <= 1):
                raise ValueError('Existing annotation P must be finite within [0,1]')
            if set(roster[reference]) & set(roster[comparison]):
                raise ValueError('Reference/comparison donor identities overlap')
            row['p'] = p; row['label'] = 'NA' if p is None else _significance_stars(p)
            seen.add(comparison); references.add(reference)
        if len(references) != 1 or seen != set(roster)-references:
            raise ValueError('Existing statistics must cover every non-reference condition exactly once')
        caption += ' Symbols use supplied native raw P values: ns >0.05; * ≤0.05; ** ≤0.01; *** ≤0.001. No statistical test is rerun.'
        if any(row['p'] is None for row in stats):
            caption += ' NA denotes an unavailable native test, not a nonsignificant result.'
    with plt.rc_context({'font.family': 'Arial', 'font.size': 11, 'pdf.fonttype': 42}):
        fig, ax = plt.subplots(figsize=(args.width, args.height))
        fig.subplots_adjust(left=.12, right=.98, bottom=.36 if marker_table else .32, top=.87)
        positions = list(range(1, len(data)+1))
        ax.boxplot([v if v else [np.nan] for v in data], positions=positions, widths=.5,
                   patch_artist=True, whis=1.5, showfliers=False,
                   boxprops={'facecolor': '#eeeeee', 'edgecolor': 'black'},
                   medianprops={'color': 'black', 'linewidth': 1.4},
                   whiskerprops={'color': 'black'}, capprops={'color': 'black'})
        for j, (c, v) in enumerate(zip(args.cohort_order, data), 1):
            point_x = np.full(len(v), j)+np.linspace(-.065, .065, len(v))
            if marker_table:
                observed = [r for r in records if r['condition'] == c and r['PSI'] is not None]
                for marker in dict.fromkeys(r['point_marker'] for r in observed):
                    mask = np.array([r['point_marker'] == marker for r in observed])
                    ax.scatter(point_x[mask], np.asarray(v)[mask], marker=marker,
                               s=42, facecolor='black', edgecolor='black', zorder=4, clip_on=False)
            else:
                ax.scatter(point_x, v, s=42, facecolor='black', edgecolor='black', zorder=4)
            ax.text(j, 1.018, 'n = '+str(len(v))+'/'+str(len(roster[c])), ha='center',
                    transform=ax.get_xaxis_transform(), fontsize=10)
            if args.show_clipped_points:
                below = sum(x < args.y_min for x in v); above = sum(x > args.y_max for x in v)
                for count, y, marker in [(below,args.y_min,'v'),(above,args.y_max,'^')]:
                    if count:
                        ax.scatter(np.full(count,j)+np.linspace(-.12,.12,count),np.full(count,y),
                                   marker=marker,s=28,facecolor='white',edgecolor='black',clip_on=False,zorder=6)
                if below:
                    ax.text(j,.03,str(below)+' below '+str(args.y_min),transform=ax.get_xaxis_transform(),
                            fontsize=8,ha='center')
        ax.set_xlim(.5, len(data)+.5)
        for row in stats:
            j = args.cohort_order.index(row['comparison'])+1
            label = (('P unavailable\nNA' if row['p'] is None else
                      'P = '+format(row['p'],'.3g')+'\n'+row['label'])
                     if args.annotations_table else row['label'])
            ax.text(j, .975, label, transform=ax.get_xaxis_transform(),
                    ha='center', va='top', fontsize=args.annotation_fontsize, fontweight='bold')
        ax.set_ylim(args.y_min, args.y_max)
        # Matplotlib expands the view to include requested ticks. A fixed tick
        # at 1 silently defeated --y-max < 1 even though the manifest kept it.
        ax.set_yticks(np.linspace(args.y_min, min(1., args.y_max), 6))
        ax.set_xticks(positions, [s.replace('\\n', '\n') for s in labels])
        ax.set_ylabel(args.ylabel.replace('\\n', '\n'), labelpad=args.ylabel_pad)
        ax.grid(axis='y', color='#dddddd', linewidth=.5)
        ax.set_axisbelow(True)
        ax.spines[['top', 'right']].set_visible(False)
        fig.suptitle(args.title, fontsize=14, y=.97)
        if marker_table:
            from matplotlib.lines import Line2D
            handles = [Line2D([], [], marker=marker, linestyle='none', color='black',
                              markersize=5, label=label) for label, marker in marker_legend.items()]
            fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.55, .29),
                       ncol=len(handles), frameon=False, fontsize=8)
        fig.text(.055, .015, textwrap.fill(caption, width=args.caption_wrap), fontsize=8.1,
                 ha='left', va='bottom', linespacing=1.25)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out)
        plt.close(fig)
    with outputs[1].open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)
    if stats:
        with outputs[3].open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(stats[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(stats)
        if statistical_units:
            with outputs[4].open('w') as handle:
                writer=csv.DictWriter(handle,fieldnames=list(statistical_units[0]),delimiter='\t',lineterminator='\n')
                writer.writeheader();writer.writerows(statistical_units)
    manifest = dict(status='COMPLETE_DONOR_PSI_PLOT', records=records,
                    inclusion_isoforms=args.inclusion_isoforms, exclusion_isoforms=args.exclusion_isoforms,
                    donor_roster=roster, cohort_order=args.cohort_order, y_limits=[args.y_min, args.y_max],
                    counts_sha256=None if args.psi_table else source_hash,
                    native_PSI_table_sha256=source_hash if args.psi_table else None, event_uid=args.event_uid,
                    donor_roster_sha256=sha256(args.donor_roster),
                    caption_sha256=sha256(args.caption_file) if args.caption_file else None,
                    code_sha256=sha256(__file__), PDF_sha256=sha256(out), caption=caption,
                    statistics=stats,
                    existing_statistics_sha256=sha256(args.annotations_table) if args.annotations_table else None,
                    display_unit=args.unit_label,
                    statistical_unit='mean saved PSI per biological donor' if args.stat_donor_column else args.unit_label,
                    clipped_points_marked=args.show_clipped_points,
                    figure_inches=[args.width,args.height],
                    all_defined_points_in_view=all(r['PSI'] is None or args.y_min <= r['PSI'] <= args.y_max
                                                   for r in records))
    if marker_table:
        manifest['marker_table_sha256'] = sha256(marker_table)
        manifest['marker_legend'] = marker_legend
    outputs[2].write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    return run(parser.parse_args(argv))


if __name__ == '__main__':
    main()
