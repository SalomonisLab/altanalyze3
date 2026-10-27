"""One-page donor PSI panels using the existing native PSI boxplot renderer."""
import argparse
import csv
import json
import math
from pathlib import Path
import textwrap

from . import psi_boxplot


def add_arguments(parser):
    parser.add_argument('--panel-manifest', required=True,
                        help='JSON with complete ordered panels: gene, event_label and native_argv')
    parser.add_argument('--manifest-sha256')
    parser.add_argument('--columns', type=int, default=4)
    parser.add_argument('--width', type=float, default=18)
    parser.add_argument('--height', type=float, default=20)
    parser.add_argument('--panel-range-step', type=float,
                        help='Enclose each panel observed PSI range by outward rounding to this step')
    parser.add_argument('--out', required=True)


def observed_psi_limits(records, step):
    """Round outward without changing, dropping or imputing any measurement."""
    if not math.isfinite(step) or not 0 < step <= 1:
        raise ValueError('PSI range step must be finite and in (0,1]')
    values=[r['PSI'] for r in records if r['PSI'] is not None]
    if not values:
        raise ValueError('No observed PSI defines an axis range; all donor roles are retained')
    minimum,maximum=min(values),max(values)
    scaled_low,scaled_high=minimum/step,maximum/step
    # Source H5AD PSI is commonly float32: e.g. nominal0.8 is stored as
    # 0.8000000119. Keep that exact point within the axis without expanding
    # its displayed upper tick to0.9 solely because of representation noise.
    low=(min(minimum,round(scaled_low)*step)
         if math.isclose(scaled_low,round(scaled_low),rel_tol=0,abs_tol=1e-6)
         else math.floor(scaled_low)*step)
    high=(max(maximum,round(scaled_high)*step)
          if math.isclose(scaled_high,round(scaled_high),rel_tol=0,abs_tol=1e-6)
          else math.ceil(scaled_high)*step)
    low,high=max(0.,low),min(1.,high)
    if math.isclose(low,high,rel_tol=0,abs_tol=step*1e-6):
        if high<1:high=min(1.,high+step)
        else:low=max(0.,low-step)
    return low,high


def run(args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.lines import Line2D

    path = Path(args.panel_manifest)
    source_hash = psi_boxplot.sha256(path)
    if args.manifest_sha256 and source_hash != args.manifest_sha256:
        raise ValueError('Panel manifest SHA256 mismatch')
    specification = json.loads(path.read_text())
    panels = specification['panels']
    genes = [p['gene'] for p in panels]
    if not panels or len(set(genes)) != len(genes):
        raise ValueError('Require a nonempty complete panel list with unique gene identities')
    if specification.get('required_genes') != genes:
        raise ValueError('Panel order/identities must match the independently required roster')
    if args.columns < 1 or args.width <= 0 or args.height <= 0:
        raise ValueError('Positive grid dimensions are required')
    out = Path(args.out)
    outputs = [out, out.with_suffix('.manifest.json'), out.with_suffix('.donor_PSI.tsv')]
    if any(p.exists() for p in outputs):
        raise FileExistsError('Existing composite outputs retained')
    parser = argparse.ArgumentParser()
    psi_boxplot.add_arguments(parser)
    results = []
    rows = math.ceil(len(panels)/args.columns)
    with plt.rc_context({'font.family': 'Arial', 'font.size': 10, 'pdf.fonttype': 42}):
        fig, axes = plt.subplots(rows, args.columns, figsize=(args.width,args.height), squeeze=False)
        spacing=1.15 if args.panel_range_step else .92
        height_units=rows+(rows-1)*spacing
        header_space=.8 if args.panel_range_step else .4
        top=min(.925,(.955*height_units+header_space*.12)/(height_units+header_space))
        fig.subplots_adjust(left=.055,right=.985,top=top,bottom=.12,wspace=.42,hspace=spacing)
        try:
            for index, (panel, ax) in enumerate(zip(panels, axes.flat)):
                native = parser.parse_args(panel['native_argv'])
                if not native.psi_table or native.reference_condition:
                    raise ValueError('Composite accepts literal native PSI and existing statistics only')
                if native.event_uid != panel['event_uid']:
                    raise ValueError('Panel/native event identities differ')
                if args.panel_range_step:
                    roster=json.loads(Path(native.donor_roster).read_text())
                    original=psi_boxplot.native_psi_records(native.psi_table,roster,native.cohort_order,
                                                            native.event_uid,native.donor_column,native.stat_donor_column)
                    native.y_min,native.y_max=observed_psi_limits(original,args.panel_range_step)
                result = psi_boxplot.run(native, panel_ax=ax)
                if not result['all_defined_points_in_view']:
                    raise ValueError('Composite must show every defined source PSI value')
                if results and (result['donor_roster'] != results[0]['donor_roster'] or
                                result['cohort_order'] != results[0]['cohort_order']):
                    raise ValueError('Every panel must retain the same independent cohort/donor roster')
                gene_id=native.event_uid.split(':',1)[0]
                expected_label=panel['gene']+':'+native.event_uid.split(':',1)[1].replace('|'+gene_id+':','|')
                if panel['event_label'] != expected_label:
                    raise ValueError('Every panel requires its full symbol:junction|junction header')
                # Titles are individual panel text, not page-level suptitles.
                event_lines = panel['event_label'].replace('|','|\n',1)
                if args.panel_range_step:
                    # Statistics stay above the data axes so the requested
                    # observed-value limits need no artificial headroom.
                    for label in ax.texts:
                        if label.get_text().startswith(('P = ','P unavailable')):
                            label.set_y(1.27)
                        elif label.get_text().startswith('n = '):
                            label.set_y(1.05)
                            label.set_fontsize(9)
                    ticks=np.arange(round(native.y_min/args.panel_range_step),
                                    round(native.y_max/args.panel_range_step)+1)*args.panel_range_step
                    if len(ticks)>6:
                        ticks=np.unique(np.r_[ticks[::math.ceil((len(ticks)-1)/5)],ticks[-1]])
                    ax.set_yticks(ticks)
                    from matplotlib.ticker import FormatStrFormatter
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f' if args.panel_range_step>=.1 else '%.2g'))
                title=panel['gene']+(' · '+panel['state'] if panel.get('state') else '')
                ax.text(.5,1.68 if args.panel_range_step else 1.32,title,ha='center',va='bottom',transform=ax.transAxes,
                        fontsize=11.5 if panel.get('state') else 13,fontweight='bold')
                ax.text(.5,1.43 if args.panel_range_step else 1.13,event_lines,ha='center',va='bottom',transform=ax.transAxes,
                        fontsize=8.3,linespacing=1.2)
                ax.text(-.13,1.68 if args.panel_range_step else 1.32,chr(65+index),transform=ax.transAxes,ha='left',va='bottom',
                        fontsize=12,fontweight='bold')
                result.update(gene=panel['gene'],event_label=panel['event_label'],
                              state=panel.get('state'),
                              source_selection=panel.get('source_selection'))
                results.append(result)
            for ax in list(axes.flat)[len(panels):]:
                ax.set_visible(False)
            fig.suptitle(specification['title'],fontsize=18,fontweight='bold',y=.985)
            handles = [Line2D([],[],marker=m,linestyle='none',color='black',markersize=6,label=l)
                       for l,m in results[0]['marker_legend'].items()]
            fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.086),
                       ncol=max(1,len(handles)),frameon=False,fontsize=11)
            fig.text(.045,.018,textwrap.fill(specification['caption'],width=228),
                     fontsize=9,va='bottom',ha='left',linespacing=1.3)
            fig.canvas.draw()
            # Reject text clipped outside the page; densely laid-out panels must
            # keep their event IDs, donor counts and annotations reviewable.
            renderer=fig.canvas.get_renderer()
            for text in fig.findobj(matplotlib.text.Text):
                if not text.get_visible() or not text.get_text():
                    continue
                bounds=text.get_window_extent(renderer)
                if bounds.x0 < -1 or bounds.y0 < -1 or bounds.x1 > fig.bbox.x1+1 or bounds.y1 > fig.bbox.y1+1:
                    raise ValueError('Composite text extends outside page: '+text.get_text())
            out.parent.mkdir(parents=True,exist_ok=True)
            fig.savefig(out)
        finally:
            plt.close(fig)
    records=[dict(r,gene=p['gene'])
             for p in results for r in (dict(x) for x in p['records'])]
    # event_uid is already present in the native records.
    with outputs[2].open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(records[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(records)
    result=dict(status='COMPLETE_SINGLE_PAGE_NATIVE_PSI_COMPOSITE',panels=results,pages=1,
                required_genes=genes,source_manifest_sha256=source_hash,
                code_sha256=psi_boxplot.sha256(__file__),
                native_renderer_sha256=psi_boxplot.sha256(psi_boxplot.__file__),
                PDF_sha256=psi_boxplot.sha256(out),figure_inches=[args.width,args.height],
                panel_range_step=args.panel_range_step,
                separate_panel_PDFs_written=False,caption=specification['caption'])
    outputs[1].write_text(json.dumps(result,indent=2)+'\n')
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    return run(parser.parse_args(argv))


if __name__=='__main__':
    main()
