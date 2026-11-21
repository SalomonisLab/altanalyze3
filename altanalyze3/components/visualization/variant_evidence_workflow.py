"""Approved KINNEX evidence-ring workflow from the complete step-04 tables. Labels may be selected without removing any calls. LOY is a separately labeled short-read source."""
import sys, math, textwrap
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Polygon, Rectangle
from .figure_style import apply_rcparams, strip_masks
from altanalyze3.components.visualization.palettes import categorical_colors
from pathlib import Path
import json, hashlib

def render(root, output, evidence_set='A_tool_thresholds', top_labels=None, label_items=None, combine=True, include_loy=False):
    R = str(root)
    SET = evidence_set
    TOPN = top_labels
    CENTER = True
    FULL1 = None
    COMBINE = combine
    ADD_LOY = include_loy
    from ..bam.clonal_inference.detection_plotting import require_arial
    require_arial()
    LOYNAME = 'LOY (loss of chromosome Y)‡'
    V = pd.read_csv(R + '/tables/variant_summary_for_figure__%s.tsv' % SET, sep='\t')
    U = pd.read_csv(R + '/tables/variant_units_by_sample__%s.tsv' % SET, sep='\t')
    LONG = pd.read_csv(R + '/tables/variant_cells_long__%s.tsv.gz' % SET, sep='\t')
    EV = pd.read_csv(R + '/tables/evidence_by_unit.tsv', sep='\t')
    apply_rcparams(6)
    vcol = dict(zip(V.sort_values(['variant_class', 'cells_detected'], ascending=[True, False]).variant, categorical_colors(len(V))))
    N_READ, N_INF = (int(V.cells_detected.sum()), int(V.cells_inferred.sum()))
    LOYNOTE = ''
    if ADD_LOY:
        _l = pd.read_csv(R + '/tables/LOY_by_patient_for_figure.tsv', sep='\t')
        _m = _l[_l.inferred_sex == 'male']
        _row = {c: np.nan for c in V.columns}
        _row.update(variant_class='chromosome_loss', variant=LOYNAME, genes='chrY', n_patients=int((_m.LOY_cells > 0).sum()), n_samples=int(_m.libraries_with_LOY_cells.sum()), clonal=False, CH=False, cells_detected=int(_m.LOY_cells.sum()), cells_inferred=0, cells_detected_plus_inferred=int(_m.LOY_cells.sum()), **{c: 0 for c in V.columns if c.startswith('inferred_')})
        V = pd.concat([V, pd.DataFrame([_row])], ignore_index=True)
        for c in ['n_patients', 'n_samples', 'cells_detected', 'cells_inferred', 'cells_detected_plus_inferred']:
            V[c] = V[c].astype(int)
        V['clonal'] = V.clonal.astype(bool)
        V['CH'] = V.CH.astype(bool)
        vcol[LOYNAME] = '#5f6368'
        LOYNOTE = '‡ LOY = cells that fastCNV assigns to its chrY-loss clone from short-read expression; %s LOY cells of %s short-read cells in %d male patients (%d of %d with at least one LOY cell; no minimum LOY fraction applied; %d female patients and 1 patient without a prediction not counted). These are short-read cells, not the KINNEX cells of the other arcs; LOY is drawn in rings 1 and 2 and has no ring-3 evidence. ' % (format(int(_m.LOY_cells.sum()), ','), format(int(_m.cells_tested.sum()), ','), len(_m), int((_m.LOY_cells > 0).sum()), len(_m), int((_l.inferred_sex == 'female').sum()))
    SHELLS = [('PV', '#d1352b', 'nuclear passenger variant'), ('MPV', '#2f6fb7', 'mitochondrial passenger variant'), ('LOH', '#2e9a62', 'chr21 loss of heterozygosity'), ('CNV', '#e08a1e', 'inferCNVpy clone'), ('DEG', '#7b4ea3', 'gene-expression classifier'), ('DEI', '#1f9aa6', 'isoform-expression classifier')]
    CODES = [s[0] for s in SHELLS]
    known = set(CODES)
    combos_seen = set((c for e in LONG.evidence[~LONG.is_read].unique() for c in e.split('+')))
    if not combos_seen <= known:
        raise ValueError('evidence code without a shell: %s' % (combos_seen - known))
    R1I, R1O = (0.36, 0.52) if FULL1 else (0.4, 0.58)
    R2I, R2O = (0.76, 0.94)
    R3I, SW, SG = (0.99, 0.052, 0.012)
    R3O = R3I + len(SHELLS) * SW + (len(SHELLS) - 1) * SG
    YMAX, STEP = (3.15, 0.24)
    XMAXS = [3.93, 3.61, 3.22]

    def lighten(hexc, f):
        r, g, b = (int(hexc[i:i + 2], 16) for i in (1, 3, 5))
        return '#%02x%02x%02x' % tuple((int(c + (255 - c) * f) for c in (r, g, b)))

    def gray_order(combos, k=0, rev=True):
        """Reflected ordering of membership tuples: the cells of shell 1 form one run, and the cells of shell 2 form one
        run that spans the edge of the shell-1 run, so coincident calls of two methods sit radially above each other."""
        if k == len(CODES) or len(combos) <= 1:
            return list(combos)
        ones = [c for c in combos if c[k]]
        zeros = [c for c in combos if not c[k]]
        if rev:
            return gray_order(ones, k + 1, False) + gray_order(zeros, k + 1, True)
        return gray_order(zeros, k + 1, False) + gray_order(ones, k + 1, True)

    def place(labels):
        out = {}
        sides = {i: 1 if math.cos(math.radians(a)) >= 0 else -1 for i, (a, _) in enumerate(labels)}
        maxn = int((2 * YMAX - 0.3) / STEP)
        for sd in (1, -1):
            while sum((1 for v in sides.values() if v == sd)) > maxn:
                j = min((i for i, v in sides.items() if v == sd), key=lambda i: abs(math.cos(math.radians(labels[i][0]))))
                sides[j] = -sd
        for side in (1, -1):
            idx = sorted((i for i in sides if sides[i] == side), key=lambda i: math.sin(math.radians(labels[i][0])))
            n = len(idx)
            if not n:
                continue
            ys = [(R3O + 0.35) * math.sin(math.radians(labels[i][0])) for i in idx]
            for k in range(1, n):
                ys[k] = max(ys[k], ys[k - 1] + STEP)
            ys[-1] = min(ys[-1], YMAX - 0.15)
            for k in range(n - 2, -1, -1):
                ys[k] = min(ys[k], ys[k + 1] - STEP)
            if not ys[0] >= -YMAX:
                raise ValueError('too many labels on one side for the panel height')
            for i, y in zip(idx, ys):
                out[i] = (side, y)
        return out

    def xy(r, a):
        return (r * math.cos(math.radians(a)), r * math.sin(math.radians(a)))

    def wedge(ax, r0, r1, a_lo, a_hi, colour):
        thin = a_hi - a_lo < 0.8
        ax.add_patch(Wedge((0, 0), r1, a_lo, a_hi, width=r1 - r0, facecolor=colour, edgecolor=colour if thin else '#ffffff', linewidth=0.25))

    def rings(ax, d, L, XMAX):
        """d: variant rows of one class in drawing order; L: the cell-variant calls of that class."""
        tot = int(d.cells_detected_plus_inferred.sum())
        deg = 360.0 / tot
        ax.set_xlim(-XMAX, XMAX)
        ax.set_ylim(-YMAX, YMAX)
        ax.set_aspect('equal')
        ax.axis('off')
        small = sum((deg * r.cells_detected_plus_inferred for r in d.itertuples() if r.cells_detected_plus_inferred / tot < 0.03))
        START = 180.0 + small / 2.0 if len(d) > 3 else 90.0
        ang = START
        order = list(d.itertuples())
        order = [r for r in order if r.cells_detected_plus_inferred / tot < 0.03][::-1] + [r for r in order if r.cells_detected_plus_inferred / tot >= 0.03]
        labels = []
        shell_cells = {c: 0 for c in CODES}
        ring1 = {}
        if FULL1:
            rdtot = float(sum((r.cells_detected for r in order)))
            d1 = 360.0 / rdtot
            m2, m1, w1 = ([], [], [])
            a_ = START
            c_ = 0.0
            for r in order:
                T_ = deg * r.cells_detected_plus_inferred
                W_ = d1 * r.cells_detected
                m2.append(a_ - T_ / 2.0)
                m1.append(-(c_ + W_ / 2.0))
                w1.append(r.cells_detected)
                a_ -= T_
                c_ += W_
            m2 = np.array(m2)
            m1 = np.array(m1)
            w1 = np.array(w1, dtype=float)
            offs = np.arange(0.0, 360.0, 0.5)
            cost = [float((w1 * ((m1 + o - m2 + 180.0) % 360.0 - 180.0) ** 2).sum()) for o in offs]
            off1 = float(offs[int(np.argmin(cost))])
            c_ = 0.0
            for r, mm in zip(order, m2):
                W_ = d1 * r.cells_detected
                hi = off1 - c_
                k_ = round((mm - (hi - W_ / 2.0)) / 360.0)
                ring1[r.variant] = (hi + 360.0 * k_ - W_, hi + 360.0 * k_)
                c_ += W_
        for r in order:
            col = vcol[r.variant]
            T = deg * r.cells_detected_plus_inferred
            Rw = deg * r.cells_detected
            a1 = ang
            F = (T - Rw) / 2.0 if CENTER else 0.0
            h1 = a1 - F
            l1 = h1 - Rw
            if FULL1:
                l1, h1 = ring1[r.variant]
            wedge(ax, R1I, R1O, l1, h1, col)
            wedge(ax, R2I, R2O, a1 - T, a1, col)
            ts = np.linspace(0, 1, 24)
            sm = 3 * ts ** 2 - 2 * ts ** 3
            if FULL1 == 'full_sweep':
                sm = ts
            lead = [xy(R1O + (R2I - R1O) * t, h1 + (a1 - h1) * u) for t, u in zip(ts, sm)]
            trail = [xy(R1O + (R2I - R1O) * t, l1 + (a1 - T - l1) * u) for t, u in zip(ts, sm)]
            if FULL1 == 'full_straight':
                lead = [xy(R1O, h1), xy(R2I, a1)]
                trail = [xy(R1O, l1), xy(R2I, a1 - T)]
            inner = [xy(R1O, a) for a in np.linspace(h1, l1, max(2, int((h1 - l1) / 2) + 2))]
            outer = [xy(R2I, a) for a in np.linspace(a1 - T, a1, max(2, int(T / 2) + 2))]
            ax.add_patch(Polygon(inner + trail[1:-1] + outer + lead[::-1][1:-1], closed=True, facecolor=lighten(col, 0.84), edgecolor='none'))
            for edge, flat in ((lead, F < 1e-09 and (not FULL1)), (trail, T - Rw - F < 1e-09 and (not FULL1))):
                if flat or len(edge) == 2:
                    ax.plot(*zip(edge[0], edge[-1]), color=col, lw=0.35, solid_capstyle='butt')
                else:
                    ax.plot(*zip(*edge), color=col, lw=0.35, solid_capstyle='butt')
            inf = L[(L.variant == r.variant) & ~L.is_read]
            if not len(inf) == r.cells_inferred:
                raise ValueError('inferred calls %d != summary %d for %s' % (len(inf), r.cells_inferred, r.variant))
            cnt = inf.evidence.value_counts()
            memb = {tuple((c in e.split('+') for c in CODES)): n for e, n in cnt.items()}
            off = 0.0
            runs = {c: [] for c in CODES}
            for m in gray_order(list(memb)):
                w = deg * memb[m]
                for k, c in enumerate(CODES):
                    if m[k]:
                        shell_cells[c] += memb[m]
                        if runs[c] and abs(runs[c][-1][1] - off) < 1e-09:
                            runs[c][-1][1] = off + w
                        else:
                            runs[c].append([off, off + w])
                off += w
            if not abs(off - (T - Rw)) < 1e-06:
                raise ValueError('Integrity check failed: abs(off - (T - Rw)) < 1e-06')
            for k, (c, scol, _) in enumerate(SHELLS):
                r0 = R3I + k * (SW + SG)
                for o0, o1 in runs[c]:
                    pieces = []
                    if o0 < F:
                        pieces.append((a1 - min(o1, F), a1 - o0))
                    if o1 > F:
                        pieces.append((a1 - Rw - o1, a1 - Rw - max(o0, F)))
                    for lo, hi in pieces:
                        ax.add_patch(Wedge((0, 0), r0 + SW, lo, hi, width=SW, facecolor=scol, edgecolor='none'))
            labels.append((a1 - T / 2.0, r, r.cells_detected_plus_inferred))
            ang = a1 - T
        if not abs(ang + 360.0 - START) < 1e-06:
            raise ValueError('Integrity check failed: abs(ang + 360.0 - START) < 1e-06')
        for k in range(len(SHELLS) + 1):
            rr = R3I + k * (SW + SG) - SG / 2.0
            ax.add_patch(plt.Circle((0, 0), rr, facecolor='none', edgecolor='#d5d8dc', linewidth=0.2))
        if TOPN is not None and len(labels) > TOPN:
            keep = sorted(range(len(labels)), key=lambda i: -labels[i][2])[:TOPN]
            labels = [labels[i] for i in sorted(keep)]
        lab = [(a, r) for a, r, _ in labels if label_items is None or r.variant in label_items]
        pos = place(lab)
        for i, (a, r) in enumerate(lab):
            side, y = pos[i]
            x0, y0 = xy(R2O if CENTER else R3O + 0.02, a)
            x1 = side * (R3O + 0.42)
            lk = dict(color=vcol[r.variant], lw=1.2, alpha=0.5) if CENTER else dict(color='#9aa0a6', lw=0.3)
            ax.plot([x0, x1], [y0, y], solid_capstyle='butt', **lk)
            ax.plot([x1, side * (R3O + 0.52)], [y, y], solid_capstyle='butt', **lk)
            dag = '†' if isinstance(r.caller_note, str) and r.caller_note else ''
            ev = ' (' + r.evidence_types + ')' if isinstance(r.evidence_types, str) and r.evidence_types else ''
            t = 'patients = %d' % r.n_patients + (' (clonal)' if r.clonal else '') + (' (CH)' if r.CH else '')
            n = '%s cells' % format(r.cells_detected, ',') if r.cells_inferred == 0 else '%s → %s cells' % (format(r.cells_detected, ','), format(r.cells_detected_plus_inferred, ','))
            ax.add_patch(Rectangle((side * (R3O + 0.56) - (0.0 if side == 1 else 0.07), y + 0.025), 0.07, 0.07, facecolor=vcol[r.variant], edgecolor='none'))
            ax.text(side * (R3O + 0.67), y, '%s%s%s\n%s, %s' % (r.variant, dag, ev, t, n), ha='left' if side == 1 else 'right', va='center', fontsize=4.8, color='#202124', linespacing=1.08)
        rd = int(d.cells_detected.sum())
        for k, (line, fs) in enumerate([(format(rd, ','), 6.2), ('read-called', 4.8), (format(tot, ','), 6.2), ('with imputed', 4.8)]):
            ax.text(0, 0.21 - 0.14 * k, line, ha='center', va='center', fontsize=fs, color='#202124')
        return (shell_cells, rd, tot)
    CLASSES = [('substitution', 'Substitutions'), ('indel', 'Small insertions and deletions'), ('rearrangement', 'Rearrangements')]
    MEMBERS = {vc: [vc] for vc, _ in CLASSES}
    FIGW = 17.0
    if COMBINE:
        CLASSES = [('substitution', 'Substitutions'), ('indel_and_rearrangement', 'Small insertions and deletions, and rearrangements')]
        MEMBERS = {'substitution': ['substitution'], 'indel_and_rearrangement': ['indel', 'rearrangement']}
        if ADD_LOY:
            CLASSES[1] = ('indel_rearrangement_LOY', 'Small insertions and deletions, rearrangements, and loss of chromosome Y')
            MEMBERS = {'substitution': ['substitution'], 'indel_rearrangement_LOY': ['indel', 'rearrangement', 'chromosome_loss']}
        FIGW = 17.0 * (XMAXS[0] + XMAXS[0]) / sum(XMAXS)
        XMAXS = [XMAXS[0], XMAXS[0]]
    fig = plt.figure(figsize=(FIGW, 6.8))
    gs = fig.add_gridspec(1, len(CLASSES), left=0.005, right=0.995, top=0.965, bottom=0.225, wspace=0.0, width_ratios=XMAXS)
    stats = {}
    for c, (vc, title) in enumerate(CLASSES):
        ax = fig.add_subplot(gs[0, c])
        d = V[V.variant_class.isin(MEMBERS[vc])].sort_values('cells_detected_plus_inferred', ascending=False)
        stats[vc] = rings(ax, d, LONG[LONG.variant_class.isin(MEMBERS[vc])], XMAXS[c])
        ax.set_title('%s   %s (%d variants)' % ('abc'[c], title, len(d)), loc='left', fontsize=7, fontweight='bold', pad=2)
    kax = fig.add_axes([0.02, 0.172, 0.96, 0.05])
    kax.set_xlim(0, 100)
    kax.set_ylim(0, 2)
    kax.axis('off')
    kax.text(0, 1.5, 'Ring 1 (inner): cells called from reads.   Ring 2 (middle): read-called plus imputed cells; a band joins the two arcs of a variant.   Ring 3 (outer): source of the imputed cells, one shell per method, white = not inferred by that method.   a → b cells = read-called → with imputed.', ha='left', va='center', fontsize=5.8, color='#202124')
    kax.text(0, 0.5, 'Shells, inside to outside:' if COMBINE else 'Ring 3 shells, inside to outside:', ha='left', va='center', fontsize=5.8, color='#202124')
    KS = 17.0 / FIGW
    x = 13.5 if COMBINE else 12.0
    alls = {c: sum((stats[vc][0][c] for vc in stats)) for c in CODES}
    for i, (c, scol, name) in enumerate(SHELLS):
        kax.add_patch(Rectangle((x, 0.25), 0.7 * KS, 0.5, facecolor=scol, edgecolor='none'))
        t = '%d  %s  %s%s' % (i + 1, c, name, ' (not computed yet)' if c == 'DEI' and (not (EV.evidence == 'DEI').any()) else '')
        kax.text(x + 0.95 * KS, 0.5, t, ha='left', va='center', fontsize=5.6, color='#202124')
        x += ((1.6 if COMBINE else 2.2) + 0.262 * len(t)) * KS
    if not x < 101:
        raise ValueError('shell key overflows the figure width')
    rd_all = LONG[LONG.is_read]
    inf_all = LONG[~LONG.is_read]
    nA = int((EV.evidence != 'MPV_weak_tier').sum())
    nT = int(EV.testable.sum())
    nB = int(EV.in_set_B_association_reproduced.sum())
    rule = 'A surrogate enters rings 2 and 3 only if its association with the read genotype holds in this object: one-sided Fisher exact test on carrier by MUT/WT among read-genotyped cells, Benjamini-Hochberg q < 0.05 over %d testable sample-variant-surrogate tests; %d of %d surrogates pass.' % (nT, nB, nA) if SET.startswith('B') else 'Rings 2 and 3 keep every surrogate at the threshold of the tool that nominated it (%d sample-variant-surrogate rows); no test inside this object is applied.' % nA
    deg = EV[(EV.evidence == 'DEG') & EV.testable]
    legend = '\n'.join(['Known MDS/AML variants called from KINNEX long reads mark %s cell-variant calls in 265,783 annotated cells, and surrogate evidence adds %s inferred calls (%s).' % (format(len(rd_all), ','), format(len(inf_all), ','), 'association reproduced' if SET.startswith('B') else 'tool thresholds'), 'Each panel: one colour per variant. Ring 1 = cells with >= 1 read carrying the variant allele (MUT), summed over the samples in which the variant is called. Ring 2 = those cells plus cells without a covering read (UNK) that carry an associated surrogate (imputed). One degree is the same number of cell-variant calls in every ring of a panel, and the full circle is the ring-2 total, so ring 1 leaves a gap where a variant gains imputed cells; the pale band joins the two arcs of a variant. Ring 3 = six shells over the imputed cells only, coloured where that method supports the call and white where it does not; colour in two shells at one angle = the same cells inferred by both methods. Imputed cells are predictions, not measurements. Imputed calls per method, all panels (a call supported by two methods counts in both): ' + ', '.join(('%s %s' % (c, format(alls[c], ',')) for c in CODES)) + '.', "Calls: %d sample-variant units, %d variants; substitutions and indels by the cohort rule of code/aggregate.py (alt reads >= 3, VAF >= max(0.05, healthy-panel VAF + 0.05, 3 x healthy-panel VAF); 9 healthy donors as the background panel); rearrangements in the documented lesion samples; CH = hotspot call in a healthy donor by a leave-one-out rule (VAF >= 0.02). %s (clonal) = within one patient the share of the variant among all variant reads of a sample changes >= 4-fold between samples (wild-type cells ignored; >= 50 variant reads in one sample; samples with >= 100 variant reads). Inherited variants are not shown: a variant in >= 95 %% of reads of the patient that is also in >= 95 %% of reads of mature lymphocytes, or of a sample in which the patient's other variants are below 5 %%, is inherited (TET2 c.86C>G and GATA1 c.212A>G in one patient). patients = patients with the variant in any sample; a xenograft, culture, therapy or transformation sample is a sample of its patient." % (len(U), len(V) - int(ADD_LOY), rule), ('Callouts name the %d most frequent variants of each panel; every variant is in tables/variant_summary_for_figure__%s.tsv. ' % (TOPN, SET) if TOPN else '') + 'Scope: 39 samples, 50 captures, cells restricted to the annotated cells of isoform_final_20261016. † the caller cannot validate this panel entry (no indel class in the label). TP53 R248R is a synonymous change. DEG: %d of %d testable units reproduce in this object; DEI: %s. One read per cell is typical, so a WT call can be allelic dropout and inferred counts for heterozygous variants are upper bounds on confirmed mutant cells.' % (int(deg.in_set_B_association_reproduced.sum()), len(deg), 'isoform classifier not finished, shell 6 is empty' if not (EV.evidence == 'DEI').any() else '%d units' % int((EV.evidence == 'DEI').sum()))])
    if LOYNOTE:
        legend += '\n' + LOYNOTE + 'The centre totals of panel b include the LOY cells.'
    if CENTER:
        a_ = 'so ring 1 leaves a gap where a variant gains imputed cells; '
        if not a_ in legend:
            raise ValueError('Integrity check failed: a_ in legend')
        legend = legend.replace(a_, 'so ring 1 leaves a gap where a variant gains imputed cells; the ring-1 arc is centred under its ring-2 arc and the imputed cells flank it on both sides; ')
    if FULL1:
        i_ = legend.index('One degree is the same number')
        j_ = legend.index('the pale band joins')
        legend = legend[:i_] + 'Ring 1 and ring 2 are each a full circle on their own scale: an arc of ring 1 is the share of the variant among the read-called cells of the panel, an arc of ring 2 its share among read-called plus imputed cells, so the width of a band changes with the share and not with the cell count; ' + legend[j_:]
    fig.text(0.02, 0.162, '\n'.join((textwrap.fill(p, int(400 / KS * (1.17 if ADD_LOY else 1.0))) for p in legend.split('\n'))), ha='left', va='top', fontsize=5.6, color='#202124', linespacing=1.25)
    strip_masks(fig)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    for ext in ['svg', 'pdf', 'png']:
        fig.savefig(out / ('variant_evidence_rings.' + ext), dpi=300)
    plt.close(fig)
    (out / 'figure_legend.txt').write_text(legend + '\n')
    counts = pd.DataFrame([dict(variant_class=vc, read_called=stats[vc][1], read_called_plus_inferred=stats[vc][2], **{'shell_' + c: stats[vc][0][c] for c in CODES}) for vc in stats])
    counts.to_csv(out / 'ring_figure_counts.tsv', sep='\t', index=False)
    for vc in stats:
        d = V[V.variant_class.isin(MEMBERS[vc])]
        if not (int(d.cells_detected.sum()) == stats[vc][1] and int(d.cells_detected_plus_inferred.sum()) == stats[vc][2]):
            raise ValueError('Integrity check failed: int(d.cells_detected.sum()) == stats[vc][1] and int(d.cells_detected_plus_inferred.sum()) == stats[vc][2]')
        for c in CODES:
            if not int(d['inferred_' + c].sum()) == stats[vc][0][c]:
                raise ValueError("Integrity check failed: int(d['inferred_' + c].sum()) == stats[vc][0][c]")
    source_files = [Path(R) / 'tables' / ('variant_summary_for_figure__' + SET + '.tsv'), Path(R) / 'tables' / ('variant_units_by_sample__' + SET + '.tsv'), Path(R) / 'tables' / ('variant_cells_long__' + SET + '.tsv.gz'), Path(R) / 'tables/evidence_by_unit.tsv']
    if ADD_LOY:
        source_files.append(Path(R) / 'tables/LOY_by_patient_for_figure.tsv')
    report = {'set': SET, 'label_items': label_items, 'all_variants_retained': True, 'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}, 'drawn_counts_equal_source': True, 'include_LOY': ADD_LOY}
    (out / 'validation.json').write_text(json.dumps(report, indent=2))
    return report

def add_arguments(p):
    p.add_argument('--root', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--evidence-set', default='A_tool_thresholds')
    p.add_argument('--top-labels', type=int)
    p.add_argument('--label-items', nargs='+')
    p.add_argument('--separate-classes', action='store_true')
    p.add_argument('--include-loy', action='store_true')

def run(a):
    return render(a.root, a.output, a.evidence_set, a.top_labels, a.label_items, not a.separate_classes, a.include_loy)

def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    add_arguments(p)
    return run(p.parse_args(argv))
if __name__ == '__main__':
    main()
