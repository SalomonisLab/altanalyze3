"""Catalog flow summary: linear versus chimeric reads, per-source iso-structures, collapsed units,
UMI>=3 classes and protein-coding calls, written as standard text tables and the flow diagram.

Definitions (Nathan Salomonis, 2026-11-11):
  UMI total          native-eligible transcript_associations rows of every source, plus chimeric
                     reads that only the split-alignment detectors see (each read counted once).
  chimeric reads     circRNA reads (back-splice joints whose joint_class is in
                     catalog_flow.circrna_classes, status kept), inter-chromosomal trans-splicing reads
                     (classification class trans_splicing_candidate, status kept) and reads assigned to
                     intra-chromosomal trans-splicing catalog isoforms (distinct-locus spliced
                     candidates of the native geometry classifier, novel, >= min_total reads).
  linear reads       every other native-eligible read.
  Iso-Structures     per-source distinct (gene, exon string) keys, summed over sources (the native
                     per-source structure definition). Chimeric: intra-trans-splicing read structures
                     plus per-source circRNA and inter-chromosomal junctions.
  Collapsed          linear: collapsed isoforms before the support filter (the native before_min_total
                     count, which includes zero-read reference exemplars) minus those whose reads are all
                     chimeric (FINAL_structure_to_exemplar_prefilter.tsv.gz); chimeric: intra-trans-splicing
                     isoforms plus junctions with >= 1 kept read.
  UMI>=3             >= min_total reads (linear reads for Known/Novel).
  Protein coding     NMD/Potential-NMD -> likely NMD; length < 25% of the gene's longest protein ->
                     truncated; otherwise protein_coding (the native summary rule).

Oncofusion/PTD/ITD events are reported in REARRANGEMENT_EVENTS.tsv and are not part of the flow.
Outputs: <output>/integrated/summary/.
"""
from __future__ import annotations
import collections
import csv
import gzip
import json
import sqlite3
import sys
from pathlib import Path
from .manifest import atomic_json, receipt, verify_receipt, write_tsv

csv.field_size_limit(sys.maxsize)
DEFAULT_CIRC_CLASSES = ('backsplice_registered_canonical',)
# combined.gff holds GTF-style novel rows (transcript_id "X") and GFF3-style GENCODE rows (transcript_id=X).
TRANSCRIPT_ID = __import__('re').compile(r'transcript_id(?:=| ")([^";]+)')
TRUNCATED_FRACTION = 0.25


def _open(path):
    path = Path(path)
    if path.suffix == '.gz':
        return gzip.open(path, 'rt', newline='')
    if not path.exists() and Path(str(path) + '.gz').exists():
        return gzip.open(str(path) + '.gz', 'rt', newline='')
    return open(path, newline='')


def _core(name):
    name = str(name).strip()
    return name.split('molecule/', 1)[1] if name.startswith('molecule/') else name


def eligible(row):
    """Native source-scan eligibility (source_stream.scan_sample): exon string present, no UNK gene."""
    return bool(row[2]) and 'UNK' not in row[0]


def protein_classes(path, fraction=TRUNCATED_FRACTION):
    """(gene, isoform) -> class, and gene -> longest protein, from the native protein_summary table."""
    classes, longest = {}, {}
    with _open(path) as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)
        if header[:4] != ['Gene ID', 'Transcript ID', 'Protein Length', 'NMD Status']:
            raise ValueError(f'{path}: unexpected protein_summary header {header[:6]}')
        for p in reader:
            gene, iso, nmd = p[0].split(';', 1)[0], p[1], p[3].strip()
            try:
                length, top = int(p[2]), int(p[5])
            except ValueError:
                classes[(gene, iso)] = 'unclassified'
                continue
            if nmd != 'Not-NMD' and 'NMD' in nmd:
                cls = 'likely_NMD'
            elif top > 0 and length < fraction * top:
                cls = 'truncated_lt_25pct_longest'
            else:
                cls = 'protein_coding'
            if (gene, iso) in classes and classes[(gene, iso)] != cls:
                raise ValueError(f'Conflicting protein class rows for {gene} {iso}')
            classes[(gene, iso)] = cls
            if nmd == 'Not-NMD':
                if longest.get(gene, top) != top:
                    raise ValueError(f'Gene longest protein differs between Not-NMD rows: {gene}')
                longest[gene] = top
    return classes, longest


def read_catalog(path):
    catalog = {}
    with _open(path) as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            key = (r['gene'], r['final_isoform_id'])
            if key in catalog:
                raise ValueError(f'Duplicate catalog isoform {key}')
            catalog[key] = (r['known'], int(r['total_reads']))
    return catalog


def catalog_geometry(gff_output, exon_model, reference_gff):
    """Distinct-locus geometry of every catalog isoform whose structure names a second gene."""
    from altanalyze3.components.annotation import trans_splicing as TS
    gff_output = Path(gff_output)
    structures = {}
    with _open(gff_output / 'transcript_associations.txt.gz') as handle:
        for p in csv.reader(handle, delimiter='\t'):
            gene, structure, iso = p[0], p[2], p[3]
            if any(':' in t and TS.token_parts(t, gene)[0] != TS.gene_id(gene) for t in structure.split('|') if t):
                structures[(gene, iso)] = structure
    exons = collections.defaultdict(list)
    wanted = {iso for _, iso in structures}
    with _open(gff_output / 'combined.gff.gz') as handle:
        for line in handle:
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'exon':
                continue
            match = TRANSCRIPT_ID.search(p[8])
            if match is None:
                raise ValueError(f'combined.gff exon row without transcript_id: {line[:200]}')
            iso = match.group(1)
            if iso in wanted:
                exons[iso].append((p[0], p[6], int(p[3]), int(p[4])))
    reference = TS.reference(exon_model, reference_gff)
    status = {}
    for (gene, iso), structure in structures.items():
        if iso not in exons:
            status[(gene, iso)] = 'no_exon_rows_in_combined_gff'
            continue
        trans = TS.characterize(gene, structure, exons[iso], reference)[2]
        status[(gene, iso)] = trans['status']
    return status


def circrna_reads(evidence_dir, classes):
    """{(source, molecule)} and per-event kept reads for back-splice joints in the requested classes."""
    evidence_dir = Path(evidence_dir)
    event_class = {}
    with _open(evidence_dir / 'classes.tsv.gz') as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            if r['event'].startswith('BSJ|'):
                previous = event_class.setdefault(r['event'], r['joint_class'])
                if previous != r['joint_class']:
                    raise ValueError(f"Event with more than one joint_class: {r['event']}")
    reads, per_event, per_source = set(), collections.defaultdict(set), collections.defaultdict(set)
    with _open(evidence_dir / 'read_support.tsv.gz') as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            e = r['event']
            if not e.startswith('BSJ|') or r['status'] != 'kept' or event_class[e] not in classes:
                continue
            key = (r['library'], _core(r['molecule_id']))
            reads.add(key); per_event[e].add(key); per_source[r['library']].add(e)
    return reads, per_event, per_source


def interchromosomal_reads(db, classification):
    """Kept reads of trans_splicing_candidate junctions, keyed to the classification rows."""
    rows = {}
    with _open(classification) as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            if r['class'] == 'trans_splicing_candidate':
                rows[(r['gene_a'], r['gene_b'], r['joint'])] = r
    reads, per_junction, per_source = set(), collections.defaultdict(set), collections.defaultdict(set)
    query = '''SELECT s.source_id, s.query_name, s.status, s.evidence_json, r.molecule_id
               FROM support s JOIN events e USING(event_id) JOIN reads r USING(source_id, query_name)
               WHERE e.modality='trans-splicing' AND s.label='trans_splicing_candidate' '''
    any_status = collections.defaultdict(set)
    for sid, q, status, evidence, molecule in db.execute(query):
        r = json.loads(evidence)
        ga = 'UNK' if r['gene_a'].startswith('UNK') else r['gene_a']
        gb = 'UNK' if r['gene_b'].startswith('UNK') else r['gene_b']
        key = (ga, gb, f"{r['a_chrom']}:{r['donor']}:{r['a_strand']} > {r['b_chrom']}:{r['acceptor']}:{r['b_strand']}")
        if key not in rows:
            raise ValueError(f'Ledger trans_splicing_candidate read without a classification row: {key}')
        any_status[key].add((sid, q))
        if status == 'kept':
            k = (sid, _core(molecule))
            reads.add(k); per_junction[key].add(k); per_source[sid].add(key)
    for key, r in rows.items():
        if len(any_status.get(key, ())) != int(r['molecules']) or len(per_junction.get(key, ())) != int(r['kept_molecules']):
            raise ValueError(f'Ledger read counts differ from classification for {key}: '
                             f"{len(any_status.get(key, ()))}/{r['molecules']} molecules, "
                             f"{len(per_junction.get(key, ()))}/{r['kept_molecules']} kept")
    return reads, per_junction, per_source, rows


def exemplar_protein(row, longest, fraction=TRUNCATED_FRACTION):
    """Protein class of the junction's most-supported catalog exemplar (classification catalog_exemplars)."""
    text = (row.get('catalog_exemplars') or '').strip()
    if not text:
        return 'no_catalog_isoform'
    first = text.split(';')[0]
    inside = first[first.index('(') + 1:first.index(')')]
    aa, nmd = [x.strip() for x in inside.split(',', 1)]
    if nmd != 'Not-NMD':
        return 'likely_NMD'
    length = int(aa.split()[0])
    top = longest.get(row['gene_a'], 0)
    return 'truncated_lt_25pct_longest' if top > 0 and length < fraction * top else 'protein_coding'


def trace_source(sid, ta_path, mol2struct, fmap, pmap, ts_isoforms, chimeric_split):
    """Per-read trace of one source. Returns per-isoform totals and per-source structure flags."""
    fin = collections.defaultdict(lambda: [0, 0]); pre = collections.defaultdict(lambda: [0, 0])
    ts_reads = set(); n = 0; split_seen = 0; ts_and_split = 0
    with open(mol2struct) as handle:
        header = handle.readline().rstrip('\n').split('\t')
        if header != ['molecule_id', 'gene', 'structure']:
            raise ValueError(f'{mol2struct}: unexpected header {header}')
        for line in handle:
            mol, gene, structure = line.rstrip('\n').split('\t')
            n += 1
            fi, pi = fmap.get((gene, structure)), pmap.get((gene, structure))
            if pi is None:
                raise ValueError(f'{sid}: structure class absent from the prefilter map: {gene} {structure}')
            is_ts = fi is not None and (gene, fi) in ts_isoforms
            is_split = mol in chimeric_split
            split_seen += is_split
            chim = is_ts or is_split
            if is_ts:
                ts_reads.add(mol)
                ts_and_split += is_split
            if fi is not None:
                a = fin[(gene, fi)]; a[0] += 1; a[1] += not chim
            b = pre[(gene, pi)]; b[0] += 1; b[1] += not chim
    structures = {}; elig = 0
    with _open(ta_path) as handle:
        for p in csv.reader(handle, delimiter='\t'):
            if not eligible(p):
                continue
            elig += 1
            key = (p[0], p[2]); m = p[3]
            flags = structures.setdefault(key, [False, False])
            flags[0] |= not (m in chimeric_split or m in ts_reads)
            flags[1] |= m in ts_reads
    if elig != n:
        raise ValueError(f'{sid}: {n} mol2struct reads differ from {elig} eligible annotation rows')
    if split_seen != len(chimeric_split):
        raise ValueError(f'{sid}: {len(chimeric_split) - split_seen} chimeric eligible reads absent from mol2struct')
    return dict(fin=fin, pre=pre, ts_reads=len(ts_reads), ts_and_split=ts_and_split, eligible=elig,
                linear_structures=sum(1 for f in structures.values() if f[0]),
                intra_structures=sum(1 for f in structures.values() if f[1]),
                structures=len(structures))


def compute(config):
    root = Path(config['output'])
    params = dict(config.get('catalog_flow') or {})
    circ_classes = tuple(params.get('circrna_classes', DEFAULT_CIRC_CLASSES))
    min_total = int(config['standard_parameters']['min_total'])
    gff = root / 'compatibility/gff-output'
    catalog = read_catalog(gff / 'FINAL_isoform_catalog.tsv.gz')
    pclass, longest = protein_classes(gff / 'protein_summary.txt.gz')
    missing_class = [k for k in catalog if k not in pclass]
    geometry = catalog_geometry(gff, config['exon_model'], config['reference_gff'])
    ts = {k for k, s in geometry.items()
          if s == 'distinct_locus_spliced_candidate' and k in catalog and catalog[k][0] == 'novel' and catalog[k][1] >= min_total}
    known_ts = {k for k, s in geometry.items() if s == 'distinct_locus_spliced_candidate' and k in catalog and catalog[k][0] == 'known'}

    circ_reads, circ_events, circ_by_source = circrna_reads(root / 'circrna/evidence', circ_classes)
    db = sqlite3.connect(f"file:{root / 'reconciliation/evidence.sqlite'}?mode=ro", uri=True)
    inter_reads, inter_events, inter_by_source, inter_rows = interchromosomal_reads(
        db, root / 'trans_splicing/classification/joint_classification.tsv.gz')
    eligible_reads = {(sid, _core(m)) for sid, m in db.execute(
        "SELECT source_id, molecule_id FROM reads WHERE legacy_status='native_eligible'")}
    split = circ_reads | inter_reads
    split_in = split & eligible_reads

    fmap, pmap, pre_total = {}, {}, {}
    for path, target, col in ((gff / 'FINAL_structure_to_exemplar.tsv.gz', fmap, 'final_isoform_id'),
                              (gff / 'FINAL_structure_to_exemplar_prefilter.tsv.gz', pmap, 'prefilter_isoform_id')):
        if not path.exists():
            raise FileNotFoundError(f'{path} is required (the collapse writes it); rerun the standard stage')
        with _open(path) as handle:
            for r in csv.DictReader(handle, delimiter='\t'):
                target[(r['gene'], r['structure'])] = r[col]
                if col == 'prefilter_isoform_id':
                    pre_total[(r['gene'], r[col])] = int(r['isoform_total_reads'])
    fin = collections.defaultdict(lambda: [0, 0]); pre = collections.defaultdict(lambda: [0, 0])
    totals = collections.Counter()
    for source in config['sources']:
        sid = source['source_id']
        res = trace_source(sid, root / 'compatibility/samples' / sid / 'gff-output/transcript_associations.txt.gz',
                           gff / 'tier1' / f'{sid}.mol2struct.tsv', fmap, pmap, ts,
                           {m for s, m in split_in if s == sid})
        for k, v in res['fin'].items():
            fin[k][0] += v[0]; fin[k][1] += v[1]
        for k, v in res['pre'].items():
            pre[k][0] += v[0]; pre[k][1] += v[1]
        for k in ('ts_reads', 'ts_and_split', 'eligible', 'linear_structures', 'intra_structures', 'structures'):
            totals[k] += res[k]
    if len(eligible_reads) != totals['eligible']:
        raise ValueError(f"Ledger native-eligible reads ({len(eligible_reads)}) differ from annotation rows ({totals['eligible']})")
    bad = [k for k in catalog if fin[k][0] != catalog[k][1]]
    if bad or set(fin) - set(catalog):
        raise ValueError(f'Traced reads differ from catalog total_reads for {len(bad)} isoforms, e.g. {bad[:3]}')
    # Collapsed units before the support filter include reference (ENST) exemplars with zero reads; the
    # trace never reaches them. Every traced unit must carry exactly its map total.
    untraced = set(pre_total) - set(pre)
    if any(pre_total[k] != 0 for k in untraced) or set(pre) - set(pre_total):
        raise ValueError('Prefilter units without traced reads must have zero reads in the map')
    bad_pre = [k for k, v in pre.items() if v[0] != pre_total[k]]
    if bad_pre:
        raise ValueError(f'Traced reads differ from prefilter totals for {len(bad_pre)} units, e.g. {bad_pre[:3]}')

    c = collections.Counter()
    for k, (known, total) in catalog.items():
        linear = fin[k][1]
        if k in ts:
            continue
        group = 'known' if known == 'known' else 'novel'
        if linear >= min_total:
            c[group] += 1
            c[group + '_protein_coding'] += pclass.get(k) == 'protein_coding'
        elif total >= min_total:
            c[group + '_below_support_after_chimeric_removal'] += 1
    inter_pc = collections.Counter(exemplar_protein(inter_rows[k], longest) for k, v in inter_events.items() if len(v) >= min_total)
    values = collections.OrderedDict()
    values['umi_label'] = params.get('umi_label', f"UMI in n={len(config['sources'])}")
    values['umi_total'] = totals['eligible'] + len(split - eligible_reads)
    values['linear_structures'] = totals['linear_structures']
    lost_units = sum(1 for v in pre.values() if v[0] > 0 and v[1] == 0)
    values['linear_collapsed'] = len(pre_total) - lost_units
    values['known'] = c['known']; values['novel'] = c['novel']
    values['known_protein_coding'] = c['known_protein_coding']; values['novel_protein_coding'] = c['novel_protein_coding']
    values['chimeric_structures'] = (totals['intra_structures'] + sum(len(v) for v in circ_by_source.values())
                                     + sum(len(v) for v in inter_by_source.values()))
    values['chimeric_collapsed'] = len(ts) + len(circ_events) + len(inter_events)
    values['intra_trans_splice'] = len(ts)
    values['inter_trans_splice'] = sum(1 for v in inter_events.values() if len(v) >= min_total)
    values['circrna'] = sum(1 for v in circ_events.values() if len(v) >= min_total)
    values['intra_trans_splice_protein_coding'] = sum(1 for k in ts if pclass.get(k) == 'protein_coding')
    values['inter_trans_splice_protein_coding'] = inter_pc['protein_coding']
    values['support_label'] = f'UMI≥{min_total}'
    detail = dict(circrna_classes=list(circ_classes), min_total=min_total,
                  eligible_reads=totals['eligible'], split_alignment_reads=len(split),
                  split_alignment_reads_in_eligible=len(split_in), split_alignment_reads_new=len(split - eligible_reads),
                  circrna_reads=len(circ_reads), interchromosomal_reads=len(inter_reads),
                  intra_trans_splice_reads=totals['ts_reads'], intra_trans_splice_and_split_reads=totals['ts_and_split'],
                  chimeric_reads=len(split) + totals['ts_reads'] - totals['ts_and_split'],
                  structures_all=totals['structures'],
                  intra_trans_splice_structures=totals['intra_structures'],
                  circrna_source_junctions=sum(len(v) for v in circ_by_source.values()),
                  interchromosomal_source_junctions=sum(len(v) for v in inter_by_source.values()),
                  circrna_junctions=len(circ_events), interchromosomal_junctions=len(inter_events),
                  interchromosomal_candidates_in_classification=len(inter_rows),
                  interchromosomal_protein_classes_umi3=dict(inter_pc),
                  prefilter_isoforms=len(pre_total), prefilter_isoforms_zero_reads=len(untraced),
                  prefilter_isoforms_losing_all_reads_to_chimeric=lost_units,
                  known_distinct_locus_isoforms_kept_in_known=len(known_ts),
                  catalog_isoforms_without_protein_class=len(missing_class),
                  geometry_status_counts=dict(collections.Counter(geometry.values())),
                  known_below_support_after_chimeric_removal=c['known_below_support_after_chimeric_removal'],
                  novel_below_support_after_chimeric_removal=c['novel_below_support_after_chimeric_removal'])
    detail['linear_reads'] = values['umi_total'] - detail['chimeric_reads']
    if detail['linear_reads'] != totals['eligible'] - (len(split_in) + totals['ts_reads'] - totals['ts_and_split']):
        raise AssertionError('Linear/chimeric read partition does not add up to the UMI total')
    events = []
    for e, v in circ_events.items():
        events.append(dict(category='circRNA', event=e, unit='back-splice junction', kept_reads=len(v),
                           sources=len({s for s, _ in v}), umi_ge_min=int(len(v) >= min_total), protein_class='not_assessed'))
    for k, v in inter_events.items():
        events.append(dict(category='inter_trans_splice', event=f'{k[0]}|{k[1]}|{k[2]}', unit='inter-chromosomal junction',
                           kept_reads=len(v), sources=len({s for s, _ in v}), umi_ge_min=int(len(v) >= min_total),
                           protein_class=exemplar_protein(inter_rows[k], longest)))
    for k in sorted(ts):
        events.append(dict(category='intra_trans_splice', event=f'{k[0]}|{k[1]}', unit='catalog isoform',
                           kept_reads=catalog[k][1], sources='', umi_ge_min=1, protein_class=pclass.get(k, 'unclassified')))
    return values, detail, events, db


def rearrangement_events(db):
    """Every oncofusion, PTD and ITD event of the ledger with distinct supporting reads per source."""
    rows = []
    query = '''SELECT e.modality, e.event_id, e.geometry_json, s.label, s.source_id, COUNT(DISTINCT s.query_name)
               FROM events e JOIN support s USING(event_id) WHERE e.modality IN ('oncofusion','ptd','itd')
               GROUP BY e.modality, e.event_id, s.label, s.source_id ORDER BY e.modality, e.event_id, s.source_id'''
    for modality, eid, geometry, label, sid, n in db.execute(query):
        rows.append(dict(modality=modality, event_id=eid, geometry=geometry, label=label, source_id=sid, reads=n))
    return rows


def run(config):
    root = Path(config['output']); out = root / 'integrated/summary'; out.mkdir(parents=True, exist_ok=True)
    done = out / 'COMPLETE.json'
    if done.exists():
        return verify_receipt(done, config)
    for proof in ('compatibility/STANDARD_COMPLETE.json', 'reconciliation/COMPLETE.json',
                  'circrna/AGGREGATE_COMPLETE.json', 'trans_splicing/CLASSIFICATION_COMPLETE.json'):
        verify_receipt(root / proof, config)
    values, detail, events, db = compute(config)
    write_tsv(out / 'CATALOG_FLOW_VALUES.tsv', ['box', 'value'], ({'box': k, 'value': v} for k, v in values.items()))
    write_tsv(out / 'CHIMERIC_EVENTS.tsv', ['category', 'event', 'unit', 'kept_reads', 'sources', 'umi_ge_min', 'protein_class'], events)
    rearr = rearrangement_events(db)
    write_tsv(out / 'REARRANGEMENT_EVENTS.tsv', ['modality', 'event_id', 'geometry', 'label', 'source_id', 'reads'], rearr)
    atomic_json(out / 'CATALOG_FLOW.json', dict(values=values, detail=detail,
                definitions=__doc__, oncofusion_ptd_itd='reported in REARRANGEMENT_EVENTS.tsv; not part of the flow'))
    from altanalyze3.components.visualization.catalog_flow import render
    figure = render(dict(values), out / 'catalog_flow')
    files = [out / n for n in ('CATALOG_FLOW_VALUES.tsv', 'CHIMERIC_EVENTS.tsv', 'REARRANGEMENT_EVENTS.tsv', 'CATALOG_FLOW.json')]
    files += [Path(figure['pdf']), Path(figure['png'])]
    return receipt(done, config, {}, files, status='COMPLETE')


if __name__ == '__main__':
    from .manifest import load_config
    run(load_config(sys.argv[1]))
