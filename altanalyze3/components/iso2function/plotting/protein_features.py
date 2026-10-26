"""Standard vector protein-feature figures from completed iso2function impact audits.

No gene-specific rules, manual coordinates, network requests, functional assays or
new domain predictions. All source features and all selected native pairs remain
in exported audits. Geometry uses the existing amino-acid projection method.
"""
import hashlib
import itertools
import json
import math
import re
import textwrap
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO

from .. import config
from ..protein_impact import dump, sha, table, write
from ..protein_isoform_audit import align_mapping, project_feature

SCHEMA = 'iso2function-protein-feature-figures-v3'
DEFAULT_TYPES = ('Domain', 'DNA binding', 'Zinc finger', 'Motif', 'Coiled coil', 'Region')
PRIMARY_SOURCES = ('reviewed_uniprot_snapshot', 'literature_construct_snapshot')
# Essential feature classes (transmembrane helices, signal peptides, kinase and phosphatase
# domains, transcription-factor DNA-binding domains, and the other families listed in the
# whitelist) are ALWAYS displayed, whatever --annotation-source and --feature-type select,
# are labelled by their class (TM1, TM2, Kinase, DNA-binding), take the lane nearest the
# protein backbone and are drawn above any annotation they overlap. Rule of record set by
# Nathan Salomonis on 2026-10-23 after a CLDND1 panel showed Pfam family hits and none of the
# four UniProt transmembrane helices. The rules live in a data file beside this module so
# they can be reviewed and extended without touching code.
ESSENTIAL_WHITELIST_PATH = Path(__file__).with_name('essential_feature_whitelist.tsv')
ESSENTIAL_RULE_COLUMNS = ('class', 'label', 'priority', 'color', 'provider', 'name_pattern', 'accession_pattern', 'description_pattern', 'note')


def load_essential_whitelist(path=None):
    """Read the essential-feature rules; every regex must compile and every row must be complete."""
    path = Path(path) if path else ESSENTIAL_WHITELIST_PATH
    rows = table(path)
    if not rows:
        raise ValueError('Essential feature whitelist is empty: ' + str(path))
    rules = []
    for i, r in enumerate(rows, 2):
        if set(ESSENTIAL_RULE_COLUMNS) - set(r):
            raise ValueError('Essential whitelist row ' + str(i) + ' lacks columns: ' + repr(sorted(set(ESSENTIAL_RULE_COLUMNS) - set(r))))
        if not (r['class'] and r['label'] and r['provider'] and re.fullmatch(r'#[0-9A-Fa-f]{6}', r['color'] or '')):
            raise ValueError('Essential whitelist row ' + str(i) + ' needs class, label, provider and a hex colour')
        patterns = {}
        for field in ('name', 'accession', 'description'):
            value = (r.get(field + '_pattern') or '').strip()
            if value:
                try:
                    patterns[field] = re.compile(value, re.IGNORECASE)
                except re.error as exc:
                    raise ValueError('Essential whitelist row ' + str(i) + ' has an invalid ' + field + ' regex: ' + str(exc))
        if not patterns:
            raise ValueError('Essential whitelist row ' + str(i) + ' has no pattern')
        rules.append(dict(cls=r['class'], label=r['label'], priority=int(r['priority']), color=r['color'],
                          provider=r['provider'], patterns=patterns, note=r.get('note', ''), row=i))
    return dict(path=str(path), sha256=sha(path), rules=rules)


def essential_rule(row, whitelist):
    """The highest-priority rule a source feature row satisfies, or None."""
    if not whitelist:
        return None
    best = None
    for rule in whitelist['rules']:
        if rule['provider'] != row.get('provider'):
            continue
        if all(p.search(str(row.get(field) or '')) for field, p in rule['patterns'].items()):
            if best is None or rule['priority'] < best['priority']:
                best = rule
    if best is not None:
        return best
    # A view may select a few Pfam identifiers, but that must not hide a named
    # UniProt domain (e.g. FES F-BAR or RHOT2 EF hands). Keep the literal source
    # name: historical snapshots can contain truncated descriptions. Different
    # names need different classes so coalescing never merges unrelated domains.
    provider = row.get('provider')
    description = str(row.get('description') or '')
    label = (description[len('DOMAIN-'):] if provider == 'uniprot_features'
             and description.startswith('DOMAIN-') else description
             if provider == 'reviewed_uniprot_snapshot' and row.get('name') == 'Domain' else '')
    if label:
        return dict(cls='NAMED_DOMAIN:' + label.casefold(), label=label,
                    priority=100, color='#4477AA', provider=provider,
                    note='Mandatory named UniProt domain; literal source label, no new prediction')
    motif = (description[len('MOTIF-'):] if provider == 'uniprot_features'
             and description.startswith('MOTIF-') else description
             if provider == 'reviewed_uniprot_snapshot' and row.get('name') == 'Motif' else '')
    if motif:
        return dict(cls='NAMED_MOTIF:' + motif.casefold(), label=motif,
                    priority=101, color='#AA4499', provider=provider,
                    note='Mandatory named UniProt motif; literal source label, no sequence inference')
    return None


def _coalesce_overlapping_same_annotation(display, audit_by_id):
    """One display bar per annotation when the same provider/type/description overlaps itself.

    The EnsMart100 Ensembl table repeats a family domain once per reference protein, with
    coordinates that differ by a few residues (CLDND1: PMP22/EMP/MP20/Claudin on 7 proteins).
    The longest placed interval represents the group; the others stay in the display audit
    as coalesced_overlapping_same_annotation with their reference proteins recorded.
    No coordinate is invented: the representative keeps its own source interval.
    """
    groups = defaultdict(list)
    for f in display:
        if f['A_start'] == '':
            # Unplaced copies of one annotation at the same reference coordinates collapse to one.
            groups[('__unplaced__', f['provider'], f['feature_type'], f['description'], str(f['reference_start']), str(f['reference_end']))].append(f)
        else:
            groups[(f['provider'], f['feature_type'], f['description'])].append(f)
    kept = []
    for key, items in groups.items():
        if len(items) == 1:
            kept.extend(items)
            continue
        if key[0] == '__unplaced__':
            rep = items[0]
            for f in items[1:]:
                rep['source_feature_ids'].extend(f['source_feature_ids'])
                rep.setdefault('coalesced_reference_proteins', []).append(f['reference_protein'])
                for fid in f['source_feature_ids']:
                    audit_by_id[fid]['display_role'] = 'coalesced_identical_display_annotation'
            kept.append(rep)
            continue
        items.sort(key=lambda f: (int(f['A_start']), -int(f['A_end'])))
        chains = []
        for f in items:
            if chains and int(f['A_start']) <= int(chains[-1][-1]['A_end']) + 1:
                chains[-1].append(f)
            else:
                chains.append([f])
        for chain in chains:
            rep = max(chain, key=lambda f: (int(f['A_end']) - int(f['A_start']), -int(f['A_start'])))
            for f in chain:
                if f is rep:
                    continue
                rep['source_feature_ids'].extend(f['source_feature_ids'])
                rep.setdefault('coalesced_reference_proteins', []).append(f['reference_protein'])
                for fid in f['source_feature_ids']:
                    audit_by_id[fid]['display_role'] = 'coalesced_overlapping_same_annotation'
            kept.append(rep)
    return _coalesce_same_essential_class(kept, audit_by_id)


def _coalesce_same_essential_class(display, audit_by_id):
    """One bar per essential class where its rows overlap, across providers.

    A kinase domain is reported by Pfam, InterPro and UniProt on nearly the same interval, and a
    zinc finger by UniProt and Pfam. The figure shows the class once, on the longest interval;
    the other rows stay in the audit as coalesced_same_essential_class with their providers.
    Non-overlapping members of a class (the four TM helices) stay separate.
    """
    groups = defaultdict(list)
    kept = []
    for f in display:
        if f.get('essential_class') and f['A_start'] != '':
            groups[f['essential_class']].append(f)
        else:
            kept.append(f)
    for cls, items in groups.items():
        items.sort(key=lambda f: (int(f['A_start']), -int(f['A_end'])))
        chains, chain_end = [], None
        for f in items:
            if chains and int(f['A_start']) <= chain_end:
                chains[-1].append(f); chain_end = max(chain_end, int(f['A_end']))
            else:
                chains.append([f]); chain_end = int(f['A_end'])
        for chain in chains:
            # The provider that resolves the chain into the most separate intervals keeps its
            # bars (two GATA zinc fingers from UniProt beat one InterPro row spanning both);
            # ties go to curated UniProt coordinates, then to the longest interval.
            by_provider = defaultdict(list)
            for f in chain:
                by_provider[f['provider']].append(f)
            rank = {'reviewed_uniprot_snapshot': 0, 'uniprot_features': 1, 'literature_construct_snapshot': 2, 'protein_features': 3}
            best = max(by_provider, key=lambda prov: (len(_disjoint(by_provider[prov])), -rank.get(prov, 9)))
            reps = _disjoint(by_provider[best])
            for f in chain:
                if f in reps:
                    continue
                target = max(reps, key=lambda r: min(int(r['A_end']), int(f['A_end'])) - max(int(r['A_start']), int(f['A_start'])))
                target['source_feature_ids'].extend(f['source_feature_ids'])
                target.setdefault('coalesced_reference_proteins', []).append(f['reference_protein'])
                target.setdefault('coalesced_providers', []).append(f['provider'] + ':' + f['description'])
                for fid in f['source_feature_ids']:
                    audit_by_id[fid]['display_role'] = 'coalesced_same_essential_class'
            kept.extend(reps)
    # Second pass on the B side: an essential row that has no A coordinate but sits on B inside a kept row
    # of the same class (a PROSITE C2H2 profile next to the curated UniProt zinc finger of the other isoform,
    # GFI1B 2026-10-24) is the same feature and is folded into it; otherwise it would be drawn as an extra
    # unnumbered bar on B and listed as "not in this protein" on A.
    placed_by_class = defaultdict(list)
    for f in kept:
        if f.get('essential_class') and f['A_start'] != '' and f.get('B_native_start') not in ('', None):
            placed_by_class[f['essential_class']].append(f)
    final = []
    for f in kept:
        if f.get('essential_class') and f['A_start'] == '' and f.get('B_native_start') not in ('', None):
            lo, hi = int(f['B_native_start']), int(f['B_native_end'])
            hosts = [g for g in placed_by_class[f['essential_class']]
                     if int(g['B_native_start']) <= hi and lo <= int(g['B_native_end'])]
            if hosts:
                target = max(hosts, key=lambda g: min(int(g['B_native_end']), hi) - max(int(g['B_native_start']), lo))
                target['source_feature_ids'].extend(f['source_feature_ids'])
                target.setdefault('coalesced_reference_proteins', []).append(f['reference_protein'])
                target.setdefault('coalesced_providers', []).append(f['provider'] + ':' + f['description'])
                for fid in f['source_feature_ids']:
                    audit_by_id[fid]['display_role'] = 'coalesced_same_essential_class'
                continue
        final.append(f)
    return final


def _disjoint(items):
    """Greedy set of non-overlapping members, longest first, in A order."""
    chosen = []
    for f in sorted(items, key=lambda f: (-(int(f['A_end']) - int(f['A_start'])), int(f['A_start']))):
        if all(int(f['A_end']) < int(g['A_start']) or int(f['A_start']) > int(g['A_end']) for g in chosen):
            chosen.append(f)
    return sorted(chosen, key=lambda f: int(f['A_start']))


def _number_essential(display):
    """TM1..TMn in N-to-C order on A; a single member keeps the bare class label."""
    by_class = defaultdict(list)
    for f in display:
        if f.get('essential_class'):
            by_class[f['essential_class']].append(f)
    for cls, items in by_class.items():
        placed = sorted((f for f in items if f['A_start'] != ''), key=lambda f: (int(f['A_start']), int(f['A_end'])))
        for i, f in enumerate(placed, 1):
            f['class_display'] = f['class_label'] + (str(i) if len(placed) > 1 else '')
        for f in items:
            if f['A_start'] == '':
                f['class_display'] = f['class_label']
COLORS = dict(reference='#4477AA', identical='#228833', missing='#BCC3CA',
              altered='#EE7733', unresolved='#AA3377', text='#192B3B')


def runs(values):
    """Closed, 1-based intervals, preserving singleton residues."""
    out=[]
    for _, group in itertools.groupby(enumerate(sorted(set(values))), lambda x:x[1]-x[0]):
        seq=[x[1] for x in group];out.append((seq[0],seq[-1]))
    return out


def alignment_segments(reference, query, mapping):
    """Each reference residue has one explicit state; query insertions stay separate."""
    by_state=defaultdict(list)
    for i, aa in enumerate(reference):
        state='deleted' if i not in mapping else 'identical' if aa==query[mapping[i]] else 'substituted'
        by_state[state].append(i+1)
    rows=[]
    for state, residues in by_state.items():
        for lo,hi in runs(residues):
            # Adjacent reference positions can flank an insertion. Split these
            # events so their query coordinates never imply an ungapped match.
            points=list(range(lo-1,hi))
            chunks=[];chunk=[]
            for i in points:
                if chunk and state!='deleted' and mapping[i]!=mapping[chunk[-1]]+1:
                    chunks.append(chunk);chunk=[]
                chunk.append(i)
            if chunk:chunks.append(chunk)
            for chunk in chunks:
                rows.append(dict(state=state,A_start=chunk[0]+1,A_end=chunk[-1]+1,
                                 B_start='' if state=='deleted' else mapping[chunk[0]]+1,
                                 B_end='' if state=='deleted' else mapping[chunk[-1]]+1))
    mapped_query=set(mapping.values())
    for lo,hi in runs(q+1 for q in range(len(query)) if q not in mapped_query):
        preceding=[i+1 for i,q in mapping.items() if q+1<lo]
        rows.append(dict(state='inserted',A_start='',A_end='',B_start=lo,B_end=hi,
                         A_anchor=max(preceding)+.5 if preceding else .5))
    return sorted(rows,key=lambda r:(r['A_start'] if r['A_start']!='' else len(reference)+1,r['B_start'] or 0,r['state']))


def load_audit(path):
    base=Path(config.require(str(path),'completed protein-impact audit'))
    manifest=json.loads((base/'impact_manifest.json').read_text())
    if manifest.get('status')!='complete_with_explicit_annotation_limits':
        raise ValueError('Protein-impact audit is not complete')
    needed=['native_model_audit.tsv','pair_summary.tsv','feature_projections.tsv','native_proteins.fasta']
    for name in needed:
        if name not in manifest['outputs'] or sha(base/name)!=manifest['outputs'][name]:
            raise ValueError('Changed or unpinned impact result: '+name)
    models=table(base/'native_model_audit.tsv')
    keys=[(r['gene_symbol'],r['isoform_id']) for r in models]
    expected=[(r['gene_symbol'],r['isoform_id']) for r in manifest['required_roster']]
    if keys!=expected or len(set(keys))!=len(keys):
        raise ValueError('Complete independent native roster was not preserved')
    sequences={}
    for rec in SeqIO.parse(base/'native_proteins.fasta','fasta'):
        if rec.id in sequences:raise ValueError('Duplicate native protein: '+rec.id)
        sequences[rec.id]=str(rec.seq)
    if sequences.keys()!={r['isoform_id'] for r in models}:
        raise ValueError('Native FASTA does not cover the complete impact roster')
    for r in models:
        seq=sequences[r['isoform_id']]
        if len(seq)!=int(r['actual_aa']) or hashlib.sha256(seq.encode()).hexdigest()!=r['protein_sha256']:
            raise ValueError('Native sequence/length differs from audited protein: '+r['isoform_id'])
    return base,manifest,models,sequences,table(base/'pair_summary.tsv'),table(base/'feature_projections.tsv')


def selected_pairs(pairs, models, genes=None, pair_file=None):
    available={r['gene_symbol'] for r in models}
    required=None
    if pair_file:
        supplied=table(pair_file)
        required=[(r.get('gene_symbol') or r.get('gene'),r.get('isoform_A'),r.get('isoform_B')) for r in supplied]
        if not required or len(set(required))!=len(required) or any(not all(k) for k in required):
            raise ValueError('Pair TSV requires unique gene_symbol,isoform_A,isoform_B')
    wanted=(genes.split(',') if genes else list(dict.fromkeys(k[0] for k in required)) if required
            else sorted({p['gene_symbol'] for p in pairs}))
    wanted=[g.strip() for g in wanted]
    if not wanted or len(wanted)!=len(set(wanted)) or any(not g for g in wanted):
        raise ValueError('Empty or duplicate requested genes')
    if set(wanted)-available:raise ValueError('Unresolved requested genes: '+repr(sorted(set(wanted)-available)))
    order={g:i for i,g in enumerate(wanted)}
    kept=[p for p in pairs if p['gene_symbol'] in order]
    if required:
        index={(p['gene_symbol'],p['isoform_A'],p['isoform_B']):p for p in kept}
        if set(required)-index.keys():raise ValueError('Requested pair is not an audited native comparison')
        kept=[index[k] for k in required]
    missing=set(wanted)-{p['gene_symbol'] for p in kept}
    if missing:raise ValueError('Requested genes without an audited pair: '+repr(sorted(missing)))
    return sorted(kept,key=lambda p:(order[p['gene_symbol']],p['isoform_A'],p['isoform_B']))


def feature_reference_selection(values,pairs,feature_rows):
    """Explicit figure-view selection; never alter the underlying impact audit."""
    genes={p['gene_symbol'] for p in pairs};selected={}
    for value in values or []:
        gene,separator,protein=value.partition('=')
        if not separator or not gene or not protein or gene not in genes:
            raise ValueError('--feature-reference requires a selected GENE=REFERENCE_PROTEIN')
        if not any(r['gene_symbol']==gene and r['reference_protein']==protein for r in feature_rows):
            raise ValueError('Requested reference protein absent from audited gene: '+value)
        refs=selected.setdefault(gene,[])
        if protein in refs:raise ValueError('Duplicate feature reference: '+value)
        refs.append(protein)
    return selected


def feature_plan(pair, feature_rows, sequences, feature_types=DEFAULT_TYPES, sources=PRIMARY_SOURCES, reference_proteins=None,
                 essential=None):
    """One explicit display policy; source rows are never silently discarded.

    Functional reviewed UniProt and literature annotations are the default view.
    Other sources, types and unresolved placements remain in the display audit.
    Essential classes (``essential``: the loaded whitelist; None loads the packaged file)
    are displayed whatever the source and type selection says.
    Identical displayed intervals/annotations are coalesced only for presentation, and so
    is the same annotation repeated on overlapping intervals from several reference proteins.
    """
    if essential is None:essential=load_essential_whitelist()
    a,b=pair['isoform_A'],pair['isoform_B'];seq_a,seq_b=sequences[a],sequences[b]
    mapping=align_mapping(seq_a,seq_b)
    relevant=[r for r in feature_rows if r['gene_symbol']==pair['gene_symbol'] and r['isoform_id'] in (a,b)]
    by_key={(r['isoform_id'],r['provider'],r['feature_id']):r for r in relevant}
    if len(by_key)!=len(relevant):raise ValueError('Duplicate native feature instance in impact audit')
    audit=[];display=[];seen={};audit_by_id={}
    for r in relevant:
        rule=essential_rule(r,essential)
        reason=('paired_B_record' if r['isoform_id']==b else 'displayed' if rule else 'source_not_in_view' if r['provider'] not in sources
                else 'type_not_in_view' if r['name'] not in feature_types else 'reference_not_in_view'
                if reference_proteins is not None and r['reference_protein'] not in reference_proteins else 'displayed')
        entry=dict(gene_symbol=pair['gene_symbol'],isoform_A=a,isoform_B=b,source_isoform=r['isoform_id'],provider=r['provider'],
                   feature_id=r['feature_id'],feature_type=r['name'],description=r['description'],source_status=r['status'],display_role=reason,
                   essential_class=rule['cls'] if rule else '')
        audit.append(entry)
        if r['isoform_id']!=b:audit_by_id[r['feature_id']]=entry
        if reason!='displayed':continue
        essential_fields=dict(essential_class=rule['cls'],class_label=rule['label'],class_priority=rule['priority'],
                              class_color=rule['color'],class_display=rule['label']) if rule else dict(
                              essential_class='',class_label='',class_priority='',class_color='',class_display='')
        rb=by_key.get((b,r['provider'],r['feature_id']))
        if rb is None:raise ValueError('Missing required paired feature projection: '+r['feature_id'])
        if r['status']!='intact_sequence':
            # Preserve functional features that cannot be placed on A, including
            # intact B-only sequences. Never fabricate an A-axis position.
            entry['display_role']='displayed_without_A_coordinate'
            item=dict(gene_symbol=pair['gene_symbol'],isoform_A=a,isoform_B=b,provider=r['provider'],feature_id=r['feature_id'],
                      feature_type=r['name'],description=r['description'],reference_protein=r['reference_protein'],
                      reference_start=r['start'],reference_end=r['end'],A_start='',A_end='',
                      B_native_start=rb['query_start'],B_native_end=rb['query_end'],
                      feature_length=int(r.get('reference_length') or 0),identical_residues='',altered_residues='',absent_residues='',
                      inserted_residues='',insertion_anchors=[],identical_segments=[],altered_segments=[],absent_segments=[],B_native_segments=[],
                      status='intact_only_in_B_unplaced' if rb['status']=='intact_sequence' and r.get('sequence_present') in ('False',False) else 'unplaced_annotation',
                      source_A_status=r['status'],source_B_status=rb['status'],pair_projection_status='not_computed_without_A_coordinate',
                      projection_disagreement=False,evidence=r.get('evidence',''),source_feature_ids=[r['feature_id']],**essential_fields)
            display.append(item)
            continue
        start,end=int(r['query_start']),int(r['query_end'])
        if not 1<=start<=end<=len(seq_a):raise ValueError('Invalid native A feature placement')
        projection=project_feature(seq_a,seq_b,start,end,mapping)
        mapped=[i for i in range(start-1,end) if i in mapping]
        identical=[i+1 for i in mapped if seq_a[i]==seq_b[mapping[i]]]
        altered=[i+1 for i in mapped if seq_a[i]!=seq_b[mapping[i]]]
        absent=[i+1 for i in range(start-1,end) if i not in mapping]
        mapped_q=[mapping[i]+1 for i in mapped]
        insertions=[e for e in alignment_segments(seq_a,seq_b,mapping) if e['state']=='inserted' and mapped_q
                    and min(mapped_q)<e['B_start'] and e['B_end']<max(mapped_q)]
        n_inserted=sum(e['B_end']-e['B_start']+1 for e in insertions)
        exact_span=(projection['status']=='intact_sequence' and len(mapped)==end-start+1 and
                    len(identical)==end-start+1 and runs(mapped_q)==[(projection['query_start'],projection['query_end'])])
        status=('retained' if exact_span else 'ambiguous_placement' if projection['status'] in ('ambiguous_repeated_sequence','intact_sequence')
                else 'absent' if not mapped else 'partial' if absent else 'altered')
        # Independent ref-to-B audit can differ from pair-to-B context when the
        # reference is not A. Never hide that disagreement or assert a loss.
        disagreement=projection['status']!=rb['status']
        if disagreement:status='projection_disagreement'
        item=dict(gene_symbol=pair['gene_symbol'],isoform_A=a,isoform_B=b,provider=r['provider'],feature_id=r['feature_id'],
                  feature_type=r['name'],description=r['description'],reference_protein=r['reference_protein'],
                  reference_start=r['start'],reference_end=r['end'],A_start=start,A_end=end,
                  B_native_start=min(mapped_q) if mapped_q else '',B_native_end=max(mapped_q) if mapped_q else '',
                  feature_length=end-start+1,identical_residues=len(identical),altered_residues=len(altered),absent_residues=len(absent),
                  inserted_residues=n_inserted,insertion_anchors=[(e['A_anchor'],e['B_end']-e['B_start']+1) for e in insertions],
                  identical_segments=runs(identical),altered_segments=runs(altered),absent_segments=runs(absent),
                  B_native_segments=runs(mapped_q),status=status,source_B_status=rb['status'],pair_projection_status=projection['status'],
                  projection_disagreement=disagreement,evidence=r.get('evidence',''),source_feature_ids=[r['feature_id']],**essential_fields)
        key=(start,end,r['name'],r['description'],r['provider'])
        if key in seen:
            seen[key]['source_feature_ids'].append(r['feature_id']);entry['display_role']='coalesced_identical_display_annotation'
        else:seen[key]=item;display.append(item)
    display=_coalesce_overlapping_same_annotation(display,audit_by_id)
    _number_essential(display)
    display.sort(key=lambda f:(f['A_start'] if f['A_start']!='' else math.inf,f['A_end'] if f['A_end']!='' else math.inf,f['feature_type'],f['description'],f['provider']))
    for i,f in enumerate(display,1):f['display_id']='F'+str(i)
    return display,audit,alignment_segments(seq_a,seq_b,mapping)


def donor_context(base,manifest):
    selection=manifest.get('read_selection')
    if not selection:return 'RNA read support was not requested.',{}
    libraries=set(selection['source_libraries'])
    candidates=[p for p in manifest['inputs'] if Path(p).name=='metadata.tsv']
    details=dict(states=selection['states'],conditions=selection['conditions'],source_libraries=len(libraries),source_state_observations=len(selection['observations']))
    if len(candidates)==1:
        p=Path(config.require(candidates[0],'audited metadata'))
        if sha(p)!=manifest['inputs'][candidates[0]]:raise ValueError('Audited metadata changed')
        selected=[r for r in table(p) if r['library'] in libraries]
        if {r['library'] for r in selected}!=libraries:raise ValueError('Incomplete RNA-support source roster')
        details['donors']=len({r['original_library'] for r in selected})
    unit=str(details['donors'])+' donor libraries' if 'donors' in details else str(len(libraries))+' source libraries'
    text=('RNA support: '+', '.join(selection['states'])+'; '+', '.join(selection['conditions'] or ['all conditions'])+
          '; pooled assigned reads from '+unit+' ('+str(len(selection['observations']))+' source/state rows).')
    return text,details


def _mpl():
    import matplotlib
    import logging
    matplotlib.use('Agg')
    matplotlib.rcdefaults()
    logging.getLogger('fontTools.subset').setLevel(logging.WARNING)
    matplotlib.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'ps.fonttype':42,
        'svg.fonttype':'none','svg.hashsalt':SCHEMA,'axes.spines.top':False,'axes.spines.right':False})
    import matplotlib.pyplot as plt
    return plt


def _clean_label(value):
    return re.sub(r'\s*\(literature construct\)\s*','',value).strip() or 'Unnamed feature'


def _text(ax,x,y,value,**kwargs):
    role=kwargs.pop('role','label')
    obj=ax.text(x,y,value,**kwargs);obj.set_gid(role);return obj


def _bar(ax,lo,hi,y,color,height=.12,**kwargs):
    from matplotlib.patches import Rectangle
    return ax.add_patch(Rectangle((lo-.5,y-height/2),hi-lo+1,height,facecolor=color,**kwargs))


def _state_text(f):
    n=f['feature_length'];k=f['identical_residues'];status=f['status']
    if status=='retained':label='Retained: '+str(k)+'/'+str(n)+' aa'
    elif status=='absent':label='Absent: 0/'+str(n)+' aa'
    elif status=='partial':label='Partial: '+str(k)+'/'+str(n)+' identical aa'
    elif status=='altered':label=('Insertion: '+str(f['inserted_residues'])+' aa' if f['inserted_residues'] else 'Altered: '+str(f['altered_residues'])+' AA change(s)')
    elif status=='intact_only_in_B_unplaced':label='Intact sequence only in B'
    elif status=='unplaced_annotation':label='Annotation not placed on A'
    else:label='Uncertain placement'
    coords=('B: '+str(f['B_native_start'])+'-'+str(f['B_native_end'])+' aa') if f['B_native_start']!='' else 'No aligned B residues'
    return label+'\n'+coords


def render_pair(pair, features, events, sequences, models, context, page_index=1, page_count=1):
    """One template in inches, with fixed row height and no gene-specific styling."""
    plt=_mpl()
    from matplotlib.patches import Patch
    from matplotlib.ticker import MaxNLocator
    count=len(features);height=4.35+.63*max(count,3)
    fig=plt.figure(figsize=(12,height),facecolor='white')
    gene=pair['gene_symbol'];a,b=pair['isoform_A'],pair['isoform_B'];length=len(sequences[a]);margin=max(length*.025,2)
    top=height-.32
    fig.text(.045,top/height,gene+' | native protein features',fontsize=17,fontweight='bold',va='top',color=COLORS['text'])
    fig.text(.045,(top-.36)/height,'Sequence-based prediction'+(' | page '+str(page_index)+'/'+str(page_count) if page_count>1 else ''),fontsize=10.5,color='#526373',va='top')
    # Full literal native IDs remain legible; no accession or length substitution.
    for side,iso,offset in [('A',a,.79),('B',b,1.11)]:
        m=models[iso];reads=m['assigned_reads'];reads=str(int(float(reads))) if reads!='' else 'not requested'
        text=side+': '+iso+'  |  '+str(len(sequences[iso]))+' aa  |  '+reads+' assigned reads  |  '+m['nmd']
        fig.text(.045,(top-offset)/height,text,fontsize=9.5,va='top',color=COLORS['text'])
    plot_left=.305;plot_width=.49
    ax=fig.add_axes([plot_left,(height-2.5)/height,plot_width,.66/height])
    ax.set_xlim(.5-margin,length+.5+margin);ax.set_ylim(-.5,1.45)
    _bar(ax,1,length,1,COLORS['reference'],height=.18)
    for e in events:
        if e['state']=='inserted':
            ax.vlines(e['A_anchor'],-.22,.22,color=COLORS['altered'],linewidth=1.3)
            continue
        lo,hi=e['A_start'],e['A_end']
        if e['state']=='deleted':_bar(ax,lo,hi,0,'none',height=.18,edgecolor=COLORS['missing'],linestyle=(0,(3,2)),linewidth=.8)
        else:_bar(ax,lo,hi,0,COLORS['identical' if e['state']=='identical' else 'altered'],height=.18)
    ax.set_yticks([1,0],['A','B']);ax.tick_params(axis='y',length=0,pad=7)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6,integer=True));ax.tick_params(axis='x',labelsize=9,length=3)
    ax.set_xlabel('Aligned amino-acid position in native A',fontsize=10,labelpad=4)
    for spine in ['left','right','top']:ax.spines[spine].set_visible(False)
    fig.text(.045,(height-2.03)/height,'Whole-protein alignment',fontsize=10,fontweight='bold',va='top')
    added=sum(e['B_end']-e['B_start']+1 for e in events if e['state']=='inserted')
    fig.text(.825,(height-1.9)/height,'Native B: '+str(len(sequences[b]))+' aa\n'+str(added)+' inserted / unaligned aa',fontsize=9.3,va='top')
    feature_top=height-3.05
    fig.text(.045,(feature_top+.23)/height,'Annotated feature',fontweight='bold',fontsize=10)
    fig.text(.825,(feature_top+.23)/height,'Native B sequence',fontweight='bold',fontsize=10)
    bounds=[]
    for row,f in enumerate(features):
        center=feature_top-row*.63-.22
        track=fig.add_axes([plot_left,(center-.21)/height,plot_width,.42/height])
        track.set_xlim(ax.get_xlim());track.set_ylim(-.4,1.4);track.axis('off')
        lo,hi=f['A_start'],f['A_end']
        if lo=='':
            _text(track,.5,.5,'No A-coordinate placement',transform=track.transAxes,ha='center',va='center',fontsize=9,color=COLORS['unresolved'])
        else:_bar(track,lo,hi,1,COLORS['reference'],height=.37,edgecolor='white',linewidth=.3)
        if lo=='':pass
        elif f['status'] in ('ambiguous_placement','projection_disagreement'):
            _bar(track,lo,hi,0,'none',height=.37,edgecolor=COLORS['unresolved'],hatch='///',linewidth=.7)
        else:
            for x,y in f['absent_segments']:_bar(track,x,y,0,'none',height=.37,edgecolor=COLORS['missing'],linestyle=(0,(3,2)),linewidth=.8)
            for x,y in f['identical_segments']:_bar(track,x,y,0,COLORS['identical'],height=.37,edgecolor='white',linewidth=.3)
            for x,y in f['altered_segments']:_bar(track,x,y,0,COLORS['altered'],height=.37,edgecolor='white',linewidth=.3)
            for anchor,n in f['insertion_anchors']:track.vlines(anchor,-.3,.3,color=COLORS['altered'],linewidth=1.2)
        # Fixed paired A/B rows preserve color-independent interpretability.
        _text(track,-.022,1,'A',transform=track.get_yaxis_transform(),fontsize=8.5,ha='right',va='center',role='pair_row_A')
        _text(track,-.022,0,'B',transform=track.get_yaxis_transform(),fontsize=8.5,ha='right',va='center',role='pair_row_B')
        desc=_clean_label(f['description'])
        prefix=f['display_id']+'  '+f['feature_type']
        label=prefix+'\n'+'\n'.join(textwrap.wrap(desc,width=34))
        if f['provider']=='literature_construct_snapshot':label+=' [literature]'
        left=fig.text(.045,center/height,label,fontsize=9.1,va='center',linespacing=1.1,color=COLORS['text']);left.set_gid('feature_label_'+f['display_id'])
        coords='A: '+str(lo)+'-'+str(hi)+' aa' if lo!='' else 'A: no intact placement'
        right=fig.text(.825,center/height,_state_text(f)+'\n'+coords,fontsize=9.1,va='center',linespacing=1.1,color=COLORS['text']);right.set_gid('feature_status_'+f['display_id'])
        bounds.extend([left,right])
    legend_y=1.13/height
    handles=[Patch(facecolor=COLORS['reference'],label='A annotation'),Patch(facecolor=COLORS['identical'],label='Identical B residues'),
             Patch(facecolor='none',edgecolor=COLORS['missing'],linestyle='--',label='Sequence absent'),
             Patch(facecolor=COLORS['altered'],label='AA change'),Patch(facecolor='none',edgecolor=COLORS['unresolved'],hatch='///',label='Uncertain placement')]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.04,legend_y),ncol=5,frameon=False,fontsize=9,handlelength=1.5,columnspacing=1.4)
    sources=sorted({f['reference_protein'] for f in features})
    caption=('A is the longer native model; B is the alternative. Both rows use A coordinates; B native coordinates are listed at right. '
             'Feature widths are to scale; orange ticks indicate B insertions (tick width is not insertion length). Annotations: '+', '.join(sources)+'. '+context+' '
             'Domains, DNA binding, zinc fingers, motifs, coiled coils and regions are shown; other source records remain in the display audit. '
             'No statistical test or confidence interval: deterministic sequence projection. RNA support does not establish translation or altered activity in vivo.')
    if any(f['provider']=='literature_construct_snapshot' for f in features):
        citation=[]
        for f in features:
            if f['provider']=='literature_construct_snapshot':
                for ev in json.loads(f['evidence'] or '[]'):
                    url=ev.get('url','');citation.append(url)
        caption+=' Literature annotation: '+', '.join(sorted(set(citation)))+'.'
    foot=fig.text(.045,.105/height,textwrap.fill(caption,width=180),fontsize=8.4,va='bottom',linespacing=1.12,color='#455564');foot.set_gid('caption')
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    layout=[];width_px,height_px=fig.canvas.get_width_height()
    for obj in [*bounds,foot]:
        box=obj.get_window_extent(renderer)
        if box.x0<0 or box.y0<0 or box.x1>width_px+1 or box.y1>height_px+1:
            raise ValueError('Figure text falls outside page: '+obj.get_gid())
        layout.append(dict(role=obj.get_gid(),bounds_px=list(box.bounds)))
    # Explicitly reject adjacent labels/statuses colliding in the standard lanes.
    for group in ('feature_label_','feature_status_'):
        texts=[v for v in bounds if v.get_gid().startswith(group)]
        for first,second in zip(texts,texts[1:]):
            if first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer)):
                raise ValueError('Adjacent feature labels overlap; full source text must fit standard rows')
    return fig,layout


def architecture_features(features, side):
    """Use existing placements, never mark substitutions as lost domains.

    Both proteins use their own native residue coordinates. This is an
    architecture view, not an amino-acid difference/alignment display.
    """
    placed=[];unplaced=[]
    for f in features:
        if side=='A':
            segments=[] if f['A_start']=='' else [(f['A_start'],f['A_end'])]
            uncertain=f['A_start']==''
        else:
            uncertain=f['status'] in ('ambiguous_placement','projection_disagreement','unplaced_annotation')
            mapped=f['B_native_segments']
            # A domain architecture shows the best mapped native interval,
            # including small insertions; it is not a residue difference map.
            segments=[(mapped[0][0],mapped[-1][1])] if mapped else []
            if f['status']=='intact_only_in_B_unplaced':
                segments=[(f['B_native_start'],f['B_native_end'])];uncertain=False
        if uncertain:
            unplaced.append(f['display_id']);continue
        if segments:
            segments=[(int(lo),int(hi)) for lo,hi in segments]
            placed.append(dict(feature=f,segments=segments,partial=(side=='B' and f['status']=='partial')))
    return placed,unplaced


def reference_architecture(features):
    """Essential classes that NO drawn isoform can place, returned on their own reference protein.

    Nathan, 2026-10-25 (IKZF1 Pro-B-1): a panel whose DNA-binding or TM annotations are unplaceable on
    both isoforms drew no domain at all and stated the loss only in a text line. The essential-class rule
    is that such a class is always SHOWN, so it is drawn here on a reference row carrying the reference
    protein's own coordinates. The row is never a measured isoform: its extent is the furthest annotated
    reference residue in the panel, which is the only reference length pinned in the impact results.

    A class that is placed on A or B is left alone, so panels that already draw their essential bars are
    unchanged. One bar per distinct reference interval per class, with the providers kept in the audit.
    """
    placed_classes=set()
    for side in ('A','B'):
        for item in architecture_features(features,side)[0]:
            if item['feature'].get('essential_class'):placed_classes.add(item['feature']['essential_class'])
    by_reference=defaultdict(dict)
    for f in features:
        cls=f.get('essential_class')
        if not cls or cls in placed_classes or f['A_start']!='':continue
        if f.get('reference_start') in ('',None) or f.get('reference_end') in ('',None):continue
        start,end=int(f['reference_start']),int(f['reference_end'])
        if end<start:raise ValueError('Reference annotation ends before it starts: '+str(f.get('feature_id')))
        by_reference[f['reference_protein']].setdefault((cls,start,end),f)
    # One row per essential class, on the MOST COMPLETE reference protein that carries it (furthest
    # annotated residue, then bar count): the row exists to show the architecture the isoforms lack, so a
    # partial reference must not replace the full-length one. The same zinc fingers repeated over ten
    # reference proteins is noise; every row stays in the display audit.
    chosen={}
    for cls in {c for bars in by_reference.values() for c,_,_ in bars}:
        counts=[(max(e for _,_,e in bars),sum(1 for c,_,_ in bars if c==cls),reference)
                for reference,bars in by_reference.items() if any(c==cls for c,_,_ in bars)]
        chosen[cls]=max(counts)[2]
    rows=[]
    for reference,bars in sorted(by_reference.items()):
        kept={k:f for k,f in bars.items() if chosen[k[0]]==reference}
        if not kept:continue
        items=[dict(feature=f,segments=[(start,end)],partial=False) for (cls,start,end),f in sorted(kept.items())]
        rows.append(dict(reference_protein=reference,extent=max(e for _,_,e in kept),features=items))
    return rows


def unplaced_summary(features, unplaced_ids):
    """Readable text for annotations that cannot be placed on a native protein: names with counts,
    never bare running ids (``Kinase x3, TM x2, Ig_I-set``)."""
    by_id={f['display_id']:f for f in features}
    counts=defaultdict(int);order=[]
    for fid in unplaced_ids:
        f=by_id.get(fid)
        name=(f['class_label'] if f and f.get('essential_class') else architecture_label(f) if f else fid)
        if name not in counts:order.append(name)
        counts[name]+=1
    return ', '.join(n+(' x'+str(counts[n]) if counts[n]>1 else '') for n in order)


def architecture_label(feature):
    """Conventional feature abbreviations; never a gene-specific label rule.

    An essential-class feature is labelled by its class (TM1, Kinase, DNA-binding).
    """
    if feature.get('essential_class'):return feature.get('class_display') or feature['class_label']
    desc=_clean_label(feature['description'])
    # EnsMart100 UniProt rows read "DOMAIN-Ig-like C2-type": the key is the feature type, the rest the name.
    if feature.get('provider')=='uniprot_features' and re.match(r'^[A-Z_]{3,}-',desc):
        desc=desc.split('-',1)[1].strip() or desc
    signal=re.search(r'\((NLS\d*|NES\d*)\)',desc)
    if signal:return signal[1]
    common={'Disordered':'IDR','HMG box':'HMG','NR LBD':'LBD','NR C4-type':'C4',
            'Antp-type hexapeptide':'Hexapeptide','Nuclear localization signal':'NLS',
            'Homeobox; degenerate':'HD-like','Fork-head':'Forkhead',
            'Required for mitochondrial import':'Mito. import',
            'Mediates interaction with CHUK/IKKA and IKBKB/IKKB':'IKK-binding',
            'Required for transcriptional activation':'Activation','Interaction with WWP1':'WWP1-binding',
            'CRF-bd':'CRF-binding','Sperma_CUB_dom_sf':'CUB fold',
            'C2H2-type 1':'C2H2-1','C2H2-type 2':'C2H2-2','C2H2-type 3':'C2H2-3'}
    if desc in common:return common[desc]
    if feature['feature_type']=='DNA binding' and desc=='Nuclear receptor':return 'DBD'
    if desc.startswith('CBP/p300-binding'):return 'CBP/p300'
    # Always the source name, never a bare running id: a bar must be readable without the key.
    # Past 40 characters the name is cut with an explicit ellipsis; the full text stays in the key and tooltip.
    return desc if len(desc)<=40 else desc[:39].rstrip()+'\u2026'


def render_compact_pair(pair, features, events, sequences, models, context, page_index=1, page_count=1, figure=None, figure_width=12):
    """Conventional native-length protein backbones with keyed feature bars.

    Colors identify annotations, never amino-acid substitutions. Missing features
    leave no mark on the shorter protein; unplaced annotations remain explicit.
    No narrative caption is drawn. Full scientific provenance stays in sidecars.
    """
    if not math.isfinite(figure_width) or figure_width<=0:raise ValueError('Protein figure width must be positive and finite')
    plt=_mpl()
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import MaxNLocator
    palette=['#4477AA','#EE6677','#228833','#AA3377','#66CCEE','#CCBB44',
             '#009988','#BB5566','#332288','#AA4499','#44AA99','#999933']
    feature_key=lambda f:(f['feature_type'],f['description'],f['provider'])
    unique=list({feature_key(f):f for f in features}.values())
    labels_by_key={feature_key(f):architecture_label(f) for f in unique}
    # Use the existing short feature identifier when a key cannot fit its
    # column. Keep the complete source description in the adjacent column.
    from matplotlib.backends.backend_agg import RendererAgg
    from matplotlib.font_manager import FontProperties
    dpi=plt.rcParams['figure.dpi']
    key_renderer=RendererAgg(1,1,dpi)
    # A key name that is wider than its column is wrapped onto more lines, never replaced by the running
    # id (F1, F2 ...), which Nathan rejected on 2026-10-23/24 as meaningless.
    key_font=FontProperties(family='DejaVu Sans',weight='bold',size=8.5)
    key_width_px=.085*figure_width*dpi-4
    for f in unique:
        key=feature_key(f)
        text=labels_by_key[key]
        if key_renderer.get_text_width_height_descent(text,key_font,False)[0] > key_width_px:
            for chars in range(max(len(text)-1,4),3,-1):
                lines=textwrap.wrap(text,width=chars,break_long_words=True,break_on_hyphens=True)
                if all(key_renderer.get_text_width_height_descent(l,key_font,False)[0] <= key_width_px for l in lines):
                    labels_by_key[key]='\n'.join(lines);break
            else:
                labels_by_key[key]='\n'.join(textwrap.wrap(text,width=4,break_long_words=True))
    color_by_key={feature_key(f):palette[i % len(palette)] for i,f in enumerate(unique)}
    for f in unique:
        if f.get('essential_class'):color_by_key[feature_key(f)]=f['class_color']
    colors={f['display_id']:color_by_key[feature_key(f)] for f in features}
    is_essential=lambda f:bool(f.get('essential_class'))
    bar_label=lambda f:(f.get('class_display') or f['class_label']) if is_essential(f) else architecture_label(f)
    a,b=pair['isoform_A'],pair['isoform_B']
    # An essential class that NEITHER isoform can place is drawn on its reference protein instead of
    # being reduced to a text line (Nathan, 2026-10-25: the IKZF1 Pro-B-1 panel drew no DNA-binding bar).
    references=reference_architecture(features)
    length=max([len(sequences[a]),len(sequences[b])]+[r['extent'] for r in references])
    # Pack annotation bars using both their true interval and the width of the
    # conventional label, including partial-domain qualifiers.
    rows=[]
    for side,iso,reference in [('A',a,None),('B',b,None)]+[('REF'+str(i+1),None,r) for i,r in enumerate(references)]:
        placed,unplaced=(architecture_features(features,side) if reference is None
                         else ([dict(item) for item in reference['features']],[]))
        lanes=[]
        # Essential classes pack first, so they take the lanes nearest the backbone.
        for item in sorted(placed,key=lambda x:(0 if is_essential(x['feature']) else 1,x['segments'][0][0],x['segments'][-1][1],x['feature']['display_id'])):
            start,end=item['segments'][0][0],item['segments'][-1][1]
            label=bar_label(item['feature'])+(' (partial)' if item['partial'] else '')
            center=(start+end)/2
            half_label=len(label)*.055/(.56*figure_width)*length/2
            left=min(start,center-half_label);right=max(end,center+half_label)
            lane=next((i for i,last in enumerate(lanes) if left>last+length*.025),len(lanes))
            if lane==len(lanes):lanes.append(right)
            else:lanes[lane]=right
            item['lane']=lane
        rows.append(dict(side=side,iso=iso,reference=reference,features=placed,unplaced=unplaced,lanes=max(len(lanes),1)))
    # A row that lists unplaced annotations needs the extra text lines below its backbone.
    unplaced_lines={r['side']:(textwrap.wrap('Unplaced annotations: '+unplaced_summary(features,r['unplaced']),width=max(40,int(figure_width*7.5))) if r['unplaced'] else []) for r in rows}
    # Reserve a key only for coding fragments whose labels may not fit. Their
    # scientific widths stay exact; the full literal IDs and AA offsets remain.
    key_width=max(40,int(figure_width*7.5))
    small_key_reserve={}
    for r in rows:
        entries=[e['exon_label']+' ['+format(e['aa_start0'],'g')+','+format(e['aa_end0'],'g')+')'
                 for e in (models[r['iso']].get('coding_exons',[]) if r['iso'] else [])
                 if e.get('exon_label') and (e['aa_end0']-e['aa_start0'])/max(length,1)*.56*figure_width*72<36]
        small_key_reserve[r['side']]=(textwrap.wrap('Small coding segments (AA offsets): '+'; '.join(entries),
                                                  width=key_width) if entries else [])
    heights=[.80+.20*r['lanes']+.15*(len(unplaced_lines[r['side']])+len(small_key_reserve[r['side']])) for r in rows]
    key_lines=[]
    for f in unique:
        label=_clean_label(f['description'])
        if f['provider']=='literature_construct_snapshot':label+=' [literature]'
        if is_essential(f):
            if not f['essential_class'].startswith('NAMED_DOMAIN:'):
                label=f['class_label']+' (essential class): '+label
            label+=' ['+f['provider']+']'
        lines=textwrap.wrap(label,width=27) or [label]
        key_lines.append((f,lines))
    partial_key=any(item['partial'] for row in rows for item in row['features'])
    key_height=sum(.16*max(len(lines),labels_by_key[feature_key(f)].count('\n')+1)+.08
                   for f,lines in key_lines)+.35+(.16 if partial_key else 0)
    height=max(2.85,.65+sum(heights),.85+key_height)
    fig=figure if figure is not None else plt.figure(figsize=(figure_width,height),facecolor='white')
    fig.text(.025,1-.13/height,pair['gene_symbol'],fontsize=14,fontweight='bold',va='top')
    fig.text(.70,1-.16/height,'Protein features (predicted)',fontsize=10,va='top')
    if any(models[i].get('coding_exons') for i in (a,b)):
        fig.text(.70,1-.36/height,'Backbone colors: coding exons',fontsize=8,va='top')
        if any(e.get('exon_label') for i in (a,b) for e in models[i].get('coding_exons',[])):
            fig.text(.70,1-.51/height,'AltAnalyze IDs: E exon; I intron',fontsize=8,va='top')
    if partial_key:
        fig.text(.70,1-.67/height,'Partial: sequence alignment only',fontsize=8,va='top')
    ax=fig.add_axes([.14,.25/height,.56,(height-.65)/height])
    ax.set_xlim(-length*.015,length*1.025);ax.set_ylim(0,height-.65);ax.axis('off')
    top=height-.78
    for row,row_height in zip(rows,heights):
        iso=row['iso'];reference=row['reference']
        n=reference['extent'] if reference else len(sequences[iso])
        base=top-.20*row['lanes']-.13
        y_fig=(base+.25)/height
        fig.text(.025,y_fig,'reference' if reference else pair['gene_symbol']+'-'+('1' if row['side']=='A' else '2'),
                 fontsize=11,va='center')
        m={} if reference else models[iso]
        length_text=str(n)+(' aa annotated' if reference else ' aa')
        if m.get('reported_aa') not in (None,'') and int(m['reported_aa'])!=n:
            length_text+=' (table: '+str(m['reported_aa'])+')'
        fig.text(.025,y_fig-.23/height,length_text,fontsize=9,va='center')
        body_height=.23 if m.get('coding_exons') else .14
        small_exon_labels=[]
        ax.add_patch(Rectangle((.5,base-body_height/2),n,body_height,facecolor='#D6D6D6',edgecolor='#777777',linewidth=.65))
        if m.get('coding_exons'):
            for exon_index,exon in enumerate(m['coding_exons'],1):
                patch=Rectangle((.5+exon['aa_start0'],base-body_height/2),exon['aa_end0']-exon['aa_start0'],body_height,
                    facecolor=exon['color'],edgecolor='#555555',linewidth=.25,gid='coding_exon_'+row['side'])
                ax.add_patch(patch)
                if 'exon_label_status' in exon and exon['exon_label_status']!='assigned':
                    raise ValueError('Unresolved AltAnalyze exon label for plotted native model: '+iso)
                if exon.get('exon_label'):
                    rgb=tuple(int(exon['color'][i:i+2],16)/255 for i in (1,3,5))
                    linear=[c/12.92 if c<=.04045 else ((c+.055)/1.055)**2.4 for c in rgb]
                    luminance=sum(c*w for c,w in zip(linear,(.2126,.7152,.0722)))
                    color='black' if (luminance+.05)/.05>=1.05/(luminance+.05) else 'white'
                    label=ax.text(.5+(exon['aa_start0']+exon['aa_end0'])/2,base,exon['exon_label'],
                                  ha='center',va='center',fontsize=8,color=color,
                                  gid='exon_label_'+row['side']+'_'+str(exon_index))
                    # Small coding fragments retain an inside label, rotated only
                    # when necessary. Never silently omit an ID or distort width.
                    fitted=False
                    for rotation,size in [(0,8),(0,7),(90,7),(90,6.5)]:
                        label.set_rotation(rotation);label.set_fontsize(size)
                        fig.canvas.draw();renderer=fig.canvas.get_renderer()
                        text_box=label.get_window_extent(renderer);block_box=patch.get_window_extent(renderer)
                        if text_box.width<=block_box.width-.3 and text_box.height<=block_box.height-1:
                            fitted=True;break
                    if not fitted:
                        label.remove()
                        small_exon_labels.append(exon['exon_label']+' ['+format(exon['aa_start0'],'g')+','+
                                                 format(exon['aa_end0'],'g')+')')
        # Essential classes are drawn last, so they lie above any annotation they overlap.
        for item in sorted(row['features'],key=lambda x:(1 if is_essential(x['feature']) else 0)):
            f=item['feature'];y=base+.17+.20*item['lane'];color=colors[f['display_id']]
            for lo,hi in item['segments']:
                ax.add_patch(Rectangle((lo-.5,y),hi-lo+1,.065,facecolor=color,edgecolor=color,linewidth=.65,
                                       zorder=3 if is_essential(f) else 2,gid='feature_'+row['side']+'_'+f['display_id']))
            center=(item['segments'][0][0]+item['segments'][-1][1])/2
            label=bar_label(f)+(' (partial)' if item['partial'] else '')
            ax.text(center,y+.085,label,ha='center',va='bottom',fontsize=8,color='#222222',
                    fontweight='bold' if is_essential(f) else 'normal',zorder=4 if is_essential(f) else 3,
                    gid='feature_label_'+row['side']+'_'+f['display_id'])
        for x,text in [(1,'1'),(n,str(n))]:
            ax.plot([x,x],[base-body_height/2,base-body_height/2-.05],color='black',linewidth=.65)
            ax.text(x,base-body_height/2-.08,text,ha='center',va='top',fontsize=8)
        ax.text(0,base-.42,reference['reference_protein'] if reference else iso,fontsize=7.5,ha='left',va='top',gid='native_id_'+row['side'])
        small_lines=(textwrap.wrap('Small coding segments (AA offsets): '+'; '.join(small_exon_labels),width=key_width)
                     if small_exon_labels else [])
        if len(small_lines)>len(small_key_reserve[row['side']]):
            raise ValueError('Small coding exon key exceeds reserved native layout space: '+str(iso))
        for k,line in enumerate(small_lines+unplaced_lines[row['side']]):
            ax.text(0,base-.57-.15*k,line,fontsize=8,color='#555555',va='top',gid='unplaced_'+row['side']+('' if k==0 else '_'+str(k+1)))
        if m.get('nmd') not in ('','Not-NMD',None):
            ax.text(length*1.02,base,m['nmd'],fontsize=8,ha='right',va='center',gid='nmd_'+row['side'])
        top-=row_height
    key_y=height-(.83 if any(e.get('exon_label') for i in (a,b) for e in models[i].get('coding_exons',[])) else .65)
    if partial_key:key_y=min(key_y,height-.99)
    for f,lines in key_lines:
        color=colors[f['display_id']]
        fig.add_artist(Rectangle((.70,(key_y-.07)/height),.013,.055/height,
                                 transform=fig.transSubfigure,facecolor=color,edgecolor=color))
        key_label=labels_by_key[feature_key(f)]
        fig.text(.725,key_y/height,key_label,fontsize=8.5,fontweight='bold',va='top',color='#222222',linespacing=1.1)
        fig.text(.81,key_y/height,'\n'.join(lines),fontsize=8.5,va='top',linespacing=1.1)
        key_y-=.16*max(len(lines),key_label.count('\n')+1)+.08
    # Check all labels, including native identifiers, rather than clipping them.
    fig.canvas.draw();renderer=fig.canvas.get_renderer();panel_box=fig.bbox
    labels=[*fig.texts,*ax.texts];layout=[]
    for obj in labels:
        box=obj.get_window_extent(renderer)
        if box.x0<panel_box.x0-1 or box.y0<panel_box.y0-1 or box.x1>panel_box.x1+1 or box.y1>panel_box.y1+1:
            raise ValueError('Compact figure text falls outside page: '+obj.get_text())
        layout.append(dict(role=obj.get_gid() or obj.get_text(),bounds_px=list(box.bounds)))
    for first,second in itertools.combinations(ax.texts,2):
        if first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer)):
            raise ValueError('Compact figure labels overlap: '+first.get_text()+' / '+second.get_text())
    for first,second in itertools.combinations(fig.texts,2):
        if first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer)):
            raise ValueError('Compact figure key labels overlap: '+first.get_text()+' / '+second.get_text())
    return fig,layout


def render_compact_figure(panels,sequences,models,context,panel_heights=None,figure_width=12):
    """Stack complete standard protein panels on one vector figure, without scaling.

    Reuse the same renderer and physical panel dimensions as individual exports.
    No PDF rasterization, source rescans, annotation selection or coordinate changes.
    """
    if not panels:raise ValueError('A combined protein figure requires at least one panel')
    plt=_mpl()
    if panel_heights is None:
        panel_heights=[]
        for panel in panels:
            preview,_=render_compact_pair(panel['pair'],panel['features'],panel['events'],sequences,models,context,figure_width=figure_width)
            panel_heights.append(preview.get_figheight());plt.close(preview)
    if len(panel_heights)!=len(panels) or any(not math.isfinite(h) or h<=0 for h in panel_heights):
        raise ValueError('Combined figure requires one positive height per complete panel')
    fig=plt.figure(figsize=(figure_width,sum(panel_heights)),facecolor='white')
    subfigures=fig.subfigures(len(panels),1,height_ratios=panel_heights,hspace=0,squeeze=False)
    layouts=[]
    for panel,subfigure in zip(panels,subfigures[:,0]):
        _,layout=render_compact_pair(panel['pair'],panel['features'],panel['events'],sequences,models,context,figure=subfigure,figure_width=figure_width)
        layouts.append(dict(gene_symbol=panel['pair']['gene_symbol'],bounds_px=list(subfigure.bbox.bounds),labels=layout))
    return fig,layouts


def run_protein_impact_plots(args):
    if not args.impact_results or not args.artifacts_dir:
        raise ValueError('protein-impact-plots requires --impact-results --artifacts-dir')
    out=Path(args.artifacts_dir);out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise ValueError('Use an empty figure output directory')
    manifest=dict(schema=SCHEMA,status='incomplete',generative_model_used=False,analysis_changed=False)
    dump(out/'figure_manifest.json',manifest)
    base,impact,models,sequences,pairs,feature_rows=load_audit(args.impact_results)
    pairs=selected_pairs(pairs,models,args.genes,args.pairs)
    types=tuple(args.feature_type or DEFAULT_TYPES)
    if len(types)!=len(set(types)) or not types:raise ValueError('Duplicate or empty feature-type selection')
    sources=tuple(args.annotation_source or PRIMARY_SOURCES)
    if len(sources)!=len(set(sources)) or not sources:raise ValueError('Duplicate or empty annotation-source selection')
    available_sources={r['provider'] for r in feature_rows}
    if set(sources)-available_sources:raise ValueError('Requested annotation source not present in the pinned audit')
    feature_references=feature_reference_selection(getattr(args,'feature_reference',None),pairs,feature_rows)
    essential=load_essential_whitelist(getattr(args,'essential_whitelist',None))
    context,units=donor_context(base,impact)
    model_index={r['isoform_id']:r for r in models}
    style=getattr(args,'protein_plot_style','compact')
    figure_layout=getattr(args,'protein_figure_layout','pages')
    figure_width=getattr(args,'protein_figure_width',12)
    if not math.isfinite(figure_width) or figure_width<=0:raise ValueError('Protein figure width must be positive and finite')
    if figure_layout=='stacked' and style!='compact':
        raise ValueError('--protein-figure-layout stacked requires --protein-plot-style compact')
    exon_root=None
    if getattr(args,'coding_exon_map',None):
        if style!='compact':raise ValueError('--coding-exon-map requires the compact protein architecture view')
        from ..coding_exons import load_coding_exon_map
        exon_models,exon_root,exon_manifest=load_coding_exon_map(args.coding_exon_map,base,models)
        model_index={k:dict(v,coding_exons=exon_models[k]['segments']) for k,v in model_index.items()}
    plans=[];all_audit=[];all_events=[];flat=[];names=defaultdict(int)
    for pair in pairs:
        feats,audit,events=feature_plan(pair,feature_rows,sequences,types,sources,feature_references.get(pair['gene_symbol']),essential)
        names[pair['gene_symbol']]+=1
        stem=re.sub(r'[^A-Za-z0-9_.-]','_',pair['gene_symbol'])+'_pair'+str(names[pair['gene_symbol']]).zfill(2)+'_protein_features'
        plans.append(dict(pair=pair,features=feats,events=events,stem=stem))
        all_audit.extend(audit)
        all_events.extend(dict(gene_symbol=pair['gene_symbol'],isoform_A=pair['isoform_A'],isoform_B=pair['isoform_B'],**e) for e in events)
        for f in feats:
            flat.append({k:json.dumps(v,separators=(',',':')) if isinstance(v,(list,dict)) else v for k,v in f.items()})
    write(out/'annotation_display_audit.tsv',all_audit)
    if any(not p['features'] for p in plans):raise ValueError('A selected pair has no annotations in the declared view; see complete display audit')
    write(out/'alignment_events.tsv',all_events)
    write(out/'plotted_feature_coordinates.tsv',flat)
    dump(out/'plot_data.json',plans)
    if style=='compact':
        native_rows=[]
        for plan in plans:
            for side in ('A','B'):
                placed,unplaced=architecture_features(plan['features'],side)
                index={v['feature']['display_id']:v for v in placed}
                for f in plan['features']:
                    item=index.get(f['display_id']);segments=item['segments'] if item else []
                    native_rows.append(dict(gene_symbol=plan['pair']['gene_symbol'],side=side,
                        isoform_id=plan['pair']['isoform_'+side],feature_id=f['feature_id'],
                        feature_type=f['feature_type'],description=f['description'],provider=f['provider'],
                        native_start=segments[0][0] if segments else '',native_end=segments[-1][1] if segments else '',
                        placement='unplaced' if f['display_id'] in unplaced else 'partial' if item and item['partial'] else 'mapped' if item else 'absent'))
        write(out/'native_architecture_coordinates.tsv',native_rows)
        if exon_root:
            selected={p['pair']['isoform_'+s] for p in plans for s in ('A','B')}
            write(out/'coding_exon_coordinates.tsv',[r for r in table(exon_root/'coding_exon_coordinates.tsv') if r['isoform_id'] in selected])
    plt=_mpl()
    render=render_compact_pair if style=='compact' else render_pair
    from matplotlib.backends.backend_pdf import PdfPages
    epoch=datetime(2000,1,1,tzinfo=timezone.utc)
    pdf_metadata=dict(Title='Native protein domain and motif comparisons',Author='AltAnalyze3 iso2function',Creator=SCHEMA,
                      CreationDate=epoch,ModDate=epoch)
    layouts=[];combined_panels=[];panel_heights=[]
    with PdfPages(out/'protein_feature_comparisons.pdf',metadata=pdf_metadata) as combined:
        for plan in plans:
            chunks=[plan['features'][i:i+12] for i in range(0,len(plan['features']),12)]
            for i,chunk in enumerate(chunks,1):
                stem=plan['stem']+('_page'+str(i) if len(chunks)>1 else '')
                fig,layout=render(plan['pair'],chunk,plan['events'],sequences,model_index,context,i,len(chunks),
                                  **({'figure_width':figure_width} if style=='compact' else {}))
                if figure_layout=='pages':combined.savefig(fig)
                else:
                    combined_panels.append(dict(pair=plan['pair'],features=chunk,events=plan['events']))
                    panel_heights.append(fig.get_figheight())
                fig.savefig(out/(stem+'.pdf'),metadata=pdf_metadata)
                fig.savefig(out/(stem+'.svg'),metadata={'Date':'2000-01-01','Creator':SCHEMA})
                fig.savefig(out/(stem+'.png'),dpi=160,metadata={'Software':SCHEMA})
                layouts.append(dict(stem=stem,labels=layout))
                plt.close(fig)
        if figure_layout=='stacked':
            fig,combined_layout=render_compact_figure(combined_panels,sequences,model_index,context,panel_heights,figure_width)
            combined.savefig(fig)
            fig.savefig(out/'protein_feature_comparisons.svg',metadata={'Date':'2000-01-01','Creator':SCHEMA})
            fig.savefig(out/'protein_feature_comparisons.png',dpi=160,metadata={'Software':SCHEMA})
            dump(out/'combined_layout_audit.json',combined_layout)
            plt.close(fig)
    dump(out/'label_layout_audit.json',layouts)
    import Bio,matplotlib,platform
    from matplotlib import font_manager,ft2font
    font_path=font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans'))
    manifest.update(status='complete',plot_style=style,figure_layout=figure_layout,figure_width_inches=figure_width,pairs=len(plans),
                    panels=len(layouts),pages=1 if figure_layout=='stacked' else len(layouts),genes=[p['pair']['gene_symbol'] for p in plans],
                    coding_exon_map=(dict(path=str(exon_root),manifest_sha256=sha(exon_root/'exon_map_manifest.json'),
                                         source_manifest_sha256=exon_manifest['source_manifest_sha256']) if exon_root else None),
                    source_manifest_sha256=sha(base/'impact_manifest.json'),
                    inputs={str(base/p):sha(base/p) for p in ['impact_manifest.json','native_model_audit.tsv','pair_summary.tsv','feature_projections.tsv','native_proteins.fasta']},
                    code_sha256={str(Path(__file__)):sha(__file__),str(Path(__file__).parents[1]/'protein_isoform_audit.py'):sha(Path(__file__).parents[1]/'protein_isoform_audit.py')},
                    selection=dict(feature_types=list(types),sources=list(sources),functional_annotation_rows=len(flat),
                                   feature_references=feature_references,
                                   selected_source_rows=len(all_audit),all_other_annotations_preserved_in_display_audit=True,
                                   essential_whitelist=dict(path=essential['path'],sha256=essential['sha256'],rules=len(essential['rules'])),
                                   named_UniProt_domains_always_displayed=True,
                                   named_domain_labels='Literal source labels unless a pinned specific whitelist rule supplies the class name',
                                   essential_rows_A=sum(1 for e in all_audit if e['essential_class'] and e['source_isoform']==e['isoform_A']),
                                   essential_rows_A_shown=sum(1 for e in all_audit if e['essential_class'] and e['source_isoform']==e['isoform_A']
                                                              and (e['display_role'].startswith('displayed') or e['display_role'].startswith('coalesced'))),
                                   essential_rows_A_hidden_placeable=sum(1 for e in all_audit if e['essential_class'] and e['source_isoform']==e['isoform_A']
                                                              and e['source_status']=='intact_sequence'
                                                              and not (e['display_role'].startswith('displayed') or e['display_role'].startswith('coalesced'))),
                                   essential_classes_displayed=sorted({f['essential_class'] for p in plans for f in p['features'] if f.get('essential_class')})),
                    coordinate_convention=('1-based inclusive native protein coordinates on each compact track; alignment coordinates in audit tables'
                                           if style=='compact' else '1-based inclusive, native A aligned coordinates; native B positions explicitly tabulated'),
                    RNA_support=units,projection_disagreements=sum(f['projection_disagreement'] for p in plans for f in p['features']),
                    software=dict(python=platform.python_version(),biopython=Bio.__version__,matplotlib=matplotlib.__version__,
                                  freetype=ft2font.__freetype_version__,font_file=font_path,font_sha256=sha(font_path)),
                    presentation=dict(amino_acid_substitutions_highlighted=style=='detailed',
                                      narrative_footer=style=='detailed',
                                      domain_colors='annotation identity' if style=='compact' else 'residue alignment state'),
                    limits=['These are sequence-based predictions, not direct native protein assays.',
                            'Missing full feature sequence is not a measurement of lost protein activity.',
                            'Inserted B residues are counted separately and not stretched onto A coordinates.',
                            'Repeated or conflicting placements are uncertain, never silently called absent.'],
                    outputs={p.name:sha(p) for p in sorted(out.iterdir()) if p.name!='figure_manifest.json'})
    dump(out/'figure_manifest.json',manifest)
    return 0
