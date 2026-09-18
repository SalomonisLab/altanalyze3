"""Cross-modality pathway representation from retained markers or matched contrasts.

This is descriptive coverage, not an enrichment test. Counts retain assay feature
identity; duplicate pathway nodes never inflate them. No differential is rerun.
"""
from functools import lru_cache
import math
import re

from ..discover_integration import _diagrams, _diagram_synonyms, _best_class, _best_group
from .modality_markers import marker_tables, _state_matches
from .cross_modal import adt_map, _adt_key

COLORS = dict(rna='#77aadd', adt='#ee8866', metabolite='#ddcc77', lipid='#88ccaa',
              lipids='#88ccaa', grn_tf='#bb99cc', grn='#aaaadd')
SUPPORTED = tuple(COLORS)


def requested(question):
    return bool(re.search(r'\bpathways?\b', question, re.I) and re.search(
        r'cross[- ]?modal|multi[- ]?modal|multi[- ]?omic|modalities|representation', question, re.I))


def norm(value):
    # Preserve stereochemical and isomer distinctions; do not strip punctuation.
    return re.sub(r'\s+', ' ', str(value or '').strip()).casefold()


@lru_cache(maxsize=1)
def pathway_index():
    result = {}
    synonyms = _diagram_synonyms()
    for pid, diagram in _diagrams()['diagrams'].items():
        genes, metabolites, classes = {}, {}, {}
        for i, node in enumerate(diagram['nodes']):
            if node['type'] in ('GeneProduct', 'Protein', 'Rna'):
                for name in (node['label'], node.get('ensembl')):
                    if name: genes.setdefault(norm(name), set()).add(i)
            elif node['type'] == 'Metabolite':
                if node['label']: metabolites.setdefault(norm(node['label']), set()).add(i)
                cls = _best_class(node['label'], synonyms)
                group = _best_group(node['label'])[1]
                for c in ([cls] if cls else group):
                    classes.setdefault(c.upper(), set()).add(i)
        result[pid] = (genes, metabolites, classes)
    return result


@lru_cache(maxsize=4)
def antibody_mapping(species):
    return adt_map(species)


def feature_nodes(feature, modality, index, species='human'):
    genes, metabolites, classes = index
    if modality in ('rna', 'grn_tf'):
        return genes.get(norm(feature), set())
    if modality == 'adt':
        key = _adt_key(feature)
        mapping = antibody_mapping(species)
        partners = mapping.get(key, mapping.get(key.split('_')[0], []))
        return set().union(*(genes.get(norm(g), set()) for g in partners))
    if modality == 'grn':
        parts = [genes.get(norm(g), set()) for g in feature.split('|') if g]
        # A composite edge counts once, only if every regulator and target maps.
        return set().union(*parts) if len(parts) > 1 and all(parts) else set()
    if modality in ('lipid', 'lipids'):
        head = re.match(r'^([A-Za-z0-9-]+)', feature)
        cls = head.group(1).upper() if head else ''
        if cls not in classes:
            cls = _best_class(feature, _diagram_synonyms()).upper()
        return classes.get(cls, set())
    if modality == 'metabolite' and not norm(feature).startswith('unknown'):
        return metabolites.get(norm(feature), set())
    return set()


def rank(features, labels, species='human'):
    """features[modality][original feature ID] = statistic dictionary."""
    mods = [m for m in SUPPORTED if m in labels]
    rows, mapped = [], {m: set() for m in mods}
    diagrams = _diagrams()['diagrams']
    for pid, index in pathway_index().items():
        hits = {m: sorted(f for f in features.get(m, {}) if feature_nodes(f, m, index, species)) for m in mods}
        for m in mods: mapped[m].update(hits[m])
        row = dict(id=pid, pathway=diagrams[pid]['name'], hits=hits,
                   modality_count=sum(bool(hits[m]) for m in mods))
        row.update({m: len(hits[m]) for m in mods})
        rows.append(row)
    denominators = {m: max((r[m] for r in rows), default=0) for m in mods}
    # Modalities with no mapped features cannot contribute; report them explicitly.
    evaluated = [m for m in mods if denominators[m]]
    for row in rows:
        row['balanced_coverage'] = sum(row[m] / denominators[m] for m in evaluated) / max(1, len(evaluated))
        row['combined_score'] = row['modality_count'] + row['balanced_coverage']
        row['total_hits'] = sum(row[m] for m in mods)
    rows = [r for r in rows if r['modality_count']]
    rows.sort(key=lambda r: (-r['combined_score'], r['pathway'], r['id']))
    return dict(rows=rows, modalities=[dict(id=m, label=labels[m], color=COLORS[m],
                  mapped_features=len(mapped[m]), input_features=len(features.get(m, {})),
                  normalization_max=denominators[m]) for m in mods], evaluated_pathways=len(diagrams))


def collect(app, meta, state, source='marker', contrast=''):
    labels = {m['id']: m.get('label', m['id']) for m in (meta.get('modalities') or {}).get('available', []) if m['id'] in SUPPORTED}
    labels.setdefault('rna', 'RNA')
    features, coverage, comparison = {}, [], ''
    if source == 'marker':
        frames, coverage = marker_tables(meta)
        for mod, frame in frames.items():
            if mod not in SUPPORTED: continue
            labels.setdefault(mod, mod)
            subset = frame.loc[(frame.cluster == state) & (frame.rho > 0)]
            features[mod] = {r.Gene: dict(statistic='MarkerFinder r', value=float(r.rho)) for r in subset.itertuples()}
    else:
        from .integration_data import integration_data
        ds = integration_data(app, meta)
        contrast = contrast or ds.current_contrast
        # Validate a run identity before asking the adapter for same-contrast siblings.
        manifest = ds.deg_manifest().get('comparisons', []) if hasattr(ds, 'deg_manifest') else ds.ds.deg_manifest().get('comparisons', [])
        entry = next((c for c in manifest if c.get('id') == contrast), None)
        if not entry:
            raise ValueError('Choose a completed comparison before requesting cross-modality contrast pathways.')
        comparison = entry.get('comparison') or entry.get('label') or contrast
        for mod in labels:
            available = ds.available(contrast, mod, state)
            coverage.append(dict(modality=mod, status='available' if available else 'no matching completed comparison for this cell state'))
            if not available: continue
            features[mod] = {}
            for row in ds.differentials(contrast, state, mod):
                feature = row.get('feature_key') or row.get('gene')
                value = row.get('log2fc')
                if feature and value is not None and math.isfinite(float(value)):
                    features[mod][str(feature)] = dict(statistic='reported log2FC', value=float(value), fdr=row.get('fdr'))
    result = rank(features, labels, meta.get('species', 'human'))
    result.update(features=features, coverage=coverage, cell_state=state, source=source, species=meta.get('species', 'human'),
                  contrast=contrast if source == 'differential' else '', comparison=comparison)
    return result


def diagram(result, pid):
    row = next((r for r in result['rows'] if r['id'] == pid), None)
    if not row: return dict(available=False, nodes=[], note='No retained features map to this pathway in the selected context.')
    raw = _diagrams()['diagrams'][pid]
    nodes = [dict(n, hits=[]) for n in raw['nodes']]
    for mod, features in row['hits'].items():
        for feature in features:
            for i in feature_nodes(feature, mod, pathway_index()[pid], result.get('species', 'human')):
                nodes[i]['hits'].append(dict(modality=mod, feature=feature, **result['features'][mod][feature]))
    modalities = []
    for mod in result['modalities']:
        hits = [result['features'][mod['id']][feature] for feature in row['hits'][mod['id']]]
        values = [h['value'] for h in hits if math.isfinite(h['value'])]
        span = max(map(abs, values), default=0.)
        modalities.append(dict(mod, scale=dict(
            available=bool(values), minimum=-span if result['source'] == 'differential' else 0.,
            maximum=span, statistic='log2FC' if result['source'] == 'differential' else 'MarkerFinder r',
            n_features=len(values))))
    for node in nodes:
        # Keep modalities separate. For multiple species/genes on one node, color
        # the strongest retained score; preserve every individual score in hover.
        node['modality_values'] = {}
        for hit in node['hits']:
            old = node['modality_values'].get(hit['modality'])
            if old is None or abs(hit['value']) > abs(old['value']):
                node['modality_values'][hit['modality']] = hit
        if node['hits']:
            node.update(link_feature=node['hits'][0]['feature'], link_modality=node['hits'][0]['modality'])
    return dict(raw, nodes=nodes, available=True, cross_modal=True, modalities=modalities,
                cell_state=result['cell_state'], comparison=result['comparison'], source=result['source'],
                note='Unique feature support by modality; repeated diagram nodes are counted once per pathway. '
                     'Each modality has its own numeric color scale. Each node stripe uses the strongest absolute retained score in that modality; hover for all feature scores.')


def answer_if_requested(app, meta, question):
    if not requested(question): return None
    source = 'differential' if re.search(r'contrast|compariso|differential|\bversus\b|\bvs\b|regulated|changed', question, re.I) else 'marker'
    result = dict(intent='cross_modality_pathways', question=question)
    frames, _ = marker_tables(meta)
    states = sorted({s for f in frames.values() for s in f.cluster.unique()})
    if not states or source == 'differential':
        from .integration_data import integration_data
        states = list(integration_data(app, meta).states)
    names = _state_matches(question, states)
    if len(names) != 1:
        return dict(result, status='clarify', answer='Specify one cell state for the cross-modality pathway analysis.', choices=dict(states=states))
    contrast = (meta.get('differential') or {}).get('run_id', '')
    if source == 'differential':
        from .integration_data import integration_data
        adapter = integration_data(app, meta)
        manifest_ds = getattr(adapter, 'ds', adapter)
        entries = manifest_ds.deg_manifest().get('comparisons', [])
        named = [c for c in entries if norm(c.get('comparison')) and norm(c['comparison']) in norm(question)]
        if named:
            contrast = named[0]['id']
        elif re.search(r'\bversus\b|\bvs\b', question, re.I):
            return dict(result, status='clarify', answer='Select a saved comparison or use its complete comparison label.',
                        choices=dict(comparisons=list(dict.fromkeys(c.get('comparison', c['id']) for c in entries))))
    try:
        data = collect(app, meta, names[0], source, contrast)
    except ValueError as exc:
        return dict(result, status='not_covered', answer=str(exc))
    if not data['rows']:
        return dict(result, status='not_covered', answer=(
            f"No retained {'differential calls' if source == 'differential' else 'positive markers'} map to the bundled pathways for {names[0]}"
            + (f" in {data['comparison']}" if data['comparison'] else '') + '. No results from another cell state or comparison were substituted.'),
            coverage=data['coverage'], reading=dict(source=source, cell_state=names[0], contrast=data['contrast']))
    mods = data['modalities']
    missing = [m['label'] for m in mods if not m['mapped_features']]
    answer = (f"{len(data['rows'])} pathways contain retained {'positive cell-state markers' if source == 'marker' else 'differential calls'} "
              f"for {names[0]}{(' — ' + data['comparison']) if data['comparison'] else ''}; {data['evaluated_pathways']} bundled WikiPathways diagrams evaluated. "
              'Ranked by modality breadth, then equally weighted coverage: each count is divided by that modality’s maximum pathway count across the full evaluated library; those values are averaged. '
              'Combined score = modality count + balanced coverage. Filters keep these denominators fixed. '
              'Counts are unique assay features within each modality; lipid species count individually. '
              'Checked modalities are all required to have at least one hit. This is representation, not statistical enrichment. '
              'Unidentified metabolites remain unmapped. GRN edges require both the regulators and target in the pathway.')
    if missing: answer += ' No mapped retained hits for: ' + ', '.join(missing) + '.'
    return dict(result, status='ok', answer=answer,
                reading=dict(intent='cross_modality_pathways', source=source, cell_state=names[0]),
                table=dict(columns=['pathway', 'modality_count', *[m['id'] for m in mods], 'balanced_coverage', 'combined_score'], rows=data['rows']),
                column_labels={**{m['id']:m['label'] for m in mods}, 'modality_count':'Modalities', 'balanced_coverage':'Balanced coverage', 'combined_score':'Combined score'},
                result_controls=dict(sort_by='combined_score', sort_direction='desc'),
                plot=dict(kind='cross_pathways', cell_state=names[0], source=source, contrast=data['contrast'], modalities=mods),
                provenance=dict(evaluated_pathways=data['evaluated_pathways'], coverage=data['coverage'],
                                normalization='count / maximum count per modality across all evaluated pathways; mean over modalities with mapped hits',
                                mapping='Exact gene/Ensembl or metabolite name; curated ADT-to-gene; Discover lipid class; complete GRN endpoints'))
