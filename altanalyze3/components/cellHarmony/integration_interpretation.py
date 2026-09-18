"""Describe the exact displayed evidence, following Discover's lookup-based Interpret.

No model-generated statistics or unrelated comparison lookups. Pathway regions
come from the published diagram's labels, not inferred reaction assignments.
"""
import math
from collections import Counter, defaultdict


def _distinct(nodes):
    best = {}
    for node in nodes:
        value = node.get('log2fc')
        name = str(node.get('label') or node.get('id') or '').strip()
        if value is None or not name:
            continue
        key = name.lower().rstrip('s')
        if key not in best or abs(value) > abs(best[key]['log2fc']):
            best[key] = node
    return sorted(best.values(), key=lambda n: -abs(n['log2fc']))


def _changes(nodes, label):
    rows = _distinct(nodes)
    if not rows:
        return f'No {label} on this figure has a reported differential result.'
    parts = []
    for direction, sign in [('increased', 1), ('decreased', -1)]:
        selected = [n for n in rows if n['log2fc'] * sign > 0]
        if selected:
            names = ', '.join(f"{n.get('label') or n['id']} ({n['log2fc']:+.2f})" for n in selected[:3])
            parts.append(f'{len(selected)} {direction}, including {names}')
    return f"Among the {len(rows)} reported {label}, " + '; '.join(parts) + ' (log2 fold change).'


def interpret_network(result):
    edges = result.get('edges') or []
    if not edges:
        return result.get('note') or 'No edges satisfy the selected criteria; no network is drawn.'
    nodes = result.get('nodes') or []
    factors = {e['source'] for e in edges}
    targets = {e['target'] for e in edges}
    lines = [f"In {result.get('cell_state', '')}, this network contains {len(edges)} retained differential edges "
             f"from {len(factors)} factors to {len(targets)} targets for {result.get('comparison') or result.get('contrast', '')}."]
    counts = Counter(e['source'] for e in edges)
    lines.append('Factors with the most displayed targets: ' + ', '.join(
        f'{name} ({count})' for name, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:4]) + '.')
    up = sum(e.get('log2fc') is not None and e['log2fc'] > 0 for e in edges)
    down = sum(e.get('log2fc') is not None and e['log2fc'] < 0 for e in edges)
    lines.append(f'{up} edge scores increased and {down} decreased in the case arm relative to the control arm.')
    strongest = sorted([e for e in edges if e.get('log2fc') is not None], key=lambda e: -abs(e['log2fc']))[:3]
    if strongest:
        lines.append('Largest displayed edge changes: ' + '; '.join(
            f"{e['source']} → {e['target']} ({e['log2fc']:+.2f} log2 fold change"
            + (f", {e.get('significance', 'FDR')} {e['fdr']:.3g}" if e.get('fdr') is not None else '') + ')'
            for e in strongest) + '.')
    lines.append(_changes([n for n in nodes if n['id'] in targets], 'target genes'))
    lines.append(_changes([n for n in nodes if n['id'] in factors], 'TF expression results'))
    activity = [dict(n, log2fc=n.get('activity_log2fc')) for n in nodes if n['id'] in factors]
    lines.append(_changes(activity, 'TF activity results'))
    lines.append('Node colour uses each gene’s own expression change, with TF activity as a fallback. '
                 'An edge change describes predicted regulatory activity; its sign does not establish activation or repression of the target.')
    return ' '.join(lines)


def interpret_pathway(result):
    nodes = result.get('nodes') or []
    genes = [n for n in nodes if n.get('measured_as') == 'gene']
    metabolites = [n for n in nodes if str(n.get('measured_as', '')).startswith(('lipid', 'metabolite'))]
    lipid = any(str(n.get('measured_as', '')).startswith('lipid') for n in metabolites)
    lines = [f"{result.get('name', '')} in {result.get('cell_state', '')}, {result.get('contrast', '')}. "
             f"The dataset measures {result.get('n_in_atlas', 0)} of this diagram’s {len(nodes)} nodes. "
             'The overlay uses the stored differential calls without an additional significance filter.']
    lines.append(_changes(genes, 'genes'))
    lines.append(_changes(metabolites, 'lipid classes/groups' if lipid else 'metabolites'))
    regions = defaultdict(list)
    for node in _distinct(genes + metabolites):
        cx = float(node.get('x') or 0) + float(node.get('w') or 0) / 2
        cy = float(node.get('y') or 0) + float(node.get('h') or 0) / 2
        closest, distance = '', 260.
        for label in result.get('labels') or []:
            text = str(label.get('text') or '').strip()
            if len(text) < 4:
                continue
            delta = math.hypot(cx - float(label.get('x') or 0), cy - float(label.get('y') or 0))
            if delta < distance:
                closest, distance = text, delta
        if closest:
            regions[closest].append(node['log2fc'])
    if regions:
        lines.append('By nearest region label on the published diagram: ' + '; '.join(
            f"{name} ({sum(v > 0 for v in values)} increased, {sum(v < 0 for v in values)} decreased)"
            for name, values in sorted(regions.items(), key=lambda item: -len(item[1]))[:4]) + '.')
    absent = sum(n.get('in_atlas') and n.get('log2fc') is None for n in nodes)
    lines.append(f'{absent} measured nodes have no reported differential call; an uncoloured node is not evidence of unchanged abundance.')
    if lipid:
        lines.append('Lipid class colours average the regulated species in each class; gene colours use their own log2 fold changes. These scales are not interchangeable.')
    lines.append('These changes locate the measured response on the pathway; they do not establish metabolic flux.')
    return ' '.join(lines)
