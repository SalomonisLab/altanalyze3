"""Compare retained MarkerFinder scores across modalities without rerunning analysis."""
from functools import lru_cache
from pathlib import Path
import re

import numpy as np
import pandas as pd


def requested(question):
    q = question.lower()
    if re.search(r'correlat|discordan|concordan|\bnetwork\b|\bpathway\b', q):
        return False
    return bool(re.search(r'\bmarkers?\b', q) and re.search(
        r'modalit|multi[- ]?omic|\badts?\b|lipid|metabolite|tf activity|transcription factor', q))


@lru_cache(maxsize=24)
def _read_table(path, mtime_ns, size):
    # File signatures invalidate cached tables when a saved analysis is replaced.
    return pd.read_csv(path, sep='\t')


def marker_tables(meta):
    from .app import _modality_marker_analysis
    available = {m['id']: m.get('label', m['id']) for m in (meta.get('modalities') or {}).get('available', [])}
    for mod in meta.get('marker_analysis_by_modality') or {}:
        available.setdefault(mod, mod)
    available.setdefault('rna', 'RNA')
    frames, coverage = {}, []
    for mod, label in available.items():
        analysis = _modality_marker_analysis(meta, mod)
        paths = []
        primary = analysis.get('markers_tsv')
        redundant = analysis.get('redundant_markers_tsv')
        if primary:
            p = Path(primary)
            redundant = redundant or p.with_name(p.stem.replace('_markers', '_redundant_markers') + '.tsv')
        for value in (redundant, primary):
            if value and Path(value).is_file() and str(value) not in paths:
                paths.append(str(value))
        parts = []
        for path in paths:
            stat = Path(path).stat()
            frame = _read_table(path, stat.st_mtime_ns, stat.st_size)
            if {'Gene', 'cluster', 'rho'}.issubset(frame.columns):
                parts.append(frame)
        if parts:
            frame = pd.concat(parts, ignore_index=True).copy()
            frame['rho'] = pd.to_numeric(frame['rho'], errors='coerce')
            frame = frame.loc[np.isfinite(frame.rho)]
            frame['Gene'] = frame.Gene.astype(str)
            frame['cluster'] = frame.cluster.astype(str)
            frames[mod] = frame.drop_duplicates(['cluster', 'Gene'], keep='first')
        coverage.append(dict(modality=mod, label=label, source_files=paths,
                             status='available' if mod in frames else 'no retained MarkerFinder correlation scores'))
    return frames, coverage


def _state_matches(question, states):
    normalize = lambda x: re.sub(r'[-_\s]+', ' ', str(x).lower()).strip()
    q = normalize(question)
    matches = []
    for state in sorted(states, key=len, reverse=True):
        for match in re.finditer(r'(?<!\w)' + re.escape(normalize(state)) + r'(?!\w)', q):
            if not any(a <= match.start() and match.end() <= b for a, b, _ in matches):
                matches.append((match.start(), match.end(), state))
    return list(dict.fromkeys(s for _, _, s in sorted(matches)))


def answer_if_requested(app, meta, question):
    if not requested(question):
        return None
    frames, coverage = marker_tables(meta)
    states = sorted({s for frame in frames.values() for s in frame.cluster.unique()})
    names = _state_matches(question, states)
    result = dict(question=question, intent='modality_markers', reading=dict(intent='modality_markers', source='saved markers'))
    if not frames:
        return dict(result, status='not_covered', answer='No retained MarkerFinder correlation scores are available for this dataset.', coverage=coverage)
    if len(names) != 1:
        return dict(result, status='clarify', answer='Specify one cell state present in the saved marker tables.', choices=dict(states=states))
    state = names[0]
    result['reading']['cell_state'] = state
    q = question.lower()
    requested_ids = set()
    if not re.search(r'modalit|multi[- ]?omic', q):
        for pattern, mod in [(r'\brna\b', 'rna'), (r'\badts?\b', 'adt'), (r'lipid', 'lipid'),
                             (r'metabolite', 'metabolite'), (r'tf activity|transcription factor', 'grn_tf')]:
            if re.search(pattern, q): requested_ids.add(mod)
    named_only = bool(re.search(r'\b(named|identified|known)\b', q))
    rows = []
    labels = {item['modality']: item['label'] for item in coverage}
    for mod, frame in frames.items():
        if requested_ids and mod not in requested_ids:
            continue
        subset = frame.loc[(frame.cluster == state) & (frame.rho > 0)].sort_values(['rho', 'Gene'], ascending=[False, True])
        for _, row in subset.iterrows():
            unidentified = bool(re.match(r'^unknown', row.Gene, re.I))
            if named_only and unidentified:
                continue
            def number(column):
                value = pd.to_numeric(row.get(column), errors='coerce')
                return float(value) if pd.notna(value) and np.isfinite(value) else None
            rows.append(dict(modality=labels[mod], modality_id=mod, feature=row.Gene,
                             marker_r=float(row.rho), annotation='Unidentified' if unidentified else 'Named',
                             state_mean=number('Query Exp'), other_states_mean=number('Ref Exp')))
    rows.sort(key=lambda r: (-r['marker_r'], r['modality'], r['feature']))
    best_by_modality = {}
    for row in rows:
        best_by_modality.setdefault(row['modality_id'], row)
    for item in coverage:
        item['n_positive_markers'] = sum(row['modality_id'] == item['modality'] for row in rows)
    if not rows:
        return dict(result, status='not_covered', answer=f'No positive {"named " if named_only else ""}markers with retained correlation scores were found for {state} in the requested modalities.', coverage=coverage)
    best = rows[0]
    leaders = '; '.join(f'{r["modality"]}: {r["feature"]} (r={r["marker_r"]:.4f})' for r in best_by_modality.values())
    missing = [item['label'] for item in coverage if not item['n_positive_markers'] and (not requested_ids or item['modality'] in requested_ids)]
    answer = (f'Among the retained marker results for {state}, the strongest {"named " if named_only else ""}marker is '
              f'{best["feature"]} ({best["modality"]}, MarkerFinder r={best["marker_r"]:.4f}). '
              f'Best per modality: {leaders}. '
              'Ranked by positive Pearson correlation with cell-state membership (that state versus all others). '
              'Means retain each modality’s own scale. This is a marker-specificity ranking, not independent validation of imputed modalities. '
              f'The table includes all {len(rows)} retained positive markers; it is not a new genome-wide marker scan.')
    if any(r['annotation'] == 'Unidentified' for r in best_by_modality.values()):
        answer += ' Unknown-labelled features remain unidentified; search Named to view the named features.'
    if missing:
        answer += ' No comparable retained positive marker scores for: ' + ', '.join(missing) + '.'
    return dict(result, status='ok', answer=answer, coverage=coverage, best=best, best_by_modality=list(best_by_modality.values()),
                table=dict(columns=['modality', 'feature', 'marker_r', 'annotation', 'state_mean', 'other_states_mean'], rows=rows),
                result_controls=dict(sort_by='marker_r', sort_direction='desc'),
                plot=dict(kind='modality_markers', cell_state=state),
                provenance=dict(source='retained marker tables', statistic='MarkerFinder Pearson r with 0/1 cell-state membership', differential_required=False),
                follow_ups=[f'What is the best named modality marker of {state}?'])
