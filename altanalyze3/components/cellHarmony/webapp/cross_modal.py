"""Descriptive, matched Expression correlations; no differential run is required."""
from __future__ import annotations

from collections import OrderedDict
import re
from threading import RLock

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

_CACHE = OrderedDict()
_LOCK = RLock()


def requested_modality(question):
    q = question.lower()
    if not re.search(r'correlat|discordan|concordan|agree|disagree|near.zero', q):
        return None
    if re.search(r'\badts?\b|antibody|surface.*gene|protein.*expression', q):
        return 'adt'
    if re.search(r'\btfs?\b|transcription factor|tf activity', q):
        return 'grn_tf'
    return None


def examples(modalities):
    ids = set(modalities)
    result = []
    if 'grn_tf' in ids:
        result += ['Which TFs have activity discordant with gene expression across cell states?',
                   'Correlate TF activity with matching gene expression across cell states',
                   'Which TFs have activity uncorrelated with gene expression across cell states?']
    if 'adt' in ids:
        result += ['Correlate ADT abundance with its cell-surface gene expression across cell states']
    return result


def _adt_key(name):
    return re.sub(r'^(hu|mm|ms|mouse)\.', '', str(name), flags=re.I).casefold()


def adt_map(species='human'):
    if str(species).lower() in ('mouse', 'mus musculus', 'mm'):
        from ...rna2adt.mouse.adt_mgi_map import load_curated_adt_mgi_map
        mapping = load_curated_adt_mgi_map()
    else:
        from ...rna2adt.adt_rna_map import load_curated_adt_rna_map
        from ...rna2adt.lung.adt_rna_map import load_curated_adt_rna_map as lung_map
        mapping = {**lung_map(), **load_curated_adt_rna_map()}
    return {_adt_key(key): list(value) for key, value in mapping.items()}


def paired_features(rna, other, modality, species='human'):
    aliases = {}
    names = list(map(str, rna['var_names']))
    for i, name in enumerate(names):
        aliases.setdefault(name.casefold(), set()).add(i)
    var = getattr(rna['adata'], 'var', pd.DataFrame())
    for column in ('gene_symbols', 'gene_symbol', 'symbol', 'features', 'gene_name'):
        if column in var:
            for i, value in enumerate(var[column]):
                if pd.notna(value):
                    aliases.setdefault(str(value).casefold(), set()).add(i)
    mapping = adt_map(species) if modality == 'adt' else {}
    other_var = getattr(other['adata'], 'var', pd.DataFrame())
    pairs, skipped = [], []
    for j, feature in enumerate(map(str, other['var_names'])):
        if modality == 'adt':
            candidates = [feature]
            if 'raw_feature_name' in other_var:
                candidates.insert(0, str(other_var.iloc[j]['raw_feature_name']))
            genes = None
            for candidate in candidates:
                key = _adt_key(candidate)
                if key in mapping:
                    genes = mapping[key]
                    break
            if genes is None:
                # A clone suffix may be removed only to find an existing curated entry.
                genes = next((mapping[_adt_key(c).split('_')[0]] for c in candidates
                              if _adt_key(c).split('_')[0] in mapping), [])
        else:
            genes = [feature]
        found = False
        for gene in genes:
            indices = aliases.get(str(gene).casefold(), set())
            if len(indices) != 1:
                continue
            i = next(iter(indices))
            pairs.append((feature, str(gene), i, j))
            found = True
        if not found:
            skipped.append(feature)
    return pairs, skipped


def _group_means(cache, rows, columns, grouper):
    """Read only matched columns, in bounded chunks; support published bundles."""
    adata = cache['adata']
    output = np.empty((grouper.shape[0], len(columns)), dtype=np.float64)
    for start in range(0, len(columns), 32):
        cols = columns[start:start + 32]
        if hasattr(adata, 'X'):
            values = adata[:, cols].X
            values = values[rows].astype(np.float64)
        else:
            values = np.column_stack([np.asarray(adata[:, str(cache['var_names'][c])].X).ravel()[rows]
                                      for c in cols])
        means = grouper @ values
        output[:, start:start + len(cols)] = means.toarray() if sparse.issparse(means) else np.asarray(means)
    # Restore the stored precision after stable float64 accumulation. Otherwise
    # averaging repeated float32 predictions can manufacture rank differences
    # of ~1e-15 between biologically identical groups.
    dtype = getattr(getattr(adata, 'X', None), 'dtype', np.dtype('float32'))
    if np.dtype(dtype).kind == 'f' and np.dtype(dtype).itemsize <= 4:
        output = output.astype(dtype).astype(np.float64)
    return output


def compute(rna, other, modality, species='human', state=''):
    rnames = pd.Index(rna['obs_names']).astype(str)
    onames = pd.Index(other['obs_names']).astype(str)
    if not rnames.is_unique or not onames.is_unique:
        raise ValueError('Cell identifiers must be unique to pair expression measurements.')
    other_rows = onames.get_indexer(rnames)
    keep = other_rows >= 0
    populations = np.asarray(rna['populations']).astype(str)
    if state:
        keep &= populations == state
        field = next((c for c in ('donor', 'Donor', 'donor_id', 'sample_id', 'Library', 'sample')
                      if c in rna['adata'].obs), '')
        if not field:
            raise ValueError('Within-state correlations require a donor or sample annotation.')
        labels = rna['adata'].obs[field].astype(object).fillna('').astype(str).to_numpy()
        unit = f'{field} averages within {state}'
    else:
        labels, unit = populations, 'cell-state averages'
    keep &= ~np.isin(labels, ['', 'nan', 'None'])
    rows = np.flatnonzero(keep)
    codes, levels = pd.factorize(labels[rows], sort=True)
    counts = np.bincount(codes, minlength=len(levels))
    retained = counts >= 5
    grouper = sparse.csr_matrix((1.0 / counts[codes], (codes, np.arange(len(rows)))),
                               shape=(len(levels), len(rows)))[retained]
    labels = np.asarray(levels)[retained].tolist()
    counts = counts[retained].tolist()
    if len(labels) < 3:
        raise ValueError('At least three groups with five matched observations each are required.')
    pairs, skipped = paired_features(rna, other, modality, species)
    if not pairs:
        raise ValueError('No unambiguous RNA partners were found for the available features.')
    rcols = sorted({p[2] for p in pairs})
    ocols = sorted({p[3] for p in pairs})
    rm = _group_means(rna, rows, rcols, grouper)
    om = _group_means(other, other_rows[rows], ocols, grouper)
    ri, oi = {c: i for i, c in enumerate(rcols)}, {c: i for i, c in enumerate(ocols)}
    results, constant = [], 0
    for feature, gene, i, j in pairs:
        x, y = rm[:, ri[i]], om[:, oi[j]]
        valid = np.isfinite(x) & np.isfinite(y)
        if valid.sum() < 3 or np.ptp(x[valid]) == 0 or np.ptp(y[valid]) == 0:
            constant += 1
            continue
        rho = float(spearmanr(x[valid], y[valid]).statistic)
        if not np.isfinite(rho):
            continue
        results.append(dict(feature=feature, gene=gene, rho=rho, n_units=int(valid.sum()),
                            x=x[valid].tolist(), y=y[valid].tolist(),
                            labels=np.asarray(labels)[valid].tolist(),
                            n_cells=np.asarray(counts)[valid].tolist()))
    return dict(pairs=results, skipped=skipped, constant=constant, unit=unit,
                matched_observations=int(len(rows)), n_units=len(labels))


def answer_if_requested(app, meta, question):
    modality = requested_modality(question)
    if modality is None:
        return None
    from .app import _get_expression_cache
    result = dict(question=question, intent='cross_modal_correlation', status='ok',
                  reading=dict(intent='cross_modal_correlation', modality=modality, source='expression'))
    try:
        rna = _get_expression_cache(app, meta, 'rna')
        other = _get_expression_cache(app, meta, modality)
        states = sorted(set(map(str, rna['populations'])))
        named = [s for s in sorted(states, key=len, reverse=True)
                 if re.search(r'(?<![\w])' + re.escape(s) + r'(?![\w])', question, re.I)]
        state = named[0] if named else ''
        species = meta.get('species') or 'human'
        # Cache only group means/correlations, not extra full expression matrices.
        key = (meta.get('job_id'), id(rna['adata']), id(other['adata']), rna.get('source_stamp'), other.get('source_stamp'), species, state)
        with _LOCK:
            if key not in _CACHE:
                _CACHE[key] = compute(rna, other, modality, species, state)
                while len(_CACHE) > 8:
                    _CACHE.popitem(last=False)
            data = _CACHE[key]
            _CACHE.move_to_end(key)
    except (ValueError, KeyError, FileNotFoundError) as exc:
        result.update(status='not_covered', answer=str(exc))
        return result
    pairs = data['pairs']
    named_pairs = [p for p in pairs if any(re.search(r'(?<![\w])' + re.escape(name) + r'(?![\w])', question, re.I)
                                          for name in (p['feature'], p['gene']))]
    if named_pairs:
        pairs = named_pairs
    discordant = bool(re.search('discordan|disagree', question, re.I))
    weak = bool(re.search(r'uncorrelat|not correlated|poorly correlat|weak\w* correlat|least correlat|near.zero', question, re.I))
    relation = 'near_zero' if weak else 'negative' if discordant else 'all'
    pairs = sorted(pairs, key=lambda p: abs(p['rho']) if weak else p['rho'] if discordant else -abs(p['rho']))
    rows = [dict(feature=p['feature'], gene=p['gene'], rho=p['rho'], abs_rho=abs(p['rho']), n_units=p['n_units']) for p in pairs]
    matched = sum(abs(p['rho']) <= .2 if weak else p['rho'] < 0 if discordant else True for p in pairs)
    selection = (f'{matched} near-zero pairs (|ρ| ≤ 0.2; adjustable), ranked closest to zero. '
                 'Near-zero rank correlation does not establish statistical independence. ' if weak else
                 f'{matched} negatively correlated pairs selected. ' if discordant else '')
    label = 'TF activity' if modality == 'grn_tf' else 'ADT abundance'
    result.update(answer=(f'{len(pairs)} available '
                          f'{label}–RNA pairs across {data["n_units"]} {data["unit"]}. ' + selection +
                          'Each point is a paired average of stored Expression values; no differential results are needed. '
                          'These are descriptive correlations, not tests across independent biological replicates. '
                          'Imputed measurements derive from RNA and are not independent validation. '
                          + (f'{len(data["skipped"])} features had no unambiguous RNA partner; ' if data['skipped'] else '')
                          + f'{data["constant"]} constant or insufficient pairs were omitted.'
                          + (' ADT subunits are separate partners; gene-level RNA cannot resolve CD45 protein isoforms.' if modality == 'adt' else '')),
                  table=dict(columns=['feature', 'gene', 'rho', 'abs_rho', 'n_units'], rows=rows),
                  result_controls=dict(relationship=relation, near_zero_threshold=.2,
                                       sort_by='rho' if discordant and not weak else 'abs_rho',
                                       sort_direction='asc' if weak or discordant else 'desc'),
                  provenance=dict(source='expression', modality=modality, unit=data['unit'],
                                  mapping='curated ADT–RNA map' if modality == 'adt' else 'gene symbol',
                                  unmapped_features=data['skipped'], matched_observations=data['matched_observations']),
                  plot=dict(kind='cross_modal', pairs=pairs, y_label=label, unit=data['unit']))
    return result
