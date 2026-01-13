#!/usr/bin/env python3
import argparse
import os
import sys
from collections import defaultdict, OrderedDict

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['figure.facecolor'] = 'white'

sys.path.insert(1, os.path.join(os.path.dirname(__file__), '..'))
from long_read import isoform_matrix as iso

GENE_ID_PREFIXES = ('ENS',)


def looks_like_gene_id(value):
    if not value:
        return False
    if not value.startswith(GENE_ID_PREFIXES):
        return False
    # Basic heuristic for Ensembl-like IDs.
    return any(char.isdigit() for char in value)


def parse_feature_identifier(feature, gene_hint=None):
    if ':' in feature:
        left, right = feature.split(':', 1)
        if looks_like_gene_id(left):
            return left, right
        if looks_like_gene_id(right):
            return right, left
        return (gene_hint if looks_like_gene_id(gene_hint) else None), left
    return (gene_hint if looks_like_gene_id(gene_hint) else None), feature


def read_gene_model(path):
    gene_segments = defaultdict(list)
    exon_lookup = {}

    with open(path, 'r') as handle:
        first = handle.readline()
        if not first:
            return gene_segments, exon_lookup
        header = first.rstrip('\n').split('\t')
        has_header = header and header[0] == 'gene'
        if has_header:
            header_map = {name: idx for idx, name in enumerate(header)}
            idx_gene = header_map.get('gene', 0)
            idx_exon = header_map.get('exon-id', 1)
            idx_chr = header_map.get('chromosome', 2)
            idx_strand = header_map.get('strand', 3)
            idx_start = header_map.get('exon-region-start(s)', 4)
            idx_end = header_map.get('exon-region-stop(s)', 5)
        else:
            idx_gene, idx_exon, idx_chr, idx_strand, idx_start, idx_end = 0, 1, 2, 3, 4, 5
            parts = header
            _ingest_gene_model_row(parts, gene_segments, exon_lookup,
                                   idx_gene, idx_exon, idx_chr, idx_strand, idx_start, idx_end)

        for line in handle:
            parts = line.rstrip('\n').split('\t')
            _ingest_gene_model_row(parts, gene_segments, exon_lookup,
                                   idx_gene, idx_exon, idx_chr, idx_strand, idx_start, idx_end)

    return gene_segments, exon_lookup


def _ingest_gene_model_row(parts, gene_segments, exon_lookup,
                           idx_gene, idx_exon, idx_chr, idx_strand, idx_start, idx_end):
    if len(parts) <= max(idx_gene, idx_exon, idx_chr, idx_strand, idx_start, idx_end):
        return
    gene = parts[idx_gene].strip()
    exon_id = parts[idx_exon].strip()
    chrom = parts[idx_chr].strip()
    strand = parts[idx_strand].strip()
    try:
        start = int(float(parts[idx_start]))
        end = int(float(parts[idx_end]))
    except ValueError:
        return
    if not gene or not exon_id:
        return
    start, end = (start, end) if start <= end else (end, start)
    feature_type = exon_id[0].upper()
    segment = {
        'gene': gene,
        'exon_id': exon_id,
        'chrom': chrom,
        'strand': strand,
        'start': start,
        'end': end,
        'type': feature_type
    }
    gene_segments[gene].append(segment)
    exon_lookup[(gene, exon_id)] = segment


def build_gene_maps(gene_segments, intron_scale):
    gene_maps = {}
    for gene, segments in gene_segments.items():
        ordered = sorted(segments, key=lambda s: s['start'])
        display_pos = 0.0
        mapped_segments = []
        for seg in ordered:
            length = max(1, abs(seg['end'] - seg['start']) + 1)
            scale = intron_scale if seg['type'] == 'I' else 1.0
            display_start = display_pos
            display_end = display_pos + (length * scale)
            mapped_segments.append({
                **seg,
                'display_start': display_start,
                'display_end': display_end,
                'scale': scale
            })
            display_pos = display_end
        gene_maps[gene] = {
            'segments': mapped_segments,
            'total_length': display_pos
        }
    return gene_maps


def map_coord(gene_map, coord):
    segments = gene_map['segments']
    if not segments:
        return coord
    for seg in segments:
        if seg['start'] <= coord <= seg['end']:
            return seg['display_start'] + (coord - seg['start']) * seg['scale']
    if coord < segments[0]['start']:
        return segments[0]['display_start']
    return segments[-1]['display_end']


def _iter_alias_tokens(value, separators=('|', ':')):
    pending = [value]
    seen = set()
    while pending:
        token = pending.pop()
        if token in seen:
            continue
        seen.add(token)
        if not token:
            continue
        for sep in separators:
            if sep in token:
                pending.extend([part for part in token.split(sep) if part])
    return [t for t in seen if t]


def load_transcript_associations(path, target_gene):
    structures = {}
    gene_strand = {}
    with open(path, 'r') as handle:
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            gene, strand, exon_sequence, transcript_id = parts[:4]
            if gene != target_gene:
                continue
            tokens = [t for t in exon_sequence.split('|') if t]
            if not tokens:
                continue
            structures[transcript_id] = tokens
            for alias in _iter_alias_tokens(transcript_id):
                if alias not in structures:
                    structures[alias] = tokens
            gene_strand[gene] = strand
    return structures, gene_strand


def sample_transcript_ids_from_associations(path, target_gene, limit=5):
    samples = []
    with open(path, 'r') as handle:
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            gene, _, _, transcript_id = parts[:4]
            if gene != target_gene:
                continue
            samples.append(transcript_id)
            if len(samples) >= limit:
                break
    return samples


def load_transcript_association_index(path, target_gene):
    transcript_ids = set()
    alias_to_transcript = {}
    with open(path, 'r') as handle:
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            gene, _, _, transcript_id = parts[:4]
            if gene != target_gene:
                continue
            transcript_ids.add(transcript_id)
            for alias in _iter_alias_tokens(transcript_id):
                if alias not in alias_to_transcript:
                    alias_to_transcript[alias] = transcript_id
    return transcript_ids, alias_to_transcript


def split_token(token, default_gene):
    gene_id = default_gene
    core = token
    if ':' in token:
        gene_id, core = token.split(':', 1)
    coord = None
    base = core
    if '_' in core:
        base, coord_str = core.rsplit('_', 1)
        try:
            coord = int(coord_str)
        except ValueError:
            base = core
            coord = None
    return gene_id, base, coord


def build_isoform_segments(tokens, exon_lookup, default_gene):
    parsed = []
    for token in tokens:
        gene_id, base, coord = split_token(token, default_gene)
        parsed.append({'gene': gene_id, 'base': base, 'coord': coord})

    grouped = []
    for item in parsed:
        if grouped and grouped[-1]['gene'] == item['gene'] and grouped[-1]['base'] == item['base']:
            if item['coord'] is not None:
                grouped[-1]['coords'].append(item['coord'])
        else:
            grouped.append({
                'gene': item['gene'],
                'base': item['base'],
                'coords': [item['coord']] if item['coord'] is not None else []
            })

    direction = None
    for i in range(len(grouped) - 1):
        seg1 = exon_lookup.get((grouped[i]['gene'], grouped[i]['base']))
        seg2 = exon_lookup.get((grouped[i + 1]['gene'], grouped[i + 1]['base']))
        if seg1 and seg2:
            if seg2['start'] > seg1['start']:
                direction = 1
            elif seg2['start'] < seg1['start']:
                direction = -1
            break
    if direction is None:
        direction = 1

    segments = []
    for idx, group in enumerate(grouped):
        seg = exon_lookup.get((group['gene'], group['base']))
        if not seg:
            continue
        start, end = seg['start'], seg['end']
        if start > end:
            start, end = end, start
        coords = [c for c in group['coords'] if c is not None]
        if coords:
            coord_min = min(coords)
            coord_max = max(coords)
            if idx == 0:
                if direction == 1:
                    start = max(start, coord_min)
                else:
                    end = min(end, coord_max)
            if idx == len(grouped) - 1:
                if direction == 1:
                    end = min(end, coord_max)
                else:
                    start = max(start, coord_min)
            if idx != 0 and idx != len(grouped) - 1 and len(coords) > 1:
                start = max(start, coord_min)
                end = min(end, coord_max)
        if end < start:
            start, end = end, start
        label = None
        if seg['type'] == 'E':
            label = normalize_token(group['base'], default_gene, 'block')
        segments.append({
            'gene': group['gene'],
            'start': start,
            'end': end,
            'type': seg['type'],
            'label': label
        })
    return segments


def normalize_token(token, default_gene, mode):
    gene_id, base, _ = split_token(token, default_gene)
    core = base
    if mode == 'full':
        return f"{gene_id}:{base}" if gene_id != default_gene else base
    if mode == 'base':
        return core
    if mode == 'block':
        if len(core) > 1 and '.' in core:
            prefix = core[0]
            number = core[1:].split('.', 1)[0]
            return f"{prefix}{number}"
        return core
    return core


def resolve_isoform_id(raw_id, transcript_structures):
    if raw_id in transcript_structures:
        return raw_id
    for candidate in _iter_alias_tokens(raw_id):
        if candidate in transcript_structures:
            return candidate
    if '_' in raw_id:
        tail = raw_id.split('_')[-1]
        if tail in transcript_structures:
            return tail
        head = raw_id.split('_')[0]
        if head in transcript_structures:
            return head
    return None


def resolve_isoform_id_fuzzy(raw_id, transcript_structures):
    resolved = resolve_isoform_id(raw_id, transcript_structures)
    if resolved:
        return resolved, []
    suffix_matches = [key for key in transcript_structures if key.endswith(raw_id)]
    if len(suffix_matches) == 1:
        return suffix_matches[0], []
    if not suffix_matches:
        contains_matches = [key for key in transcript_structures if raw_id in key]
        if len(contains_matches) == 1:
            return contains_matches[0], []
        return None, contains_matches
    return None, suffix_matches


def _load_adata(input_path):
    if input_path.endswith('.h5ad'):
        return ad.read_h5ad(input_path)
    if input_path.endswith('.mtx') or input_path.endswith('.mtx.gz') or os.path.isdir(input_path):
        mtx_dir = input_path if os.path.isdir(input_path) else os.path.dirname(input_path)
        feature = _find_mtx_file(mtx_dir, ['genes.tsv', 'features.tsv', 'genes.tsv.gz', 'features.tsv.gz'])
        barcode = _find_mtx_file(mtx_dir, ['barcodes.tsv', 'barcodes.tsv.gz'])
        matrix = _find_mtx_file(mtx_dir, ['matrix.mtx', 'matrix.mtx.gz'])
        if not feature or not barcode or not matrix:
            raise ValueError(f"Missing 10x mtx files in {mtx_dir}.")
        return iso.mtx_to_adata(
            int_folder=mtx_dir,
            gene_is_index=True,
            feature=os.path.basename(feature),
            feature_col=0,
            barcode=os.path.basename(barcode),
            barcode_col=0,
            matrix=os.path.basename(matrix),
            rev=False
        )
    raise ValueError(f"Unsupported input format: {input_path}")


def _find_mtx_file(folder, candidates):
    for name in candidates:
        path = os.path.join(folder, name)
        if os.path.exists(path):
            return path
    return None


def load_isoform_counts(input_path, target_gene, groupby=None, group_values=None):
    adata = _load_adata(input_path)
    if groupby and groupby in adata.obs:
        if group_values:
            mask = adata.obs[groupby].isin(group_values)
            adata = adata[mask, :]
    counts = np.asarray(adata.X.sum(axis=0)).ravel()
    var_names = np.asarray(adata.var_names)
    gene_col = None
    if 'gene' in adata.var.columns:
        gene_col = np.asarray(adata.var['gene'])

    records = []
    for idx, var in enumerate(var_names):
        gene_hint = gene_col[idx] if gene_col is not None else None
        gene_id, isoform_id = parse_feature_identifier(str(var), gene_hint=gene_hint)
        if gene_id is not None and gene_id != target_gene:
            continue
        records.append({
            'var': var,
            'gene_id': gene_id,
            'isoform_id': isoform_id,
            'count': counts[idx]
        })
    return records


def filter_records_by_structures(records, transcript_structures):
    filtered = []
    for rec in records:
        if resolve_isoform_id(rec['isoform_id'], transcript_structures):
            filtered.append(rec)
    return filtered


def coerce_read_count(value):
    try:
        count = float(value)
    except (TypeError, ValueError):
        return 0
    if count <= 0:
        return 0
    rounded = int(round(count))
    if abs(count - rounded) < 1e-6:
        return rounded
    return max(1, rounded)


def build_feature_set(tokens, feature_mode):
    if feature_mode == 'junctions':
        return {f"{tokens[i]}->{tokens[i + 1]}" for i in range(len(tokens) - 1)}
    if feature_mode == 'both':
        features = {f"{tokens[i]}->{tokens[i + 1]}" for i in range(len(tokens) - 1)}
        features.update(tokens)
        return features
    return set(tokens)


def strip_terminal_coords(tokens, default_gene):
    if not tokens:
        return []
    trimmed = list(tokens)
    drop_indices = set()
    for idx in (0, len(trimmed) - 1):
        gene_id, base, coord = split_token(trimmed[idx], default_gene)
        if coord is not None and base.startswith('E'):
            drop_indices.add(idx)
    if drop_indices:
        return [tok for i, tok in enumerate(trimmed) if i not in drop_indices]
    return trimmed


def filter_exon_tokens(tokens, default_gene):
    exon_tokens = []
    for tok in tokens:
        gene_id, base, _ = split_token(tok, default_gene)
        if base.startswith('E'):
            if '_' in tok:
                continue
            if gene_id != default_gene and ':' not in tok:
                tok = f"{gene_id}:{base}"
            exon_tokens.append(tok)
    return exon_tokens


def filter_exon_intron_tokens(tokens, default_gene):
    kept_tokens = []
    for tok in tokens:
        gene_id, base, _ = split_token(tok, default_gene)
        if base.startswith(('E', 'I')):
            if base.startswith('E') and '_' in tok:
                continue
            if gene_id != default_gene and ':' not in tok:
                tok = f"{gene_id}:{base}"
            kept_tokens.append(tok)
    return kept_tokens


def cluster_by_feature(items, min_split_fraction=0.05):
    def split(group):
        if len(group) <= 1:
            return [group]
        total_weight = sum(item['weight'] for item in group)
        if total_weight <= 0:
            return [group]

        feature_weights = defaultdict(float)
        for item in group:
            for feat in item['feature_set']:
                feature_weights[feat] += item['weight']

        best_feat = None
        best_score = 0
        for feat, weight in feature_weights.items():
            if weight <= 0 or weight >= total_weight:
                continue
            score = min(weight, total_weight - weight)
            if score > best_score:
                best_score = score
                best_feat = feat

        if best_feat is None:
            return [group]

        split_fraction = best_score / total_weight
        if split_fraction < min_split_fraction:
            return [group]

        present = [item for item in group if best_feat in item['feature_set']]
        absent = [item for item in group if best_feat not in item['feature_set']]
        if not present or not absent:
            return [group]

        weight_present = sum(item['weight'] for item in present)
        weight_absent = sum(item['weight'] for item in absent)
        if weight_present >= weight_absent:
            ordered = [present, absent]
        else:
            ordered = [absent, present]
        return split(ordered[0]) + split(ordered[1])

    return split(items)


def cluster_by_jaccard(items, threshold=0.4):
    clusters = []

    def jaccard(a_set, b_set):
        union = a_set | b_set
        inter = a_set & b_set
        return (len(inter) / len(union)) if union else 0.0

    for item in sorted(items, key=lambda x: x['weight'], reverse=True):
        best_idx = None
        best_score = -1.0
        for idx, cluster in enumerate(clusters):
            score = jaccard(item['cluster_set'], cluster['rep_set'])
            if score > best_score:
                best_score = score
                best_idx = idx
        if best_idx is not None and best_score >= threshold:
            clusters[best_idx]['items'].append(item)
        else:
            clusters.append({'items': [item], 'rep_set': item['cluster_set']})
    return [cluster['items'] for cluster in clusters]


def longest_common_substring_length(seq_a, seq_b):
    if not seq_a or not seq_b:
        return 0
    if len(seq_a) > len(seq_b):
        seq_a, seq_b = seq_b, seq_a
    prev = [0] * (len(seq_a) + 1)
    best = 0
    for token in seq_b:
        current = [0]
        for idx, a_token in enumerate(seq_a, start=1):
            if token == a_token:
                val = prev[idx - 1] + 1
            else:
                val = 0
            current.append(val)
            if val > best:
                best = val
        prev = current
    return best


def substring_similarity(seq_a, seq_b):
    if not seq_a or not seq_b:
        return 0.0
    short = seq_a if len(seq_a) <= len(seq_b) else seq_b
    long = seq_b if short is seq_a else seq_a
    lcs = longest_common_substring_length(short, long)
    return lcs / float(len(short))


def cluster_by_substring(items, threshold=0.4):
    clusters = []

    for item in sorted(items, key=lambda x: x['weight'], reverse=True):
        best_idx = None
        best_score = -1.0
        for idx, cluster in enumerate(clusters):
            min_score = 1.0
            for other in cluster['items']:
                score = substring_similarity(item['cluster_seq'], other['cluster_seq'])
                if score < min_score:
                    min_score = score
                if min_score < threshold:
                    break
            if min_score >= threshold and min_score > best_score:
                best_score = min_score
                best_idx = idx
        if best_idx is not None:
            clusters[best_idx]['items'].append(item)
        else:
            clusters.append({'items': [item]})
    return [cluster['items'] for cluster in clusters]


def is_subsequence(short_seq, long_seq):
    if not short_seq:
        return True
    long_iter = iter(long_seq)
    for item in short_seq:
        for candidate in long_iter:
            if candidate == item:
                break
        else:
            return False
    return True


def build_plotted_isoforms(isoforms, transcript_structures, target_gene, cluster_mode,
                           cluster_features, min_count, max_isoforms):
    filtered = [i for i in isoforms if i['count'] >= min_count]
    missing_ids = []
    resolved_items = []
    for iso in filtered:
        resolved = resolve_isoform_id(iso['isoform_id'], transcript_structures)
        if not resolved:
            missing_ids.append(iso['isoform_id'])
            continue
        resolved_items.append((iso, resolved))

    resolved_items.sort(key=lambda x: x[0]['count'], reverse=True)
    resolved_items = resolved_items[:max_isoforms]

    plotted = []
    for iso, resolved in resolved_items:
        tokens = transcript_structures[resolved]
        tokens_trimmed = strip_terminal_coords(tokens, target_gene)
        cluster_tokens = filter_exon_tokens(tokens_trimmed, target_gene)
        cluster_set = set(cluster_tokens)
        tokens_norm = [normalize_token(tok, target_gene, cluster_mode) for tok in tokens]
        plotted.append({
            **iso,
            'resolved_id': resolved,
            'tokens': tokens,
            'tokens_trimmed': tokens_trimmed,
            'cluster_tokens': cluster_tokens,
            'cluster_set': cluster_set,
            'tokens_norm': tokens_norm,
            'feature_set': build_feature_set(tokens_trimmed, cluster_features)
        })
    return plotted, len(missing_ids), missing_ids


def build_plotted_rows(plotted, cluster_label_map=None, cluster_index_map=None):
    rows = []
    seen = set()
    for iso in plotted:
        isoform_id = iso.get('isoform_id') or iso.get('resolved_id')
        resolved_id = iso.get('resolved_id') or isoform_id
        key = (isoform_id, resolved_id)
        if key in seen:
            continue
        seen.add(key)
        row = {
            'isoform_id': isoform_id,
            'resolved_id': resolved_id,
            'count': iso.get('count', 0),
            'cluster_label': (cluster_label_map or {}).get(resolved_id)
                            or (cluster_label_map or {}).get(isoform_id)
        }
        if cluster_index_map:
            row['cluster_index'] = cluster_index_map.get(resolved_id)
        rows.append(row)
    return rows


def export_plotted_isoforms(plotted, output_path, cluster_label_map=None, cluster_index_map=None):
    if not plotted:
        return None, []
    rows = build_plotted_rows(plotted, cluster_label_map, cluster_index_map)
    output_base = os.path.splitext(output_path)[0]
    isoform_path = f"{output_base}_isoform_ids.tsv"
    pd.DataFrame(rows).to_csv(isoform_path, sep='\t', index=False)
    return isoform_path, rows


def sample_isoform_ids(records, limit=5):
    samples = []
    for rec in records:
        isoform_id = rec.get('isoform_id')
        if isoform_id:
            samples.append(isoform_id)
        if len(samples) >= limit:
            break
    return samples


def build_cluster_label_map(records, transcript_structures, target_gene, cluster_mode,
                            cluster_features, cluster_strategy,
                            cluster_similarity_threshold, min_split_fraction,
                            max_isoforms, min_count, cluster_isoforms,
                            inspect_isoforms=None):
    if cluster_strategy is None:
        cluster_strategy = 'substring' if cluster_isoforms else 'jaccard'
    if cluster_strategy in ('subsequence', 'overlap'):
        cluster_strategy = 'substring'
    if cluster_similarity_threshold is None:
        cluster_similarity_threshold = 0.85 if cluster_isoforms else 0.4
    plotted, _, _ = build_plotted_isoforms(
        records,
        transcript_structures,
        target_gene,
        cluster_mode,
        cluster_features,
        min_count,
        max_isoforms
    )
    plotted_ids = set()
    for iso in plotted:
        iso_key = iso.get('isoform_id', iso.get('resolved_id'))
        if iso_key:
            plotted_ids.add(iso_key)
        if iso.get('resolved_id'):
            plotted_ids.add(iso.get('resolved_id'))
    if not plotted:
        return {}, plotted_ids
    if inspect_isoforms:
        max_isoforms = max(max_isoforms, len(records))
        plotted, _, _ = build_plotted_isoforms(
            records,
            transcript_structures,
            target_gene,
            cluster_mode,
            cluster_features,
            min_count,
            max_isoforms
        )
    grouped_structures = group_structures(
        plotted,
        target_gene,
        cluster_features,
        cluster_strategy,
        cluster_similarity_threshold,
        min_split_fraction,
        include_introns=cluster_isoforms,
        report=False
    )
    if cluster_isoforms:
        for group in grouped_structures:
            group.sort(key=structure_max_length, reverse=True)
            for structure in group:
                structure['items'].sort(
                    key=lambda iso: (isoform_length(iso), iso.get('count', 0)),
                    reverse=True
                )
        grouped_structures.sort(
            key=lambda group: (group_max_length(group), sum(item['weight'] for item in group)),
            reverse=True
        )
    label_map = {}
    for group in grouped_structures:
        label = select_cluster_label(group) if cluster_isoforms else None
        for structure in group:
            for iso in structure['items']:
                iso_key = iso.get('isoform_id', iso.get('resolved_id'))
                if cluster_isoforms:
                    label_map[iso_key] = label
                    label_map[iso.get('resolved_id')] = label
                else:
                    label_map[iso_key] = iso.get('resolved_id')
    return label_map, plotted_ids


def report_isoform_structures(records, transcript_structures, target_gene, isoform_ids,
                              include_introns=False, cluster_label_map=None,
                              plotted_ids=None, plot_cluster_label_map=None,
                              plot_cluster_index_map=None):
    if not isoform_ids:
        return
    record_lookup = {}
    for rec in records:
        isoform_id = rec.get('isoform_id')
        var = rec.get('var')
        if isoform_id:
            record_lookup.setdefault(isoform_id, rec)
        if var:
            record_lookup.setdefault(var, rec)
    resolved_ids = []
    cluster_seqs = []
    print("Isoform structure inspection:")
    for raw_id in isoform_ids:
        resolved, matches = resolve_isoform_id_fuzzy(raw_id, transcript_structures)
        if not resolved:
            if matches:
                preview = ", ".join(matches[:3])
                print(f"  {raw_id}: multiple matches ({len(matches)}) {preview}")
            else:
                print(f"  {raw_id}: no match in transcript_associations")
            continue
        tokens = transcript_structures[resolved]
        trimmed = strip_terminal_coords(list(tokens), target_gene)
        if include_introns:
            cluster_seq = filter_exon_intron_tokens(list(trimmed), target_gene)
        else:
            cluster_seq = filter_exon_tokens(list(trimmed), target_gene)
        resolved_ids.append(resolved)
        cluster_seqs.append(cluster_seq)
        rec = record_lookup.get(raw_id) or record_lookup.get(resolved)
        label = resolved if raw_id == resolved else f"{raw_id} -> {resolved}"
        print(f"  {label}")
        if cluster_label_map:
            display_label = cluster_label_map.get(resolved) or cluster_label_map.get(raw_id)
            if display_label:
                print(f"    cluster_label={display_label}")
        if plot_cluster_label_map:
            plot_label = plot_cluster_label_map.get(resolved) or plot_cluster_label_map.get(raw_id)
            if plot_label:
                print(f"    plot_cluster_label={plot_label}")
        if plot_cluster_index_map:
            plot_index = plot_cluster_index_map.get(resolved) or plot_cluster_index_map.get(raw_id)
            if plot_index is not None:
                print(f"    plot_cluster_index={plot_index}")
        if plotted_ids is not None:
            plotted_flag = (resolved in plotted_ids) or (raw_id in plotted_ids)
            print(f"    plotted={str(plotted_flag).lower()}")
        if rec and rec.get('count') is not None:
            print(f"    count={rec.get('count')}")
        if rec and rec.get('var') and rec.get('var') not in (raw_id, resolved):
            print(f"    var={rec.get('var')}")
        print(f"    tokens={'|'.join(tokens)}")
        print(f"    trimmed={'|'.join(trimmed)}")
        print(f"    cluster_seq={'|'.join(cluster_seq)}")
    if len(resolved_ids) == 2:
        similarity = substring_similarity(cluster_seqs[0], cluster_seqs[1])
        print(f"  substring_similarity={similarity:.3f}")


def group_structures(plotted, target_gene, cluster_features, cluster_strategy,
                     cluster_similarity_threshold, min_split_fraction,
                     include_introns=False, report=False):
    structure_groups = OrderedDict()
    for iso in plotted:
        key = tuple(iso['tokens'])
        trimmed_key = tuple(strip_terminal_coords(list(key), target_gene))
        if include_introns:
            cluster_tokens = filter_exon_intron_tokens(list(trimmed_key), target_gene)
        else:
            cluster_tokens = filter_exon_tokens(list(trimmed_key), target_gene)
        cluster_set = set(cluster_tokens)
        if key not in structure_groups:
            structure_groups[key] = {
                'structure_key': key,
                'structure_key_trimmed': trimmed_key,
                'items': [],
                'weight': 0,
                'feature_set': build_feature_set(list(trimmed_key), cluster_features),
                'cluster_set': cluster_set,
                'cluster_seq': cluster_tokens
            }
        structure_groups[key]['items'].append(iso)
        structure_groups[key]['weight'] += iso['count']

    structure_items = list(structure_groups.values())
    if cluster_strategy == 'jaccard':
        grouped_structures = cluster_by_jaccard(
            structure_items,
            threshold=cluster_similarity_threshold
        )
    elif cluster_strategy in ('substring', 'subsequence', 'overlap'):
        grouped_structures = cluster_by_substring(
            structure_items,
            threshold=cluster_similarity_threshold
        )
    elif cluster_strategy == 'feature':
        grouped_structures = cluster_by_feature(
            structure_items,
            min_split_fraction=min_split_fraction
        )
    elif cluster_strategy == 'structure':
        grouped_structures = [[item] for item in structure_items]
    else:
        grouped_structures = [[item] for item in sorted(structure_items, key=lambda x: x['weight'], reverse=True)]

    if report:
        total_structures = len(structure_items)
        cluster_count = len(grouped_structures)
        merged_count = total_structures - cluster_count
        cluster_sizes = [len(group) for group in grouped_structures]
        mean_size = (sum(cluster_sizes) / len(cluster_sizes)) if cluster_sizes else 0
        max_size = max(cluster_sizes) if cluster_sizes else 0
        summary_parts = [
            "Cluster summary",
            f"strategy={cluster_strategy}",
            f"structures={total_structures}",
            f"clusters={cluster_count}",
            f"merged={merged_count}",
            f"mean_size={mean_size:.2f}",
            f"max_size={max_size}"
        ]
        if cluster_strategy in ('substring', 'subsequence', 'overlap', 'jaccard'):
            summary_parts.append(f"threshold={cluster_similarity_threshold}")
        else:
            summary_parts.append(f"min_split_fraction={min_split_fraction}")
        print(" ".join(summary_parts))

    for group in grouped_structures:
        rep = max(group, key=lambda x: x['weight'])
        rep_set = rep['cluster_set']

        def similarity_key(item):
            subset = item['cluster_set'].issubset(rep_set)
            union = rep_set | item['cluster_set']
            inter = rep_set & item['cluster_set']
            jaccard = len(inter) / len(union) if union else 0.0
            overlap = (len(inter) / min(len(rep_set), len(item['cluster_set']))) if rep_set and item['cluster_set'] else 0.0
            rep_rank = 0 if item is rep else 1
            return (rep_rank, -int(subset), -overlap, -jaccard, -item['weight'])

        group.sort(key=similarity_key)
        for entry in group:
            entry['items'].sort(key=lambda x: x['count'], reverse=True)

    return grouped_structures


def assign_isoform_colors(grouped_structures, palette):
    isoform_colors = {}
    color_idx = 0
    for group in grouped_structures:
        for structure in group:
            for iso in structure['items']:
                iso_key = iso.get('isoform_id', iso['resolved_id'])
                if iso_key not in isoform_colors:
                    isoform_colors[iso_key] = palette[color_idx % len(palette)]
                    color_idx += 1
    return isoform_colors


def assign_cluster_colors(grouped_structures, palette):
    isoform_colors = {}
    for group_idx, group in enumerate(grouped_structures):
        color = palette[group_idx % len(palette)]
        for structure in group:
            for iso in structure['items']:
                iso_key = iso.get('isoform_id', iso['resolved_id'])
                if iso_key not in isoform_colors:
                    isoform_colors[iso_key] = color
    return isoform_colors


def select_cluster_label(group):
    best_label = None
    best_score = None
    for structure in group:
        for iso in structure['items']:
            label = iso.get('resolved_id') or iso.get('isoform_id')
            if not label:
                continue
            tokens = iso.get('tokens_trimmed') or iso.get('tokens') or []
            score = (len(tokens), iso.get('count', 0), len(str(label)))
            if best_score is None or score > best_score:
                best_score = score
                best_label = label
    return best_label


def isoform_length(iso):
    tokens = iso.get('tokens_trimmed') or iso.get('tokens') or []
    return len(tokens)


def structure_max_length(structure):
    return max((isoform_length(iso) for iso in structure['items']), default=0)


def group_max_length(group):
    return max((structure_max_length(structure) for structure in group), default=0)


# Parameters:
# - isoforms: list of dicts with isoform IDs and read counts.
# - transcript_structures: transcript_id -> list of exon/intron tokens.
# - exon_lookup: (gene, exon_id) -> coordinate metadata.
# - gene_maps: gene -> display coordinate map (scaled introns).
# - cluster_mode: token normalization mode (full/base/block).
# - target_gene: gene ID used for filtering and token defaults.
# - intron_scale: display scale factor for intron lengths.
# - output_path: file path for the saved plot.
# - max_isoforms: cap on number of isoforms to render.
# - min_count: minimum read count to keep an isoform.
# - row_height: height of each read row in plot units.
# - row_gap: vertical gap between read rows.
# - group_gap: vertical gap between clustered isoform groups.
# - isoform_gap: vertical gap between isoform stacks (defaults to group_gap).
# - gene_gap: horizontal gap between genes in trans-splicing plots.
# - cluster_strategy: ordering method for isoform structures.
# - cluster_features: features used for clustering (tokens/junctions/both).
# - label_mode: label display mode (first/none/all).
# - min_split_fraction: minimum weighted split for feature clustering.
# - cluster_similarity_threshold: minimum similarity threshold for clustering.
# - grouped_structures: precomputed cluster order to reuse across files.
# - isoform_colors: precomputed isoform-to-color mapping.
# - isoform_counts: per-file isoform counts for consistent ordering.
# - cluster_isoforms: color isoforms by cluster and label with a representative isoform ID.
# - debug_transcripts: transcript IDs to print segment details for.
def plot_isoform_structures(isoforms, transcript_structures, exon_lookup, gene_maps,
                            cluster_mode, target_gene, intron_scale, output_path,
                            max_isoforms=50, min_count=1, row_height=0.0125,
                            row_gap=0.0, group_gap=0.0, isoform_gap=None, gene_gap=1000,
                            cluster_strategy='jaccard', cluster_features='tokens',
                            label_mode='first', min_split_fraction=0.05,
                            cluster_similarity_threshold=0.4,
                            debug_isoforms=None,
                            debug_transcripts=None,
                            grouped_structures=None,
                            isoform_colors=None,
                            isoform_counts=None,
                            transcript_associations_path=None,
                            cluster_isoforms=False,
                            return_cluster_labels=False):
    """Render stacked isoform read tracks for a single gene.

    isoforms: list of dicts with isoform IDs and read counts.
    transcript_structures: transcript_id -> list of exon/intron tokens.
    exon_lookup: (gene, exon_id) -> coordinate metadata.
    gene_maps: gene -> display coordinate map (scaled introns).
    cluster_mode: token normalization mode (full/base/block).
    target_gene: gene ID used for filtering and token defaults.
    intron_scale: display scale factor for intron lengths.
    output_path: file path for the saved plot.
    max_isoforms: cap on number of isoforms to render.
    min_count: minimum read count to keep an isoform.
    row_height: height of each read row in plot units.
    row_gap: vertical gap between read rows.
    group_gap: vertical gap between clustered isoform groups.
    isoform_gap: vertical gap between isoform stacks (defaults to group_gap).
    gene_gap: horizontal gap between genes in trans-splicing plots.
    cluster_strategy: ordering method for isoform structures.
    cluster_features: features used for clustering (tokens/junctions/both).
    label_mode: label display mode (first/none/all).
    min_split_fraction: minimum weighted split for feature clustering.
    cluster_similarity_threshold: minimum similarity threshold for clustering.
    grouped_structures: precomputed cluster order to reuse across files.
    isoform_colors: precomputed isoform-to-color mapping.
    isoform_counts: per-file isoform counts for consistent ordering.
    debug_transcripts: transcript IDs to print segment details for.
    """
    if cluster_strategy is None:
        cluster_strategy = 'substring' if cluster_isoforms else 'jaccard'
    if cluster_strategy in ('subsequence', 'overlap'):
        cluster_strategy = 'substring'
    if cluster_similarity_threshold is None:
        cluster_similarity_threshold = 0.85 if cluster_isoforms else 0.4
    if grouped_structures is None:
        plotted, missing_structures, missing_ids = build_plotted_isoforms(
            isoforms,
            transcript_structures,
            target_gene,
            cluster_mode,
            cluster_features,
            min_count,
            max_isoforms
        )
        if not plotted:
            print("No isoforms with transcript structures available to plot.")
            missing_sample = sorted(set(missing_ids))[:5]
            if missing_sample:
                print("Sample isoform IDs missing from transcript_associations:", missing_sample)
            if transcript_associations_path:
                transcript_ids, alias_to_transcript = load_transcript_association_index(
                    transcript_associations_path, target_gene
                )
                matrix_ids = {rec['isoform_id'] for rec in isoforms if rec.get('isoform_id')}
                matched_map = {}
                for matrix_id in matrix_ids:
                    resolved = resolve_isoform_id(matrix_id, transcript_structures)
                    if not resolved:
                        continue
                    matched_map[matrix_id] = alias_to_transcript.get(resolved, resolved)
                matched_matrix_ids = set(matched_map.keys())
                matched_transcript_ids = set(matched_map.values())
                matrix_only = sorted(matrix_ids - matched_matrix_ids)
                transcript_only = sorted(transcript_ids - matched_transcript_ids)
                overlap_examples = sorted(matched_matrix_ids)[:5]
                print(f"QC matrix isoforms: {len(matrix_ids)}")
                print(f"QC transcript_associations isoforms: {len(transcript_ids)}")
                print(f"QC overlap (matrix->transcript): {len(matched_matrix_ids)}")
                print(f"QC matrix-only: {len(matrix_only)}")
                print(f"QC transcript-only: {len(transcript_only)}")
                if overlap_examples:
                    print("QC overlap examples:", overlap_examples)
                if matrix_only:
                    print("QC matrix-only examples:", matrix_only[:5])
                if transcript_only:
                    print("QC transcript-only examples:", transcript_only[:5])
            raise ValueError("No isoforms with transcript structures available to plot.")
        grouped_structures = group_structures(
            plotted,
            target_gene,
            cluster_features,
            cluster_strategy,
            cluster_similarity_threshold,
            min_split_fraction,
            include_introns=cluster_isoforms,
            report=cluster_isoforms
        )
    else:
        missing_structures = 0
        plotted = []
        seen = set()
        for group in grouped_structures:
            for structure in group:
                for iso in structure['items']:
                    iso_key = iso.get('isoform_id', iso.get('resolved_id'))
                    if iso_key in seen:
                        continue
                    seen.add(iso_key)
                    plotted.append(iso)

    if cluster_isoforms:
        for group in grouped_structures:
            group.sort(key=structure_max_length, reverse=True)
            for structure in group:
                structure['items'].sort(
                    key=lambda iso: (isoform_length(iso), iso.get('count', 0)),
                    reverse=True
                )
        grouped_structures.sort(
            key=lambda group: (group_max_length(group), sum(item['weight'] for item in group)),
            reverse=True
        )

    group_labels = {}
    if cluster_isoforms:
        for group_idx, group in enumerate(grouped_structures):
            group_labels[group_idx] = select_cluster_label(group)

    if debug_isoforms:
        debug_set = set(debug_isoforms)
        debug_hits = []
        debug_lookup = {}
        for structure_idx, group in enumerate(grouped_structures):
            for item in group:
                for iso in item['items']:
                    debug_lookup[iso['resolved_id']] = iso
                    if iso.get('var'):
                        debug_lookup[iso['var']] = iso
                    if iso['resolved_id'] in debug_set or iso.get('var') in debug_set:
                        debug_hits.append({
                            'iso_id': iso['resolved_id'],
                            'var': iso.get('var'),
                            'structure_key': item['structure_key'],
                            'group_index': structure_idx,
                            'weight': item['weight']
                        })
        if debug_hits:
            print("Debug isoform clustering:")
            for hit in debug_hits:
                print(
                    f"  {hit['iso_id']} (var={hit['var']}) group={hit['group_index']} "
                    f"weight={hit['weight']} key={'|'.join(hit['structure_key'])}"
                )
            if len(debug_isoforms) >= 2:
                iso_a = debug_lookup.get(debug_isoforms[0])
                iso_b = debug_lookup.get(debug_isoforms[1])
                if iso_a and iso_b:
                    def exon_only(tokens):
                        return [t for t in tokens if t.startswith('E')]

                    def junctions(tokens):
                        return [f"{tokens[i]}->{tokens[i + 1]}" for i in range(len(tokens) - 1)]

                    raw_a = iso_a['tokens']
                    raw_b = iso_b['tokens']
                    trim_a = strip_terminal_coords(raw_a, target_gene)
                    trim_b = strip_terminal_coords(raw_b, target_gene)
                    cluster_a = filter_exon_tokens(trim_a, target_gene)
                    cluster_b = filter_exon_tokens(trim_b, target_gene)
                    norm_a = iso_a['tokens_norm']
                    norm_b = iso_b['tokens_norm']
                    exon_a = exon_only(trim_a)
                    exon_b = exon_only(trim_b)
                    junc_a = junctions(trim_a)
                    junc_b = junctions(trim_b)

                    def jaccard(a_set, b_set):
                        union = a_set | b_set
                        inter = a_set & b_set
                        return (len(inter) / len(union)) if union else 0.0
                    def overlap(a_set, b_set):
                        if not a_set or not b_set:
                            return 0.0
                        inter = a_set & b_set
                        return len(inter) / min(len(a_set), len(b_set))

                    print("  Debug comparison:")
                    print(f"    A={iso_a['resolved_id']} raw={'|'.join(raw_a)}")
                    print(f"    B={iso_b['resolved_id']} raw={'|'.join(raw_b)}")
                    print(f"    A_trim={'|'.join(trim_a)}")
                    print(f"    B_trim={'|'.join(trim_b)}")
                    print(f"    A_cluster={'|'.join(cluster_a)}")
                    print(f"    B_cluster={'|'.join(cluster_b)}")
                    print(f"    A_exons={'|'.join(exon_a)}")
                    print(f"    B_exons={'|'.join(exon_b)}")
                    print(f"    A_juncs={'|'.join(junc_a)}")
                    print(f"    B_juncs={'|'.join(junc_b)}")
                    print(
                        "    Similarity:",
                        f"cluster_jaccard={jaccard(set(cluster_a), set(cluster_b)):.3f}",
                        f"cluster_overlap={overlap(set(cluster_a), set(cluster_b)):.3f}",
                        f"exon_jaccard={jaccard(set(exon_a), set(exon_b)):.3f}",
                        f"junc_jaccard={jaccard(set(junc_a), set(junc_b)):.3f}",
                        f"a_sub_b={is_subsequence(cluster_a, cluster_b)}",
                        f"b_sub_a={is_subsequence(cluster_b, cluster_a)}",
                    )
                else:
                    missing = [iso for iso in debug_isoforms[:2] if iso not in debug_lookup]
                    print(f"  Debug comparison skipped; missing: {missing}")

    # Track gene order for trans-splicing and compute panel offsets.
    genes_in_plot = OrderedDict()
    genes_in_plot[target_gene] = None
    for iso in plotted:
        for token in iso['tokens']:
            gene_id, _, _ = split_token(token, target_gene)
            if gene_id not in genes_in_plot:
                genes_in_plot[gene_id] = None

    gene_offsets = {}
    x_offset = 0.0
    for gene in genes_in_plot.keys():
        gene_map = gene_maps.get(gene)
        if not gene_map or not gene_map['segments']:
            continue
        gene_offsets[gene] = x_offset
        x_offset += gene_map['total_length'] + gene_gap
    total_width = max(1.0, x_offset - gene_gap)

    # Default isoform spacing to the group spacing when not explicitly provided.
    if isoform_gap is None:
        isoform_gap = group_gap
    if cluster_isoforms:
        isoform_gap = 0.0

    # Expand isoforms into per-read rows with per-isoform colors.
    colors = plt.get_cmap('tab20').colors
    if isoform_colors is None:
        if cluster_isoforms:
            isoform_colors = assign_cluster_colors(grouped_structures, colors)
        else:
            isoform_colors = assign_isoform_colors(grouped_structures, colors)
    rows = []
    y_positions = []
    y = 0.0
    for group_idx, group in enumerate(grouped_structures):
        for structure in group:
            for iso in structure['items']:
                iso_key = iso.get('isoform_id', iso['resolved_id'])
                iso_count = iso['count']
                if isoform_counts is not None:
                    iso_count = isoform_counts.get(iso_key, 0)
                read_count = coerce_read_count(iso_count)
                if read_count <= 0:
                    continue
                for row_idx in range(read_count):
                    rows.append({
                        'isoform': iso,
                        'color': isoform_colors.get(iso_key, colors[0]),
                        'group_idx': group_idx,
                        'group_label': group_labels.get(group_idx),
                        'row_idx': row_idx,
                        'read_count': read_count
                    })
                    y_positions.append(y)
                    y += row_height + row_gap
                y += isoform_gap
        y += max(0.0, group_gap - isoform_gap)
    if rows:
        y -= max(0.0, group_gap - isoform_gap)
        y -= isoform_gap
    total_height = max(1.0, y)

    # Scale figure height to the total stack height for consistent spacing.
    fig_height = max(1.5, total_height * 0.35)
    fig_width = max(10.0, total_width / 3000.0)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    row_index = 0
    group_label_y = {}
    if cluster_isoforms and rows:
        group_min = {}
        group_max = {}
        for idx, row in enumerate(rows):
            center = y_positions[idx] + (row_height / 2.0)
            group_idx = row['group_idx']
            group_min[group_idx] = min(center, group_min.get(group_idx, center))
            group_max[group_idx] = max(center, group_max.get(group_idx, center))
        for group_idx in group_min:
            group_label_y[group_idx] = (group_min[group_idx] + group_max[group_idx]) / 2.0
    isoform_first_row = {}
    group_first_row = {}
    exon_labels_drawn = False
    debug_transcript_set = set(debug_transcripts or [])
    debug_transcript_seen = set()
    effective_label_mode = label_mode
    if cluster_isoforms and label_mode == 'all':
        effective_label_mode = 'first'
    for row in rows:
        iso = row['isoform']
        is_first_row = row_index == 0
        # Center the read within its row height.
        y_center = y_positions[row_index] + (row_height / 2.0)
        row_index += 1
        segments = build_isoform_segments(iso['tokens'], exon_lookup, target_gene)
        if iso['resolved_id'] in debug_transcript_set and iso['resolved_id'] not in debug_transcript_seen:
            debug_transcript_seen.add(iso['resolved_id'])
            print(f"Debug transcript {iso['resolved_id']}: tokens={'|'.join(iso['tokens'])}")
            for seg in segments:
                print(
                    "  segment",
                    seg['gene'],
                    seg['type'],
                    seg['start'],
                    seg['end'],
                    seg.get('label')
                )
        span_min = None
        span_max = None
        intron_spans = []
        exon_spans = []
        for seg in segments:
            gene = seg['gene']
            gene_map = gene_maps.get(gene)
            if not gene_map:
                continue
            x0 = map_coord(gene_map, seg['start']) + gene_offsets.get(gene, 0.0)
            x1 = map_coord(gene_map, seg['end']) + gene_offsets.get(gene, 0.0)
            if x1 < x0:
                x0, x1 = x1, x0
            if span_min is None or x0 < span_min:
                span_min = x0
            if span_max is None or x1 > span_max:
                span_max = x1
            if seg['type'] == 'I':
                intron_spans.append({
                    'x0': x0,
                    'x1': x1
                })
            else:
                exon_spans.append({
                    'x0': x0,
                    'x1': x1,
                    'label': seg['label']
                })
        # Draw a thin backbone across the read span.
        if span_min is not None and span_max is not None:
            ax.plot(
                [span_min, span_max],
                [y_center, y_center],
                color='lightgrey',
                linewidth=0.01,
                zorder=1
            )
        # Draw intron retention boxes behind exons.
        ir_face = row['color']
        ir_alpha = None
        for span in intron_spans:
            ax.add_patch(
                plt.Rectangle(
                    (span['x0'], y_center - (row_height / 2.0)),
                    max(1e-6, span['x1'] - span['x0']),
                    row_height,
                    facecolor=ir_face,
                    edgecolor='none',
                    linewidth=0,
                    alpha=ir_alpha,
                    zorder=1.5
                )
            )
        # Draw exon boxes on top of the backbone.
        for span in exon_spans:
            ax.add_patch(
                plt.Rectangle(
                    (span['x0'], y_center - (row_height / 2.0)),
                    max(1e-6, span['x1'] - span['x0']),
                    row_height,
                    facecolor=row['color'],
                    edgecolor='white',
                    linewidth=0.000000,
                    zorder=2
                )
            )
        # Add exon labels only once for the first isoform row shown.
        if is_first_row and not exon_labels_drawn:
            label_y = y_center - row_height * 0.7
            seen_labels = set()
            for span in exon_spans:
                label = span['label']
                if label and label not in seen_labels:
                    seen_labels.add(label)
                    ax.text(
                        (span['x0'] + span['x1']) / 2.0,
                        label_y,
                        label,
                        va='bottom',
                        ha='center',
                        fontsize=5,
                        color='black',
                        zorder=3
                    )
            exon_labels_drawn = True
        if effective_label_mode != 'none':
            if effective_label_mode == 'all':
                show_label = True
                label_y = y_center
            else:
                if cluster_isoforms:
                    show_label = row['row_idx'] == 0 and group_first_row.get(row['group_idx'], True)
                    if show_label:
                        label_y = group_label_y.get(row['group_idx'], y_center)
                else:
                    show_label = row['row_idx'] == 0 and isoform_first_row.get(iso['resolved_id'], True)
                    if show_label:
                        block_height = row_height * row['read_count'] + row_gap * (row['read_count'] - 1)
                        label_y = y_positions[row_index - 1] + (block_height / 2.0)
            if show_label:
                label_text = iso['resolved_id']
                if cluster_isoforms:
                    label_text = row.get('group_label') or label_text
                # Center the isoform label over its read stack.
                ax.text(-total_width * 0.01, label_y, label_text,
                        va='center', ha='right', fontsize=8, color=row['color'])
                if cluster_isoforms:
                    group_first_row[row['group_idx']] = False
                else:
                    isoform_first_row[iso['resolved_id']] = False

    for gene, offset in gene_offsets.items():
        gene_map = gene_maps.get(gene)
        if not gene_map:
            continue
        center = offset + gene_map['total_length'] / 2.0
        ax.text(center, -row_height * 4.0, gene, ha='center', va='bottom', fontsize=9)

    # Tight vertical framing around the read stack.
    ax.set_ylim(-row_height * 5.0, total_height)
    ax.set_xlim(-total_width * 0.1, total_width + total_width * 0.02)
    ax.invert_yaxis()
    ax.axis('off')

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

    if missing_structures:
        print(f"Skipped {missing_structures} isoforms with no transcript_associations mapping.")
    if not return_cluster_labels:
        return plotted
    cluster_label_map = {}
    cluster_index_map = {}
    for group_idx, group in enumerate(grouped_structures):
        label = group_labels.get(group_idx)
        for structure in group:
            for iso in structure['items']:
                iso_key = iso.get('isoform_id', iso.get('resolved_id'))
                resolved_id = iso.get('resolved_id')
                if cluster_isoforms:
                    if iso_key:
                        cluster_label_map[iso_key] = label
                        cluster_index_map[iso_key] = group_idx
                    if resolved_id:
                        cluster_label_map[resolved_id] = label
                        cluster_index_map[resolved_id] = group_idx
                else:
                    if iso_key and resolved_id:
                        cluster_label_map[iso_key] = resolved_id
                        cluster_index_map[iso_key] = group_idx
    return plotted, cluster_label_map, cluster_index_map


def main():
    parser = argparse.ArgumentParser(
        description="Plot single-gene long-read isoform structures from an isoform matrix."
    )
    parser.add_argument('--h5ad', nargs='+', required=True,
                        help='Isoform h5ad file(s) or 10x mtx folder(s)')
    parser.add_argument('--gene', required=True, help='Target gene (e.g., ENSG...)')
    parser.add_argument('--transcript-associations', required=True,
                        help='Path to gff-output/transcript_associations.txt')
    parser.add_argument('--gene-model', required=True,
                        help='Path to gene_model TSV reference')
    parser.add_argument('--out', default='isoform_structures.pdf', help='Output figure path')
    parser.add_argument('--groupby', default=None,
                        help='Cell annotation column for filtering')
    parser.add_argument('--group', nargs='*', default=None,
                        help='Optional values in groupby column to keep')
    parser.add_argument('--max-isoforms', type=int, default=200,
                        help='Maximum isoforms to display')
    parser.add_argument('--min-count', type=float, default=1,
                        help='Minimum summed count to include an isoform')
    parser.add_argument('--cluster-mode', choices=['full', 'base', 'block'], default='block',
                        help='How to group isoform structures for clustering')
    parser.add_argument('--cluster-strategy', choices=['substring', 'subsequence', 'overlap', 'jaccard', 'feature', 'structure', 'count'], default='subsequence',
                        help='How to order isoforms (substring/jaccard similarity, feature splits, structure grouping, or count)')
    parser.add_argument('--cluster-similarity-threshold', '--cluster-overlap-threshold',
                        dest='cluster_similarity_threshold', type=float, default=0.85,
                        help='Minimum similarity for clustering (subsequence when --cluster-isoforms)')
    parser.add_argument('--cluster-features', choices=['tokens', 'junctions', 'both'], default='tokens',
                        help='Features to use for clustering')
    parser.add_argument('--min-split-fraction', type=float, default=0.05,
                        help='Minimum fraction for a structural split (feature clustering)')
    parser.add_argument('--cluster-isoforms', action='store_true',
                        help='Color isoforms by cluster and label using the longest isoform in each cluster')
    parser.add_argument('--no-cluster-isoforms', action='store_false', dest='cluster_isoforms',
                        help='Disable isoform clustering and display isoforms independently')
    parser.add_argument('--inspect-isoforms', nargs='*', default=None,
                        help='Isoform IDs to print structure details for')
    parser.add_argument('--intron-scale', type=float, default=0.2,
                        help='Scale factor for intron lengths (0-1)')
    parser.add_argument('--row-height', type=float, default=0.0125,
                        help='Height of each read row')
    parser.add_argument('--row-gap', type=float, default=0.0,
                        help='Gap between read rows')
    parser.add_argument('--group-gap', type=float, default=0.1,
                        help='Gap between clustered groups')
    parser.add_argument('--isoform-gap', type=float, default=None,
                        help='Gap between isoform stacks (defaults to group-gap)')
    parser.add_argument('--label-mode', choices=['first', 'none', 'all'], default='first',
                        help='When to display isoform labels')
    parser.add_argument('--debug-transcripts', nargs='*', default=None,
                        help='Transcript IDs to print segment debugging for')
    parser.add_argument('--debug-isoforms', nargs='*', default=None,
                        help='Isoform IDs to print clustering diagnostics for')
    parser.set_defaults(cluster_isoforms=True)
    args = parser.parse_args()
    if args.cluster_strategy is None:
        args.cluster_strategy = 'substring' if args.cluster_isoforms else 'jaccard'
    if args.cluster_strategy in ('subsequence', 'overlap'):
        args.cluster_strategy = 'substring'
    if args.cluster_similarity_threshold is None:
        args.cluster_similarity_threshold = 0.85 if args.cluster_isoforms else 0.4

    gene_segments, exon_lookup = read_gene_model(args.gene_model)
    gene_maps = build_gene_maps(gene_segments, args.intron_scale)
    transcript_structures, _ = load_transcript_associations(
        args.transcript_associations, args.gene
    )
    if not transcript_structures:
        raise ValueError(f"No transcript structures found for gene {args.gene}.")

    input_paths = args.h5ad
    if len(input_paths) == 1:
        isoform_records = load_isoform_counts(
            input_paths[0], args.gene, groupby=args.groupby, group_values=args.group
        )
        isoform_records = filter_records_by_structures(isoform_records, transcript_structures)
        if not isoform_records:
            raise ValueError("No isoform counts found for the target gene.")

        cluster_label_map, _plotted_ids = build_cluster_label_map(
            isoform_records,
            transcript_structures,
            args.gene,
            args.cluster_mode,
            args.cluster_features,
            args.cluster_strategy,
            args.cluster_similarity_threshold,
            args.min_split_fraction,
            args.max_isoforms,
            args.min_count,
            args.cluster_isoforms,
            inspect_isoforms=args.inspect_isoforms
        )
        output_path = args.out
        input_path = input_paths[0]
        if os.path.isdir(input_path) or input_path.endswith('.mtx') or input_path.endswith('.mtx.gz'):
            out_dir = args.out
            if args.out.endswith('.pdf'):
                out_dir = os.path.dirname(args.out) or "."
            os.makedirs(out_dir, exist_ok=True)
            if os.path.isdir(input_path):
                basename = os.path.basename(input_path)
            else:
                basename = os.path.basename(os.path.dirname(input_path))
            output_path = os.path.join(out_dir, f"{basename}_isoform_structures.pdf")

        plotted, plot_cluster_label_map, plot_cluster_index_map = plot_isoform_structures(
            isoform_records,
            transcript_structures,
            exon_lookup,
            gene_maps,
            args.cluster_mode,
            args.gene,
            args.intron_scale,
            output_path,
            max_isoforms=args.max_isoforms,
            min_count=args.min_count,
            row_height=args.row_height,
            row_gap=args.row_gap,
            group_gap=args.group_gap,
            isoform_gap=args.isoform_gap,
            cluster_strategy=args.cluster_strategy,
            cluster_similarity_threshold=args.cluster_similarity_threshold,
            cluster_features=args.cluster_features,
            label_mode=args.label_mode,
            min_split_fraction=args.min_split_fraction,
            debug_isoforms=args.debug_isoforms,
            debug_transcripts=args.debug_transcripts,
            transcript_associations_path=args.transcript_associations,
            cluster_isoforms=args.cluster_isoforms,
            return_cluster_labels=True
        )
        plotted_ids = set(plot_cluster_index_map.keys())
        report_isoform_structures(
            isoform_records,
            transcript_structures,
            args.gene,
            args.inspect_isoforms,
            include_introns=args.cluster_isoforms,
            cluster_label_map=cluster_label_map,
            plotted_ids=plotted_ids,
            plot_cluster_label_map=plot_cluster_label_map,
            plot_cluster_index_map=plot_cluster_index_map
        )
        print(f"Saved isoform structure plot to: {output_path}")
        isoform_path, _isoform_rows = export_plotted_isoforms(
            plotted,
            output_path,
            plot_cluster_label_map,
            plot_cluster_index_map
        )
        if isoform_path:
            print(f"Saved plotted isoform IDs to: {isoform_path}")
    else:
        combined_counts = defaultdict(float)
        record_by_isoform = {}
        per_file_counts = {}
        for path in input_paths:
            records = load_isoform_counts(
                path, args.gene, groupby=args.groupby, group_values=args.group
            )
            records = filter_records_by_structures(records, transcript_structures)
            per_file_counts[path] = {rec['isoform_id']: rec['count'] for rec in records}
            for rec in records:
                key = rec['isoform_id']
                combined_counts[key] += rec['count']
                if key not in record_by_isoform:
                    record_by_isoform[key] = rec

        combined_records = []
        for key, rec in record_by_isoform.items():
            merged = dict(rec)
            merged['count'] = combined_counts[key]
            combined_records.append(merged)

        cluster_label_map, _plotted_ids = build_cluster_label_map(
            combined_records,
            transcript_structures,
            args.gene,
            args.cluster_mode,
            args.cluster_features,
            args.cluster_strategy,
            args.cluster_similarity_threshold,
            args.min_split_fraction,
            args.max_isoforms,
            args.min_count,
            args.cluster_isoforms,
            inspect_isoforms=args.inspect_isoforms
        )
        grouped_structures = None
        isoform_colors = None
        if combined_records:
            plotted, missing_structures, _missing_ids = build_plotted_isoforms(
                combined_records,
                transcript_structures,
                args.gene,
                args.cluster_mode,
                args.cluster_features,
                args.min_count,
                args.max_isoforms
            )
            if not plotted:
                raise ValueError("No isoform counts found for the target gene.")
            grouped_structures = group_structures(
                plotted,
                args.gene,
                args.cluster_features,
                args.cluster_strategy,
                args.cluster_similarity_threshold,
                args.min_split_fraction,
                include_introns=args.cluster_isoforms,
                report=args.cluster_isoforms
            )
            if args.cluster_isoforms:
                isoform_colors = assign_cluster_colors(grouped_structures, plt.get_cmap('tab20').colors)
            else:
                isoform_colors = assign_isoform_colors(grouped_structures, plt.get_cmap('tab20').colors)

        out_base = args.out
        out_dir = out_base
        if out_base.endswith('.pdf'):
            out_dir = os.path.dirname(out_base) or "."
        os.makedirs(out_dir, exist_ok=True)

        for path in input_paths:
            if os.path.isdir(path):
                basename = os.path.basename(path)
            elif path.endswith('.mtx') or path.endswith('.mtx.gz'):
                basename = os.path.basename(os.path.dirname(path))
            else:
                basename = os.path.splitext(os.path.basename(path))[0]
            out_path = os.path.join(out_dir, f"{basename}_isoform_structures.pdf")
            plotted, plot_cluster_label_map, plot_cluster_index_map = plot_isoform_structures(
                record_by_isoform.values(),
                transcript_structures,
                exon_lookup,
                gene_maps,
                args.cluster_mode,
                args.gene,
                args.intron_scale,
                out_path,
                max_isoforms=args.max_isoforms,
                min_count=args.min_count,
                row_height=args.row_height,
                row_gap=args.row_gap,
                group_gap=args.group_gap,
                isoform_gap=args.isoform_gap,
                cluster_strategy=args.cluster_strategy,
                cluster_similarity_threshold=args.cluster_similarity_threshold,
                cluster_features=args.cluster_features,
                label_mode=args.label_mode,
                min_split_fraction=args.min_split_fraction,
                debug_isoforms=args.debug_isoforms,
                debug_transcripts=args.debug_transcripts,
                grouped_structures=grouped_structures,
                isoform_colors=isoform_colors,
                isoform_counts=per_file_counts.get(path, {}),
                transcript_associations_path=args.transcript_associations,
                cluster_isoforms=args.cluster_isoforms,
                return_cluster_labels=True
            )
            plotted_ids = set(plot_cluster_index_map.keys())
            report_isoform_structures(
                combined_records,
                transcript_structures,
                args.gene,
                args.inspect_isoforms,
                include_introns=args.cluster_isoforms,
                cluster_label_map=cluster_label_map,
                plotted_ids=plotted_ids,
                plot_cluster_label_map=plot_cluster_label_map,
                plot_cluster_index_map=plot_cluster_index_map
            )
            print(f"Saved isoform structure plot to: {out_path}")
            isoform_path, _isoform_rows = export_plotted_isoforms(
                plotted,
                out_path,
                plot_cluster_label_map,
                plot_cluster_index_map
            )
            if isoform_path:
                print(f"Saved plotted isoform IDs to: {isoform_path}")


if __name__ == '__main__':
    main()
