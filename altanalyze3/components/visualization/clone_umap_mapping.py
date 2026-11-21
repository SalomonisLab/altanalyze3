"""Map clones to saved coordinates using exact library and reverse-complement barcode; preserve both complete rosters."""
from pathlib import Path
import pandas as pd, numpy as np, json, hashlib, shutil

def prepare(coordinates, figures, config_path, output):
    I = Path(output)
    I.mkdir(parents=True, exist_ok=True)
    F = Path(figures)
    config = json.loads(Path(config_path).read_text())
    u = pd.read_csv(Path(coordinates), sep='\t', header=None, names=['umap_id', 'UMAP1', 'UMAP2'])
    u['source_order'] = np.arange(len(u))
    if not (u.umap_id.is_unique and np.isfinite(u[['UMAP1', 'UMAP2']]).all().all()):
        raise ValueError("Integrity check failed: u.umap_id.is_unique and np.isfinite(u[['UMAP1', 'UMAP2']]).all().all()")
    p = u.umap_id.str.extract('^(?P<illumina_barcode>[ACGT]{16})-1\\.(?P<library>.+)$')
    if not p.notna().all().all():
        raise ValueError('Integrity check failed: p.notna().all().all()')
    u = pd.concat([u, p], axis=1)
    c = pd.read_csv(F / 'cell_assignments.tsv.gz', sep='\t')
    o = pd.read_csv(F / 'cell_order.tsv.gz', sep='\t')
    k = pd.read_csv(F / 'clone_key.tsv', sep='\t')
    if not (c.cell_id.is_unique and o.cell_id.is_unique and (set(c.cell_id) == set(o.cell_id))):
        raise ValueError('Integrity check failed: c.cell_id.is_unique and o.cell_id.is_unique and (set(c.cell_id) == set(o.cell_id))')
    c = c.merge(o[['cell_id', 'clone']], on='cell_id', validate='one_to_one', sort=False)
    if not c.clone.eq(c.node.map(k.set_index('node').clone)).all():
        raise ValueError("Integrity check failed: c.clone.eq(c.node.map(k.set_index('node').clone)).all()")
    if not c.barcode_rc.eq(c.barcode.map(lambda x: x.translate(str.maketrans('ACGT', 'TGCA'))[::-1])).all():
        raise ValueError("Integrity check failed: c.barcode_rc.eq(c.barcode.map(lambda x: x.translate(str.maketrans('ACGT', 'TGCA'))[::-1])).all()")
    if not not c.duplicated(['library', 'barcode_rc']).any():
        raise ValueError("Integrity check failed: not c.duplicated(['library', 'barcode_rc']).any()")
    cols = ['cell_id', 'sample', 'library', 'barcode', 'barcode_rc', 'clone', 'stage', 'accepted_MUT_any', 'max_reporting_membership', 'cell_state']
    m = u.merge(c[cols], left_on=['library', 'illumina_barcode'], right_on=['library', 'barcode_rc'], how='left', validate='one_to_one', sort=False)
    if not (m.umap_id.tolist() == u.umap_id.tolist() and np.array_equal(m[['UMAP1', 'UMAP2']], u[['UMAP1', 'UMAP2']])):
        raise ValueError("Integrity check failed: m.umap_id.tolist() == u.umap_id.tolist() and np.array_equal(m[['UMAP1', 'UMAP2']], u[['UMAP1', 'UMAP2']])")
    if not m.cell_id.dropna().is_unique:
        raise ValueError('Integrity check failed: m.cell_id.dropna().is_unique')
    healthy = m.library.isin(config.get('healthy_libraries', []))
    matched = m.cell_id.notna()
    if not not (healthy & matched).any():
        raise ValueError('Healthy library also has patient clone assignments')
    m['match_status'] = np.select([matched, healthy, m.library.isin(c.library)], ['matched', 'healthy_control_no_patient_assignment', 'same_library_no_KINNEX_assignment'], default='library_not_in_KINNEX')
    m['plot_clone'] = m.clone.fillna('Unmatched')
    m.loc[healthy, 'plot_clone'] = 'Aged healthy'
    stage_map = config.get('library_timepoints', {})
    m['umap_timepoint'] = m.library.map(stage_map).fillna('Unmapped library')
    if config.get('require_all_library_timepoints', True):
        if not set(m.library) <= set(stage_map):
            raise ValueError('Missing library timepoint mapping')
    stage_normal = config.get('stage_normalization', {'sAML CR': 'sAML-CR'})
    if not m.loc[matched, 'stage'].replace(stage_normal).eq(m.loc[matched, 'umap_timepoint']).all():
        raise ValueError("Integrity check failed: m.loc[matched, 'stage'].replace(stage_normal).eq(m.loc[matched, 'umap_timepoint']).all()")
    summary = m.groupby(['library', 'umap_timepoint'], sort=False).agg(umap_cells=('umap_id', 'size'), matched_cells=('cell_id', 'count')).reset_index()
    summary['unmatched_cells'] = summary.umap_cells - summary.matched_cells
    summary['matched_percent'] = 100 * summary.matched_cells / summary.umap_cells
    status = m.groupby('match_status').size().to_dict()
    reverse = c.merge(m.loc[matched, ['cell_id', 'umap_id', 'UMAP1', 'UMAP2']], on='cell_id', how='left', validate='one_to_one', sort=False)
    if not (len(reverse) == len(c) and set(reverse.cell_id) == set(c.cell_id)):
        raise ValueError('Integrity check failed: len(reverse) == len(c) and set(reverse.cell_id) == set(c.cell_id)')
    reverse['has_saved_UMAP_coordinate'] = reverse.umap_id.notna()
    reverse_summary = reverse.groupby(['sample', 'library', 'stage'], sort=False).agg(KINNEX_cells=('cell_id', 'size'), cells_on_saved_UMAP=('has_saved_UMAP_coordinate', 'sum')).reset_index()
    counts = pd.crosstab(m.loc[matched, 'umap_timepoint'], m.loc[matched, 'clone']).reindex(columns=k.clone, fill_value=0)
    pct = counts.div(counts.sum(axis=1), axis=0) * 100
    m.to_csv(I / 'mapped_umap_cells.tsv.gz', sep='\t', index=False)
    reverse.to_csv(I / 'all_KINNEX_mapping.tsv.gz', sep='\t', index=False)
    summary.to_csv(I / 'matching_by_library.tsv', sep='\t', index=False)
    reverse_summary.to_csv(I / 'matching_by_KINNEX_sample.tsv', sep='\t', index=False)
    counts.to_csv(I / 'clone_counts_among_matched_UMAP_cells.tsv', sep='\t')
    pct.to_csv(I / 'clone_percent_among_matched_UMAP_cells.tsv', sep='\t')
    k.to_csv(I / 'clone_key.tsv', sep='\t', index=False)
    labels = k.set_index('clone').mutation_combination.to_dict()
    levels = ['Unmatched', 'Aged healthy'] + k.clone.tolist()
    group = pd.DataFrame({'cell': m.umap_id, 'group': m.plot_clone.map({x: i + 1 for i, x in enumerate(levels)}), 'label': m.plot_clone.map(lambda x: x + ' = ' + labels[x] if x in labels else x)})
    if not len(group) == len(u):
        raise ValueError('Integrity check failed: len(group) == len(u)')
    group.to_csv(I / 'groups.KINNEX_clones.txt', sep='\t', header=False, index=False)
    files = [Path(coordinates), Path(config_path), F / 'cell_assignments.tsv.gz', F / 'cell_order.tsv.gz', F / 'clone_key.tsv']
    manifest = {'sources': [{'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in files], 'matching_rule': 'Exact library + exact reverse-complement barcode; no fuzzy or cross-library matching, reclustering or UMAP refitting.', 'UMAP_cells': len(u), 'matched_cells': int(matched.sum()), 'KINNEX_cells_retained': len(c), 'status_counts': {x: int(v) for x, v in status.items()}, 'all_coordinates_and_order_unchanged': True, 'all_clone_assignments_unchanged': True, 'all_input_cells_retained': True, 'render_inputs_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in I.iterdir() if p.is_file() and p.name != 'mapping_validation.json'}}
    (I / 'mapping_validation.json').write_text(json.dumps(manifest, indent=2))
    return manifest

def add_arguments(p):
    p.add_argument('--coordinates', required=True)
    p.add_argument('--figures', required=True)
    p.add_argument('--config', required=True)
    p.add_argument('--output', required=True)

def run(a):
    return prepare(a.coordinates, a.figures, a.config, a.output)

def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    add_arguments(p)
    return run(p.parse_args(argv))
if __name__ == '__main__':
    main()
