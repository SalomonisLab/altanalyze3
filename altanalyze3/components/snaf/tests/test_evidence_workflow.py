"""Integration of ranked isoforms with SNAF-B surface checks and reports."""
import csv
import json

import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.snaf.surface import evidence_isoform as E
from altanalyze3.components.snaf.surface import evidence_workflow as W
from altanalyze3.components.snaf.surface import main as M


@pytest.fixture
def reference(monkeypatch):
    gene = 'ENSG00000000001'
    uid = gene + ':I1.1_1000-E2.1'
    coords = {gene: {'E1.1': ('chr1', '+', 2001, 2101),
                    'I1.1': ('chr1', '+', 2102, 2499),
                    'E2.1': ('chr1', '+', 1500, 1649)}}
    # Deliberately use a downstream reference first exon: the selected novel donor
    # cannot inherit a 5' boundary and must get 250 bases.
    coords[gene]['E1.1'] = ('chr1', '+', 100, 200)
    coords[gene]['I1.1'] = ('chr1', '+', 201, 1499)
    exonlist = pd.DataFrame([[gene, 'ENST1', 'ENSP1', 'E1.1|E2.1']],
                           columns=['EnsGID', 'EnsTID', 'EnsPID', 'Exons'])
    mrna = 'C' * 200 + 'ATG' + 'GCC' * 49 + 'TAA' + 'C' * 47
    sequence = list('C' * 2000)
    sequence[750:1000] = mrna[:250]
    sequence[1499:1649] = mrna[250:]
    # origin=start-2000=1, so the synthetic reference starts at genomic base 1.
    for name, value in dict(df_exonlist=exonlist, dict_exonCoords=coords,
            dict_fa={gene: ['chr1', '2001', '2001', ''.join(sequence)]},
            dict_biotype={gene: {'ENSP1': 'protein_coding'}},
            dict_uni_fa={gene: {'REF': 'M' + 'A' * 39}}).items():
        monkeypatch.setattr(M, name, value, raising=False)
    matrix = pd.DataFrame([[10, 0]], index=[uid], columns=['A', 'B'])
    tuples = [(uid, 0.01, None, None, 0.5)]
    return gene, uid, tuples, matrix, mrna


@pytest.mark.parametrize('method,depth', [('learned', 4), ('evidence', 2)])
def test_surface_run_keeps_ranked_orfs_and_exports_boundaries(reference, tmp_path, method, depth):
    gene, uid, tuples, matrix, mrna = reference
    result = M.run(tuples, str(tmp_path), junction_counts=matrix,
                   isoform_method=method, first_exon_junctions=[uid])
    sa = result[0]
    assert sa.orfp == ['M' + 'A' * 49]
    assert sa.full_length == [mrna]
    assert sa.nmd == ['#'] and sa.translatability == ['#'] and sa.alignment == [True]
    assert sa.synthetic_predictions[0].chain == [(751, 1000), (1500, 1649)]
    with (tmp_path / 'synthetic_isoforms/predictions.tsv').open() as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    assert rows[0]['supporting_samples'] == 'A'
    assert rows[0]['first_exon_boundary_source'] == 'assumed_250nt'
    manifest = json.loads((tmp_path / 'synthetic_isoforms/manifest.json').read_text())
    assert manifest['method'] == method and manifest['max_edits'] == depth
    assert not manifest['long_reads_used_for_inference']
    # Real surface filtering + reporting must not index candidates as reference ENST rows.
    cc, cf = M.generate_results(str(tmp_path / 'surface_antigen_sr.p'), outdir=str(tmp_path))
    assert cc == 1 and cf == 0
    freq = tmp_path / 'frequency.tsv'
    pd.DataFrame({'tumor_specificity_mean': [0.0], 'tumor_specificity_mle': [0.0]},
                 index=[uid + ',' + uid]).to_csv(freq, sep='\t')
    M.report_candidates(str(tmp_path / 'surface_antigen_sr.p'),
        str(tmp_path / 'candidates_3_sr_None_False.txt'),
        str(tmp_path / 'validation_3_sr_None_False.txt'), str(freq), mode='short_read',
        outdir=str(tmp_path / 'reports'), name='report.tsv')
    with (tmp_path / 'reports/report.tsv').open() as fh:
        report = next(csv.DictReader(fh, delimiter='\t'))
    assert report['mode'] == 'synthetic_' + method
    assert report['evidence'] == rows[0]['artifact_id']
    assert report['mRNA_sequence'] == mrna


def test_python_workflow_defaults_to_learned_and_enforces_sample(reference, tmp_path):
    _, uid, tuples, matrix, _ = reference
    result = M.run(tuples, str(tmp_path), junction_counts=matrix,
                   first_exon_junctions=[uid], isoform_sample='B')
    assert result[0].isoform_method == 'learned'
    assert result[0].synthetic_predictions == []


def test_missing_reference_sequence_leaves_event_unresolved(reference, tmp_path):
    gene, uid, tuples, matrix, _ = reference
    M.dict_fa[gene][3] = M.dict_fa[gene][3][:900]
    result = M.run(tuples, str(tmp_path), junction_counts=matrix,
                   first_exon_junctions=[uid])
    assert result[0].synthetic_predictions == []


def test_reference_blocks_merge_and_all_evidence_rows_are_retained(reference):
    gene, uid, _, matrix, _ = reference
    coords = M.dict_exonCoords
    coords[gene]['E1.2'] = ('chr1', '+', 200, 240)
    exonlist = M.df_exonlist.copy()
    exonlist.loc[0, 'Exons'] = 'E1.1|E1.2|E2.1'
    models, _, _ = W.reference_models(exonlist, coords, {gene})
    assert models[gene]['ENST1'][1] == [(100, 240), (1500, 1649)]
    # A background junction need not have passed tumor-specificity sifting.
    matrix.loc[gene + ':E1.2-E2.1'] = [12, 20]
    matrix.loc[uid + '=chr1:1000-1500'] = [8, 5]
    counts, _, skipped = W.matrix_evidence(matrix, coords, {gene})
    assert len(counts[gene]) == 2 and not skipped
    np.testing.assert_array_equal(counts[gene][(1000, 1500)], [10, 5])
    assert W.resolve_junction(gene + ':E1.1-ENSG2:E1.1', coords) is None


def test_gene_fasta_minus_strand_is_exposed_in_forward_orientation():
    assert W._gene_fetch(['chr1', 2001, 2004, 'AAGC'], '-')(1, 4) == 'GCTT'
    assert W._gene_fetch(['chr1', 2001, 2004, 'AAGC'], '-')(0, 3) == ''


def test_bundled_model_and_schema_validation(tmp_path):
    model = E.load_ranker()
    assert np.isfinite(model.predict(np.zeros((2, len(E.FEATURE_NAMES))))).all()
    path = tmp_path / 'bad.json'
    path.write_text(json.dumps({'feature_names': ['wrong']}))
    with pytest.raises(ValueError, match='schema'):
        E.load_ranker(path)


def test_no_validation_catalog_never_calls_long_read_support(monkeypatch, tmp_path):
    calls = []
    def generate(**kwargs):
        calls.append(kwargs)
        return 0, 0
    def forbidden(*args, **kwargs):
        pytest.fail('No long-read validation should run without an independent catalog')
    monkeypatch.setattr(M, 'generate_results', generate)
    monkeypatch.setattr(M, 'prewarm_support_cache', forbidden)
    M.generate_full_results(str(tmp_path), 'unused.tsv', 'short_read', None)
    assert len(calls) == 6
    assert all(c['strigency'] == 3 and c['gtf'] is None for c in calls)
