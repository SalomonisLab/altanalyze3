"""The SNAF-T collapsed summary must reproduce the validated Python-2 reference.

`altanalyze3/components/snaf/t_summary.py` ports `updatedSNAFTMerge` from AltAnalyze's
`stats_scripts/SNAFintegration.py`. These tests pin the two things a port can quietly get
wrong: the per-(HLA, sample) de-duplication the reference gets for free from its dicts, and
the EventAnnotation join read by column POSITION. They also pin the direction of each added
metric, because a "best binder" computed as a maximum would be silently backwards.
"""
from __future__ import annotations

import os
import sys

import pytest

pytest.importorskip('pandas')
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from altanalyze3.components.snaf.t_summary import (  # noqa: E402
    EA_POS, import_event_annotations, summarize_t_candidates)

CANDIDATE_HEADER = ['sample', 'peptide', 'uid', 'junction_count', 'phase', 'evidences', 'hla',
                    'binding_affinity', 'immunogenicity', 'tumor_specificity_mean',
                    'tumor_specificity_mle', 'n_sample', 'coord', 'symbol', 'in_db']

J1 = 'ENSG00000000419:E6.1-E7.1'
EV = "((130611533, 'ENST00000409602', '+'),)"


def _row(sample, peptide, uid, hla, binding, immuno, jc='50', n_sample='3'):
    return [sample, peptide, uid, jc, '1', EV, hla, binding, immuno,
            '0.26', '0.13', n_sample, 'chr20:50942031-50941209', 'DPM1', 'True']


def _write_candidates(tmp_path, rows):
    d = tmp_path / 'T_candidates'
    d.mkdir(parents=True, exist_ok=True)
    p = d / 'T_antigen_candidates_all.txt'
    with open(p, 'w') as fh:
        fh.write('\t'.join(CANDIDATE_HEADER) + '\n')
        for r in rows:
            fh.write('\t'.join(r) + '\n')
    return str(d)


def _write_event_annotation(tmp_path):
    p = tmp_path / 'ea.txt'
    head = ['Symbol', 'Description', 'Examined-Junction', 'Background-Major-Junction',
            'AltExons', 'ProteinPredictions', 'dPSI', 'ClusterID', 'UID', 'Coordinates',
            'EventAnnotation', 'S1.bed', 'S2.bed']
    row = ['DPM1', 'dolichyl-phosphate', J1, 'ENSG00000000419:E6.1-E8.2',
           'ENSG00000000419:E7.1', '(+)alt-coding', '0.2154', 'clu_4',
           'DPM1:' + J1, 'chr20:50942031-50941209', 'cassette-exon', '0.1', '0.2']
    with open(p, 'w') as fh:
        fh.write('\t'.join(head) + '\n')
        fh.write('\t'.join(row) + '\n')
    return str(p)


def _blank(v):
    """An unfilled cell: the TSV holds an empty field, which pandas reads back as NaN."""
    return pd.isna(v) or str(v).strip() == ''


def _summary(tmp_path, rows, **kw):
    d = _write_candidates(tmp_path, rows)
    out = summarize_t_candidates(d, **kw)
    return pd.read_csv(out, sep='\t', dtype=str)


# --- the reference's aggregation -----------------------------------------------------

def test_one_row_per_peptide_and_junction(tmp_path):
    rows = [_row('S1.bed', 'SNPEQDLKL', J1, 'HLA-C*03:04', '1.0', '0.70'),
            _row('S2.bed', 'SNPEQDLKL', J1, 'HLA-C*03:04', '3.0', '0.90'),
            _row('S1.bed', 'OTHERPEPT', J1, 'HLA-C*03:04', '2.0', '0.50')]
    df = _summary(tmp_path, rows)
    assert len(df) == 2
    r = df[df['Peptide'] == 'SNPEQDLKL'].iloc[0]
    assert r['Num Samples'] == '2'
    assert set(r['Samples'].split('|')) == {'S1', 'S2'}, 'the .bed suffix must be stripped'
    assert float(r['Mean Binding']) == pytest.approx(2.0)
    assert float(r['Mean Immunogenicity']) == pytest.approx(0.80)


def test_a_repeated_hla_sample_pair_is_counted_once(tmp_path):
    """The reference keys its dicts on (hla, sample), so a duplicate row contributes once.

    Grouping the raw rows instead would weight the duplicate twice and move the mean.
    """
    rows = [_row('S1.bed', 'PEPTIDEAA', J1, 'HLA-A*02:01', '1.0', '0.10'),
            _row('S1.bed', 'PEPTIDEAA', J1, 'HLA-A*02:01', '1.0', '0.10'),
            _row('S1.bed', 'PEPTIDEAA', J1, 'HLA-B*07:02', '5.0', '0.50')]
    r = _summary(tmp_path, rows).iloc[0]
    assert r['Num Observations'] == '2', 'the duplicate pair must collapse'
    assert float(r['Mean Binding']) == pytest.approx(3.0)   # (1+5)/2, not (1+1+5)/3
    assert r['Num Samples'] == '1'


def test_in_frame_follows_the_reference_rule(tmp_path):
    """The reference sets In-Frame from the LENGTH of the evidence string, not the phase."""
    long_ev = _row('S1.bed', 'PEPTIDEAA', J1, 'HLA-A*02:01', '1.0', '0.1')
    short_ev = _row('S2.bed', 'PEPTIDEBB', J1, 'HLA-A*02:01', '1.0', '0.1')
    short_ev[5] = '()'
    df = _summary(tmp_path, [long_ev, short_ev])
    assert df.set_index('Peptide').loc['PEPTIDEAA', 'In-Frame'] == 'True'
    assert df.set_index('Peptide').loc['PEPTIDEBB', 'In-Frame'] == 'UNK'


def test_the_peptide_id_matches_the_reference_transform(tmp_path):
    r = _summary(tmp_path, [_row('S1.bed', 'SNPEQDLKL', J1, 'HLA-A*02:01', '1.0', '0.1')]).iloc[0]
    assert r['peptideID'] == 'SNPEQDLKL.ENSG00000000419.E6.1.E7.1'


def test_every_reference_column_is_present_and_in_order(tmp_path):
    from altanalyze3.components.snaf.t_summary import REFERENCE_COLUMNS
    df = _summary(tmp_path, [_row('S1.bed', 'PEPTIDEAA', J1, 'HLA-A*02:01', '1.0', '0.1')])
    assert list(df.columns)[:len(REFERENCE_COLUMNS)] == REFERENCE_COLUMNS


# --- the added metrics ---------------------------------------------------------------

def test_best_binding_is_the_minimum_and_names_its_allele(tmp_path):
    """binding_affinity is a percentile: LOWER binds more strongly."""
    rows = [_row('S1.bed', 'PEPTIDEAA', J1, 'HLA-A*02:01', '9.0', '0.20'),
            _row('S1.bed', 'PEPTIDEAA', J1, 'HLA-B*07:02', '0.5', '0.80'),
            _row('S1.bed', 'PEPTIDEAA', J1, 'HLA-C*03:04', '4.0', '0.40')]
    r = _summary(tmp_path, rows).iloc[0]
    assert float(r['Best Binding']) == pytest.approx(0.5)
    assert r['Best Binding HLA'] == 'HLA-B*07:02'
    assert float(r['Mean Binding']) == pytest.approx(4.5)
    assert float(r['Max Immunogenicity']) == pytest.approx(0.80)
    assert r['Num HLA Alleles'] == '3'
    assert r['HLA Alleles'] == 'HLA-A*02:01|HLA-B*07:02|HLA-C*03:04'
    assert r['Peptide Length'] == '9'


def test_sample_frequency_uses_the_cohort_denominator(tmp_path):
    rows = [_row('S1.bed', 'PEPTIDEAA', J1, 'HLA-A*02:01', '1.0', '0.1'),
            _row('S2.bed', 'PEPTIDEAA', J1, 'HLA-A*02:01', '1.0', '0.1'),
            _row('S3.bed', 'OTHERPEPT', J1, 'HLA-A*02:01', '1.0', '0.1')]
    r = _summary(tmp_path, rows).set_index('Peptide').loc['PEPTIDEAA']
    assert float(r['Sample Frequency']) == pytest.approx(2 / 3, abs=1e-4)  # written at 4 dp
    d = _write_candidates(tmp_path / 'b', rows)
    r2 = pd.read_csv(summarize_t_candidates(d, cohort_size=100), sep='\t',
                     dtype=str).set_index('Peptide').loc['PEPTIDEAA']
    assert float(r2['Sample Frequency']) == pytest.approx(0.02, abs=1e-4)


def test_max_junction_count_beats_the_arbitrary_first(tmp_path):
    rows = [_row('S1.bed', 'PEPTIDEAA', J1, 'HLA-A*02:01', '1.0', '0.1', jc='10'),
            _row('S2.bed', 'PEPTIDEAA', J1, 'HLA-A*02:01', '1.0', '0.1', jc='90')]
    r = _summary(tmp_path, rows).iloc[0]
    assert r['JunctionCount'] == '10', 'the reference keeps whichever row came first'
    assert float(r['Max Junction Count']) == pytest.approx(90)


def test_rows_are_ordered_by_sharing_then_by_best_binder(tmp_path):
    rows = [_row('S1.bed', 'RAREPEPTI', J1, 'HLA-A*02:01', '0.1', '0.1'),
            _row('S1.bed', 'COMMONPEP', J1, 'HLA-A*02:01', '5.0', '0.1'),
            _row('S2.bed', 'COMMONPEP', J1, 'HLA-A*02:01', '5.0', '0.1')]
    assert list(_summary(tmp_path, rows)['Peptide']) == ['COMMONPEP', 'RAREPEPTI']


# --- the EventAnnotation join --------------------------------------------------------

def test_event_annotation_is_read_by_position(tmp_path):
    ea = import_event_annotations(_write_event_annotation(tmp_path))
    assert set(ea) == {J1}
    assert ea[J1]['event_annotation'] == 'cassette-exon'
    assert ea[J1]['dpsi'] == '0.2154'
    assert ea[J1]['protein_predictions'] == '(+)alt-coding'
    assert EA_POS['event_annotation'] == 10, 'position is the contract; the header grows per sample'


def test_the_splice_event_reaches_the_summary(tmp_path):
    rows = [_row('S1.bed', 'PEPTIDEAA', J1, 'HLA-A*02:01', '1.0', '0.1')]
    r = _summary(tmp_path, rows, event_annotation=_write_event_annotation(tmp_path)).iloc[0]
    assert r['EventAnnotation'] == 'cassette-exon'
    assert r['ClusterID'] == 'clu_4'
    assert r['UID'] == 'DPM1:' + J1
    # the four the reference parses and then drops
    assert r['ProteinPredictions'] == '(+)alt-coding'
    assert r['dPSI'] == '0.2154'
    assert r['AltExons'] == 'ENSG00000000419:E7.1'
    assert r['Description'] == 'dolichyl-phosphate'


def test_a_junction_absent_from_the_annotation_is_blank_not_dropped(tmp_path):
    rows = [_row('S1.bed', 'PEPTIDEAA', 'ENSG00000999999:E1.1-E2.1', 'HLA-A*02:01', '1.0', '0.1')]
    df = _summary(tmp_path, rows, event_annotation=_write_event_annotation(tmp_path))
    assert len(df) == 1
    assert _blank(df.iloc[0]['EventAnnotation'])
    assert df.iloc[0]['Symbol'] == 'DPM1', 'the candidate symbol survives'


def test_without_an_annotation_the_run_still_writes_every_candidate(tmp_path):
    df = _summary(tmp_path, [_row('S1.bed', 'PEPTIDEAA', J1, 'HLA-A*02:01', '1.0', '0.1')])
    assert len(df) == 1 and _blank(df.iloc[0]['EventAnnotation'])


# --- guards --------------------------------------------------------------------------

def test_a_missing_candidate_directory_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        summarize_t_candidates(str(tmp_path / 'nothing'))


def test_a_table_missing_the_snaf_columns_raises(tmp_path):
    d = tmp_path / 'T_candidates'
    d.mkdir()
    (d / 'T_antigen_candidates_all.txt').write_text('a\tb\n1\t2\n')
    with pytest.raises(ValueError, match='candidate columns'):
        summarize_t_candidates(str(d))


def test_the_sample_filter_restricts_the_cohort(tmp_path):
    rows = [_row('S1.bed', 'PEPTIDEAA', J1, 'HLA-A*02:01', '1.0', '0.1'),
            _row('S2.bed', 'PEPTIDEAA', J1, 'HLA-A*02:01', '1.0', '0.1')]
    r = _summary(tmp_path, rows, sample_filter=['S1']).iloc[0]
    assert r['Num Samples'] == '1' and r['Samples'] == 'S1'
