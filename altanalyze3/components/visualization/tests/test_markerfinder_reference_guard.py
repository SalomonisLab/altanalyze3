from copy import deepcopy
import pytest
from altanalyze3.components.visualization.markerfinder_reference_guard import (
    ReferenceMismatch, file_sha256, validate_preflight, validate_artifacts,
)


def fixture(tmp_path):
    p = tmp_path / 'verified.py'
    p.write_text('original code\n')
    pins = {str(p): file_sha256(p)}
    c = dict(unresolved=[], historical_code_verified=True,
             source_samples=['a', 'b'], source_features=['e1', 'e2'],
             patterns=['up', 'down'], callout_config={'shared_exons': True, 'color': 'red'},
             parameters={'rho': .4, 'top_n': 100, 'missing': 'observed', 'coverage': .7},
             authorized_changes={'missing': 'row_median', 'coverage': 0},
             reference_files=pins, execution_code=pins)
    run = {k: deepcopy(c[k]) for k in ('source_samples', 'source_features', 'patterns', 'callout_config')}
    run.update(parameters={**c['parameters'], **c['authorized_changes']}, code=pins)
    return c, run, p


def test_authorized_only(tmp_path):
    c, run, _ = fixture(tmp_path)
    validate_preflight(c, run)


@pytest.mark.parametrize('key,value', [
    ('source_samples', ['a', 'wrong']), ('source_features', ['e1', 'wrong']),
    ('patterns', ['up']), ('patterns', ['down', 'up']),
    ('callout_config', {'shared_exons': False, 'color': 'red'}),
])
def test_same_count_substitutions_and_dropped_patterns_fail(tmp_path, key, value):
    c, run, _ = fixture(tmp_path)
    run[key] = value
    with pytest.raises(ReferenceMismatch): validate_preflight(c, run)


@pytest.mark.parametrize('key,value', [('rho', .5), ('top_n', 20), ('unexpected_filter', True)])
def test_unrequested_setting_fails(tmp_path, key, value):
    c, run, _ = fixture(tmp_path)
    run['parameters'][key] = value
    with pytest.raises(ReferenceMismatch): validate_preflight(c, run)


def test_changed_code_fails(tmp_path):
    c, run, p = fixture(tmp_path)
    p.write_text('different code\n')
    with pytest.raises(ReferenceMismatch): validate_preflight(c, run)


def test_unresolved_provenance_blocks_even_if_hashes_match(tmp_path):
    c, run, _ = fixture(tmp_path)
    c['unresolved'] = ['Historical native module not recovered']
    with pytest.raises(ReferenceMismatch): validate_preflight(c, run)


def test_missing_callout_and_same_count_wrong_event_fail():
    expected = dict(samples=['a', 'b'], events=['e1', 'e2'], event_patterns=['up', 'down'],
                    callout_rows=[{'event': 'e2', 'label': 'CEP57', 'color': 'red'}],
                    render_settings={'range': 1.5})
    validate_artifacts(expected, deepcopy(expected))
    for key, value in [('callout_rows', []), ('events', ['e1', 'wrong']),
                       ('event_patterns', ['down', 'up']), ('render_settings', {'range': .2})]:
        actual = deepcopy(expected)
        actual[key] = value
        with pytest.raises(ReferenceMismatch): validate_artifacts(expected, actual)
