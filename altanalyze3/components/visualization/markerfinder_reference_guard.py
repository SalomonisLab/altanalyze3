"""Fail-closed checks for reproductions of explicitly named MarkerFinder figures.

This module performs no analysis and writes no scientific outputs. Call preflight
before fitting/rendering; call validate_artifacts on staged output before publishing.
A captured hash alone is not proof that a file generated a historical figure:
unresolved provenance must be recorded in ``unresolved`` and blocks execution.
"""
from hashlib import sha256
from pathlib import Path


class ReferenceMismatch(RuntimeError):
    """The requested figure's verified reference contract was not satisfied."""


def file_sha256(path):
    digest = sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def require_equal(name, actual, expected):
    if actual != expected:
        raise ReferenceMismatch(f'{name} differs from the verified reference')


def validate_preflight(contract, proposed):
    """Reject unverified code, changed identities/order, parameters or annotations.

    ``authorized_changes`` must contain only changes explicitly requested by the
    user. It does not exempt code: reviewed updated code needs its own pinned hash.
    Required keys are deliberately accessed directly: absent evidence is an error.
    """
    if contract['unresolved']:
        raise ReferenceMismatch('Unresolved reference provenance: ' + '; '.join(contract['unresolved']))
    if not contract['historical_code_verified']:
        raise ReferenceMismatch('Historical generating code has not been verified')
    for key in ('source_samples', 'source_features', 'patterns', 'callout_config'):
        require_equal(key, proposed[key], contract[key])
    expected = dict(contract['parameters'])
    expected.update(contract['authorized_changes'])
    require_equal('parameters (including authorized changes only)', proposed['parameters'], expected)
    require_equal('execution code paths and hashes', proposed['code'], contract['execution_code'])
    for key in ('reference_files', 'execution_code'):
        if not contract[key]:
            raise ReferenceMismatch(f'No pinned {key}')
        for path, expected_hash in contract[key].items():
            if not Path(path).is_file():
                raise ReferenceMismatch(f'Missing required {key}: {path}')
            require_equal(str(path), file_sha256(path), expected_hash)


def validate_artifacts(expected, actual):
    """Require exact staged row identities/order, pattern assignments and callouts.

    For a rerun with authorized changed scoring, expected selection must come from
    that original selection algorithm on the new scores, not the historical top-N
    list. No count-only checks or implicit intersections are permitted.
    """
    for key in ('samples', 'events', 'event_patterns', 'callout_rows', 'render_settings'):
        require_equal(key, actual[key], expected[key])
    if len(set(actual['samples'])) != len(actual['samples']):
        raise ReferenceMismatch('Duplicate output sample identity')
    if len(set(actual['events'])) != len(actual['events']):
        raise ReferenceMismatch('Duplicate output event identity')
