"""Explicit native/corrected BnpC adapters; caller-provided pinned source stays intact."""
import difflib
import json
import pickle
import shutil
import sys
from copy import deepcopy
from pathlib import Path
from types import MethodType
import numpy as np
import pandas as pd
from .integrity import require, sha256

PINNED_COMMIT = 'edb1ddb1beb1bf977d44aede4c50e6bf985e06a8'
REPLACEMENTS = (
    ('gamma_fct(*self.DP_a_gamma)',
     'gamma_fct(self.DP_a_gamma[0], loc=0, scale=1 / self.DP_a_gamma[1])'),
    ('self.DP_a_gamma[0] + k, self.DP_a_gamma[1] - np.log(eta)',
     'self.DP_a_gamma[0] + k, 1 / (self.DP_a_gamma[1] - np.log(eta))'),
    ('self.DP_a_gamma[0] + k - 1, self.DP_a_gamma[1] - np.log(eta)',
     'self.DP_a_gamma[0] + k - 1, 1 / (self.DP_a_gamma[1] - np.log(eta))'),
)


def verify_source(source):
    source = Path(source).resolve()
    manifest = json.loads(Path(__file__).with_name('bnpc_source_manifest.json').read_text())
    for relative, expected in manifest['files'].items():
        require((source / relative).is_file(), f'Missing pinned BnpC source: {relative}')
        require(sha256(source / relative) == expected, f'Unsupported/modified BnpC source: {relative}')
    return manifest


def prepare_source(source, destination, corrected):
    source, destination = Path(source).resolve(), Path(destination)
    manifest = verify_source(source)
    require(not destination.exists(), 'Refuse to overwrite source snapshot')
    shutil.copytree(source, destination, ignore=shutil.ignore_patterns('.git', '__pycache__', '*.pyc'))
    old = (source / 'libs/CRP.py').read_text()
    new = old
    if corrected:
        for a, b in REPLACEMENTS:
            require(new.count(a) == 1, 'Pinned gamma expression not unique')
            new = new.replace(a, b)
        (destination / 'libs/CRP.py').write_text(new)
    (destination / 'ALTANALYZE_PATCH.diff').write_text(''.join(difflib.unified_diff(
        old.splitlines(True), new.splitlines(True), fromfile='native/CRP.py', tofile='selected/CRP.py')))
    # All copied code except the three explicitly chosen expressions must be identical.
    for relative, expected in manifest['files'].items():
        require(sha256(source / relative) == expected, 'Native source changed during preparation')
        if relative != 'libs/CRP.py':
            require(sha256(destination / relative) == expected, f'Unexpected source change: {relative}')
    return destination


def summarize_sampling(X, result, point_genotypes, output):
    """Report all saved native draws without applying burn-in twice to params."""
    burn = int(result['burn_in'])
    assignments = result['assignments'][burn:]
    parameters = result['params']  # Native get_results already removes parameter burn-in.
    require(len(assignments) == len(parameters) and len(assignments) > 1,
            'Saved assignment/parameter states are not aligned')
    mean = np.zeros(X.shape)
    fraction = np.zeros(X.shape)
    prediction = np.zeros(X.shape)
    point_support = np.zeros(len(X))
    max_error = 0.
    for index, assignment in enumerate(assignments):
        ids, inverse = np.unique(assignment, return_inverse=True)
        gp = parameters[index, :len(ids)][inverse].astype(float)
        require(gp.shape == X.shape and np.isfinite(gp).all() and ((gp >= 0) & (gp <= 1)).all(),
                'Invalid saved genotype parameters')
        fp, fn = result['FP'][burn + index], result['FN'][burn + index]
        require(fp >= 0 and fn >= 0 and fp + fn < 1, 'Invalid saved error rates')
        p = gp * (1 - fn) + (1 - gp) * fp
        ll = np.where(X >= 0, np.where(X == 1, np.log(p), np.log1p(-p)), 0).sum()
        saved = result['ML'][burn + index]
        error = abs(ll - saved)
        require(error < max(.1, abs(saved) * 1e-6), 'Saved state fails independent likelihood check')
        max_error = max(max_error, float(error))
        mean += gp
        fraction += gp >= .5
        prediction += p
        point_support += ((gp >= .5) == (point_genotypes >= .5)).all(1)
    states = len(assignments)
    np.savez_compressed(Path(output) / 'sampling_summary.npz',
                        parameter_mean=mean / states,
                        threshold_positive_fraction=fraction / states,
                        predictive_call_probability=prediction / states,
                        MAP_signature_state_fraction=point_support / states,
                        informative=(X >= 0).any(1))
    return dict(retained_states=states, max_saved_state_likelihood_difference=max_error,
                support_semantics='Fraction of saved states supporting the thresholded MAP genotype; not calibrated cell posterior or native MPEAR')


def fit(X, output, source, corrected, seed, steps, burn_fraction, initialization,
        concentration_shape=None, concentration_rate=1.):
    require(steps >= 10 and 0 < burn_fraction < 1, 'Invalid sampling budget/burn-in')
    require(initialization in ['together', 'random'], 'Unsupported initialization')
    require(concentration_rate > 0, 'Concentration rate must be positive')
    require(concentration_shape is None or concentration_shape > 0, 'Concentration shape must be positive')
    require(not any(k == 'libs' or k.startswith('libs.') for k in sys.modules),
            'Run BnpC in a fresh Python process to prevent mixing native/corrected imports')
    prepared = prepare_source(source, Path(output) / 'bnpc_source', corrected)
    sys.path.insert(0, str(prepared))
    from libs.CRP_learning_errors import CRP_errors_learning
    from libs.MCMC import MCMC, Chain
    from libs.utils import get_latents_point
    data = np.asarray(X, dtype=float).copy()
    data[data < 0] = np.nan
    alpha = [-1, -1] if concentration_shape is None else [concentration_shape, concentration_rate]
    model = CRP_errors_learning(data, DP_alpha=alpha, param_beta=[.75, .75],
                                FP_mean=.01, FP_sd=.01, FN_mean=.2, FN_sd=.1)
    if initialization == 'together':
        def together(self, mode='random', assign=False):
            require(not assign, 'Fixed assignments are not accepted')
            CRP_errors_learning.init(self, mode='together', assign=assign)
            require(len(self.cells_per_cluster) == 1 and np.all(self.assignment == 0), 'Initialization failed')
        model.init = MethodType(together, model)
    sampler = MCMC(model, sm_prob=.33, dpa_prob=.25, error_prob=.25,
                   sm_ratios=[.75, .25], sm_steps=3)
    probe = deepcopy(model)
    require(probe.__class__ is CRP_errors_learning, 'Native class identity changed')
    if initialization == 'together':
        require(probe.init.__self__ is probe, 'Deepcopy initializer bound to wrong object')
    require(Chain(probe, sampler.params, 1).learning_errors, 'Native error learning disabled')
    sampler.run((steps, int(steps * burn_fraction)), seed, n=1, verbosity=2, debug=True)
    require(sampler.chains[0].learning_errors and not sampler.chains[0].fix_assign,
            'Native learning/assignment behavior changed')
    result = sampler.get_results()
    with (Path(output) / 'native_results.pkl').open('wb') as stream:
        pickle.dump(result, stream)
    point = get_latents_point(result, 'MAP', data, single_chains=True)[0]
    membership = np.asarray(point['assignment']).astype(str)
    gp = np.asarray(point['genotypes']).T
    require(gp.shape == X.shape, 'Native genotype export shape mismatch')
    require(np.isfinite(gp).all() and ((gp >= 0) & (gp <= 1)).all(), 'Invalid native genotype parameters')
    probability = gp * (1 - point['FN']) + (1 - gp) * point['FP']
    ll = np.where(X >= 0, np.where(X == 1, np.log(probability), np.log1p(-probability)), 0).sum()
    saved = result[0]['ML'][int(point['step'])]
    require(abs(ll - saved) < max(.1, abs(saved) * 1e-6), 'Exported MAP state fails independent likelihood check')
    ids, first, inverse = np.unique(membership, return_index=True, return_inverse=True)
    require(np.allclose(gp, gp[first][inverse]),
            'Cell/cluster genotype export mismatch')
    pd.DataFrame({key: result[0][key] for key in ['ML', 'MAP', 'FP', 'FN', 'DP_alpha']}).to_csv(
        Path(output) / 'trace.tsv', sep='\t', index=False)
    summary = {key: float(point[key]) for key in ['step', 'a', 'FP', 'FN']}
    summary.update(states=len(result[0]['ML']), burn_in=int(result[0]['burn_in']),
                   clusters=len(ids), MAP_likelihood_difference=float(ll - saved),
                   genotype_parameter_semantics='Native MAP cluster Bernoulli parameters; not cell posterior probabilities',
                   convergence='Not established by reaching requested iteration count',
                   native_alpha_floor='Native concentration clipping at >=1 retained')
    summary.update(summarize_sampling(X, result[0], gp, output))
    (Path(output) / 'bnpc_summary.json').write_text(json.dumps(summary, indent=2))
    return membership, pd.DataFrame(gp[first], index=ids), summary
