"""Clone inference with complete input contracts and explicit uncertainty reporting."""
import argparse
import json
import os
import platform
from pathlib import Path
import time
import numpy as np
import pandas as pd
from .integrity import load_inputs, numerical_gate, require, sha256
from .reporting import export_point_assignments


def algorithm_spec(args):
    if args.method == 'detection-tree':
        from .detection_workflow import specification
        return specification(args.seed, args.sex, args.pooling_strength, args.rna_concentration)
    common = dict(method=args.method, seed=args.seed, unknown=-1, features='complete contract roster', source_precedence='already consolidated upstream; unchanged')
    if args.method == 'supported-combinations':
        common.update(grouping='exact accepted MUT set; WT/unknown retained per cell', additional_imputation=False, population_extrapolation=False, no_positive_categories=['WT', 'WTobs', 'NE'], selection='none; deterministic', annotation_only='uncertain inherited status retained, not asserted somatic')
    elif args.method == 'exhaustive':
        common.update(FP=args.false_positive_rate, FN=args.false_negative_rate, sample_weights='independent',
                      candidate_space='all 2**J configurations', convergence_gap=1e-5)
    elif args.method == 'shared-genotype':
        common.update(FP_cap=args.fp_cap, FN_cap=args.fn_cap, starts=args.starts, K_grid=args.k_grid,
                      max_iterations=args.max_iterations, validation_seed=args.validation_seed,
                      validation_fraction=.2, inner_mask_seed=args.validation_seed+1,
                      weight_Dirichlet=1.5, error_Beta=[1.5,1.5], EM_tolerance=1e-7,
                      initialization='25 soft EM steps then binary coordinate EM',
                      sample_weights='specimen-specific; shared discrete genotypes',
                      selection='minimum mean inner log loss; best full training objective across starts',
                      global_unobserved_reporting='retain event, report unknown, sum equivalent component posterior mass')
    else:
        common.update(iterations=args.iterations, burn_fraction=args.burn_fraction,
                      initialization=args.initialization, concentration_shape=args.concentration_shape,
                      concentration_rate=args.concentration_rate, param_beta=[.75,.75],
                      FP_prior=[.01,.01], FN_prior=[.2,.1], sm_prob=.33, dpa_prob=.25, error_prob=.25,
                      sm_ratios=[.75,.25], sm_steps=3, native_source_commit='edb1ddb1beb1bf977d44aede4c50e6bf985e06a8')
    return common


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix', required=True, help='NPZ: X, cell_ids, features; unknown=-1')
    p.add_argument('--cells', required=True, help='TSV: cell_id, sample; optional stage/cell_state')
    p.add_argument('--events', required=True, help='TSV: feature plus optional provenance/labels')
    p.add_argument('--contract', required=True, help='Independent complete roster, provenance and input hash contract')
    p.add_argument('--output', required=True)
    p.add_argument('--method', default='detection-tree', choices=['detection-tree', 'supported-combinations', 'shared-genotype', 'bnpc-native', 'bnpc-gamma-corrected', 'exhaustive'], help='Default: reviewed detection-aware gain-only tree; cells without accepted MUT calls are reported WT. Legacy engines remain explicit comparisons.')
    p.add_argument('--fp-cap', type=float, default=.10)
    p.add_argument('--fn-cap', type=float, default=.30)
    p.add_argument('--starts', type=int, default=3)
    p.add_argument('--k-grid', type=lambda s:[int(v) for v in s.split(',')], default=[2,4,6,8,12,16,24,32,48,64,96,128])
    p.add_argument('--max-iterations', type=int, default=10000)
    p.add_argument('--validation-seed', type=int, default=20261114)
    p.add_argument('--seed', type=int, default=None)
    p.add_argument('--detector', help='Verified NPZ: Z, available, cell_ids, features')
    p.add_argument('--sex', choices=['M','F','U'])
    p.add_argument('--pooling-strength', type=float, default=2.)
    p.add_argument('--rna-concentration', type=float, default=2.)
    p.add_argument('--bnpc-source', help='Unmodified pinned BnpC checkout; never changed in place')
    p.add_argument('--iterations', type=int, default=1000)
    p.add_argument('--burn-fraction', type=float, default=.33)
    p.add_argument('--initialization', choices=['random', 'together'], default='together')
    p.add_argument('--concentration-shape', type=float)
    p.add_argument('--concentration-rate', type=float, default=1.)
    p.add_argument('--false-positive-rate', type=float)
    p.add_argument('--false-negative-rate', type=float)
    p.add_argument('--families', help='Reporting-only JSON: family -> {exact feature ID: 0 or 1}')
    args = p.parse_args(argv)
    if args.seed is None: args.seed = 222 if args.method == 'detection-tree' else 0
    if args.method == 'detection-tree':
        require(args.bnpc_source is None and args.false_positive_rate is None and args.false_negative_rate is None and args.concentration_shape is None and args.concentration_rate == 1. and args.iterations == 1000 and args.burn_fraction == .33 and args.initialization == 'together' and args.fp_cap == .10 and args.fn_cap == .30 and args.starts == 3 and args.max_iterations == 10000 and args.validation_seed == 20261114 and args.k_grid == [2,4,6,8,12,16,24,32,48,64,96,128] and args.families is None, 'Legacy engine options cannot silently change detection-tree')
        from .detection_workflow import run
        run(args)
        return
    require(args.detector is None and args.sex is None and args.pooling_strength == 2. and args.rna_concentration == 2., 'Detection-tree options cannot silently change legacy engines')
    X, cells, events, contract = load_inputs(args.matrix, args.cells, args.events, args.contract, args.method)
    require(contract.get('algorithm') == algorithm_spec(args), 'Algorithm differs from explicit contract')
    numerical_gate()
    families = json.loads(Path(args.families).read_text()) if args.families else {}
    for name, definition in families.items():
        require(bool(definition), f'Empty family: {name}')
        require(set(definition) <= set(events.feature), f'Unresolved family features: {name}')
        require(set(definition.values()) <= {0, 1}, f'Invalid family states: {name}')
    if args.method != 'shared-genotype':
        require(args.fp_cap == .10 and args.fn_cap == .30 and args.starts == 3
                and args.k_grid == [2,4,6,8,12,16,24,32,48,64,96,128]
                and args.max_iterations == 10000 and args.validation_seed == 20261114,
                'Shared-genotype-only options cannot alter another engine; remove them')
    if args.method == 'exhaustive':
        require(args.bnpc_source is None and args.concentration_shape is None and args.concentration_rate == 1.
                and args.iterations == 1000 and args.burn_fraction == .33 and args.initialization == 'together',
                'BnpC-only options cannot alter an exhaustive fit; remove them')
        require(args.false_positive_rate is not None and args.false_negative_rate is not None,
                'Explicit FP and FN required; there are no calibrated universal defaults')
        require(0 < args.false_positive_rate < .5 and 0 < args.false_negative_rate < .5, 'FP/FN outside supported range')
        require(len(events) <= 16, 'Exhaustive enumeration exceeds 16 events; no automatic feature reduction')
    elif args.method == 'shared-genotype':
        require(args.bnpc_source is None and args.false_positive_rate is None and args.false_negative_rate is None
                and args.concentration_shape is None and args.concentration_rate == 1.
                and args.iterations == 1000 and args.burn_fraction == .33 and args.initialization == 'together',
                'BnpC/exhaustive-only options cannot change shared-genotype; use explicit error caps')
        require(.001 <= args.fp_cap < .5 and .001 <= args.fn_cap < .5, 'Invalid learned-error caps')
        require(args.starts >= 2 and args.max_iterations > 25 and args.seed >= 0, 'Invalid starts, iterations or seed')
        require(bool(args.k_grid) and all(k >= 1 for k in args.k_grid), 'Invalid K grid')
    elif args.method == 'supported-combinations':
        require(args.bnpc_source is None and args.false_positive_rate is None and args.false_negative_rate is None
                and args.concentration_shape is None and args.concentration_rate == 1.
                and args.iterations == 1000 and args.burn_fraction == .33 and args.initialization == 'together',
                'Sampler/error options do not apply to supported combinations')
    else:
        require(args.false_positive_rate is None and args.false_negative_rate is None,
                'Fixed FP/FN options apply only to exhaustive; this BnpC adapter learns error rates')
        require(args.bnpc_source is not None, '--bnpc-source required')
        require(args.seed > 0, 'BnpC seed must be positive for reproducibility')
        require(args.concentration_shape is not None or args.concentration_rate == 1,
                'A nondefault concentration rate requires an explicit shape')
        from .bnpc import verify_source
        verify_source(args.bnpc_source)
    output = Path(args.output).resolve()
    require(not output.exists(), 'Refuse to overwrite an existing analysis directory')
    output.mkdir(parents=True)
    started = time.time()
    manifest = dict(arguments=vars(args), contract=contract, python=platform.python_version(),
                    numpy=np.__version__, OPENBLAS_CORETYPE=os.environ.get('OPENBLAS_CORETYPE'),
                    cells=len(cells), events=len(events), samples=cells['sample'].nunique(),
                    code_sha256={f.name: sha256(f) for f in Path(__file__).parent.glob('*.py')},
                    interpretation='Experimental inference; reaching an iteration limit is not convergence or biological validation')
    (output / 'run_manifest.json').write_text(json.dumps(manifest, indent=2))
    cells.to_csv(output / 'cell_roster.tsv', sep='\t', index=False)
    events.to_csv(output / 'event_roster.tsv', sep='\t', index=False)
    if args.method.startswith('bnpc'):
        from .bnpc import fit
        membership, G, summary = fit(X, output, args.bnpc_source, args.method == 'bnpc-gamma-corrected',
                                     args.seed, args.iterations, args.burn_fraction, args.initialization,
                                     args.concentration_shape, args.concentration_rate)
        with np.load(output / 'sampling_summary.npz', allow_pickle=False) as saved:
            sampling = {key: saved[key] for key in saved.files}
        np.savez_compressed(output / 'sampling_summary.npz', **sampling,
                            cell_ids=cells.cell_id.to_numpy(dtype=str), features=events.feature.to_numpy(dtype=str))
        assignments = export_point_assignments(X, cells, events, membership, G, output)
        G.columns = events.feature
        rows = []
        for sample, ix in cells.groupby('sample', sort=False).groups.items():
            ix = np.asarray(ix)
            informative = (X[ix] >= 0).any(1)
            gp = G.loc[membership[ix]].to_numpy() >= .5
            for name, definition in families.items():
                cols = [list(events.feature).index(f) for f in definition]
                count = int(((gp[:, cols] == np.array(list(definition.values()))).all(1) & informative).sum())
                rows.append(dict(sample=sample, family=name, cells=count, all_cells=len(ix), percent_all=100*count/len(ix),
                                 uncertainty='MAP point assignment; not a posterior probability'))
        if rows:
            pd.DataFrame(rows).to_csv(output / 'family_frequencies.tsv', sep='\t', index=False)
    elif args.method == 'supported-combinations':
        from .supported import run_supported
        run_supported(X, cells, events, output)
    elif args.method == 'shared-genotype':
        from .shared_workflow import run_shared
        run_shared(X, cells, events, output, args, families)
    else:
        from .exhaustive import fit, equivalent_bounds
        summaries, frequencies, rows = [], [], []
        assignment = cells.copy()
        assignment['native_config'] = -1
        assignment['max_configuration_probability'] = np.nan
        assignment['observed_events'] = (X >= 0).sum(1)
        event_p = np.zeros(X.shape)
        for number, (sample, ix) in enumerate(cells.groupby('sample', sort=False).groups.items()):
            ix = np.asarray(ix)
            f = fit(X[ix], args.false_positive_rate, args.false_negative_rate, args.seed)
            G, P, inv, w = f['G'], f['P'], f['inverse'], f['w']
            sample_name = f'sample_{number:03d}'  # User sample names never become paths.
            np.savez_compressed(output / (sample_name + '.npz'), **{k: f[k] for k in ['G','patterns','inverse','counts','w','P','offset']},
                                cell_ids=cells.cell_id.to_numpy(dtype=str)[ix], features=events.feature.to_numpy(dtype=str))
            assignment.loc[ix, 'native_config'] = P.argmax(1)[inv]
            assignment.loc[ix, 'max_configuration_probability'] = P.max(1)[inv]
            gp = P @ G
            reference = np.column_stack([(P * G[:, j]).sum(1) for j in range(len(events))])
            require(np.allclose(gp, reference, atol=1e-10, rtol=0), 'Posterior numerical gate failed')
            event_p[ix] = gp[inv]
            known = (f['patterns'] >= 0).any(1)
            counts = f['counts'] * known
            soft_counts = counts @ P
            for k in range(len(w)):
                frequencies.append(dict(sample=sample, config=k, genotype=''.join(map(str,G[k])),
                                        population_weight=w[k], expected_informative_cells=soft_counts[k],
                                        all_cells=len(ix), expected_percent_all=100*soft_counts[k]/len(ix)))
            # Family definitions only enter after fitting and never constrain configurations.
            for name, definition in families.items():
                cols = [list(events.feature).index(e) for e in definition]
                mask = (G[:, cols] == np.array(list(definition.values()))).all(1)
                endpoints = equivalent_bounds(f, mask)
                n = float(soft_counts[mask].sum())
                rows.append(dict(sample=sample, family=name, population_fraction=float(w[mask].sum()),
                                 likelihood_equivalent_min=endpoints[0]['value'], likelihood_equivalent_max=endpoints[1]['value'],
                                 expected_informative_cells=n, all_cells=len(ix), expected_percent_all=100*n/len(ix),
                                 uncertainty='Likelihood-equivalence range under fixed error model; not a confidence interval'))
            summaries.append(dict(sample=sample, artifact=sample_name+'.npz', loglik=f['loglik'], dual_gap=f['gap'],
                                  iterations=f['iterations'], polish_iterations=f['polish_iterations'],
                                  no_evidence=int((X[ix] < 0).all(1).sum()), cells=len(ix)))
        assignment['display_config'] = np.where(assignment.observed_events > 0, assignment.native_config.astype(str), 'Unresolved_no_evidence')
        assignment.to_csv(output / 'cell_assignments.tsv.gz', sep='\t', index=False)
        np.savez_compressed(output / 'cell_event_posteriors.npz', genotype_probability=event_p,
                            cell_ids=cells.cell_id.to_numpy(dtype=str), features=events.feature.to_numpy(dtype=str))
        pd.DataFrame(summaries).to_csv(output / 'fit_summary.tsv', sep='\t', index=False)
        pd.DataFrame(frequencies).to_csv(output / 'configuration_frequencies.tsv', sep='\t', index=False)
        if rows:
            pd.DataFrame(rows).to_csv(output / 'family_frequencies.tsv', sep='\t', index=False)
    (output / 'COMPLETE.json').write_text(json.dumps(dict(cells=len(cells), events=len(events), elapsed_seconds=time.time()-started,
                                                        scope='Computation completed; scientific validation must be assessed separately'), indent=2))


if __name__ == '__main__':
    main()
