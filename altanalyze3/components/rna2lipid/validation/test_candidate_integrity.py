"""Executable failure-mode tests; no real model fit or user approval is made."""
import copy
import importlib
import inspect
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import candidate_integrity as gate


class IntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'source.txt'
        self.source.write_text('Synthetic unit-test evidence, not real experimental measurements.')
        self.record = {'path': str(self.source), 'sha256': gate.sha256(self.source)}
        self.contract = {
            'schema_version': 1, 'samples': ['D018_END', 'D071_END'],
            'lipids': ['L1', 'L2', 'L3'], 'genes': ['G1', 'G2'],
            'production_bundle': self.record, 'RNA_source': self.record,
            'pinned_code': {k: self.record for k in ['trainer', 'api', 'delivered_trainer']},
            'method': {'estimator': 'ElasticNetCV', 'cv': 3, 'alpha_grid': [0.01, 1],
                       'preprocessing': 'baseline', 'target_scale': 'log2', 'seed': 1,
                       'differential_policy': 'full-panel BH'},
        }
        self.x = pd.DataFrame([[1., 2.], [3., 4.]], index=self.contract['samples'], columns=self.contract['genes'])
        self.y = pd.DataFrame([[1., 2., 3.], [4., 5., 6.]], index=self.contract['samples'], columns=self.contract['lipids'])
        self.provenance = {'schema_version': 1, 'data_role': 'corrected_training_targets', 'unresolved': [],
            'transformation_description': 'Synthetic test only', 'imputed_entries': [],
            'target_sha256': self.record['sha256']}
        for key, ids in [('sample_sources', self.contract['samples']), ('feature_sources', self.contract['lipids'])]:
            self.provenance[key] = [{'id': i, 'status': 'verified', 'evidence': 'Synthetic test only',
                'source_identifier': i, 'source': self.record} for i in ids]

    def test_exact_coverage_passes(self):
        gate.validate_contract(self.contract)
        gate.validate_tables(self.x, self.y, self.contract)
        gate.validate_provenance(self.provenance, self.contract)

    def test_full_metacell_roster_check_finishes_without_quadratic_work(self):
        ids = [f'metacell_{i}' for i in range(230057)]
        started = time.monotonic()
        gate.exact_ids(ids, list(ids), 'Full metacell roster')
        self.assertLess(time.monotonic() - started, 5.0)

    def test_reduced_lipid_panel_rejected(self):
        with self.assertRaisesRegex(gate.IntegrityError, 'L3'):
            gate.validate_tables(self.x, self.y.iloc[:, :2], self.contract)

    def test_missing_D071_rejected(self):
        with self.assertRaisesRegex(gate.IntegrityError, 'D071_END'):
            gate.validate_tables(self.x.iloc[:1], self.y.iloc[:1], self.contract)

    def test_same_count_wrong_identity_rejected(self):
        with self.assertRaisesRegex(gate.IntegrityError, 'missing=.*L3'):
            gate.validate_tables(self.x, self.y.rename(columns={'L3': 'OTHER'}), self.contract)

    def test_same_count_duplicate_identity_rejected(self):
        y = self.y.copy(); y.columns = ['L1', 'L2', 'L2']
        with self.assertRaisesRegex(gate.IntegrityError, 'duplicates'):
            gate.validate_tables(self.x, y, self.contract)

    def test_row_order_mismatch_rejected(self):
        with self.assertRaisesRegex(gate.IntegrityError, 'order'):
            gate.validate_tables(self.x, self.y.iloc[::-1], self.contract)

    def test_input_gene_loss_rejected(self):
        with self.assertRaisesRegex(gate.IntegrityError, 'G2'):
            gate.validate_tables(self.x.iloc[:, :1], self.y, self.contract)

    def test_nonfinite_cannot_be_silently_filled(self):
        for value in [float('nan'), float('inf'), 'unmatched']:
            with self.subTest(value=value):
                y = self.y.astype(object); y.iloc[0, 0] = value
                with self.assertRaisesRegex(gate.IntegrityError, 'Nonfinite'):
                    gate.validate_tables(self.x, y, self.contract)

    def test_algorithm_and_parameter_changes_rejected(self):
        changes = {'estimator': 'KernelRidge', 'cv': 5, 'alpha_grid': [.1],
                   'preprocessing': 'study_alignment', 'target_scale': 'log10',
                   'seed': 8, 'differential_policy': '47-panel BH'}
        for key, value in changes.items():
            with self.subTest(key=key):
                proposed = dict(self.contract['method'], **{key: value})
                with self.assertRaisesRegex(gate.IntegrityError, key):
                    gate.require_method(self.contract['method'], proposed)

    def test_missing_provenance_record_rejected(self):
        self.provenance['sample_sources'].pop()
        with self.assertRaisesRegex(gate.IntegrityError, 'D071_END'):
            gate.validate_provenance(self.provenance, self.contract)

    def test_unresolved_provenance_rejected(self):
        self.provenance['unresolved'] = ['D071 upstream source']
        with self.assertRaisesRegex(gate.IntegrityError, 'unresolved'):
            gate.validate_provenance(self.provenance, self.contract)

    def test_reconstructed_legacy_not_accepted_as_corrected(self):
        self.provenance['data_role'] = 'legacy_reconstruction'
        with self.assertRaisesRegex(gate.IntegrityError, 'legacy reconstruction'):
            gate.validate_provenance(self.provenance, self.contract)

    def test_unverified_mapping_rejected(self):
        self.provenance['feature_sources'][0]['status'] = 'hypothesis'
        with self.assertRaisesRegex(gate.IntegrityError, 'unresolved source'):
            gate.validate_provenance(self.provenance, self.contract)

    def test_changed_source_rejected(self):
        self.source.write_text('Changed after review')
        with self.assertRaisesRegex(gate.IntegrityError, 'hash changed'):
            gate.validate_contract(self.contract)

    def test_missing_source_requires_question(self):
        self.source.unlink()
        with self.assertRaisesRegex(gate.IntegrityError, 'Ask the user where it is'):
            gate.validate_provenance(self.provenance, self.contract)

    def test_new_target_imputation_requires_decision(self):
        self.provenance['imputed_entries'] = [{'sample': 'D071_END', 'lipid': 'L1'}]
        with self.assertRaisesRegex(gate.IntegrityError, 'imputation requires'):
            gate.validate_provenance(self.provenance, self.contract)

    def synthetic_review(self):
        return {'schema_version': 1, 'analysis_status': 'approved', 'questions': [],
            'approval': {'manifest_sha256': 'fixture', 'user_message': 'SYNTHETIC TEST FIXTURE ONLY',
                         'conversation_reference': 'unit-test-not-a-user-approval', 'granted_by': 'user'}}

    def test_pending_question_blocks_even_with_approval_label(self):
        state = self.synthetic_review()
        state['questions'] = [{'status': 'pending', 'question': 'Where is D071?'}]
        with self.assertRaisesRegex(gate.IntegrityError, 'Where is D071'):
            gate.require_review(state, 'fixture')

    def test_silence_or_missing_approval_never_passes(self):
        for approval in [None, {}, {'granted_by': 'user'}]:
            state = self.synthetic_review(); state['approval'] = approval
            with self.assertRaises(gate.IntegrityError):
                gate.require_review(state, 'fixture')

    def test_agent_cannot_self_authorize(self):
        state = self.synthetic_review(); state['approval']['granted_by'] = 'assistant'
        with self.assertRaisesRegex(gate.IntegrityError, 'Only the user'):
            gate.require_review(state, 'fixture')

    def test_approval_not_transferable_to_another_run(self):
        with self.assertRaisesRegex(gate.IntegrityError, 'exact run manifest'):
            gate.require_review(self.synthetic_review(), 'different-manifest')

    def test_resolved_label_without_answer_rejected(self):
        state = self.synthetic_review(); state['questions'] = [{'status': 'resolved'}]
        with self.assertRaisesRegex(gate.IntegrityError, 'lacks its user answer'):
            gate.require_review(state, 'fixture')

    def test_approved_D071_scope_retains_all_features_and_other_samples(self):
        c = gate.read_json(gate.CONTRACT)
        s = gate.read_json(gate.STATE)
        result = gate.candidate_contract(c, {'approved_scope_change': s['approved_scope_change']}, s)
        self.assertEqual(len(result['samples']), 45)
        self.assertEqual(result['lipids'], c['lipids'])
        self.assertEqual(result['genes'], c['genes'])
        self.assertEqual(result['samples'], [k for k in c['samples'] if not k.startswith('D071_')])
        self.assertEqual(len(c['samples']), 50)

    def test_scope_cannot_add_another_exclusion_or_change_panel(self):
        c = gate.read_json(gate.CONTRACT)
        s = gate.read_json(gate.STATE)
        for key, value in [('excluded_profiles', s['approved_scope_change']['excluded_profiles'] + ['D018_MIC']),
                           ('required_lipid_outputs', 47), ('required_RNA_inputs', 1302)]:
            altered = copy.deepcopy(s['approved_scope_change']); altered[key] = value
            with self.subTest(key=key), self.assertRaises(gate.IntegrityError):
                gate.candidate_contract(c, {'approved_scope_change': altered}, s)

    def test_phase_approval_cannot_authorize_pending_abundance_calibration(self):
        s = self.synthetic_review()
        s['approval']['analysis_phase'] = 'native_fold_validation'
        s['questions'] = [{'status': 'pending', 'question': 'Missing MS1 identification',
                           'required_for': ['MS1_abundance']}]
        gate.require_review(s, 'fixture', 'native_fold_validation')
        with self.assertRaisesRegex(gate.IntegrityError, 'Missing MS1 identification'):
            gate.require_review(s, 'fixture', 'MS1_abundance')
        with self.assertRaisesRegex(gate.IntegrityError, 'does not cover'):
            gate.require_review(s, 'fixture', 'unapproved_other_analysis')

    def test_missing_state_fails_closed(self):
        with patch.object(gate, 'STATE', self.root/'missing-state.json'):
            with self.assertRaisesRegex(gate.IntegrityError, 'Missing or unreadable'):
                gate.require_diagnostic_authorization('full_panel_source_recovery')

    def test_narrowed_contract_cannot_redefine_expected_panel(self):
        contract = self.root/'contract.json'; contract.write_text(json.dumps(self.contract))
        state = self.root/'state.json'; state.write_text(json.dumps(self.synthetic_review()))
        with patch.object(gate, 'CONTRACT', contract), patch.object(gate, 'STATE', state):
            with self.assertRaisesRegex(gate.IntegrityError, 'Baseline contract changed'):
                gate.validate_run_manifest(self.root/'irrelevant-manifest.json')

    def test_gateway_checks_actual_tables_before_estimator(self):
        import sklearn
        sys.path.insert(0, str(HERE.parents[2]))
        from altanalyze3.components.rna2lipid import training
        contract = copy.deepcopy(self.contract)
        contract['method']['sklearn_version'] = sklearn.__version__
        xp, yp = self.root/'rna.csv', self.root/'targets.csv'
        self.x.T.to_csv(xp)
        self.y.iloc[:, :2].to_csv(yp)
        with patch.object(gate, 'validate_run_manifest', return_value=(contract, {'RNA':xp, 'corrected_targets':yp})), \
                patch.object(training, 'fit_sparse_lipidwise_elasticnet') as fit:
            with self.assertRaisesRegex(gate.IntegrityError, 'Lipid output features mismatch'):
                gate.fit_verified_candidate(self.root/'synthetic-manifest.json')
            fit.assert_not_called()

    def test_complete_manifest_validation_without_training(self):
        contract_path = self.root/'contract.json'
        contract_path.write_text(json.dumps(self.contract))
        provenance_path = self.root/'provenance.json'; provenance_path.write_text(json.dumps(self.provenance))
        manifest = {'schema_version': 1, 'baseline_contract_sha256': gate.sha256(contract_path),
            'method': self.contract['method'], 'inputs': {'RNA': self.record,
            'corrected_targets': dict(self.record, data_role='corrected_training_targets'),
            'provenance': {'path': str(provenance_path), 'sha256': gate.sha256(provenance_path)}}}
        manifest_path = self.root/'manifest.json'; manifest_path.write_text(json.dumps(manifest))
        state = self.synthetic_review(); state['approval']['manifest_sha256'] = gate.sha256(manifest_path)
        state_path = self.root/'state.json'; state_path.write_text(json.dumps(state))
        with patch.object(gate, 'CONTRACT', contract_path), patch.object(gate, 'STATE', state_path), \
                patch.object(gate, 'BASELINE_CONTRACT_SHA256', gate.sha256(contract_path)):
            c, paths = gate.validate_run_manifest(manifest_path)
            self.assertEqual(c['lipids'], ['L1','L2','L3'])
            self.assertEqual(paths['RNA'], self.source)


class RealWorkflowGuards(unittest.TestCase):
    def test_current_production_contract_and_paused_state(self):
        c = gate.read_json(gate.CONTRACT)
        gate.validate_contract(c)
        self.assertEqual((len(c['samples']),len(c['lipids']),len(c['genes'])), (50,202,1303))
        with self.assertRaisesRegex(gate.IntegrityError, 'ANALYSIS BLOCKED'):
            gate.require_review(gate.read_json(gate.STATE), '')

    def test_actual_reduced_candidate_is_rejected(self):
        root = HERE/'artifacts/LungMAP_ElasticNet_MS1_candidate_20261005/elasticnet_MS1_47'
        x = pd.read_csv(root/'candidate_training_RNA.csv',index_col=0)
        y = pd.read_csv(root/'candidate_training_lipids_log2.csv',index_col=0)
        with self.assertRaisesRegex(gate.IntegrityError, 'D071'):
            gate.validate_tables(x,y,gate.read_json(gate.CONTRACT))

    def test_reconstruction_diagnostic_stays_paused(self):
        with self.assertRaisesRegex(gate.IntegrityError, 'DIAGNOSTIC BLOCKED'):
            gate.require_diagnostic_authorization('full_panel_source_recovery')

    def test_candidate_gateway_blocks_before_loading_tables(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory)/'manifest.json'; manifest.write_text('{}')
            with patch('pandas.read_csv', side_effect=AssertionError('Loaded scientific data before approval')):
                with self.assertRaisesRegex(gate.IntegrityError, 'ANALYSIS BLOCKED'):
                    gate.fit_verified_candidate(manifest)

    def test_all_retired_CLIs_block_without_writes_or_bypass_flags(self):
        drivers = gate.read_json(HERE/'integrity/protected_entrypoints.json')
        with tempfile.TemporaryDirectory() as directory:
            marker = Path(directory)/'sentinel'; marker.write_text('preserve')
            for name in drivers:
                with self.subTest(driver=name):
                    result = subprocess.run([sys.executable, str(HERE/name), '--out', directory,
                        '--force', '--skip-integrity'], cwd=directory, capture_output=True, text=True)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn('RETIRED WORKFLOW BLOCKED',result.stderr)
                    self.assertEqual(list(Path(directory).iterdir()), [marker])
                    self.assertEqual(marker.read_text(), 'preserve')

    def test_imported_orchestrators_also_block(self):
        drivers = gate.read_json(HERE/'integrity/protected_entrypoints.json')
        for name, functions in drivers.items():
            module = importlib.import_module(Path(name).stem)
            for name in functions:
                with self.subTest(module=module.__name__,function=name):
                    fn = getattr(module,name)
                    kwargs = {k:None for k,p in inspect.signature(fn).parameters.items()
                              if p.default is inspect.Parameter.empty}
                    with self.assertRaisesRegex(gate.IntegrityError,'RETIRED WORKFLOW BLOCKED'):
                        fn(**kwargs)


if __name__ == '__main__':
    unittest.main(verbosity=2)
