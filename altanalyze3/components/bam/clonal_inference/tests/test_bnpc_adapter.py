"""Optional pinned-source integration test, using a tiny synthetic dataset."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
from altanalyze3.components.bam.clonal_inference.bnpc import summarize_sampling


class SavedStateTest(unittest.TestCase):
    def test_postburn_alignment_and_uncertainty(self):
        X=np.array([[1,0],[-1,1]])
        # Parameters contain two retained states; assignments also contain the initial state.
        result=dict(burn_in=1,assignments=np.zeros((3,2),dtype=int),
                    params=np.array([[[.8,.2]],[[.2,.8]]]),FP=np.full(3,.1),FN=np.full(3,.3),
                    ML=np.array([0.,np.log(.58*.78*.22),np.log(.22*.42*.58)]))
        with tempfile.TemporaryDirectory() as d:
            summary=summarize_sampling(X,result,np.array([[.8,.2],[.8,.2]]),d)
            q=np.load(Path(d)/'sampling_summary.npz')
            np.testing.assert_allclose(q['parameter_mean'],.5)
            np.testing.assert_allclose(q['threshold_positive_fraction'],.5)
            np.testing.assert_allclose(q['predictive_call_probability'],.4)
            np.testing.assert_allclose(q['MAP_signature_state_fraction'],.5)
            self.assertEqual(summary['retained_states'],2)
            bad=dict(result,params=np.repeat(result['params'][:1],3,axis=0))
            with self.assertRaisesRegex(ValueError,'not aligned'):summarize_sampling(X,bad,np.zeros(X.shape),d)


@unittest.skipUnless(os.environ.get('BNPC_SOURCE_PATH'), 'Set BNPC_SOURCE_PATH to the audited native checkout')
class NativeAdapterTest(unittest.TestCase):
    def test_native_class_sampling_and_source_preservation(self):
        script = r'''
import json,sys
from pathlib import Path
import numpy as np
from altanalyze3.components.bam.clonal_inference.bnpc import fit,verify_source
source=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir()
before=verify_source(source)
X=np.array([[0,0,-1],[0,0,0],[1,1,0],[1,1,1],[1,-1,1],[0,0,0]]*3)
labels,G,summary=fit(X,out,source,True,5801,30,.33,'together')
assert len(labels)==18 and G.shape[1]==3
assert summary['states']>=30 and abs(summary['MAP_likelihood_difference'])<.1
assert summary['retained_states']==summary['states']-summary['burn_in']
q=np.load(out/'sampling_summary.npz');assert q['parameter_mean'].shape==X.shape
assert np.all((q['MAP_signature_state_fraction']>=0)&(q['MAP_signature_state_fraction']<=1))
assert verify_source(source)==before
patch=(out/'bnpc_source/ALTANALYZE_PATCH.diff').read_text()
assert 'scale=1 / self.DP_a_gamma[1]' in patch
assert (out/'bnpc_source/LICENSE').read_text()==(source/'LICENSE').read_text()
assert (out/'native_results.pkl').is_file()
print('PASS native class identity, error-learning dispatch, sampling, MAP export and preserved source')
'''
        with tempfile.TemporaryDirectory() as d:
            r=subprocess.run([sys.executable,'-c',script,os.environ['BNPC_SOURCE_PATH'],str(Path(d)/'fit')],text=True,capture_output=True)
            self.assertEqual(r.returncode,0,r.stdout+'\n'+r.stderr)


if __name__=='__main__':unittest.main()
