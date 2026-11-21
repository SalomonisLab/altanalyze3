import unittest
import numpy as np
from altanalyze3.components.bam.clonal_inference.wt_evidence_guard import protect_wt_evidence


class WTEvidenceGuardTests(unittest.TestCase):
    def test_overwhelming_prior_cannot_override_wt_only(self):
        X = np.array([[0,-1,-1],[-1,-1,0],[-1,-1,-1],[1,0,0],[0,1,-1]])
        P = np.tile([1e-12, 1-1e-12], (5,1))
        original_X, original_P = X.copy(), P.copy()
        Q, groups = protect_wt_evidence(X,P,[True,True,False])
        np.testing.assert_array_equal(Q[:3,2:],np.eye(3))
        np.testing.assert_array_equal(Q[3:,:2],P[3:])
        np.testing.assert_allclose(Q.sum(1),1)
        np.testing.assert_array_equal(X, original_X)
        np.testing.assert_array_equal(P, original_P)
        self.assertEqual(list(groups),['WT_variant_evidence','LOH_negative_only',
                                     'no_call','accepted_mutation','accepted_mutation'])

    def test_no_loh_feature_and_all_unknown(self):
        Q, group = protect_wt_evidence(np.array([[0,-1],[-1,-1]]),
                                      np.array([[.3,.7],[.9,.1]]),[True,True])
        np.testing.assert_array_equal(Q[:,2:],[[1,0,0],[0,0,1]])

    def test_reject_invalid_inputs(self):
        for P in [np.array([[.5,np.nan]]),np.array([[.8,.8]]),np.array([[-.1,1.1]])]:
            with self.assertRaises(ValueError):
                protect_wt_evidence(np.array([[0,-1]]),P,[True,False])
        with self.assertRaises(ValueError):
            protect_wt_evidence(np.array([[0,2]]),np.array([[.3,.7]]),[True,False])
        with self.assertRaises(ValueError):
            protect_wt_evidence(np.array([[0,-1]]),np.array([[.3,.7]]),[True])
