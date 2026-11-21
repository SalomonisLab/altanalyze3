import unittest
import numpy as np
from altanalyze3.components.bam.clonal_inference.exhaustive import configurations,likelihood,fit,equivalent_bounds


class LikelihoodTests(unittest.TestCase):
    def test_known_mixture_recovery_and_unknown_invariance(self):
        G=configurations(2);L,o=likelihood(G,G,.1,.2);w=np.array([.4,.1,.1,.4])
        counts=np.rint((L*np.exp(o[:,None])@w)*10000).astype(int)
        X=np.repeat(G,counts,axis=0);a=fit(X,.1,.2);b=fit(np.vstack([X,np.full((20,2),-1)]),.1,.2)
        np.testing.assert_allclose(a['w'],w,atol=1e-5);np.testing.assert_allclose(a['w'],b['w'],atol=1e-5)
        self.assertAlmostEqual(a['loglik'],b['loglik'],places=5)

    def test_missing_linkage_preserves_wide_joint_bounds(self):
        X=np.array([[0,-1],[1,-1],[-1,0],[-1,1]]*100)
        f=fit(X,.1,.1);b=equivalent_bounds(f,(f['G']==1).all(1))
        self.assertLess(b[0]['value'],1e-5);self.assertAlmostEqual(b[1]['value'],.5,places=5)

    def test_full_panel_likelihood_against_product(self):
        G=configurations(11);x=G[:337].copy();x[::3,::2]=-1
        L,o=likelihood(x,G,.01,.3);p=.01+.69*G
        expected=np.ones(L.shape)
        for j in range(11):expected*=np.where(x[:,j,None]==1,p[:,j],np.where(x[:,j,None]==0,1-p[:,j],1.))
        np.testing.assert_allclose(L*np.exp(o[:,None]),expected,atol=1e-12)

    def test_no_implicit_feature_reduction(self):
        with self.assertRaisesRegex(ValueError,'no implicit feature reduction'):configurations(17)


if __name__=='__main__':unittest.main()
