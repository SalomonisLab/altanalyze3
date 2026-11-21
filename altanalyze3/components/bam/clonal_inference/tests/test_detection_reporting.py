import unittest
import numpy as np
import pandas as pd
from altanalyze3.components.bam.clonal_inference.detection_reporting import memberships, tables, POLICY


class DefaultReportingTests(unittest.TestCase):
    def test_no_mutation_includes_unknown_and_loh_negative(self):
        X=np.array([[-1,-1],[-1,0],[0,0],[1,0],[0,1]])
        P=np.array([[.01,.49,.5],[.01,.49,.5],[.01,.49,.5],[.8,.05,.15],[.4,.4,.2]])
        active=np.array([0,1,2]);G=np.array([[0,0],[1,0],[1,1]])
        oldX,oldP=X.copy(),P.copy()
        Q=memberships(X,P,active,G)
        np.testing.assert_array_equal(Q[:3],[[1,0,0]]*3)
        np.testing.assert_allclose(Q[3:],[[0,.25,.75],[0,2/3,1/3]])
        np.testing.assert_array_equal(X,oldX);np.testing.assert_array_equal(P,oldP)
        # Node order is an identity, not an assumed matrix position.
        np.testing.assert_allclose(memberships(X,P[:,[2,0,1]],active[[2,0,1]],G),Q[:,[2,0,1]])

    def test_zero_malignant_support_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'lacks malignant'):
            memberships(np.array([[1]]),np.array([[1.,0.]]),np.array([0,1]),np.array([[0],[1]]))

    def test_tables_keep_all_cells_and_specimens(self):
        X=np.array([[-1],[0],[1],[1]])
        c=pd.DataFrame({'cell_id':['a','b','c','d'],'sample':['s','s','s','graft']})
        e=pd.DataFrame({'feature':['f'],'event_label':['M']})
        m=dict(cell_ids=np.array(c.cell_id),features=np.array(e.feature),G=np.array([[0],[1]]),
               active=np.array([0,1]),posterior=np.array([[.01,.99]]*4))
        Q,A,F,audit=tables(X,c,e,m)
        self.assertEqual(list(A.node),[0,0,1,1]);self.assertEqual(set(A.reporting_policy),{POLICY})
        self.assertAlmostEqual(F.query('sample=="s" and node==0').reported_percent_all.iloc[0],200/3)
        self.assertAlmostEqual(F.query('sample=="graft" and node==1').reported_percent_all.iloc[0],100)
        self.assertTrue(A.all_events_unknown.iloc[0]);self.assertEqual(len(A),4)
        with self.assertRaisesRegex(ValueError,'cell identity'):
            tables(X,c.iloc[::-1],e,m)
