import numpy as np
import pytest
from altanalyze3.components.visualization.pseudobulk_depth_qc import expected_detected_features


def test_without_replacement_expectation_and_unavailable_target():
    c=np.array([2,1,0]);original=c.copy()
    assert np.isclose(expected_detected_features(c,2),5/3)
    assert expected_detected_features(c,0)==0
    assert expected_detected_features(c,3)==2
    assert np.isnan(expected_detected_features(c,4))
    assert np.array_equal(c,original)
    with pytest.raises(ValueError):expected_detected_features(np.array([.5,1.5]),1)
