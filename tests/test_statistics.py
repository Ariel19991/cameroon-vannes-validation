import numpy as np

from vannes_validation.statistics import compare_to_reference


def test_comparison_statistics_have_known_values():
    result = compare_to_reference([1.0, 2.0, 3.0], [0.0, 2.0, 4.0], "candidate")
    assert result["bias_mgal"] == 0.0
    np.testing.assert_allclose(result["mae_mgal"], 2.0 / 3.0)
    np.testing.assert_allclose(result["rmse_mgal"], np.sqrt(2.0 / 3.0))
    assert result["negative_count"] == 0


def test_comparison_statistics_reject_shape_mismatch():
    with np.testing.assert_raises(ValueError):
        compare_to_reference([1.0], [1.0, 2.0], "invalid")
