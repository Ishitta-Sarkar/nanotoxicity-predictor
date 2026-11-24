import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import LeaveOneOut, cross_val_predict

from nanotoxicity.model import leave_one_out, leave_one_out_mean


def test_leave_one_out_mean_holds_out_the_row():
    y = np.array([1.0, 2.0, 3.0])
    baseline = leave_one_out_mean(y)
    np.testing.assert_allclose(baseline, [2.5, 2.0, 1.5])


def test_perfect_line_has_no_leave_one_out_error():
    x = np.array([[0.0], [1.0], [2.0], [3.0], [4.0]])
    y = 3.0 * x[:, 0] - 1.0
    fit = leave_one_out(x, y)
    np.testing.assert_allclose(fit.y_pred, y, atol=1e-8)
    assert fit.model.rmse < 1e-8
    assert fit.baseline.rmse > fit.model.rmse


def test_predictions_match_unscaled_least_squares():
    rng = np.random.default_rng(0)
    x = rng.normal(size=(12, 2)) * np.array([1000.0, 0.01])
    y = 0.4 + x @ np.array([-0.002, 3.0]) + rng.normal(scale=0.05, size=12)
    fit = leave_one_out(x, y)
    unscaled = cross_val_predict(LinearRegression(), x, y, cv=LeaveOneOut())
    np.testing.assert_allclose(fit.y_pred, unscaled, atol=1e-6)
