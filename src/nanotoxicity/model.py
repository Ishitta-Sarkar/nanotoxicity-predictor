"""Leave-one-out linear regression and a mean baseline.

Standardization is fit inside each training fold. For ordinary least squares
it does not change the predictions; it puts the full-sample coefficients in
standard-deviation units.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


@dataclass(frozen=True)
class ErrorMetrics:
    rmse: float
    mae: float
    r2: float


@dataclass(frozen=True)
class LoocvResult:
    y_pred: np.ndarray
    y_baseline: np.ndarray
    model: ErrorMetrics
    baseline: ErrorMetrics
    intercept: float
    coefficients: tuple[float, ...]


def leave_one_out(x: np.ndarray, y: np.ndarray) -> LoocvResult:
    """Predict each row from a linear model fit on the other rows."""

    if len(y) < 3:
        raise ValueError("Leave-one-out regression needs at least three rows.")
    if x.shape[0] != len(y):
        raise ValueError("x and y have different numbers of rows.")
    if x.ndim != 2 or x.shape[1] < 1:
        raise ValueError("x must have at least one descriptor column.")

    pipeline = Pipeline(
        [
            ("scale", StandardScaler()),
            ("model", LinearRegression()),
        ]
    )
    predicted = cross_val_predict(pipeline, x, y, cv=LeaveOneOut())
    baseline = leave_one_out_mean(y)
    pipeline.fit(x, y)
    coefficients = tuple(float(value) for value in pipeline.named_steps["model"].coef_)
    intercept = float(pipeline.named_steps["model"].intercept_)
    return LoocvResult(
        y_pred=np.asarray(predicted, dtype=float),
        y_baseline=baseline,
        model=_metrics(y, predicted),
        baseline=_metrics(y, baseline),
        intercept=intercept,
        coefficients=coefficients,
    )


def leave_one_out_mean(y: np.ndarray) -> np.ndarray:
    """For each row, the mean of the endpoint on the other rows."""

    y = np.asarray(y, dtype=float)
    total = float(y.sum())
    return (total - y) / (len(y) - 1)


def _metrics(y_true: np.ndarray, y_hat: np.ndarray) -> ErrorMetrics:
    rmse = float(np.sqrt(mean_squared_error(y_true, y_hat)))
    mae = float(mean_absolute_error(y_true, y_hat))
    r2 = float(r2_score(y_true, y_hat))
    return ErrorMetrics(rmse=rmse, mae=mae, r2=r2)
