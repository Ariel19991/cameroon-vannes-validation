"""Comparison statistics used by both scripts and tests."""

from __future__ import annotations

import numpy as np
import pandas as pd


def compare_to_reference(values, reference, label: str) -> dict[str, float | int | str]:
    values = np.asarray(values, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    if values.shape != reference.shape or values.size == 0:
        raise ValueError("values and reference must be non-empty arrays of equal shape")
    difference = values - reference
    return {
        "comparison": label,
        "bias_mgal": float(np.mean(difference)),
        "mae_mgal": float(np.mean(np.abs(difference))),
        "rmse_mgal": float(np.sqrt(np.mean(difference**2))),
        "median_difference_mgal": float(np.median(difference)),
        "p95_absolute_difference_mgal": float(np.quantile(np.abs(difference), 0.95)),
        "maximum_absolute_difference_mgal": float(np.max(np.abs(difference))),
        "pearson_r": float(pd.Series(values).corr(pd.Series(reference))),
        "negative_count": int(np.count_nonzero(values < 0.0)),
    }
