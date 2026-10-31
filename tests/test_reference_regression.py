from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_checked_in_full_station_result_passes_acceptance_gates():
    results = pd.read_csv(ROOT / "outputs/reference/station_comparison.csv")
    summary = pd.read_csv(ROOT / "outputs/reference/pairwise_summary.csv")
    row = summary.loc[summary.comparison == "Station-centred SCRC - Nagy prism"].iloc[0]
    assert len(results) == 1991
    assert abs(row.bias_mgal) <= 0.01
    assert row.rmse_mgal <= 0.05
    assert row.p95_absolute_difference_mgal <= 0.05
    assert row.pearson_r >= 0.999
    assert row.negative_count == 0
