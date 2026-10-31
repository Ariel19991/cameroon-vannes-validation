from pathlib import Path

import numpy as np

from vannes_validation.io.grid import load_xyz_grid
from vannes_validation.io.stations import audit_stations, load_stations
from vannes_validation.io.surfer import read_surfer7_binary

ROOT = Path(__file__).resolve().parents[1]


def test_supplied_station_count_and_known_quality_flags():
    audit = audit_stations(load_stations(ROOT / "data/raw/original_stations.xlsx"))
    assert audit["count"] == 1991
    assert audit["zero_altitude_count"] == 147
    assert audit["duplicate_coordinate_count"] == 12


def test_dense_grid_is_only_a_factor_three_resampling():
    x0, y0, z0 = read_surfer7_binary(ROOT / "data/reference/original_grid_200x200.grd")
    x1, y1, z1 = load_xyz_grid(
        ROOT / "data/reference/interpolated_grid_598x598.xyz",
        stride=3,
        clip_negative=False,
    )
    assert z0.shape == (200, 200)
    assert z1.shape == z0.shape
    np.testing.assert_allclose(x1, x0, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(y1, y0, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(z1, z0, rtol=0.0, atol=1e-9)
