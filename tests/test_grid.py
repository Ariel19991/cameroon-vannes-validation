import numpy as np

from vannes_validation.io.grid import cells_for_station


def test_cells_for_station_preserves_shape_and_positive_measures():
    longitude = np.array([9.0, 9.01])
    latitude = np.array([4.0, 4.01])
    height = np.array([[100.0, 120.0], [80.0, 140.0]])
    cells = cells_for_station(longitude, latitude, height, 9.005, 4.005)
    assert len(cells) == 4
    assert np.all(cells.area_m2 > 0.0)
    assert np.all(cells.solid_angle_sr > 0.0)
    np.testing.assert_allclose(cells.height_m, height.ravel())


def test_cells_for_station_rejects_shape_mismatch():
    with np.testing.assert_raises(ValueError):
        cells_for_station(np.array([0.0, 1.0]), np.array([0.0, 1.0]), np.ones((3, 2)), 0, 0)
