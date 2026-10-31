import numpy as np

from vannes_validation.io.terrain import TerrainCells
from vannes_validation.methods.point_mass import point_mass_correction


def _cell(height_m):
    one = np.array([1.0])
    return TerrainCells(
        longitude_deg=np.array([0.01]),
        latitude_deg=np.array([0.0]),
        height_m=np.array([height_m]),
        east_m=np.array([1000.0]),
        north_m=np.array([0.0]),
        width_east_m=one,
        width_north_m=one,
        area_m2=np.array([10_000.0]),
        distance_m=np.array([1000.0]),
        azimuth_rad=np.array([np.pi / 2]),
        solid_angle_sr=np.array([1e-10]),
    )


def test_point_mass_is_zero_for_level_terrain_and_positive_on_both_sides():
    assert point_mass_correction(_cell(100.0), 100.0) == 0.0
    assert point_mass_correction(_cell(200.0), 100.0) > 0.0
    assert point_mass_correction(_cell(0.0), 100.0) > 0.0
