import numpy as np

from vannes_validation.constants import EARTH_RADIUS
from vannes_validation.io.terrain import TerrainCells
from vannes_validation.methods.numerical_reference import spherical_cell_correction
from vannes_validation.methods.vannes import vannes_station_centred_correction


def test_midpoint_spherical_cell_matches_independent_3d_quadrature_when_far_and_small():
    distance = 20_000.0
    width = 100.0
    centre_lon = distance / EARTH_RADIUS
    half_angle = 0.5 * width / EARTH_RADIUS
    solid_angle = (2.0 * half_angle) * (2.0 * np.sin(half_angle))
    cells = TerrainCells(
        longitude_deg=np.array([np.rad2deg(centre_lon)]),
        latitude_deg=np.array([0.0]),
        height_m=np.array([800.0]),
        east_m=np.array([distance]),
        north_m=np.array([0.0]),
        width_east_m=np.array([width]),
        width_north_m=np.array([width]),
        area_m2=np.array([EARTH_RADIUS**2 * solid_angle]),
        distance_m=np.array([distance]),
        azimuth_rad=np.array([np.pi / 2.0]),
        solid_angle_sr=np.array([solid_angle]),
    )
    production = vannes_station_centred_correction(
        cells,
        0.0,
        0.0,
        100.0,
        near_field_radius_m=0.0,
        quadrature_order=8,
    )
    reference = spherical_cell_correction(
        0.0,
        0.0,
        100.0,
        (centre_lon - half_angle, centre_lon + half_angle),
        (-half_angle, half_angle),
        800.0,
        order=8,
    )
    np.testing.assert_allclose(production, reference, rtol=2e-5, atol=1e-12)
