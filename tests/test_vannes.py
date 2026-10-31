import numpy as np

from vannes_validation.io.terrain import TerrainCells
from vannes_validation.methods.prism import nagy_prism_correction
from vannes_validation.methods.vannes import (
    dnewtc,
    dnewtc_components,
    vannes_direct_effect,
    vannes_station_centred_correction,
)


def _single_cell(distance_m: float, terrain_height_m: float, width_m: float = 30.0):
    radius = 6_371_000.0
    longitude = np.rad2deg(distance_m / radius)
    solid_angle = (width_m / radius) ** 2
    return TerrainCells(
        longitude_deg=np.array([longitude]),
        latitude_deg=np.array([0.0]),
        height_m=np.array([terrain_height_m]),
        east_m=np.array([distance_m]),
        north_m=np.array([0.0]),
        width_east_m=np.array([width_m]),
        width_north_m=np.array([width_m]),
        area_m2=np.array([width_m**2]),
        distance_m=np.array([distance_m]),
        azimuth_rad=np.array([np.pi / 2.0]),
        solid_angle_sr=np.array([solid_angle]),
    )


def test_dnewtc_vector_is_finite_away_from_self_point():
    x = np.cos(np.deg2rad(np.array([0.01, 0.1, 1.0])))
    result = dnewtc(6_371_500.0, x, 6_371_000.0, np.array([6_371_100.0] * 3), np.array([100.0] * 3))
    assert np.isfinite(result).all()


def test_dnewtc_components_reconstruct_original_kernel():
    x = np.cos(np.deg2rad(np.array([0.01, 0.1, 1.0])))
    upper = np.array([6_371_100.0, 6_371_500.0, 6_372_000.0])
    height = upper - 6_371_000.0
    topographic, condensed = dnewtc_components(6_371_250.0, x, 6_371_000.0, upper, height)
    np.testing.assert_allclose(
        dnewtc(6_371_250.0, x, 6_371_000.0, upper, height),
        topographic - condensed,
        rtol=0.0,
        atol=0.0,
    )


def test_station_centred_vannes_is_zero_for_level_terrain():
    cells = _single_cell(1_000.0, 250.0)
    result = vannes_station_centred_correction(cells, 0.0, 0.0, 250.0, near_field_radius_m=0.0)
    assert result == 0.0


def test_station_centred_vannes_hill_and_valley_are_positive():
    hill = _single_cell(1_000.0, 300.0)
    valley = _single_cell(1_000.0, 100.0)
    assert vannes_station_centred_correction(hill, 0.0, 0.0, 200.0) > 0.0
    assert vannes_station_centred_correction(valley, 0.0, 0.0, 200.0) > 0.0


def test_station_centred_vannes_matches_nagy_in_local_planar_limit():
    cells = _single_cell(10_000.0, 300.0)
    station_height = 100.0
    adapted = vannes_station_centred_correction(
        cells, 0.0, 0.0, station_height, near_field_radius_m=0.0, quadrature_order=6
    )
    reference = nagy_prism_correction(cells, station_height)
    np.testing.assert_allclose(adapted, reference, rtol=3e-3, atol=1e-10)


def test_station_centred_vannes_rejects_invalid_quadrature_order():
    cells = _single_cell(1_000.0, 300.0)
    with np.testing.assert_raises(ValueError):
        vannes_station_centred_correction(cells, 0.0, 0.0, 100.0, quadrature_order=1)


def test_station_centred_vannes_scales_linearly_with_density():
    cells = _single_cell(10_000.0, 450.0, width_m=100.0)
    low = vannes_station_centred_correction(
        cells, 0.0, 0.0, 100.0, density_kg_m3=2000.0, near_field_radius_m=0.0
    )
    high = vannes_station_centred_correction(
        cells, 0.0, 0.0, 100.0, density_kg_m3=3000.0, near_field_radius_m=0.0
    )
    np.testing.assert_allclose(high / low, 1.5, rtol=1e-12)


def test_radial_quadrature_is_converged_for_a_far_cell():
    cells = _single_cell(25_000.0, 1_500.0, width_m=250.0)
    order_6 = vannes_station_centred_correction(
        cells, 0.0, 0.0, 100.0, near_field_radius_m=0.0, quadrature_order=6
    )
    order_8 = vannes_station_centred_correction(
        cells, 0.0, 0.0, 100.0, near_field_radius_m=0.0, quadrature_order=8
    )
    np.testing.assert_allclose(order_6, order_8, rtol=1e-9, atol=1e-12)


def test_near_field_branch_is_exactly_the_nagy_reference():
    cells = _single_cell(1_000.0, 450.0, width_m=250.0)
    corrected = vannes_station_centred_correction(
        cells, 0.0, 0.0, 100.0, near_field_radius_m=2_000.0
    )
    nagy = nagy_prism_correction(cells, 100.0)
    np.testing.assert_allclose(corrected, nagy, rtol=0.0, atol=1e-14)


def test_original_vannes_baseline_remains_finite_and_density_linear():
    cells = _single_cell(10_000.0, 450.0, width_m=100.0)
    low = vannes_direct_effect(cells, 0.0, 0.0, 100.0, density_kg_m3=2000.0)
    high = vannes_direct_effect(cells, 0.0, 0.0, 100.0, density_kg_m3=3000.0)
    assert np.isfinite(low)
    np.testing.assert_allclose(high / low, 1.5, rtol=1e-12)
