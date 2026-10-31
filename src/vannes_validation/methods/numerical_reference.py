"""Independent 3-D spherical quadrature for synthetic validation cases."""

from __future__ import annotations

import numpy as np

from ..constants import DEFAULT_DENSITY, EARTH_RADIUS, SI_TO_MGAL, G


def spherical_cell_correction(
    station_lon_rad: float,
    station_lat_rad: float,
    station_height_m: float,
    cell_lon_bounds_rad: tuple[float, float],
    cell_lat_bounds_rad: tuple[float, float],
    terrain_height_m: float,
    density_kg_m3: float = DEFAULT_DENSITY,
    geoid_radius_m: float = EARTH_RADIUS,
    order: int = 8,
) -> float:
    """Integrate one residual spherical voxel and its station shell directly.

    Unlike the production estimator, this function integrates longitude,
    latitude and radius. It is intentionally slower and is used as an
    independent synthetic benchmark, not as a production method.
    """
    if order < 2:
        raise ValueError("order must be at least two")
    nodes, weights = np.polynomial.legendre.leggauss(order)

    def mapped(bounds):
        midpoint = 0.5 * (bounds[0] + bounds[1])
        half = 0.5 * (bounds[1] - bounds[0])
        return midpoint + half * nodes, half * weights

    longitudes, lon_weights = mapped(cell_lon_bounds_rad)
    latitudes, lat_weights = mapped(cell_lat_bounds_rad)
    station_radius = geoid_radius_m + station_height_m
    lower = geoid_radius_m + min(station_height_m, terrain_height_m)
    upper = geoid_radius_m + max(station_height_m, terrain_height_m)
    sign = float(np.sign(terrain_height_m - station_height_m))
    radii, radial_weights = mapped((lower, upper))
    signed_volume_per_sr = sign * (upper**3 - lower**3) / 3.0

    station_unit = np.array(
        [
            np.cos(station_lat_rad) * np.cos(station_lon_rad),
            np.cos(station_lat_rad) * np.sin(station_lon_rad),
            np.sin(station_lat_rad),
        ]
    )
    station_position = station_radius * station_unit
    integral = 0.0
    for latitude, latitude_weight in zip(latitudes, lat_weights):
        cosine_latitude = np.cos(latitude)
        for longitude, longitude_weight in zip(longitudes, lon_weights):
            source_unit = np.array(
                [
                    cosine_latitude * np.cos(longitude),
                    cosine_latitude * np.sin(longitude),
                    np.sin(latitude),
                ]
            )
            direct = 0.0
            for radius, radial_weight in zip(radii, radial_weights):
                displacement = radius * source_unit - station_position
                distance = np.linalg.norm(displacement)
                direct += (
                    radial_weight * radius**2 * np.dot(displacement, station_unit) / distance**3
                )
            shell_displacement = station_radius * (source_unit - station_unit)
            shell_distance = np.linalg.norm(shell_displacement)
            condensed = (
                signed_volume_per_sr * np.dot(shell_displacement, station_unit) / shell_distance**3
            )
            integral += (
                (sign * direct - condensed) * cosine_latitude * latitude_weight * longitude_weight
            )
    return float(G * density_kg_m3 * integral * SI_TO_MGAL)
