"""Original Vannes kernel and the station-centred corrected formulation.

The original transcription is retained to make the scientific difference
auditable. ``vannes_station_centred_correction`` is a derived method and must
not be presented as an unchanged implementation of Vannes (2011).
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from ..constants import DEFAULT_DENSITY, DEG_TO_RAD, EARTH_RADIUS, SI_TO_MGAL, G
from ..io.terrain import TerrainCells
from .prism import nagy_prism_correction


def dnewtc_components(
    computation_radius_m: float,
    cosine_psi: np.ndarray,
    lower_radius_m: float,
    upper_radius_m: np.ndarray,
    height_m: np.ndarray,
    geoid_radius_m: float = EARTH_RADIUS,
) -> tuple[np.ndarray, np.ndarray]:
    """Return Vannes' topographic and condensed-mass kernel components."""
    r = float(computation_radius_m)
    x = np.asarray(cosine_psi, dtype=np.float64)
    lower = float(lower_radius_m)
    upper = np.asarray(upper_radius_m, dtype=np.float64)
    height = np.asarray(height_m, dtype=np.float64)
    t1 = lower / r
    t2 = upper / r
    d1 = np.sqrt(np.maximum(1.0 - 2.0 * t1 * x + t1 * t1, 0.0))
    d2 = np.sqrt(np.maximum(1.0 - 2.0 * t2 * x + t2 * t2, 0.0))
    numerator = t2 - x + d2
    denominator = t1 - x + d1
    ratio = np.divide(
        numerator,
        denominator,
        out=np.ones_like(numerator),
        where=np.abs(denominator) > np.finfo(float).eps,
    )
    log_term = (3.0 * x * x - 1.0) * np.log(np.maximum(np.abs(ratio), np.finfo(float).tiny))
    p1 = (x * t1 * t1 + 3.0 * x + (1.0 - 6.0 * x * x) * t1) / d1
    p2 = (x * t2 * t2 + 3.0 * x + (1.0 - 6.0 * x * x) * t2) / d2
    topographic = r * (p2 - p1 + log_term)

    distance = np.sqrt(np.maximum(r * r + geoid_radius_m**2 - 2.0 * r * geoid_radius_m * x, 0.0))
    sigma = height * (1.0 + height / geoid_radius_m + height * height / (3.0 * geoid_radius_m**2))
    condensed = np.divide(
        -(geoid_radius_m**2) * sigma * (r - geoid_radius_m * x),
        distance**3,
        out=np.zeros_like(distance),
        where=distance > 0.0,
    )
    return topographic, condensed


def dnewtc(
    computation_radius_m: float,
    cosine_psi: np.ndarray,
    lower_radius_m: float,
    upper_radius_m: np.ndarray,
    height_m: np.ndarray,
    geoid_radius_m: float = EARTH_RADIUS,
) -> np.ndarray:
    """First derivative of the residual Newton kernel from Vannes Appendix A.7.2."""
    topographic, condensed = dnewtc_components(
        computation_radius_m,
        cosine_psi,
        lower_radius_m,
        upper_radius_m,
        height_m,
        geoid_radius_m,
    )
    return topographic - condensed


@lru_cache(maxsize=8)
def _gauss_legendre(order: int) -> tuple[np.ndarray, np.ndarray]:
    if not 2 <= order <= 16:
        raise ValueError("quadrature_order must be between 2 and 16")
    return np.polynomial.legendre.leggauss(order)


def _masked_cells(cells: TerrainCells, mask: np.ndarray) -> TerrainCells:
    return TerrainCells(
        **{name: getattr(cells, name)[mask] for name in TerrainCells.__dataclass_fields__}
    )


def vannes_station_centred_correction(
    cells: TerrainCells,
    station_lon_deg: float,
    station_lat_deg: float,
    station_height_m: float,
    density_kg_m3: float = DEFAULT_DENSITY,
    geoid_radius_m: float = EARTH_RADIUS,
    near_field_radius_m: float = 5_000.0,
    quadrature_order: int = 4,
) -> float:
    """Return a station-centred spherical terrain correction in mGal.

    The original Vannes residual condenses absolute topographic masses onto the
    geoid. Here, the signed column between station and terrain height is
    integrated on the sphere and the same residual mass is condensed onto the
    sphere through the station. This removes attraction caused solely by Earth
    curvature. Near cells use exact Nagy prisms to regularize the self cell and
    enforce the local planar limit.
    """
    if near_field_radius_m < 0.0:
        raise ValueError("near_field_radius_m must be non-negative")
    nodes, weights = _gauss_legendre(quadrature_order)
    if len(cells) == 0:
        return 0.0

    station_radius = geoid_radius_m + float(station_height_m)
    if station_radius <= 0.0:
        raise ValueError("station radius must be positive")
    source_lat = cells.latitude_deg * DEG_TO_RAD
    source_lon = cells.longitude_deg * DEG_TO_RAD
    station_lat = station_lat_deg * DEG_TO_RAD
    station_lon = station_lon_deg * DEG_TO_RAD
    cosine_psi = np.sin(source_lat) * np.sin(station_lat) + np.cos(source_lat) * np.cos(
        station_lat
    ) * np.cos(source_lon - station_lon)
    cosine_psi = np.clip(cosine_psi, -1.0, 1.0)

    local_scale = 2.0 * float(
        np.max(np.maximum(cells.width_east_m, cells.width_north_m), initial=0.0)
    )
    near = cells.distance_m < max(float(near_field_radius_m), local_scale)
    total_mgal = 0.0
    if np.any(near):
        total_mgal += nagy_prism_correction(
            _masked_cells(cells, near), station_height_m, density_kg_m3
        )

    far = ~near
    if not np.any(far):
        return float(total_mgal)
    terrain_height = cells.height_m[far]
    lower_radius = geoid_radius_m + np.minimum(terrain_height, station_height_m)
    upper_radius = geoid_radius_m + np.maximum(terrain_height, station_height_m)
    sign = np.sign(terrain_height - station_height_m)
    if np.any(lower_radius <= 0.0):
        raise ValueError("terrain radius must be positive")
    midpoint = 0.5 * (lower_radius + upper_radius)
    half_range = 0.5 * (upper_radius - lower_radius)
    x = cosine_psi[far]
    direct = np.zeros_like(x)
    for node, weight in zip(nodes, weights):
        source_radius = midpoint + half_range * node
        distance_squared = np.maximum(
            station_radius**2 + source_radius**2 - 2.0 * station_radius * source_radius * x,
            np.finfo(float).tiny,
        )
        direct += (
            weight * source_radius**2 * (source_radius * x - station_radius) / distance_squared**1.5
        )
    direct *= sign * half_range

    signed_volume_per_sr = sign * (upper_radius**3 - lower_radius**3) / 3.0
    reference_distance_squared = np.maximum(
        2.0 * station_radius**2 * (1.0 - x), np.finfo(float).tiny
    )
    condensed = signed_volume_per_sr * station_radius * (x - 1.0) / reference_distance_squared**1.5
    residual_integral = np.sum((direct - condensed) * cells.solid_angle_sr[far], dtype=np.float64)
    total_mgal += G * density_kg_m3 * residual_integral * SI_TO_MGAL
    return float(total_mgal)


# Backward-compatible scientific alias used by the first controlled run.
vannes_adapted_correction = vannes_station_centred_correction


def vannes_direct_effect(
    cells: TerrainCells,
    station_lon_deg: float,
    station_lat_deg: float,
    station_height_m: float,
    density_kg_m3: float = DEFAULT_DENSITY,
    geoid_radius_m: float = EARTH_RADIUS,
) -> float:
    """Compute Vannes' signed direct residual effect on surface gravity in mGal."""
    source_height = np.maximum(cells.height_m, 0.0)
    station_height = max(float(station_height_m), 0.0)
    source_lat = cells.latitude_deg * DEG_TO_RAD
    source_lon = cells.longitude_deg * DEG_TO_RAD
    station_lat = station_lat_deg * DEG_TO_RAD
    station_lon = station_lon_deg * DEG_TO_RAD
    cosine_psi = np.sin(source_lat) * np.sin(station_lat) + np.cos(source_lat) * np.cos(
        station_lat
    ) * np.cos(source_lon - station_lon)
    cosine_psi = np.clip(cosine_psi, -1.0, 1.0)
    nonself = cosine_psi < 1.0 - 1e-15
    if not np.any(nonself):
        return 0.0
    x = cosine_psi[nonself]
    r = geoid_radius_m + station_height
    source_kernel = dnewtc(
        r,
        x,
        geoid_radius_m,
        geoid_radius_m + source_height[nonself],
        source_height[nonself],
        geoid_radius_m,
    )
    reference_height = np.full_like(x, station_height)
    reference_kernel = dnewtc(
        r,
        x,
        geoid_radius_m,
        geoid_radius_m + reference_height,
        reference_height,
        geoid_radius_m,
    )
    integral = np.sum(
        density_kg_m3 * (source_kernel - reference_kernel) * cells.solid_angle_sr[nonself],
        dtype=np.float64,
    )
    return float(G * integral * SI_TO_MGAL)
