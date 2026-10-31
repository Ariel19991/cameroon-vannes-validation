"""Load the supplied legacy XYZ grid and construct common terrain cells."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..constants import DEG_TO_RAD, EARTH_RADIUS
from .terrain import TerrainCells


def load_xyz_grid(path: str | Path, stride: int = 1, clip_negative: bool = True):
    """Return longitude, latitude and height arrays from a regular XYZ grid."""
    if stride < 1:
        raise ValueError("stride must be a positive integer")
    data = np.loadtxt(path, dtype=np.float64)
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError("Expected an ASCII XYZ grid with at least three columns")
    longitude = np.unique(data[:, 0])
    latitude = np.unique(data[:, 1])
    if longitude.size * latitude.size != data.shape[0]:
        raise ValueError("XYZ points do not form a complete regular grid")
    ordered = data[np.lexsort((data[:, 0], data[:, 1]))]
    height = ordered[:, 2].reshape(latitude.size, longitude.size)
    longitude = longitude[::stride]
    latitude = latitude[::stride]
    height = height[::stride, ::stride]
    if clip_negative:
        height = np.maximum(height, 0.0)
    return longitude, latitude, height


def cells_for_station(
    longitude_axis: np.ndarray,
    latitude_axis: np.ndarray,
    height_grid: np.ndarray,
    station_lon_deg: float,
    station_lat_deg: float,
) -> TerrainCells:
    """Build one common cell collection for every compared method."""
    longitude_axis = np.asarray(longitude_axis, dtype=np.float64)
    latitude_axis = np.asarray(latitude_axis, dtype=np.float64)
    height_grid = np.asarray(height_grid, dtype=np.float64)
    if height_grid.shape != (latitude_axis.size, longitude_axis.size):
        raise ValueError("height_grid shape must be (latitude, longitude)")
    if longitude_axis.size < 2 or latitude_axis.size < 2:
        raise ValueError("At least two nodes are required on each axis")
    if np.any(np.diff(longitude_axis) <= 0) or np.any(np.diff(latitude_axis) <= 0):
        raise ValueError("Grid axes must be strictly increasing")

    lon_grid, lat_grid = np.meshgrid(longitude_axis, latitude_axis)
    longitude = lon_grid.ravel()
    latitude = lat_grid.ravel()
    height = height_grid.ravel()
    step_lon_rad = float(np.median(np.diff(longitude_axis))) * DEG_TO_RAD
    step_lat_rad = float(np.median(np.diff(latitude_axis))) * DEG_TO_RAD
    dlon_rad = (longitude - station_lon_deg) * DEG_TO_RAD
    dlat_rad = (latitude - station_lat_deg) * DEG_TO_RAD
    mean_lat_rad = 0.5 * (latitude + station_lat_deg) * DEG_TO_RAD
    east = EARTH_RADIUS * np.cos(mean_lat_rad) * dlon_rad
    north = EARTH_RADIUS * dlat_rad
    distance = np.hypot(east, north)
    lat_rad = latitude * DEG_TO_RAD
    solid_angle = step_lon_rad * (
        np.sin(lat_rad + 0.5 * step_lat_rad) - np.sin(lat_rad - 0.5 * step_lat_rad)
    )
    return TerrainCells(
        longitude_deg=longitude,
        latitude_deg=latitude,
        height_m=height,
        east_m=east,
        north_m=north,
        width_east_m=EARTH_RADIUS * np.cos(lat_rad) * step_lon_rad,
        width_north_m=np.full_like(height, EARTH_RADIUS * step_lat_rad),
        area_m2=EARTH_RADIUS**2 * solid_angle,
        distance_m=distance,
        azimuth_rad=np.mod(np.arctan2(east, north), 2.0 * np.pi),
        solid_angle_sr=solid_angle,
    )
