"""Common terrain-cell representation shared by every numerical method."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TerrainCells:
    """Terrain geometry with both local Cartesian and spherical measures."""

    longitude_deg: np.ndarray
    latitude_deg: np.ndarray
    height_m: np.ndarray
    east_m: np.ndarray
    north_m: np.ndarray
    width_east_m: np.ndarray
    width_north_m: np.ndarray
    area_m2: np.ndarray
    distance_m: np.ndarray
    azimuth_rad: np.ndarray
    solid_angle_sr: np.ndarray

    def __len__(self) -> int:
        return int(self.height_m.size)

