"""Readers for the legacy Surfer grids used only as a reference baseline."""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def read_surfer7_binary(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read the regular GRID and DATA sections of a Surfer 7 binary grid."""
    raw = Path(path).read_bytes()
    if raw[:4] != b"DSRB":
        raise ValueError("Expected a Surfer 7 DSRB file")
    offset = 12
    grid = None
    data = None
    while offset + 8 <= len(raw):
        tag = raw[offset : offset + 4]
        size = struct.unpack_from("<I", raw, offset + 4)[0]
        payload = offset + 8
        if payload + size > len(raw):
            raise ValueError("Truncated Surfer section")
        if tag == b"GRID":
            if size < 72:
                raise ValueError("Invalid Surfer GRID section")
            nrow, ncol = struct.unpack_from("<II", raw, payload)
            values = struct.unpack_from("<8d", raw, payload + 8)
            grid = (nrow, ncol, *values)
        elif tag == b"DATA":
            data = np.frombuffer(raw, dtype="<f8", count=size // 8, offset=payload).copy()
        offset = payload + size
    if grid is None or data is None:
        raise ValueError("Missing GRID or DATA section")
    nrow, ncol, x_min, y_min, dx, dy, _z_min, _z_max, _rotation, blank = grid
    if data.size != nrow * ncol:
        raise ValueError("Surfer DATA size does not match GRID dimensions")
    height = data.reshape(nrow, ncol)
    height[height >= blank] = np.nan
    longitude = x_min + np.arange(ncol, dtype=np.float64) * dx
    latitude = y_min + np.arange(nrow, dtype=np.float64) * dy
    return longitude, latitude, height


def sample_xyz_grid(path: str | Path, longitude: np.ndarray, latitude: np.ndarray) -> np.ndarray:
    data = np.loadtxt(path, dtype=np.float64)
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError("Expected an ASCII XYZ grid with at least three columns")
    x = np.unique(data[:, 0])
    y = np.unique(data[:, 1])
    if x.size * y.size != data.shape[0]:
        raise ValueError("XYZ points do not form a complete regular grid")
    ordered = data[np.lexsort((data[:, 0], data[:, 1]))]
    z = ordered[:, 2].reshape(y.size, x.size)
    interpolator = RegularGridInterpolator((y, x), z, bounds_error=True)
    return interpolator(np.column_stack([latitude, longitude]))
