#!/usr/bin/env python3
"""Measure dependence on near-field radius, quadrature order and grid stride."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from vannes_validation.io.grid import cells_for_station, load_xyz_grid
from vannes_validation.io.stations import load_stations
from vannes_validation.methods.prism import nagy_prism_correction
from vannes_validation.methods.vannes import vannes_station_centred_correction
from vannes_validation.statistics import compare_to_reference


def stratified_indices(count: int, sample_size: int) -> np.ndarray:
    if not 1 <= sample_size <= count:
        raise ValueError("sample_size must be between one and station count")
    return np.unique(np.linspace(0, count - 1, sample_size, dtype=int))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stations", type=Path, default=Path("data/raw/original_stations.xlsx"))
    parser.add_argument(
        "--grid", type=Path, default=Path("data/reference/interpolated_grid_598x598.xyz")
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/sensitivity"))
    parser.add_argument("--sample-size", type=int, default=80)
    parser.add_argument("--density", type=float, default=2670.0)
    args = parser.parse_args()

    stations = load_stations(args.stations)
    selected = stratified_indices(len(stations), args.sample_size)
    configurations = [
        (3, 0.0, 4),
        (3, 2.0, 4),
        (3, 5.0, 2),
        (3, 5.0, 4),
        (3, 5.0, 6),
        (3, 5.0, 8),
        (3, 10.0, 4),
        (6, 5.0, 4),
        (12, 5.0, 4),
    ]
    rows = []
    for stride, near_km, order in configurations:
        lon_axis, lat_axis, height = load_xyz_grid(args.grid, stride=stride)
        nagy_values = []
        corrected_values = []
        for index in selected:
            lon = float(stations.longitude_deg[index])
            lat = float(stations.latitude_deg[index])
            altitude = float(stations.altitude_m[index])
            cells = cells_for_station(lon_axis, lat_axis, height, lon, lat)
            nagy_values.append(nagy_prism_correction(cells, altitude, args.density))
            corrected_values.append(
                vannes_station_centred_correction(
                    cells,
                    lon,
                    lat,
                    altitude,
                    args.density,
                    near_field_radius_m=near_km * 1000.0,
                    quadrature_order=order,
                )
            )
        label = f"stride={stride}; near={near_km:g} km; order={order}"
        record = compare_to_reference(corrected_values, nagy_values, label)
        record.update(
            {
                "stride": stride,
                "near_field_km": near_km,
                "quadrature_order": order,
                "station_count": len(selected),
                "grid_nodes": int(height.size),
            }
        )
        rows.append(record)
        print(label, flush=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(rows)
    summary.to_csv(args.output_dir / "sensitivity_summary.csv", index=False)

    subset = summary[(summary.stride == 3) & (summary.quadrature_order == 4)]
    figure, axis = plt.subplots(figsize=(7.5, 4.8))
    axis.plot(subset.near_field_km, subset.rmse_mgal, marker="o")
    axis.set_xlabel("Near-field radius (km)")
    axis.set_ylabel("RMSE from Nagy (mGal)")
    axis.set_title("Sensitivity to the exact-prism near field")
    axis.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(args.output_dir / "near_field_sensitivity.png", dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
