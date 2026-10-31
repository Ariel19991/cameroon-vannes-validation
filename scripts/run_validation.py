#!/usr/bin/env python3
"""Evaluate original and station-centred Vannes variants on a common legacy grid."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from vannes_validation.constants import DEG_TO_RAD, EARTH_RADIUS
from vannes_validation.io.stations import load_stations
from vannes_validation.io.terrain import TerrainCells
from vannes_validation.methods.hammer import hammer_sector_correction, logarithmic_edges
from vannes_validation.methods.point_mass import point_mass_correction
from vannes_validation.methods.prism import nagy_prism_correction
from vannes_validation.methods.vannes import (
    vannes_direct_effect,
    vannes_station_centred_correction,
)

METHODS = {
    "hammer_mgal": "Hammer-style",
    "nagy_prism_mgal": "Nagy prism",
    "point_mass_mgal": "Point mass",
    "vannes_direct_mgal": "Vannes original residual",
    "vannes_station_centred_mgal": "Station-centred SCRC",
}


def load_information_equivalent_grid(path: Path, factor: int = 3) -> tuple[np.ndarray, ...]:
    """Read the 598x598 grid and retain its information-equivalent 200x200 nodes."""
    data = np.loadtxt(path, dtype=np.float64)
    longitude = np.unique(data[:, 0])
    latitude = np.unique(data[:, 1])
    if longitude.size * latitude.size != data.shape[0]:
        raise ValueError("Reference XYZ file is not a complete regular grid")
    ordered = data[np.lexsort((data[:, 0], data[:, 1]))]
    height = ordered[:, 2].reshape(latitude.size, longitude.size)
    if (longitude.size - 1) % factor or (latitude.size - 1) % factor:
        raise ValueError("Grid dimensions are incompatible with the requested factor")
    longitude = longitude[::factor]
    latitude = latitude[::factor]
    # Vannes clips negative topography to sea level. Apply that same documented
    # convention to every method in this controlled comparison.
    height = np.maximum(height[::factor, ::factor], 0.0)
    return longitude, latitude, height


def cells_for_station(
    longitude_axis: np.ndarray,
    latitude_axis: np.ndarray,
    height_grid: np.ndarray,
    station_lon: float,
    station_lat: float,
) -> TerrainCells:
    lon_grid, lat_grid = np.meshgrid(longitude_axis, latitude_axis)
    longitude = lon_grid.ravel()
    latitude = lat_grid.ravel()
    height = height_grid.ravel()
    step_lon = float(np.median(np.diff(longitude_axis)))
    step_lat = float(np.median(np.diff(latitude_axis)))
    dlon_rad = (longitude - station_lon) * DEG_TO_RAD
    dlat_rad = (latitude - station_lat) * DEG_TO_RAD
    mean_lat_rad = 0.5 * (latitude + station_lat) * DEG_TO_RAD
    east = EARTH_RADIUS * np.cos(mean_lat_rad) * dlon_rad
    north = EARTH_RADIUS * dlat_rad
    distance = np.hypot(east, north)
    step_lon_rad = step_lon * DEG_TO_RAD
    step_lat_rad = step_lat * DEG_TO_RAD
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


def pairwise_summary(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    reference = data["nagy_prism_mgal"]
    for column, label in METHODS.items():
        if column == "nagy_prism_mgal":
            continue
        difference = data[column] - reference
        rows.append(
            {
                "comparison": f"{label} - Nagy prism",
                "bias_mgal": difference.mean(),
                "mae_mgal": difference.abs().mean(),
                "rmse_mgal": np.sqrt(np.mean(difference**2)),
                "median_difference_mgal": difference.median(),
                "p95_absolute_difference_mgal": difference.abs().quantile(0.95),
                "maximum_absolute_difference_mgal": difference.abs().max(),
                "pearson_r": data[column].corr(reference),
                "negative_count": int((data[column] < 0.0).sum()),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stations", type=Path, default=Path("data/raw/original_stations.xlsx"))
    parser.add_argument(
        "--grid", type=Path, default=Path("data/reference/interpolated_grid_598x598.xyz")
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/reference"))
    parser.add_argument("--density", type=float, default=2670.0)
    parser.add_argument("--near-field-km", type=float, default=5.0)
    parser.add_argument("--quadrature-order", type=int, default=4)
    args = parser.parse_args()

    stations = load_stations(args.stations)
    longitude_axis, latitude_axis, height_grid = load_information_equivalent_grid(args.grid)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir = args.output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, float | int]] = []
    for index in range(len(stations)):
        lon = float(stations.longitude_deg[index])
        lat = float(stations.latitude_deg[index])
        altitude = float(stations.altitude_m[index])
        cells = cells_for_station(longitude_axis, latitude_axis, height_grid, lon, lat)
        edges = logarithmic_edges(float(cells.distance_m.max()) * 1.000001)
        rows.append(
            {
                "excel_row": int(stations.excel_row[index]),
                "longitude_deg": lon,
                "latitude_deg": lat,
                "station_altitude_m": altitude,
                "cell_count": len(cells),
                "maximum_distance_km": float(cells.distance_m.max() / 1000.0),
                "hammer_mgal": hammer_sector_correction(
                    cells, altitude, edges, density_kg_m3=args.density
                ),
                "nagy_prism_mgal": nagy_prism_correction(cells, altitude, args.density),
                "point_mass_mgal": point_mass_correction(cells, altitude, args.density),
                "vannes_direct_mgal": vannes_direct_effect(cells, lon, lat, altitude, args.density),
                "vannes_station_centred_mgal": vannes_station_centred_correction(
                    cells,
                    lon,
                    lat,
                    altitude,
                    args.density,
                    near_field_radius_m=args.near_field_km * 1000.0,
                    quadrature_order=args.quadrature_order,
                ),
            }
        )
        if (index + 1) % 100 == 0 or index + 1 == len(stations):
            print(f"Completed {index + 1}/{len(stations)} stations", flush=True)

    results = pd.DataFrame(rows)
    results["original_minus_nagy_mgal"] = results["vannes_direct_mgal"] - results["nagy_prism_mgal"]
    results["station_centred_minus_nagy_mgal"] = (
        results["vannes_station_centred_mgal"] - results["nagy_prism_mgal"]
    )
    results.to_csv(args.output_dir / "station_comparison.csv", index=False)
    summary = pairwise_summary(results)
    summary.to_csv(args.output_dir / "pairwise_summary.csv", index=False)

    figure, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharex=True, sharey=True)
    limits = [
        min(results.nagy_prism_mgal.min(), results.vannes_direct_mgal.min()),
        max(results.nagy_prism_mgal.max(), results.vannes_direct_mgal.max()),
    ]
    axes[0].scatter(results.nagy_prism_mgal, results.vannes_direct_mgal, s=8, alpha=0.55)
    axes[0].plot(limits, limits, "k--", linewidth=1)
    axes[0].set_title("Original Vannes residual")
    axes[1].scatter(results.nagy_prism_mgal, results.vannes_station_centred_mgal, s=8, alpha=0.55)
    axes[1].plot(limits, limits, "k--", linewidth=1)
    axes[1].set_title("Station-centred Vannes adaptation")
    for axis in axes:
        axis.set_xlabel("Nagy prism (mGal)")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("Vannes result (mGal)")
    figure.tight_layout()
    figure.savefig(figure_dir / "vannes_before_after.png", dpi=180)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(8, 5.2))
    axis.scatter(
        results.station_altitude_m,
        results.original_minus_nagy_mgal,
        s=8,
        alpha=0.45,
        label="Original - Nagy",
    )
    axis.scatter(
        results.station_altitude_m,
        results.station_centred_minus_nagy_mgal,
        s=8,
        alpha=0.45,
        label="SCRC - Nagy",
    )
    axis.axhline(0.0, color="black", linewidth=1)
    axis.set_xlabel("Station altitude (m)")
    axis.set_ylabel("Difference from Nagy (mGal)")
    axis.set_title("Removal of altitude-dependent Vannes bias")
    axis.grid(alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(figure_dir / "difference_vs_altitude.png", dpi=180)
    plt.close(figure)

    metadata = {
        "station_count": len(stations),
        "grid_shape": list(height_grid.shape),
        "grid_node_count": int(height_grid.size),
        "grid_policy": "Every third node of 598x598 Surfer grid; negative heights clipped to zero for every method",
        "domain_policy": "Complete rectangular legacy grid for every station; no radial truncation",
        "density_kg_m3": args.density,
        "near_field_radius_km": args.near_field_km,
        "quadrature_order": args.quadrature_order,
        "important_limitation": "Controlled common-grid comparison, not the final hybrid-FABDEM 200 km production run",
    }
    (args.output_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    corrected = summary[summary.comparison.str.startswith("Station-centred SCRC")].iloc[0]
    original = summary[summary.comparison.str.startswith("Vannes original")].iloc[0]
    within_001 = float((results.station_centred_minus_nagy_mgal.abs() <= 0.01).mean())
    within_005 = float((results.station_centred_minus_nagy_mgal.abs() <= 0.05).mean())
    within_010 = float((results.station_centred_minus_nagy_mgal.abs() <= 0.10).mean())
    report = f"""# Validation of the station-centred Vannes adaptation

This controlled calculation uses all {len(stations):,} stations and the complete
information-equivalent 200 x 200 legacy terrain grid. Every method receives the
same cells, density, station elevation and treatment of negative terrain.

| Variant compared with Nagy | Bias (mGal) | MAE (mGal) | RMSE (mGal) | P95 absolute (mGal) | Maximum absolute (mGal) | Pearson r | Negative values |
|---|---:|---:|---:|---:|---:|---:|---:|
| Original Vannes residual | {original.bias_mgal:.6f} | {original.mae_mgal:.6f} | {original.rmse_mgal:.6f} | {original.p95_absolute_difference_mgal:.6f} | {original.maximum_absolute_difference_mgal:.6f} | {original.pearson_r:.6f} | {int(original.negative_count)} |
| Station-centred adaptation | {corrected.bias_mgal:.6f} | {corrected.mae_mgal:.6f} | {corrected.rmse_mgal:.6f} | {corrected.p95_absolute_difference_mgal:.6f} | {corrected.maximum_absolute_difference_mgal:.6f} | {corrected.pearson_r:.6f} | {int(corrected.negative_count)} |

For the adapted result, {100 * within_001:.2f}% of stations are within 0.01 mGal,
{100 * within_005:.2f}% within 0.05 mGal and {100 * within_010:.2f}% within 0.10 mGal
of the Nagy reference.

The adaptation is not a regression or empirical rescaling. It changes the
reference condensation surface from the geoid to the sphere through the station,
so that the residual mass and reference mass represent the same column. Near-field
cells use exact Nagy prisms to regularize the self cell.

This run validates the mathematical direction of the correction on the supplied
legacy grid. Final publication values require rerunning the normal hybrid-FABDEM
pipeline and its radial-convergence analysis with `vannes_station_centred_mgal` enabled.
"""
    (args.output_dir / "RESULTS.md").write_text(report, encoding="utf-8")
    report_fr = f"""# Validation de l'adaptation de Vannes centrée sur la station

Ce calcul contrôlé utilise les {len(stations):,} stations et la totalité de la
grille historique informationnellement équivalente de 200 x 200 nœuds. Toutes
les méthodes reçoivent les mêmes cellules, la même densité, la même altitude de
station et le même traitement des altitudes négatives.

| Variante comparée à Nagy | Biais (mGal) | MAE (mGal) | RMSE (mGal) | P95 absolu (mGal) | Maximum absolu (mGal) | Pearson r | Valeurs négatives |
|---|---:|---:|---:|---:|---:|---:|---:|
| Vannes résiduelle originale | {original.bias_mgal:.6f} | {original.mae_mgal:.6f} | {original.rmse_mgal:.6f} | {original.p95_absolute_difference_mgal:.6f} | {original.maximum_absolute_difference_mgal:.6f} | {original.pearson_r:.6f} | {int(original.negative_count)} |
| Adaptation centrée sur la station | {corrected.bias_mgal:.6f} | {corrected.mae_mgal:.6f} | {corrected.rmse_mgal:.6f} | {corrected.p95_absolute_difference_mgal:.6f} | {corrected.maximum_absolute_difference_mgal:.6f} | {corrected.pearson_r:.6f} | {int(corrected.negative_count)} |

Pour la méthode adaptée, {100 * within_001:.2f} % des stations diffèrent de Nagy
de moins de 0,01 mGal, {100 * within_005:.2f} % de moins de 0,05 mGal et
{100 * within_010:.2f} % de moins de 0,10 mGal.

L'adaptation n'est ni une régression ni un réétalonnage empirique. Elle déplace
la surface de condensation de référence du géoïde vers la sphère passant par la
station, afin que la masse résiduelle et la masse de référence représentent la
même colonne. Les cellules proches sont calculées par les prismes exacts de Nagy
afin de régulariser la cellule contenant la station.

Ce résultat valide l'orientation mathématique de la correction sur la grille
historique fournie. Les valeurs destinées à une publication nécessitent encore
un nouveau calcul de la chaîne MNT hybride/FABDEM et de sa convergence radiale,
avec la colonne `vannes_station_centred_mgal` activée.
"""
    (args.output_dir / "RESULTS_FR.md").write_text(report_fr, encoding="utf-8")


if __name__ == "__main__":
    main()
