# Data provenance and limitations

| File | Role | Verified property |
|---|---|---|
| `data/raw/original_stations.xlsx` | Publication-safe station table | 1,991 records; longitude, latitude and elevation only |
| `data/reference/original_grid_200x200.grd` | Surfer 7 binary grid | 40,000 nodes |
| `data/reference/interpolated_grid_598x598.xyz` | ASCII XYZ grid | 357,604 nodes |

The source equations from Vannes (2011) are cited bibliographically but the
supplied PDF excerpt is not redistributed because no separate permission from
its copyright holder has been established.

Every third node of the 598 × 598 XYZ grid reproduces the 200 × 200 grid to
floating-point precision. The dense file is therefore a resampling, not a DEM
with nine times more independent information. Negative resampling overshoots
are clipped to sea level for all methods in the controlled comparison because
the Vannes program applies that rule. This is a declared convention, not a
claim that bathymetry is absent.

## Publication rights

Cameroon Mining Group S.A. authorization 00352/10/09/26/CMG permits scientific
publication and permanent archiving of the pseudonymized station coordinates,
station elevations, terrain grid, calculated terrain corrections, figures and
statistics used in this study. The authorization excludes raw gravity readings,
client or prospect names, commercial interpretations and operational metadata.
Accordingly, the internal free-air anomaly column has been removed from both
the public station workbook and the station-level output table. Ownership of
the source data remains with Cameroon Mining Group S.A.; see
`DATA_USE_NOTICE.md`.

