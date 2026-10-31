# Credibility assessment

## Verdict

The unchanged Vannes residual output is **not credible as a direct substitute
for the classical station terrain correction** on the supplied data. Its sign,
reference surface and altitude dependence describe a different residual
observable.

The station-centred CRCS/SCRC estimator is **credible as a physically derived,
Nagy-equivalent hybrid method under the controlled common-grid experiment**.
It passes the predefined bias, RMSE, P95, correlation and sign gates.

## What has been demonstrated

- The correction follows from Newtonian attraction and conservation of the
  same signed residual mass; no fitted coefficient is used.
- The station reference surface removes curvature-only attraction from a level
  surface through the observation point.
- The local planar limit agrees with exact rectangular prisms.
- Agreement is observed at all 1,991 stations, not at selected examples.

## What remains unresolved

- Near-field agreement is partly imposed by the Nagy branch.
- The 598 × 598 grid is a densification of the 200 × 200 Surfer grid and does
  not add independent terrain information.
- Kriging parameters and uncertainty were not supplied.
- The grid has no exterior buffer; edge stations lack symmetric terrain.
- Constant density 2,670 kg/m³ is conventional, not locally calibrated.
- Eshagh and Bagherbandi implementations require exact source equations and
  conventions before they can be treated as reproducible reference methods.
- Independent tesseroid/direct-volume and field benchmarks remain necessary
  before claiming definitive physical validation.

