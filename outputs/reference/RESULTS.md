# Validation of the station-centred Vannes adaptation

This controlled calculation uses all 1,991 stations and the complete
information-equivalent 200 x 200 legacy terrain grid. Every method receives the
same cells, density, station elevation and treatment of negative terrain.

| Variant compared with Nagy | Bias (mGal) | MAE (mGal) | RMSE (mGal) | P95 absolute (mGal) | Maximum absolute (mGal) | Pearson r | Negative values |
|---|---:|---:|---:|---:|---:|---:|---:|
| Original Vannes residual | -1.866119 | 2.323545 | 5.098439 | 9.448417 | 73.950194 | -0.274040 | 604 |
| Station-centred adaptation | -0.002618 | 0.002619 | 0.006768 | 0.011113 | 0.127825 | 0.999999 | 0 |

For the adapted result, 93.97% of stations are within 0.01 mGal,
99.75% within 0.05 mGal and 99.95% within 0.10 mGal
of the Nagy reference.

The adaptation is not a regression or empirical rescaling. It changes the
reference condensation surface from the geoid to the sphere through the station,
so that the residual mass and reference mass represent the same column. Near-field
cells use exact Nagy prisms to regularize the self cell.

This run validates the mathematical direction of the correction on the supplied
legacy grid. Final publication values require rerunning the normal hybrid-FABDEM
pipeline and its radial-convergence analysis with `vannes_station_centred_mgal` enabled.
