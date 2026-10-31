# Cameroon Vannes validation

This repository corrects and tests the point-by-point formulation transcribed
from Nathalie Vannes (2011) for station terrain corrections. It preserves the
original algorithm as a baseline and introduces a clearly named derivative:
the **station-centred residual-condensation method** (SCRC).

The target observable is the positive classical terrain correction. Nagy's
exact rectangular-prism solution is the conventional numerical reference.
Hammer compartments and centroid point masses are diagnostic approximations;
they are not interchangeable reference solutions.

## Main result on the supplied data

All 1,991 stations were evaluated on the same information-equivalent 200 × 200
terrain grid, with density 2,670 kg/m³, negative heights clipped to zero, the
same complete rectangular domain and no station-dependent truncation.

| Method minus Nagy | Bias (mGal) | MAE | RMSE | P95 absolute | Maximum absolute | r | Negative results |
|---|---:|---:|---:|---:|---:|---:|---:|
| Original Vannes residual | -1.866119 | 2.323545 | 5.098439 | 9.448417 | 73.950194 | -0.274040 | 604 |
| Station-centred SCRC | -0.002618 | 0.002619 | 0.006768 | 0.011113 | 0.127825 | 0.999999 | 0 |

SCRC is within 0.01 mGal of Nagy at 93.97% of stations and within 0.05 mGal at
99.75%. This is strong evidence of numerical equivalence under the controlled
conditions, but not yet fully independent validation: the default first 5 km
use Nagy prisms to regularize the self-cell and near field.

## Reproduce

```bash
conda env create -f environment.yml
conda activate cameroon-vannes-validation
python -m pip install -e .
pytest
python scripts/run_validation.py
python scripts/run_sensitivity.py
```

The full run writes station-level CSV files and figures under `outputs/`.

## Data and rights

The software is released under the MIT License. The data are not covered by
that software license. Cameroon Mining Group S.A. authorized the publication
and permanent scientific archiving of the pseudonymized station coordinates,
station elevations, terrain grid, calculated terrain corrections, figures and
statistics for this study (authorization 00352/10/09/26/CMG). Ownership of the
source data remains with Cameroon Mining Group S.A. Raw gravity measurements,
client or prospect names, commercial interpretations and operational metadata
are excluded. See [DATA_USE_NOTICE.md](DATA_USE_NOTICE.md).

The public station workbook therefore contains only longitude, latitude and
station elevation. The free-air anomaly column present in an internal working
file is deliberately excluded from this release.

## Scientific safeguards

- The original Vannes kernel is never silently replaced.
- Every comparison uses identical terrain cells, density, station height,
  sea-level convention and integration domain.
- No regression or empirical calibration is applied to force agreement.
- Sensitivity to quadrature order, near-field radius and grid stride is
  reported separately.
- Credibility is conditional until an independent tesseroid implementation,
  an external DEM with adequate margins and field benchmarks are included.

See [the derivation](docs/derivation.md), [validation protocol](docs/validation_protocol.md),
[credibility verdict](docs/credibility.md), and [data provenance](docs/data.md).

## Author

- Ariel Bounou Guimking — Faculty of Physics and Technology, Novosibirsk State
  Technical University; ORCID: [0009-0003-1847-2436](https://orcid.org/0009-0003-1847-2436)

## Citation

Please cite the numbered GitHub/Zenodo release described in `CITATION.cff`.

