# Changelog

## [1.0.0] - 2026-11-24

### Added

- Download of the QsarDB archives for the Gajewicz et al. (2014) HaCaT table and the Puzyn et al. (2011) *Escherichia coli* table, with SHA-256 checks.
- Leave-one-out ordinary least squares on the deposited descriptors, compared with a leave-one-out mean baseline.
- `results/loocv_observed_vs_predicted.png`, `results/metrics.json`, and `results/loocv_predictions.csv`.
- Pinned `requirements.txt`, pytest, and a GitHub Actions workflow that lints and tests.

### Removed

- The 8-row illustrative table in `data/nanomaterial_dataset.csv`.
- The weighted risk score. That score used cell viability, which is an outcome, as an input.
