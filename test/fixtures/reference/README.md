# References

## Approved references

The directory `test/integration/known_cases/` contains one CSV file per case. Each file
records the analytical Z and Y at the 11 FEM frequencies, together with the available
FEM and PSCAD references. Its header identifies the sources and conductor order.

Each case's deviation curves, statuses and explanations are recorded in `status.toml`
beside the CSV files. `test/integration/known_cases.jl` pins the analytical values and
checks FEM deviations at the 11 FEM frequencies.

FEM decides: discrepancies against PSCAD lose to FEM. PSCAD curves are advisory and
do not need an explanation. The analytical model collapses insulators and semicons
into one generic insulation bundle. Against FEM, this is the same type of discrepancy
as the lay (helical path) correction of stranded cores.

FEM discrepancies are judged as curves. A continuous, smooth, monotonic deviation that
follows FEM's trend is generally accepted. Within 5 % of FEM a curve is accepted. A
larger deviation known from this state is advisory. The test fails when a FEM deviation
rises at any frequency. The analytical pin (rtol 1e-12) keeps the engine from changing
against itself.

## Retired inherited references

The inherited numerical payloads have all been removed. No file in this directory is
an approved initial-release reference. Necessary scenarios and independent expectation
recipes belong beside their maintained test families.

Do not regenerate rejected expectations after a failure.
