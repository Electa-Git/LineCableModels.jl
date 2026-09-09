```@meta
EditURL = "../literate/gauntlet.jl"
```

# Gauntlet

Gauntlet runs cases manually and retains numerical results. This page is its
documentation overview. Full tables and matrix-cell overlays are available at
the REPL from explicitly configured, persisted benchmarks.
Documentation generation never starts PSCAD, FEM, analytical calculations or
uncertainty propagation. Full-catalogue and UQ results remain stored for
analysis through the ordinary result, observation and plotting APIs.

## Recorded comparisons

Both errors are calculated element-wise by
[`compare`](@ref LineCableModels.Engine.compare), with **A = the benchmark's
reference** and **B = its candidate**. Neither role is inferred from a backend.

```math
\mathrm{NRMSE}=100\sqrt{\frac{\sum_k|B_k-A_k|^2}{\sum_k|A_k|^2}},
\qquad
\mathrm{RMS}_{\mathrm{pointwise}}=100\sqrt{\frac1N\sum_k\left|\frac{B_k-A_k}{A_k}\right|^2}.
```

**One row is one benchmark.** Reference and candidate columns identify the
backend and retained formulation record. Every recorded quantity is included;
each overview cell contains **maximum error in percent (response, excitation)**.
The maximum is across matrix entries, not a whole-matrix error. The two
normalizations can attain their maxima at different entries. Expand the
benchmark identities below each full-band table to see the terminal order.

**The full-band summary comes first.** Frequency slices follow in separate
sections, with the same column layout and their actual stored ranges.
No interpolation or additional simulations are performed.

!!! note "No recorded benchmarks selected"
    Set `LINECABLEMODELS_GAUNTLET_RESULTS` to a saved benchmark directory.
    Use `lcm gauntlet compare` to bind completed calculations explicitly.
    No calculations run during documentation generation.


## Comparison conventions

No denominator floor or solver-data modification is applied. Absolute RMS
retains the measured difference. If the reference trace lies within the recorded
observable tolerance, relative RMS is unavailable with a per-cell reason.
Pointwise normalization is also unavailable when any selected reference sample
is numerically zero; no samples are omitted. Defaults are 1e-10 Ω/m for R,
1e-15 H/m for L, 1e-12 S/m for G and 1e-16 F/m for C. Z and Y thresholds follow
R + 2πfL and G + 2πfC. Override them with `compare(...; atol=(G=..., C=...))`.
Unsupported observables and empty bands retain their separate reasons.

Saved comparisons retain the settings used when they were calculated. Historical
tables can therefore contain zero or infinite ratios under the previous settings;
recompute comparisons from their retained raw operands to apply the current
settings. The report renderer does not change stored metrics.

Deterministic values and UQ mean/std comparisons have separate sections. The
benchmark's recorded settings supply quantities, bands and normalizations.
There is no backend truth ranking or requirement for different models to agree.
Inputs, terminal order, basis and frequencies must be comparable before errors
are calculated. Report generation does not reinterpret the saved settings.

## Inspect results at the REPL

Read a completed benchmark directory, or pass a particular `snapshot.jld2` with
`load_results=true`. Both forms load the explicit reference/candidate pair and
saved analyses. Relative operand paths resolve from the snapshot directory.

```julia
using LineCableModels, DataFrames
using LineCableModels.ReportBuilder: BenchmarkTableDefinition
include(joinpath(pkgdir(LineCableModels), "gauntlet", "Gauntlet.jl"))
using .Gauntlet

benchmark = read_benchmark("/path/to/campaign/benchmark_id")
# Alternatively: read_benchmark("/path/to/snapshot.jld2"; load_results=true)
tables = report(BenchmarkTableDefinition(false), benchmark).table

tables.calculations  # IDs, actual formulation records, controls, axes and hashes
tables.comparisons   # One row per quantity/statistic/band/normalization/point pair
tables.terms         # Every matrix entry, with terminal names, units and reasons

z = filter(row -> row.quantity === :Z && row.statistic === :value &&
    row.band === :all && row.normalization === :reference_rms &&
    row.reference_point == 1 && row.candidate_point == 1, tables.comparisons)
only(z.absolute_rms)          # Complete matrix, in the units shown in absolute_unit
only(z.relative_rms_percent)  # Complete matrix; missing stays missing
only(z.reason)               # Per-entry unavailable reasons

g = filter(row -> row.quantity === :G && row.band === :all, tables.terms)
show(g; allrows=true, allcols=true)

using GLMakie
plots = LineCableModels.plot(benchmark, (Z, Y); pair=(1, 1))
```

`Z` produces the existing R/X views, and `Y` produces G/B. Each quantity has its
own matrix layout, with reference and candidate overlaid in every corresponding
cell, including both off-diagonal entries. `(R, L, G, C)` selects those quantities
directly. The existing plot controls, unit keywords and SVG export remain available.
Use `length_unit=:base` to show the native per-metre units used by the RMS tables;
the ordinary plotting default displays per-kilometre quantities.

For a single compared point pair, `pair` can be omitted. With multiple points,
choose a pair retained in the analysis; zipped and product result spaces keep
their declared pairing. `tables.calculations.axes` retains the axis descriptions
needed to identify the formulations and problems behind each point. Plotting does
not choose a replacement reference or assemble a different comparison.

`BenchmarkTableDefinition(false)` disables display clipping, and the benchmark
plot overload also defaults to `clip=false`. Absolute RMS remains available when
relative RMS cannot be normalized, such as a numerical-zero G reference. Tables
include all recorded quantities, statuses, tolerances, actual sample indices and
requested/actual frequency bounds. Stored UQ mean and standard-deviation RMS
matrices appear as separate `statistic` values; the line-parameter overlay method
does not reinterpret moment products as line-parameter results.

The engine currently selects ordinary band endpoints by the nearest stored
frequencies. Inspect `requested_bounds_Hz`, `actual_bounds_Hz` and `sample_indices`
when using a coarse grid: a band name or a 50 Hz endpoint does not establish that
a 50 Hz sample exists. Reporting displays the saved settings and values; it does not
silently recompute comparisons under another band-selection rule.

Tables load without Makie and no solver runs during inspection. The documentation
page provides a maximum-relative-error overview across every recorded quantity;
the REPL tables expose the full matrices and absolute errors.

## Run and select stored results

```bash
lcm gauntlet run --definition gauntlet/benchmarks/examples/compare_soil.jl \
  --directory /path/to/campaign
lcm gauntlet status --directory /path/to/campaign
lcm gauntlet resume --directory /path/to/campaign

lcm gauntlet compare --definition /path/to/benchmarks.toml \
  --output /path/to/comparisons
LINECABLEMODELS_GAUNTLET_RESULTS=/path/to/comparisons \
  julia --project=docs docs/make.jl
```

Catalogue cases and benchmark examples default to 101 logarithmically spaced
frequencies from 0.1 Hz to 10 MHz: `10.0 .^ range(-1, 7; length=101)`.
These are 100 increments including both endpoints, accepted by the PSCAD scan.
Explicit frequency overrides remain authoritative. Pass `--frequencies FILE.toml` with a
`frequencies = [...]` vector, or declare a grid using
`--grid log --count 101 --bounds 0.1,1e7`. The runner passes those values to
the benchmark constructor before materialization; it does not replace them
with a backend-specific grid.

A native base frequency such as 50 Hz controls physical conversion; it does
not insert a comparison sample. Reports retain the actual selected band bounds;
the engine's current endpoint selection uses the nearest stored frequencies.
Retained archives keep the frequency vectors on which they were computed.

Formula grids, material selections and uncertainty runs are documented in the
[Gauntlet CLI guide](https://github.com/Electa-Git/LineCableModels.jl/blob/main/gauntlet/README.md).

A benchmark definition names exactly one reference and one candidate artifact
(paths and SHA-256 checksums), plus quantities, bands and normalizations. See the
CLI guide for the TOML format. Comparing saved files never reruns their solvers.
Missing operands are errors; they never trigger a replacement reference.
Multiple comparison directories can be selected with the platform path-list
separator (`:` on Unix, `;` on Windows); identical benchmark records appear once.
Counts, frequency ranges and timing scopes come from the selected records.
No simulations are resumed by this page. Batch elapsed-at-completion values
are not per-selection cold or warmed execution timings.

## Numerical references for CI

Stored Gauntlet results are not automatically approved CI references. The
separate [numerical-reference gate](https://github.com/Electa-Git/LineCableModels.jl/tree/main/test/numerical)
requires reviewed results, explicit tolerances and artifact bindings. It never
starts a Gauntlet campaign or refreshes a reference. Approval remains separate
from collecting results or displaying this summary.

## Benchmark data API

```@docs
LineCableModels.Engine.RMSError
LineCableModels.Engine.LineParametersBenchmark
LineCableModels.Engine.compare
LineCableModels.Engine.absolute_error
LineCableModels.Engine.relative_error
```
