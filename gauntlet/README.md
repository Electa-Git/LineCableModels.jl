# Gauntlet: declarations, calculations and retained comparisons

Gauntlet is the application in `gauntlet/`. It ingests cases, materializes benchmark
definitions, schedules calculations, compares saved operands and locks artifacts.
`test/gauntlet/` contains its tests. `LineCableModels.PSCAD` owns the PSCAD adapter
in `ext/LineCableModelsPSCADExt/`; it loads with LineCableModels and launches the
external tool only through explicit computation or station identification.

## Declare the calculations

All catalogue cases and benchmark examples default to 101 logarithmically spaced
frequencies from 0.1 Hz to 10 MHz (100 increments, including both endpoints):

```julia
10.0 .^ range(-1, 7; length = 101)
```

Explicit `frequencies` overrides remain authoritative. The defaults contain no
extra 50 Hz sample. Saved comparisons retain their actual selected samples;
inspect requested and actual band bounds when working on a coarse grid.
Existing retained results keep their original frequency axes.

```julia
using LineCableModels
include(joinpath(pkgdir(LineCableModels), "gauntlet", "Gauntlet.jl"))
using .Gauntlet

model = load_case(:two_insulated_wires;
    variation=ExactOverrides(frequencies=[1.0, 10.0, 37.0, 1e3], temperature=60.0))
problem = model.problem # nominal_problem retains the original case baseline.
reference = BenchmarkCalculation(:default, problem, Formulation())
candidate = BenchmarkCalculation(:pollaczek, problem,
    Formulation(earth_impedance=formula(:Pollaczek1926)))
definition = benchmark_definition(:soil_comparison, model.id, :manual,
    @__FILE__, model, reference, candidate,
    (; quantities=(:Z,:Y,:R,:L,:G,:C)), (;))
result = run_benchmark(definition; directory="/path/to/calculation")
```

Use a real declaration file for `source_file` (the `@__FILE__` argument above).
For a REPL, pass the path of the file that declares the benchmark. The two
`BenchmarkCalculation` objects carry their own problem, formulation and options.
No backend owner enum, reference-grid conversion, formula substitution or
automatic reduction is applied by the runner.

Comparison settings are a named tuple. The constructor fills omitted defaults
once; `validate(definition)` checks them before either operand runs. For example:

```julia
settings = (quantities=(:Z, :G), bands=(:all, (10.0, 1e5)),
    normalizations=(:reference_rms,), atol=(G=1e-12,))
```

`run_benchmark(...).comparison` contains one row per requested quantity, band,
normalization and point pair. Each row carries its `RMSError` in `error`. Saving
uses those same comparisons. `compare_saved` explicitly requests another analysis
of completed operands. Reading, tabulating and plotting those benchmark records never recalculate RMS.

Moment comparisons select `(quantities=(:R, :L, :C, :G), statistics=(:mean, :std))`.
Their current implementation supports the full frequency range and reference-RMS
normalization; other moment settings fail before calculation. The line-parameter
statistics default is `(:value,)`.

A comparison reference defines direction and normalization. PSCAD, FEM and LCM
remain different models. Cross-model RMS values are observations. Incomparable
coordinates are rejected; numerical-zero reference traces retain absolute RMS
and report unavailable relative RMS with the recorded reason and tolerance.

Catalogue declarations are ordinary constructors under `gauntlet/benchmarks/`:

```julia
definition = benchmark_definition(:benchmark_two_insulated_wires_pscad;
    frequencies=10.0.^range(-1,7;length=101),
    reference_options=(remote=station,work_root=station.local_root))
```

Construction does not compute or discover station configuration. Frequency
overrides use the same explicit case-variation mechanism as other parameters.
Existing `Formulation(...=Grid(...))`, product/zip composition, `ParametricProblem`,
`LinearError` and `MonteCarlo` supply calculation choices. Higher-order execution
options remain at their existing formulation owner; unsupported options fail at
the public compute method. Performance repetition occurs only when explicitly
declared and records compute wall time, samples and whether native work was reused.

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

## Standalone PSCAD

```julia
using LineCableModels
station = PSCAD.RemoteConfig("pscad-host",
    raw"Z:\shared\lcm", raw"C:\LCM\scratch",
    raw"C:\Julia\bin\julia.exe", raw"C:\Python\python.exe";
    local_root="/path/to/shared/lcm")
selected = Formulation(:pscad; options=(base_frequency=60.0,))
result = compute(problem, selected;
    options=(remote=station, work_root=joinpath(station.local_root,"my-run"),
        verbosity=(default=0,PSCAD=2)))
```

`local_root` and `shared_root` name the same exchanged files from the caller and
station. `remote_root` is station scratch. `work_root` must be within `local_root`
and defaults to it. Directory nesting and system identifiers do not define a
Gauntlet-specific filesystem grammar. `gauntlet/local.example` illustrates an
explicit CLI configuration; the backend never searches for or loads that file.

The supported PSCAD physical model has air and one infinite horizontal soil
half-space. Indexed formula validation covers air, earth and mixed self/mutual
interactions. EHEM, unavailable material laws, unsupported formula hooks and
incompatible native switches fail before native execution. Conductor temperature
correction uses the selected constitutive law. `base_frequency` defaults to 50 Hz
and controls the native physical conversion; it is separate from the requested
frequency vector. Requested/exported dielectric loss tangents and the cap of ten
are retained in `details(result)`.

The minimum calculation frequency is 0.1 Hz. The verified phase-scan route accepts
101, 201, 501 or 1001 logarithmic samples, retaining the original request and raw
native frequencies. Roundoff-equivalent log grids are accepted. Other grids are
an explicit adapter limitation; they are never clamped, interpolated or replaced.
The native phase scan exposes [100, 200, 500 or 1000 increments](https://www.pscad.com/webhelp-pscad-v5.1.0-ol/Master_Library_Models/Transmission_Lines_Cables/Distributed_Line_Models/fd_phase_options.htm).
Passivity sampling is a different control and is not a Z/Y frequency-grid selector.

`details(result).files` enumerates retained backend evidence by relative name,
original file path and SHA-256. Gauntlet copies those files into the calculation
bundle; original workstation paths remain calculation records. FEM supplies this evidence
when its declaration retains its native run directory.

## CLI and recovery

```sh
cli/lcm gauntlet case list
cli/lcm gauntlet case import --id my_case --source case.jl --project . --dry-run
cli/lcm gauntlet run --definition benchmarks.jl --directory runs/example
cli/lcm gauntlet run --definition gauntlet/benchmarks/pscad/benchmark_two_insulated_wires_pscad.jl \
  --configuration gauntlet/local.jl --grid log --count 101 --bounds 0.1,1e7 --directory runs/pscad
cli/lcm gauntlet status --directory runs/example
cli/lcm gauntlet resume --directory runs/example
cli/lcm gauntlet lock --directory runs/example --output retained/example
```

A declaration file returns one definition, a vector, or a constructor accepting
explicit keywords. A configuration file returns those keywords as a named tuple.
`--frequencies frequencies.toml` accepts `frequencies=[...]`; alternatively a
constructed grid requires `--grid`, `--count` and `--bounds` together. Omitting
these flags preserves the definition's frequencies.

Campaign execution passes the materialized calculations to `compute`. It stores
atomic completion records and checksummed result payloads. Reuse checks actual
numerical inputs, formula options and captured runtime sources. Dirty source bytes
and dependency environments are retained, not just a Git revision. Failed operands
retain diagnostics; resume preserves completed calculations. Changing RMS bands
or tolerances creates another analysis over those retained operands.

PSCAD separately owns native completion, installation identity, raw-output checks,
readbacks and `resume_run_directory`. A campaign-level reused result is historical
completed evidence, not a new native timing sample. To require a fresh calculation,
declare a new campaign and omit native resume. Use explicit performance declarations
for repeated solves.

`lock_campaign` copies the completed campaign into an immutable bundle. Numerical
operands, analysis bindings, raw files and source evidence are checksummed. Read it
with `read_campaign` after moving it to another directory; no original work directory
or native solver is needed. Changes require a new bundle identity. Locking does not
publish anything or certify one numerical model as truth.

## Navigation and checks

- `definitions.jl`, `benchmarks/`: concrete benchmark composition and catalogue.
- `cases.jl`, `cases/`: model inputs, variations and source assets.
- `benchmarks.jl`, `campaigns.jl`: execution, scheduling and recovery.
- `records.jl`, `fingerprints.jl`: source capture and numerical identities.
- `comparisons/`: saved-operand binding and the package comparison contract.
- `artifacts.jl`, `read.jl`: packaging, immutable bundles and retained records.
- `docs/gauntlet_report.jl`: rendering of saved comparison objects.

```sh
julia --project=gauntlet -e 'using Pkg; Pkg.instantiate()'
julia --project=gauntlet test/gauntlet/runtests.jl
julia --project=test test/runtests.jl pscad
```

Ordinary tests run without coverage instrumentation. Coverage collection retains
the LCOV report and removes source-adjacent `.cov` traces after workers have exited.


## Historical records

Portable schema-1 calculation files remain readable with `read_calculation`; the
reader verifies their original checksums and does not change published hashes.
Old collection snapshots also retain their finite schema-2 reader. They are
historical observations, not scientific acceptance gates.

For typed `problem`/`result` JLD2 payloads from the checkpoint before extraction,
use `migrate_typed.jl` once in a detached checkout of
`2d694a2444f52942e9d14cd5ff835760295ff532` and its `test/gauntlet` environment:

```sh
julia --project=/path/to/pinned/test/gauntlet gauntlet/migrate_typed.jl \
  old-result.jld2 new/calculation.jld2 declaration.toml
```

The declaration explicitly supplies `case_id`, `backend` and `port_order` in
stored matrix order; `result_key` and `problem_key` default to `result` and
`problem`. The bridge requires the pinned implementation to be unchanged. It
copies the original payload, records its hash and the old source bytes, and
exports plain numerical fields. It runs no solver and modifies no original file.
The new reader needs neither the old module namespace nor its checkout.

Fresh-process campaign resume loads the recorded Julia dependencies and verifies
the case builder sources before reading typed execution declarations. The CLI
also records and reloads the explicit declaration/configuration source files.
Changed numerical source bytes reject execution reuse; reading a completed bundle
and making new comparisons remain independent of those execution dependencies.

All built-in cases can be checked without a solver using
`cli/lcm gauntlet case validate`. A benchmark factory may accept explicit
frequencies, station options or variations. There is no implicit backend catalogue
that replaces a declared formulation or suppresses an unavailable selection.
