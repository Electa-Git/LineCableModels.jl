# Gauntlet: declarations, calculations and retained comparisons

Gauntlet is the application in `gauntlet/`. It ingests cases, materializes benchmark
definitions, schedules calculations, compares saved operands and locks artifacts.
`test/gauntlet/` contains its tests. `LineCableModels.PSCAD` owns the PSCAD adapter
in `ext/LineCableModelsPSCADExt/`; it loads with LineCableModels and launches the
external tool only through explicit computation or station identification.

Use `lcm gauntlet ...` through the [package CLI](../cli/README.md), or
`./cli/lcm gauntlet ...` from the repository root before global installation.
`lcm --paths` identifies the checkout supplying each application. Gauntlet owns
its Julia environment and command arguments; global routing starts no solver.

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
using LineCableModels.ReportBuilder: BenchmarkTableDefinition

model = load_case(:two_insulated_wires;
    variation=ExactOverrides(frequencies=[1.0, 10.0, 37.0, 1e3], temperature=60.0))
problem = model.problem # nominal_problem retains the original case baseline.
reference = Formulation()
formulations = Formulation(earth_impedance=Grid((:default, :Pollaczek1926)))
definition = benchmark_definition(model; id=:soil_comparison,
    source_file=@__FILE__, reference, formulations, collection=:manual,
    report=BenchmarkTableDefinition())
result = run_benchmark(definition; directory="/path/to/calculation")
```

Spectral controls remain on each formula selection, including when it is an
entry in a formulation Gridspace. For trapz and CIM, `samples=nothing` selects
adaptive construction; an integer caps construction kernel evaluations per scalar
integral. These are separate from the problem's frequency samples and timing
repetitions. Saving and reporting retain the requested budget and the normalized
per-interaction controls. An exhausted budget fails explicitly.

Use a real declaration file for `source_file` (the `@__FILE__` argument above).
For a REPL, pass the path of the file that declares the benchmark. An explicit
`BenchmarkCalculation` carries its own problem, formulation and execution options;
it can also be passed as `reference` to the convenience constructor.
No backend owner enum, reference-grid conversion, formula substitution or
automatic reduction is applied by the runner.

`BenchmarkTableDefinition` owns comparison settings. The constructor fills omitted defaults
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

The shared report works with live results independently of Gauntlet:

```julia
using LineCableModels
using LineCableModels.ReportBuilder: BenchmarkTableDefinition

selected = Formulation(earth_impedance=Grid((:default, :Pollaczek1926)))
reference = compute(problem, Formulation())
candidates = compute(problem, selected)
artifact = report(BenchmarkTableDefinition(), (; reference, candidate=candidates))
artifact.table.summary
artifact.table.formulations
artifact.table.maxima
artifact.table.terms

using GLMakie
plots = LineCableModels.plot(artifact; ydata=(Z, Y))
# Directly from completed results:
plots = LineCableModels.plot(candidates; ydata=(Z, Y), reference)
```

A scalar reference is compared once with each candidate. Problem and formulation
axes retain the ordinary `results[problem, formulation]` ordering. Two result
spaces require explicit `pairing=[(reference_index, candidate_index), ...]`, with
every candidate named once. Equal lengths do not establish correspondence.
Native results carry actual output terminal identities, including reductions.
Scalar results also retain complete requested and resolved formula records in
`details(result).formulations`: `requested` contains the supplied parameters and
hooks, and `methods` contains the resolved owner records. This is the same
`NamedTuple(formulation)` representation used by formulation-space reports.
Changing a physical reference or hook therefore remains identifiable after
saving and reopening a scalar result.
External results without terminal identities require explicit context:

```julia
operand = (result=external_result,
    metadata=(port_order=["cable:1:core", "cable:2:core"],
        formulation=(backend=:external, requested=(solver=:chosen,)), axes=nothing))
artifact = report(BenchmarkTableDefinition(), (reference=operand, candidate=candidates))
```

Read saved results using the same report:

```julia
benchmark = read_benchmark("/path/to/staging/benchmark_id")
# One particular analysis:
# benchmark = read_benchmark("/path/to/snapshot.jld2"; load_results=true)
artifact = report(BenchmarkTableDefinition(), benchmark)
tables = artifact.table

tables.summary      # Compact maxima, grouped by problem/formulation and band
tables.maxima       # Numerical per-term maxima, units, coordinates and availability
tables.formulations # Stable keys, labels and complete selection records
tables.calculations # Inputs, output coordinates, execution records and file hashes
tables.comparisons  # Full RMS matrices per quantity, band and formulation
tables.terms        # Every matrix entry, absolute/relative RMS, status and reason

z = filter(row -> row.quantity === :Z && row.band === :all &&
    row.formulation_index == 2, tables.comparisons)
only(z.absolute_rms)
only(z.relative_rms_percent)
only(z.reason)

# Select products already retained; no comparison is recalculated.
report(BenchmarkTableDefinition(quantities=(G,), bands=(:all,)), benchmark)

using GLMakie
plots = LineCableModels.plot(benchmark; ydata=(Z, Y))
selected_plots = LineCableModels.plot(artifact; ydata=(R, L), formulations=[2])
export_svg(first(plots); path="impedance.svg", open_file=false)
```

Explicitly requesting one problem overlays its reference and every candidate in
each matrix cell, including both off-diagonals. `Z` gives R/X pages and `Y` gives
G/B pages. `(R, L, G, C)` requests those quantities directly. Multiple problems
require `problem=2` or an explicit list; each receives its own pages. Formula
filtering preserves original indices and colors. Equal curves remain separate
selections. Scalar and air/earth/mixed choices, parameters, hooks and numerical
options remain in the full records. Legends show differing fields; common choices
remain in `tables.formulations.record`.

The existing unit controls, legends, zoom and SVG export remain available.
Benchmark plots also retain the log-y toggle for signed matrix entries, using a
sign-preserving pseudo-log scale on those panels. Use
`length_unit=:base` to display native per-metre quantities. `band=:dc` uses saved
sample indices. A band or numerical setting that was not retained requires explicit
reanalysis of the raw operands, through `report` or `compare_saved`. Loading,
tabulating, filtering and plotting a retained analysis never calculate RMS again.

ReportBuilder owns the default order: entire range, near DC (0.1–100 Hz), harmonic
range (50–2500 Hz by default), narrowband (1 kHz–1 MHz), wideband (>1 MHz). Groups
overlap. Ordinary endpoints use the engine's existing nearest-sample selection.
Tables retain requested/actual bounds and sample indices; no interpolation occurs.
`clip=false` is the default. Numerical-zero G retains absolute RMS and unavailable
relative RMS with its reason. UQ mean/std products remain separate statistics under
their existing full-band contract.

`artifact.published` holds the unformatted scientific products. Ordinary display
shows the compact summary. No figure is constructed unless requested by `plot`
or `BenchmarkTableDefinition(illustration=true, plot_options=(...))`.

## What the documentation publishes

The standard Gauntlet page contains completion counts, compact formulation keys
and five ordered band summaries. Each quantity cell shows the maximum of per-term
RMS discrepancies and its terminal pair, with unavailable-term counts. These are
neither whole-matrix RMS values nor averages across formulations. When relative
RMS is unavailable, the page shows absolute RMS with units and the reason.

The page has no default matrix figures, thumbnails or interactive plot payloads.
A build reads retained summaries and explicitly selected image files; it does not
solve, calculate RMS or render benchmark illustrations. Detailed plots remain an
explicit REPL request for a chosen benchmark/problem.

`docs/gauntlet.toml` lists version-specific artifact bindings, for example
`artifacts = ["gauntlet_soil_v1_0_0"]`. An empty list reports that no published
results are selected. For a local preview, explicitly set
`LINECABLEMODELS_GAUNTLET_RESULTS=/path/to/staging` (or a vault/package artifact
directory). The page labels that source as a local preview. It never selects the
latest staging directory implicitly.

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

## CLI, drafts and acceptance

```sh
cli/lcm gauntlet case list
cli/lcm gauntlet case import --id my_case --source case.jl --project . --dry-run
cli/lcm gauntlet run --definition benchmarks.jl --directory runs/example
cli/lcm gauntlet status --directory runs/example
cli/lcm gauntlet resume --directory runs/example
cli/lcm gauntlet lock --directory runs/example --benchmark soil_comparison \
  --expected INSPECTED_SNAPSHOT --output vault/soil-reviewed --note "Inspected at the REPL"
```

A declaration file returns one definition, a vector, or a constructor accepting
explicit keywords. A configuration file returns those keywords as a named tuple.
`--frequencies frequencies.toml` accepts `frequencies=[...]`; a generated grid
requires `--grid`, `--count` and `--bounds` together. Omitted flags preserve the
definition's frequencies. Backend constraints remain the backend's responsibility.

`run` starts a new attempt for each selected benchmark. On success it replaces
that benchmark's previous draft. Unselected drafts and vault bundles remain intact.
If the replacement fails, the previous complete result remains accessible with
`read_benchmark(path; previous=true)`; ordinary reads report the failed latest
attempt. `status` reports this distinction and the current inspected identity.
Every selected declaration is saved before the first solver starts. `resume`
checks numerical inputs and saved-file integrity, and reuses matching completed
calculations. It does not compare the live source tree. A scalar problem with a
formulation `Gridspace` checkpoints each completed point, so an interruption need
not repeat earlier points.
The low-level `run_benchmark` retains and
reuses calculations in an explicitly supplied attempt directory; use `run_campaign`
or the CLI for replaceable drafts.

Each invocation records its Julia version, loaded package versions, Git revision
and dirty flag once. Reused operands retain their original session metadata; a
resumed campaign can therefore contain results from different execution sessions.
This records provenance, not a guarantee of a frozen or reproducible environment.
Use a pinned checkout/environment when that guarantee is required.

Numerical operands are saved before reporting; reports are saved before optional
timing repetitions. A reporting or timing failure leaves completed work available.
Changing report bands requires reanalysis, not another solver run.

`resume --recover-solvers` also asks FEM and PSCAD to recover compatible native
run directories, after the backend checks its actual inputs and saved outputs.
For older campaigns with unsaved queued declarations, supply the original factory
and options with `run ... --resume --recover-solvers`. Existing attempts must have
matching numerical declarations; missing queue entries are saved before execution.

PSCAD independently owns native reuse and installation/readback verification.
A campaign retry is not a new cold timing sample. Explicit performance definitions
own repeated timing measurements. Reanalysis over saved operands creates another
checksummed analysis, with its own snapshot identity, without running a backend.

`lock` is the acceptance action. It copies only the selected complete benchmarks
into a self-contained vault directory. Omit `--benchmark` to accept the declared
whole campaign; incomplete members then reject the operation. `--expected` checks
that a single selected benchmark is still the snapshot inspected by the user.
Locking neither applies an error threshold nor claims agreement with another model.

A vault retains calculations, all saved analyses, captured declarations and session metadata,
backend files and relative operand bindings. It remains readable after staging is
moved or deleted. Missing, changed and unexpected files reject reads. Writers reject
vault destinations; corrections require a new bundle. Locking does not upload,
change Git history or update published artifact bindings.

## Live progress and performance

Campaigns display progress on stderr by default. `--progress auto` uses an updating
bar in a terminal and throttled plain output when redirected; `plain` always emits
plain output, and `off` disables monitoring and its disposable snapshots. REPL
calls accept the same choices as `progress=:auto`, `:plain` or `:off`.

```bash
./gauntlet/lcm gauntlet run --definition gauntlet/benchmark_all_references.jl \
  --directory gauntlet/.work/all-references --configuration gauntlet/local.jl \
  --grid log --count 101 --bounds 0.1,1e7 --on-error continue --progress auto
./gauntlet/lcm gauntlet status --directory gauntlet/.work/all-references --watch
```

The first command starts fresh attempts. Add `--resume --recover-solvers` only when
you want verified completed work reused. PSCAD excludes the bare-wire case. The
local configuration supplies your station and diagnostic verbosity independently
of progress; set `verbosity=(default=0, PSCAD=0)` there for quiet native chatter.

The display separates selected benchmarks, calculation jobs, native worker jobs
and MC accepted trials/attempts/rejections. Failed, skipped and interrupted work
are not counted as successful. It shows outer point/formulation and backend stages,
not the frequency currently inside a solver. Warnings remain visible. Explicit
Gmsh console diagnostics switch the bar to append-only plain output. `status
--watch` only reads metadata and active-session snapshots; Ctrl-C stops the watcher,
not the campaign. Old runs without snapshots still show coarse campaign status.

Each reference/candidate completion leaves a persistent summary above the live
display, using the same labels for Owned, FEM, PSCAD, Monte Carlo and LEP runs.
It shows execution wall time (through required result persistence), compute-call
wall time, completed/reused calculation jobs and whether recovery was involved.
Compute-call time covers calls actually made in this invocation; a saved operand
instead says `not run (saved result)`, without displaying its historical compute
time as a new measurement. Failed/interrupted operands show their elapsed wall
time and leave unavailable measurements explicit. These are ordinary execution
observations, not controlled performance samples. Summaries are emitted by the
renderer, never from numerical loops, and are disabled by `--progress off`.
Routine PSCAD messages honor the same 0/1/2 verbosity choices as FEM. With explicit
verbosity, PSCAD's completion diagnostic labels the compile-call duration and its
scope; it does not advertise that duration as elapsed scan time.

Stage ETAs use bounded observed throughput. FEM estimates use MSH 4.1 ASCII node
counts as a workload proxy, calibrated from at least three completed workers;
remaining terminal columns, worker concurrency and the final serial tail are
included. Unsupported mesh formats retain the identical-mesh throughput fallback.
Campaign estimates learn from completed calculations during the current invocation,
grouped by backend, propagation method and execution settings, with frequency,
matrix and trial counts as workload proxies. Exact compatible history takes
precedence. PSCAD uses comparable whole-call durations, including native setup.
Reused/recovered calculations do not teach fresh computation costs.
Finalization has its own observations, including reporting, persistence and any
requested performance pass. Approximate estimates use `~`; when some workloads
are unknown, the display reports the estimated portion separately. An untouched
backend remains unknown; PSCAD observations cannot predict FEM or Monte Carlo.
Overdue work returns to `estimating` until new evidence arrives. Heartbeats establish
recent observation, not numerical progress. Estimates are advisory and never feed
performance comparisons or numerical reuse decisions.

Timing records distinguish execution wall, compute-call wall, and backend-native
time. `sessions/SESSION.timing.toml` retains invocation wall time even with progress
off. `timing.toml` retains operand wall/compute durations; calculation metadata
retains compute policy and source timings. GetDP phase sums are accumulated worker
time, **not** elapsed scan time. PSCAD source time covers remote `line.compile()`;
output readiness and transfer belong to broader adapter time. Historical elapsed
fields keep their legacy meaning. Verified recovery retains its original timings
and is not a new zero-second solver measurement.

Declared performance checks warm owned Julia calls, then time complete compute
calls without progress, optional diagnostics or `on_result` callbacks. Required
numerical/checkpoint IO remains included. MC samples include draws, reconstruction,
retries and aggregation. The display says `Performance sample … — display paused`:
there are no live trial updates, redraws or snapshot writes inside that sample.
Normal execution still delivers user callbacks. Julia allocation counts cover
Julia, not native worker memory. Native repetitions are not enabled implicitly.

Operational runs include monitoring overhead and system load. Quiet samples remove
the observer's active work, not unrelated machine contention. Speedup decisions
require matching timing scope, measurement policy and environment, with no recovered
native work or coverage/allocation instrumentation. One sample is limited evidence.

## Package and publish accepted bundles

Create an explicit release definition, with bundle paths relative to this file:

```toml
# release.toml
collection = "soil"
version = "1.0.0"
description = "Reviewed soil formulation comparisons"
bundles = ["vault/soil-reviewed"]
```

```sh
cli/lcm gauntlet package --definition release.toml --output packages/soil-1.0.0
# Upload the exact generated .tar.gz using the chosen host's existing tool.
cli/lcm gauntlet bind --package packages/soil-1.0.0 \
  --url https://your-host.example/soil-1.0.0.tar.gz
```

Packaging consumes locked bundles only. It verifies their contents, rejects
conflicting case/benchmark entries, and produces a standalone Julia artifact
archive with its tree hash and archive SHA-256. Repeating the same release into
the same destination verifies and reuses it. A different definition requires a
new version/destination; there is no `--force` replacement. No calculation,
comparison, plot, upload or binding update occurs during packaging.

`bind` downloads the served archive, verifies both its bytes and extracted artifact
tree, then writes a version-specific binding. Failed verification leaves the
binding unchanged. A published version cannot be rebound to different bytes.
`--current` explicitly updates a convenience binding; documentation pins the
version-specific name. `--artifacts FILE.toml` selects a separate binding file.

Static illustrations are opt-in. Export the selected plot through PlotBuilder,
then pass `lock --illustrations illustrations.toml`. Every entry declares the
accepted benchmark, existing file, caption and display selection:

```toml
[[illustrations]]
benchmark = "soil_comparison"
path = "figures/resistance.svg"
caption = "Self and mutual resistance for the selected soil formulations"
[illustrations.selection]
problem = 1
quantities = ["R"]
formulations = [1, 2]
length_unit = "base"
```

The lock records the image bytes and selection. Packaging copies them unchanged;
the documentation embeds only those retained files. Neither operation generates
additional figures.

## Navigation and checks

- `definitions.jl`, `benchmarks/`: concrete benchmark composition and catalogue.
- `cases.jl`, `cases/`: model inputs, variations and source assets.
- `benchmarks.jl`, `campaigns.jl`: execution, scheduling and recovery.
- `records.jl`, `fingerprints.jl`: source capture and numerical identities.
- `comparisons/`: saved-operand binding and the package comparison contract.
- `artifacts.jl`, `read.jl`: packaging, immutable bundles and retained records.
- `docs/gauntlet_report.jl`: summary-only publication of retained report products.

```sh
julia --project=gauntlet -e 'using Pkg; Pkg.instantiate()'
julia --project=gauntlet test/gauntlet/runtests.jl
julia --project=test test/runtests.jl pscad
```

Pass file or test-name fragments to select only the affected Gauntlet tests:

```sh
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --project=gauntlet \
  --startup-file=no --compiled-modules=existing \
  test/gauntlet/runtests.jl spectral_reporting_tests
```

Multiple selectors are alternatives. Omitting selectors retains the full suite.
For upstream spectral-control reconciliation, the selected regression covers
formula Gridspace computation, saved results, reporting and budget exhaustion.
The wider numerical, FEM and documentation suites can run in CI when local power
or time is limited.

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

Fresh-process resume loads the recorded Julia dependencies from the current
environment and restores captured case/declaration/configuration code before
reading typed execution declarations. Original declaration files may be edited or
deleted. Relative case assets declared in `assets` are restored beside the case;
configuration files should be self-contained or explicitly capture their dependencies.
Opaque execution objects use Julia's native serialization inside the checkpoint;
portable numerical fields and source evidence remain independently readable JLD2 data.
Resume does not reinstall old package versions. Reading completed numerical bundles
and making new comparisons remain independent of the original execution sources.

All built-in cases can be checked without a solver using
`cli/lcm gauntlet case validate`. A benchmark factory may accept explicit
frequencies, station options or variations. There is no implicit backend catalogue
that replaces a declared formulation or suppresses an unavailable selection.
