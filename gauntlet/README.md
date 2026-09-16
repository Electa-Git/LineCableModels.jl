# Gauntlet: declarations, calculations and retained comparisons

Gauntlet is the application in `gauntlet/`. It ingests cases, materializes benchmark
definitions, schedules calculations, compares saved operands and locks artifacts.
Its role is **run, compare and report**. Numerical differences, unavailable relative
RMS values and timing ratios are observations for the user to interpret. Gauntlet
assigns no scientific or performance acceptance verdicts.
`test/gauntlet/` contains its tests. `LineCableModels.PSCAD` owns the PSCAD adapter
in `ext/LineCableModelsPSCADExt/`; it loads with LineCableModels and launches the
external tool only through explicit computation or station identification.

Use `lcm gauntlet ...` through the [package CLI](../cli/README.md), or
`./cli/lcm gauntlet ...` from the repository root before global installation.
`lcm --paths` identifies the checkout supplying each application. Gauntlet owns
its Julia environment and command arguments; global routing starts no solver.

## Declare the calculations

All catalog cases and benchmark examples default to 101 logarithmically spaced
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
formulations = Formulation(earth_impedance=Grid((:default, :pollaczek1926)))
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

Statistical comparisons select ordinary requests such as
`BenchmarkTableDefinition(((statistics, R, mean), (statistics, R, std)); bands=(:all, :dc, :harmonic, :narrow, :wide))`.
The `quantities`/`statistics` shorthand remains accepted and is expanded once.
MC and LEP retain their native results; each outer point remains a separate population.
Mean and standard-deviation comparisons use the same physical cutoffs, bands and
normalizations as deterministic quantities. Retained MC percentiles remain available
without inventing a LEP output distribution.

`artifact.table.features` contains numeric formula-by-band absolute and relative
DataFrames. `statistics` and `sampling` expose retained UQ summaries and sampling
precision. `execution`, `source_timings`, `performance`, `performance_samples`,
`performance_environment` and `performance_policy` describe recorded measurements;
they never collect new timings. Worker-time sums and PSCAD compile-only times retain
their original scopes. A reference/candidate timing ratio describes those workloads,
not equal accuracy or cost per MC trial.

Read a saved benchmark with `read_benchmark(path; load_results=true)` and render its
historical analysis with `report(BenchmarkTableDefinition(), benchmark)`. To reanalyse,
pass its `reference` and `candidate` operands explicitly to a new report definition.
Neither operation solves a model. `evidence=:numerical` permits numerical inspection
when auxiliary solver files changed; missing evidence is reported and strict
acceptance/locking is not relaxed.

A comparison reference defines direction and normalization. PSCAD, FEM and LCM
remain different models. Cross-model RMS values are observations. Incomparable
coordinates are rejected; numerical-zero reference traces retain absolute RMS
and report unavailable relative RMS with the recorded reason and tolerance.
That numerical-resolution tolerance decides whether a ratio is meaningful; it is
not an acceptance limit or an absolute-error fallback verdict. UQ declarations
request timing samples with `uq_timing_settings()`; they set no moment-error or
minimum-speedup thresholds. `performance_comparison(reference, current)` reports
ratios and whether the measurement environments/scopes are comparable.

Historical declarations retain their serialized `tolerances` field so saved work
orders remain readable. Only its performance sampling budget is used; old numerical
limits and minimum-speedup thresholds have no effect on execution or reporting.

Catalog declarations are ordinary constructors under `gauntlet/benchmarks/`:

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

## Catalog uncertainty law

New catalog LEP/Monte Carlo declarations use
`bounded_correlated_geometry_v1`, a **synthetic correlated-geometry study**:

- All selected cross-section lengths share one uniform scale of mean 1 and
  standard deviation 0.1, bounded by `1 ± sqrt(3)*0.1`.
- Lay ratios and the Milliken fillet factor have separate independent uniform
  inputs with 10% standard uncertainty. The explicit solid-core area is an
  independent area input, not a length.
- Counts, formation controls and unselected properties remain exact unless
  explicitly varied. The existing reuse of a design across cable positions
  retains its shared inputs.

Role tags `:length`, `:area` and `:dimensionless` make that declaration explicit.
Finite exact/grid overrides become base values before applying the law; an
already uncertain geometric source conflicts with the default law. Use a
complete custom joint declaration instead of layering incompatible assumptions.
Nominal/base construction is checked before execution, and the Milliken case
retains a conservative contact/inventory support certificate.

The new law preserves the specified marginal length means and standard
deviations, **not** independence or the old normal distribution. It is not a
claim about measured manufacturing correlations. A derived area scales as the
square of the common scale and has mean `1.01A₀`; first-order LEP does not include
that second-order shift. Moment differences are not removed by relaxing tolerances.

`JointParameterGrids(source; record)` is an application ingestion adapter for an
ordinary `Gridspace{NamedTuple{names}}`. Its record describes independent factor
inputs, supports, base points, roles and output dependencies. The executable
source is passed once to the case builder; `LoadedCase.sources` projections are
only marginal views. Later independent overrides of joint-owned fields fail.
Neither the runner nor the solver implements a separate correlated sampler.

Use a **fresh `run` declaration**, without `--resume`, to adopt this law. Resume
continues the saved declaration and must not relabel an old independent-normal
study. Existing attempts/results and solver caches need not be deleted.

The 320/380 kV armored cases now explicitly declare their previously resolved
9.376520171 mm bedding (a reconstructed benchmark dimension, not a verified
manufacturer value). The 525 kV/1600 mm² case removes its artificial buffer:
bedding is 3 mm and outer radius is 79.04885 mm instead of 94.85862 mm. That case
requires matching new FEM/PSCAD references. A nonzero
`armor_packing_clearance_ratio` remains an explicit design override, never an
automatic repair. Other PSCAD references remain candidates for native recovery
when actual exported inputs and solver settings match; source edits alone do
not establish compatibility.

## Inspect results at the REPL

The shared report works with live results independently of Gauntlet:

```julia
using LineCableModels
using LineCableModels.ReportBuilder: BenchmarkTableDefinition

selected = Formulation(earth_impedance=Grid((:default, :pollaczek1926)))
reference = compute(problem, Formulation())
candidates = compute(problem, selected)
artifact = report(BenchmarkTableDefinition(), (; reference, candidate=candidates))
display(artifact)
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
`details(result).data.formulations`: `requested` contains the supplied parameters and
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

tables.features     # Per-quantity/statistic DataFrames: unique formulas × bands
tables.overview     # Compact coverage, timing and MC sampling tables
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
filtering selects stored array positions and preserves source colors. Candidates
with the same quantity-relevant selection share one curve and feature-table row;
different formulas remain separate even when their curves agree. Scalar and
air/earth/mixed choices, parameters, hooks and numerical
options remain in the full published operand metadata. `tables.formulations`
contains readable labels; `tables.formula_details` contains owner-dispatched
scientific explanations. Legends use only the choices relevant to the plotted
quantity. References are named `Reference · FEM`, `Reference · PSCAD`, or
`Reference · Monte Carlo`; candidate labels have no numbering, for example `LEP`.
Every relevant route and control participates in equality; descriptions do not.
Reordering with `formulations=[3,1]` applies before deduplication. Conflicting
observations under the same selection raise an error; saved results stay intact.

For UQ, mean and standard-deviation comparisons have separate band tables.
`tables.sampling` contains scalar trial counts and CDF precision, while
`tables.mean_sampling_precision` identifies the standard error of each sampled
mean by quantity, point, terminals and frequency. Neither is a precision estimate
for the sample standard deviation. Timing/allocation tables contain scalar
measurements, original scopes and actual versus requested timed repetitions;
they do not contain whole workload/session objects.

`terms` and `maxima` expose scalar errors, method labels and terminal coordinates;
frequency bounds have separate lower/upper Hz columns. `comparisons` remains the
explicit structured audit view. Unformatted maxima are retained in native
DataFrame metadata for the existing saved-summary writer, not reconstructed
from display strings. Descriptions use retained consumed selections when present;
a declaration alone does not establish which geometry-dependent routes executed.

The existing unit controls, legends, zoom and SVG export remain available.
Benchmark plots also retain the log-y toggle for signed matrix entries, using a
sign-preserving pseudo-log scale on those panels. Use
`length_unit=:base` to display native per-meter quantities. `band=:dc` uses saved
sample indices. A band or numerical setting that was not retained requires explicit
reanalysis of the raw operands, through `report` or `compare_saved`. Loading,
tabulating, filtering and plotting a retained analysis never calculate RMS again.

ReportBuilder owns the default order: entire range, near DC (0.1–100 Hz), harmonic
range (50–2500 Hz by default), narrowband (1 kHz–1 MHz), wideband (>1 MHz). Groups
overlap. Ordinary endpoints use the engine's existing nearest-sample selection.
Tables retain requested/actual bounds and sample indices; no interpolation occurs.
`clip=false` is the default. Numerical-zero G retains absolute RMS and unavailable
relative RMS with its reason. UQ statistics use the same selectable bands and
two-sided relative eligibility, without pooling populations or matrix terms.

`artifact.published` holds the unformatted scientific products. Ordinary display
shows the compact summary. No figure is constructed unless requested by `plot`
or `BenchmarkTableDefinition(illustration=true, plot_options=(...))`.

## What the documentation publishes

The standard Gauntlet page contains completion counts and one numeric DataFrame
per physical quantity: unique relevant formulations are rows, retained bands are
adjacent columns (`all`, `dc`, `harmonic`, `narrow`, `wide`). UQ mean and standard
deviation have separate tables. Cells show maximum eligible per-term relative RMS
[%], not whole-matrix RMS or averages across formulations. Missing stays missing;
absolute errors, winning pairs and unavailable-term counts remain in detailed
`features`, `maxima` and `terms` tables. References are comparison methods, not truth.

`artifact.table.overview` owns compact frequency coverage, recorded whole-workload
timings, controlled median seconds, timed-call counts, cumulative Julia allocations
in MiB, and the recorded reference/candidate time ratio with its comparability flag.
MC trial counts, input distribution and CDF precision are separate tables: the CDF
bound is not relative error in the mean or standard deviation. Unrecorded timings
are not estimated, a single timing cannot establish variability, and a batch timing
is not a per-formula measurement. REPL and HTML displays use these same views.

The summary page contains no plots, including retained illustrations. A build reads
saved comparisons; it does not solve, sample, calculate RMS, collect timings or copy
images. Detailed plots remain an explicit REPL request for a chosen benchmark/problem.
The disposable `dev/inspect_saved_gauntlet.jl` and `dev/inspect_saved_gauntlet_mc.jl`
scripts show the shared compact summary and expose its DataFrames in the IDE, plus
optional detail tables and plots. Unlike publication, they explicitly reanalyse the
saved operands using the selected quantities and bands, without rerunning solvers.

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

`details(result).data.files` enumerates retained backend evidence by relative name,
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

`run` reconciles the supplied definitions with saved work. Identical completed
benchmarks are skipped, identical unfinished declarations continue, and changed
declarations start new attempts. Matching saved reference/candidate operands and
formulation points are reused before entering any backend. `--force` explicitly
requests new calculations for the selection, including after an interrupted forced
attempt is resumed. Historical attempts, unselected drafts and vault bundles remain intact.
When only report settings change, matching controlled timing observations also
retain their original execution session. Changing the declared timing workload
or timing settings requests new measurements.
If the replacement fails, the previous complete result remains accessible with
`read_benchmark(path; previous=true)`; ordinary reads report the failed latest
attempt. `status` reports this distinction and the current inspected identity.
Every selected declaration is saved before the first solver starts. `resume`
continues the saved work order, including its input law; it does not adopt current
catalog defaults. Completed entries skip before restoring their executable
declarations or loading numerical results. A skip verifies declaration, calculation
and report payload checksums, without traversing native evidence or regenerating
comparisons. Public saved-artifact readers and verified status retain the full
evidence checks. Neither route compares the live source tree. A scalar problem with a
formulation `Gridspace` checkpoints each completed point, so an interruption need
not repeat earlier points.
The low-level `run_benchmark` retains and
reuses calculations in an explicitly supplied attempt directory; use `run_campaign`
or the CLI for replaceable drafts.

Select exact IDs with `--benchmark ID[,ID]` on `run` or `resume`. Unknown or repeated
IDs are rejected before campaign writes. Preview per-operand scheduling decisions
with `run ... --dry-run`; it performs no campaign writes, solver calls or full
artifact integrity check. In Julia these controls are `benchmark`, `dry_run` and
`force` keywords. A skipped outcome has `state=:complete`, `skipped=true` and
`result=nothing`; load results explicitly with `read_benchmark` when needed.
Use `--force` for a numerical implementation correction whose declared inputs have
not changed. It cannot be combined with strict `run --resume`.

Each invocation records its Julia version, loaded package versions, Git revision
and dirty flag once. Reused operands retain their original session metadata; a
resumed campaign can therefore contain results from different execution sessions.
To repeat a calculation in the same environment, pin the checkout and dependency versions.

Numerical operands are saved before reporting; reports are saved before optional
timing repetitions. A reporting or timing failure leaves completed work available.
Changing report bands requires reanalysis, not another solver run.

`--recover-solvers` asks FEM and PSCAD to recover compatible native run directories
only when there is no matching Gauntlet calculation. The backend then checks its
actual inputs and saved outputs; PSCAD may inspect its remote installation.
A whole-benchmark skip or saved-operand reuse never contacts that station.
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

## Campaign progress and performance

`--progress auto` publishes lightweight snapshots and prints an exact-session watch
command plus a final summary. It installs no live display in the calculation
process. `plain` adds single-line status at outer boundaries, throttled to about ten
seconds. `off` disables optional observation, estimation, publication, and progress
output. Julia entrypoints accept the same values as `progress=:auto`, `:plain`, or
`:off`. Diagnostics retain their separately declared verbosity and callbacks.

Run the printed watcher command in a separate terminal:

```bash
./gauntlet/lcm gauntlet run --definition FILE.jl --directory DIR --progress auto
./gauntlet/lcm gauntlet resume --directory DIR --progress auto
./gauntlet/lcm gauntlet status --directory DIR --watch --session SESSION
```

The startup hint contains the actual safely quoted directory and session. The
optional session selector pins the invocation across benchmark/attempt transitions
and after completion. Directory-only watching remains valid and uses conservative
metadata discovery; ambiguous sessions produce an inventory-only display without a
combined ETA. Session identity does not establish solver liveness. Ctrl-C or a
broken output pipe ends only watching. One-shot status remains available.

At ordinary terminal sizes, the watcher owns six rows:

```text
Benchmarks [====................] 12/59 finished
Complete 11 | Failed 1 | Running 1 | Pending 46
13/59 NA2XS2Y trefoil | reference / Monte Carlo
Elapsed 00:42:18 | ETA ~00:08:40
Current run / (96/512 scans) | sampling; rejected 3
Last work observation 2s ago
```

A scan is an accepted complete frequency-grid result for one realized problem and
core formulation. MC counts accepted trials, including an inferred target once it
is known; rejections consume time without advancing scans. Batched traversals count
accepted complete outputs without changing batching. Nested FEM contributes
validated completed-frequency detail while MC retains the primary scan counter.
Unknown totals use `?`. Recovered scans satisfy work without becoming fresh solver
timing evidence; partial native recovery remains incomplete until its result is
accepted. Validation, persistence, reports, and requested timing measurements must
finish before the benchmark becomes terminal.

Subset invocations say `Selected benchmarks`. Failed, skipped, and interrupted
outcomes remain distinct. `Complete` means the requested execution, comparison and
report persistence finished. `Failed` means an execution error, including errors
in calculation, comparison inputs, timing collection or required persistence.
Differences between operands and unavailable relative RMS do not fail a job.
Older progress snapshots are displayed using their recorded execution state;
saved comparisons and timing records are left intact. If the invocation
exhausts its selected execution, ETA says `done`, including when some jobs failed.
An interruption, fail-fast termination, or other early abort says `stopped`. Valid
accepted-scan durations remain useful regardless of numerical differences.

The ASCII spinner advances every 250 ms in the watcher using cached state. Elapsed,
freshness, ETA, and snapshot polling advance about once per second. The spinner
means the viewer is active; an opaque computation can publish no new observation
for a long time. Heartbeats do not reset work freshness or teach duration estimates.
Final elapsed freezes and the spinner stops. Missing, malformed, legacy, or stale
snapshots produce a compact status; staleness never cancels work or implies failure.
The watcher reads metadata and snapshots, never executable declarations, numerical
results, or worker artifacts. Narrow/short terminals and redirected output receive
plain summaries without cursor escapes, color, or spinner frames. Color honors the
IO setting and `NO_COLOR`.

There is one approximate remaining campaign duration, including unfinished
calculation and required overhead/performance work. The collector smooths eligible
observed durations with an EWMA coefficient of 0.25. It uses fresh repeated-operation
throughput, compatible operation timings, coarse backend/mode timing evidence, and
provisional pooled assumptions, keeping scan, whole-operation, controlled-call, and
whole-benchmark timing scopes separate. Batch throughput cannot silently predict
independent scalar calls. Mesh sizes, frequency detail, geometry, and matrix sizes
are not workload models. Unresolved resumed definitions use a provisional whole-
benchmark observation or forecast; a lone scan duration cannot establish seconds
per unresolved benchmark, so `ETA --` is legitimate during bootstrap.

Unmeasured required overhead gets one provisional allowance of at least one second
or ten percent of the whole calculation forecast. The allowance does not shrink
with the calculation countdown and is not added to work already covered by a whole
budget. For active residual budget B and elapsed time A since a genuine estimate
anchor, remaining active time is `max(B-A, 1 + 0.25*max(A-B, 0))` seconds. It becomes
revising only after A exceeds B. Queued work is added once. These allowances and the
overdue extension are assumptions, never measured samples. ETA can increase after
new duration evidence and never enters numerical reuse or performance pass/fail.

Execution publishes at most once per second at ordinary outer boundaries and
forces lifecycle and controlled-call boundaries. Same-directory rename replaces
snapshots; optional output failure disables that output with one diagnostic.
Required result/checkpoint IO failures retain their normal behavior. There is no
producer rendering timer, progress-driven yield, or handshake with the watcher.

Timing records still distinguish invocation wall, operand wall, compute-call wall,
and backend-native scopes. `sessions/SESSION.timing.toml` exists even with progress
off. Operand `timing.toml` and calculation metadata retain durations, callback policy,
and source timings. GetDP phase/process sums are accumulated worker durations, not
elapsed scan time. PSCAD source `timing.txt` measures remote `line.compile()` only;
output readiness and transfer belong to broader adapter time. Recovery retains its
original timing provenance and is not a new zero-second solver measurement.

Before a controlled call, execution publishes sample identity and suspended
observation together, then enters the task-scoped quiet context. The watcher shows
`performance sample 2/3; observation suspended` and can keep animating its cache.
No optional reporting, publication, monitoring transport, or progress-driven yield
occurs within the measured call, including nested MC/backend work. After timing,
execution records the measurement and publishes restored observation with the actual
outcome. Exceptions and interrupts restore scoped state. Normal correctness callbacks
remain enabled; controlled samples retain their recorded callback-free policy and
required numerical/checkpoint IO. Native warmup/repeats remain absent unless already
declared. Allocation counts describe Julia allocations, not native worker memory.

Ordinary publication has measurable overhead. The [validation report](progress-validation.md) compares off,
publisher-only, and publisher-plus-watcher workloads, with variability and the
approximately one-percent median target. This is a measurement target, not a noisy
wall-clock regression assertion. Quiet samples exclude optional execution-side
observation; they do not isolate computation from unrelated machine load.

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

- `definitions.jl`, `benchmarks/`: concrete benchmark composition and catalog.
- `cases.jl`, `cases/`: model inputs, variations and source assets.
- `benchmarks.jl`, `campaigns.jl`: execution, scheduling and recovery.
- `records.jl`, `fingerprints.jl`: source capture and numerical identities.
- `comparisons/`: saved-operand binding and the package comparison rules.
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

Multiple file/name selectors are alternatives, as are multiple `tag:` selectors.
If both groups are supplied, both must match. Append `--list` to inspect the
selection without executing it; `tag:uq` and `tag:pscad` select the owned families.
Omitting selectors retains the full suite. See [test commands](../test/README.md)
for release scopes and unresolved outcomes.
For upstream spectral-control reconciliation, the selected contract covers
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
frequencies, station options or variations. There is no implicit backend catalog
that replaces a declared formulation or suppresses an unavailable selection.
