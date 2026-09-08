# Validation gauntlet

## Manual campaigns

Run indexed cases with the ordinary formulation Gridspace and backend batch
dispatch:

```bash
./cli/lcm gauntlet run --directory /path/to/new-campaign \
  --cases two_bare_wires,cable_18kv_1000mm2_trefoil \
  --backends coaxial,fem,pscad --formulas catalogue
./cli/lcm gauntlet status --directory /path/to/new-campaign
./cli/lcm gauntlet resume --directory /path/to/new-campaign
```

Omit `--cases` for the indexed catalogue. `--formulas default` runs only the
contextual default; `catalogue` additionally varies each registered earth-return
formula independently, keeping other selections at their defaults. The coaxial
backend records inapplicable selections and reasons. FEM records the same
requested identifiers but uses its fixed field equations; identical effective
FEM inputs reuse a solve. PSCAD selects the native methods supported for the
case placement. `--dielectric Ametani2004` explicitly requests lossy insulation
and semicons; the default is lossless.

Catalogue applicability comes from the formula owners' `validate` methods, not
an author allowlist. The planner prepares the resolved assembly geometry once
per case, then checks each selection's placement, earth layers and pair
restrictions without evaluating numerical kernels. Rejected selections retain
the validator's reason. Exceptions during numerical execution remain failures;
they are not retroactively labelled inapplicable. Catalogue descriptions contain
no case-independent applicability verdict, and catalogue size is not fixed.

To vary other formula slots, specify the axes explicitly. Slot names and
identifiers are the same as for `Formulation`; Gauntlet passes them to its normal
`Gridspace` constructor, without a separate combination algorithm:

```bash
./cli/lcm gauntlet run --directory /path/to/dielectric-comparison \
  --cases cable_18kv_1000mm2_trefoil --backends coaxial,fem,pscad \
  --formulas default \
  --select insulation_admittance=default,Ametani2004 \
  --select semicon_admittance=default,Ametani2004 --combine zip
```

This requests two paired selections. `--combine product` instead requests all
four combinations. Singleton axes broadcast in zip mode, following the public
API. All nine formula slots are admitted, including internal/insulation/earth
impedance, insulation/semicon/earth admittance, FrequencyDependent soil (`earth_properties`), EquivalentHomogeneous
and `pipe_impedance`. Equivalent-earth reductions are nested in their consuming
external formula definitions. Adding an identifier to an existing
owner does not require adding it to a CLI allowlist. Explicit selections are
not silently skipped: an unsupported backend/topology fails with its own error.
In particular, selecting pipe formulas does not implement the coaxial pipe
backend. Unspecified slots keep their defaults; an explicit dielectric slot
overrides `--dielectric` for that slot only. The CLI accepts identifiers, not
arbitrary Julia expressions or custom route functions.

Uncertainty propagation uses the same formula selections and case catalogue:

```bash
./cli/lcm gauntlet run --directory /path/to/uq-campaign \
  --cases two_bare_wires --backends coaxial --formulas default \
  --propagation linear_error,monte_carlo \
  --uncertainty-percent 1 --uncertainty-tags geometry,cable_layer \
  --select insulation_admittance=default,Ametani2004 \
  --trials 512 --seed 1234
```

`RelativeStandardUncertainty` supplies the parameter variation: each selected
parameter ID is one uncertain primitive, and repeated uses remain correlated.
Both the percentage and matching tags are explicit. `LinearError` uses
Measurements derivatives; `MonteCarlo` uses the existing normal sampler with
an explicit seed. The same root seed is used for every formula of a case, so
comparisons use the same random draws. The default trial count is 512, not a
claim of convergence for every case. Sampling errors fail the calculation;
there is no implicit rejection/resampling policy.

Linear propagation is supported only by the coaxial backend here. Monte Carlo
can request coaxial or PSCAD calculations, but is explicitly forbidden on FEM.
These checks run before creating a campaign directory or contacting a solver.
UQ operates on each case's existing uncertain problem Gridspace. Each complete
propagation result is checkpointed before the next selection; an interrupted
Monte Carlo calculation restarts that incomplete calculation from its saved
seed. Completed moment artifacts are not rerun. Trial counts, seeds,
uncertainty/correlation declarations and implementation fingerprints are stored.
There are no automatic RMS or performance acceptance judgements.

Add `deterministic` to `--propagation` to also retain nominal phase matrices.
New moment artifacts are single-calculation records; paired UQ snapshots
remain readable through `read_moments`. Neither format implies CI approval.
The documentation summary includes only deterministic default selections.

Use `--frequency-range 0.1,1e6` to run every case on the same 101-point
logarithmic frequency grid. The bounds are recorded in the campaign manifest
and retained by `resume`; original case declarations remain unchanged.
Without this option, each case keeps its frequency range, raised to 0.1 Hz when
necessary, with 101 logarithmic samples. Campaigns do not run FEM Monte Carlo.
FEM uses the package-owned GetDP artifact by default; `LINECABLEMODELS_GETDP`
is the external solver override. PSCAD uses the existing `local.jl`
configuration (or `LINECABLEMODELS_GAUNTLET_CONFIG`).

The manifest fixes the requested selections. A completion callback writes each
result atomically, without replacement, before the next formulation runs. If a
later selection fails, already saved results survive and resume skips them.
Each attempt keeps its own completion/failure record under the job's `attempts/`
directory. Artifact timings record elapsed time at that result's completion;
they are not independent per-formulation performance measurements. Successful
attempt records contain the full batch duration.

Resume verifies numerical input and selected
implementation fingerprints and checks stored numerical payloads before reuse;
changed inputs require a new directory. A Unix process lock prevents concurrent
writers and releases automatically on process exit, including power loss.
`status` distinguishes an interrupted ledger entry from a live locked process.

FEM resumption can also read a compatible completed run for a selection that
was not checkpointed before interruption. It verifies the resolved inputs,
adapter/solver sources, GetDP executable/build identity and raw-result checksums,
then reapplies reduction and output-basis selection without changing the saved
run. Backend-owned Gmsh sessions ignore user configuration files; caller-owned
sessions, UI operation and explicit remeshing are not completed-run reuse paths.
Changed GetDP executable bytes also invalidate campaign checkpoint reuse.

The GetDP `.pro` equations remain parametric. Frequency, excitation and domain
dimensions are supplied through native `-setnumber`/`-setstring` arguments,
while generated includes supply the resolved materials and physical groups.
Preserving the solver source used by a run does not freeze these parameters;
it prevents a working-tree edit from changing the equations halfway through a
calculation or its resumption.

PSCAD uses the same `resume_run_directory` computation option: `nothing` runs a
fresh calculation (the direct API default), `:latest` searches the case's
completed runs, and a path requires that particular run to match. Campaigns
select `:latest`. Reuse compares the exported project, frequency samples, native
solver setting, frozen remote scripts, and the station's PSCAD, `tline.exe` and
master-library identities. Source and raw-output hashes are checked before the
four matrix files are parsed again. Output-basis conversion is reapplied locally;
requested author labels do not prevent reuse of identical effective inputs.
Old runs lacking completion/solver evidence remain readable as historical data
but are not assumed reusable.

PSCAD retains the saved profile, including its license configuration, and reads
the actual selected LCP executable. The station identity is checked before and
after execution, and recorded in campaign
signatures. A change prevents reuse of that campaign's checkpoints. An optional
`solver_identity=PSCADBenchmarks.identify(remote)` pins the same identity across
several direct calls. Reused results carry `execution.reused=true`, zero new
solver execution time, and a separate `source_elapsed_seconds`; this is not a
new solver performance measurement. Remote checking, parsing and conversion
still take time.

PSCAD's `elapsed_seconds` measures the native `line.compile()` call, not the
whole remote operation. `elapsed_scope` records that boundary explicitly;
`source_elapsed_scope` describes the original timing when a result is reused.
Output-readiness polling and transfer are outside this timer. The automation
manual does not promise that returning from
[`compile()`](https://www.pscad.com/webhelp-v502-al/reference/component.html#mhi.pscad.UserCmp.compile)
means every detailed matrix file is ready, so this duration is not advertised
as total simulation time. Older artifacts without a recorded timing scope
remain labelled as unrecorded.

Local warmed performance comparisons require matching Julia, OS, architecture,
CPU model, CPU-thread inventory, Julia/BLAS thread counts and BLAS configuration.
Missing or different environment fields make a comparison diagnostic only.
Coverage/allocation-instrumented executions cannot establish timing regressions.
Matching metadata is necessary, not proof of an idle or frequency-stable machine;
wall-time comparisons still require a controlled runner. The core suite uses
inference and allocation checks independently of wall-clock timing.

The required `mhi.pscad 3.1.2` settings decoder fails on this 5.1 station during
unrelated Fortran-compiler discovery. The identity reader therefore uses that
version's underlying read-only settings call, solely to resolve `file_lcp`.
It neither changes application settings nor selects a license. This narrow
workaround must be revisited when the automation dependency is updated.

These commands execute scientific calculations manually, not through a test
runner or CI. A completed artifact is explicitly **unreviewed** as a numerical
reference. Publication and CI-reference approval remain separate actions.
New campaign records preserve the serialized nominal problem and complete
physical formulation declaration alongside the matrices or UQ moments. The
separate [numerical-reference gate](../numerical/README.md) can replay reviewed
phase records without importing the current case catalogue. Legacy records
remain readable as stored data; missing input declarations are not fabricated.

### Documentation summary

Gauntlet has no standalone report command or HTML/plot generator. Its only report
is the compact defaults comparison in the documentation:

```bash
lcm gauntlet compare --definition /path/to/benchmarks.toml --output /path/to/comparisons
LINECABLEMODELS_GAUNTLET_RESULTS=/path/to/comparisons \
  julia --project=docs docs/make.jl
```

Use the platform path-list separator (`:` on Unix, `;` on Windows) for multiple
comparison directories. The page reads already calculated benchmark errors,
verifies record/operand checksums, and displays the explicit reference and
candidate. It never infers pairings from backend names. Select only the benchmark
records intended for publication; collecting a broad campaign does not publish it.
The deterministic summary has one row per benchmark, with full-band Z/Y errors
in adjacent columns and both RMS normalizations side by side. Frequency slices
are separate sections with the same layout. Conductance G is opt-in for loss
analysis through `compare(reference, candidate, G; ...)`, not a default summary
column. Existing extra comparisons remain stored. UQ means and standard
deviations have their own sections, separate from deterministic comparisons.
The summary does not copy numerical artifacts, generate per-case HTML pages,
embed plots or dump input objects. It never runs missing calculations or selects
the newest run implicitly.

Use the ordinary `LineParameters`, `observe`, `compare` and plotting APIs to
inspect retained data interactively. No Gauntlet-specific detailed reporter is
provided. Collection archives contain snapshots and their checksums only, plus
release metadata; packaging instructions are below.

The gauntlet separates reusable physical models from the calculations used to
validate them:

```text
cases/index.toml
        │
        ▼
cases/<case-id>.jl                 physical model and parameter manifest
        │ load_case(case_id; variation=...)
        ▼
LoadedCase / Gridspace{LineParametersProblem}
        │
        ├── benchmarks/pscad/<benchmark-id>.jl
        │       PSCAD reference versus native-engine candidate
        │
        └── benchmarks/uq/<benchmark-id>.jl
                LEP reference versus Monte Carlo candidate

pscad_reference.jl                 exhaustive deterministic PSCAD references
fem_reference.jl                   one deterministic FEM reference per case
formulation_comparisons.jl         every native formulation versus both backends
```

The hard invariant is one file equals one benchmark. Every Julia file below
`benchmarks/` contains exactly one top-level gauntlet `@testitem`, loads exactly
one indexed case, and declares one benchmark ID matching its filename. Toolkit
tests enforce that rule, ID uniqueness, and the absence of case loops. A case
may be used by any number of benchmarks.

## Cases and the loader

Files in `cases/` contain all physical information from materials, cable parts,
layout, earth, temperature, length, and frequencies through the final
`LineParametersProblem`. They do not contain formulations, tolerances,
comparisons, or test assertions. Each file's final expression is one
`CaseDefinition` with:

- a stable lowercase case ID;
- named `CaseParameter` values and auditable tags;
- a builder closure that receives those named values;
- canonical terminal order.

Fixed imported cases may declare `parameters = (;)`. A `CaseParameter` exists
only when the case intentionally exposes one variable input; it is not
boilerplate required to index a materialised problem.

`cases/index.toml` is the only catalog. `load_case` resolves only indexed files
beneath `cases/`, verifies the declared ID, builds fresh model state on every
call, and records the source path and SHA-256. Unknown IDs or overrides,
unmatched tag selectors, duplicate paths, escaping paths, and ID/file
mismatches fail immediately.

### Importing existing problems

The repository CLI can import a trusted Julia file whose final expression is
one concrete `LineParametersProblem`:

```bash
./cli/lcm gauntlet case import \
  --id two_wire_variant \
  --source /path/to/problem.jl \
  --description "Two buried conductors used for ..."
```

`--description` supplies the documentation-page title. It defaults to the case
ID when omitted.

The source runs in a fresh Julia process with `--startup-file=no`. The importer
normalises only `system_id`, verifies that the coaxial backend can flatten every
design, orders ports by positive phase ID, and writes a versioned JSON problem,
a small `CaseDefinition` wrapper, and the index entry. `--dry-run` validates
without writing; `--force` replaces the same ID. The stored numbers are the
fully materialised values, not a reconstruction of the source expressions.
Pass `--project DIR` when the trusted source belongs to another Julia project;
that project must provide `LineCableModels`.

Use `./cli/lcm gauntlet case list`, `show`, or `validate` to inspect the
catalogue. `case catalogue --output FILE` emits the detached manifest consumed
by the documentation build; `--check` rejects a stale manifest.

Variations are applied before the builder runs:

```julia
model = load_case(
    :cable_132kv_630mm2_flathor;
    variation = RelativeStandardUncertainty(
        10.0;
        tags = (:geometry, :cable_layer),
    ),
)
```

Available policies are `NoVariation`, `ExactOverrides`, `ParameterGrids`,
`RelativeStandardUncertainty`, and `compose_variations`. They use the existing
`Grid` and `Gridspace` semantics. One parameter ID becomes one uncertain
primitive even when the builder reuses it in several layers; different IDs are
independent. This preserves intended covariance instead of wrapping each leaf
constructor argument independently.

The 10% policy means standard uncertainty, not a hard bound:

```text
sigma = 0.1 * abs(nominal)
source = Grid(nominal, 10.0)
```

The layer-geometry selector includes continuous diameters, radii, thicknesses,
strip widths, and nonzero lay ratios. It excludes integer topology, materials,
earth, frequency, temperature, length, cable placement, and fixed design
constraints. All singleton descriptors form one uncertainty-bearing Gridspace
point. LEP materializes that point with Measurements; Monte Carlo repeatedly
realizes the same unresolved point.

Feasibility belongs to the case definition, not to either UQ technique. The
trefoil center spacing is therefore derived as `2.2 * outer_radius`, leaving a
clearance equal to 20% of the realized outer radius for every realization. The
525 kV fixed-count armor case similarly declares a fixed 20% packing-clearance
ratio relative to its unbuffered outer radius. Its compliant bedding absorbs an
extreme residual packing shortfall rather than constructing overlapping armor
wires. Neither clearance ratio is an uncertain manufacturing dimension.

## Benchmark definitions

A benchmark owns its benchmark and case IDs, collection, source digest,
reference and candidate `BenchmarkCalculation`s, comparison policy, tolerances,
and execution options. `run_benchmark` dispatches through the calculations and
comparison policy, so the same runner can compare UQ techniques today and two
Engine formulas or option sets later.

### Compare completed calculations without rerunning them

`benchmark_definition`, `BenchmarkCalculation`, `LineParametersPolicy` and
`UQMomentPolicy` are shared by live and saved calculations. An external owner is
valid on either side. The saved-file CLI accepts explicit scalar bindings; each
Gridspace selection retains its own saved identity. No catalogue or solver is
loaded by `compare`.

```toml
schema_version = 1
collection = "manual"

[comparison]
kind = "line_parameters"
quantities = ["Z", "Y"]
bands = ["all", "dc", "harmonic", "narrow", "wide"]
normalizations = ["reference_rms", "pointwise"]
fundamental = 50.0
harmonics = 50
# Optional quantity-specific zero tolerances, in the stored result basis:
# [comparison.atol]
# G = 1e-12
# C = 1e-16

[[benchmarks]]
id = "trefoil_default_fem_reference"
case = "cable_18kv_1000mm2_trefoil"
description = "18 kV trefoil"
reference = {path = "fem/0001.jld2", sha256 = "<64 hex digits>", owner = "external"}
candidate = {path = "coaxial/0001.jld2", sha256 = "<64 hex digits>", owner = "engine"}
```

Paths are relative to the definition file (absolute paths are also accepted).
Each benchmark may override the complete `comparison` table. For saved UQ
moments, use `kind = "uq_moments"`: R/L/C/G means and standard deviations remain
separate and retain the existing moment comparison semantics. No moments are
mistaken for deterministic Z/Y. Use `owner = "uq"` for those calculations.

`compare_saved` uses the existing `compare`/`RMSError` API. For deterministic
quantities, `:reference_rms` means RMS difference divided by reference RMS;
`:pointwise` means RMS of sample-wise relative differences. Both are stored as
fractions, displayed as percentages. Exact 0/0 contributes zero; nonzero/0 is
infinite. The configured numerical-zero policy applies before either metric,
without modifying arrays or dropping samples. Every maximum identifies its
matrix entry. Full-band errors remain primary; sub-bands use stored samples.

Records retain the ordered checksummed operands, formula selections and
implementation records, input/terminal/basis/frequency identities, comparison
settings, element-wise errors, and timing scope. Missing inputs never trigger a
fallback reference. A source mismatch is an error, not permission to interpolate
or relabel. Existing output records are not overwritten: choose a new directory
when changing comparison settings. Publication and CI approval remain separate.

Campaign `elapsed_at_completion_seconds` measures time since the pending batch
started. It is not an independent cold call or warmed median. These historical
records cannot support a per-selection speedup claim; the documentation says so.

The seven migrated PSCAD benchmarks retain their external formulations,
mappings, and numerical gates. They use the legacy external-reference runner
but now receive a neutral `LoadedCase`; their work and artifact identities are
benchmark IDs, not physical case IDs. A numerical case-input change
intentionally invalidates its PSCAD reference. Comments, descriptions, and
unrelated repository changes do not.

### Exhaustive deterministic formulation comparisons

The cross-backend catalogue covers every case in `cases/index.toml`. It is
separate from the seven legacy benchmark files and does not run uncertainty or
Monte Carlo calculations.

`reference_case` normalizes every external solve to exactly 101 logarithmically
spaced frequencies. The lower bound is `max(first(case.frequencies), 0.1)` Hz,
because PSCAD cannot calculate below 0.1 Hz, and the upper bound is the case's
declared maximum. PSCAD and FEM therefore produce directly comparable tensors
on the same frequency axis.

Run the stages in order:

```bash
julia --project=test/gauntlet --startup-file=no test/gauntlet/pscad_reference.jl

julia --project=test/gauntlet --startup-file=no test/gauntlet/fem_reference.jl

julia --project=test/gauntlet --startup-file=no \
  test/gauntlet/formulation_comparisons.jl
```

The PSCAD runner executes every field formulation applicable to the case's
overhead, underground, or mixed placement and records honest skips for the
other fields. The FEM runner executes each case once. The comparison runner
then computes every applicable registered native formulation and compares it
with every available PSCAD reference and the FEM reference. Its locked metric
is element-wise absolute and reference-normalized RMS error across frequency
for every entry of `Z` and `Y`.

Run only one `fem_reference.jl` process against its shared output collection at
a time. The backend now assigns a separate ONELAB socket to each run, so
independent processes with distinct run directories do not share that socket.
This does not make Gmsh's process-global API safe for concurrent Julia tasks,
or retrofit socket isolation into an already-running older process.

All three runners are resumable. Their generated records live under
`.linecablemodels/`. Numerical reuse is keyed by the materialised problem and
its resolved geometry, terminal ownership and longitudinal paths,
selected formulation, exact Git blob identities of the selected formulas and
shared numerical implementation (including flattening equivalences), and reference bytes. Adding an unrelated
formula, plot, extension, or library entry does not invalidate prior results.
Every newly written artefact also records the full repository commit and dirty
state for historical provenance; official publication requires a clean tree.

FEM also fingerprints the resolved material domains, boundary shapes and mesh
settings actually handed to Gmsh. A declaration-only fingerprint from an older
schema is insufficient for automatic reuse. Such records remain readable for
stored results; they are not rewritten, relabelled or automatically recomputed.

### LEP versus Monte Carlo

The legacy uncertainty collection has one owned LEP-versus-Monte-Carlo
benchmark for each of its seven cable and overhead-line models. New indexed
cases are covered by the exhaustive deterministic catalogue above without
inventing uncertainty experiments for them.

Both calculations consume the same single-point `ParametricProblem`, the same
inner native `Formulation`, and the same execution options:

- reference: `LinearError(inner)`;
- candidate: `MonteCarlo(inner; trials=<fixed>, seed=<fixed>,
  distribution=:normal)`;
- retained samples and histograms: disabled.

The trial count is explicit; automatic DKW sizing is not used. Five cases use
512 trials. The 525 kV armor and two-bare-wire cases use 2,048 trials because a
fixed-seed convergence check showed that their 512-trial standard-deviation
estimates were not stable enough for the common gate. This is convergence by
increasing the sample, not seed selection or rejection/resampling of physical
trials.

“Meaningful” relative error excludes a ratio when that term's absolute RMS is
already below its quantity/statistic-specific floor. The locked gates are 5%
for means and 10% for standard deviations, with explicit absolute floors. At
10% input uncertainty, the gate represents practical engineering equivalence,
not numerical identity: LEP is first-order and local, whereas Monte Carlo also
contains nonlinear propagation and finite-sample noise.

The accepted local recording produced:

| Case | Trials | Maximum meaningful mean difference | Maximum meaningful uncertainty difference | Monte Carlo / LEP median |
|:--|--:|--:|--:|--:|
| 132 kV 630 mm² flat horizontal | 512 | 1.23% | 3.58% | 23.37× |
| 18 kV 1000 mm² trefoil | 512 | 1.88% | 8.51% | 6.25× |
| 380 kV 2000 mm² flat vertical | 512 | 1.68% | 9.16% | 24.99× |
| 525 kV 1600 mm² bipole | 2,048 | 1.31% | 9.25% | 80.76× |
| 640 kV 2000 mm² bipole | 512 | 0.86% | 6.68% | 21.28× |
| Solid 1000 mm² single phase | 512 | 0.68% | 7.33% | 34.28× |
| Two bare wires | 2,048 | 1.75% | 6.24% | 165.28× |

Across this family, the largest meaningful mean difference is 1.88%, the
largest propagated-uncertainty difference is 9.25%, and the smallest observed
speedup is 6.25×. Timing values are machine-specific; their environment is
stored in the artifact.

The scientific KPIs are mean and standard deviation for `R`, `L`, `C`, and `G`.
LEP moments come from `nominal` and `uncertainty`; Monte Carlo moments
come from `SampleSummary.mean` and `.std`. Before numerical comparison, the
adapter requires exact agreement in quantity set, tensor shape, terminal order,
frequency samples, domain, and basis. It then computes absolute and
reference-normalized RMS errors for each matrix term across frequency. A term
passes when its absolute or relative limit passes. This comparison gates every
mode.

Timing is a separate warmed benchmark of the same two calculations. Each side
requests up to three `BenchmarkTools` samples within a 20-second sampling
budget after the scientific executions; calculations longer than that budget
still produce one complete warmed sample. The reported figure is
`Monte Carlo median / LEP median`. A comparable local run must show at least a
2× LEP advantage. Coverage or allocation-instrumented runs retain the timing
diagnostics but do not gate performance.

## Modes

`LINECABLEMODELS_GAUNTLET_MODE` accepts `snapshot`, `live`, or `record`.

| Mode | PSCAD reference | Owned reference | Accepted artifact | Writes staging |
|:--|:--:|:--:|:--:|:--:|
| `snapshot` | Loaded, never executed | Executed locally | Required | No |
| `live` | Executed | Executed locally | No | No |
| `record` | Executed | Executed locally | Published collection when available | Yes |

Snapshot mode for PSCAD does not load host configuration, initialize Python,
start a process, or contact a network. Owned benchmarks still execute both
current calculations in snapshot mode, then regress both moment products
against their accepted artifact. CI permits only snapshot mode.

Instantiate the isolated environment from the repository root:

```bash
julia --project=test/gauntlet --startup-file=no -e \
  'using Pkg; Pkg.instantiate()'
```

Run the complete gauntlet:

```bash
LINECABLEMODELS_GAUNTLET_MODE=snapshot \
julia --project=test/gauntlet --startup-file=no test/gauntlet/runtests.jl
```

Run one family while developing. The PSCAD command below is explicitly live;
without the mode setting it would use the default snapshot mode:

```bash
julia --project=test/gauntlet --startup-file=no -e '
using TestItemRunner
TestItemRunner.run_tests(joinpath(pwd(), "test");
          filter=ti -> :uq in ti.tags, verbose=true)
'
LINECABLEMODELS_GAUNTLET_MODE=live \
julia --project=test/gauntlet --startup-file=no -e '
using TestItemRunner
TestItemRunner.run_tests(joinpath(pwd(), "test");
          filter=ti -> :pscad in ti.tags, verbose=true)
'
```

Run the reusable infrastructure checks with:

```bash
julia --project=test/gauntlet --startup-file=no -e \
  'using TestItemRunner; TestItemRunner.run_tests(joinpath(pwd(), "test"); filter=ti -> :gauntlet_toolkit in ti.tags, verbose=true)'
```

Record collections only through the dedicated runner:

```bash
LINECABLEMODELS_GAUNTLET_PERSIST=true \
LINECABLEMODELS_GAUNTLET_MODE=record \
LINECABLEMODELS_GAUNTLET_STAGE_FORCE=true \
LINECABLEMODELS_GAUNTLET_CLEANUP=true \
julia --project=test/gauntlet --startup-file=no test/gauntlet/runtests.jl
```

Set `LINECABLEMODELS_GAUNTLET_STAGE_FORCE=true` to replace the complete
unversioned staging area. This setting never replaces a release package or an
`Artifacts.toml` binding. Set
`LINECABLEMODELS_GAUNTLET_CLEANUP=true` to remove
`test/gauntlet/benchmarks/.work/` after a fully successful run. Failed runs
retain diagnostics.

## PSCAD live setup

Copy `local.example` to the ignored `local.jl` and configure a `RemoteConfig`.
Its shared Windows root must point to
`test/gauntlet/benchmarks/.work`. The host requires PSCAD 5.1.0, a working
license, Julia 1.12, PythonCall 0.9, and `mhi.pscad` 3.1.2. Instantiate the
remote project once:

```powershell
$env:JULIA_PYTHONCALL_EXE = "C:\Python311\python.exe"
julia --project=test/gauntlet/pscad/remote --startup-file=no -e \
  'using Pkg; Pkg.instantiate(); using PythonCall'
```

`transport=:ssh` uses system SSH. A custom transport defines
`remote_command(::Val{:name}, config, powershell)` in `local.jl`; the supplied
Tailscale example routes to the local libvirt guest. Authentication, VM setup,
licenses, and tunnels remain outside the gauntlet.

Generated projects and logs live below:

```text
test/gauntlet/benchmarks/.work/pscad/<benchmark-id>/reference/
```

The supervisor records the exact remote Julia PID and runner path, targets only
that process tree on cancellation, and preserves the shared and Windows scratch
directories after failure. `verbosity=(default=0, PSCAD=2)` streams milestones;
PSCAD's blocking `compile()` call cannot stream intermediate project messages.

## Numerical artifacts

Artifact lifecycle and release versioning are intentionally separate. Julia
owns benchmark execution, validation, snapshot schema 2, and unversioned local
staging. It does not inspect Git tags, choose a version, create a release, or
upload an archive. The external packaging command owns those release concerns
and operates on one collection at a time.

This is the locked layout:

```text
test/gauntlet/.artifacts/
├── staging/
│   ├── pscad/
│   │   └── benchmarks/<benchmark-id>/snapshot.{jld2,sha256}
│   └── uq/
│       └── benchmarks/<benchmark-id>/snapshot.{jld2,sha256}
└── releases/<collection>/vX.Y.Z/
    ├── benchmarks-<collection>-vX.Y.Z.tar.gz
    └── package.toml
```

Staging and release packages are ignored working data. `Artifacts.toml` is the
tracked runtime registry and uses stable keys such as `gauntlet_pscad` and
`gauntlet_uq`; advancing a collection updates its stable binding. Collections
version independently through `gauntlet-pscad-vX.Y.Z`,
`gauntlet-uq-vX.Y.Z`, and corresponding future tag families. A snapshot schema
change is independent of any collection release version.

Each snapshot stores separate case and benchmark IDs and SHA-256 values, the
parameter manifest, applied variation, parameter-identity correlation record,
reference/candidate calculation records and options, tolerances, terminal and
frequency metadata, results or plain moment products, comparison, seed/trials,
timings, environment, and timestamp. Snapshot digests and stored comparisons
are checked by `read_collection` and before packaging. Finalizing a collection
validates its snapshots; it does not write a report or duplicate results.

After recording and reviewing the numerical results, commit the definitive source tree.
Then package one collection by supplying the exact next version and a reason:

```bash
python .github/scripts/package_gauntlet.py \
  --collection uq \
  --version 1.0.0 \
  --reason "Initial LEP versus Monte Carlo baseline"
```

The script requires a clean worktree (ignored staging remains available), reads only
`gauntlet-<collection>-v*` Git tags, and accepts `1.0.0` for an unreleased
collection or exactly one patch, minor, or major successor of its latest tag.
It passes the explicit version, reason, and full Git commit to Julia's narrow
packager. Add `--create-tag` to create the validated annotated tag locally;
uploading assets and pushing tags remain explicit release/CI operations.
An explicitly dispatched CI job may call the same script only after staging has
been supplied to that job; an ordinary hosted checkout cannot see ignored local
staging.

After the archive has a real immutable URL, update the stable runtime binding:

```bash
julia --project=test/gauntlet --startup-file=no test/gauntlet/bind.jl \
  uq 1.0.0 \
  https://github.com/OWNER/REPOSITORY/releases/download/gauntlet-uq-v1.0.0/benchmarks-uq-v1.0.0.tar.gz
```

Commit the resulting `Artifacts.toml` change. There is deliberately no
pre-commit hook: ordinary Julia runs stage data, while an explicit Git-aware
release action versions and packages it. Published collection releases are
immutable; corrections use the next version of only the affected collection.
