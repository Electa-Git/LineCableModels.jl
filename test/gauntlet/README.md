# Validation gauntlet

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

`cases/index.toml` is the only catalog. `load_case` resolves only indexed files
beneath `cases/`, verifies the declared ID, builds fresh model state on every
call, and records the source path and SHA-256. Unknown IDs or overrides,
unmatched tag selectors, duplicate paths, escaping paths, and ID/file
mismatches fail immediately.

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

The seven migrated PSCAD benchmarks retain their external formulations,
mappings, and numerical gates. They use the legacy external-reference runner
but now receive a neutral `LoadedCase`; their work and artifact identities are
benchmark IDs, not physical case IDs. A physical case edit intentionally
invalidates its PSCAD snapshot digest and requires a fresh external recording.

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

LINECABLEMODELS_GETDP=/path/to/getdp \
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

Run only one `fem_reference.jl` process at a time. The current Gmsh/GetDP
ONELAB socket transport is not safe across concurrent gauntlet reference
processes.

All three runners are resumable. Their generated records live under
`.linecablemodels/` and are invalidated by case, source, formula, frequency, or
reference changes.

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

## Artifacts and reports

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
│   │   ├── benchmarks/<benchmark-id>/snapshot.{jld2,sha256}
│   │   └── report.{jld2,tsv,sha256}
│   └── uq/
│       ├── benchmarks/<benchmark-id>/snapshot.{jld2,sha256}
│       └── report.{jld2,tsv,sha256}
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
are recomputed on load and report generation.

Display a staged collection report with:

```bash
julia --project=test/gauntlet --startup-file=no test/gauntlet/report.jl uq
julia --project=test/gauntlet --startup-file=no test/gauntlet/report.jl pscad
```

The UQ report exposes maximum mean/std RMS absolute and display-safe relative
errors for every `R/L/C/G` quantity, raw relative maxima, matrix locations,
failing locations, Monte Carlo count and seed, timings, and integrity metadata.
The PSCAD report retains the established `Z/Y` RMS and performance fields while
reporting benchmark ID and physical case ID separately.

After recording and reviewing the reports, commit the definitive source tree.
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
