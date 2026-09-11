```@meta
EditURL = "../literate/gauntlet.jl"
```

# Gauntlet

Gauntlet declares a case, an explicit reference and a scalar or Gridspace of
candidate formulations. Its calculations use the ordinary `compute` methods;
ReportBuilder owns the comparison bands, grouping and tables, and PlotBuilder owns
matrix overlays. A reference defines comparison direction and normalization.
LCM, PSCAD and FEM remain different models whose discrepancies are measured.

## Retained summaries

This page reads the immutable versions selected in `docs/gauntlet.toml`. An
explicit `LINECABLEMODELS_GAUNTLET_RESULTS` directory enables a labelled local
preview. Building this page performs no solve, RMS calculation or plotting.

Each row identifies the case, problem, complete candidate formulation, reference
and comparison snapshot. The default sections are the entire range, near DC,
harmonic range, narrowband and wideband. Each quantity cell shows the **maximum
of per-term RMS discrepancies**, its terminal pair, and unavailable-term count.
This is neither a whole-matrix RMS nor an average across formulations. If relative
RMS is unavailable, the summary retains absolute RMS and its units.

The default page contains no matrix figures or plot thumbnails. Detailed plots
require an explicit request for the selected benchmark/problem. Only explicitly
exported and retained illustrations are embedded during publication.

No published benchmark artifacts are selected. Add immutable version bindings
to `docs/gauntlet.toml` after `lcm gauntlet package`, upload and `lcm gauntlet bind`.
An explicit `LINECABLEMODELS_GAUNTLET_RESULTS` directory enables a local draft preview.


## Numerical comparison

With **A = reference** and **B = candidate**, each matrix term uses
[`compare`](@ref LineCableModels.Engine.compare):

```math
\mathrm{NRMSE}=100\sqrt{\frac{\sum_k|B_k-A_k|^2}{\sum_k|A_k|^2}},
\qquad
\mathrm{RMS}_{\mathrm{pointwise}}=100\sqrt{\frac1N\sum_k\left|\frac{B_k-A_k}{A_k}\right|^2}.
```

Absolute RMS retains the measured difference. No denominator floor is applied.
Numerical-zero reference traces give unavailable relative RMS, with a per-term
reason and absolute RMS. Pointwise normalization is unavailable if any selected
reference sample is numerically zero; samples are never silently omitted.
Defaults are 1e-10 Ω/m for R, 1e-15 H/m for L, 1e-12 S/m for G and 1e-16 F/m
for C. Z and Y thresholds follow R + 2πfL and G + 2πfC. Explicit `atol` overrides
belong to the comparison request. Empty bands and unsupported quantities retain
their separate reasons. No error threshold is imposed by artifact acceptance.

Report bands overlap: near DC is 0.1–100 Hz, the default harmonic range is
50–2500 Hz, narrowband is 1 kHz–1 MHz, and wideband is strictly above 1 MHz.
Ordinary endpoints retain the engine's nearest-sample selection. Requested and
actual bounds and sample counts remain visible. There is no interpolation.

## Inspect a benchmark

```julia
using LineCableModels
using LineCableModels.ReportBuilder: BenchmarkTableDefinition
include(joinpath(pkgdir(LineCableModels), "gauntlet", "Gauntlet.jl"))
using .Gauntlet

benchmark = read_benchmark("/path/to/staging/benchmark_id")
# A vault can be moved and read independently of its original staging folder:
# benchmark = only(read_campaign("/path/to/vault/accepted"))
artifact = report(BenchmarkTableDefinition(), benchmark)
artifact.table.summary
artifact.table.formulations
artifact.table.maxima
artifact.table.terms

using GLMakie
plots = LineCableModels.plot(artifact; ydata=(Z, Y))
# Optional: plot(artifact, (R, L); problem=1, formulations=[2], band=:dc)
export_svg(first(plots); path="impedance.svg", open_file=false)
```

All selected formulations and the reference appear in each corresponding matrix
cell, including both off-diagonals. Z gives R/X pages; Y gives G/B pages. Original
formulation indices, complete selections and colors survive filtering and reload.
Multiple problem points require an explicit `problem` selection. Equal numerical
curves remain separate formulation choices. UQ mean/std errors remain separate
statistics under their existing full-band comparison contract.

The same report and plot APIs accept live completed results without Gauntlet:

```julia
reference = compute(problem, Formulation())
candidates = compute(problem, Formulation(earth_impedance=Grid((:default, :Pollaczek1926))))
artifact = report(BenchmarkTableDefinition(), (; reference, candidate=candidates))
```

Native output terminal identities, matrix basis/domain and frequency coordinates
must agree for RMS. Externally constructed results require explicit terminal
metadata when they lack it. Different constitutive laws and formulations are
legitimate comparisons. Two result spaces need explicit reference/candidate pairing.

`artifact.published` retains unformatted scientific products. Report display and
the standard page show compact summaries; complete RMS matrices and per-term data
remain in the detailed tables. Selecting saved products does not recalculate RMS.
An unrecorded band or changed tolerance requires explicit reanalysis of raw operands.

## Draft, inspect, lock, package, bind

```bash
lcm gauntlet run --definition gauntlet/benchmarks/examples/compare_soil.jl \
  --directory /path/to/staging
lcm gauntlet status --directory /path/to/staging
lcm gauntlet resume --directory /path/to/staging
lcm gauntlet lock --directory /path/to/staging --benchmark compare_soil \
  --expected INSPECTED_SNAPSHOT --output /path/to/vault/soil
lcm gauntlet package --definition release.toml --output /path/to/package
# Upload the generated archive with the hosting service's existing tool, then:
lcm gauntlet bind --package /path/to/package --url https://host.example/archive.tar.gz
```

A new `run` replaces only selected drafts once complete. Failed replacement leaves
the previous complete result explicitly accessible with `previous=true`. `resume`
reuses only matching completed calculations. `lock` copies selected accepted
benchmarks into an immutable, self-contained bundle; it does not publish them or
apply a numerical agreement threshold. Packaging accepts only explicit locked
bundles. Binding verifies the served archive bytes and extracted tree before
recording an immutable version-specific download.

Catalogue defaults remain 101 logarithmic samples from 0.1 Hz to 10 MHz (100
increments). Frequency overrides are authoritative and backend limitations are
validated separately. Declarations are saved before execution. Resume checks
numerical inputs and stored-file integrity, not the live source tree. Each
execution session records its environment; reused operands keep their original
provenance. Actual grids, backend evidence and timing scopes are retained.

The [Gauntlet CLI guide](https://github.com/Electa-Git/LineCableModels.jl/blob/main/gauntlet/README.md)
contains declaration, release and illustration-file examples. Standard publication
pins artifact versions. Its default illustration list is empty; explicit exports
are retained with their display selections and embedded without rerendering.

Stored comparisons are not automatically numerical references for CI. The separate
[numerical-reference gate](https://github.com/Electa-Git/LineCableModels.jl/tree/main/test/numerical)
continues to require reviewed references and explicit tolerances.

## Benchmark data API

```@docs
LineCableModels.Engine.RMSError
LineCableModels.Engine.LineParametersBenchmark
LineCableModels.Engine.compare
LineCableModels.Engine.absolute_error
LineCableModels.Engine.relative_error
```
