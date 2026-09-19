# # Gauntlet
#
# Gauntlet declares a case, an explicit reference and a scalar or Gridspace of
# candidate formulations. Its calculations use the ordinary `compute` methods;
# ReportBuilder owns the comparison bands, grouping and tables, and PlotBuilder owns
# matrix overlays. A reference defines comparison direction and normalization.
# LCM, PSCAD and FEM remain different models whose discrepancies are measured.
#
# ## Catalog uncertainty declarations
#
# New catalog UQ declarations record `bounded_correlated_geometry_v1`: one
# bounded uniform 10%-standard-deviation scale for cross-section lengths,
# independent bounded lay/fillet inputs, and an independent explicit area input.
# This synthetic correlated study is not the old independent-normal study and
# does not assert measured manufacturing correlations. Counts stay fixed over
# the declared support. Derived areas scale quadratically, so first-order LEP
# and nonlinear Monte Carlo need not have identical means.
#
# Joint sources use ordinary core Gridspace builders. Gauntlet records their
# primitive supports and output dependencies and passes each joint source once.
# A fresh declaration adopts the law; resume retains the saved study. Existing
# results and solver caches are preserved. Removing the artificial buffer in
# the 525 kV/1600 mm² case changes its nominal geometry and requires matching new
# references; the 320/380 kV nominal resolved geometry is retained explicitly.
#
# ## Retained summaries
#
# This page reads the immutable versions selected in `docs/gauntlet.toml`. An
# explicit `LINECABLEMODELS_GAUNTLET_RESULTS` directory enables a labeled local
# preview. Building this page performs no solve, RMS calculation or plotting.
#
# Each benchmark has one numeric DataFrame per physical quantity, with unique
# relevant formulations as rows and retained bands side by side: `all`, `dc`,
# `harmonic`, `narrow`, `wide`. UQ mean and standard deviation have separate tables.
# Each cell is the **maximum of eligible per-term relative RMS discrepancies [%]**,
# not a whole-matrix RMS or an average across formulations. Missing stays missing;
# absolute errors, winning terminal pairs and availability counts remain in the
# detailed report tables. References identify comparison methods, not ground truth.
#
# Compact performance tables show recorded whole-workload timings, timed calls,
# cumulative Julia allocations and the reference/candidate time ratio with its
# comparability flag. MC workload and CDF precision are separate from physical spread.
# Missing measurements are not manufactured. No plots, thumbnails or retained
# illustrations are embedded; detailed plots remain an explicit inspector request.
#
# <!-- GAUNTLET_REPORT -->
#
# ## Numerical comparison
#
# With **A = reference** and **B = candidate**, each matrix term uses
# [`compare`](@ref LineCableModels.Engine.compare):
#
# ```math
# \mathrm{NRMSE}=100\sqrt{\frac{\sum_k|B_k-A_k|^2}{\sum_k|A_k|^2}},
# \qquad
# \mathrm{RMS}_{\mathrm{pointwise}}=100\sqrt{\frac1N\sum_k\left|\frac{B_k-A_k}{A_k}\right|^2}.
# ```
#
# Absolute RMS retains the measured difference. No denominator floor is applied.
# Both operands must exceed the reporting resolution at every selected sample.
# Otherwise relative RMS is missing for either normalization, with a per-term
# reason and absolute RMS; samples are never silently omitted.
# Defaults are 1e-10 Ω/m for R, 1e-15 H/m for L, 1e-12 S/m for G and 1e-16 F/m
# for C. X/B thresholds follow 2πfL/2πfC; Z/Y follow R + 2πfL and G + 2πfC. Explicit `atol` overrides
# belong to the comparison request. Empty bands and unsupported quantities retain
# their separate reasons. No error threshold is imposed by artifact acceptance.
#
# These are native-unit reporting cutoffs, not certified floating-point forward-error
# bounds. The corresponding `:total` defaults use Ω, H, S and F; they do not infer a
# line length. Equivalent per-length/total comparisons require correspondingly
# scaled explicit cutoffs. Absolute RMS uses unchanged raw values even when relative
# RMS is missing. A missing whole-band KPI does not remove a significant part of
# that band's plotted curve.
#
# New analyses retain their resolution revision, effective cutoffs and per-operand
# unresolved-sample counts. The existing analysis identity includes these semantics,
# without changing any calculation identity. Historical/unversioned RMS remains
# readable as recorded; `compare_saved` explicitly creates current analysis from
# saved operands without a solver run. Reading or plotting never silently rewrites
# old errors. Report tables expose the revision; plots warn when historical RMS and
# current observation resolution differ.
#
# Report bands overlap: near DC is 0.1–100 Hz, the default harmonic range is
# 50–2500 Hz, narrowband is 1 kHz–1 MHz, and wideband is strictly above 1 MHz.
# Ordinary endpoints retain the engine's nearest-sample selection. Requested and
# actual bounds and sample counts remain visible. There is no interpolation.
#
# ## Inspect a benchmark
#
# ```julia
# using LineCableModels
# using LineCableModels.ReportBuilder: BenchmarkTableDefinition
# include(joinpath(pkgdir(LineCableModels), "gauntlet", "Gauntlet.jl"))
# using .Gauntlet
#
# benchmark = read_benchmark("/path/to/staging/benchmark_id")
# # A vault can be moved and read independently of its original staging folder:
# # benchmark = only(read_campaign("/path/to/vault/accepted"))
# artifact = report(BenchmarkTableDefinition(), benchmark)
# display(artifact)             # Compact published-summary layout in the REPL
# artifact.table.features      # Quantity/statistic DataFrames, bands side by side
# artifact.table.overview      # Compact coverage, performance and sampling tables
# artifact.table.formulations
# artifact.table.maxima
# artifact.table.terms
#
# using GLMakie
# plots = LineCableModels.plot(artifact; ydata=(Z, Y))
# # Optional: plot(artifact, (R, L); problem=1, formulations=[2], band=:dc)
# export_svg(first(plots); path="impedance.svg", open_file=false)
# ```
#
# The reference and unique quantity-relevant candidates appear in each matrix
# cell, including both off-diagonals. Z gives R/X pages; Y gives G/B pages.
# Complete saved selections and source colors survive filtering and reload.
# Multiple problem points require an explicit `problem` selection. Equal numerical
# curves remain separate formulation choices. UQ mean/std errors remain separate
# statistics, with the same selectable frequency bands and physical resolution as
# deterministic quantities. `artifact.table.features` contains numeric formula-by-band
# absolute/relative DataFrames; `statistics` retains available UQ summaries and
# `sampling` records finite-sample precision separately from physical spread.
# Execution/native timings and controlled performance samples are available in
# `execution`, `source_timings`, `performance` and `performance_samples`. These are
# recorded measurements, never new benchmark runs triggered by a report.
# `overview` supplies the scalar compact views used by the REPL, HTML publication
# and both `dev/inspect_saved_gauntlet*.jl` scripts. Detailed statistics and native
# timing records are not dumped into the default summary. In the inspectors,
# `overview_tables`, `performance_display_df`, `timing_ratio_df` and
# `band_coverage_df` remain ordinary IDE-accessible tables. Those scripts explicitly
# reanalyse saved operands using their selected bands; publication uses saved RMS.
#
# Backend and method names come from owned `description` methods. References are
# labelled `Reference · FEM`, `Reference · PSCAD`, or `Reference · Monte Carlo`;
# candidate labels are unnumbered, for example `LEP`. Common coaxial
# identity is omitted. Quantity-specific legends show the applicable equation
# choices; `formulations` lists labels and `formula_details` supplies scientific
# explanations. Complete declarations remain in the published operand metadata.
# `formulations=[3,1]` selects stored array positions before deduplication.
# Repeated quantity-relevant selections share one curve and feature-table row;
# every relevant route and control participates in equality, not the label or
# numerical curve. Conflicting repeated observations raise an error.
# Explicit composite selections retain every named branch in the relevant legend,
# even when branches are unchanged or default. Internal selections use
# `inner`/`outer`/`transfer`, and earth selections use `air`/`earth`/`mixed`.
# Formula-local overrides also remain visible when shared by every candidate.
# Short and detailed names both use the owning formula's `description` method;
# identifiers and calculation-reuse decisions do not depend on that text.
#
# `sampling` contains scalar MC counts and CDF bounds. `mean_sampling_precision`
# contains one numeric standard error per quantity, point, terminal pair and
# frequency, with its physical unit. This is precision of the sampled mean, not
# precision of the estimated standard deviation. That latter precision is not
# established by a CDF bound or small LEP/MC RMS discrepancy.
# Performance tables retain seconds, Julia allocated bytes/MiB, actual timed
# repetitions, requested repetitions and the original timing scope. Complete
# workload/session records remain in the publication rather than table cells.
# `terms` and `maxima` keep method labels, numeric errors and separately
# filterable response/excitation coordinates. Requested/actual frequency bounds
# have numeric lower/upper columns. `comparisons` is the explicit structured
# audit view; raw maxima remain in DataFrame metadata for unchanged saved-summary
# writing. Retained consumed formula IDs take precedence when available. With
# only a declaration, descriptions state its selected routes, not inferred
# geometry-dependent execution.
#
# The same report and plot APIs accept live completed results without Gauntlet:
#
# ```julia
# reference = compute(problem, Formulation())
# candidates = compute(problem, Formulation(earth_properties=Grid((:constant, :longmire1975))))
# artifact = report(BenchmarkTableDefinition(), (; reference, candidate=candidates))
# ```
#
# Native output terminal identities, matrix basis/domain and frequency coordinates
# must agree for RMS. Externally constructed results require explicit terminal
# metadata when they lack it. Different constitutive laws and formulations are
# legitimate comparisons. Two result spaces need explicit reference/candidate pairing.
#
# `artifact.published` retains unformatted scientific products. Report display and
# the standard page show compact summaries; complete RMS matrices and per-term data
# remain in the detailed tables. Selecting saved products does not recalculate RMS.
# An unrecorded band or changed tolerance requires explicit reanalysis of raw operands.
#
# ## Draft, inspect, lock, package, bind
#
# ```bash
# lcm gauntlet run --definition gauntlet/benchmarks/examples/compare_soil.jl \
#   --directory /path/to/staging
# lcm gauntlet status --directory /path/to/staging
# lcm gauntlet resume --directory /path/to/staging
# lcm gauntlet lock --directory /path/to/staging --benchmark compare_soil \
#   --expected INSPECTED_SNAPSHOT --output /path/to/vault/soil
# lcm gauntlet package --definition release.toml --output /path/to/package
# # Upload the generated archive with the hosting service's existing tool, then:
# lcm gauntlet bind --package /path/to/package --url https://host.example/archive.tar.gz
# ```
#
# `run` skips identical completed benchmarks, continues identical drafts, and creates
# new attempts for changed declarations. Matching operands and formulation points
# are reused before any backend is entered. `--force` explicitly requests fresh
# calculations; historical attempts remain retained. Failed replacement leaves
# the previous complete result explicitly accessible with `previous=true`. `resume`
# preserves the saved work order, skipping complete entries before restoring builders
# or loading results. `--benchmark ID[,ID]` selects run/resume entries; `run --dry-run`
# prints per-operand scheduling decisions without writes or solver calls.
#
# A completed skip verifies declaration and numerical/report payload checksums;
# full native-evidence verification remains in saved-artifact readers and verified
# status. Saved-result reuse does not contact PSCAD. `--recover-solvers` applies only
# to unfinished native work without a matching Gauntlet calculation. Source edits
# alone do not invalidate results; use explicit force for implementation corrections
# with unchanged numerical declarations. Skipped Julia outcomes have
# `state=:complete`, `skipped=true` and `result=nothing`.
# Report-only changes retain matching controlled timing observations and their
# original session; they do not silently repeat the numerical workload for timing.
#
# `lock` copies selected accepted
# benchmarks into an immutable, self-contained bundle; it does not publish them or
# apply a numerical agreement threshold. Packaging accepts only explicit locked
# bundles. Binding verifies the served archive bytes and extracted tree before
# recording an immutable version-specific download.
#
# Catalog defaults remain 101 logarithmic samples from 0.1 Hz to 10 MHz (100
# increments). Frequency overrides are authoritative and backend limitations are
# validated separately. Declarations are saved before execution. Resume checks
# numerical inputs and stored-file integrity, not the live source tree. Each
# execution session records its environment; reused operands keep their original
# execution records. Actual grids, solver details, and timing scopes are retained.
#
# Campaign progress is enabled by default. `run_campaign(...; progress=:auto)`
# publishes lightweight snapshots and prints an exact-session watch command and final
# summary. `:plain` adds throttled single-line execution status; `:off` disables
# optional observation, estimation, and publication. The CLI uses
# `--progress auto|plain|off`.
#
# Run the printed command in a separate terminal:
#
# ```bash
# ./gauntlet/lcm gauntlet status --directory DIR --watch --session SESSION
# ```
#
# `--session` pins the selected invocation across benchmarks and final closure.
# Directory-only watching remains valid; ambiguous sessions show inventory without a
# combined ETA. The watcher reads lightweight metadata and snapshots only. Closing it
# never affects computation, and session identity or stale observations do not prove
# solver liveness.
#
# The six-row display counts selected terminal benchmarks and accepted complete
# frequency scans. MC retains its accepted-trial counter when a child FEM call supplies
# optional validated-frequency detail. Unknown scan totals use `?`. Validation,
# persistence, reports, and declared performance work remain unfinished after the
# last scan. `Complete` means that work finished; `Failed` means an execution error.
# Numerical differences, unavailable relative RMS and timing ratios remain reported
# observations, without acceptance verdicts. Exhausting selected execution shows
# `ETA done`, including when jobs failed; early abort shows `ETA stopped`.
#
# One approximate campaign ETA combines remaining operation, overhead, and declared
# performance budgets. Scope-correct observations and explicitly provisional fallback
# forecasts replace workload scaling; an individual scan seed cannot predict an
# unresolved benchmark's multiplicity. ETA may be unavailable during bootstrap or
# increase when costs change. Backend frequency detail and heartbeats do not train it.
# Only actual work changes reset freshness. The watcher animates its cache independently
# of execution, freezes elapsed at closure, honors `NO_COLOR`, and uses plain output
# when the terminal cannot fit the panel or output is redirected.
#
# Controlled performance calls publish their sample identity and suspended state
# before timing. Nested UQ/backend work performs no optional reporting, snapshot IO,
# or monitoring transport during that call; the external watcher may keep animating.
# Actual outcome and restored observation publish afterward, including on exceptions.
# Normal callbacks remain enabled outside the separate timing pass. Execution wall,
# compute-call wall, GetDP worker-time sums, and PSCAD compile-only timings retain
# their separate scopes. Recovery is not a new cold timing sample. Solver diagnostics are independent of the progress switch.
#
# The [Gauntlet CLI guide](https://github.com/Electa-Git/LineCableModels.jl/blob/main/gauntlet/README.md)
# contains declaration, release and illustration-file examples. Standard publication
# pins artifact versions. Explicit plot exports may be retained with their display
# selections in artifacts, but the results-summary page never embeds or copies them.
#
# Stored comparisons are not automatically numerical references for CI. The separate
# [numerical-reference gate](https://github.com/Electa-Git/LineCableModels.jl/tree/main/test/numerical)
# continues to require reviewed references and explicit tolerances.
#
# ## Benchmark data API
#
# ```@docs
# LineCableModels.Engine.RMSError
# LineCableModels.Engine.LineParametersBenchmark
# LineCableModels.Engine.compare
# LineCableModels.Engine.absolute_error
# LineCableModels.Engine.relative_error
# ```
