# Changelog

All changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and versions follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `Engine.AbstractModalOperators` is the supertype of modal-to-phase operator bases.
  `ModalOperators` subtypes it and implements `size`.
- `Pose2` is callable. `pose(point)` rotates a point by the pose's orientation and
  translates it by the pose's position.
- `initialize_buffers(reduction::ReductionPlan, T, input, plan, buffers)` adds
  `reduction`, the `ReductionBuffers` of a reduction plan, to a buffer record.

### Changed

- `validate` is owned by `Commons`. `LineCableModels.validate` and
  `import LineCableModels: validate` are unchanged. Every `validate` method takes the
  checked input first and returns it unchanged or throws. Further arguments say what the
  input is checked for. These signatures changed:
  - `validate(settings, compare)` replaces `validate(compare; settings...)`.
  - `validate(workbooks, definition)` checks XLSX workbooks.
  - The dielectric and temperature law checks take the law result first and no longer
    convert it.
  - The Unified check is `validate(Γ, formula type)`.
  - The earth checks take the earth model or the resistivity first.
  - A custom equivalent-earth reduction is admitted with `validate(reduction, expression)`.
  - PSCAD checks take `PSCADFormulation` instead of `Val(:pscad)`.
  - Result-space element types are checked with `validate(T, result space)`.
- `Commons.Expression(formula, operation, selectors...)` is one expression of a formula.
  Calling it evaluates `operation(formula, selectors..., functor, workspace)`.
- `Commons.Functor(formula, input, state)` stores the input of one evaluation point and the
  state of the formula, for every formula family.
- Formula methods take the formula object first. An expression's operation reads its input
  from `functor.input`.
- `Commons.formulas(family)` returns the identifiers that a formula family registers.
- Earth formulas and equivalent-earth reductions resolve one expression per pair with
  `Expression(formula, pair)`. A missing expression throws before the frequency loop.
- One `formulation_options(formula, expressions)` projection serves every formula family and
  reports an unused option as `unused formulation options (...) for :id`.
- `Formulation(tag; kwargs...)` builds the backend formulations `:coaxial`,
  `:cable_constants`, `:fem` and `:pscad`, and throws for an unknown tag.
- A saved `LineParameters` record keeps its `details` whole and loads them as data.
- On an earth with more than two layers, an earth formula with no expression beyond layer 2
  takes the `:default` reduction.
- `InternalImpedance` takes one formula per slot.
- The ShuntModel formulas are `:equivalent` and `:boundary`. `:default` resolves to
  `:equivalent`.
- Unified's `EarthAdmittance.source_coefficients` returns the axial and potential
  coefficients of one pair for both earth families.
- `ModalAnalysis` is a child module of `Engine`. The package root exports it and its public
  names.
- The Unified Γ is validated when the formula is built. The `:unified` docstring and the
  theory page on earth-return impedance state the validity ranges of the formula.
- Normalized author-formula identifiers to lowercase main-author-year symbols
  and made their compact descriptions readable author names.
- Renamed the dielectric constitutive choices to `:lossy` and `:lossless`,
  with `:default` retained as a routing alias. Ametani (2004) remains an
  application reference rather than the formula name.
- Registered the literature identities for the Schelkunoff internal-impedance,
  Ametani insulation-impedance, and Chrysochos modal routes while retaining
  their package-owned `:default` selections. Formula hook IDs are now directly
  queryable through `formula_id`.
- Refactored the `Material`-to-`compute` path around natural Julia promotion,
  owner-local validation, explicit definitions, immutable solver input, and scoped console
  logging.
- Cable parts and `CableDesign` now represent the common material reference state.
  Operating temperature is owned by the line problem and applied once to local
  resistivity values inside `compute`.
- Earth models now store static physical layers only. Analysis frequencies
  belong to the problem. Frequency-dependent soil laws and their ephemeral
  `EarthMaterial` values belong to `EarthProps` and are selected by the
  line-parameter formulation. The default relation is an exact static
  pass-through.
- Mathematical functions use short physical names without the former `calc_` prefix.
- Material and cable JSON files use an explicit versioned schema and dispatched
  type tags.
  JLS loading remains supported for trusted matching package types only.
- Wire-pattern searches return typed `WireEstimate` results, including ranked
  best-effort candidates for feasible search inputs that cannot meet every limit.
- Consolidated shared problem, formulation, and result roots plus the common
  action generics under `LineCableModels.Commons`.
- The coaxial and FEM backends share one matrix reduction,
  `Commons.reduce_line_matrices!`. Ideal transposition now averages the
  potential-coefficient matrix ``P`` before inverting it to ``Y`` in every backend.
  The coaxial engine previously averaged ``Y`` after the inversion.
- `ideal_transposition` now defaults to `false` for line-parameter and FEM
  formulations.
- Defined the vacuum constants once as `Commons.vacuum_permittivity(T)` and
  `Commons.vacuum_permeability(T)`. The FEM input writer, PSCAD import, and boundary
  shunt model now use the shared Float64 permittivity, one ulp below the former
  `8.8541878128e-12` literal.
- Replaced the prototype parameter and uncertainty paths with typed `Grid` and
  inferred `Gridspace` construction from the public declarative builders.
- Converged computation ownership around `Commons` action generics,
  `ParametricBuilder` deterministic traversal, `UQ` uncertainty propagation,
  and complete native execution in `Engine`.
- Reduced Gridspace to explicit finite `Grid` sources, local product or zip
  composition, one internal unresolved point, and recursive materialization or
  realization through concrete callable builders. Scalar public construction
  calls now invoke the corresponding scalar action.
- Replaced the declarative PlotBuilder renderer with a compact native Makie
  layer. `plot` retains automatic `(R, X, G, B)` views, scientific formatting,
  log correction, visible-series limits, uncertainty bars, widgets, previews,
  Monte Carlo plots, reports, and live SVG export while returning caller-owned
  figures, axes, legends, colorbars, and controls. Legends and colorbars now
  accept independent positions and native attributes. Observable quantities
  now determine axes, while `layout` only arranges those axes in arbitrary
  horizontal, vertical, or grid compositions.
- Restored `PlotBuilder` as a deliberately thin owner of optional plotting
  entry points and live `UIPlot` handles. Plot request normalization, preview
  presentation, and material palettes now live exclusively in the Makie
  extension. Engine and DataModel retain only scientific results, public
  observations, physical geometry, and property ranges.
- Split material visualization into the reusable `materialcolors(property,
  range)` palette and `materialscale!(position, scheme)` single-colorbar
  primitive. The three-scale reference remains optional preview sugar.
- Added `@observe` as three-index syntax over the native `observe` protocol.
- Added typed computation details with explicit higher-order retention and no
  default per-point or per-trial record collection.
- Moved human-facing XLSX line-parameter workbooks to ReportBuilder. The
  existing `export_data(:xlsx, ...)` call now delegates to
  `XLSXReportDefinition`.
- Routed ordinary and higher-order execution through `compute`, with
  `Combinatorial(inner)`, `LinearError(inner)`, and `MonteCarlo(inner)` selected
  explicitly.
- Named the concentric backend `LineCableModelsCoaxial`, retained
  `Formulation()` as the default line-parameter method bundle, and exposed
  optional trace data through `details(result)`.
- Replaced nominal author-formula types with stable literature symbols,
  overridable leaf routes, deterministic formula discovery, and an explicit
  longitudinal propagation-constant path.
- Ported nine legacy soil-dispersion laws into discovered
  `earth/frequencydependent/formulas/mainauthorYear.jl` files per literature
  identity, with SI conversion and typed parameters in place of MATLAB flags
  and runtime coefficient files.
- Moved equivalent homogeneous-earth rules to `EarthProps.EHEM`, separated
  before-FD and after-FD composition by dispatch, and evaluated the resulting
  material per conductor pair before the shared earth-impedance and admittance
  calculations. Added the conductivity-only `:martinsbritto2020` and complex
  propagation-constant `:xue2021` recurrences as discovered formulas.
- Added the public `formula(:AuthorYear; ...)` selector. Each formulation owner
  resolves the same wrapper locally, while EHEM order is selected with
  `order=:before` or `order=:after` without exposing sequence wrappers.
- Separated modal transformations from the coaxial backend into an independent
  problem-formulation-compute workflow. Modal results retain complete
  frequency-dependent voltage and current operators for reverse transformation.
- Made `build` the complete construction action for cable designs and systems.
  Explicit `Grid` inputs materialize the same action through `Gridspace`.
- Separated unresolved `*Definition` geometry from resolved primitives carrying
  absolute poses.
- Removed `NominalData` from `CableDesign`. `CablesLibrary` now binds optional
  named-tuple catalog records beside stored designs.
- Moved Measurements.jl and Distributions.jl integrations into package
  extensions.
- Restricted radial declarations to numeric radius or thickness semantics.
- Line-parameter plots and tables now select quantities with accessor tuples,
  such as `(R, L, G, C)` or `(abs, angle)`. Direct plotting of
  `LineParameters` returns R, X, G, and B matrix-dashboard pages in that order.
  A direct complex `Z` or `Y` coordinate request expands to a paired component
  page. Matrix coordinates identify semantic subplots, while legends identify
  overlaid result containers through `series_labels`.
- The PSCAD sources moved from the extensions directory to `src/pscad/`, and its remote
  runner project to `src/pscad/remote/`. `LineCableModels.PSCAD` is unchanged.
- `Commons` owns `initialize_buffers`, which builds every buffer of a computation.
  Each formula and shared component extends `Commons.initialize_buffers(selected, T,
  input, plan, buffers)` and shapes its own buffers, usually as a named tuple of arrays.
  The fourth argument, formerly `invariants`, is the computation plan.
- The coaxial workspace fields are `input`, `plan`, `buffers` and `trace`, formerly
  `invariants` and `capture` for the second and fourth.
- `SpectralIntegral` is a formulation, a subtype of `AbstractFormulation`, and
  `integrate` takes it first: `integrate(integral, Val(:quad), controls, buffers)`. The
  buffers are the record that `initialize_buffers(SpectralIntegral, Val(:quad), T, input,
  plan, buffers)` extends with `quadrature`.

### Removed

- Removed PlotBuilder pages, recipes, rendering contexts, backend registration,
  and fixed legend and colorbar docks. Plotting backends are selected per call or
  inherited from Makie's active backend.
- Removed the former `Commons` and `Utils` utility modules, package scalar-union aliases,
  coercion macros, operating-temperature cable fields, `EMTWorkspace`,
  intermediate-storage options, file logging, and the constructor proxy types `MaxFill`
  and `WireArray`.
- Retired unversioned and legacy JSON loading. The error identifies commit
  `a71bdfe1ac832f27a0c88b1d02596194aac46ec7` as the last snapshot able to migrate those
  files.
- Removed the former parameter tuple grammar, duplicate execution entrypoints,
  specialized analysis containers, and radial proxy wrapper types.
- Removed Grid identity and binding machinery, public temporary point records,
  traversal coordinates and metadata, result-side failure dictionaries, and
  passive forwarding definitions from the parametric construction path.
- Removed the `mode=:ZY`/`:RLCG` and `coord=:cart`/`:polar` keywords from
  line-parameter presentation.
- Removed the `InputValidation` module. `Commons` owns `validate`.
- Removed `Commons.check_core_result`. Use `validate(T, AbstractResultSpace)`.
- Removed `Engine.validate_modal_operators`. Modal operators subtype
  `AbstractModalOperators` and implement `size`.
- Removed `TextDisplay.quantity`, which no package code used. `show` of a
  `Units.Quantity` displays the quantity identity.
- Removed the `ParametricBuilder.Conductor`, `Insulator` and `Semiconductor` modules and
  their root exports. Each of their constructions forwarded to a public builder with a
  `tag`: `Conductor.Solid` to `core`, `Conductor.Shell` to `sheath`, `Conductor.Wires` to
  `wires`, `Conductor.Strip` to `tape` within `terminal`, `Conductor.Tubular` to `core`
  within `Group`, `Insulator.Shell` to `insulation` and `Semiconductor.Shell` to `screen`.
- Removed `Engine.integration_workspace`. `initialize_buffers(SpectralIntegral, Val(:quad),
  T, input, plan, buffers)` builds the quadrature buffers.

### Fixed

- `uncertainty` of a complex value or an array gives the uncertainty of each part or
  element, as `nominal` does.
- `uncertainty` of a deterministic number is the zero of its type. Values that are not
  numbers throw a `MethodError`.
- The methods that `TextDisplay.@showfields` defines record their call site, so `@which`,
  `methods` and stack traces show the caller.

## [0.2.0] - 2026-08-13

### Added

- Julia package extensions for the Makie backends.
- Aqua, SciML formatting, gitlint, clean-install, and modular documentation
  checks.
- Citation and contribution metadata.
- Type-stable core and statistical result containers.
- One declarative PlotBuilder renderer with interactive legends and
  one-click, non-overwriting SVG export.

### Changed

- Makie is an optional dependency. Core loading no longer imports Makie,
  CairoMakie, GLMakie, or WGLMakie.
- Plotting requires the caller to load CairoMakie, GLMakie, or WGLMakie
  explicitly.
- Documentation examples and development conventions were consolidated.
- Line-parameter results now include an explicit `:pul` or `:total` basis,
  and use `Z`, `Y`, `R`, `X`, `L`, `G`, `B`, and `C` accessors consistently.
- `Units` now maps physical accessors to quantity, unit, label, symbol,
  and scaling semantics without extracting values from result containers.
- `preview` and statistical plots return one `UIPlot`. Line-parameter plots
  return `Vector{UIPlot}`.

### Removed

- Binder experiments and Binder-specific notebook bootstrapping.
- The accidental standalone cable-construction subsystem. Parameter-space and
  covariance work remain under `ParametricBuilder`.
- FEM/Gmsh/GetDP support and sector-shaped cable support. The final
  pre-removal snapshot is `legacy/fem-sector` at commit
  `b75dd2723f90a83ec090b20605ea42af57f4a9c3`.
- The obsolete TODO scraper and duplicate tag-release workflow.
- The old `ResultsView`, `CableDesignMC`, `LineParametersPDF`, `plotmetadata`,
  duplicated plotting UIs, and direct Makie renderers.

### Migration

Load optional integrations explicitly:

```julia
using LineCableModels

using CairoMakie
# preview(...) and plot(...) are now available

R(line_parameters, 1, 1)       # complete frequency response
R(line_parameters, 1, 1, 2:5)  # selected frequencies
```

Projects that require the removed FEM or sector APIs must pin the archived
snapshot:

```sh
julia --project=. -e "using Pkg; Pkg.add(Pkg.PackageSpec(url=ARGS[1], rev=ARGS[2]))" https://github.com/Electa-Git/LineCableModels.jl.git b75dd2723f90a83ec090b20605ea42af57f4a9c3
```

## [0.1.0] - 2025-03-29

### Added

- Initial release.

[Unreleased]: https://github.com/Electa-Git/LineCableModels.jl/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/Electa-Git/LineCableModels.jl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/Electa-Git/LineCableModels.jl/releases/tag/v0.1.0
