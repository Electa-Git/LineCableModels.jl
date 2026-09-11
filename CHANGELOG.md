# Changelog

All notable changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and versions follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Refactored the `Material`-to-`compute` path around natural Julia promotion,
  owner-local validation, explicit definitions, immutable solver input, and scoped console
  logging.
- Cable parts and `CableDesign` now represent the common material reference state.
  Operating temperature is owned by the line problem and applied once to local
  resistivity values inside `compute`.
- Earth models now store static physical layers only. Analysis frequencies
  belong to the problem; frequency-dependent soil laws and their ephemeral
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
  action generics under `LineCableModels.Grammar`.
- Replaced the prototype parameter and uncertainty paths with typed `Grid` and
  inferred `Gridspace` construction from the public declarative builders.
- Converged calculation ownership around `Grammar` action generics,
  `ParametricBuilder` deterministic traversal, `UQ` uncertainty propagation,
  and complete native execution in `Engine`.
- Reduced Gridspace to explicit finite `Grid` sources, local product or zip
  composition, one internal unresolved point, and recursive materialization or
  realisation through concrete callable builders. Scalar public construction
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
  extension; Engine and DataModel retain only scientific results, public
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
  formula-owned typed functors, overridable leaf routes, deterministic formula
  discovery, and an explicit longitudinal propagation-constant path.
- Ported nine legacy soil-dispersion laws into one discovered
  `earthprops/fd/formulas/authoryear.jl` file per literature identity, with SI
  conversion and typed assumptions in place of MATLAB flags and runtime
  coefficient files.
- Moved equivalent homogeneous-earth rules to `EarthProps.EHEM`, separated
  before-FD and after-FD composition by dispatch, and evaluated the resulting
  material per conductor pair before the shared earth-impedance/admittance
  calculations. Added the conductivity-only `:MartinsBritto2020` and complex
  propagation-constant `:Xue2021` recurrences as discovered formulas.
- Added the public `formula(:AuthorYear; ...)` selector. Each formulation owner
  resolves the same wrapper locally, while EHEM order is selected with
  `order=:before` or `order=:after` without exposing sequence wrappers.
- Separated modal transformations from the coaxial backend into an independent
  problem–formulation–compute workflow. Modal results retain complete
  frequency-dependent voltage and current operators for reverse transformation.
- Made `build` the complete construction action for cable designs and systems.
  Explicit `Grid` inputs materialize the same action through `Gridspace`.
- Separated unresolved `*Definition` geometry from resolved primitives carrying
  absolute poses.
- Removed `NominalData` from `CableDesign`; `CablesLibrary` now binds optional
  named-tuple catalogue records beside stored designs.
- Moved Measurements.jl and Distributions.jl integrations into package
  extensions.
- Restricted radial declarations to numeric radius or thickness semantics.
- Line-parameter plots and tables now select quantities with accessor tuples,
  such as `(R, L, G, C)` or `(abs, angle)`. Direct plotting of
  `LineParameters` returns R, X, G, and B matrix-dashboard pages in that order.
  A direct complex `Z` or `Y` coordinate request expands to a paired component
  page. Matrix coordinates identify semantic subplots, while legends identify
  overlaid result containers through `series_labels`.

### Removed

- Removed PlotBuilder pages, recipes, rendering contexts, backend registration,
  and fixed legend/colorbar docks. Plotting backends are selected per call or
  inherited from Makie's active backend.
- Removed `Commons`, `Utils`, package scalar-union aliases, coercion macros,
  operating-temperature cable fields, `EMTWorkspace`, intermediate-storage options, file
  logging, and the constructor proxy types `MaxFill` and `WireArray`.
- Retired unversioned and legacy JSON loading. The error identifies commit
  `a71bdfe1ac832f27a0c88b1d02596194aac46ec7` as the last snapshot able to migrate those
  files.
- Removed the former parameter tuple grammar, duplicate execution entrypoints,
  specialised analysis containers, and radial proxy wrapper types.
- Removed Grid identity and binding machinery, public temporary point records,
  traversal coordinates and metadata, result-side failure dictionaries, and
  passive forwarding definitions from the parametric construction path.
- Removed the `mode=:ZY`/`:RLCG` and `coord=:cart`/`:polar` keywords from
  line-parameter presentation.

## [0.2.0] - 2026-08-13

### Added

- Julia package extensions for the Makie backends.
- Aqua, SciML formatting, gitlint, clean-install, and modular documentation
  checks.
- Citation and contribution metadata.
- Type-stable core and statistical result containers.
- A single declarative PlotBuilder renderer with interactive legends and
  one-click, non-overwriting SVG export.

### Changed

- Makie is an optional dependency. Core loading no longer imports Makie,
  CairoMakie, GLMakie, or WGLMakie.
- Plotting requires the caller to load CairoMakie, GLMakie, or WGLMakie
  explicitly.
- Documentation examples and development conventions were consolidated.
- Line-parameter results now carry an explicit `:pul` or `:total` basis,
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
