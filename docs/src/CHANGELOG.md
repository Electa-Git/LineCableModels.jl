```@meta
EditURL = "https://github.com/Electa-Git/LineCableModels.jl/blob/master/CHANGELOG.md"
```

# Changelog

All notable changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and versions follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased](https://github.com/Electa-Git/LineCableModels.jl/compare/v0.2.0...HEAD)

### Changed

- Refactored the eager `Material`-to-`compute` path around natural Julia promotion,
  owner-local validation, explicit definitions, immutable solver input, and scoped console
  logging.
- Cable parts and `CableDesign` now represent the common material reference state.
  Operating temperature is owned by the line problem and applied once to local
  resistivity values inside `compute`.
- Earth models now store static physical layers only. Analysis frequencies belong to the
  problem and frequency-dependent earth properties belong to the formulation.
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
  and complete analytical execution in `Engine`.
- Reduced Gridspace to explicit finite `Grid` sources, local product or zip
  composition, one internal unresolved point, and recursive materialization or
  realisation through concrete callable builders. Scalar public builder calls
  now construct eager domain values.
- Converged PlotBuilder on detached definition-owned pages and one fixed core
  sequence: entitle, parse, resolve, fetch, and finish. Loaded Makie extensions
  draw those pages directly through the standard shell.
- Added `@observe` as three-index syntax over the native `observe` protocol.
- Added typed computation details with explicit higher-order retention and no
  default per-point or per-trial record collection.
- Moved human-facing XLSX line-parameter workbooks to ReportBuilder. The
  existing `export_data(:xlsx, ...)` call now delegates to `XLSXReport`.
- Routed ordinary and higher-order execution through `compute`, with
  `Combinatorial(inner)`, `LinearError(inner)`, and `MonteCarlo(inner)` selected
  explicitly.
- Renamed the built-in backend to `Formulation(:analytical)` and selected
  parameter or trace output through formulation options.
- Made the declarative builder the default modelling API while retaining strict
  materialised constructors through explicit submodule imports.
- Moved Measurements.jl and Distributions.jl integrations into package
  extensions.
- Restricted radial declarations to numeric radius or thickness semantics.
- Line-parameter plots and tables now select quantities with accessor tuples,
  such as `(R, L, G, C)` or `(abs, angle)`. Direct plotting of
  `LineParameters` produces separate Z and Y figures with real and imaginary
  parts in side-by-side panels.

### Removed

- Removed `Commons`, `Utils`, package scalar-union aliases, coercion macros,
  operating-temperature cable fields, `EMTWorkspace`, primitive-storage options, file
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
- Type-stable primitive and statistical result containers.
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
- The accidental standalone `src/cablebuilder` subsystem. The maintained
  `ParametricBuilder.CableBuilder`, covariance work, and its tests remain.
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

[0.2.0]: https://github.com/Electa-Git/LineCableModels.jl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/Electa-Git/LineCableModels.jl/releases/tag/v0.1.0
