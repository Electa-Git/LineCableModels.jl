```@meta
EditURL = "https://github.com/Electa-Git/LineCableModels.jl/blob/master/CHANGELOG.md"
```

# Changelog

All notable changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and versions follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased](https://github.com/Electa-Git/LineCableModels.jl/compare/v0.2.0...HEAD)

## [0.2.0] - 2026-08-13

### Added

- Julia package extensions for Gmsh and the Makie backends.
- A manual, checksum-verified FEM release gate.
- Aqua, SciML formatting, gitlint, clean-install, and modular documentation
  checks.
- Citation and contribution metadata.
- Type-stable `CableConstants`, `SampleSummary`, `RLCG`, `CableConstantsMC`,
  `LineParametersMC`, and `HistogramPDF` result containers.
- A single declarative PlotBuilder renderer with interactive legends and
  one-click, non-overwriting SVG export.

### Changed

- Gmsh and Makie are optional dependencies. Core loading no longer imports
  Gmsh, Makie, CairoMakie, GLMakie, WGLMakie, or GetDP.
- The legacy `Engine.FEM` API delegates to the Gmsh extension and emits
  deprecation warnings.
- The compatible GetDP frontend is now a private, licensed source snapshot.
  GetDP itself must be supplied through `GETDP_EXECUTABLE` or `PATH`.
- Plotting requires the caller to load CairoMakie, GLMakie, or WGLMakie
  explicitly.
- Documentation examples and development conventions were consolidated.
- Line-parameter results now carry an explicit `:per_length` or `:total` basis,
  and use `Z`, `Y`, `R`, `X`, `L`, `G`, `B`, and `C` accessors consistently.
- `preview` and statistical plots return one `UIPlot`; line-parameter plots
  return `Vector{UIPlot}`.

### Removed

- Binder experiments and Binder-specific notebook bootstrapping.
- The accidental standalone `src/cablebuilder` subsystem. The maintained
  `ParametricBuilder.CableBuilder`, sector support, covariance work, and
  their tests remain.
- The hard dependency on the external, unregistered GetDP.jl package.
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

using Gmsh
using LineCableModels.Engine.FEM
formulation = FormulationSet(:FEM)
```

Set `GETDP_EXECUTABLE=/absolute/path/to/getdp` or place `getdp` on
`PATH` before invoking the transitional FEM solver.

## [0.1.0] - 2025-03-29

### Added

- Initial release.

[0.2.0]: https://github.com/Electa-Git/LineCableModels.jl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/Electa-Git/LineCableModels.jl/releases/tag/v0.1.0
