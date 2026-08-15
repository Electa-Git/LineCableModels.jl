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

- Julia package extensions for the Makie backends.
- Aqua, SciML formatting, gitlint, clean-install, and modular documentation
  checks.
- Citation and contribution metadata.
- Type-stable `CableConstants`, `SampleSummary`, `RLCG`, `CableConstantsMC`,
  `LineParametersMC`, and `HistogramPDF` result containers.
- A type-stable `ParametricSweep` container for ordered deterministic cases and
  their results, ready for later GridSpace integration.
- A single declarative PlotBuilder renderer with interactive legends and
  one-click, non-overwriting SVG export.

### Changed

- Makie is an optional dependency. Core loading no longer imports Makie,
  CairoMakie, GLMakie, or WGLMakie.
- Plotting requires the caller to load CairoMakie, GLMakie, or WGLMakie
  explicitly.
- Documentation examples and development conventions were consolidated.
- Line-parameter results now carry an explicit `:per_length` or `:total` basis,
  and use `Z`, `Y`, `R`, `X`, `L`, `G`, `B`, and `C` accessors consistently.
- `UnitHandler` now maps physical accessors to quantity, unit, label, symbol,
  and scaling semantics without extracting values from result containers.
- `preview` and statistical plots return one `UIPlot`; line-parameter plots
  return `Vector{UIPlot}`.

### Removed

- Binder experiments and Binder-specific notebook bootstrapping.
- The accidental standalone `src/cablebuilder` subsystem. The maintained
  `ParametricBuilder.CableBuilder`, covariance work, and its tests remain.
- FEM/Gmsh/GetDP integration and sector-shaped cable support. The final
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
